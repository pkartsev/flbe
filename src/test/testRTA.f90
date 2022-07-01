subroutine testRTA

use omp_lib
use flbe
implicit none

integer L, problemVersion, numThreads_saved, kroPossible

logical useKroneckers

	L=4

	numThreads_saved = omp_get_max_threads()

	call omp_set_num_threads( numThreads )
	call flbe_init ! init FFTW3 inner data for this number of threads
	
	call printTwice_openFile( 'testRTA.log', ' [ test RTA...'//CHAR(10) )
	
	do problemVersion=0,2 !0,2 TODO
		call testRTAforParameters( L, problemVersion )
	end do
	call printTwice_closeFile( CHAR(13)//' test RTA: all OK ]               ' )
	
	call flbe_done
	call omp_set_num_threads( numThreads_saved )
	

end subroutine testRTA


subroutine testRTAforParameters( L, problemVersion )

use flbe
use v0simple
implicit none 

integer L, problemVersion

double precision, parameter :: deltaN=1d-5
double precision, parameter :: requiredPrecision=1d-7

double precision, parameter :: density = 1.87d0
double precision, parameter :: kT = 4.d0

class(Geometry), pointer :: sizes
class(Problem), pointer :: prob
class(BoltzmannEquation), pointer :: be1, be4
class(Occupations), pointer :: n, dn, rtaCoeffs
class(TableOfValues), pointer :: V_q
class(Perturbation), pointer :: p
class(Subsystem), pointer :: s1,s2

class(Output), pointer :: oH, oT

real*8, allocatable :: n0_symm(:,:,:,:), dn0_symm(:,:,:,:), V_q_symm(:,:,:)

double precision oldvalue
integer L2, kx,ky,kz, i, j, k, v_q_v, kro, kroPossible
double precision vPlus, vMinus, coeff0, coeff1, t, de
logical V_q_dependent, useKroneckers
class(DispersionRelation), pointer :: dr

	sizes => Geometry( L,L,L )

	L2 = L/2
	allocate( n0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( dn0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( V_q_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1) )

	do kz=-L2, L2-1
	do ky=-L2, L2-1
	do kx=-L2, L2-1
		V_q_symm( kx,ky,kz ) = 1.d0/sqrt( 1.d0 + kx*kx + ky*ky + kz*kz )
	end do
	end do
	end do
	V_q => TableOfValues( sizes % L(:) )
	call symm2asymm( L, V_q_symm, V_q % values )

	kroPossible = 0
	if (problemVersion.EQ.1) kroPossible = 1
	do kro=0,kroPossible

		do v_q_v=0, 1
			V_q_dependent = (v_q_v.EQ.0)

			oH => OutputHTML( 'testRTA.html', debugFileNumber+1 )
			call oH % markToDelete(.TRUE.)
			oT => OutputHTML( 'testRTA.tex', debugFileNumber+2 )	
			call oT % markToDelete(.TRUE.)

			select case (problemVersion)
				case(0)
					call printTwice( '  = Fermi gas =               ' )
					prob => InteractingFermiGas( sizes, .FALSE., particleDensity=density, kT=kT )
				case(1)
					call printTwice( '  = Bose gas  =               ' )
					prob => InteractingBoseGas( sizes, particleDensity=density, kT=kT )
				case(2)
					call printTwice( '  = Bose-Bose2 pair interaction =                 ' )
					prob => BoseGas( sizes, 2, particleDensity=density, kT=kT ) ! no interaction yet; 2 subsystems
					s1 => prob % hamilt % subsystems(0) % obj
					s2 => prob % hamilt % subsystems(1) % obj
					p => InteractionPairwise( s1, s2 )
					call prob % hamilt % addPerturbation( p )
			end select
			dr => prob % hamilt % subsystems(0) % obj % dispersionLaw
			de = dr % getEnergy( (/1, 0, 0/), sizes) - dr % getEnergy( (/0, 0, 0/), sizes)

			useKroneckers = (kro.GT.0)
			call printTwiceSL('   == useKroneckers: ', useKroneckers)
			prob % useKroneckers = useKroneckers


			if (V_q_dependent) then 
				call printTwice( '    === V(q) q-dependent' )
				call prob % hamilt % perturbations(0) % obj % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
			else 
				call printTwice( '    === V(0) q-independent' )
				call prob % hamilt % perturbations(0) % obj % setAmplitudeConst(1.d0)
			end if
		
!			be4 => BoltzmannEquation( prob, de, VERSION_AS_IS_v3  ) ! TODO: temporary not FFT
			be4 => BoltzmannEquation( prob, de )
			be1 => BoltzmannEquation( prob, de, VERSION_AS_IS_v1 )
		
			if (verbose) call showRtaDebug( oH, oT, prob, be1, be4)
		

			n => Occupations( sizes % L, prob % numS )
			dn => Occupations( n )
			rtaCoeffs => Occupations( n )

			! 1) take equilibrium occupations
			if (.NOT. prob % equilibriumValid) then
				call prob % updateEquilibriumOccupations
			end if
			n % values(:,:,:,:) = prob % equilibriumOccupations % values(:,:,:,:)

			! 2) calculate tau from RTA using BE calculator
			call be4 % getCoeffsForRTA( rtaCoeffs )

			do kz=0,L2
			do ky=kz,L2
			do kx=ky,L2

				coeff0 = rtaCoeffs % values(kx,ky,kz,0)
			
				t = 0.d0

				oldvalue = n % values(kx,ky,kz,0)
				n % values(kx,ky,kz,0) = oldvalue*(1 + deltaN)
				call be4 % getRightPart( t, n, dn )
				vPlus  = dn % values(kx,ky,kz,0)

				n % values(kx,ky,kz,0) = oldvalue*(1 - deltaN)
				call be4 % getRightPart( t, n, dn )
				vMinus = dn % values(kx,ky,kz,0)

				n % values(kx,ky,kz,0) = oldvalue

				coeff1 = (vPlus-vMinus)/(2*deltaN*oldvalue)

				write(debugFileNumber, '(A,3I2,A,ES14.8,A,2ES17.8,A,ES20.10)') &
&					' k=', kx,ky,kz, ' nk=', oldvalue, ' coeff(k)=', coeff0, coeff1, &
&					' rel. difference=', coeff0/coeff1-1.d0
				if (verbose) write(*, '(A,3I2,A,ES14.8,A,2ES17.8,A,ES20.10)') &
&					' k=', kx,ky,kz, ' nk=', oldvalue, ' coeff(k)=', coeff0, coeff1, &
&					' rel. difference=', coeff0/coeff1-1.d0

				if (abs(coeff0/coeff1 - 1.d0).GT.requiredPrecision) then
					call printTwice( ' FAIL' )

					! if not shown already
					if (.NOT.verbose) call showRtaDebug( oH, oT, prob, be1, be4)
				
					call oH % markToDelete(.FALSE.)
					call oT % markToDelete(.FALSE.)
					deallocate( oH )
					deallocate( oT )
				
					stop 'see testRTA.log/html/tex for details'
				end if

			end do
			end do
			end do ! kx,ky,kz
			deallocate( be1 )
			deallocate( be4 )
			deallocate( dn )
			deallocate( rtaCoeffs )
			deallocate( n )
			deallocate( prob )
			deallocate( oH )
			deallocate( oT )
		end do ! v_q
	end do ! useKroneckers
	deallocate( V_q )
	deallocate( V_q_symm )
	deallocate( sizes )

end subroutine testRTAforParameters

subroutine showRtaDebug( oH, oT, prob, be1, be4 )
use flbe
class(Output), pointer :: oH, oT
class(Problem), pointer :: prob
class(BoltzmannEquation), pointer :: be1, be4

	call prob % hamilt % show( oH )
	call be1 % qtc % show( oH )
	call be4 % qtc % show( oH )
	call prob % hamilt % show( oT )
	call be1 % qtc % show( oT )
	call be4 % qtc % show( oT )

end subroutine showRtaDebug
