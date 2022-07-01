subroutine testQbeTermsCalculator

use flbe
use omp_lib
implicit none

integer L, serial, numThreads_saved, rtaCalc
logical slow, rta

	L = 4

	numThreads_saved = omp_get_max_threads()

	open(11, file='testQbeTermsCalculator.html')
	close(11, status='delete')
	open(12, file='testQbeTermsCalculator.tex')
	close(12, status='delete')


	call printTwice_openFile('testQbeTermsCalculator.log', ' [ test QbeTermsCalculator...' )
	
!-----------------------------------
	do serial=0,1
		slow = (serial.EQ.1)
		if (slow.OR.testParallel) then
			if (slow) then
				call printTwice( '  --- serial version ---     ', .TRUE.)
				call omp_set_num_threads( 1 )
			else
				call printTwiceSIS( CHAR(10)//'  --- OpenMP version (num.threads: ', &
&					numThreads, ') ---    ')
				call omp_set_num_threads( numThreads )
			end if

			call flbe_init ! init FFTW3 inner data for this number of threads

			do rtaCalc=0,1
				rta = (rtaCalc.EQ.1)
				if (rta) then
					call printTwice( '   === calculate relaxation times ===     ', .TRUE.)
				else
					call printTwice( '   === calculate dn/dt ===     ', .TRUE.)
				end if

				call testQTCwithParams( L, rta, slow )
			end do
			call flbe_done
		end if
	end do

	call printTwice_closeFile( CHAR(13)//' test QbeTermsCalculator: all OK ]' )

	call omp_set_num_threads( numThreads_saved )
	
end subroutine testQbeTermsCalculator






subroutine testQTCwithParams( L, rta, slow )
use flbe
use v0simple
implicit none

integer L
logical rta, slow

class(Geometry), pointer :: sizes
class(DispersionRelation), pointer :: dr, dr0, dr2
type(Subsystem), pointer :: sf,sb,sf2,sb2
type(Perturbation), pointer :: p
class(Hamiltonian), pointer :: hamilt
class(TableOfValues), pointer :: V_q

class(Efactor), pointer :: smooth
class(QbeTermsCalculator), pointer :: qtc1, qtc4, qtc3, qtc2

class(Occupations), pointer :: n, dn0,dn1,dn4
real*8, allocatable :: n0_symm(:,:,:,:), dn0_symm(:,:,:,:), V_q_symm(:,:,:)

class(Output), pointer :: oH, oT

real r
double precision t, diff01, diff14, diff04
integer kro, problemVersion, L2, kx,ky,kz, kroneckersPossible, s, V_q_possible, v_q_v
logical useKroneckers, good, V_q_dependent

integer timeStart, timeEnd, rate, maxVersion, i
double precision dt0, dt1, dt4, V0, de

	CALL system_clock(count_rate=rate)

	sizes => Geometry(L,L,L)

	dr => DispersionRelationParabolic( sizes )
	call dr % setCoeff( 1.d0 )

	dr0 => DispersionRelationConst( 0.d0 )
	dr2 => DispersionRelationParabolic( sizes )
	call dr2 % setCoeff( 2.d0 )

	smooth => Efactor( Efactor_EXACT, 1.d-8 )

	t = 0.d0


	L2 = sizes % L2(1)

	allocate( n0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( dn0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( V_q_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1) )
	V0 = 1.0d0
	V_q_symm (:,:,:) = V0

	do kz=-L2, L2-1
	do ky=-L2, L2-1
	do kx=-L2, L2-1
		do s=0, 1
			call RANDOM_NUMBER(r)
			n0_symm( kx,ky,kz, s ) = r
		end do

		V_q_symm( kx,ky,kz ) = V0/sqrt( 1.d0 + kx*kx + ky*ky + kz*kz )
	end do
	end do
	end do

	V_q => TableOfValues( sizes % L(:) )
	call symm2asymm( L, V_q_symm, V_q % values )

	sf => SubsystemFermi( dr )
	sf2 => SubsystemFermi( dr2 )
	de = dr % getEnergy( (/1, 0, 0/), sizes) - dr % getEnergy( (/0, 0, 0/), sizes)

	sb => SubsystemBose( dr )
	sb2 => SubsystemBose( dr )
	
	
	maxVersion = 6
	if (slow .OR. rta) maxVersion = 5

	do problemVersion=0, maxVersion
!	do problemVersion=5, 5

		kroneckersPossible = 0
		V_q_possible = 1
		select case(problemVersion)
			case(0)
				call printTwice( '    * Fermi-Fermi pair interaction (spin-up x spin-down) ')
				p => InteractionPairwise( sf, sf ) ! same Fermions with different spin projections
			case(1)
				call printTwice( '    * Bose-Bose pair interaction ')
				p => InteractionPairwise( sb ) ! Bose gas: a+a+ a a
				kroneckersPossible = 1
			case(2)
				call printTwice( '    * Bose-Bose2 pair interaction ')
				p => InteractionPairwise( sb, sb2 ) ! two Bose systems : a+a b+b
			case(3)
				call printTwice( '    * Bose31: a+a+a+ a')
				p => Perturbation( 'c3a1') ! general case
				call p % addSubsystem( sb, 3, 1)
				kroneckersPossible = 1
				V_q_possible = 0
			case(4)
				call printTwice( '    * Bose21: a+a+ a')
				p => Perturbation( 'c2a1') ! general case
				call p % addSubsystem( sb, 2, 1)
				kroneckersPossible = 1
				V_q_possible = 0
			case(5)
				call printTwice( '    * Fermi-Fermi2 pair interaction ')
				p => InteractionPairwise( sf, sf2 ) ! two Bose systems : a+a b+b
			case(6)
				call printTwice( '    * Bose21-Bose21: A+A+ A B+B+ B')
				p => Perturbation( '2121') ! general case
				call p % addSubsystem( sb, 2, 1)
				call p % addSubsystem( sb2, 2, 1)
				V_q_possible = 0
		end select

		do v_q_v=1-V_q_possible, 1 ! 0/1 or 1

			V_q_dependent = (v_q_v.EQ.0) ! first calculate V(q), then revert to constant
			if (V_q_dependent) then
				call printTwice( '      / V(q) q-dependent  ')
				call p % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
			else
				call printTwice( '      / V_0  q-independent  ')
				call p % setAmplitudeConst( V0 ) ! revert to constant
				
			end if
			do kro=0,kroneckersPossible
				useKroneckers = (kro.GT.0)

				hamilt => Hamiltonian()
				call hamilt % addPerturbation( p )

				n => Occupations( sizes % L, hamilt % numSubsystems )	
				dn0 => Occupations( n )
				dn4 => Occupations( n )
				do s=0, hamilt%numSubsystems-1
					call symm2asymm( L, n0_symm(:,:,:,s), n % values(:,:,:,s) )
				end do
		
				call printTwiceSL( '      // useKroneckers: ', useKroneckers)

				qtc4 => QbeTermsCalculatorFFT( sizes, hamilt, smooth, useKroneckers, de )
!				qtc4 => QbeTermsCalculatorSimple( sizes, hamilt, smooth, useKroneckers, VERSION_AS_IS_v2 )

				dn4 % values(:,:,:,:) = 0.d0			
			    call SYSTEM_CLOCK( timeStart )
				call qtc4 % addValue( t, n, dn4, rta)
			    call SYSTEM_CLOCK( timeEnd )
			    dt4 = (timeEnd-timeStart+0.d0)/rate
				if (dt4 .GT. 2.d0) write(*, '(A,F8.3,A)') '       v4 in ', dt4, ' seconds       '

			    call SYSTEM_CLOCK( timeStart )
				select case(problemVersion)
					case(0)
						call v0_f_c2a2( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta )
					case(1)
						call v0_b_c2a2( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta, useKroneckers )
					case(2)
						call v0_bb2_c1a1c1a1( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta )
					case(3)
						call v0_b_c3a1( L, n0_symm, dn0_symm, V0, rta, useKroneckers )
					case(4)
						call v0_b_c2a1( L, n0_symm, dn0_symm, V0, rta, useKroneckers )
					case(5)
						call v0_ff2_c1a1c1a1( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta )
					case(6)
						call v0_bb2_c2a1c2a1( L, n0_symm, dn0_symm, V0, rta, useKroneckers ) 
				end select
			    call SYSTEM_CLOCK( timeEnd )
			    dt0 = (timeEnd-timeStart+0.d0)/rate
				if (dt0 .GT. 2.d0) write(*, '(A,F8.3,A)') '       v0 in ', dt0, ' seconds       '

				call printTwiceSFSFS('       v0/v4 in ', dt0, ' /', dt4, ' seconds       ')

				do s=0, hamilt%numSubsystems-1
					call symm2asymm( L, dn0_symm(:,:,:,s), dn0 % values(:,:,:,s) )
				end do
			
				diff04 = dn0 % compareToRelative( dn4 )
			
				good = (abs(diff04).LT.1d-14)
				if (good) then
					call printTwice( '        OK  ')
					write(debugFileNumber, *) '      difference: v0-v4: ', diff04, '               '
				else
					call printTwice( '        FAIL')
					write(*, *) '      difference: v0-v4: ', diff04, '               '

					write(debugFileNumber, *) '      dn_v0=', dn0 % values(:,0,0,:)
					if (verbose) then
						do i=0, hamilt % numSubsystems-1
							write(*, *) '      v4(subs)=', dn4 % values(:,0,0,i)
						end do
						do i=0, hamilt % numSubsystems-1
							write(*, *) '      v0(subs)=', dn0 % values(:,0,0,i)
						end do
						do i=0, hamilt % numSubsystems-1
							write(*, *) '      v0-v4(subs)=', dn0 % values(:,0,0,i)-dn4 % values(:,0,0,i)
						end do
					end if

					call printTwice( '      see equations in testQbeTermsCalculator.html and .tex')

					qtc1 => QbeTermsCalculatorSimple( sizes, hamilt, smooth, useKroneckers, VERSION_AS_IS_v1 )
					qtc2 => QbeTermsCalculatorSimple( sizes, hamilt, smooth, useKroneckers, VERSION_AS_IS_v2 )
					qtc3 => QbeTermsCalculatorSimple( sizes, hamilt, smooth, useKroneckers, VERSION_AS_IS_v3 )
	
					oH => OutputHTML( 'testQbeTermsCalculator.html', debugFileNumber+1 )	
					oT => OutputTeX( 'testQbeTermsCalculator.tex', debugFileNumber+2 )

					call hamilt % show( oH )
					call qtc1 % show( oH )
					call qtc2 % show( oH )
					call qtc3 % show( oH )
					call qtc4 % show( oH )

					call hamilt % show( oT )
					call qtc1 % show( oT )
					call qtc2 % show( oT )
					call qtc3 % show( oT )
					call qtc4 % show( oT )

					deallocate( oH )
					deallocate( oT )
					deallocate( qtc3 )
					deallocate( qtc2 )
				
					write(debugFileNumber, *) '       dn_v4=', dn4 % values(:,0,0,:)

					if (verbose) then
						! comparison failed: check also with earlier version v1
						print*, '        calculating again with v1...'

						dn1 => Occupations( n )
						dn1 % values(:,:,:,:) = 0.d0
					    call SYSTEM_CLOCK( timeStart )
						call qtc1 % addValue( t, n, dn1, rta )
				    	call SYSTEM_CLOCK( timeEnd )
					    dt1 = (timeEnd-timeStart+0.d0)/rate
						if (dt1 .GT. 2.d0) write(*, '(A,F8.3,A)') 'v1 in ', dt1, ' seconds  '
	
						diff01 = dn0 % compareToRelative( dn1 )
						diff14 = dn1 % compareToRelative( dn4 )
			    
						write(debugFileNumber, *) '       dn_v1=', dn1 % values(:,0,0,:)
						call printTwiceSFSFS( '       difference: v0-v1: ', diff01, ' v1-v4: ', diff14, '  ' )
						deallocate( dn1 )

					end if ! verbose
					deallocate( qtc1 )
				
					stop 'see details in testQbeTermsCalculator.log/html/tex'
				end if ! good
			
				deallocate( qtc4 )
				deallocate( dn4 )
				deallocate( dn0 )
				deallocate( n )
				deallocate( hamilt )
			end do ! kroneckers
		end do ! V_q
		deallocate( p )
	end do
	deallocate( sb2 )
	deallocate( sb )
	deallocate( sf2 )
	deallocate( sf )
	deallocate( V_q )
	deallocate( V_q_symm )
	deallocate( n0_symm )
	deallocate( dn0_symm )
	deallocate( smooth )
	deallocate( dr )
	deallocate( sizes )

	if (slow) call printTwice( '    * Bose21-Bose21: skipped in serial version (takes too long)' )

end subroutine testQTCwithParams
