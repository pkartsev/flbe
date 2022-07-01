include '../flbe/flbe.f90'


use flbe 
implicit none

class(DispersionRelation), pointer :: dr2
type(Subsystem), pointer ::  sf, sf2, sb, sb2, sb3
type(Perturbation), pointer :: v_bb, v_ee, v_21, v_fb, v_eph, v_21b
type(TableOfValues), pointer :: V_q 
type(Hamiltonian), pointer :: h
type(Geometry), pointer :: g
type(Problem), pointer :: prob

class(BoltzmannEquation), pointer :: be1, be2, be3, be4
class(FiniteDifferenceSolver), pointer :: fd
type(Efactor), pointer :: F
class(Output), pointer :: oH, oT

integer LX, LY, LZ
double precision dv, t, de

integer iTimes1, iTimes2, rate
real cpuTimeDelta

	CALL system_clock(count_rate=rate)

! openmp number of threads
	call omp_set_num_threads( 16 )

	call flbe_init

! sizes
!	g => Geometry(2,2,2)
!	g => Geometry(4,4,4)
!	g => Geometry(6,6,6)
!	g => Geometry(8,8,1)
!	g => Geometry(6,6,1)
!	g => Geometry(4,8,1)
!	g => Geometry(10,10,1)
	g => Geometry(4,4,1)
!	g => Geometry(8,8,8)
!	g => Geometry(10,10,10)
!	g => Geometry(12,12,12)

	LX = g % L(1)
	LY = g % L(2)
	LZ = g % L(3)

! energies
!	dr2 => DispersionRelationParabolic( g )
!	dr2 => DispersionRelationLinear( g )
	dr2 => DispersionRelationParabolic(  )
	call dr2 % setCoeff(1.0d0)

	de = dr2 % getEnergy( (/1, 0, 0/), g) - dr2 % getEnergy( (/0, 0, 0/), g)
!	print*, ' de=', de, ' energies=', dr2 % getEnergy( (/1, 0, 0/), g),  dr2 % getEnergy( (/0, 0, 0/), g)

! subsystems
	sf => SubsystemFermi( dr2 )
	sf2 => SubsystemFermi( dr2 )
	sb  => SubsystemBose( dr2 )
	sb2 => SubsystemBose( dr2 )
	sb3 => SubsystemBose( dr2 )


! perturbation
	v_ee => InteractionPairwise( sf, sf ) ! spin-up x spin-down
!	v_ee => InteractionPairwise( sf, sf, .TRUE. ) ! same spins
	v_bb => InteractionPairwise( sb ) ! same
!	v_bb => InteractionPairwise( sb, sb2 ) ! not same
	v_eph => InteractionWithPhonons( sf, sb )
!	v_eph => InteractionWithPhonons( sb2, sb )
	v_fb => InteractionPairwise( sf, sb )
!	v_fb => InteractionPairwise( sb, sf )
!	v_fb % title = 'pairFB'
	v_21 => Perturbation( '2-1') ! general case
	v_21b => Perturbation( '1-2') ! general case
!	call v_21 % addSubsystem( sf, 1, 1, .TRUE.)
!	call v_21 % addSubsystem( sf, 2, 0, .FALSE.)

!	call v_21 % addSubsystem( sb, 2, 0)
!	call v_21 % addSubsystem( sf2, 2, 0)
!	call v_21 % addSubsystem( sb, 3, 1)
!	call v_21 % addSubsystem( sb, 1, 2)
!	call v_21 % addSubsystem( sb2, 1, 2)
	call v_21 % addSubsystem( sb, 2, 1)
!	call v_21 % addSubsystem( sb2, 2, 1)
!	call v_21b % addSubsystem( sb, 1, 3)
!	call v_21 % addSubsystem( sb2, 1, 2)
 
!	call v_21 % addSubsystem( sf2, 2, 0)
!	call v_21 % addSubsystem( sb, 0, 2)
!	call v_21 % addSubsystem( sf, 0, 2)
!	call v_21 % addSubsystem( sf2, 2, 0)
	
! interaction amplitude V(q)
	V_q => TableOfValues( g % L(:) )
	V_q % values = 1.d0
!	V_q % values(   2,0,0) = 30.d0
!	V_q % values(LX-2,0,0) = 30.d0
	V_q % values(   1,0,0) = 5.d0
	V_q % values(LX-1,0,0) = 5.d0
!	V_q % values(0,   1,0) = 20.d0
!	V_q % values(0,LY-1,0) = 20.d0
!	V_q % values(0,0,   1) = 10.d0
!	V_q % values(0,0,LZ-1) = 10.d0
!	call v_eph % setAmplitudeTable( V_q, (/0, 2/), (/1, -1/) )
!	call pump % setAmplitudeTable( V_q, (/0, -1/), (/1, -1/) )
!	call v_ee % setAmplitudeTable( V_q )
!	call v_fb % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
!	call v_ee % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
!	call v_bb % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
!	call v_21 % setAmplitudeTable( V_q, (/0, 3/), (/1, 2/) )
!	call v_21 % setAmplitudeTable( V_q, (/0, 3/), (/2, 1/) )
!	call v_21 % setAmplitudeTable( V_q, (/2, -1/) )
!	print*, ' V_q size: ', SIZE(v_q % values(:,0,0)), SIZE(v_q % values(0,:,0)), SIZE(v_q % values(0,0,:))
!	print*, ' L=', g % L


! hamiltonian
	h => Hamiltonian()
!	call h % addSubsystem( sb )
!	call h % addPerturbation(v_eph)
!	call h % addPerturbation(v_ee)
	call h % addPerturbation(v_bb)
!	call h % addPerturbation(v_fb)
!	call h % addPerturbation(v_21)
!	call h % addPerturbation(v_21b)

!	print*, ' h=', TRIM( h % toString() )
!	print*, ' g=', g % toString()


	F => EfactorExact( 1.d-8 )
!	F => EfactorGauss( 0.1d0 )
!	F => EfactorLorentz( 1.d-1 )

	allocate(prob)
	call prob % init( g, h, F )
	
	call prob % setTemperature( 1.d0 )
	call prob % setChemicalPotential( -1.d0 ) ! same value for all subsystems

!	prob => ElectronGasWithPhonons( g, 0.9d0, 1.d0, 1.d0 )
!	prob => InteractingBoseGas( g, 0.9d0, 1.d0, 1.d0 )
!	prob % F =>	F

	prob % useKroneckers = .FALSE.
	
! equation
	be1 => BoltzmannEquation( prob, de, VERSION_AS_IS_v1 )
	be2 => BoltzmannEquation( prob, de, VERSION_AS_IS_v2 )
	be3 => BoltzmannEquation( prob, de, VERSION_AS_IS_v3 )
	be4 => BoltzmannEquation( prob, de, VERSION_FAST )

	oH => OutputHTML( '1.html', 11 )
!	call oH % test
	call prob % show( oH )
	call be1 % show( oH )
	call be2 % show( oH )
	call be3 % show( oH )
	call be4 % show( oH )
	deallocate( oH )

!stop

	oT => OutputTeX( '1.tex', 12 )
!	call oT % test
	call prob % show( oT )
	call be1 % show( oT )
	call be2 % show( oT )
	call be3 % show( oT )
	call be4 % show( oT )
	deallocate( oT )	
	
!	stop

! solver
	fd => FD_Euler( be4 )
!	fd => FD_DormandPrince( be4 )
	
!	t = 0.d0
!	call prob % updateEquilibriumOccupations
!	fd % n % values(:,:,:,:) = prob % equilibriumOccupations % values(:,:,:,:)
!	open(18, file='v3_bose_qbe_3d.dat', form='unformatted')
!	open(18, file='v3_fermi_qbe_3d.dat', form='unformatted')
!	write(18) fd % n % values
!	close(18)
!	stop 'dat saved'


!	call fd % qbe % getRightPart( t, fd % n, fd % k1 )
!	open(18, file='v3_bose_qbe_3d.dat', form='unformatted')
!	write(18) fd % k1 % values
!	close(18)
!	stop 'dat saved'

	fd % n % values(:,:,:,:) = 0.1d0
	fd % n % values(1,0,0,:) = 0.2d0
	fd % n % values(LX-1,0,0,:) = 0.2d0
	fd % n % values(0,1,0,:) = 0.2d0
	fd % n % values(0,LY-1,0,:) = 0.2d0
	if (fd % n % N4.GT.1) then
		fd % n % values(:,:,:,1) = 0.3d0
		fd % n % values(1,0,0,1) = 0.74d0
		fd % n % values(LX-1,0,0,1) = 0.74d0
	end if
	if (fd % n % N4.GT.2) then
!		fd % n % values(:,:,:,2) = 0.56d0
!		fd % n % values(1,0,0,2) = 0.1245d0
	end if

	print*, ' comparison of be1 and be2 getRightPart()'
	dv = be1 % compareTo( be2, 0.d0, fd % n )
!	dv = be2 % compareTo( be1, 0.d0, fd % n )
!	print*, ' comparison of be1 and be4 getRightPart()'
	write(*, *) ' relative difference = ', dv
!stop

	print*, ' comparison of be2 and be3 getRightPart()'
	dv = be2 % compareTo( be3, 0.d0, fd % n )
	write(*, *) ' relative difference = ', dv

!    call SYSTEM_CLOCK(iTimes1)
!	print*, ' comparison of be1 and be3 getRightPart()'
!	dv = be1 % compareTo( be3, 0.d0, fd % n )
	print*, ' comparison of be3 and be4 getRightPart()'
	dv = be3 % compareTo( be4, 0.d0, fd % n )
	write(*, *) ' relative difference = ', dv
!	call SYSTEM_CLOCK(iTimes2)
!	cpuTimeDelta = real(iTimes2-iTimes1)/real(rate)
!	write(*, '(A,F8.2,A)') ' DONE in ', cpuTimeDelta, ' seconds.'
	
	deallocate(fd)

	deallocate(be4)
	deallocate(be3)
	deallocate(be2)
	deallocate(be1)

	deallocate(prob)
	deallocate(F)
	deallocate(h)
	deallocate(V_q)
	deallocate(V_bb)

	deallocate(sb)
	deallocate(dr2)
	deallocate(g)

	call flbe_done

end program
