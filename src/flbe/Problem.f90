module Problem_module

	use Output_module
	use Efactor_module
	use Geometry_module
	use Subsystem_module
	use Perturbation_module
	use Occupations_module
	use Hamiltonian_module
	
	implicit none
	
	type Problem
		type(Geometry), pointer :: sizes
		type(Hamiltonian), pointer :: hamilt
		type(Efactor), pointer :: F
		integer numS ! number of subsystems
		logical useKroneckers
		logical hamiltonianCreatedInside
		
		double precision temperature
		double precision, allocatable :: chemicalPotential(:) ! for each subsystem

		! work data
		type(Occupations), pointer :: equilibriumOccupations
		logical equilibriumValid
		
	contains
		procedure :: init => Problem_init_sub
		procedure :: setTemperature => Problem_setTemperature
		procedure :: setChemicalPotential => Problem_setChemicalPotential
		procedure :: setChemicalPotentialForNumberOfParticles => Problem_setChemicalPotentialForNumberOfParticles
		procedure :: updateEquilibriumOccupations => Problem_updateEquilibriumOccupations
		procedure :: updateN0ForSubsystem => Problem_updateN0ForSubsystem
		
		procedure :: show => Problem_show
		
		final :: Problem_destructor
	end type  Problem
	
	interface Problem
		module procedure Problem_init
	end interface Problem
	interface BoseGas ! no interaction (to add later) , several Bose subsystems (if needed)
		module procedure BoseGas_init
	end interface BoseGas
	interface InteractingBoseGas ! Bose gas with contact interaction, one subsystem
		module procedure InteractingBoseGas_init
	end interface InteractingBoseGas
	interface InteractingFermiGas ! Fermi gas with contact interaction (n_up x n_down)
		module procedure InteractingFermiGas_init
	end interface InteractingFermiGas
	interface ElectronGasWithPhonons
		module procedure ElectronGasWithPhonons_init
	end interface ElectronGasWithPhonons
	
contains

!----------------------------------------------------------------

	function Problem_init( sizes, hamilt, F, kT, mu ) result (this)
	type( Problem ), pointer :: this

	type( Geometry ), pointer :: sizes
	type( Hamiltonian ), pointer :: hamilt
	type( Efactor ), pointer :: F
	double precision, optional :: kT, mu

		allocate(this)
		call this % init( sizes, hamilt, F, kT, mu )

	end function Problem_init



	subroutine Problem_init_sub( this, sizes, hamilt, F, kT, mu )
	class( Problem ), intent(INOUT) :: this

	type( Geometry ), pointer :: sizes
	type( Hamiltonian ), pointer :: hamilt
	type( Efactor ), pointer :: F
	double precision, optional :: kT, mu

		this % sizes => sizes
		this % hamilt => hamilt
		this % hamiltonianCreatedInside = .FALSE.
		this % F => F
		this % numS = hamilt % numSubsystems

		this % useKroneckers = .TRUE. ! default value
		
		this % equilibriumOccupations => Occupations( this % sizes % L(:), this % hamilt % numSubsystems )
		allocate( this % chemicalPotential( 0: this%numS-1 ) )

		if (present(kT)) then
			this % temperature = kT
		else
			this % temperature = 1.234d0 ! default value
		end if
		
		if (present(mu)) then
			this % chemicalPotential(:) = mu
			call this % updateEquilibriumOccupations()
		else
			this % chemicalPotential(:) = 0.d0 ! invalid ; expecting to set the value later
			this % equilibriumValid = .FALSE.
		end if
		
	end subroutine Problem_init_sub


	subroutine Problem_destructor( this )
	type( Problem ) :: this
	
		deallocate( this % chemicalPotential )
		deallocate( this % equilibriumOccupations )
		if (this % hamiltonianCreatedInside) then
			deallocate( this % hamilt )
		end if

	end subroutine Problem_destructor

! -- main functions

	subroutine Problem_setTemperature( this, v )
	class( Problem ), intent(INOUT) :: this
	double precision v
	
		this % temperature = v
		this % equilibriumValid = .FALSE.
		
	end subroutine Problem_setTemperature

	subroutine Problem_setChemicalPotential( this, mu, s )
	class( Problem ), intent(INOUT) :: this
	integer, optional :: s
	double precision mu

		if (present(s)) then	
			this % chemicalPotential(s) = mu
		else
			this % chemicalPotential(:) = mu
		end if
		this % equilibriumValid = .FALSE.
		
	end subroutine Problem_setChemicalPotential



	subroutine Problem_updateEquilibriumOccupations( this )
	class( Problem ), intent(INOUT) :: this
	
	integer ss
	
		do ss=0, this%numS-1
!			print*, ' updateEO: ss=', ss
			call this % updateN0ForSubsystem( ss )
		end do
		
		this % equilibriumValid = .TRUE.
		
!		print*, ' Problem_updateEquilibriumOccupations'
		
	end subroutine Problem_updateEquilibriumOccupations


	subroutine Problem_updateN0ForSubsystem( this, ss )
	class( Problem ), intent(INOUT) :: this
	integer ss	

	integer volume, k3d, k(1:3)
	type( Subsystem ), pointer :: s
	type( DispersionRelation ), pointer :: dr
	double precision e,n0,kT,mu

			kT = this % temperature
			volume = this%sizes%volume

			s => this % hamilt % subsystems(ss) % obj
!			if (.NOT. s % lifetimePresent) then
!				this % equilibriumOccupations % values(:,:,:,:) = -1024.d0 ! for testing: not expected to use
!			else
				dr => s % dispersionLaw
				mu = this % chemicalPotential(ss)
!$omp parallel do private(k3d,k,e,n0)
				do k3d=0, volume-1
					call this % sizes % toVector( k3d, k )
					e = dr % getEnergy( k, this % sizes )
					n0 = s % getEquilibriumOccupation( e, kT, mu )
					call this % equilibriumOccupations % setValue( k, ss, n0 )
				end do
!$omp end parallel do
!			end if
!            print*, ' done'
	end subroutine Problem_updateN0ForSubsystem


	subroutine Problem_setChemicalPotentialForNumberOfParticles( this, ss, aim )
	class( Problem ), intent(INOUT) :: this
	integer ss
	double precision aim
	
	type( Subsystem ), pointer :: s
		
	double precision a, b, c, fa, fb, fc
	integer iter
	
		s => this % hamilt % subsystems(ss) % obj
		if (s % boseStatistics) then
			a = -1d-12
			b = -1d4 ! expecting large enough (abs. value)
		else
			a = -1d4 ! somewhere below conduction band
			b = 1d4 ! expecting large enough
		end if
		this % chemicalPotential(ss) = a
		call this % updateN0ForSubsystem(ss)
		fa = sum( this % equilibriumOccupations % values(:,:,:,ss) ) - aim
		this % chemicalPotential(ss) = b
		call this % updateN0ForSubsystem(ss)
		fb = sum( this % equilibriumOccupations % values(:,:,:,ss) ) - aim
		if (fa*fb.GT.0) stop 'bad range for bisection'
		do iter=1, 52
			c = (a+b)/2
			this % chemicalPotential(ss) = c
			call this % updateN0ForSubsystem(ss)
			fc = sum( this % equilibriumOccupations % values(:,:,:,ss) ) - aim
			if (fa*fc.LT.0) then
				b = c
			else
				a = c
				fa = fc
			end if
		end do
		this % chemicalPotential(ss) = c
		this % equilibriumValid = .FALSE.

	end subroutine Problem_setChemicalPotentialForNumberOfParticles
	








! -- for debug

	subroutine Problem_show( this, o )
	class( Problem ) :: this
	class( Output ), pointer :: o
	
		call o % horizontalLine
		call o % newSection( 'Problem' )
		call this % hamilt % show( o )
		call this % sizes % show( o )
		call this % F % show( o )
	
	end subroutine Problem_show

!------------------------------------------------------------
	
	function InteractingBoseGas_init( sizes, chemicalPotential, particleDensity, kT, U0 ) result(this)
	type( Geometry ), pointer :: sizes
	double precision, optional :: chemicalPotential
	double precision, optional :: particleDensity ! number of particles divided by volume=L^3
	double precision, optional :: kT
	double precision, optional :: U0
	type( Problem ), pointer :: this

	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: v
	type( Hamiltonian ), pointer :: h

	double precision U0value, numParticles
	
		this => BoseGas( sizes, 1, chemicalPotential, particleDensity, kT )
		h => this % hamilt
		s => h % subsystems(0) % obj

		v => InteractionPairwise( s )
		if (present(U0)) then
			U0value = U0
		else
			U0value = 1.d0
		end if
		call v % setAmplitudeConst( U0value ) ! later can be changed

		call h % addPerturbation( v )		

	end function InteractingBoseGas_init



	function BoseGas_init( sizes, numSubsystems, chemicalPotential, particleDensity, kT ) result(this)
	type( Geometry ), pointer :: sizes
	double precision, optional :: chemicalPotential
	double precision, optional :: particleDensity ! number of particles divided by volume=L^3
	double precision, optional :: kT
	integer numSubsystems
	type( Problem ), pointer :: this

	class(DispersionRelation), pointer :: dr
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: v
	type( Efactor ), pointer :: F	
	type( Hamiltonian ), pointer :: h

	double precision numParticles
	integer i
	
!		dr => DispersionRelationParabolic( sizes )
		dr => DispersionRelationParabolic()

		h => Hamiltonian() ! no perturbations: add later
		do i=0, numSubsystems-1
			s => SubsystemBose( dr )
			call h % addSubsystem( s )
		end do
		
		F => EfactorExact( 1d-8 ) ! later can be changed

		allocate(this)
		call this % init( sizes, h, F, kT, chemicalPotential )
		this % hamiltonianCreatedInside = .TRUE. ! to call destructor later

		if (present(particleDensity)) then
			numParticles = sizes % volume * particleDensity
			call this % setChemicalPotentialForNumberOfParticles( 0, numParticles )
			do i=1, numSubsystems-1
				this % chemicalPotential(i) = this % chemicalPotential(0)
			end do
		end if
		
	end function BoseGas_init


	
	function InteractingFermiGas_init( sizes, spinRepeat, chemicalPotential, particleDensity, kT, U0_in ) result(this)
	type( Geometry ), pointer :: sizes
	logical spinRepeat ! TRUE = both spin projections have always equal occupation; FALSE = can differ
	double precision, optional :: chemicalPotential
	double precision, optional :: particleDensity ! number of particles divided by volume=L^3, from 0 to 2
	double precision, optional :: kT
	double precision, optional :: U0_in
	type( Problem ), pointer :: this

	class(DispersionRelation), pointer :: dr
	type( Subsystem ), pointer :: s, s2
	type( Perturbation ), pointer :: v
	type( Efactor ), pointer :: F	
	type( Hamiltonian ), pointer :: h

	double precision U0, numElectrons
	
		if (present(U0_in)) then
			U0 = U0_in
		else
			U0 = 1.d0
		end if

		dr => DispersionRelationParabolic( sizes )
		s => SubsystemFermi( dr )
		if (spinRepeat) then
			v => InteractionPairwise( s, s ) ! spin-up x spin-down
		else	
			s2 => SubsystemFermi( dr )
			v => InteractionPairwise( s, s2 ) ! one x two
		end if
		call v % setAmplitudeConst( U0 ) ! later can be changed

		h => Hamiltonian()
		call h % addPerturbation( v )
				
		F => EfactorExact( 1d-8 ) ! later can be changed

		allocate(this)
		call this % init( sizes, h, F, kT, 1d3 )
		this % hamiltonianCreatedInside = .TRUE. ! to call destructor later

		if (present(particleDensity)) then
			if (particleDensity.GT.2.d0) stop 'expecting electron density from 0 to 2'	
			numElectrons = sizes % volume * particleDensity
			call this % setChemicalPotentialForNumberOfParticles( 0, numElectrons/2 ) ! two spins
		end if
		if (present(chemicalPotential)) then
			this % chemicalPotential(:) = chemicalPotential
		end if
		if (.NOT.spinRepeat) then
			this % chemicalPotential(1) = this % chemicalPotential(0)
		end if
	end function InteractingFermiGas_init
	
	
	function ElectronGasWithPhonons_init( sizes, electronDensity, kT, &
&		M0, Upair0 ) result(this)
!sizes, electronDensity=2.d0, kT=10.d0, M0=1.d0, Upair=1.d0 )
	type( Geometry ), pointer :: sizes
	double precision electronDensity ! per site, i.e. from 0 do 2
	double precision kT
	double precision, optional ::  M0, Upair0
	type( Problem ), pointer :: this

	class(DispersionRelation), pointer :: dr1, dr2
	type( Subsystem ), pointer :: el, phon
	type( Perturbation ), pointer :: V_eph, V_ee
	type( Efactor ), pointer :: F
	type( Hamiltonian ), pointer :: h

	double precision M0value, numElectrons
	
		if (electronDensity.GT.2.d0) stop 'expecting electron density from 0 to 2'

		if (present(M0)) then
			M0value = M0
		else
			M0value = 1.d0
		end if

		dr2 => DispersionRelationParabolic( sizes ) ! TODO!!!!! deallocate in the end
		el => SubsystemFermi( dr2 )
		dr1 => DispersionRelationLinear( sizes )
		call dr1 % setCoeff( 0.1d0 )
		phon => SubsystemBose( dr1 )

		V_eph => InteractionWithPhonons( el, phon )
		call V_eph % setAmplitudeConst( M0value ) ! later can be changed

		h => Hamiltonian()
		call h % addSubsystem( el )
		call h % addSubsystem( phon )
		call h % addPerturbation( V_eph )
		
		if (present(Upair0)) then ! if asked to add e-e interaction
			V_ee => InteractionPairwise( el, el )
			call V_ee % setAmplitudeConst( Upair0 ) ! later can be changed
			call h % addPerturbation( V_ee )
		end if
		
		F => EfactorGauss( 1d-3 ) ! later can be changed

		allocate(this)
		call this % init( sizes, h, F, kT, 1d3 )
		this % hamiltonianCreatedInside = .TRUE. ! to call destructor later
		numElectrons = sizes % volume * electronDensity
		call this % setChemicalPotentialForNumberOfParticles( 0, numElectrons/2 ) ! two spins
		call this % setChemicalPotential( 0.d0, 1 ) ! zero for phonons
		
	end function ElectronGasWithPhonons_init
	

end module Problem_module
