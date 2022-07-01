module QbeTermsCalculator_module

	use, intrinsic :: iso_c_binding 

	use fftw3
	use Output_module
	use Efactor_module
	use Geometry_module
	use Subsystem_module
	use Perturbation_module
	use QbeTerm_module
	use Occupations_module
	use Hamiltonian_module
	use MultidimensionalCycle_module

	implicit none
	
	integer, parameter :: VERSION_AS_IS_v1 = 1
	integer, parameter :: VERSION_AS_IS_v2 = 2
	integer, parameter :: VERSION_AS_IS_v3 = 3
	integer, parameter :: VERSION_AS_IS = VERSION_AS_IS_v3
	integer, parameter :: VERSION_FAST = 4

	type QbeTermsCalculator
		type(Geometry), pointer :: sizes
		type(Hamiltonian), pointer :: hamilt
		integer numS, numP, numT ! numSubsystems, numPerturbations, numThreads
		logical useKroneckers
	contains
		procedure :: addValue => QTC_Test_addValue
		procedure :: compareTo => QTC_compareTo ! qbe2, occ
		procedure :: show => QTC_Test_show
		
		final :: QTC_Test_destructor
		
	end type  QbeTermsCalculator

	type, extends(QbeTermsCalculator) :: QbeTermsCalculatorSmooth
		type(Efactor), pointer :: F
		double precision :: factorBoundary

		type(MultidimensionalCycle_p), allocatable :: mdc(:,:,:) ! thread, half, perturbation
		logical, allocatable :: term_ps(:,:) ! subsystem, perturbation
		integer, allocatable :: numHalves(:) ! 1 or 2
		
	contains
		procedure :: initFromGHFK => QTC_Smooth_initFromGHFK ! geometry, Hamiltonian, Efactor, useKroneckers
		final :: QTC_Smooth_destructor
	end type QbeTermsCalculatorSmooth
	
	type, extends(QbeTermsCalculatorSmooth) :: QbeTermsCalculatorSimple
		type(QbeTermAsIs_p), allocatable :: term(:,:)  ! 1:2, 0:numP-1
		integer qt_version ! v1/v2/v3 = AsIs / GeneralSum / SumOfABCD
		
	contains
		procedure :: addValue => QTC_Simple_addValue
		procedure :: show => QTC_Simple_show
		final :: QTC_Simple_destructor
	end type QbeTermsCalculatorSimple



		
	type, extends(QbeTermsCalculatorSmooth) :: QbeTermsCalculatorFFT
		type(QbeTermFFT_p), allocatable :: term(:,:)   ! 1:2, 0:numP-1
		double precision de
		double precision, allocatable :: emin(:) ! for each subsystem
		integer numE
		integer*8 volume4D, volume3D
		logical prepared, cachePrepared(0:1)
		
		type(MDC_CachedList_p), allocatable :: mdc_cache(:,:,:) ! p x half x rtaIndex
		integer maxk1multiplier
				
	! FFT plans
		! for E
		type(C_PTR) :: plan1d
		type(C_PTR) :: plan1dBack

		! for X Y Z E
		type(C_PTR) :: plan4dR2C 
		type(C_PTR) :: plan4dBackC2R

		! more for U(q) : X Y Z
		type(C_PTR) :: plan3dMany
		type(C_PTR) :: plan3dManyBack

	! work arrays
		complex*16, allocatable :: lineshape_e(:)
		complex*16, allocatable :: lineshape_g(:)

		double precision, allocatable :: n_keR(:,:,:,:)	! 4D : e, k
		double precision, allocatable :: S_keR(:,:,:,:)	! 4D : e, k  ! drop later / TODO: update if energies change
		double precision, allocatable :: calc_keR(:,:,:,:)	! 4D : e, k

		complex*16, allocatable :: S_rgR(:,:,:,:,:)	! 5D : g, r, subsystem
		complex*16, allocatable :: n_rgR(:,:,:,:,:)	! 5D : g, r, subsystem
		complex*16, allocatable :: calc_rgR(:,:,:,:)	! 4D : e, k

		complex*16, allocatable :: calc_k_many(:,:,:,:)	! 3D : k x numE
		complex*16, allocatable :: calc_r_many(:,:,:,:)	! 3D : r x numE

		logical, allocatable :: minusR(:,:,:,:) ! operator, half, ss, pp

	contains
		procedure :: init2 => QTC_FFT_init2
		procedure :: addValue => QTC_FFT_addValue
		procedure :: show => QTC_FFT_show

		procedure :: latePreparation => QTC_FFT_latePreparation
		procedure :: prepareWorkArrays => QTC_FFT_prepareWorkArrays
		procedure :: prepareMoreArrays => QTC_FFT_prepareMoreArrays
		procedure :: prepareFftPlans => QTC_FFT_prepareFftPlans
		procedure :: prepareLineshape => QTC_FFT_prepareLineshape
		procedure :: make_S_rg => QTC_FFT_make_S_rg
		procedure :: make_N_rg => QTC_FFT_make_N_rg
		
		procedure :: calculateTermForRightPart => QTC_FFT_calculateTermForRightPart
		procedure :: prepareCacheForKroneckers => QTC_FFT_prepareCacheForKroneckers
		procedure :: simpleCalculationForKroneckers => QTC_FFT_simpleCalculationForKroneckers
		
		procedure :: make_calc_rg => QTC_FFT_make_calc_rg
		procedure :: general_interaction_XY => QTC_FFT_general_interaction_XY
		procedure :: constant_interaction_XYZ => QTC_FFT_constant_interaction_XYZ

		procedure :: make_calc_ke => QTC_FFT_make_calc_ke
		procedure :: make_BC_for_Q => QTC_FFT_make_BC_for_Q
		procedure :: make_calc_k_3D_many => QTC_FFT_make_calc_k_3D_many
		procedure :: make_calc_r_3D_many => QTC_FFT_make_calc_r_3D_many
		procedure :: updateDN => QTC_FFT_updateDN

		final :: QTC_FFT_destructor
		
	end type QbeTermsCalculatorFFT



	interface QbeTermsCalculatorTest
		module procedure  QTC_Test_init
	end interface  QbeTermsCalculatorTest
	interface QbeTermsCalculatorSimple
		module procedure  QTC_Simple_init
	end interface  QbeTermsCalculatorSimple
	interface QbeTermsCalculatorFFT
		module procedure  QTC_FFT_init
	end interface  QbeTermsCalculatorFFT
	
! polymorphic subroutine
	abstract interface
		subroutine showT4P(this, p, s_i, version, fd)
			import QbeTermsCalculatorSmooth
			import Perturbation
			class(QbeTermsCalculatorSmooth) :: this
			type(Perturbation), pointer :: p
		    integer s_i, version, fd
		end subroutine showT4P
	end interface
	
	
contains

!----------------------------------------------------------------

	function QTC_Test_init( sizes, hamilt ) result(this)
	type(Geometry), pointer :: sizes
	type(Hamiltonian), pointer :: hamilt

	type( QbeTermsCalculator ), pointer :: this

		allocate( this )
		
		this % sizes => sizes
		this % hamilt => hamilt

	end function QTC_Test_init



	subroutine QTC_Test_destructor( this )
	type( QbeTermsCalculator ) :: this
	
!		print*, 'QTC_Test destructor'
		
	end subroutine QTC_Test_destructor

! -- main function

	subroutine QTC_Test_addValue( this, t, n, dn, rta )
	class( QbeTermsCalculator ), intent(INOUT) :: this
	type( Occupations ), pointer :: n, dn
	double precision t
	logical rta ! FALSE = get dn/dt, TRUE = get dn/dt for relaxation time

		if (rta) then
			dn % values(:,:,:,:) = n % values(:,:,:,:)	! solution=exp(t)
!			dn % values(:,:,:,:) = -sin(t)				! solution=cos(t)
		else
			dn % values(:,:,:,:) = - 1.d0	! relaxation time = -1/dn = 1.d0  
		end if
	
	end subroutine QTC_Test_addValue

! -- for debug

!--

	subroutine QTC_Test_show( this, o )
	class( QbeTermsCalculator ) :: this
	class( Output ), pointer :: o
	
		call o % write( 'QbeTermsCalculator, Test' )
	
	end subroutine QTC_Test_show

!----------------------------------------------------------------

	double precision function QTC_compareTo( this, qtc2, t, n, rta )
	class( QbeTermsCalculator ) :: this
	class( QbeTermsCalculator ), pointer :: qtc2
	double precision t
	type( Occupations ), pointer :: n
	logical rta
	
	class( Occupations ), pointer :: dn1, dn2
	
	double precision v
	integer i, j
		
		dn1 => Occupations( n )
		dn2 => Occupations( n )

		dn1 % values = 0.d0
		dn2 % values = 0.d0
		
		call this % addValue( t, n, dn1, rta )
		call qtc2 % addValue( t, n, dn2, rta )
		v = dn1 % compareTo( dn2 )
		if (v.GT.1d-8) then
			print*, ' dn1: sum=', sum(dn1%values(:,:,:,:)) 
			print*, ' dn2: sum=', sum(dn2%values(:,:,:,:)) 
			do i=0, n%N4-1
				write(*, '(A,I0,A)', advance='no') ' (:00',i,')='
				do j=0, SIZE(dn1%values(:,0,0,i))-1
					write(*, '(F24.13)', advance='no') dn1%values(j,0,0,i)
				end do
				print*
				write(*, '(A,I0,A)', advance='no') ' (:00',i,')='
				do j=0, SIZE(dn2%values(:,0,0,i))-1
					write(*, '(F24.13)', advance='no') dn2%values(j,0,0,i)
				end do
				print*
			end do
		end if

		deallocate( dn2 )
		deallocate( dn1 )

		QTC_compareTo = v
	end function QTC_compareTo

!----------------------------------------------------------------

	subroutine QTC_Smooth_initFromGHFK( this, sizes, hamilt, F, useKroneckers )

	use omp_lib
	class( QbeTermsCalculatorSmooth ), intent(INOUT) :: this	
	type(Geometry), pointer :: sizes
	type(Hamiltonian), pointer :: hamilt
	type(Efactor), pointer :: F
	logical useKroneckers

	
	type(Perturbation), pointer :: p		
	type(Subsystem), pointer :: s
	integer pp, ss, pos, t, numHalves, half

		this % sizes => sizes
		this % hamilt => hamilt
		
		this % F => F
		this % factorBoundary = defaultFactorBoundary
				
		this % numS = hamilt % numSubsystems
		this % numP = this % hamilt % numPerturbations
		this % numT = omp_get_max_threads()

		this % useKroneckers = useKroneckers
		
		allocate( this % mdc(0:this % numT-1, 1:2, 0:this % numP-1) )
		allocate( this % term_ps(0:this % numP-1, 0:this % numS-1) )
		allocate( this % numHalves(0:this % numP-1) )

		do pp=0, this % numP - 1
			p => hamilt % perturbations(pp) % obj
			call p % prepare

			numHalves = 2
			if ( p%hermitian ) numHalves = 1
			this % numHalves(pp) = numHalves
			
			do half=1, this%numHalves(pp)
				if (half.GT.1) p => p % hc

				do t=0, this % numT - 1
					this % mdc(t,half,pp) % obj => MultidimensionalCycle( sizes, p, F )
				end do

			end do

			do ss=0, this % numS-1 
				s => this % hamilt % subsystems( ss ) % obj
				pos = p % indexOf( s )
				this % term_ps(pp,ss) = (pos.GE.0)
			end do
		end do
		
	end subroutine QTC_Smooth_initFromGHFK

	subroutine QTC_Smooth_destructor( this )
	type( QbeTermsCalculatorSmooth ) this	
	
	integer t, pp, half

!		print*, 'QTC_Smooth_destructor'
		
		deallocate( this % term_ps )
		
		do pp=this % numP - 1, 0, -1
			do half=1, this % numHalves(pp)
				do t=this % numT - 1, 0, -1
					deallocate( this % mdc( t, half, pp ) % obj )
				end do
			end do
		end do
		deallocate( this % mdc )
		deallocate( this % numHalves )

	end subroutine QTC_Smooth_destructor
	
!----------------------------------------------------------------

	include 'QbeTermsCalculatorSimple.f90'
	include 'QbeTermsCalculatorFFT.f90'

end module QbeTermsCalculator_module
