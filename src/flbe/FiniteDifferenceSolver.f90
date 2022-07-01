module FiniteDifferenceSolver_module

	use Geometry_module
	use Occupations_module
	use BoltzmannEquation_module

	implicit none
	
	type, abstract :: FiniteDifferenceSolver
		class( FD_RightPart ), pointer :: rp
		class(Occupations), pointer :: n
		type(Occupations), pointer :: k1 ! k1 for Euler
	contains
		procedure :: constructor => FD_constructor
		procedure(FD_step), deferred :: step
		procedure :: simulate => FD_simulate
		procedure :: compareTo => FD_compareTo
		
	end type FiniteDifferenceSolver

! - children classes

	type, extends (FiniteDifferenceSolver) :: FD_Euler
	contains
		procedure :: step => Euler_step
		final :: FD_Euler_destructor
	end type FD_Euler


	type, extends (FD_Euler) :: FD_RungeKutta2
		type(Occupations), pointer :: tmp ! k1, tmp for RK2
	contains
		procedure :: step => RK2_step
		final :: FD_RK2_destructor
	end type FD_RungeKutta2

	type, extends (FD_RungeKutta2) :: FD_RungeKutta4
		type(Occupations), pointer :: k2 ! k1..k4, tmp for RK4
		type(Occupations), pointer :: k3
		type(Occupations), pointer :: k4
	contains
		procedure :: step => RK4_step
		final :: FD_RK4_destructor
	end type FD_RungeKutta4

	
	type, extends (FD_Euler) :: FD_DormandPrince
		type(Occupations), pointer :: k2
		type(Occupations), pointer :: k3
		type(Occupations), pointer :: k4
		type(Occupations), pointer :: k5
		type(Occupations), pointer :: k6
		type(Occupations), pointer :: tmp
	contains
		procedure :: step => DP_step
		final :: FD_DP_destructor
	end type FD_DormandPrince


! -- constructor names

	interface FD_Euler
		module procedure FD_Euler_Init
	end interface FD_Euler

	interface FD_RungeKutta2
		module procedure FD_RK2_Init
	end interface FD_RungeKutta2

	interface FD_RungeKutta4
		module procedure FD_RK4_Init
	end interface FD_RungeKutta4

	interface FD_DormandPrince
		module procedure FD_DP_Init
	end interface FD_DormandPrince
	
! polymorphic subroutine
	abstract interface
		subroutine FD_step(this, dt, t)
			import FiniteDifferenceSolver
			class(FiniteDifferenceSolver) :: this
			double precision dt
			double precision, intent(INOUT) :: t
		end subroutine
	end interface
	
contains

	subroutine FD_constructor( this, rp )
	class(FiniteDifferenceSolver), intent(INOUT) :: this
	class( FD_RightPart ), target :: rp

		this % rp => rp
		this % n => Occupations( rp % sizes % L(:), rp % hamilt % numSubsystems )
		this % k1 => Occupations( this % n )
	
	end subroutine FD_constructor
	
! --
	
	function FD_Euler_Init( rp ) result(this)
	class( FD_Euler ), pointer :: this
	class( FD_RightPart ), target :: rp
	
		allocate( this )
		call this % constructor( rp )! creates n, k1

	end function FD_Euler_Init
	
! --	
	function FD_RK2_Init( rp ) result(this)
	type( FD_RungeKutta2 ), pointer :: this
	class( FD_RightPart ), target :: rp

		allocate( this )

		call this % constructor( rp ) ! creates n, k1
		this % tmp => Occupations( this % n )

	end function FD_RK2_Init

! --	
	function FD_RK4_Init( rp ) result(this)
	type( FD_RungeKutta4 ), pointer :: this
	class( FD_RightPart ), target :: rp

		allocate( this )

		call this % constructor( rp ) ! creates n, k1
		this % k2  => Occupations( this % n )
		this % k3  => Occupations( this % n )
		this % k4  => Occupations( this % n )
		this % tmp => Occupations( this % n )

	end function FD_RK4_Init

! --	
	function FD_DP_Init( rp ) result(this)
	type( FD_DormandPrince ), pointer :: this
	class( FD_RightPart ), target :: rp
	

		allocate( this )

		call this % constructor( rp ) ! creates n, k1
		this % k2  => Occupations( this % n )
		this % k3  => Occupations( this % n )
		this % k4  => Occupations( this % n )
		this % k5  => Occupations( this % n )
		this % k6  => Occupations( this % n )
		this % tmp => Occupations( this % n )

	end function FD_DP_Init
	
! -- destructors

	subroutine FD_Euler_destructor( this )
	type( FD_Euler ) :: this

!		print*, 'FD_Euler: destructor'
		deallocate( this % k1 )
		deallocate( this % n )

	end subroutine FD_Euler_destructor
 	
!--

	subroutine FD_RK2_destructor( this )
	type( FD_RungeKutta2 ) :: this

!		print*, 'FD_RungeKutta2: destructor'
		deallocate( this % tmp )

	end subroutine FD_RK2_destructor
 	
!--

	subroutine FD_RK4_destructor( this )
	type( FD_RungeKutta4 ) :: this

!		print*, 'FD_RungeKutta4: destructor'
		deallocate( this % k4 )
		deallocate( this % k3 )
		deallocate( this % k2 )

	end subroutine FD_RK4_destructor
 	
!--

	subroutine FD_DP_destructor( this )
	type( FD_DormandPrince ) :: this

!		print*, 'FD_DormandPrince: destructor'
		deallocate( this % tmp )
		deallocate( this % k6 )
		deallocate( this % k5 )
		deallocate( this % k4 )
		deallocate( this % k3 )
		deallocate( this % k2 )

	end subroutine FD_DP_destructor
 	
!--
	
	subroutine FD_simulate(this, tmin, tmax, numSteps, numShow_in)
	class(FiniteDifferenceSolver) :: this
	double precision tmin, tmax
	integer numSteps
	integer, optional :: numShow_in

	double precision t, dt
	integer i, numShow
	integer iTimes1, iTimes2, rate
	real cpuTimeDelta

	numShow = 1 + (numSteps-1)/10
	if (present(numShow_in)) numShow = numShow_in

	CALL system_clock(count_rate=rate)
    call SYSTEM_CLOCK(iTimes1)

	dt = (tmax-tmin)/numSteps
	t = tmin
	do i=0, numSteps-1
		call this%step(dt, t)
		if (MOD(i+1,numShow).EQ.0) then
			write(*, '(A,I0,A,I0)', advance='no') CHAR(13), i+1, '/', numSteps
		end if
	end do
    call SYSTEM_CLOCK(iTimes2)
    cpuTimeDelta = real(iTimes2-iTimes1)/real(rate)
	write(*, '(A,F8.2,A)') ' DONE in ', cpuTimeDelta, ' seconds.'

	end subroutine FD_simulate
	

!--------

	subroutine Euler_step( this, dt, t )
	class(FD_Euler) :: this
	double precision, intent(INOUT) :: t
	double precision dt
	
! k1 = f ( t, n )
		call this % rp % getRightPart( t, this % n, this % k1 )
! result = n  +  dt * k1
		call this % n % add1( this % k1, this % n, dt )
		
		t = t + dt
	
	end subroutine Euler_step

	subroutine RK2_step( this, dt, t )
	class(FD_RungeKutta2) :: this
	double precision, intent(INOUT) :: t
	double precision dt
	
! k1 = f ( t, n )
		call this % rp % getRightPart( t, this % n, this % k1 )
		call this % n % add1( this % k1, this % tmp, dt/2 )
				
! k2 = f( t+dt/2,  n  +  dt/2 * k1 )
		call this % rp % getRightPart( t+dt/2, this % tmp, this % k1 )

! result =  n  +  dt * k2
		call this % n % add1( this % k1, this % n, dt ) ! new = old + dt*k2
				
		t = t + dt
	
	end subroutine RK2_step

	subroutine RK4_step( this, dt, t )
	class(FD_RungeKutta4) :: this
	double precision, intent(INOUT) :: t
	double precision dt
	
! k1 = f ( t, n )
		call this % rp % getRightPart( t, this % n, this % k1 )
		call this % n % add1( this % k1, this % tmp, dt/2 )
				
! k2 = f ( t+dt/2,  n + dt/2 * k1 )
		call this % rp % getRightPart( t+dt/2, this % tmp, this % k2 )
		call this % n % add1( this % k2, this % tmp, dt/2 )
				
! k3 = f ( t+dt/2,  n + dt/2 * k2 )
		call this % rp % getRightPart( t+dt/2, this % tmp, this % k3 )
		call this % n % add1( this % k3, this % tmp, dt )
				
! k4 = f ( t+dt,  n + dt * k3 )
		call this % rp % getRightPart( t+dt, this % tmp, this % k4 )

! result =  n + dt/6 * ( k1 + 2*k2 + 2*k3 + k4 )
		call this % n % add4( this%k1, this%k2, this%k3, this % k4, this % n, &
&			1.d0, 2.d0, 2.d0, 1.d0, dt/6 )
		
		t = t + dt
	
	end subroutine RK4_step

	subroutine DP_step( this, dt, t )
	class(FD_DormandPrince) :: this
	double precision, intent(INOUT) :: t
	double precision dt
	
! k1 = f ( t, n )
		call this % rp % getRightPart( t, this % n, this % k1 )

! k2 = f ( t+dt/5,  n + dt/5 * k1 )
		call this % n % add1( this % k1,   this % tmp, dt/5 )
		call this % rp % getRightPart( t+dt/5, this % tmp, this % k2 )
				
! k3 = f ( t+(3/10)dt, n + dt * ( 3/40 * k1 + 9/40 * k2) )
		call this % n % add2( this % k1, this % k2, &
&			this % tmp, &
&			3.d0/40, 9.d0/40, dt )
		call this % rp % getRightPart( t+0.3d0*dt, this % tmp, this % k3 )				

! k4 = f ( t + (4/5)dt, n + dt * ( 44/45 * k1 - 56/15 * k2 + 32/9 * k3 ) )
		call this % n % add3( this % k1, this % k2, this % k3,  &
&			this % tmp, &
&			44.d0/45, -56.d0/15, 32.d0/9,  dt )
		call this % rp % getRightPart( t+0.8d0*dt, this % tmp, this % k4 )

! k5 = f ( t+(8/9)dt, n + dt * ( 19372/6561 * k1 - 25360/2187 * k2 + 64448/6561 * k3  - 212/729 * k4 ) )
		call this % n % add4( this % k1, this % k2, this % k3, this % k4, &
&			this % tmp, &
&  			19372.d0/6561,  -25360.d0/2187, 64448.d0/6561, -212.d0/729,    dt )
		call this % rp % getRightPart( t+(8.d0/9)*dt, this % tmp, this % k5 )

! k6 = f ( t+dt, n + dt * ( 9017/3168 * k1 - 355/33 * k2 + 46732/5247 * k3  + 49/176 * k4  - 5103/18656 ( k5 ) )
		call this % n % add5( this % k1, this % k2, this % k3, &
&			this % k4, this % k5,  &
&			this % tmp, &
& 			9017.d0/3168,  -355.d0/33,  46732.d0/5247, 49.d0/176,  -5103.d0/18656,     dt )
		call this % rp % getRightPart( t+dt, this % tmp, this % k6 )

! result = n + ((35/384)*k1 + (500/1113)*k3 + (125/192)*k4 - (2187/6784)*k5 + (11/84)*k6 ) * dt
		call this % n % add5( this % k1, this % k3, this % k4, &
&			this % k5, this % k6, &
&			this % n, &
&			35.d0/384, 500.d0/1113, 125.d0/192, -2187.d0/6784, 11.d0/84,    dt)
		
		t = t + dt
	
	end subroutine DP_step
	
!----------------------------------------

	double precision function FD_compareTo( this, fd2, dt, t, numSteps_in )
	class(FiniteDifferenceSolver) :: this
	class(FiniteDifferenceSolver), pointer :: fd2
	double precision t
	double precision dt
	integer, optional :: numSteps_in
	
	integer :: numSteps, i
	double precision t1, t2
		
		numSteps = 1
		if (present(numSteps_in)) numSteps = numSteps_in

		t1 = t
		t2 = t

		do i=0, numSteps-1
			call this % step( dt, t )
			call fd2  % step( dt, t )
		end do
		FD_compareTo = this % n % compareTo( fd2 % n )
		
	end function FD_compareTo
	
end module FiniteDifferenceSolver_module
