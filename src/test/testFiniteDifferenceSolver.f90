subroutine testFiniteDifferenceSolver
! FD version Euler/RK2/RK4/DormandPrince
! problem version: const/exp/sin

use omp_lib
implicit none

integer problemVersion, fdVersion, serial, numThreads, numThreads_saved, L, numSteps
double precision dt
	
	numThreads_saved = omp_get_max_threads()
	numThreads = 16

	numSteps = 1000
	dt = 1d-4
	call printTwice_openFile('testFiniteDifferenceSolver.log', ' [ test FiniteDifferenceSolver...')
	
	call printTwiceSIS( '  num.steps= ', numSteps, '')
	
	call omp_set_num_threads( numThreads )
	if (testParallel) then
		call printTwiceSIS( '  = OpenMP version (num. threads: ', numThreads, ') = ')
	end if
	do problemVersion=1,3
		call testFDwithParams( problemVersion, numSteps, dt )
	end do
	call printTwice_closeFile( ' OK ]' )
	
	call omp_set_num_threads( numThreads_saved )

end subroutine testFiniteDifferenceSolver




subroutine testFDwithParams( problemVersion, numSteps, dt )
use flbe
implicit none

integer problemVersion, numSteps
logical verbose
double precision dt

class(Geometry), pointer :: sizes
class(Problem), pointer :: prob

class( FD_RightPart ), pointer :: be
class(FiniteDifferenceSolver), pointer :: euler, rk2, rk4, dp
integer i, L
double precision t, v1,v2,v4,v6, n0, density, kT
double precision diff01, diff02, diff04, diff06
double precision v0 ! analytical result
logical good
real r

	L = 4
	density = 0.5d0
	kT = 1.d0
	
	sizes => Geometry( L, L, L )
	prob => BoseGas( sizes, 1, density, kT )

	be => FD_RightPart_Test( problemVersion, prob )

	euler => FD_Euler( be )
	rk2 => FD_RungeKutta2( be )
	rk4 => FD_RungeKutta4( be )
	dp => FD_DormandPrince( be )

	call RANDOM_NUMBER(r)

	select case( problemVersion )
		case(1)	! n' = 0,  n(0) = R
			n0 = r
			v0 = r	! solution n(t) = R
			call printTwice( '  == test problem #1: n(t)=const')
		case(2)		! n'=n,  n(0) = R
			n0 = r
			v0 = r*exp( dt*numSteps )	! solution n(t) = R*exp(t)
			call printTwice( '  == test problem #2: n(t)=C*exp(t)')
		case(3)		! n' = cos(t), n(0) = R
			n0 = r
			v0 = sin( dt*numSteps ) + r	! solution n(t) = sin(t) + R
			call printTwice( '  == test problem #3: n(t)=sin(t) + C')
	end select

	t=0.d0
	euler % n % values = n0
	do i=1, numSteps
		call euler % step( dt, t )
	end do
	v1 = euler % n % values(0,0,0,0)
	
	t=0.d0
	rk2 % n % values = n0
	do i=1, numSteps
		call rk2 % step( dt, t )
	end do
	v2 = rk2 % n % values(0,0,0,0)
	
	t=0.d0
	rk4 % n % values = n0
	do i=1, numSteps
		call rk4 % step( dt, t )
	end do
	v4 = rk4 % n % values(0,0,0,0)
	
	t=0.d0
	dp % n % values = n0
	do i=1, numSteps
		call dp % step( dt, t )
	end do
	v6 = dp % n % values(0,0,0,0)
	
	diff01 = abs( v0-v1 )
	diff02 = abs( v0-v2 )
	diff04 = abs( v0-v4 )
	diff06 = abs( v0-v6 )
	
	good = ((abs(diff01).LT.1d-5) .AND. (abs(diff02).LT.1d-9) .AND. (abs(diff04).LT.1d-14) .AND. (abs(diff06).LT.1d-14))
	if (good) then
		call printTwice( '       OK  ')
	else
		call printTwice( '       FAIL')
		call printTwiceSFSFS( ' difference from analytics: Euler=', diff01, '            RK2=', diff02, ' ')
		call printTwiceSFSFS( ' difference from analytics:   RK4=', diff04, ' Dormand-Prince=', diff06, ' ')
		stop 'see details in testFiniteDifferenceSolver.log'
	end if
	
	deallocate( dp )
	deallocate( rk4 )
	deallocate( rk2 )
	deallocate( euler )
	deallocate( be )
	deallocate( prob )
	deallocate( sizes )
	
	
end subroutine testFDwithParams
