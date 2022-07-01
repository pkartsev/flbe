include 'flbe.f90'

use flbe
implicit none
integer, parameter :: L=4

class(Problem), pointer :: prob
class(Geometry), pointer :: sizes
class(BoltzmannEquation), pointer :: equation
class(FiniteDifferenceSolver), pointer :: solver
class(Output), pointer :: o

	call omp_set_num_threads( 4 )
	call flbe_init

	sizes => Geometry( L, L, L )
	prob => InteractingBoseGas( sizes, particleDensity=2.d0, kT=10.d0 )

	equation => BoltzmannEquation( prob )
	solver => FD_Euler( equation )	
	solver % n % values(:,:,:,:) = prob % equilibriumOccupations % values(:,:,:,:)
	solver % n % values( 1, 1, 1,:) = prob % equilibriumOccupations % values( 1, 1, 1,:) + 0.1d0

	print*, ' n(111)=', solver % n % values( 1, 1, 1,:)
	call solver % simulate( tmin=0.d0, tmax=1.0d-5, numSteps=1000 )
	print*, ' n(111)=', solver % n % values( 1, 1, 1,:)
	
	call flbe_done

end

