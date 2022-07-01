module fftw3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module

include 'Output.f90'
include 'TableOfValues.f90'
include 'Lifetime.f90'
include 'Geometry.f90'
include 'DispersionRelation.f90'
include 'Subsystem.f90'
include 'Perturbation.f90'
include 'Hamiltonian.f90'
include 'Occupations.f90'
include 'QbeTerm.f90'
include 'Efactor.f90'
include 'MultidimensionalCycle.f90'
include 'QbeTermsCalculator.f90'
include 'Problem.f90'
include 'BoltzmannEquation.f90'
include 'FiniteDifferenceSolver.f90'


module flbe

use Lifetime_module
use DispersionRelation_module
use Subsystem_module
use Perturbation_module
use Hamiltonian_module
use QbeTerm_module
use Occupations_module
use Geometry_module
use MultidimensionalCycle_module
use BoltzmannEquation_module
use Problem_module
use FiniteDifferenceSolver_module



contains


	subroutine flbe_init

	use fftw3
	use omp_lib
	integer void

		void = fftw_init_threads()
		if (void.EQ.0) stop 'flbe init: fftw3 init fail'
		call fftw_plan_with_nthreads( omp_get_max_threads() )
!		call fftw_set_timelimit(20.d0) ! stop estimations at 20 seconds
		
	end subroutine flbe_init



	subroutine flbe_done
		call fftw_cleanup_threads()
		subsystem_count = 0 ! reset subsystem IDs
	end subroutine flbe_done


end module flbe
