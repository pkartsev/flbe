include 'flbe.f90'

include 'v0simple.f90'

module tests

	integer, parameter :: debugFileNumber = 10
	logical, parameter :: verbose = .FALSE.
!	logical, parameter :: verbose = .TRUE.

	integer, parameter :: numThreads = 4
	logical, parameter :: testParallel = (numThreads.GT.1)
	
	! levels to test: from 1 to 4
	integer, parameter :: minLevel = 1
	integer, parameter :: maxLevel = 4
	
contains

	include 'test_printTwice.f90'

	include 'testGeometry.f90'
	include 'testDispersionRelation.f90'
	include 'testOccupations.f90'
	include 'testQbeTerms.f90'
	include 'testFiniteDifferenceSolver.f90'
	include 'testQbeTermsCalculator.f90'
	include 'testRTA.f90'

end module

!-----------------------------------------------

use tests

integer level

	if (verbose) then
		write(*,'(A)') 'Testing FLBE functionality...'
	else
		write(*,'(A)') 'Testing FLBE functionality... (set verbose=.TRUE. to leave files with debug info)'
	end if
	do level=minLevel, maxLevel
	
		select case(level)
			case(1)
				print*, '--- Level 1: basic support classes ---'
				call testGeometry
				call testDispersionRelation
				call testOccupations
			case(2)
				print*, '--- Level 2: finite-difference solver ---'
				call testFiniteDifferenceSolver
			case(3)
				print*, '--- Level 3: term calculations ---'
				call testQbeTerms
				call testQbeTermsCalculator
			case(4)
				print*, '--- Level 4: RTA vs dn/dt at thermal equilibrium ---'
				call testRTA
		end select
	end do
	write(*,'(A)') 'Testing FLBE functionality... ALL OK                                     '
	
end program
