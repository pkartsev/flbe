
!	call benchmark( 16,16,16, 100, 8 )
!	call benchmark_byThreads( 4, 1000, 64 )
!	call benchmark_byThreads( 8, 500, 16 )
!	call benchmark_byThreads( 16, 10, 16 )
!	call benchmark_byThreads( 24, 5, 16 )
!	call benchmark_byThreads( 32, 1, 16 )
!	call benchmark_byThreads( 48, 3, 5, 16 )
!	call benchmark_byThreads( 16, 2, 50, 16 )
!	call benchmark_byThreads( 4, 1000, 32 )
	call benchmark_bySize( 8, 10, 2, 8, 32, 4 )
!	call benchmark_bySize( 8, 5, 3, 40, 52, 4 )
!	call benchmark_bySize( 8, 5, 3, 56, 64, 4 )
!	call benchmark_bySize( 8, 5, 2, 56, 96, 8 )
!	call benchmark_bySize( 8, 1, 2, 104, 128, 8 )
!	call benchmark_bySize( 8, 1, 2, 8, 128, 8 )
!	call benchmark_bySize( 8, 1, 4, 60, 4 )
!	call benchmark_bySize( 16, 10, 4, 32, 4 )
!	call benchmark_bySize( 8, 1, 32, 48, 4 )

end program

!-----------------------------------

include 'flbe.f90'

!-----------------------------------

subroutine benchmark_byThreads( L, d, numSteps, numThreadsMax )
use flbe
implicit none
integer L, d, numSteps, numThreadsMax

integer nt

	nt = 1
	do while (nt.LE.numThreadsMax)
		call benchmark( L, d, numSteps, nt )
		nt = nt*2
	end do

end subroutine benchmark_byThreads




subroutine benchmark_bySize( numThreads, numSteps, d, Lmin, Lmax, dL )
use flbe
implicit none
integer numThreads, numSteps, Lmin, Lmax, dL, d

integer L

	

	do L = Lmin, Lmax, dL
		call benchmark( L, d, numSteps, numThreads )
	end do

end subroutine benchmark_bySize


!-----------------------------------


subroutine benchmark(L, d, numSteps, numThreads)

use flbe
implicit none
integer LX, LY, LZ, L, d, numSteps, numThreads

type(Geometry), pointer :: g
type(Problem), pointer :: prob
class(BoltzmannEquation), pointer :: qbe
class(FD_RightPart), pointer :: rp
class(FiniteDifferenceSolver), pointer :: fd
type(Efactor), pointer :: F

integer i, numShow
double precision dt, t, dv, de

integer iTimes1, iTimes2, rate
real cpuTimeDelta

	LX = L
	LY = L
	LZ = L
	if (d.EQ.2) LZ=1
	if (d.EQ.1) LY=1
	

	write(*, '(A)', advance='no') 'Benchmark: '
	write(*, '(A,I0,A,I0,A,I0)', advance='no') 'Size ',LX,'x',LY,'x',LZ
	write(*, '(A,I0,A,I0)') ' num.threads: ', numThreads, ' num.steps: ', numSteps

	CALL system_clock(count_rate=rate)

! openmp number of threads
	call omp_set_num_threads( numThreads )
	call flbe_init

! sizes
	g => Geometry(LX, LY, LZ )

! equation
	prob => InteractingBoseGas( g, 3.d0, 1.d0 )
	prob % useKroneckers = .FALSE. 
	de = 1.d0
!	qbe => BoltzmannEquation( prob, de, VERSION_FAST )
	qbe => BoltzmannEquation( prob, de, VERSION_AS_IS_v1 )
	
! solver
	rp => qbe
	fd => FD_Euler( rp )
	
	fd % n % values(:,:,:,:) = 0.1d0
	fd % n % values(1,0,0,:) = 0.4d0
	if (subsystem_count.GT.1) then
		fd % n % values(:,:,:,1) = 0.3d0
		fd % n % values(1,0,0,1) = 0.713d0
	end if
	if (subsystem_count.GT.2) then
		fd % n % values(:,:,:,2) = 0.56d0
		fd % n % values(1,0,0,2) = 0.1245d0
	end if

! simulation and measuring time
	dt = 1d-5
	t = 0.d0

!	call qbe % latePreparation

    numShow = numSteps/10
    if (numShow.EQ.0) numShow=1

	! warm-up for FFTW3 late initialization
	write(*, '(A)', advance='no') 'warm-up...'
	call fd%step(dt, t)
	print*, 'done'
	t = 0.d0
	
    call SYSTEM_CLOCK(iTimes1)
	do i=0, numSteps-1
		call fd%step(dt, t)
		if (MOD(i+1,numShow).EQ.0) then
			write(*, '(A,I0,A,I0)', advance='no') CHAR(13), i+1, '/', numSteps
		end if
	end do

    call SYSTEM_CLOCK(iTimes2)
    cpuTimeDelta = real(iTimes2-iTimes1)/real(rate)
	write(*, '(A,F8.2,A)') ' DONE in ', cpuTimeDelta, ' seconds.'
	print*	
	print*

	deallocate(fd)
	deallocate(qbe)
	deallocate(prob)
	deallocate(g)

	call flbe_done
	open(1, file='benchmark_result.txt', access='append')
	write(1, *) LX, LY, LZ, numThreads, cpuTimeDelta/numSteps
	close(1)

end subroutine benchmark

