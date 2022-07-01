
!	call benchmark( 16,16,16, 100, 8 )
!	call benchmark_v0_byThreads( 4, 1000, 64 )
!	call benchmark_v0_byThreads( 8, 500, 16 )
!	call benchmark_v0_byThreads( 16, 10, 16 )
!	call benchmark_v0_byThreads( 24, 5, 16 )
!	call benchmark_v0_byThreads( 32, 1, 16 )
!	call benchmark_v0_byThreads( 48, 3, 5, 16 )
!	call benchmark_v0_byThreads( 16, 2, 50, 16 )
!	call benchmark_v0_byThreads( 4, 1000, 32 )
!	call benchmark_v0_bySize( 8, 10, 2, 16, 32, 4 )
	call benchmark_v0_bySize( 8, 5, 3, 4, 12, 4 )
!	call benchmark_v0_bySize( 8, 5, 3, 56, 64, 4 )
!	call benchmark_v0_bySize( 8, 5, 2, 56, 96, 8 )
!	call benchmark_v0_bySize( 8, 1, 2, 104, 128, 8 )
!	call benchmark_v0_bySize( 8, 1, 2, 8, 128, 8 )
!	call benchmark_v0_bySize( 8, 1, 4, 60, 4 )
!	call benchmark_v0_bySize( 16, 10, 4, 32, 4 )
!	call benchmark_v0_bySize( 8, 1, 32, 48, 4 )

end program

!-----------------------------------

include 'flbe.f90'

include 'v0simple.f90'

!-----------------------------------

subroutine benchmark_v0_byThreads( L, d, numSteps, numThreadsMax )
use flbe
implicit none
integer L, d, numSteps, numThreadsMax

integer nt

	nt = 1
	do while (nt.LE.numThreadsMax)
		call benchmark( L, d, numSteps, nt )
		nt = nt*2
	end do

end subroutine benchmark_v0_byThreads




subroutine benchmark_v0_bySize( numThreads, numSteps, d, Lmin, Lmax, dL )
use flbe
implicit none
integer numThreads, numSteps, Lmin, Lmax, dL, d

integer L

	

	do L = Lmin, Lmax, dL
		call benchmark( L, d, numSteps, numThreads )
	end do

end subroutine benchmark_v0_bySize


!-----------------------------------


subroutine benchmark(L, d, numSteps, numThreads)

use flbe
use v0simple
implicit none
integer LX, LY, LZ, L, d, numSteps, numThreads

type(Geometry), pointer :: g
type(Problem), pointer :: prob
class(BoltzmannEquation), pointer :: qbe
class(FD_RightPart), pointer :: rp
type(Efactor), pointer :: F

class(TableOfValues), pointer :: V_q

class(Occupations), pointer :: n, dn0,dn1,dn4
real*8, allocatable :: n0_symm(:,:,:,:), dn0_symm(:,:,:,:), V_q_symm(:,:,:)

integer i, numShow, kx_asymm, ky_asymm, kz_asymm, kx,ky,kz, L2

double precision dt, dv, V0

integer iTimes1, iTimes2, rate
real cpuTimeDelta, r

logical useKroneckers, rta, V_q_dependent

	LX = L
	LY = L
	LZ = L
	if (d.EQ.2) stop 'v0 does not support d=2'
	if (d.EQ.1) stop 'v0 does not support d=1'
	

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
	call prob % updateEquilibriumOccupations()
	
	n => Occupations( g % L(:), 1 )

	n % values(:,:,:,:) = 0.1d0
	n % values(1,0,0,:) = 0.4d0

	L2 = g % L2(1)

	allocate( n0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( dn0_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1, 0:1) )
	allocate( V_q_symm(-L2:L2-1, -L2:L2-1, -L2:L2-1) )
	V0 = 1.0d0
	V_q_symm (:,:,:) = V0
	V_q_dependent = .FALSE.
	rta = .FALSE.
	useKroneckers = .FALSE.

	do kz=-L2, L2-1
		kz_asymm = MOD(kz + L, L)
		do ky=-L2, L2-1
			ky_asymm = MOD(ky + L, L)
			do kx=-L2, L2-1
				kx_asymm = MOD(kx + L, L)

				n0_symm( kx,ky,kz, 0 ) = prob % equilibriumOccupations % values( kx_asymm, ky_asymm, kz_asymm, 0 )

			end do
		end do
	end do

    numShow = numSteps/10
    if (numShow.EQ.0) numShow=1

	! warm-up for FFTW3 late initialization
	write(*, '(A)', advance='no') 'warm-up...'
	call v0_b_c2a2( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta, useKroneckers )
	print*, 'done'
	
    call SYSTEM_CLOCK(iTimes1)
	do i=0, numSteps-1
		call v0_b_c2a2( L, n0_symm, dn0_symm, V_q_symm, V_q_dependent, rta, useKroneckers )
		if (MOD(i+1,numShow).EQ.0) then
			write(*, '(A,I0,A,I0)', advance='no') CHAR(13), i+1, '/', numSteps
		end if
	end do

    call SYSTEM_CLOCK(iTimes2)
    cpuTimeDelta = real(iTimes2-iTimes1)/real(rate)
	write(*, '(A,F8.2,A)') ' DONE in ', cpuTimeDelta, ' seconds.'
	print*	
	print*

	deallocate( V_q_symm )
	deallocate( dn0_symm )
	deallocate( n0_symm )
	deallocate( n )
	deallocate( prob )
	deallocate( g)

	call flbe_done
	open(1, file='benchmark_v0_result.txt', access='append')
	write(1, *) LX, LY, LZ, numThreads, cpuTimeDelta/numSteps
	close(1)

end subroutine benchmark

