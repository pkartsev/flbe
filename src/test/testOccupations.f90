subroutine testOccupations
use flbe
implicit none

integer dims, subs
	write(*, '(A)', advance='no') ' [ test Occupations...'
	
	do dims=2,3
		do subs=1,3
!			print*, dims, subs
			call testOccupationsDimsSubs( dims, subs )
		end do
	end do
	print*, ' OK ]'
	
end subroutine testOccupations
	
subroutine testOccupationsDimsSubs( dims, numSubsystems )

use flbe
implicit none
integer dims, numSubsystems

type(Geometry), pointer :: geom
type(Occupations), pointer :: n1,n2,n3
type(Subsystem), pointer :: s1, s2, s3
type(Hamiltonian), pointer :: h
class(DispersionRelation), pointer :: dr

integer kx,ky,kz, LX, i, volume
integer L(1:3), k(1:3), ksym(1:3)
double precision dv
real r
logical bad
	
	bad = .FALSE.

	LX=8
	select case (dims)
		case(2)
			geom => Geometry(LX,LX,1)
		case(3)
			geom => Geometry(LX,LX,LX)
	end select
	L(:) = geom % L(:)
	
	dr => DispersionRelationParabolic(  )

	h => Hamiltonian()
!	do i=0, numSubsystems-1
		s1 => SubsystemFermi( dr )
		call h % addSubsystem(s1)
!	print*, '-1-'
	if (numSubsystems.GT.1) then
!	print*, '-2-'
		s2 => SubsystemFermi( dr )
!	print*, '-3-'
		call h % addSubsystem(s2)
	end if
!	print*, '-4-'
	if (numSubsystems.GT.2) then
		s3 => SubsystemFermi( dr )
		call h % addSubsystem(s3)
	end if
!	end do

	volume = L(1) * L(2) * L(3) * numSubsystems
	
	n1 => Occupations( L, h%numSubsystems )	! type 1 of constructor
	n2 => Occupations( L, numSubsystems )
	n3 => Occupations( n2 ) 				! type 2 of constructor

	if (n1 % N4 .NE. numSubsystems) stop 'n1 % N4'
	
	do i=0, numSubsystems-1
		if (dims.EQ.2) then
			do ky=0, LX-1
			do kx=0, LX-1
				call RANDOM_NUMBER(r)
				n1 % values(kx,ky,0,i) = r + 2
			end do
			end do
		end if
		if (dims.EQ.3) then
			do kz=0, LX-1
			do ky=0, LX-1
			do kx=0, LX-1
				call RANDOM_NUMBER(r)
				n1 % values(kx,ky,kz,i) = r + 2
			end do					
			end do					
			end do					
		end if
	end do
!	print*, '-10-'
	
	n2 % values(:,:,:,:) = n1 % values(:,:,:,:)
	n3 % values(:,:,:,:) = n1 % values(:,:,:,:)

	dv = n1 % compareTo(n2)
	if (dv.GT.1d-12) stop 'equal arrays compareTo() gives nonzero value'

	call RANDOM_NUMBER(r)
	r = r + 1.0
	n2 % values(0,0,0,0) = 0.d0
	n3 % values(0,0,0,0) = r
	dv = n2 % compareTo(n3)
!	print*, ' dv=', dv, '  (r**2)=',  (r**2)
	dv = dv - (r**2) 
!	print*, ' dv=', dv
	if (abs(dv).GT.1d-4) stop 'changed single value compareTo() gives wrong distance'

!	print*, '-20-'
	call RANDOM_NUMBER(r)
	r = r + 1.0
	n2 % values(:,:,:,:) = r
	n3 % values(:,:,:,:) = 0.d0
	dv = n2 % compareTo(n3)
!	print*, ' dv=', dv, '  (r**2)*volume=',  (r**2)*volume, ' volume=', volume
	dv = dv - (r**2)*volume
!	print*, ' dv=', dv
	if (abs(dv).GT.1d-2) stop 'probably not all volume compared'	
	
!	print*, '-21-'

	deallocate(h)
!	print*, '-22-'

	if (numSubsystems.GT.2) deallocate(s3)
	if (numSubsystems.GT.1) deallocate(s2)

!	print*, '-23-'
	deallocate(s1)
!	print*, '-24-'
	deallocate(dr)
!	print*, '-25-'
	deallocate(geom)
!	print*, '-26-'

	if (bad) stop 'FAILED'

end subroutine testOccupationsDimsSubs
