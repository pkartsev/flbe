subroutine testGeometry

use flbe
implicit none

type(Geometry), pointer :: g

integer, dimension(1:3) :: a, b, c, apb, bmc, z, minusb, anorm
integer LX, LY, LZ, i, L
logical bad, good

	write(*, '(A)', advance='no') ' [ test Geometry...'

	LX = 4
	LY = 4
	LZ = 4

	g => Geometry( LX, LY, LZ )

	a = (/ 1,3,2 /)
	b = (/ 2,1,3 /)
	c = (/ 1,3,3 /)

	call g % plus1 ( a, b, apb )
	call g % minus1( b, c, bmc )
	call g % zero( z )

!	print*, 'a=', a
!	print*, 'b=', b
!	print*, 'c=', c
!	print*, 'a+b=', apb
!	print*, 'bmc=', bmc
	

	minusb(:) = b(:)
	call g % minus( minusb )
	
	bad = .FALSE.
	do i=1,3
		L = g%L(i)
		! apb = a+b mod L
		if ( apb(i) .LT. 0 ) bad = .TRUE.
		if ( apb(i) .GE. L ) bad = .TRUE.
		if ( MOD( 2*L + apb(i) - a(i) - b(i), L ) .NE. 0 ) bad = .TRUE.

		! bmc = b-c mod L   in range 0..L-1
		if ( bmc(i) .LT. 0 ) bad = .TRUE.
		if ( bmc(i) .GE. L ) bad = .TRUE.
		if ( MOD( 2*L + bmc(i) - b(i) + c(i), L ) .NE. 0 ) bad = .TRUE.

		! zero
		if ( z(i) .NE. 0) bad = .TRUE.
		! -b -> mod(L-b, L)
		if ( MOD(minusb(i)+b(i), L) .NE. 0) bad = .TRUE.
		
	end do	

	L = LX
	do i=0, L-1 ! symmetrize: 0...L-1 -> -L/2 .. L/2-1
		a(:) = i
		call g % symmetrize( a, anorm )
		if ( i .LT. L/2 ) then
!			if ( anorm(1) .NE. a(1) ) bad = .TRUE.
		else
!			if ( anorm(1) .NE. a(1)-L ) bad = .TRUE.
		end if
	end do

	deallocate(g)

	if (bad) stop 'FAILED'
	print*, ' OK ]'
	
end subroutine testGeometry
