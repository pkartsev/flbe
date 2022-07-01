module Occupations_module

	implicit none
	
	type Occupations
		integer, dimension(1:3) :: L
		integer N4
		double precision, pointer :: values(:,:,:,:) ! x,y,z, index
	contains
		procedure :: getValue => Occupations_getValue
		procedure :: setValue => Occupations_setValue
		procedure :: incValue => Occupations_incValue
		procedure :: add1 => Occupations_add1
		procedure :: add2 => Occupations_add2
		procedure :: add3 => Occupations_add3
		procedure :: add4 => Occupations_add4
		procedure :: add5 => Occupations_add5
		procedure :: write2D => Occupations_write2D
		procedure :: write3D => Occupations_write3D

		procedure :: compareTo => Occupations_compareTo
		procedure :: compareToRelative => Occupations_compareToRelative

		final :: Occupations_destructor
		
	end type Occupations
	
	interface Occupations
		module procedure Occupations_init4D
		module procedure Occupations_initSimilar
	end interface Occupations
	
	
contains

	subroutine Occupations_destructor( this )
	type( Occupations ) :: this

!		print*, ' Occupations: destructor'
		deallocate( this % values )
		
		
	end subroutine Occupations_destructor

! --
	
	function Occupations_init4D( L, N4 ) result(this)
	integer L(1:3), N4
	type( Occupations ), pointer :: this

	integer LX, LY, LZ

		allocate( this )

		this % N4 = N4
		this % L(:) = L(:)
				
		LX = L(1)
		LY = L(2)
		LZ = L(3)

		allocate( this % values( 0:LX-1,  0:LY-1,  0:LZ-1, 0:N4-1 ) )
		
	end function Occupations_init4D
	
! --

	function Occupations_initSimilar( occ0 ) result(this)
	type( Occupations ), pointer :: occ0
	type( Occupations ), pointer :: this

		this => Occupations_init4D( occ0 % L(:), occ0 % N4 )
		
	end function Occupations_initSimilar

! --

	subroutine Occupations_setValue( this, k, n4, v )
	class( Occupations ), intent(INOUT) :: this
	integer, dimension(1:3) :: k
	integer n4
	double precision v
	
		this % values( k(1), k(2), k(3), n4 ) = v
		
	end subroutine Occupations_setValue

! --

	subroutine Occupations_incValue( this, k, n4, v )
	class( Occupations ), intent(INOUT) :: this
	integer, dimension(1:3) :: k
	integer n4
	double precision v
	
		this % values( k(1), k(2), k(3), n4 ) = this % values( k(1), k(2), k(3), n4 ) + v

	end subroutine Occupations_incValue

! --

	double precision function Occupations_getValue( this, k, n4 )
	class( Occupations ) :: this
	integer, dimension(1:3) :: k
	integer n4
	
		Occupations_getValue = this % values( k(1), k(2), k(3), n4 )
		
	end function Occupations_getValue

!---------------------------------------------------


	subroutine Occupations_add1( n_from, dn, n_to, factor )
	class( Occupations ) :: n_from
	type( Occupations ), pointer :: n_to, dn
	double precision :: factor

	integer x,y,z,s
		
if (1.EQ.1) then

		n_to % values(:,:,:,:) =  n_from % values(:,:,:,:) + &
&			dn % values(:,:,:,:) * factor

else
do s=0, n_from%N4-1
!$omp parallel do private(x,y,z)
	do z=0, n_from%L(3)-1
	do y=0, n_from%L(2)-1
	do x=0, n_from%L(1)-1
		n_to % values(x,y,z,s) =  n_from % values(x,y,z,s) + &
&			dn % values(x,y,z,s) * factor
	end do
	end do
	end do
!$omp end parallel do
	end do
end if
	
	end subroutine Occupations_add1
	
	subroutine Occupations_add2( n_from, dn1, dn2, n_to, factor1, factor2, dt )
	class( Occupations ) :: n_from
	type( Occupations ), pointer :: dn1, dn2, n_to
	double precision :: factor1, factor2, dt
	
		n_to % values(:,:,:,:) =  n_from % values(:,:,:,:) + &
&			( &
&				  dn1 % values(:,:,:,:) * factor1  &
&				+ dn2 % values(:,:,:,:) * factor2  &
&			) * dt
	end subroutine Occupations_add2

	subroutine Occupations_add3( n_from, dn1, dn2, dn3, n_to, factor1, factor2, factor3, dt )
	class( Occupations ) :: n_from
	type( Occupations ), pointer :: dn1, dn2, dn3, n_to
	double precision :: factor1, factor2, factor3, dt
	
		n_to % values(:,:,:,:) =  n_from % values(:,:,:,:) + &
&			( &
&				  dn1 % values(:,:,:,:) * factor1  &
&				+ dn2 % values(:,:,:,:) * factor2  &
&				+ dn3 % values(:,:,:,:) * factor3  &
&			) * dt
	end subroutine Occupations_add3
	
	subroutine Occupations_add4( n_from, dn1, dn2, dn3, dn4, &
&		n_to, factor1, factor2, factor3, factor4, dt )
	class( Occupations ) :: n_from
	type( Occupations ), pointer :: dn1, dn2, dn3, dn4, n_to
	double precision :: factor1, factor2, factor3, factor4, dt
	
		n_to % values(:,:,:,:) =  n_from % values(:,:,:,:) + &
&			( &
&				  dn1 % values(:,:,:,:) * factor1  &
&				+ dn2 % values(:,:,:,:) * factor2  &
&				+ dn3 % values(:,:,:,:) * factor3  &
&				+ dn4 % values(:,:,:,:) * factor4  &
&			) * dt
	end subroutine Occupations_add4
	
	subroutine Occupations_add5( n_from, dn1, dn2, dn3, dn4, dn5, &
&		n_to, factor1, factor2, factor3, factor4, factor5, dt )
	class( Occupations ) :: n_from
	type( Occupations ), pointer :: dn1, dn2, dn3, dn4, dn5, n_to
	double precision :: factor1, factor2, factor3, factor4, factor5, dt
	
		n_to % values(:,:,:,:) =  n_from % values(:,:,:,:) + &
&			( &
&				  dn1 % values(:,:,:,:) * factor1  &
&				+ dn2 % values(:,:,:,:) * factor2  &
&				+ dn3 % values(:,:,:,:) * factor3  &
&				+ dn4 % values(:,:,:,:) * factor4  &
&				+ dn5 % values(:,:,:,:) * factor5  &
&			) * dt
	end subroutine Occupations_add5
	
	subroutine Occupations_write2D( this, fd, zs, s )	
	class( Occupations ) :: this
	integer fd, zs, s
	integer x,y, xs,ys, z, Lx, Ly, Lz, Lx2, Ly2, Lz2

		Lx = this % L(1)
		Ly = this % L(2)
		Lz = this % L(3)
		Lx2 = Lx / 2
		Ly2 = Ly / 2
		Lz2 = Lz / 2

		z = MOD(zs+Lz, Lz)

		do ys=-Ly2,Ly2-1 ! L=8: -4..3 ; L=1 : 0..-1 :(
			y = MOD(ys+Ly, Ly)
			write( fd, *) (this % values(MOD(xs+Lx, Lx),y,z,s), xs=-Lx2,Lx2-1 )
		end do

	end subroutine Occupations_write2D
	
	
	subroutine Occupations_write3D( this, filename)	
	class( Occupations ) :: this
	character(len=*) :: filename

	integer x,y,z, xs,ys,zs, Lx, Ly, Lz, Lx2, Ly2, Lz2, s

	Lx = this % L(1)
	Ly = this % L(2)
	Lz = this % L(3)
	Lx2 = Lx / 2
	Ly2 = Ly / 2
	Lz2 = Lz / 2
	
	open(1, file=filename)
	do s=0, this%N4-1
	do xs=-Lx2,Lx2-1 ! L=8: -4..3 ; L=1 : 0..-1 :(
		x = MOD(xs+Lx, Lx)
		do ys=-Ly2,Ly2-1
			y = MOD(ys+Ly, Ly)
			do zs=-Lz2,Lz2-1
				z = MOD(zs+Lz, Lz)
				write(1, *) s, xs,ys,zs, this % values(x,y,z,s)
			end do
		end do
	end do
	end do
	close(1)

	end subroutine Occupations_write3D
	
	
!-------------

	double precision function Occupations_compareTo( this, occ2 )
	class( Occupations ) :: this
	type( Occupations ), pointer :: occ2
	
	double precision s
	
		s = sum( (this % values(:,:,:,:) - occ2 % values(:,:,:,:) ) **2 )
		Occupations_compareTo = s
		
	end function Occupations_compareTo



	double precision function Occupations_compareToRelative( this, occ2 )
	class( Occupations ) :: this
	type( Occupations ), pointer :: occ2
	
	double precision numer, denom
	
		numer = this % compareTo( occ2 )
		denom = sum( this % values(:,:,:,:)**2) + sum( occ2 % values(:,:,:,:)**2 )
		Occupations_compareToRelative = numer / ( denom + 1.0d-10 )
		
	end function Occupations_compareToRelative

end module Occupations_module
