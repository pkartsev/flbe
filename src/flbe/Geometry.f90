module Geometry_module

	use Output_module
	implicit none
	
	type Geometry
		integer, dimension(1:3) :: L, L2
		integer volume
	contains
		procedure :: toString => Geometry_toString
		procedure :: show => Geometry_show
		procedure :: toVector => Geometry_toVector
		procedure :: fromVector => Geometry_fromVector

		procedure :: zero => Geometry_zero
		procedure :: minus => Geometry_minus
		procedure :: add1 => Geometry_add1
		procedure :: sub1 => Geometry_sub1
		procedure :: minus1 => Geometry_minus1
		procedure :: plus1 => Geometry_plus1

		procedure :: add1sub1 => Geometry_add1sub1

		procedure :: check => Geometry_check
		procedure :: symmetrize => Geometry_symmetrize ! convert 0..L-1 to -L/2..L/2-1

		final :: Geometry_destructor

	end type Geometry
	
	interface Geometry
		module procedure Geometry_init
	end interface Geometry
	
	
contains

	subroutine Geometry_destructor( this )
	type( Geometry ) :: this

!		print*, ' Geometry: destructor'
		
	end subroutine Geometry_destructor



	function Geometry_init(LX, LY, LZ) result(this)
	integer LX, LY, LZ
	
	type(Geometry), pointer :: this
	
		allocate(this)
	
		this % L(1) = LX
		this % L(2) = LY
		this % L(3) = LZ
		this % L2(1) = LX/2
		this % L2(2) = LY/2
		this % L2(3) = LZ/2
		this % volume = LX*LY*LZ
		
	end function Geometry_init
		
!--

	character*100 function Geometry_toString( this )
	class(Geometry) :: this
	
	character*10 s1,s2,s3
	
		write(s1, '(I3)') this % L(1)
		write(s2, '(I3)') this % L(2)
		write(s3, '(I3)') this % L(3)
		Geometry_toString = 'Geometry : ' // TRIM(s1) // 'x' //TRIM(s2) // 'x' //TRIM(s3)
		
	end function Geometry_toString
		
!--

	subroutine Geometry_toVector(this, i, v )
	class(Geometry) :: this
	integer, intent(IN) :: i
	integer, dimension(1:3), intent(OUT) :: v

		v(1) = MOD( i, this % L(1) )
		v(2) = MOD( i / this % L(1), this % L(2) )
		v(3) = i / (this % L(1) * this % L(2) )

	end subroutine Geometry_toVector

	integer function Geometry_fromVector(this, v )
	class(Geometry) :: this
	integer, dimension(1:3), intent(IN) :: v

		Geometry_fromVector = v(1) + this % L(1) * ( v(2) + this % L(2) * v(3) )

	end function Geometry_fromVector

!--

	subroutine Geometry_print( v )
	integer, dimension(1:3), intent(IN) :: v
		write(*, '(A,I0,A,I0,A,I0,A)', advance='no')  '(', v(1),',',v(2),',',v(3),')'
	end subroutine

!--

	subroutine Geometry_check(this, v )
	class(Geometry) :: this
	integer, dimension(1:3), intent(IN) :: v

	integer i
	
		do i=1, 3
			if (v(i).GE. this % L(i)) stop 'v GE L'
			if (v(i).LT. 0) stop 'v LT 0'
		end do

	end subroutine Geometry_check

!--

	subroutine Geometry_symmetrize(this, v_in, v_out ) ! convert from 0..L-1 to -L/2..L/2-1
	class(Geometry) :: this
	integer, dimension(1:3), intent(IN) :: v_in
	integer, dimension(1:3), intent(OUT) :: v_out

	integer i, L, L2
	
		v_out(:) = MOD( v_in(:) + this%L2(:), this%L(:)) - this%L2(:)
		do i=1, 3
			L = this % L(i)
			L2 = L/2
			if ((L2.GT.0).AND.(v_out(i).GE.L2)) v_out(i) = v_out(i) - L
			if (v_out(i) .NE. (MOD(v_in(i)+L2, L)-L2) ) stop 'symmetrize fail'
		end do
 
	end subroutine Geometry_symmetrize

!--

	subroutine Geometry_add1( this, v, dv )
	class(Geometry) :: this
	integer, dimension(1:3), intent(INOUT) :: v
	integer, dimension(1:3), intent(IN) :: dv
	
		v(:) = MOD( v(:) + dv(:), this % L(:) )
		
	end subroutine Geometry_add1

	subroutine Geometry_plus1( this, v, dv, res )
	class(Geometry) :: this
	integer, dimension(1:3), intent(IN) :: v,dv
	integer, dimension(1:3), intent(OUT) :: res
	
		res(:) = MOD( v(:) + dv(:), this % L(:) )

	end subroutine Geometry_plus1

!--

	subroutine Geometry_sub1( this, v, dv )
	class(Geometry) :: this
	integer, dimension(1:3), intent(INOUT) :: v
	integer, dimension(1:3), intent(IN) :: dv
	
		v(:) = MOD( v(:) + this % L(:) - dv(:), this % L(:) )
		
	end subroutine Geometry_sub1

	subroutine Geometry_minus1( this, v, dv, res )
	class(Geometry) :: this
	integer, dimension(1:3), intent(IN) :: v,dv
	integer, dimension(1:3), intent(OUT) :: res
	
		res(:) = MOD( v(:) + this % L(:) - dv(:), this % L(:) )
		
	end subroutine Geometry_minus1

!--

	subroutine Geometry_add1sub1( this, v, a1, s1 )
	class(Geometry) :: this
	integer, dimension(1:3), intent(INOUT) :: v
	integer, dimension(1:3), intent(IN) :: a1, s1
	
		v(:) = MOD( v(:) + this % L(:) + a1(:) - s1(:), this%L(:) )
		
	end subroutine Geometry_add1sub1

!--
	subroutine Geometry_minus( this, v )
	class(Geometry) :: this
	integer, dimension(1:3), intent(INOUT) :: v
	
		v(:) = MOD( this %L(:)-v(:), this%L(:) )

	end subroutine Geometry_minus
	
!--
	
	subroutine Geometry_zero( this, v )
	class(Geometry) :: this
	integer, dimension(1:3), intent(OUT) :: v
	
		v(:) = 0
		
	end subroutine Geometry_zero
	
! -- for debug

	subroutine Geometry_show( this, o )
	class(Geometry) :: this
	class( Output ), pointer :: o

	character*100 info
	
		write(info, '(I0,A,I0,A,I0,A)') &
&			this%L(1), TRIM(o%times), &
&			this%L(2), TRIM(o%times), &
&			this%L(3)

		call o % write ( 'System size: ' )
		call o % beginEquation()
		call o % formulaItem ( info )
		call o % endEquation()
		
	end subroutine Geometry_show
	
	
end module Geometry_module
