module TableOfValues_module
	implicit none

	type TableOfValues
		integer, dimension(1:3) :: L ! sizes
		double precision, allocatable :: values(:,:,:)
	contains
		procedure :: getValue => ToV_getValue

		final :: ToV_destructor
	end type TableOfValues

! constructors

	interface TableOfValues
		module procedure TableOfValues_initEmpty
!		module procedure TableOfValues_initWithValues
	end interface TableOfValues


contains

! constructor

	function TableOfValues_initEmpty( L ) result(this)
	integer, dimension(1:3) :: L ! sizes
	integer LX, LY, LZ
	
	type(TableOfValues), pointer :: this
		
		allocate(this)
		
		this % L(:) = L(:)
		LX = L(1)
		LY = L(2)
		LZ = L(3)
		allocate( this % values(0:LX-1, 0:LY-1, 0:LZ-1) )
		
	end function TableOfValues_initEmpty

!-- destructor

	subroutine ToV_destructor( this )
	type( TableOfValues ) :: this
	
!		print*, 'TableOfValues: destructor'
		deallocate( this % values )
		
	end subroutine ToV_destructor

!--

	double precision function ToV_getValue( this, k )
	class( TableOfValues ) :: this
	integer, dimension(1:3) :: k

		if ((k(1).GT.100).OR.(k(2).GT.100).OR.(k(3).GT.100)) then
			print*, ' ToV.getValue for k=', k
			stop
		end if

!		print*, ' ToV_getValue (',  k(1), k(2), k(3) , ')'
		ToV_getValue = this % values( k(1), k(2), k(3) )
		
	end function ToV_getValue

end module TableOfValues_module


