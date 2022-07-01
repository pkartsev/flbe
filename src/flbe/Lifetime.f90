module Lifetime_module

	use Output_module
	use TableOfValues_module

	implicit none

	integer :: lifetimes_count = 0

! modes, instead of polymorphism, can change type on the fly
	integer, parameter :: LIFETIME_CONST = 0
	integer, parameter :: LIFETIME_TABLE = 3

	type Lifetime
		integer id
		double precision :: value
		integer mode
		type( TableOfValues ), pointer :: table
		character*40 :: title
		
	contains
		procedure :: setValue => Lifetime_setValue
		procedure :: getValue => Lifetime_getValue
		procedure :: toString => Lifetime_toString
		procedure :: show => Lifetime_show
	end type Lifetime

! constructors

	interface LifetimeConst
		module procedure LifetimeConst_init
	end interface LifetimeConst

	interface LifetimeTable
		module procedure LifetimeTable_init
	end interface LifetimeTable


contains

	integer function Lifetime_nextId() ! single-threaded!  not reenterable,  not atomic

		Lifetime_nextId = lifetimes_count ! starting from 0
		lifetimes_count = lifetimes_count + 1

	end function Lifetime_nextId


! constructors: bodies

	function LifetimeConst_init(c) result(this)
	double precision c
	type(Lifetime), pointer :: this

		allocate(this)

		this % id = Lifetime_nextId()
		this % mode = LIFETIME_CONST
		this % value = c
		this % title = 'Lifetime(Const)'

	end function LifetimeConst_init

	function LifetimeTable_init( table ) result(this)
	class( TableOfValues ), pointer :: table
	type(Lifetime), pointer :: this
	
		allocate(this)
	
		this % mode = LIFETIME_TABLE
		this % value = 1.d0
		this % table = table
		this % title = 'Lifetime(Table)'
		
	end function LifetimeTable_init


! setter for coefficient/value

	subroutine Lifetime_setValue( this, c )
	class(Lifetime), intent(INOUT) :: this
	double precision :: c

		if (this%mode .NE. LIFETIME_CONST) stop 'Lifetime_setValue wrong mode'
		
		this % value = c

	end subroutine Lifetime_setValue



! toString() useful for information and debug

	character*80 function Lifetime_toString( this )
	class(Lifetime) :: this
	character*21 :: s
		
		select case (this % mode) 
			case (LIFETIME_CONST)
				write(s, '(ES21.15)' ) this % value
				Lifetime_toString = TRIM(this % title) // ' [value=' // TRIM(s) // ']'
			case (LIFETIME_TABLE)
				Lifetime_toString = TRIM(this % title) 
		end select

	end function Lifetime_toString

!--

	subroutine Lifetime_show(this, version, fd, superscript, subscript)
    class( Lifetime ) :: this
    integer version
    integer fd
    character superscript, subscript

! nothing

	end subroutine Lifetime_show    
	
!--
		
	double precision function Lifetime_getValue( this, k )
	class(Lifetime) :: this
	integer, dimension(1:3) :: k
	double precision v
	
		select case (this % mode) 
			case (LIFETIME_CONST)
				v = this % value
			case (LIFETIME_TABLE)
				v = this % table % getValue( k )
		end select
		
		Lifetime_getValue = v
		return
		
	end function Lifetime_getValue

!--

end module Lifetime_module


