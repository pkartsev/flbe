module Subsystem_module

	use Output_module
	use DispersionRelation_module
	use Lifetime_module
	
	implicit none

	public		:: Subsystem

	integer :: subsystem_count = 0
	integer, parameter :: subsystem_numVersions = 6
	character*6, parameter :: operatorLettersFermi='abcdfg'
	character*6, parameter :: operatorLettersBose ='ABCDFG'
	character*6, parameter :: indexLetters='kplnrsuv'

	type Subsystem
		logical		:: boseStatistics
		integer		:: id
		integer		:: occupationIndex
		character*40 :: title
		character :: operatorLetter, indexLetter
		
		class(DispersionRelation), pointer :: dispersionLaw
		type(Lifetime), pointer :: lifetime ! use only if present
		logical lifetimePresent
	contains
		procedure :: setBoseStatistics => Subsystem_setBoseStatistics
		procedure :: setLifetime => Subsystem_setLifetime
		procedure :: setOccupationIndex => Subsystem_setOccupationIndex
		procedure :: setOperatorLetter => Subsystem_setOperatorLetter
		procedure :: getEquilibriumOccupation => Subsystem_getEquilibriumOccupation

		procedure :: toString => Subsystem_toString
		procedure :: show => Subsystem_show

		procedure :: show_N_Operator => Subsystem_show_N_Operator

		final :: Subsystem_destructor
	end type Subsystem

	type Subsystem_p
		class(Subsystem), pointer :: obj
	end type Subsystem_p

! constructors
	interface Subsystem
		module procedure Subsystem_init
	end interface Subsystem
	interface SubsystemFermi
		module procedure SubsystemFermi_init
	end interface SubsystemFermi
	interface SubsystemBose
		module procedure SubsystemBose_init
	end interface SubsystemBose




contains

	integer function Subsystem_nextId() ! single-threaded!  not reenterable,  not atomic

		Subsystem_nextId = subsystem_count ! starting from 0
		subsystem_count = subsystem_count + 1

	end function Subsystem_nextId


! constructors: bodies

	function Subsystem_init( dr ) result(this)
	class(DispersionRelation), pointer :: dr
	integer id
	type(Subsystem), pointer :: this
	
		allocate(this)
		
		id = Subsystem_nextId()
		this % id = id ! to compare

		this % indexLetter = indexLetters( 1:1 )
		
		this % title = 'Subsystem(General)'
		this % boseStatistics = .TRUE.
		this % dispersionLaw => dr
		this % lifetimePresent = .FALSE.
		this % occupationIndex = -1 ! for calculation and show

	end function Subsystem_init
	
	function SubsystemFermi_init( dr ) result(this)
	class(DispersionRelation), pointer :: dr
	type(Subsystem), pointer :: this

		this => Subsystem_init( dr )
		this % operatorLetter = 'x'
		call this % setBoseStatistics( .FALSE. )

	end function SubsystemFermi_init

	function SubsystemBose_init( dr ) result(this)
	class(DispersionRelation), pointer :: dr
	type(Subsystem), pointer :: this
	
		this => Subsystem_init( dr )
		call this % setBoseStatistics( .TRUE. )

	end function SubsystemBose_init

	subroutine Subsystem_destructor(this)
	type(Subsystem) :: this

!		print*, ' Subsystem destructor: ', loc(this)
	
	end subroutine Subsystem_destructor
	


! ----- setters -------

	subroutine Subsystem_setBoseStatistics( this, boseStatistics )
	class(Subsystem), intent(INOUT) :: this
	logical boseStatistics 
	integer id
	
		this % boseStatistics = boseStatistics

		id = MOD(this % id, subsystem_numVersions)

		this % operatorLetter = '?'
		if (boseStatistics) then
			this % title = 'Subsystem(Bose)'
			this % operatorLetter = operatorLettersBose( id+1:id+1 )
		else
			this % title = 'Subsystem(Fermi)'
			this % operatorLetter = operatorLettersFermi( id+1:id+1 )
		end if
		
	end subroutine Subsystem_setBoseStatistics
	
	subroutine Subsystem_setOccupationIndex( this, ind )
	class(Subsystem), intent(INOUT) :: this
	integer ind
	
		this % occupationIndex = ind
		ind = MOD(ind, subsystem_numVersions)

		this % operatorLetter = '?'
		if (this%boseStatistics) then
			this % operatorLetter = operatorLettersBose( ind+1:ind+1 )
		else
			this % operatorLetter = operatorLettersFermi( ind+1:ind+1 )
		end if
		
	end subroutine Subsystem_setOccupationIndex
	
! --	

	subroutine Subsystem_setLifetime( this, lifetime1 )
	class(Subsystem), intent(INOUT) :: this
	type(Lifetime), pointer :: lifetime1

		this % lifetime => lifetime1
		this % lifetimePresent = .TRUE.

	end subroutine Subsystem_setLifetime
	
! --	

	subroutine Subsystem_setOperatorLetter( this, operatorLetter )
	class(Subsystem), intent(INOUT) :: this
	character operatorLetter
	
		this % operatorLetter = operatorLetter
		
	end subroutine Subsystem_setOperatorLetter

!----------------------

	function Subsystem_getEquilibriumOccupation( this, e, T, mu ) result(n0)
	class(Subsystem) :: this
	double precision e, T, mu, n0
	
	double precision arg
	
		n0 = 0.d0
		arg = (e-mu)/T
		if (arg.GT.100) then
			n0 = 0.d0
		else
			if( this % boseStatistics ) then
!				n0 = 1.d0 / (exp( (e-mu)/T ) - 1.d0 )
				if (arg.GT.-100) then
					n0 = 1.d0 / (exp( arg ) - 1.d0 )
				else
					n0 = 1d100
				end if
			else
!				n0 = 1.d0 / (exp( (e-mu)/T ) + 1.d0 )
				if (arg.GT.-100) then
					n0 = 1.d0 / (exp( arg ) + 1.d0 )
				else
					n0 = 1.d0
				end if
			end if
		end if
	
	end function Subsystem_getEquilibriumOccupation


! show() for HTML and TeX output

!--

	subroutine Subsystem_show(this, o )
    class( Subsystem ) :: this
	class( Output ), pointer :: o

    character a, k
    character*10 kbold, ahat, superscript
    
		a = this % operatorLetter
		if (subsystem_count.EQ.1) then
			superscript = ''
		else
			superscript = '(' // a // ')'
		end if		
	    k = this % indexLetter
		call o % hat( a, ahat )
    	call o % bold( k, kbold )
	
		call o % summation()
		call o % beginLimits()
		call o % subscript( kbold )
		call o % endLimits()

		call this % dispersionLaw % showTitle( superscript, o )
		call o % subscript( kbold )

		call o % formulaItem( ahat )
		call o % upperLowerIndex( o % dagger, kbold )
		
		call o % formulaItem( ahat )
		call o % upperLowerIndex( ' ', kbold )
		
	
	end subroutine Subsystem_show

!-- more output -------------------------------------

	subroutine Subsystem_show_N_operator(this, o )
    class( Subsystem ) :: this
	class( Output ), pointer :: o

	character a, k
	character*10 ahat, kbold, a_brackets
	
		a = this % operatorLetter
		k = this % indexLetter
		call o % bold( k, kbold )
		call o % hat( a, ahat )
		a_brackets = '(' // a // ')'
		if (subsystem_count.EQ.1) a_brackets=''
		
    	call o % formulaItem('n' )
    	call o % upperLowerIndex( a_brackets, kbold )
    	call o % formulaItem('=' )

 		call o % openBracket( BRACKET_ANGLE )
    	
    	call o % formulaItem( ahat )
    	call o % upperLowerIndex( o % dagger, kbold )

    	call o % formulaItem( ahat )
    	call o % upperLowerIndex( ' ', kbold )

		call o % closeBracket
    	
		
	end subroutine Subsystem_show_N_operator


	
! toString() useful for information and debug

	character*100 function Subsystem_toString( this )
	class(Subsystem) :: this
	character*20 :: sloc
	character*100 :: ans

		write(sloc, '(I0)' ) loc(this)
		ans = 'Sybsystem@' // TRIM(sloc) // ': ' // TRIM(this % title)
		if ( this % lifetimePresent) then
			ans = trim(ans) // ' {lifetime=' // TRIM(this % lifetime % toString()) // '}'
		end if
		Subsystem_toString  = ans
		return

	end function Subsystem_toString

end module Subsystem_module


