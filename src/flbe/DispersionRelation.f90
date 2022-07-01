module DispersionRelation_module

	use Output_module
	use TableOfValues_module
	use Geometry_module

	implicit none

	integer :: DR_count = 0

	integer, parameter :: DR_CONST = 0
	integer, parameter :: DR_LINEAR = 1
	integer, parameter :: DR_PARABOLIC = 2
	integer, parameter :: DR_TABLE = 3

	type DispersionRelation
		double precision :: coeff ! or value, for CONST
		integer id
		integer mode
		type( TableOfValues ), pointer :: table
		character*40 :: title
		
	contains
		procedure :: constructor => DR_constructor
		procedure :: setCoeff => DR_setCoeff
		procedure :: toString => DR_toString
		procedure :: getEnergy => DR_getEnergy ! beware geometry! convert 0..L-1 -> -L/2 .. L/2
		procedure :: getMinEnergy => DR_getMinEnergy
		procedure :: getMaxEnergy => DR_getMaxEnergy

		procedure :: show => DR_show
		procedure :: showTitle => DR_showTitle
		procedure :: showBody => DR_showBody
		
	end type DispersionRelation

	type DispersionRelation_p
		class(DispersionRelation), pointer :: obj
	end type DispersionRelation_p

! constructors

	interface DispersionRelationConst
		module procedure DR_init_Const
	end interface DispersionRelationConst

	interface DispersionRelationLinear
		module procedure DR_init_LinearWithFormula
		module procedure DR_init_Linear
	end interface DispersionRelationLinear

	interface DispersionRelationParabolic
		module procedure DR_init_ParabolicWithFormula
		module procedure DR_init_Parabolic
	end interface DispersionRelationParabolic

	interface DispersionRelationTable
		module procedure DR_init_TableFromToV
		module procedure DR_init_TableFromDR
	end interface DispersionRelationTable


contains

	integer function DR_nextId() ! single-threaded!  not reenterable,  not atomic

		DR_nextId = DR_count ! starting from 0
		DR_count = DR_count + 1

	end function DR_nextId




! constructors: bodies

	subroutine DR_constructor( this, mode )
	class(DispersionRelation), intent(INOUT) :: this
	integer mode
	
		this % id = DR_nextId()
		this % mode = mode
		
	end subroutine DR_constructor

!--

	function DR_init_Const( c ) result(this)
	double precision c	
	type(DispersionRelation), pointer :: this
	
		allocate(this)
		call this % constructor( DR_CONST )
		this % coeff = c
		this % title = 'DispersionRelation(Const)'

	end function DR_init_Const

!--

	function DR_init_LinearWithFormula() result(this)
	type(DispersionRelation), pointer :: this
	
		allocate(this)
		call this % constructor( DR_LINEAR )
		this % coeff = 1.d0
		this % title = 'DispersionRelation(Linear)'

	end function DR_init_LinearWithFormula

	function DR_init_Linear( sizes, coeff ) result(this)
	type(DispersionRelation), pointer :: this, dr
	type( Geometry ), pointer :: sizes	
	double precision, optional :: coeff
	
		dr => DR_init_LinearWithFormula()
		if (present(coeff)) call dr % setCoeff( coeff )
		this => DispersionRelationTable( dr, sizes )
		deallocate(dr)
!		print*, ' size=', this%table%L(:)
!		print*, ' energies=', this%table%values(:,:,:)
!@stop

	end function DR_init_Linear

!--

	function DR_init_ParabolicWithFormula() result(this)
	type(DispersionRelation), pointer :: this
	
		allocate(this)
		call this % constructor( DR_PARABOLIC )

		this % coeff = 1.d0
		this % title = 'DispersionRelation(Parabolic)'

	end function DR_init_ParabolicWithFormula

!--

	function DR_init_Parabolic( sizes,  coeff  ) result(this)
	type(DispersionRelation), pointer :: this, dr
	type( Geometry ), pointer :: sizes	
	double precision, optional :: coeff
	
		dr => DR_init_ParabolicWithFormula()
		if (present(coeff)) call dr % setCoeff( coeff )
		this => DispersionRelationTable( dr, sizes )
		deallocate(dr)

	end function DR_init_Parabolic

!--

	function DR_init_TableFromToV( tov ) result(this)
	class( TableOfValues ), pointer :: tov
	type(DispersionRelation), pointer :: this
	
		allocate(this)
		call this % constructor( DR_TABLE )
		this % coeff = 1.d0 ! not used though
		this % title = 'DispersionRelation(Table)'
		this % table => tov
		
	end function DR_init_TableFromToV

	function DR_init_TableFromDR( dr, sizes ) result(this)
	type(DispersionRelation), pointer :: this, dr
	type(Geometry), pointer :: sizes
	
	integer k(1:3), kx, ky, kz, LX, LY, LZ

		allocate(this)
		call this % constructor( DR_TABLE )
		this % coeff = 1.d0 ! not used though
		this % title = 'DispersionRelation(Table) from ' // TRIM(dr%title)
		this % table => TableOfValues( sizes % L )

		LX = sizes % L(1)
		LY = sizes % L(2)
		LZ = sizes % L(3)
		do kz=0, LZ-1
			k(3) = kz
			do ky=0, LY-1
				k(2) = ky
				do kx=0, LX-1
					k(1) = kx
!					call this % table % setValue( k, dr % getEnergy( k, sizes ))
					this % table % values( kx, ky, kz ) = dr % getEnergy( k, sizes )
				end do		
			end do		
		end do		
	end function DR_init_TableFromDR

! setter for coefficient/value

	subroutine DR_setCoeff( this, c )
	class(DispersionRelation), intent(INOUT) :: this
	double precision :: c

		this % coeff = c

	end subroutine DR_setCoeff



! toString() useful for information and debug

	character*80 function DR_toString( this )
	class(DispersionRelation), intent(INOUT) :: this
	character*21 :: s

		write(s, '(ES21.15)' ) this % coeff
		DR_toString = TRIM(this % title) // ' (coeff=' // TRIM(s) // ')'

	end function DR_toString

! main function : getEnergy (kx, ky, kz)

	double precision function DR_getEnergy( this, k_in, sizes )
	class(DispersionRelation), intent(INOUT) :: this
	integer k_in(1:3)
	type(Geometry), pointer :: sizes
	
	integer k(1:3)
	double precision v

!		print*, ' DR.getEnergy: k=', k	
		select case (this % mode) 
			case (DR_CONST)
				v = 1.d0
			case (DR_LINEAR)
				call sizes % symmetrize( k_in, k )
				v = sqrt( 0.d0 + sum( k(:)**2 ) )
			case (DR_PARABOLIC)
				call sizes % symmetrize( k_in, k )
				v = ( 0.d0 + sum( k(:)**2 ) )
			case (DR_TABLE)
				v = this % table % getValue( k_in )
		end select
		DR_getEnergy = v * this % coeff
		return
	end function DR_getEnergy

!--

	double precision function DR_getMaxEnergy( this, sizes )
	class(DispersionRelation), intent(INOUT) :: this
	type(Geometry), pointer :: sizes
	
	integer kmax(1:3)
	double precision v

		select case (this % mode) 
			case (DR_CONST)
				v = 1.d0
			case (DR_LINEAR)
				if (this%coeff.GT.0) then
					kmax(:) = sizes % L2(:)
					v = sqrt( 0.d0 + sum( kmax(:)**2 ) )
				else
					v = 0.d0
				end if
			case (DR_PARABOLIC)
				if (this%coeff.GT.0) then
					kmax(:) = sizes % L2(:)
					v = 0.d0 + sum( kmax(:)**2 ) 
				else
					v = 0.d0
				end if
			case (DR_TABLE)
				v = MAXVAL( this % table % values(:,:,:) )
		end select
		DR_getMaxEnergy = v * this % coeff
		return
	end function DR_getMaxEnergy
	
	double precision function DR_getMinEnergy( this, sizes )
	class(DispersionRelation), intent(INOUT) :: this
	type(Geometry), pointer :: sizes
	
	integer kmax(1:3)
	double precision v

		select case (this % mode) 
			case (DR_CONST)
				v = 1.d0
			case (DR_LINEAR)
				if (this%coeff.LT.0) then
					kmax(:) = sizes % L2(:)
					v = sqrt( 0.d0 + sum( kmax(:)**2 ) )
				else
					v = 0.d0
				end if
			case (DR_PARABOLIC)
				if (this%coeff.LT.0) then
					kmax(:) = sizes % L2(:)
					v = 0.d0 + sum( kmax(:)**2 ) 
				else
					v = 0.d0
				end if
			case (DR_TABLE)
				v = MINVAL( this % table % values(:,:,:) )
		end select
		DR_getMinEnergy = v * this % coeff
		return
	end function DR_getMinEnergy

!--

	subroutine DR_show( this, superscript, o )
	class(DispersionRelation) :: this
	character(len=*) superscript
	class( Output ), pointer :: o
	character k
	character*10 kbold
	
		k = 'k'
		call o % bold( k, kbold )

		call o % formulaItem( o % varepsilon )
		call o % upperLowerIndex( superscript, kbold )
		call o % formulaItem( '=' )
		call this % showBody( 'k', o ) 

	end subroutine DR_show


	subroutine DR_showTitle( this, superscript, o )
	class(DispersionRelation) :: this
	character(len=*) superscript
	class( Output ), pointer :: o
	
		call o % formulaItem( o % varepsilon )
		call o % superscript( superscript )
	
	end subroutine DR_showTitle
	
	
	subroutine DR_showBody( this, k, o )
	class(DispersionRelation) :: this
	character k
	class( Output ), pointer :: o
	
	character*10 id_str, kbold
	character*40 coeff_str
	
		write(id_str, '(I0)') this % id
		write(coeff_str, '(F10.1)') this % coeff
		call o % bold( k, kbold )
		
		select case(this % mode) 
			case (DR_CONST)
				call o % formulaItem( coeff_str )
			case (DR_LINEAR)
				call o % formulaItem( 'c' )
				call o % upperLowerIndex( ' ', id_str )
				call o % formulaItem( '|' // TRIM(kbold) // '|' )
			case (DR_PARABOLIC)
				call o % formulaItem( 'c' )
				call o % upperLowerIndex( ' ', id_str )
				call o % formulaItem( kbold )
				call o % upperLowerIndex( '2', ' ' )
			case (DR_TABLE)
				call o % formulaItem( o % varepsilon )
				call o % upperLowerIndex( ' ', id_str )
				call o % formulaItem( '(' // TRIM(kbold) // ')' )
		end select

	
	end subroutine DR_showBody

end module DispersionRelation_module


