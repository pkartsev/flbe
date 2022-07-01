module Perturbation_module

use Subsystem_module
use Output_module
implicit none

! operator types
integer, parameter :: CREATION = 0
integer, parameter :: ANNIHILATION = 1

!interaction amplitude version
integer, parameter :: Q_INDEPENDENT = 0
integer, parameter :: Q_DEPENDENT = 1

integer, parameter :: max_num_subsystems = 4
integer, parameter :: max_num_operators = max_num_subsystems * 4


integer ops_count

type OperatorInfo ! pre-created info for show(), qbe terms
	 ! 0/1 HTML/TeX
	 ! op = operator (a, b, ...)   /   n = occupation number a(+) a
	character*64, pointer :: subscript(:)
	character*64, pointer :: superscript(:)
	integer :: pm ! +1 or -1 (or 0?)
	integer subsystemIndex ! index in this perturbation % subsystems(:)  
	character a, k
	integer indexPos
	logical dagger, bose, firstOperatorForSubsystem
	logical free ! not fixed like k4 in formula delta( k1 + k2 = k3 + k4 )
	! a+a+ a a -> 2/2, 2/2 ;  a+ a b -> 1/1, 1/1, 1/1 ; a+ a+ a -> 2/1, 2/1, 1/2
	
	type(Subsystem), pointer :: subsystem ! to get ID for calculation
	integer numPossibleKroneckers

end type

type Perturbation ! Interaction term

    integer id
    character(80) title
    integer numSubsystems
    integer currentSubsystem
	
	integer mode ! amplitude: V(q) or V0
	double precision amplitudeValue	
	type( TableOfValues ), pointer :: amplitudeTable
	
	type( Subsystem_p ), allocatable :: subsystems(:)
	
	integer numOperators(0:1, 0:max_num_subsystems-1) ! 1/0 = creation/annihilation
	integer, allocatable :: firstOperatorEntry(:) ! subsystemIndex
	integer, allocatable :: lastOperatorEntry(:)  ! subsystemIndex
	logical operatorSpinUp(0:max_num_subsystems-1) ! T/F for fermi, ignored for Bose
	type(OperatorInfo), allocatable :: info(:)
	character*300, pointer :: lastOperatorSubscript(:)
	character*10, allocatable :: pm(:,:) ! (0:1,-1:1) ! HTML/TeX , 1=plus, -1=minus, 0=empty

	integer, dimension(0:1, 0:1) :: operatorPositionForQ ! first/second operator, first/second version
	logical, dimension(0:1) :: v_q_double
	
	logical prepared
	logical hermitian
	logical master ! or slave when it is H.c. of his master
	type( Perturbation ), pointer :: hc	
	
	integer numAllOperators

	! for Kroneckers:
	integer, allocatable :: deltaPositions(:)
		! 12 34 with delta12, delta34 -> 2, 4
		! 123 456 with delta12/23/13, delta45/46/56 -> 2,3, 5,6
	integer, allocatable :: deltaNumCounterparts(:)
		! 12 34 with delta12, delta34 -> 1, 3 -> counts=1,1
		! 123 456 with delta12/23/13, delta45/46/56 -> 1,2, 3,4 -> counts=1,2,1,2
	integer numDeltaPositions
		! 12 34 -> 2,4 -> 2 ; 123 456 -> 2,3,5,6 -> 4
		! i.e. sum( numCre-1 + numAnn-1 )

	integer numDeltaCombinations
		! 12 3 -> 2 ; 123 4 -> 2*3 = 6
		! 12 34 -> 2*2=4 ; 123 456 -> 2*3*2*3 = 36
		! i.e. prod( numCre! * numAnn! )

	integer maxNumEqualMomenta
	! delta_12 delta_13 delta_14 -> up to 4 equal momenta ; delta12 -> up to 2
	

contains
    procedure :: addSubsystem => Perturbation_addSubsystem
    procedure :: setAmplitudeConst => Perturbation_setAmplitudeConst
    procedure :: setAmplitudeTable => Perturbation_setAmplitudeTable

    procedure :: prepare => Perturbation_prepare
    procedure :: dropPreparations => Perturbation_dropPreparations
	
	procedure :: show => Perturbation_show

	procedure :: showQbeTermTitle => Perturbation_showQbeTermTitle
	procedure :: showSumHead => Perturbation_showSumHead
	procedure :: showSumTail => Perturbation_showSumTail

    procedure :: toString => Perturbation_toString
    
    procedure :: containsSubsystem => Perturbation_containsSub
    procedure :: indexOf => Perturbation_indexOf

	final :: Perturbation_destructor

end type Perturbation

type Perturbation_p
	class(Perturbation), pointer :: obj
end type Perturbation_p



interface Perturbation
	module procedure Perturbation_init
end interface
interface InteractionPairwise
	module procedure InteractionPairwise12_init
	module procedure InteractionPairwise1_init
end interface
interface InteractionWithPhonons
	module procedure InteractionWithPhonons_init
end interface
interface OpticalPumpingBose
	module procedure OpticalPumpingBose_init
end interface
interface OpticalPumpingFermi
	module procedure OpticalPumpingFermi_init 
end interface

contains

! ------------ constructors -------

	function Perturbation_init( title ) result(this)
    character(len=*) title
	type(Perturbation), pointer :: this

		allocate( this )

		this % title = title
		this % numSubsystems = 0
		this % currentSubsystem = 0
		
		this % mode = Q_INDEPENDENT
		this % amplitudeValue = 1.d0
		this % numOperators(:,:) = 0 ! ready for adding
		this % operatorSpinUp(:) = .TRUE. ! ready for adding
    	this % numAllOperators = 0
		
		allocate( this % subsystems(0:max_num_subsystems-1) )
		
		allocate( this % pm(0:1, -1:1) ) ! text lines for debug output
		
		this % hermitian = .TRUE.  ! but actually unknown
		this % prepared = .FALSE.
		this % master = .TRUE.
		
    end function Perturbation_init

!--

    function InteractionPairwise12_init( subsystem1, subsystem2, sameSpin ) result(this)
    type(Subsystem), pointer :: subsystem1, subsystem2
	logical, optional :: sameSpin
	logical secondSpinUp
    
	type(Perturbation), pointer :: this

		secondSpinUp = .FALSE.
		if (present(sameSpin)) then
			secondSpinUp = sameSpin
		end if

		this => Perturbation_init( 'e-e' )
		if (secondSpinUp) then
	    	call this % addSubsystem( subsystem1, 2, 2 )
		else
			call this % addSubsystem( subsystem1, 1, 1,  .TRUE. )
		   	call this % addSubsystem( subsystem2, 1, 1,  .FALSE. )
		end if


    end function InteractionPairwise12_init

!--

	function InteractionPairwise1_init( subsystem1 ) result(this)
    type(Subsystem), pointer :: subsystem1
    
	type(Perturbation), pointer :: this

		this => Perturbation_init( 'pair' )
    	call this % addSubsystem( subsystem1, 2, 2 )

    end function InteractionPairwise1_init

!--

    function InteractionWithPhonons_init( subsystem1, subsystem2 ) result(this)
    type(Subsystem), pointer :: subsystem1, subsystem2
    
	type(Perturbation), pointer :: this

		if (.NOT.subsystem2 % boseStatistics) stop 'InteractionWithPhonons : subsystem2 not bosons'

    	this => Perturbation_init( 'e-ph' )
    	call this % addSubsystem( subsystem1, 1, 1 )
!		call this % addSubsystem( subsystem2, 1, 0 )
		call this % addSubsystem( subsystem2, 0, 1 )
		
    end function InteractionWithPhonons_init

!--

	function OpticalPumpingBose_init( subsystem1 ) result(this)
	type(Subsystem), pointer :: subsystem1
    
	type(Perturbation), pointer :: this

		stop 'OpticalPumpingBose : not ready (nonzero k?)'
		if (.NOT.subsystem1 % boseStatistics) stop 'OpticalPumpingBose : subsystem1 not bosons'

    	this => Perturbation_init( 'OpticalPumpingBose' )
    	call this % addSubsystem( subsystem1, 1, 0 )
    	
    end function OpticalPumpingBose_init

	function OpticalPumpingFermi_init( subsystem1 ) result(this)
	type(Subsystem), pointer :: subsystem1
    
	type(Perturbation), pointer :: this

		stop 'OpticalPumpingFermi : not ready (h omega...)'
		if (subsystem1 % boseStatistics) stop 'OpticalPumpingFermi : subsystem1 not fermions'

    	this => Perturbation_init( 'OpticalPumpingFermi' )
    	call this % addSubsystem( subsystem1, 2, 0 )
    	
    end function OpticalPumpingFermi_init

	
! ----------------

	subroutine Perturbation_destructor( this )
	type(Perturbation) :: this
	
!		print*, ' Perturbation: destructor'
		call this % dropPreparations
		deallocate( this % pm )
		deallocate( this % subsystems )
		
	end subroutine Perturbation_destructor	
	
	
! ------------ setters ------------------

	subroutine Perturbation_addSubsystem( this, s, numCreation, numAnnihilation , operatorSpinUp )
    class( Perturbation ), intent(INOUT) :: this
	type( Subsystem ), pointer :: s
	integer numCreation, numAnnihilation
	logical, optional :: operatorSpinUp
	
	integer pos
	
		pos = this % numSubsystems
		if (pos.GE.max_num_subsystems) then
			stop 'exceeded max_num_subsystems'
		end if

    	this % numOperators( CREATION, pos ) = numCreation
    	this % numOperators( ANNIHILATION, pos ) = numAnnihilation
    	if (present(operatorSpinUp)) this % operatorSpinUp( pos ) = operatorSpinUp ! or TRUE from constructor
		this % subsystems(pos) % obj => s

    	this % numAllOperators = this % numAllOperators + numCreation + numAnnihilation

		this % numSubsystems = pos + 1
		call this % dropPreparations
				
    end subroutine Perturbation_addSubsystem


    subroutine Perturbation_dropPreparations( this )
    class( Perturbation ), intent(INOUT) :: this

	integer i

		if (this % prepared) then
			deallocate( this % deltaNumCounterparts )
			deallocate( this % deltaPositions )
			do i=this % numAllOperators-1, 0, -1
				deallocate( this % info(i) % subscript )
				deallocate( this % info(i) % superscript )
			end do
			deallocate( this % info )
		end if
		this % prepared = .FALSE.

    end subroutine Perturbation_dropPreparations

!--

    subroutine Perturbation_setAmplitudeConst( this, v )
    class( Perturbation ), intent(INOUT) :: this
	double precision v
	
		this % mode = Q_INDEPENDENT
		this % amplitudeValue = v
	
	end subroutine Perturbation_setAmplitudeConst

!--

    subroutine Perturbation_setAmplitudeTable( this, ToV, operators1forQ, operators2ForQ )
    class( Perturbation ), intent(INOUT) :: this
	type( TableOfValues ), pointer :: ToV
	integer :: operators1forQ(0:1) ! -1 if single operator
	integer, optional :: operators2forQ(0:1)
	
		this % mode = Q_DEPENDENT
		this % amplitudeValue = -666.d0
		this % amplitudeTable => ToV

! V(q) source for q : A(+) A(+) A A -> k2-k3 or k1-k4 ;   A(+) B(+) A => k2 or k1-k3
		! 2 versions of terms, to use one not used in first index of QBE term 
	
		this % operatorPositionForQ(:,0) = operators1forQ(:)
		this % V_q_double(0) = (operators1forQ(1).GE.0)
		if (present(operators2forQ)) then
			this % operatorPositionForQ(:,1) = operators2forQ(:)
			this % V_q_double(1) = (operators2forQ(1).GE.0)
		else
			this % operatorPositionForQ(:,1) = -1
		end if


	
	end subroutine Perturbation_setAmplitudeTable


! ------------ more methods -----------------

	subroutine Perturbation_show(this, o )
    class( Perturbation ) :: this
    class( Output ), pointer :: o
    
    integer i, pos, sId
    character a
    character*200 subscript
    character*20 spin_subscript
    character*10 ahat, qbold

		call this % prepare

		if ( this%mode .EQ. Q_INDEPENDENT ) then
			call o % formulaItem( 'V' )
			call o % upperLowerIndex( ' ', this % title, .FALSE., .TRUE. )
		end if

		call o % summation()
		call o % beginLimits()

		subscript = ''
		do i=0, this%numAllOperators-2
			subscript = TRIM(subscript) // TRIM( this % info(i) % subscript( o % version ) )
		end do
		call o % subscript( subscript )

		call o % endLimits()

		if (.NOT.this % hermitian) then
			call o % openBracket( BRACKET_SQUARE )
		end if
		if ( this%mode .EQ. Q_DEPENDENT ) then
			call o % formulaItem( 'V' )
			call o % upperLowerIndex( ' ', this % title, .FALSE., .TRUE. )
			call o%bold( 'q', qbold)
			call o % formulaItem( '(' // qbold // '=' )
			do i=0, 1
				pos = this % operatorPositionForQ(0,i)
				if ( pos.GE.0 ) then
					if (i.GT.0) call o % formulaItem( ',' )
					call o % formulaItem( this % info(pos) % subscript(o%version) )
					pos = this % operatorPositionForQ(1,i)
					if ( pos.GE.0 ) then
						call o % formulaItem( TRIM( o%minus ) // this % info(pos) % subscript(o%version) )
					end if
				end if
			end do
			call o % formulaItem( ')' )
		end if
		do i=0, this%numAllOperators-1
			a = this % info(i) % a
			call o % hat( a, ahat )
			call o % formulaItem( ahat )
			if (this % info(i) % dagger) then
				call o % superscript( o % dagger )
			else
				call o % superscript( ' ' )
			end if

			
			subscript = this % info(i) % subscript( o % version )
			if (i .EQ. this%numAllOperators-1 ) then
				subscript = TRIM(subscript) // '=' // TRIM(this % lastOperatorSubscript( o % version ))
			end if

			spin_subscript = ''
			if (.NOT.this % info(i) % bose) then
				sId = this % info(i) % subsystemIndex
				if (this%operatorSpinUp(sId)) then
					spin_subscript = o % up_arrow
				else
					spin_subscript = o % down_arrow
				end if
				if (i .EQ. this%numAllOperators-1 ) then
					spin_subscript = ',' // TRIM(spin_subscript)
				end if
			end if

			call o % subscript( TRIM(subscript) // TRIM(spin_subscript) )
			
		end do

		if (.NOT.this % hermitian) then
			call o % FormulaItem( '+ H.c.' )
			call o % closeBracket
		end if

    
	end subroutine Perturbation_show


!------- show QBE formulas

	subroutine Perturbation_showQbeTermTitle(this, s_index, o )
    class( Perturbation ) :: this
	integer s_index
    class( Output ), pointer :: o
    
	type( Subsystem ), pointer :: s
	character a, k
	character*20 acomma, kbold
    integer i
	
		call this % prepare

		s => this % subsystems(s_index) % obj
		a = s % operatorLetter
		acomma = ',' // a
		if (subsystem_count.EQ.1) acomma=''

		i = this % firstOperatorEntry(s_index)
		call o % formulaItem( 'J' )
		call o % upperLowerIndex( '(' // &
&			TRIM(this % title) // TRIM(acomma) // ')' , &
&			this % info(i) % subscript( o%version )	 )

	end subroutine Perturbation_showQbeTermTitle


	subroutine Perturbation_showSumHead( this, o, skipIndex )
	class( Perturbation ) :: this
	class( Output ), pointer :: o
	logical :: skipIndex(0:this%numAllOperators-1)

	integer i
	character*200 s

		call o % summation()

		call o % beginLimits()

		s = ''
		do i=0, this % numAllOperators-1
			if (.NOT.skipIndex(i)) then
				s = TRIM(s) //  TRIM(this % info(i) % subscript( o % version ))
			end if
		end do

		call o % subscript( s )

		call o % endLimits()
		call o % openBracket( BRACKET_CURLY )

	end subroutine Perturbation_showSumHead

	subroutine Perturbation_showSumTail( this, o )
	class( Perturbation ) :: this
	class( Output ), pointer :: o

		call o % closeBracket
		
	end subroutine Perturbation_showSumTail



!-------

	character*80 function Perturbation_toString(this)
    class( Perturbation ) :: this

	character*10 :: n1
	character*20 ::ploc

		write(ploc, '(I0)') loc(this)
		write(n1, '(I5)') this % numSubsystems

		Perturbation_toString = 'Perturbation@' // TRIM(ploc) // ': ' &
&		 // TRIM(this % title) // ' num. subsystems: ' // TRIM(n1)

	end function Perturbation_toString

!--

	logical function Perturbation_containsSub( this, s )
    class( Perturbation ) :: this
	type( Subsystem ), pointer :: s, s1

	integer i
	logical ans
		ans = .FALSE.
		do i=0, this % numSubsystems - 1
			s1 => this % subsystems(i) % obj
			if (associated(s1 , s )) ans = .TRUE.
		end do
		Perturbation_containsSub = ans
	end function Perturbation_containsSub

	integer function Perturbation_indexOf( this, s )
    class( Perturbation ) :: this
	type( Subsystem ), pointer :: s
	
	type( Subsystem ), pointer :: s1
	integer i, ans
	
		ans = -1
		do i=0, this % numSubsystems - 1
			if (ans.LT.0) then
				s1 => this % subsystems(i) % obj
				if (associated(s1 , s )) ans = i
			end if
		end do
		Perturbation_indexOf = ans

	end function Perturbation_indexOf

!--

	subroutine Perturbation_prepare(this)
	class( Perturbation ), intent(INOUT) :: this
	integer i, pos, numCre, numAnn, j, sId, pm, npk, numAllOp
	type(Subsystem), pointer :: s
	logical lastDagger, minus
	character a, k
	character*2 s1, s2
	character*300 texSubs, htmlSubs, sigma_k

	  if ( .NOT. this % prepared ) then
		numAllOp = this%numAllOperators

		allocate( this % firstOperatorEntry( 0:this%numSubsystems-1) )
		allocate( this %  lastOperatorEntry( 0:this%numSubsystems-1) )
		this % firstOperatorEntry(:) = -1
		this %  lastOperatorEntry(:) = -1
		
		allocate( this % info(0:numAllOp-1) )

		do i=0, this % numAllOperators-1
			allocate( this % info(i) % superscript(0:1) )
			allocate( this % info(i) % subscript(0:1) )
			this % info(i) % dagger = .FALSE.
			this % info(i) % free = ( i .LT. numAllOp-1 ) ! last index is not free : sum=0
			this % info(i) % firstOperatorForSubsystem = .FALSE.
		end do

		this % numDeltaCombinations = 1 ! including one without delta : (1+delta)*(1+delta) -> 4 terms
		this % numDeltaPositions = 0
		allocate( this % deltaPositions(0: numAllOp-1 ) )
		allocate( this % deltaNumCounterparts(0: numAllOp-1 ) )
		this % deltaPositions(:) = -1
		this % deltaNumCounterparts(:) = 0

		this % hermitian = .TRUE.

! ordering		
		pos = 0
		do i=0, this % numSubsystems-1
			s => this % subsystems(i) % obj
			numCre = this % numOperators( CREATION, i )
			numAnn = this % numOperators( ANNIHILATION, i )
			if (numCre.NE.numAnn) this % hermitian = .FALSE.

			do j=0, numCre-1
				this % info(pos) % dagger = .TRUE. ! creation
				this % info(pos) % subsystemIndex = i
				this % info(pos) % subsystem => s
				this % info(pos) % indexPos = j
				npk = j
				this % info(pos) % numPossibleKroneckers = npk
				if (npk.GT.0) then
					this % numDeltaCombinations = this % numDeltaCombinations * (npk+1)
					this % deltaPositions( this % numDeltaPositions ) = pos
					this % deltaNumCounterparts( this % numDeltaPositions ) = npk+1
					this % numDeltaPositions = this % numDeltaPositions + 1
				end if
				pos = pos + 1
			end do
		end do

		do i=this % numSubsystems-1, 0, -1
			s => this % subsystems(i) % obj
			numCre = this % numOperators( CREATION, i )
			numAnn = this % numOperators( ANNIHILATION, i )
			do j=numCre, numCre+numAnn-1
				this % info(pos) % dagger = .FALSE. ! annihilation
				this % info(pos) % subsystemIndex = i
				this % info(pos) % subsystem => s
				this % info(pos) % indexPos = j
				npk = j - numCre
				this % info(pos) % numPossibleKroneckers = npk
				if (npk.GT.0) then
					this % numDeltaCombinations = this % numDeltaCombinations * (npk+1)
					this % deltaPositions( this % numDeltaPositions ) = pos
					this % deltaNumCounterparts( this % numDeltaPositions ) = npk+1
					this % numDeltaPositions = this % numDeltaPositions + 1
				end if
				pos = pos + 1
			end do
		end do

		do j=0, numAllOp-1
			i = this % info(j) % subsystemIndex
			if (this % firstOperatorEntry(i).LT.0) then
				this % firstOperatorEntry(i) = j
				this % info (i) % firstOperatorForSubsystem = .TRUE.
			end if
		end do


! description texts		
		lastDagger = this % info(this % numAllOperators-1) % dagger
		texSubs = ''
		htmlSubs = ''
		do i=0, numAllOp-1
			sId = this % info(i) % subsystemIndex
			s => this % subsystems(sId) % obj
			a = s % operatorLetter
			k = indexLetters(sId+1:sId+1)
			
			write(s1, '(I0)')  this % info(i) % indexPos + 1

			this % info(i) % a = a
			this % info(i) % k = k
			this % info(i) % bose = s % boseStatistics
			this % info(i) % superscript(HTML) = '<sup>(' // a // ')</sup>'
			this % info(i) % superscript(TeX) = '\mathrm{(' // a // ')}'
			if (subsystem_count.EQ.1) then
				this % info(i) % superscript(HTML) = ''
				this % info(i) % superscript(TeX) = ''
			end if
			this % info(i) % subscript(HTML) = '<b>' // k // '</b><sub>' // TRIM(s1) // '</sub>'
			this % info(i) % subscript(TeX) = '{{\bf ' // k // '}_{' // TRIM(s1) // '}}'
			
			if (this % info(i) % free) then
				minus = (this%info(i)%dagger .EQV. lastDagger) 
				if (minus .OR. (i.GT.0)) then
					if (minus) then
						texSubs  = TRIM(texSubs)  // '-'
						htmlSubs = TRIM(htmlSubs) // HTML_MINUS
					else
						texSubs  = TRIM(texSubs)  // '+'
						htmlSubs = TRIM(htmlSubs) // '+'
					end if
				end if
				texSubs = TRIM(texSubs) // this % info(i) % subscript(TeX)
				htmlSubs = TRIM(htmlSubs) // this % info(i) % subscript(HTML)
			end if

			minus = (this%info(i)%dagger .NEQV. this%info(0)%dagger) 
			if (minus) then
				pm  = -1
			else
				pm  = +1
			end if
			this % info(i) % pm  = pm

		end do		
		this % pm(TeX,  +1) = '+'
		this % pm(HTML, +1) = '+'
		this % pm(TeX,  -1) = '-'
		this % pm(HTML, -1) = HTML_MINUS
		this % pm(TeX,   0) = 'EMPTY'
		this % pm(HTML,  0) = '&nbsp;'
		allocate( this % lastOperatorSubscript(0:1) )
		this % lastOperatorSubscript(TeX)  = texSubs
		this % lastOperatorSubscript(HTML) = htmlSubs

		this % maxNumEqualMomenta = maxval( this % deltaNumCounterparts( 0:this % numDeltaPositions-1) )
		this % prepared = .TRUE.

		if (.NOT.this%hermitian .AND. this%master) then
			this % hc => Perturbation( 'H.c.(' // TRIM(this%title) // ')' )
			this % hc % master = .FALSE.
			do i=0,this%numSubsystems-1
				s => this % subsystems(i) % obj
				numCre = this % numOperators( CREATION, i )
				numAnn = this % numOperators( ANNIHILATION, i )
				call this % hc % addSubsystem( s, numAnn, numCre )
			end do
			if (this % mode.EQ. Q_DEPENDENT) then
				call this % hc % setAmplitudeTable( this % amplitudeTable, &
&					this % operatorPositionForQ(:,0), this % operatorPositionForQ(:,1) )
			else
				call this % hc % setAmplitudeConst( this % amplitudeValue )
			end if
			call this % hc % prepare
		end if

	  end if



	end subroutine Perturbation_prepare



end module Perturbation_module

