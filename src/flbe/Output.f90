module Output_module

	implicit none
	
	type textIndex_p
		character(len=:), allocatable :: s
		logical text
	end type
	integer, parameter :: max_num_subscripts = 200 ! or something else

! -- OutputStream
	
	integer, parameter :: max_bracket_level=20
	
! output versions
	integer, parameter :: HTML = 0 ! output versions can be used as array index
	integer, parameter :: TeX = 1

! bracket versions
	integer, parameter :: BRACKET_INVISIBLE = -1
	integer, parameter :: BRACKET_ROUND		= 0 ! ()
	integer, parameter :: BRACKET_SQUARE	= 1 ! []
	integer, parameter :: BRACKET_CURLY		= 2 ! {}
	integer, parameter :: BRACKET_ANGLE		= 3 ! <>
	
	type, abstract :: Output
		integer version ! 0 or 1
		integer f ! file
		integer superscriptTextLevel ! how many times switched on : >0 means YES
		integer subscriptTextLevel
		integer, allocatable :: bracketValues(:)
		integer bracketLevel
		integer bracketNextVersion
		logical :: deleteWhenClose

	! text constants
		character*20 :: dagger ! 20 due to longer HTML code
		character*10 :: times
		character*10 :: minus
		character*10 :: gamma
		character*10 :: sigma
		character*10 :: tau
		character*20 :: varepsilon
		character*10 :: pi
		character*10 :: hbar
		character*20 :: twopi_hbar
		character*10 :: deltasmall
		character*20 :: right_arrow
		character*20 :: left_arrow
		character*20 :: up_arrow
		character*20 :: down_arrow
		character*10 :: empty_space

	contains
	
! level-0 :
		procedure :: markToDelete => Output_proc_markToDelete

		procedure(OutputFunction1s1s), deferred :: bold
		procedure(OutputFunction1s1s), deferred :: hat

! level-1 : write command
		procedure :: write => Output_proc_writeS123
		procedure :: writeln => Output_proc_writelnS123
		procedure :: plusminus => Output_proc_plusminus
		procedure :: writeInteger => Output_proc_writeInteger
		procedure :: writeIntegerExcludingPM1 => Output_proc_writeIntegerExcludingPM1
		procedure :: writeIntegerExcludingPM1WithPlusSign => Output_proc_writeIntegerExcludingPM1WithPlusSign
		
		procedure(OutputSubroutine0), deferred :: newLine
		procedure(OutputSubroutine0), deferred :: horizontalLine
		procedure(OutputSubroutine1), deferred :: newSection
		
		procedure :: openBracket  => Output_proc_openBracket
		procedure :: closeBracket => Output_proc_closeBracket
		procedure(OutputSubroutineShowBracket), deferred :: showBracket
		procedure :: closeBrackets  => Output_proc_closeBrackets
		procedure :: openPseudoBrackets  => Output_proc_openPseudoBrackets
		procedure :: closePseudoBrackets => Output_proc_closePseudoBrackets

		procedure(OutputSubroutine0), deferred :: beginEquation
		procedure(OutputSubroutine0_INOUT), deferred :: endEquation

		procedure(OutputSubroutine0), deferred :: summation

		procedure(OutputSubroutine0), deferred :: beginLimits
		procedure(OutputSubroutine0), deferred :: endLimits

		procedure(OutputSubroutine0), deferred :: beginFormulaItem
		procedure(OutputSubroutine0), deferred :: endFormulaItem

		procedure(OutputSubroutine0), deferred :: beginGreater
		procedure(OutputSubroutine0), deferred :: endGreater

		procedure(OutputSubroutine0), deferred :: beginText
		procedure(OutputSubroutine0), deferred :: endText

		procedure(OutputSubroutine0), deferred :: beginSuperscript
		procedure(OutputSubroutine0), deferred :: endSuperscript

		procedure :: switchSuperscriptText => Output_proc_switchSuperscriptText
		procedure :: switchSubscriptText => Output_proc_switchSubscriptText

		procedure(OutputSubroutine0), deferred :: beginSubscript
		procedure(OutputSubroutine0), deferred :: endSubscript

		procedure(OutputSubroutine0), deferred :: beginFrac
		procedure(OutputSubroutine0), deferred :: middleFrac
		procedure(OutputSubroutine0), deferred :: endFrac

! level-2 : combine several commands
		procedure :: formulaItem => Output_proc_formulaItem
		procedure :: formulaItemGreater => Output_proc_formulaItemGreater
		procedure :: formulaItemText => Output_proc_formulaItemText
		procedure :: upperLowerIndex => Output_proc_upperLowerIndex

!		procedure :: test => Output_proc_test
		
		procedure :: superscript => Output_proc_superscript
		procedure :: superscriptText => Output_proc_superscriptText
		procedure(OutputSubroutine1_INOUT), deferred :: subscript
		
		procedure(OutputSubroutine0_INOUT), deferred :: equationNewLine

	end type Output



	type, extends(Output) :: OutputTeX
	contains
! level-0 : modify text
		procedure :: bold => TeX_proc_bold
		procedure :: hat => TeX_proc_hat

! level-1 : write command
		procedure :: showBracket => TeX_proc_showBracket

		procedure :: newLine => TeX_proc_newLine
		procedure :: horizontalLine => TeX_proc_horizontalLine
		procedure :: newSection => TeX_proc_newSection

		procedure :: beginEquation => TeX_proc_beginEquation
		procedure :: endEquation => TeX_proc_endEquation

		procedure :: summation => TeX_proc_summation

		procedure :: beginLimits => TeX_do_nothing
		procedure :: endLimits => TeX_do_nothing

		procedure :: beginFormulaItem => TeX_do_nothing
		procedure :: endFormulaItem => TeX_do_nothing

		procedure :: beginGreater => TeX_do_nothing
		procedure :: endGreater => TeX_do_nothing

		procedure :: beginText => TeX_proc_beginText
		procedure :: endText => TeX_proc_endText

		procedure :: beginSuperscript => TeX_proc_beginSuperscript
		procedure :: endSuperscript => TeX_proc_endSuperscript

		procedure :: beginSubscript => TeX_proc_beginSubscript
		procedure :: endSubscript => TeX_proc_endSubscript

		procedure :: beginFrac => TeX_proc_beginFrac
		procedure :: middleFrac => TeX_proc_middleFrac
		procedure :: endFrac => TeX_proc_endFrac

! level-2 : combine several commands
		procedure :: subscript => TeX_proc_subscript
		procedure :: equationNewLine => TeX_proc_equationNewLine

		final :: TeX_proc_endDocument
	end type OutputTeX

	type, extends(Output) :: OutputHTML
		integer numSubscripts
		type(textIndex_p), allocatable :: subscriptStorage(:)
	contains
		procedure :: dropSubscripts => HTML_proc_dropSubscripts
		procedure :: showLateSubscripts => HTML_proc_showLateSubscripts

! level-0 : modify text
		procedure :: bold => HTML_proc_bold
		procedure :: hat => HTML_proc_hat

! level-1 : write command
		procedure :: showBracket => HTML_proc_showBracket

		procedure :: newLine => HTML_proc_newLine
		procedure :: horizontalLine => HTML_proc_horizontalLine
		procedure :: newSection => HTML_proc_newSection

		procedure :: beginEquation => HTML_proc_beginEquation
		procedure :: endEquation => HTML_proc_endEquation

		procedure :: summation => HTML_proc_summation

		procedure :: beginLimits => HTML_do_nothing
		procedure :: endLimits => HTML_do_nothing

		procedure :: beginFormulaItem => HTML_proc_beginFormulaItem
		procedure :: endFormulaItem => HTML_proc_endFormulaItem

		procedure :: beginGreater => HTML_proc_beginGreater
		procedure :: endGreater => HTML_proc_endGreater

		procedure :: beginText => HTML_do_nothing
		procedure :: endText => HTML_do_nothing

		procedure :: beginSuperscript => HTML_proc_beginSuperscript
		procedure :: endSuperscript => HTML_proc_endSuperscript

		procedure :: beginSubscript => HTML_proc_beginSubscript
		procedure :: endSubscript => HTML_proc_endSubscript

		procedure :: beginFrac => HTML_proc_beginFrac
		procedure :: middleFrac => HTML_proc_middleFrac
		procedure :: endFrac => HTML_proc_endFrac

! level-2 : combine several commands
		procedure :: subscript => HTML_proc_subscript
		procedure :: equationNewLine => HTML_proc_equationNewLine
		
		final :: HTML_proc_endDocument
	end type OutputHTML

! polymorphic subroutine
	abstract interface
		subroutine OutputSubroutine0(this) ! 0 arguments
			import Output
			class(Output), intent(IN) :: this
		end subroutine
		subroutine OutputSubroutine0_INOUT(this) ! 0 arguments, this INOUT
			import Output
			class(Output), intent(INOUT) :: this
		end subroutine
		subroutine OutputSubroutine1(this, s) ! 1 string in
			import Output
			class(Output), intent(IN) :: this
			character(len=*) :: s
		end subroutine
		subroutine OutputSubroutine1I(this, n) ! 1 integer in
			import Output
			class(Output), intent(IN) :: this
			integer n
		end subroutine
		subroutine OutputSubroutine1_INOUT(this, s) ! 1 string in, this INOUT
			import Output
			class(Output), intent(INOUT) :: this
			character(len=*) :: s
		end subroutine
		subroutine OutputFunction0s1s(this, s_out) ! 1 string out
			import Output
			class(Output), intent(IN) :: this
			character(len=*), intent(OUT) :: s_out
		end subroutine
		subroutine OutputFunction1s1s(this, s_in, s_out) ! 1 string is, 1 string out
			import Output
			class(Output), intent(IN) :: this
			character(len=*), intent(IN) :: s_in
			character(len=*), intent(OUT) :: s_out
		end subroutine
		subroutine OutputSubroutineL(this, b) ! logical
			import Output
			class(Output), intent(IN) :: this
			logical b
		end subroutine
		subroutine OutputSubroutineShowBracket(this, version, opening) ! integer, logical
			import Output
			class(Output), intent(IN) :: this
			integer version
			logical opening
		end subroutine
	end interface

	interface OutputTeX
		module procedure TeX_beginDocument
	end interface OutputTeX
	interface OutputHTML
		module procedure HTML_beginDocument
	end interface OutputHTML


! HTML symbols

	character*4, parameter :: HTML_NEW_LINE = '<br>'

	character*31, parameter :: HTML_TABLE_HEAD =  '<table border=0 class=a><tbody>123'
	character*16, parameter :: HTML_TABLE_TAIL =  '</tbody></table>'

	character*7, parameter :: HTML_MINUS = '&minus;'
	character*7, parameter :: HTML_SUM = '&#8721;'
	character*6, parameter :: HTML_HAT = '&#770;'
	character*7, parameter :: HTML_LEFT_LT   = '&#9001;'
	character*7, parameter :: HTML_RIGHT_GT = '&#9002;'
	
	character*7, parameter :: HTML_SIGMA = '&sigma;'
	character*8, parameter :: TeX_SIGMA = '{\sigma}'

	
! TeX symbols

	character*12, parameter :: TeX_NEW_LINE = '\\ \nonumber'

!-------------
	
contains

!--------------------------------

	include 'OutputStream.f90'

!--------------------------------


end module Output_module
