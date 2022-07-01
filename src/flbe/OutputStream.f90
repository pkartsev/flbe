
! constructors : filename + file number

	function HTML_beginDocument( filename, fn ) result(this)
	class(OutputHTML), pointer :: this
	character(len=*) :: filename
	integer, optional :: fn ! file number

	integer f
		allocate(this)
		this % version = HTML

		this % superscriptTextLevel = 0
		this % subscriptTextLevel = 0

		allocate( this % bracketValues( 0:max_bracket_level-1 ) )
		this % bracketLevel = 0
		this % bracketNextVersion = BRACKET_ROUND
		
		this % numSubscripts = 0
		allocate( this % subscriptStorage( 0:max_num_subscripts-1) )
		
		this % dagger = '<sup>&dagger;</sup>'
		this % minus  = '&minus;'
		this % times  = '&times;'
		this % gamma  = '&gamma;'
		this % sigma  = '&sigma;'
		this % tau    = '&tau;'
		this % varepsilon    = '&epsilon;'
		this % deltasmall = '&#948;'
		this % pi  = '&pi;'
		this % hbar = '&#x210F;'
		this % twopi_hbar = '(2&pi;/&#x210F;)'
		this % right_arrow = '&rarr;'
		this % left_arrow = '&larr;'
		this % up_arrow = '&uarr;'
		this % down_arrow = '&darr;'
		this % empty_space= '&nbsp;'

		f = 12
		if (present(fn)) f = fn
	
		open(f, file=filename)
		this % f = f
		this % deleteWhenClose = .FALSE.
		
		write( this % f, '(A)') '<!doctype html><html><head><meta charset="utf-8"><style>'
		write( this % f, '(A)') 'body{font-size: 18pt} td{font-style:italic} .a{font-size: 18pt}'
		write( this % f, '(A)') '.greater {font-size: 32pt} .smaller {font-size: 14pt}'
		write( this % f, '(A)') 'td.normal {font-style: normal}</style></head><body>'
		
	end function HTML_beginDocument
	
	function TeX_beginDocument( filename, fn ) result(this)
	class(OutputTeX), pointer :: this
	character(len=*) :: filename
	integer, optional :: fn ! file number

	integer f
	
		allocate(this)
		this % version = TeX

		this % superscriptTextLevel = 0
		this % subscriptTextLevel = 0

		allocate( this % bracketValues( 0:max_bracket_level-1 ) )
		this % bracketLevel = 0
		this % bracketNextVersion = BRACKET_ROUND

		this % dagger = '{\dagger}'
		this % times  = '{\times}'
		this % minus = '-'
		this % gamma = '{\gamma}'
		this % sigma = '{\sigma}'
		this % tau = '{\tau}'
		this % varepsilon    = '{\varepsilon}'
		this % deltasmall = '{\delta}'
		this % pi  = '{\pi}'
		this % hbar = '{\hbar}'
		this % twopi_hbar = '\frac{2 \pi}{\hbar}'
		this % right_arrow = '{\rightarrow}'
		this % left_arrow = '{\leftarrow}'
		this % up_arrow = '{\uparrow}'
		this % down_arrow = '{\downarrow}'
		this % empty_space= ' '

		f = 13
		if (present(fn)) f = fn
	
		open(f, file=filename)
		this % f = f
		this % deleteWhenClose = .FALSE.
		
		write( this % f, '(A)') '\documentclass[12pt]{article}'
		write( this % f, '(A)') '\begin{document}'
		
	end function TeX_beginDocument

	subroutine Output_proc_markToDelete( this, v )
	class(Output), intent(INOUT) :: this	
	logical, optional :: v

		if (present(v)) then
			this % deleteWhenClose = v
		else
			this % deleteWhenClose = .TRUE.
		end if

	end subroutine Output_proc_markToDelete

! destructors

	subroutine HTML_proc_endDocument( this )
	type(OutputHTML) :: this
	
		if (this % deleteWhenClose) then
			close( this % f, status='delete' )
		else
			write( this % f, '(A)') '</body></html>'
			close( this % f )
		end if
		deallocate( this % subscriptStorage )
		deallocate( this % bracketValues )
		
	end subroutine HTML_proc_endDocument

	subroutine TeX_proc_endDocument( this )
	type(OutputTeX) :: this
	
		if (this % deleteWhenClose) then
			close( this % f, status='delete' )
		else
			write( this % f, '(A)') '\end{document}'
			close( this % f )
		end if
		deallocate( this % bracketValues )
		
	end subroutine TeX_proc_endDocument

!--------------------------------------------------

	subroutine Output_proc_writeS123( this, s1, s2, s3 )
	class(Output), intent(IN) :: this	
	character(LEN=*), optional :: s1
	character(LEN=*), optional :: s2, s3
	
		write( this % f, '(A)') TRIM(s1)
		if (present(s2)) write( this % f, '(A)') TRIM(s2)
		if (present(s3)) write( this % f, '(A)') TRIM(s3)
		
	end subroutine Output_proc_writeS123
	
	subroutine Output_proc_writelnS123( this, s1, s2, s3 )
	class(Output), intent(IN) :: this	
	character(LEN=*), optional :: s1
	character(LEN=*), optional :: s2, s3
	
		write( this % f, '(A)') TRIM(s1)
		if (present(s2)) write( this % f, '(A)') TRIM(s2)
		if (present(s3)) write( this % f, '(A)') TRIM(s3)
		write( this % f, *) ''
!		call this % newLine
				
	end subroutine Output_proc_writelnS123
	
	subroutine Output_proc_writeSIS( this, s1, i, s2 )
	class(Output), intent(IN) :: this	
	character(LEN=*) s1, s2
	integer i

		write( this % f, '(A)',  advance='no') s1
		write( this % f, '(I0)', advance='no') i
		write( this % f, '(A)',  advance='no') s2
		
	end subroutine Output_proc_writeSIS

	subroutine Output_proc_upperLowerIndex( this, super, sub, superText, subText )
	class(Output), intent(INOUT) :: this	
	character(LEN=*) super, sub
	logical, optional, intent(IN) :: superText, subText
	
		if (present(superText)) then
			if (superText) call this % switchSuperscriptText( .TRUE. )
		end if
		if (present(subText)) then
			if (subText) call this % switchSubscriptText( .TRUE. )
		end if
		
		call this % superscript( super )
		call this % subscript( sub )

		if (present(superText)) then
			if (superText) call this % switchSuperscriptText( .FALSE. )
		end if
		if (present(subText)) then
			if (subText) call this % switchSubscriptText( .FALSE. )
		end if
		
	end subroutine Output_proc_upperLowerIndex

!---------------------------------------------------------

	subroutine HTML_proc_dropSubscripts( this )
	class(OutputHTML), intent(INOUT) :: this
	integer i
	
		if (this%numSubscripts.GT.0) then
			do i=this%numSubscripts-1, 0, -1
				deallocate( this % subscriptStorage(i) % s )
			end do
		end if
		this%numSubscripts = 0
		
	end subroutine HTML_proc_dropSubscripts
	

!------------ level 2 --------------

	subroutine HTML_proc_subscript( this, s )
	class(OutputHTML), intent(INOUT) :: this
	character(LEN=*) s

		this % subscriptStorage( this % numSubscripts ) % s = s
		this % subscriptStorage( this % numSubscripts ) % text = &
&			(this % subscriptTextLevel .GT. 0)
		this % numSubscripts = this % numSubscripts + 1
		
	end subroutine HTML_proc_subscript
 
	subroutine TeX_proc_subscript( this, s )
	class(OutputTeX), intent(INOUT) :: this
	character(LEN=*) s
	
		call this % beginSubscript()
		write( this % f, '(A)', advance='no') TRIM(s)
		call this % endSubscript()
		
	end subroutine TeX_proc_subscript




	subroutine Output_proc_superscript( this, s )
	class(Output), intent(IN) :: this
	character(LEN=*) s
	
		call this % beginSuperscript()
		if (LEN(TRIM(s)).GT.0) then
			write( this % f, '(A)', advance='no') TRIM(s)
		else
			write( this % f, '(A)', advance='no') TRIM(this % empty_space)
		end if
		call this % endSuperscript()
		
	end subroutine Output_proc_superscript


	subroutine Output_proc_superscriptText( this, s )
	class(Output), intent(INOUT) :: this
	character(LEN=*) s
	
		call this % switchSuperscriptText( .TRUE. )
		call this % superscript( s )
		call this % switchSuperscriptText( .FALSE. )
		
	end subroutine Output_proc_superscriptText





	subroutine Output_proc_formulaItem( this, s )
	class(Output), intent(IN) :: this		
	character(LEN=*) s
	
		call this % beginFormulaItem()
		write( this % f, '(A)', advance='no') TRIM(s)
		call this % endFormulaItem()
		
	end subroutine Output_proc_formulaItem



	subroutine Output_proc_formulaItemGreater( this, s )
	class(Output), intent(IN) :: this	
	character(LEN=*) s
	
		call this % beginFormulaItem()
		call this % beginGreater()
		write( this % f, '(A)', advance='no') s
		call this % endGreater()
		call this % endFormulaItem()
		
	end subroutine Output_proc_formulaItemGreater

	subroutine Output_proc_formulaItemText( this, s )
	class(Output), intent(IN) :: this	
	character(LEN=*) s
	
		call this % beginFormulaItem()
		call this % beginText()
		write( this % f, '(A)', advance='no') s
		call this % endText()
		call this % endFormulaItem()
		
	end subroutine Output_proc_formulaItemText




!------------




	
	subroutine HTML_proc_bold( this, s_in, s_out )
	class(OutputHTML), intent(IN) :: this
	character(LEN=*), intent(IN) :: s_in
	character(LEN=*), intent(OUT) :: s_out

		s_out = '<b>' // TRIM(s_in) // '</b>'
	
	end subroutine HTML_proc_bold

	subroutine TeX_proc_bold( this, s_in, s_out )
	class(OutputTeX), intent(IN) :: this
	character(LEN=*), intent(IN) :: s_in
	character(LEN=*), intent(OUT) :: s_out

		s_out = '{\bf ' // TRIM(s_in) // '}'
	
	end subroutine TeX_proc_bold




	subroutine HTML_proc_hat( this, s_in, s_out )
	class(OutputHTML), intent(IN) :: this
	character(LEN=*), intent(IN) :: s_in
	character(LEN=*), intent(OUT) :: s_out

		s_out = s_in // HTML_HAT

	end subroutine HTML_proc_hat

	subroutine TeX_proc_hat( this, s_in, s_out )
	class(OutputTeX), intent(IN) :: this
	character(LEN=*), intent(IN) :: s_in
	character(LEN=*), intent(OUT) :: s_out

		s_out = '\hat{' // s_in // '}'

	end subroutine TeX_proc_hat

!------------ level 1 ---------------

	subroutine HTML_proc_NewLine( this )
	class(OutputHTML), intent(IN) :: this	
	
		write( this % f, *) HTML_NEW_LINE
		
	end subroutine HTML_proc_NewLine

	subroutine TeX_proc_NewLine( this )
	class(OutputTeX), intent(IN) :: this	
	
		write( this % f, *) TeX_NEW_LINE
		
	end subroutine TeX_proc_NewLine




	subroutine HTML_proc_horizontalLine( this )
	class(OutputHTML), intent(IN) :: this	
	
		write( this % f, '(A)') '<hr>'
		
	end subroutine HTML_proc_horizontalLine

	subroutine TeX_proc_horizontalLine( this )
	class(OutputTeX), intent(IN) :: this	
	
		write( this % f, *) ''
		write( this % f, '(A)') '\vspace{5pt}'
		write( this % f, '(A)') '\hrulefill'
		write( this % f, *) ''
	
	end subroutine TeX_proc_horizontalLine




	subroutine HTML_proc_newSection( this, s )
	class(OutputHTML), intent(IN) :: this	
	character(len=*) :: s	
	
		write( this % f, '(A)') '<h2>' // TRIM(s) // '</h2>'
		
	end subroutine HTML_proc_newSection

	subroutine TeX_proc_newSection( this, s )
	class(OutputTeX), intent(IN) :: this	
	character(len=*) :: s	
	
		write( this % f, '(A)') '\section{' // TRIM(s) // '}'
	
	end subroutine TeX_proc_newSection






	subroutine Output_proc_plusminus( this, b )
	class(Output), intent(IN) :: this
	logical b
	
		if (b) then
			write( this % f, '(A)') '+'
		else
			write( this % f, '(A)') this % minus
		end if
		
	end subroutine Output_proc_plusminus




	subroutine TeX_proc_equationNewLine( this )
	class(OutputTeX), intent(INOUT) :: this	
	
		call this % closePseudoBrackets
		call this % newLine()
		call this % openPseudoBrackets
		
	end subroutine TeX_proc_equationNewLine

	subroutine HTML_proc_equationNewLine( this )
	class(OutputHTML), intent(INOUT) :: this	
	
		call this % closePseudoBrackets
		call this % showLateSubscripts()
		write( this % f, '(A)') HTML_TABLE_TAIL
		write( this % f, '(A)') HTML_TABLE_HEAD
		call this % openPseudoBrackets
		
	end subroutine HTML_proc_equationNewLine



	subroutine Output_proc_writeInteger( this, n )
	class(Output), intent(IN) :: this	
	integer n
	character*10 n_str

		if (n.LT.0) then
			write(n_str, '(I0)') -n
			call this % formulaItem( TRIM(this%minus) // n_str )
		else
			write(n_str, '(I0)') n
			call this % formulaItem( n_str )
		end if
	
	end subroutine Output_proc_writeInteger

	subroutine Output_proc_writeIntegerExcludingPM1( this, n )
	class(Output), intent(IN) :: this	
	integer n
	character*10 n_str

		if (n.LT.0) then
			if (n.NE.-1) then
				write(n_str, '(I0)') -n
				call this % formulaItem( TRIM(this%minus) // n_str )
			else
				call this % formulaItem( this%minus )
			end if
		else
			if (n.NE.1) then
				write(n_str, '(I0)') n
				call this % formulaItem( n_str )
			end if
		end if
	
	end subroutine Output_proc_writeIntegerExcludingPM1

	subroutine Output_proc_writeIntegerExcludingPM1WithPlusSign( this, n, withPlusSign )
	class(Output), intent(IN) :: this	
	integer n
	logical withPlusSign
	character*10 n_str

		if (n.LT.0) then
			if (n.NE.-1) then
				write(n_str, '(I0)') -n
				call this % formulaItem( TRIM(this%minus) // n_str )
			else
				call this % formulaItem( this%minus )
			end if
		else
			if (withPlusSign) then
				call this % formulaItem( '+' )
			end if
			if (n.NE.1) then
				write(n_str, '(I0)') n
				call this % formulaItem( n_str )
			end if
		end if
	
	end subroutine Output_proc_writeIntegerExcludingPM1WithPlusSign


	subroutine HTML_proc_beginEquation( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)') HTML_TABLE_HEAD
		
	end subroutine HTML_proc_beginEquation

	subroutine TeX_proc_beginEquation( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)') '\begin{eqnarray}'
		write( this % f, '(A)') '\nonumber'
		
	end subroutine TeX_proc_beginEquation



	subroutine HTML_proc_endEquation( this )
	class(OutputHTML), intent(INOUT) :: this
	
		call this % showLateSubscripts()
		call this % closeBrackets()
		write( this % f, '(A)') HTML_TABLE_TAIL
		this % superscriptTextLevel = 0
		this % subscriptTextLevel = 0
		
	end subroutine HTML_proc_endEquation

	subroutine TeX_proc_endEquation( this )
	class(OutputTeX), intent(INOUT) :: this
	
		call this % closeBrackets()
		write( this % f, '(A)') '\end{eqnarray}'
		this % superscriptTextLevel = 0
		this % subscriptTextLevel = 0
		
	end subroutine TeX_proc_endEquation



	subroutine HTML_proc_beginFrac( this )
	class(OutputHTML), intent(IN) :: this	

		call this % formulaItemGreater( '(' )

	end subroutine HTML_proc_beginFrac

	subroutine TeX_proc_beginFrac( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '\frac{'
		
	end subroutine TeX_proc_beginFrac




	subroutine HTML_proc_middleFrac( this )
	class(OutputHTML), intent(IN) :: this
	
		call this % formulaItemGreater( '/' )
		
	end subroutine HTML_proc_middleFrac

	subroutine TeX_proc_middleFrac( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '}{'
		
	end subroutine TeX_proc_middleFrac




	subroutine HTML_proc_endFrac( this )
	class(OutputHTML), intent(IN) :: this
	
		call this % formulaItemGreater( ')' )
		
	end subroutine HTML_proc_endFrac

	subroutine TeX_proc_endFrac( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)') '}'
		
	end subroutine TeX_proc_endFrac



	subroutine HTML_proc_summation( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)') '<td rowspan=3 class=greater>', HTML_SUM,'</td>'
		write( this % f, '(A)') '<td rowspan=2>&nbsp;</td>' ! rowspan=1 for limits later
		
	end subroutine HTML_proc_summation

	subroutine TeX_proc_summation( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)') '\sum'
		
	end subroutine TeX_proc_summation


	subroutine Output_do_nothing( this ) ! 
	class(Output), intent(IN) :: this	
! do nothing
	end subroutine Output_do_nothing

	subroutine HTML_do_nothing( this ) ! 
	class(OutputHTML), intent(IN) :: this	
! do nothing
	end subroutine HTML_do_nothing

	subroutine TeX_do_nothing( this ) ! 
	class(OutputTeX), intent(IN) :: this	
! do nothing
	end subroutine TeX_do_nothing

	subroutine TeX_proc_beginLimits( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '\limits'
		
	end subroutine TeX_proc_beginLimits




	subroutine HTML_proc_beginFormulaItem( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '<td rowspan=3>'
		
	end subroutine HTML_proc_beginFormulaItem



	subroutine HTML_proc_endFormulaItem( this )
	class(OutputHTML), intent(IN) :: this

		write( this % f, '(A)') '</td>'
	
	end subroutine HTML_proc_endFormulaItem




	subroutine HTML_proc_beginGreater( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '<span class=greater>'

	end subroutine HTML_proc_beginGreater

	subroutine HTML_proc_endGreater( this )
	class(OutputHTML), intent(IN) :: this

		write( this % f, '(A)') '</span>'
	
	end subroutine HTML_proc_endGreater




	subroutine TeX_proc_beginText( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '\mathrm{'

	end subroutine TeX_proc_beginText

	subroutine TeX_proc_endText( this )
	class(OutputTeX), intent(IN) :: this

		write( this % f, '(A)') '}'	
		
	end subroutine TeX_proc_endText




	subroutine HTML_proc_beginSuperscript( this )
	class(OutputHTML), intent(IN) :: this
	
		if (this % superscriptTextLevel .GT. 0) then
			write( this % f, '(A)', advance='no') '<td rowspan=2 class="smaller normal">'
		else
			write( this % f, '(A)', advance='no') '<td rowspan=2 class=smaller>'
		end if
	end subroutine HTML_proc_beginSuperscript

	subroutine TeX_proc_beginSuperscript( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '^{'
		if (this % superscriptTextLevel .GT. 0) write( this % f, '(A)', advance='no') '\mathrm{'
		
	end subroutine TeX_proc_beginSuperscript



	subroutine Output_proc_switchSuperscriptText( this, b )
	class(Output), intent(INOUT) :: this
	logical b
	
		if (b) then
			this % superscriptTextLevel = this % superscriptTextLevel + 1
		else
			this % superscriptTextLevel = this % superscriptTextLevel - 1
		end if
		
	end subroutine Output_proc_switchSuperscriptText

	subroutine Output_proc_switchSubscriptText( this, b )
	class(Output), intent(INOUT) :: this
	logical b
	
		if (b) then
			this % subscriptTextLevel = this % subscriptTextLevel + 1
		else
			this % subscriptTextLevel = this % subscriptTextLevel - 1
		end if
		
	end subroutine Output_proc_switchSubscriptText




	subroutine HTML_proc_beginSubscript( this )
	class(OutputHTML), intent(IN) :: this
	
		if (this % subscriptTextLevel .GT. 0) then
			write( this % f, '(A)', advance='no') '<td class="smaller normal"><sub>'
		else
			write( this % f, '(A)', advance='no') '<td class=smaller><sub>'
		end if
	end subroutine HTML_proc_beginSubscript

	subroutine TeX_proc_beginSubscript( this )
	class(OutputTeX), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '_{'
		if (this % subscriptTextLevel .GT. 0) write( this % f, '(A)', advance='no') '\mathrm{'
		
	end subroutine TeX_proc_beginSubscript




	subroutine HTML_proc_endSuperscript( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)', advance='no') '</td>'
		
	end subroutine HTML_proc_endSuperscript

	subroutine TeX_proc_endSuperscript( this )
	class(OutputTeX), intent(IN) :: this
	
		if (this % superscriptTextLevel .GT. 0) write( this % f, '(A)', advance='no') '}'
		write( this % f, '(A)', advance='no') '}'
		
	end subroutine TeX_proc_endSuperscript




	subroutine HTML_proc_endSubscript( this )
	class(OutputHTML), intent(IN) :: this
	
		write( this % f, '(A)') '</sub></td>'
		
	end subroutine HTML_proc_endSubscript

	subroutine TeX_proc_endSubscript( this )
	class(OutputTeX), intent(IN) :: this
	
		if (this % subscriptTextLevel .GT. 0) write( this % f, '(A)', advance='no') '}'
		write( this % f, '(A)') '}'
		
	end subroutine TeX_proc_endSubscript




	subroutine HTML_proc_showLateSubscripts( this )
	class(OutputHTML), intent(INOUT) :: this
	
	integer i, tmp
	
		write( this % f, '(A)') '</tr><tr></tr><tr>'
		tmp = this % subscriptTextLevel
		do i=0, this % numSubscripts-1
			
			if (this % subscriptStorage(i) % text) then
				this % subscriptTextLevel = 1
			else
				this % subscriptTextLevel = 0
			end if
			call this % beginSubscript()
			write( this % f, '(A)', advance='no') TRIM( this % subscriptStorage(i) % s )
			call this % endSubscript()
		end do
		call this % dropSubscripts
		this % subscriptTextLevel = tmp
		
	end subroutine HTML_proc_showLateSubscripts






!---------------------------------- brackets

	subroutine Output_proc_openBracket( this, i )
	class(Output), intent(INOUT) :: this
	integer, optional :: i
	integer version, pos

		if (present(i)) then
			version = i
		else
			version = this % bracketNextVersion ! 3 versions: ( [ { ; manually also <
		end if

		if (version.LT.3) then
			this % bracketNextVersion = MOD( version + 3 - 1, 3 ) ! open in backwards order: { [ ( 2 - 1 - 0
		else
			this % bracketNextVersion = 0 ! reset levels
		end if

		pos = this % bracketLevel
		this % bracketValues( pos ) = version
		pos = pos + 1
		this % bracketLevel = pos
		
		call this % showBracket( version, .TRUE. )
	end subroutine Output_proc_openBracket

	subroutine Output_proc_closeBracket( this )
	class(Output), intent(INOUT) :: this
	integer version, pos

		pos = this % bracketLevel
		if (pos.GT.0) then
			pos = pos - 1
			this % bracketLevel = pos
			version = this % bracketValues( pos )
			call this % showBracket( version, .FALSE. )
		end if
		
	end subroutine Output_proc_closeBracket
	
	subroutine Output_proc_closePseudoBrackets( this ) ! at new line
	class(Output), intent(INOUT) :: this
	integer pos

		if (this % bracketLevel .GT. 0) then
			do pos = this % bracketLevel-1, 0, -1
				call this % showBracket( BRACKET_INVISIBLE, .FALSE. )
			end do
		end if
		
	end subroutine Output_proc_closePseudoBrackets

	subroutine Output_proc_closeBrackets( this ) ! at new line
	class(Output), intent(INOUT) :: this
	integer pos

		if (this % bracketLevel .GT. 0) then
			do pos = this % bracketLevel-1, 0, -1
				call this % showBracket( this % bracketValues(pos), .FALSE. )
			end do
		end if
		this % bracketLevel = 0 ! close forever
		this % bracketNextVersion = BRACKET_ROUND
		
	end subroutine Output_proc_closeBrackets

	subroutine Output_proc_openPseudoBrackets( this )
	class(Output), intent(INOUT) :: this
	integer pos

		if (this % bracketLevel .GT. 0) then
			do pos = 0, this % bracketLevel-1
				call this % showBracket( BRACKET_INVISIBLE, .TRUE. )
			end do
		end if
		
	end subroutine Output_proc_openPseudoBrackets


	
	subroutine HTML_proc_showBracket( this, version, opening )
	class(OutputHTML), intent(IN) :: this
	integer version
	logical opening
	
		if (opening) then ! opening
			select case(version)
				case( BRACKET_ROUND )
					call this % FormulaItemGreater( '(' )
				case( BRACKET_SQUARE )
					call this % FormulaItemGreater( '[' )
				case( BRACKET_CURLY )
					call this % FormulaItemGreater( '{' )
			! special cases
				case( BRACKET_ANGLE )
					call this % FormulaItemGreater( HTML_LEFT_LT )
!				case( BRACKET_STRAIGHT )
!					call this % FormulaItemGreater( '|' )
!				case( BRACKET_INVISIBLE )
					! do nothing
			end select
		else ! closing
			select case(version) ! version
				case( BRACKET_ROUND )
					call this % FormulaItemGreater( ')' )
				case( BRACKET_SQUARE )
					call this % FormulaItemGreater( ']' )
				case( BRACKET_CURLY )
					call this % FormulaItemGreater( '}' )
			! special cases
				case( BRACKET_ANGLE )
					call this % FormulaItemGreater( HTML_RIGHT_GT )
!				case( BRACKET_STRAIGHT )
!					call this % FormulaItemGreater( '|' )
!				case( BRACKET_INVISIBLE )
					! do nothing
			end select
		end if
	end subroutine HTML_proc_showBracket
			

	subroutine TeX_proc_showBracket( this, version, opening )
	class(OutputTeX), intent(IN) :: this
	integer version
	logical opening
	
		if (opening) then ! opening
			select case(version)
				case( BRACKET_ROUND )
					call this % FormulaItem( '\left(' )
				case( BRACKET_SQUARE )
					call this % FormulaItem( '\left[' )
				case( BRACKET_CURLY )
					call this % FormulaItem( '\left\{' )
			! special cases
				case( BRACKET_ANGLE)
					call this % FormulaItem( '\left<' )
!				case( BRACKET_STRAIGHT )
!					call this % FormulaItem( '\left|' )
				case( BRACKET_INVISIBLE )
					call this % FormulaItem( '\left.' )
			end select
		else ! closing
			select case(version) ! version
				case( BRACKET_ROUND )
					call this % FormulaItem( '\right)' )
				case( BRACKET_SQUARE )
					call this % FormulaItem( '\right]' )
				case( BRACKET_CURLY )
					call this % FormulaItem( '\right\}' )
			! special cases
				case( BRACKET_ANGLE )
					call this % FormulaItem( '\right>' )
!				case( BRACKET_STRAIGHT )
!					call this % FormulaItem( '\right|' )
				case( BRACKET_INVISIBLE )
					call this % FormulaItem( '\right.' )
			end select
		end if
!		call this % writeln( '' )
		
	end subroutine TeX_proc_showBracket



