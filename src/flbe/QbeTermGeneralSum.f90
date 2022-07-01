	
	function QTGS_init(p, useKroneckers) result(this)
	class( QbeTermGeneralSum ), pointer :: this
	type(Perturbation), pointer :: p
	logical useKroneckers
		
	integer i, i0, s_index
	logical hasKron, current_dagger

		allocate(this)
		this % title = 'general sum'
		
		this % currentSubsystemIndex = -1 ! not set
		this % useKroneckers = useKroneckers
		
		call this % setPerturbation(p)

		allocate(this % A_b(0:2, 0:this % order - 1, 0:1) )
		allocate(this % A_b_desc(0:this % order - 1, 0:1) )
		this % A_b(:,:,:) = 0
		this % A_b_desc(:,:) = -1

		allocate( this % desc(0:2) )
		this % desc(0) = 'n'
		this % desc(1) = '(n+1)'
		this % desc(2) = '(1-n)'

!		print*, ' -------------'
			do i=0, this % order - 1
				hasKron = (p%info(i)%numPossibleKroneckers.GT.0)
				s_index = p%info(i)% subsystemIndex
				! dagger => (n+1) Bose, (1-n) Fermi
				
				current_dagger = this % partDagger(s_index)
				call this % setTwoCoeffs( 0, i, this % dagger(i), p%info(i)%bose, hasKron) 
				call this % setTwoCoeffs( 1, i, .NOT.this % dagger(i), p%info(i)%bose, hasKron) 
			end do
!		end do
!		stop
	    		
	end function QTGS_init


	subroutine QTGS_destructor( this )
	type( QbeTermGeneralSum ) :: this

!		print*, 'QTGS: destructor'
		if (associated( this % p )) then
			deallocate( this % desc )
			deallocate( this % A_b_desc )
			deallocate( this % A_b )
		end if

	end subroutine QTGS_destructor

!--
	
	subroutine QTGS_setTwoCoeffs(this, half, i, dagger, bose, hasKroneckers)
	class( QbeTermGeneralSum ), intent(INOUT) :: this
	integer half, i
	logical dagger, bose, hasKroneckers
	
	integer c1, cn, cdelta, desc
	
		if (dagger) then
			c1 = 1
			if (bose) then
				cn = 1			! Bose: n+1 or n+1+delta
				desc = 1
				cdelta = 1
			else
				cn = -1			! Fermi: 1-n or 1-n-delta
				desc = 2
				cdelta = -1
			end if
		else
			c1 = 0
			cn = 1 
			if (bose) then
				cdelta = -1		! Bose: n or n-delta 
			else
				cdelta = -1		! Fermi: n or n-delta 
			end if
			desc = 0
		end if
		if (.NOT.hasKroneckers) cdelta = 0
!		print*, ' QTGS: setTwoCoeffs: i half part=', i, half, hcPart, ' c1 cn=', c1,cn
		this % A_b(0, i, half ) = c1
		this % A_b(1, i, half ) = cn
		this % A_b(2, i, half ) = cdelta
		this % A_b_desc(i, half) = desc
	
	end subroutine QTGS_setTwoCoeffs

!--

	subroutine QTGS_showSingleTerm( this, o, current_dagger )
	class( QbeTermGeneralSum ) :: this
	class( Output ), pointer :: o
	logical current_dagger

	integer i, desc, s_index, iter, v, cdelta, j, half, part
	logical as_current, first
	type( OperatorInfo ) :: info
	logical hasKroneckers
	integer cnt

		s_index = this % currentSubsystemIndex
		v = o % version
		cnt = 0
		do half=0,1	
			if (half.EQ.1) then
				call o % equationNewLine
				call o % formulaItem( o % minus )
!				cnt = 0
			end if

			first = (half.EQ.0)
			do i=0, this % order - 1

				part = half
				if (.NOT.current_dagger) part = 1-part

				desc = this % A_b_desc(i, part)
				info = this % p % info(i)

				cdelta = this % A_b(2, i, part)

				hasKroneckers = (cdelta.NE.0) .AND. this%useKroneckers
				if (hasKroneckers) then
					if (cnt.GT.10) then
						call o % equationNewLine
						call o % formulaItem( o % times )
						cnt = cnt - 10
					end if
					call o % openBracket( BRACKET_ROUND )
				end if
				call o % formulaItem( this % desc( desc ) )
				cnt = cnt + 2
				call o % upperLowerIndex( info % superscript( v ) ,  info % subscript( v ) )
				if (hasKroneckers) then
					do j=1, info%numPossibleKroneckers
						call o % beginFormulaItem
						call o % plusminus ( this % A_b(2, i, part ).GT.0 ) 
						call o % write( o % deltasmall )
						call o % endFormulaItem
						call o % superscript( ' ' )
						call o % subscript( &
&							TRIM( this % p % info(i-j) % subscript( o % version )) // &
&							TRIM( info % subscript( o % version ))  )
						cnt = cnt + 2
					end do
					call o % closeBracket
				end if
			end do
		end do
					
	end subroutine QTGS_showSingleTerm

!--

	double precision function QTGS_getValue(this, n, k, excludePosition)
	implicit none
	class( QbeTermGeneralSum ) :: this
	double precision n(0:this%order-1)
	integer k(0:this%order-1)
	integer excludePosition ! .GE.0 for RTA ; -1 for getRightPart

	double precision :: term, nvalue, factor, ans
	
	integer i, half, part, s_index, c1, cn,  cdelta, j, delta12
	logical as_current, first
	logical current_dagger
	    
!		print*, ' GS.getValue()'
		ans = 0.d0	    
		s_index = this % currentSubsystemIndex

	    current_dagger = this % partDagger( s_index )
!		print*, 'GS input n(:)=', n(:), 'k=', k(:)

		do half=0,1 ! 0 and 1    mean   (nnn - nnn)
			term = 1.d0		
!			write(*, '(A,I0,A,I0,A)', advance='no') 'getValue: hc=', &
!&				hcIndex, ' half #', half, ' is product of '

			first = (half.EQ.0)
			do i=0, this % order - 1
				
				nvalue = n(i)
				part = half
				if (.NOT.current_dagger) part = 1-part

				c1 = this % A_b(0, i, part )
				cn = this % A_b(1, i, part )
				if (i.NE.excludePosition) then
					factor = nvalue*cn + c1
					cdelta = this % A_b(2, i, part )
					delta12 = 0
					if (this % useKroneckers) then
						do j=1, this % p % info(i) % numPossibleKroneckers
							if (k(i).EQ.k(i-j)) delta12 = delta12 + 1
						end do
					end if
					factor = factor + cdelta*delta12
				else ! RTA: n+dn-> dn,  so we take   n->1, n+1->1, 1-n->-1
					factor = cn
				end if
!				select case(this % A_b_desc(i, part))
!					case(0)
!						factor = nvalue
!					case(1)
!						factor = nvalue+1.d0
!					case(2)
!						factor = 1.d0-nvalue
!				end select

!				write(*, '(F5.2,A)', advance='no') factor, ' '
!				if (delta12.NE.0) then
!					write(*, '(A,I0,A,I0,A)', advance='no') '(delta=', delta12, '*', cdelta, ') '
!				end if
				term = term * factor
			end do
!			write(*, '(A,F15.10)') '=', term 
			if (half.EQ.1) term = -term ! (nnn - nnn)
			ans = ans + term
		end do

		ans = ans	 * this % partMultiplier( s_index )	
!		write(*, '(A,F15.10)') 'result=', ans
		QTGS_getValue = ans
		
	end function QTGS_getValue


