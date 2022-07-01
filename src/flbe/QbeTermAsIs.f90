
! constructors: bodies

	subroutine QTAI_setPerturbation(this, p)
	class( QbeTermAsIs ), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	
	integer i, numParts, s_index
	integer numCre, numAnn

		this % p => p
		call p % prepare		

		this % order = p % numAllOperators
		
		allocate(this % subsystemIndex(0:this % order - 1) )
		allocate(this % subsystemHasDaggers(0:p%numSubsystems - 1) )
		allocate(this % bose(0:this % order - 1) )
		allocate(this % dagger(0:this % order - 1) )
		this % subsystemHasDaggers(:) = .FALSE.
		do i=0, this % order - 1
			s_index = p % info(i) % subsystemIndex
			this % subsystemIndex(i) = s_index
			if (p % numOperators(CREATION, s_index).GT.0) then
				this % subsystemHasDaggers(s_index) = .TRUE.
			end if
			this % bose(i) = p % info(i) % bose
			this % dagger(i) = p % info(i) % dagger
		end do

!		allocate(this % numParts(0:p % numSubsystems - 1) )
		allocate(this % partMultiplier(0:p % numSubsystems - 1) )
		allocate(this % partDagger(0:p % numSubsystems - 1) )
		do s_index=0, p%numSubsystems-1
			numCre = p % numOperators(CREATION, s_index)
			numAnn = p % numOperators(ANNIHILATION, s_index)
			if (numCre.GT.0) then
				this % partDagger(s_index) = .TRUE.
				this % partMultiplier(s_index) = numCre
			else
				this % partDagger(s_index) = .FALSE.
				this % partMultiplier(s_index) = numCre
			end if

		end do
!		if (this%useKroneckers) then
			this % kroVolume = p % numDeltaCombinations
!		else
!			this % kroVolume = 1
!		end if
		
	end subroutine QTAI_setPerturbation
	
!--

	function QTAI_init(p, useKroneckers) result(this)
	class( QbeTermAsIs ), pointer :: this
	type( Perturbation ), pointer :: p
	logical useKroneckers
	
		allocate(this)
		this % title = 'as-is'
		this % currentSubsystemIndex = -1 ! not set
		this % useKroneckers = useKroneckers
		call this % setPerturbation( p )
		
		
	end function QTAI_init
	
	subroutine QTAI_destructor(this)
	type( QbeTermAsIs ) :: this
	
!		print*, 'QTAI: destructor'
		if (associated( this % p )) then
			
			deallocate( this % partMultiplier )
			deallocate( this % partDagger )
			deallocate( this % dagger )
			deallocate( this % bose )
			deallocate( this % subsystemHasDaggers )
			deallocate( this % subsystemIndex )

		end if
	end subroutine QTAI_destructor
	
!--

	subroutine QTAI_setCurrentSubsystemIndex(this, currentSubsystemIndex)
	class( QbeTermAsIs ), intent(INOUT) :: this
	integer currentSubsystemIndex

		this % currentSubsystemIndex = currentSubsystemIndex

	end subroutine QTAI_setCurrentSubsystemIndex

!--

	subroutine QTAI_show( this, o )
	class( QbeTermAsIs ) :: this
	class( Output ), pointer :: o
		
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	integer s_index
	logical current_dagger
	    
		    p => this % p
			s_index = this % currentSubsystemIndex
!		print*, ' s_index=', s_index
!			call o % writeIntegerExcludingPM1( this % partMultiplier( s_index) )
			s => p % subsystems( s_index ) % obj
		    current_dagger = this % dagger( s_index )

			call o % formulaItem( o % twopi_hbar )
			if (this % p % mode .EQ.Q_INDEPENDENT) then
				call o % formulaItem( '|V' )
				call o % upperLowerIndex( ' ', this % p % title, .FALSE., .TRUE. )
				call o % formulaItem( '|' )
				call o % upperLowerIndex( '2', ' ' )
			end if

			call this % showHead( o )
!		print*, ' numParts=', this % numParts(s_index)
		    current_dagger = this % partDagger( s_index )
!print*, ' QTAI_show: part=', part, ' current_dagger=', current_dagger
!			call o % writeIntegerExcludingPM1( this % partMultiplier( s_index) )
			call o % openBracket( BRACKET_SQUARE )
			call this % showSingleTerm( o, current_dagger )
			call o % closeBracket()
			call this % showVFDeltaKDeltaE( o )
			call p % showSumTail( o )		
	
	end subroutine QTAI_show	


	subroutine QTAI_showVFDeltaKDeltaE( this, o )
	class( QbeTermAsIs ) :: this
	class( Output ), pointer :: o
	
	integer i, pm, s_index, versionForQ, pos(0:1)
	type(OperatorInfo), pointer :: info
	type( Subsystem ), pointer :: s
	character a
	character*10 askobki
	
		s_index = this % currentSubsystemIndex

		versionForQ = 0
		if (s_index.EQ.0) versionForQ = 1
!		call o % equationNewLine()
!		call o % formulaItem( o % times )

		if ( this % p % mode .EQ. Q_DEPENDENT) then
			call o % formulaItem( '|V' )
			call o % upperLowerIndex( ' ', this % p % title, .FALSE., .TRUE. )
			call o % formulaItem( '(' )
			versionForQ=-1
			do i=1,0,-1
				pos(:) = this % p % operatorPositionForQ(:,i)
				if ( (pos(0).NE.s_index).AND.(pos(1).NE.s_index) ) versionForQ=i
			end do
			if (versionForQ.GE.0) then
				pos(:) = this % p % operatorPositionForQ(:,versionForQ)
				if ( pos(0).GE.0 ) then
					if (i.GT.0) call o % formulaItem( ',' )
					call o % formulaItem( this % p % info(pos(0)) % subscript(o%version) )
					if ( pos(1).GE.0 ) then
call o % formulaItem( TRIM( o%minus ) // this % p % info(pos(1)) % subscript(o%version) )
					end if
				end if
			else
				stop 'versionForQ.LT.0'
			end if
			call o % formulaItem( ')' )
			call o % formulaItem( '|' )
			call o % upperLowerIndex( '2', ' ' )
		end if
		call o % beginFormulaItem()
		call o % write( TRIM(o % deltaSmall) // '(' )
												
		do i=0, this % order - 1
			info => this % p % info(i)
			pm = info % pm
			if ((i.GT.0) .OR. (pm.LT.0)) then
				call o % plusminus( pm.GT.0 )
			end if
			call o % write( info % subscript( o % version ) )
		end do
		call o % write( ')' )
		call o % endFormulaItem()
		call o % formulaItem( 'F(' )
		do i=0, this % order - 1
			info => this % p % info(i) ! TODO: =>
			pm = info % pm
			a = info % a
			if ((i.GT.0) .OR. (pm.LT.0)) then
				call o % formulaItem( this % p % pm(o % version, pm) )
			end if

			s => this % p % subsystems( info % subsystemIndex ) % obj
			askobki = '(' // a // ')'
			if (subsystem_count.EQ.1) askobki = ''
			call s % dispersionLaw % showTitle( askobki, o )
			call o % subscript( info % subscript( o % version ) )
		end do
		call o % formulaItem( ')' )
				
	end subroutine QTAI_showVFDeltaKDeltaE


	subroutine QTAI_showHead( this, o )
	class( QbeTermAsIs ) :: this
	class( Output ), pointer :: o

	logical :: skipIndex(0:this % order-1)
	integer i, s_index

		s_index = this % currentSubsystemIndex

		skipIndex(:) = .FALSE.
		do i=0, this % order - 1
			if (( s_index .EQ. this%p%info(i)%subsystemIndex) .AND. ( this%p%info(i)%indexPos.EQ.0)) then
				skipIndex(i) = .TRUE.
			end if
		end do
		
		call this % p % showSumHead( o, skipIndex )

	end subroutine QTAI_showHead


	subroutine QTAI_showTail( this, o )
	class( QbeTermAsIs ) :: this
	class( Output ), pointer :: o

!		call o % closeBracket
		call this % p % showSumTail( o )
		
	end subroutine QTAI_showTail


	subroutine QTAI_showSingleTerm( this, o, current_dagger )
	class( QbeTermAsIs ) :: this
	class( Output ), pointer :: o
	logical current_dagger
	
	integer i, s_index, iter, j
	logical as_current, first
	type(OperatorInfo), pointer :: info
	logical dagger, hasKroneckers, externalIndexShown
	integer cnt
		
		s_index = this % currentSubsystemIndex

		cnt = 0
		do iter=0,1	
			if (iter.EQ.1) then
				call o % equationNewLine
				call o % formulaItem( o % minus )
				cnt = 0
			end if
	
			externalIndexShown = .FALSE.
			do i=0, this % order - 1
				info => this % p % info(i)
				call o % beginFormulaItem()
				as_current = this % dagger(i).EQV.current_dagger
				if (iter.EQ.1) then
					 as_current = .NOT.as_current
				end if
				dagger = as_current
				hasKroneckers = (info%numPossibleKroneckers .GT. 0) .AND. this%useKroneckers
				if (dagger .OR. hasKroneckers) then
					if (cnt.GT.10) then
!						call o % endFormulaItem()
!						call o % equationNewLine
!!						call o % formulaItem( o % times )
!						cnt = cnt - 10
!						call o % beginFormulaItem()
					end if
					call o % write( '(' )
					if (.NOT. info % bose) then
						call o % write( '1' // o % minus )
					end if
					cnt = cnt + 2
				end if
				call o % write( 'n' )
				call o % endFormulaItem()
				call o % superscript( TRIM( info % superscript( o % version )) )
!				if ((.NOT.externalIndexShown) .AND. (s_index.EQ.info%subsystemIndex) ) then
!					call o % subscript( 'Y' )
!					externalIndexShown = .TRUE.
!				else
					call o % subscript( TRIM( info % subscript( o % version )) )
!				end if
				cnt = cnt + 2
				if (dagger .OR. hasKroneckers) then
					call o % beginFormulaItem()
					if (dagger) then
						if (info % bose) call o % write( '+1' )
					end if
					if (hasKroneckers) then
						do j=1, info%numPossibleKroneckers
							call o % plusminus( dagger .AND. info % bose )
							call o % write( o % deltasmall )
							call o % endFormulaItem()
							call o % upperLowerIndex( ' ', &
&								TRIM( this % p % info(i-j) % subscript( o % version )) // &
&								TRIM( info % subscript( o % version ))  )
							cnt = cnt + 2
							call o % beginFormulaItem()
						end do
					end if
					call o % write( ')' )
				end if
			end do
		end do
		
	end subroutine QTAI_showSingleTerm

!--

	function QTAI_addValueForRTA(this, n, k, occupationIndex) result(ans)
	class( QbeTermAsIs ) :: this
	double precision n(0:this%order-1)
	integer k(0:this%order-1)

	integer occupationIndex, i, first
	
	double precision ans
	
		ans = 0.d0
		
		first = this % p % firstOperatorEntry( occupationIndex ) 

		do i=0, this % order-1
			if (this % subsystemIndex(i) .EQ. this % currentSubsystemIndex) then
				if (k(i).EQ.k(first)) then
					ans = ans + this % getValue( n, k, i )
				end if

			end if
		end do
	end function QTAI_addValueForRTA
	
!--

	double precision function QTAI_getValue(this, n, k, excludePosition)
	class( QbeTermAsIs ) :: this
	double precision n(0:this%order-1)
	integer k(0:this%order-1)
	integer excludePosition ! .GE.0 for RTA ; -1 for getRightPart
	
	double precision :: term, nvalue, ans, factor
	
	integer i, iter, s_index, delta_12, j
	logical as_current, current_dagger
	type(OperatorInfo), pointer :: info
	logical dagger
	    
		ans = 0.d0	    

		s_index = this % currentSubsystemIndex

	    current_dagger = this % partDagger( s_index )
!		print*, 'input n(:)=', n(:), 'k=', k(:)
		do iter=0,1
			term = 1.d0		
!			write(*, '(A,I0,A,I0,A)', advance='no') 'getValue: hc=', &
!&				hcIndex, ' iter #', iter, ' is product of '
			do i=0, this % order - 1
				as_current = this % dagger(i).EQV.current_dagger
				nvalue = n(i)
				dagger = ((iter.EQ.0) .EQV. as_current)
				delta_12 = 0
				if (this % useKroneckers) then
					info => this % p % info(i)
					do j=1, info%numPossibleKroneckers
						if (k(i).EQ.k(i-j)) then
							delta_12 = delta_12 + 1
						end if
					end do
				end if
				if (i.NE.excludePosition) then				
!					if (.NOT.this % useKroneckers) then
!						delta_12 = 0
!					end if
					if (dagger) then  ! a(+) -> (n+1) / (1-n)
						if (this%bose(i)) then
							factor = nvalue + 1.d0 + delta_12	! Bose : (n+1+delta)
						else
							factor = 1.d0 - nvalue - delta_12	! Fermi : (1-n-delta)
						end if
					else
						factor = nvalue - delta_12 ! Bose and Fermi : n-delta 
					end if
				else
					factor = 1.d0  !+ delta_12
					if (dagger) then  ! a(+) -> (n+1) / (1-n)
						if (.NOT.this%bose(i)) then
							factor = -1.d0 ! Fermi: (1-n)
						end if
					end if
					! with delta12:  (n+dn+1)*(n+dn+2) = ... + dn*(2n+3)
					
				end if
!				write(*, '(F5.2,A)', advance='no') factor, ' '
!				if (delta_12.NE.0) then
!					write(*, '(A,I0,A)', advance='no') '(delta=', delta_12, ') '
!				end if
				term = term * factor
			end do
!			write(*, '(A,F15.10)') '=', term 
			if (iter.EQ.1) term = -term ! (nnn-nnn)
			ans = ans + term
		end do
		
		ans = ans	 * this % partMultiplier( s_index )
!		write(*, '(A,F15.10)') 'result=', ans
		
		QTAI_getValue = ans
		
	end function QTAI_getValue

