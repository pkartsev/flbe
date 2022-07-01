

!-------------	

	function QbeTermFFT_fromPerturbation( p, useKroneckers ) result(this) ! reuse ABCD(parent) constructor
	class( QbeTermFFT ), pointer :: this
	class( QbeTermSumOfABCD ), pointer :: tmp
	type( Perturbation ), pointer :: p
	logical useKroneckers
	
		tmp => QbeTermSumOfABCD( p, useKroneckers )

		allocate( this )
		call this % simplify( tmp )
		deallocate( tmp)
!		call this % setABCDforRTA

	end function QbeTermFFT_fromPerturbation

!----------------------




	subroutine QbeTermFFT_show( this, o )
	class( QbeTermFFT ) :: this
	class( Output ), pointer :: o
		
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	integer s_index, NS, cnt, v, kroIndex, maxk1Multiplier, k1multiplier
	integer countTerms, multip, externalIndex, i
	logical current_dagger, with_Nk
	type(OperatorInfo), pointer :: info
	    
	logical useABCD, combinedMode
	
		v = o % version
	    p => this % p
		s_index = this % currentSubsystemIndex
		externalIndex = this % p % firstOperatorEntry(s_index)
!		print*, ' QbeTermFFT_show: s_index=', s_index, ' externalIndex=', externalIndex
!		print* , 'subs='
!		do i=0, p%numAllOperators-1
!			print* , '   i=', i, ' subs=', TRIM(this%p%info(i)%subscript(v,0) )
!		end do
		s => p % subsystems( s_index ) % obj
		info => p % info( s_index )
!	    current_dagger = this % dagger( s_index )

		call o % formulaItem( '=' // o % twopi_hbar)

		combinedMode = ((p%mode .EQ. Q_DEPENDENT) .AND. (this%kroVolume.GT.1))
		if (combinedMode) then
		
			call o % openBracket( BRACKET_ROUND )
			do kroIndex = 0, this%kroVolume-1
				useABCD = ((p%mode .EQ. Q_DEPENDENT) .AND. (kroIndex.GT.0) )
				if (useABCD) then
					    current_dagger = this % partDagger( s_index )
						call o % writeIntegerExcludingPM1( this % partMultiplier( s_index ))
						call this % QbeTermSumOfABCD % showSingleTermForKronecker( &
&							o, current_dagger, kroIndex )
				else
				
					if (kroIndex.GT.0) call o % formulaItem('+')
					do NS=0,1
						with_Nk = (NS.EQ.0) 
						call o % equationNewLine()
						if (NS.GT.0) call o % formulaItem( '+' )


						if (with_Nk) then
							call o % formulaItem( 'n' )
							info => this % p % info( externalIndex )
							call o % upperLowerIndex( info % superscript( v ), info % subscript( v ) )
						end if
			
						call this % showHead( o )
						cnt = 0
			
						call o % openBracket( BRACKET_CURLY )	
						multip = this % partMultiplier( s_index) 
					    current_dagger = this % partDagger( s_index )
						call this % showTermsForKronecker( o, &
&							current_dagger, with_Nk, cnt, kroIndex, multip )

						call o % closeBracket
						call this % showVdKdEforKroneckers( o, kroIndex )
						call this % showTail( o )
					end do ! N or S
				end if ! useABCD
			end do ! kroIndex
			call o % closeBracket
		else
			call o % formulaItem( '|V' )
			call o % upperLowerIndex( ' ', this  % p % title )
			call o % formulaItem( '|' )
			call o % upperLowerIndex( '2', ' ' )
			call o % openBracket( BRACKET_ROUND )

			maxk1Multiplier = 1
			if (this%useKroneckers) then
				maxk1Multiplier = this % p % maxNumEqualMomenta
			end if

			do NS=0,1
				with_Nk = (NS.EQ.0)
				if (with_Nk) then
					call o % formulaItem( 'n' )
					info => this % p % info( externalIndex )
!					print*, ' exter=', externalIndex, ' subs=', info % subscript( v, 0 )
					call o % upperLowerIndex( info % superscript( v ), info % subscript( v ) )
					call o % openBracket( BRACKET_CURLY )	
					call o % equationNewLine()
				else
					call o % equationNewLine()
					call o % formulaItem( '+' )
				end if
				do k1multiplier=1, maxK1multiplier
					countTerms = this % countTerms( current_dagger, with_Nk, k1multiplier )
			if (countTerms.GT.0) then
					if (k1multiplier.GT.1) then
						call o % equationNewLine()
						call o % formulaItem( '+' )
					end if
					call this % showHeadWithMultiplier( o, k1multiplier, externalIndex )
					if (countTerms.GT.-1) call o % openBracket( BRACKET_CURLY )
					cnt = 0
					do kroIndex = 0, this%kroVolume-1
						if (k1multiplier .EQ. this % positionMultiplier( s_index, kroIndex ) ) then
						    current_dagger = this % partDagger( s_index )
							multip = this % partMultiplier( s_index) 
							call this % showTermsForKronecker( o, &
&								current_dagger, with_Nk, cnt, kroIndex, multip )
						end if ! good k1multiplier
					end do ! kroIndex
					if (countTerms.GT.-1) then
						call o % closeBracket
					end if
					call o % formulaItem( 'F(' // TRIM(o % gamma) // ')' )
					call this % showTail( o )
			end if
				end do ! k1multiplier
				if (with_Nk) then
					call o % closeBracket
				end if
			end do ! N or S
			call o % closeBracket
		end if
		
	
	end subroutine QbeTermFFT_show	


	function QbeTermFFT_countTerms( this, &
&		current_dagger, with_Nk, k1multiplier ) result(cnt)
	implicit none
	class( QbeTermFFT ) :: this
	logical current_dagger, with_Nk
	integer cnt, k1multiplier

	integer kroIndex, pos, coeff, s_index, b, numTerms, externalIndex
	type(ABCD_Info) :: ABCD

		cnt = 0

		s_index = this % currentSubsystemIndex
		externalIndex = this % p % firstOperatorEntry(s_index)
		do kroIndex = 0, this % kroVolume - 1
			if (k1multiplier .EQ. this % positionMultiplier( s_index, kroIndex ) ) then
				numTerms = this % numTerms(kroIndex,s_index)
				do pos=0, numTerms-1
					ABCD = this % info(pos, kroIndex,s_index)
					coeff = ABCD % coeff
					if (coeff .NE. 0) then
!						b = ABCD % NS(s_index)
						b = ABCD % NS(externalIndex) ! THIS RIGHT
						if ( (b .EQ. 1) .EQV. with_Nk) then
							cnt = cnt + 1
						end if
					end if ! coeff<>0
				end do ! numTerms
			end if ! k1multiplier
		end do ! kroIndex
!		print*, ' countTerms: ', cnt

	end function QbeTermFFT_countTerms


	
	subroutine QbeTermFFT_showTermsForKronecker( this, o, &
&		current_dagger, with_Nk, cnt, kroIndex, multip )
	class( QbeTermFFT ) :: this
	class( Output ), pointer :: o
	logical current_dagger, with_Nk
	integer, intent(INOUT) :: cnt
	integer kroIndex, multip

	integer i, s_index, coeff, pos, numTerms, b, externalIndex
	integer localCnt, round, startRound, posForQ(0:1), versionForQ, localCntBoundary
	logical as_current, nothing, good
	type(OperatorInfo), pointer :: info
	type(ABCD_Info) :: ABCD
	character*20 coeff_str
	
		s_index = this % currentSubsystemIndex

		externalIndex = this % p % firstOperatorEntry(s_index)
		
		numTerms = this % numTerms(kroIndex,s_index)
		localCnt = 0
		localCntBoundary = 12
		if (this % p % mode.EQ.Q_DEPENDENT) localCntBoundary = 1
		do pos=0, numTerms-1
				ABCD = this % info(pos, kroIndex,s_index)
				coeff = ABCD % coeff * multip
				if (coeff .NE. 0) then
					b = ABCD % NS(externalIndex)
					if ( (b .EQ. 1) .EQV. with_Nk) then
!						print*, ' term coeff=', coeff
						if (cnt.GE.localCntBoundary) then
							call o % equationNewLine()
						end if
						call o % beginFormulaItem()
						call o % plusminus( coeff.GT.0 )
						call o % endFormulaItem()
						if (cnt.GE.localCntBoundary) then
							cnt = 0
						end if
					
						nothing = .TRUE.
						if (coeff.GT.0) then
							if (coeff.NE.1)	then
								nothing = .FALSE.
								write(coeff_str, '(I0)') coeff
								call o % formulaItem( coeff_str )
							end if
						else
							if (coeff.NE.-1) then
								nothing = .FALSE.
								write(coeff_str, '(I0)') -coeff
								call o % formulaItem( coeff_str )
							end if
						end if

						startRound = 1
						if (this % p % mode.EQ.Q_DEPENDENT) then
							startRound = 0

							call o % formulaItemText( 'FFT3D' )
							call o % upperLowerIndex( TRIM(o%minus) // '1' , &
&								'q' // TRIM(o%right_arrow) // 'r' )
							call o % openBracket( BRACKET_SQUARE )

							call o % formulaItem( '|V' )
							call o % upperLowerIndex( '', this  % p % title )
							call o % formulaItem( '(q)|' )
							call o % upperLowerIndex( '2', ' ' )
							
							call o % formulaItemText( 'FFT3D' )
							call o % upperLowerIndex( '' , &
&								'r' // TRIM(o%right_arrow) // 'q' )

							call o % openBracket( BRACKET_CURLY )

							versionForQ = 0
if ( (s_index.EQ.this%p%operatorPositionForQ(0,0)) .OR. (s_index.EQ.this%p%operatorPositionForQ(1,0))) then
								versionForQ = 1
end if
							posForQ(:) = this % p % operatorPositionForQ(:,versionForQ)

						else
							posForQ(:) = -1
						end if
						do round=startRound, 1
						
						do i=0, this%order-1
							info => this % p % info(i)
							if (i.NE.externalIndex) then
								good = ((i.NE.posForQ(0)).AND.(i.NE.posForQ(1)))
								if (round.EQ.0) good = .NOT.good
								if (good) then
								  coeff = this % positionMultiplier( i, kroIndex ) 
								  if (coeff.NE.0) then
									nothing = .FALSE.

									if (ABCD % NS(i).EQ.1) then
										call o % formulaItem( 'N' )
									else
										call o % formulaItem( 'S' )
									end if
									cnt = cnt + 1

									if (coeff.GT.1) then
										write(coeff_str, '(I0)') coeff
									else
										coeff_str = ''
									end if
									call o % superscript( info % superscript( o % version ) )
									as_current = this % dagger(i).EQV.current_dagger
									if (as_current) then
										call o % subscript( TRIM(o%minus) // TRIM(coeff_str) // 'R' )
									else
										call o % subscript( TRIM(coeff_str) // 'R' )
									end if
								  end if ! not skip
								end if ! good for this round
							end if ! not s_index
						end do ! i
							if (round.EQ.0) then
								call o % closeBracket
								call o % closeBracket
							end if
						end do ! round
						if (nothing) then
							call o % formulaItem('1')
						else
							cnt = cnt + 1
							localCnt = localCnt + 1
						end if
						
					end if ! S or N
				end if ! coeff.NE.0
		end do ! pos=0..power2-1
!		print*, ' localCnt=', localCnt
		
	end subroutine QbeTermFFT_showTermsForKronecker
	
	
	subroutine QbeTermFFT_showHeadWithMultiplier( this, o, multiplier, externalIndex )
	class( QbeTermFFT ) :: this
	class( Output ), pointer :: o
	integer multiplier, externalIndex
	
	type(OperatorInfo), pointer :: info	
	character*10 multip_str
			
		if (multiplier.GT.1) then
			write(multip_str, '(I0)') multiplier
		else
			multip_str = ''
		end if
		info => this % p % info( externalIndex )
		
!		call o % writeIntegerExcludingPM1( coeff )
		call o % formulaItemText( 'FFT4D' )
		call o % upperLowerIndex( '', TRIM(multip_str) // &
&			TRIM( info % subscript( o%version) ) )
		call o % openBracket( BRACKET_SQUARE )

	end subroutine QbeTermFFT_showHeadWithMultiplier

	subroutine QbeTermFFT_showTail( this, o )
	class( QbeTermFFT ) :: this
	class( Output ), pointer :: o

		call o % closeBracket

	end subroutine QbeTermFFT_showTail

	subroutine QbeTermFFT_showVFDeltaKDeltaE( this, o )
	class( QbeTermFFT ) :: this
	class( Output ), pointer :: o
	
	integer s_index
	
		s_index = this % currentSubsystemIndex

		call o % equationNewLine()

		call o % formulaItem( o % times )

		if ( this % p % mode .EQ. Q_INDEPENDENT ) then
		end if

		call o % formulaItem( 'F(' // o % gamma // ')' )
				
	end subroutine QbeTermFFT_showVFDeltaKDeltaE

	
!--

