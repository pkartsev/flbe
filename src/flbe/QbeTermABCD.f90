! called from constructor

	function QTSumABCD_fromPerturbation( p, useKroneckers ) result(this)
	class( QbeTermSumOfABCD ), pointer :: this
	class( QbeTermSumOfABCD ), pointer :: tmp
	type( Perturbation ), pointer :: p
	logical useKroneckers

		allocate( this )
!		call this % initFromPerturbation( p, useKroneckers )
		
		allocate( tmp )
		call tmp % initFromPerturbation( p, useKroneckers )
		call this % simplify( tmp )
		deallocate( tmp)

	end function QTSumABCD_fromPerturbation

	


	subroutine QTSumABCD_simplify(this, src) 
	class( QbeTermSumOfABCD ), intent(INOUT) :: this
	class( QbeTermSumOfABCD ), pointer :: src

	class( Perturbation ), pointer :: p

	integer nu, numS, s_index, maxNumTerms
	integer numTerms, cnt, kroIndex
	type(ABCD_info) :: ABCD
	
		this % title = 'sum of ABCD, simplified'

		p => src % p
		this % useKroneckers = src % useKroneckers
		call this % setPerturbation( p )

		this % order = src % order
		maxNumTerms = src % maxNumTerms
		this % maxNumTerms = maxNumTerms

		numS = p % numSubsystems

		allocate( this % numTerms( 0:this % kroVolume-1, 0:numS-1 ) ) ! numKro x numSubsystems
		this % numTerms(:,:) = 0

		allocate( this % info( 0:maxNumTerms-1, 0:this % kroVolume-1, 0:numS-1 ) ) ! numOperators x 2(hc) x numKro x numSubsystems
		ABCD % coeff = 0
		ABCD % NS(:) = -1 ! undefined
		this % info(:,:,:) = ABCD

		do s_index = 0, numS-1

			do kroIndex = 0, this % kroVolume - 1
					numTerms = src%numTerms(kroIndex,s_index)

					cnt = 0
					do nu=0, numTerms-1
						if (src % info(nu,kroIndex,s_index) % coeff.NE. 0) then
							cnt = cnt + 1
						end if
					end do				

					this % numTerms(kroIndex,s_index) = cnt
					cnt = 0
					do nu=0, numTerms-1
						if (src % info(nu,kroIndex,s_index) % coeff.NE. 0) then
!							print*, 's_index=', s_index, ' nonzero coeff at nu=', nu, &
!&								'NS=', src % info(nu,part,kroIndex,s_index)%NS(0:3)


							this % info(cnt,kroIndex,s_index) = src % info(nu,kroIndex,s_index)
							this % info(cnt,kroIndex,s_index) % addr = cnt
							cnt = cnt + 1
						end if
					end do	! numTerms
			end do ! kroVolume
		end do	! numS


		allocate( this % positionMultiplier ( 0: this % order - 1, 0:this%kroVolume-1 ) )
		allocate( this % replacementPosition ( 0: this % order - 1, 0:this%kroVolume-1 ) )
		allocate( this % skipIndex( 0:this % order-1, 0:this%kroVolume-1 ))
		call this % fillKroneckerInfo

	end subroutine QTSumABCD_simplify

!--

	subroutine QTSumABCD_initFromGeneralSum( this, gs )
	class( QbeTermSumOfABCD ), intent(INOUT) :: this
	class(QbeTermGeneralSum), pointer :: gs
	class( Perturbation ), pointer :: p

	integer i, nu, half, numParts, numS, s_index, kroIndex, i0, nu_norm, j, NSd2
	integer NSd, maxNumTerms, numTermsHere, modX, base, kroIndexBase
	integer coeff, partForBit, term
	logical as_current, current_dagger
	type(ABCD_info) :: ABCD, ABCD_norm

		this % title = 'sum of ABCD made from general sum'
		p => gs % p
		this % useKroneckers = gs % useKroneckers
		call this % setPerturbation( p )

		maxNumTerms = 1
		do i=0, this % order - 1
			numTermsHere = 2
			if (this%useKroneckers) then ! 2 / 3+ with kronekers
				numTermsHere = 2 + p % info(i) % numPossibleKroneckers
			end if
			maxNumTerms = maxNumTerms * numTermsHere
		end do 
		this % maxNumTerms = maxNumTerms

		numS = p % numSubsystems
		! numKro x numSubsystems
		allocate( this % numTerms( 0:this % kroVolume-1, 0:numS-1 ) ) 
		this % numTerms(:,:) = 0

		! numOperators x numKro x numSubsystems
		allocate( this % info( 0:maxNumTerms-1, 0:this % kroVolume-1, 0:numS-1 ) ) 
		ABCD % coeff = 0
		ABCD % NS(:) = -1 ! undefined
		ABCD % addr = -1 ! undefined
		this % info(:,:,:) = ABCD
		
		do s_index = 0, numS-1

				current_dagger = gs % partDagger( s_index )
				do nu=0, maxNumTerms-1
					coeff = 0
					do half=0,1
						base = 1
						kroIndexBase = 1
						kroIndex = 0
						term = 1

						partForBit = half
						if (.NOT.current_dagger) partForBit = 1-partForBit

						do i=0, this%order-1
!						do i0=0, this % order - 1
!							i = i0
!							if (half.EQ.1) i = this % order - 1 - i0
							numTermsHere = 2
							if (this%useKroneckers) then ! 2 / 3+ with kronekers
								numTermsHere = 2 &
&									+ p % info(i) % numPossibleKroneckers
							end if
							modX = MOD(nu/base, numTermsHere)
							NSd = min( modX, 2) ! 3+ remain at 2
							ABCD % NS(i) = NSd

							! NSd (0=c1 / 1=cn / 2=cdelta), operator, half, hcIndex
							term = term * gs % A_b( NSd, i, partForBit) 
							if (numTermsHere.GT.2) then ! kroneckers allowed
								if (NSd.EQ.2) then
!									kroIndex = kroIndex + power2
									kroIndex = kroIndex + kroIndexBase*(modX-1)
								end if
								kroIndexBase = kroIndexBase * (numTermsHere-1) ! 2->nothing; 3->2
							end if
!							print*, ' i NTH base=', i, numTermsHere, base
							base = base * numTermsHere
						end do ! i=0..order-1
						if (half.EQ.0) then
							coeff = term
						else
							coeff = coeff - term
						end if

					end do ! halves of (n+1)(n+1) - nn

		! switching off:
		! conflict with testQbeTerm:  not summing by k3 and k4
		! expecting problems with kroneckers
		if ( 1.EQ.2) then 
!		if ( this % p % mode .EQ. Q_INDEPENDENT) then ! then k3 and k4 can be exchanged
!			print*, ' Q_independent'
!					find common terms, e.g. NS + SN = 2*NS
					nu_norm = 0
					ABCD_norm % NS(:) = ABCD % NS(:)
					base = 1
					do i=0, this%order-1
						NSd = ABCD_norm % NS(i) 
						if (.NOT.p % info(i) % firstOperatorForSubsystem) then
							if (NSd.LT.2) then ! not kronecker (maybe)
								do j=i+1, this%order-1
									if ( ((p % info(i) % dagger).EQV.(p % info(j) % dagger)) .AND. &
&									(p % info(i) % subsystemIndex .EQ. p % info(j) % subsystemIndex) ) then
										NSd2 = ABCD_norm % NS(j) 
										if (NSd2.LT.2) then ! not kronecker (maybe)
											if (NSd.LT.NSd2) then ! exchange: SN -> NS
!												print*, '   exchange ', i, j, ' NS=', NSd, NSd2
												ABCD_norm % NS(j) = NSd
												NSd = NSd2 ! to use later
												ABCD_norm % NS(i) = NSd
											end if
										end if
									end if
								end do
							end if
						end if
						numTermsHere = 2
						if (this%useKroneckers) then ! 2 / 3+ with kronekers
							numTermsHere = 2 &
&								+ p % info(i) % numPossibleKroneckers
						end if
						nu_norm = nu_norm + NSd*base
!						print*, ' i nth base NSd=', i, numTermsHere, base, NSd
						base = base * numTermsHere						
					end do 
					if (nu_norm.LT.nu) then
						print*, 'NS_ini =', ABCD % NS(0:this%order-1)
						print*, 'NS_norm=', ABCD_norm % NS(0:this%order-1)
!						print*, ' nu -> norm: ', nu, nu_norm
						this % info(nu_norm,kroIndex,s_index) % coeff = &
&							this % info(nu_norm,kroIndex,s_index) % coeff + coeff 
						coeff = 0						
					end if
end if
					this % info(nu,kroIndex,s_index) % NS(:) = ABCD % NS(:)
					this % info(nu,kroIndex,s_index) % coeff = coeff
					this % info(nu,kroIndex,s_index) % addr = nu
				end do ! nu
				this % numTerms(:,s_index) = maxNumTerms
		end do	! s_index


		allocate( this % positionMultiplier ( 0: this % order - 1, 0:this%kroVolume-1 ) )
		allocate( this % replacementPosition ( 0: this % order - 1, 0:this%kroVolume-1 ) )
!		print*, ' replacementPosision dimensions: ', this % order, this%kroVolume
		allocate( this % skipIndex( 0:this % order - 1, 0:this%kroVolume-1 ))

	end subroutine QTSumABCD_initFromGeneralSum
	
!--



	subroutine QTSumABCD_initFromPerturbation( this, p, useKroneckers )  ! first allocate(this)
	class( QbeTermSumOfABCD ), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	logical useKroneckers
	
	class(QbeTermGeneralSum), pointer :: gs
	
		gs => QbeTermGeneralSum( p, useKroneckers )
		call QTSumABCD_initFromGeneralSum( this, gs )
		deallocate(gs)

	
	end subroutine QTSumABCD_initFromPerturbation




	subroutine QTSumABCD_destructor( this )
	type( QbeTermSumOfABCD ) :: this
	
!		print*, 'QbeTermSumOfABCD: destructor'
		if ( associated( this% p )) then
			deallocate( this % info )
			deallocate( this % numTerms )
			deallocate( this % skipIndex )
			deallocate( this % replacementPosition )
			deallocate( this % positionMultiplier )
		end if
		
	end subroutine QTSumABCD_destructor
	
!----------------------

	subroutine QTSumABCD_show( this, o )
	class( QbeTermSumOfABCD ) :: this
	class( Output ), pointer :: o
		
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	integer s_index, part, cnt, kroIndex
	logical current_dagger
	    
	    p => this % p
		s_index = this % currentSubsystemIndex
		s => p % subsystems( s_index ) % obj
!	    current_dagger = p % info( s_index )%dagger
	    current_dagger = this % dagger( s_index )

		call o % formulaItem( o%twopi_hbar )
		
		if (this % p % mode .EQ.Q_INDEPENDENT ) then
			call o % formulaItem( '|V' )
			call o % upperLowerIndex( ' ', this  % p % title )
			call o % formulaItem( '|' )
			call o % upperLowerIndex( '2', ' ' )
		end if
		call o 	% openBracket( BRACKET_SQUARE )
		
		cnt = 0
		do kroIndex =0, this%kroVolume-1
			if (cnt.GE.4) then
				call o % equationNewLine()
				cnt = 0
			end if
			if (kroIndex.GT.0) then
				call o % equationNewLine()
				call o % formulaItem('+')
				cnt = 0
			end if
!			call this % useKroneckerIndex( kroIndex )

			call this % showHeadForKronecker( o, kroIndex )
			call this % showSingleTermForKronecker( o, current_dagger, kroIndex )
			cnt = cnt + 1
			call this % showTail( o )		
		end do ! kroIndex
		
		call o 	% closeBracket()

	end subroutine QTSumABCD_show	


	subroutine QTSumABCD_showHeadForKronecker( this, o, kroIndex )
	class( QbeTermSumOfABCD ) :: this
	class( Output ), pointer :: o
	integer kroIndex
	
	logical :: skipIndex(0:this%order-1)
	integer i, s_index

		s_index = this % currentSubsystemIndex

		skipIndex(:) = this % skipIndex(:, kroIndex)
		do i=0, this % order-1
			if (( s_index .EQ. this%p%info(i)%subsystemIndex) .AND. ( this%p%info(i)%indexPos.EQ.0)) then
				skipIndex(i) = .TRUE.
			end if
		end do
		
		call this % p % showSumHead( o, skipIndex )

	end subroutine QTSumABCD_showHeadForKronecker





	subroutine QTSumABCD_fillKroneckerInfo( this )
	class( QbeTermSumOfABCD ), intent(INOUT) :: this
	
	integer i, j, base, numTermsHere, order, kroIndex

		this % positionMultiplier(:,:) = 1
		this % replacementPosition(:,:) = -1
		this % skipIndex(:,:) = .FALSE.
		
		do kroIndex=0, this % kroVolume-1
			base = 1
			do order=0, this % p % numDeltaPositions-1
				i = this % p % deltaPositions(order)
				numTermsHere = this % p % deltaNumCounterparts( order )
				j = i - MOD(kroIndex/base, numTermsHere)
				if (j .LT. i ) then							
					do while ( this % replacementPosition( j, kroIndex ).GE.0) 
						j = this % replacementPosition( j, kroIndex )
					end do 
					this % positionMultiplier( j, kroIndex ) = &
&						this % positionMultiplier( j, kroIndex ) + &
&						this % positionMultiplier( i, kroIndex )
					this % positionMultiplier( i, kroIndex ) = 0
					this % skipIndex( i, kroIndex ) = .TRUE.
					this % replacementPosition( i, kroIndex ) = j
!					if (kroIndex.EQ.3) then
!						print*, ' j i=', j, i
!print*, ' posMu=', this % positionMultiplier( :, kroIndex ), ' repPos=', this % replacementPosition(:, kroIndex )
!					end if
				end if
				base = base * numTermsHere
			end do
!			if (kroIndex.EQ.3) then
!print*, 'FINALLY posMu=', this % positionMultiplier( :, kroIndex ), &
!&	' repPos=', this % replacementPosition(:, kroIndex )
!			end if
		end do
	
	end subroutine QTSumABCD_fillKroneckerInfo




	subroutine QTSumABCD_writeKroneckerFactor( this, kroIndex, o )
	class( QbeTermSumOfABCD ) :: this
	integer kroIndex
	class( Output ), pointer :: o

	integer i, j, base, numTermsHere, v, order

		v = o % version
		base = 1
		do order=0, this % p % numDeltaPositions-1
			i = this % p % deltaPositions(order)
			numTermsHere = this % p % deltaNumCounterparts( order )
			j = i - MOD(kroIndex/base, numTermsHere)
			if (j.LT.i) then
				call o % formulaItem( o % deltasmall )
				call o % upperLowerIndex( ' ', &
&	TRIM(this % p % info(j) % subscript( v )) // TRIM(this % p % info(i) % subscript( v )) )
			end if
			base = base * numTermsHere
		end do

	end subroutine QTSumABCD_writeKroneckerFactor



	function QTSumABCD_calculateKroneckerFactor( this, k, kroIndex )  result(r)
	class( QbeTermSumOfABCD ) :: this
	integer k(0:this%order-1)
	integer kroIndex

	integer r
	
	integer i, j, base, numTermsHere, order

		r = 1
		base = 1
		do order=0, this % p % numDeltaPositions-1
			i = this % p % deltaPositions(order)
			numTermsHere = this % p % deltaNumCounterparts( order )
			j = i - MOD(kroIndex/base, numTermsHere)
			if (j.LT.i) then
				if ( k(i) .NE. k(j)) r = 0 ! TODO: maby later logical TRUE/FALSE : if ((integer).EQ.0) vs if (logical)
			end if
			base = base * numTermsHere
		end do

	end function QTSumABCD_calculateKroneckerFactor




	subroutine QTSumABCD_showVdKdEforKroneckers( this, o, kroIndex )
	class( QbeTermSumOfABCD ) :: this
	class( Output ), pointer :: o
	integer kroIndex
	
	integer i, pm, s_index, versionForQ, pos(0:1), multip
	type(OperatorInfo), pointer :: info
	type( Subsystem ), pointer :: s
	character a
	character*10 askobki, multip_str
	
		s_index = this % currentSubsystemIndex
				
		versionForQ = 0
		if (s_index.EQ.0) versionForQ = 1
!		print*, 'QT as is : show : s_index=', s_index, ' versionForQ=', versionForQ
		call o % equationNewLine()
		call o % formulaItem( o % times )

		if ( this % p % mode .EQ. Q_DEPENDENT) then
			call o % formulaItem( '|V' )
			call o % upperLowerIndex( ' ', this % p % title, .FALSE., .TRUE. )
			call o % formulaItem( '(' )
!			print*, ' pos=', this % p % operatorPositionForQ(:,:)
			versionForQ=-1
			do i=1,0,-1
				pos(:) = this % p % operatorPositionForQ(:,i)
				if (pos(0).GE.0) then
					if ( (pos(0).NE.s_index).AND.(pos(1).NE.s_index) ) versionForQ=i
				end if
			end do
!			print*,' versionForQ=', versionForQ
			if (versionForQ.GE.0) then
				pos(:) = this % p % operatorPositionForQ(:,versionForQ)
				if (this % positionMultiplier(pos(0), kroIndex).EQ.0) then
					pos(0) = this % replacementPosition( pos(0), kroIndex )
				end if
!				print*, ' pos=', pos
				if (pos(1).GE.0) then
					if (this % positionMultiplier(pos(1), kroIndex).EQ.0) then
						pos(1) = this % replacementPosition( pos(1), kroIndex )
					end if
				end if
!				print*, ' versionForQ=', versionForQ, ' pos=', pos
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
												
		do i=0, this % order-1
			multip = this%positionMultiplier(i, kroIndex)
			if (multip.NE.0) then
				info => this % p % info(i)
				pm = info % pm
				if ((i.GT.0) .OR. (pm.LT.0)) then
					call o % plusminus( pm.GT.0 )
				end if
				
				if (multip.NE.1) then
					write( multip_str, '(I0)' ) multip
					call o % write( TRIM( multip_str ) )
				end if
				call o % write( info % subscript( o % version ) )
			end if
		end do
		call o % write( ')' )
		call o % endFormulaItem()
		call o % formulaItem( 'F(' )
		do i=0, this % order-1
			multip = this%positionMultiplier(i, kroIndex)
			if (multip.NE.0) then
				info => this % p % info(i)
				pm = info % pm
				a = info % a
				if ((i.GT.0) .OR. (pm.LT.0)) then
!					call o % plusminus( pm.GT.0 )
					call o % formulaItem( this % p % pm(o % version, pm) )
				end if
				s => this % p % subsystems( info % subsystemIndex ) % obj
				askobki = '(' // a // ')'
				if (subsystem_count.EQ.1) askobki = ''
				if (multip.NE.1) then
					write( multip_str, '(I0)' ) multip
					call o % formulaItem( TRIM( multip_str ) )
				end if
				call s % dispersionLaw % showTitle( askobki, o )
				call o % subscript( info % subscript( o % version ) )
			end if
		end do
		call o % formulaItem( ')' )
				
	end subroutine QTSumABCD_showVdKdEforKroneckers









	subroutine QTSumABCD_showSingleTermForKronecker( this, o, current_dagger, kroIndex )
	class( QbeTermSumOfABCD ) :: this
	class( Output ), pointer :: o
	logical current_dagger
	integer kroIndex

	integer i, s_index, coeff, pos, cnt, numTerms, v, NSd
	type(OperatorInfo) :: info
	type(ABCD_info) :: ABCD
	character*20 coeff_str
	logical nothing

		s_index = this % currentSubsystemIndex
		v = o % version
		
		cnt = 0

		numTerms = this % numTerms(kroIndex,s_index)
!		if (numTerms.GT.0) then
			call o % equationNewLine()
			call o % formulaItem( '+' )

			if (kroIndex.GT.0) then
				call o % openBracket( BRACKET_ROUND )
				call this % writeKroneckerFactor( kroIndex, o )
!			write(coeff_str, '(I0)') numTerms
!			call o % formulaItem( 'numTerms=' // TRIM(coeff_str) )
				call o % closeBracket
			end if
!			call this % showHeadForKronecker( o, kroIndex )
			call o % writeIntegerExcludingPM1( this % partMultiplier( s_index ) )
			call o 	% openBracket( BRACKET_ROUND )
		
			do pos=0, numTerms-1
				ABCD = this % info(pos,kroIndex,s_index)
				coeff = ABCD % coeff
!				print*, '     coeff=', coeff
!				if (use_hc) coeff = -coeff
				if (coeff .NE. 0) then
					if ((cnt.GE.5)) then
						call o % equationNewLine()
						cnt = 0
					end if
					cnt = cnt + 1
					if ((pos.GT.0).OR.(coeff.LT.0)) then
						call o % beginFormulaItem()
						call o % plusminus( coeff.GT.0 )
						call o % endFormulaItem()
					end if
					
					if (coeff.GT.0) then
						if (coeff.NE.1)	then
							write(coeff_str, '(I0)') coeff
							call o % formulaItem( coeff_str )
						end if
					else
						if (coeff.NE.-1) then
							write(coeff_str, '(I0)') -coeff
							call o % formulaItem( coeff_str )
						end if
					end if

					nothing = .TRUE.
					do i=0, this%order-1
						NSd = ABCD % NS(i)
						info = this % p % info(i)
						if (NSd.EQ.1) then
							call o % formulaItem( 'n' )
							call o % upperLowerIndex( info % superscript( v ), info % subscript( v ) )
							nothing = .FALSE.
						end if
					end do
					if (nothing) call o % formulaItem( '1' )
				end if
			end do
			call o % closeBracket
			call this % showVFDeltaKDeltaE( o )
!			call this % showTail( o )		
!		end if
	end subroutine QTSumABCD_showSingleTermForKronecker


!--

	double precision function QTSumABCD_getValue( this, n, k, excludePosition )
	class( QbeTermSumOfABCD ) :: this
	double precision n(0:this%order-1)
	integer k(0:this%order-1)
	integer excludePosition ! .GE.0 for RTA ; -1 for getRightPart
	
	integer kroIndex, maxKroIndex, s_index
	double precision :: ans
		
		ans = 0.d0

		s_index = this % currentSubsystemIndex

		maxKroIndex = 0
		if (this % useKroneckers) then
			maxKroIndex = this%kroVolume-1
		end if
!		print*, ' QT ABCD useKro=', this % useKroneckers,' maxKro=', maxKroIndex
		do kroIndex = 0, maxKroIndex

			if (this % calculateKroneckerFactor( k, kroIndex ).NE.0) then
				ans = ans + this % getValueForKronecker( n, kroIndex, excludePosition )
			end if
		end do
		
		ans = ans * this % partMultiplier( s_index )

		QTSumABCD_getValue = ans
		
	end function QTSumABCD_getValue



	function QTSumABCD_addValueForRTAForKronecker( &
&		this, n, k, kroIndex, occupationIndex ) result(ans)
	class( QbeTermSumOfABCD ) :: this
	double precision n(0:this%order-1)
	integer k(0:this%order-1)
	integer kroIndex
	integer occupationIndex
	
	double precision ans
	integer i, first, s_index
	
		ans = 0.d0

		s_index = this % currentSubsystemIndex

!
!	reference formula from addValue() : 
!
!		do i=0, this % order-1
!			if (this % subsystemIndex(i) .EQ. subsystemIndex) then
!				ans = ans + this % getValueForKronecker( n, kroIndex, i)
!			end if
!		end do
		first = this % p % firstOperatorEntry( occupationIndex ) 
		if (this % calculateKroneckerFactor( k, kroIndex ).NE.0) then
			do i=0, this % order-1
				if ((kroIndex.GT.0).OR.(i.NE.first)) then
					if (this % subsystemIndex(i) .EQ. this % currentSubsystemIndex) then
						if (k(i).EQ.k(first)) then
							ans = ans + this % getValueForKronecker( n, kroIndex, i)
						end if
					end if
				end if
			end do
		else 
			print*, ' k=', k, ' kroIndex=', kroIndex, ' n=', n
			print* , ' this % p % deltaPositions=', this % p % deltaPositions(:)
			stop 'unexpected QtABCD kroneckerFactor=0'
		end if

	end function QTSumABCD_addValueForRTAForKronecker
	
	function QTSumABCD_getValueForKronecker( this, n, kroIndex, excludePosition ) result(ans)
	class( QbeTermSumOfABCD ) :: this
	double precision n(0:this%order-1)
	integer kroIndex
	integer excludePosition ! .GE.0 for RTA ; -1 for getRightPart
	double precision ans
	
	integer s_index, numTerms, pos, coeff, i, NSd
	type(ABCD_info) :: ABCD
	double precision term, factor
	logical doit

		ans = 0.d0

		s_index = this % currentSubsystemIndex
		numTerms = this % numTerms(kroIndex,s_index)
		do pos=0, numTerms-1
			ABCD = this % info(pos,kroIndex,s_index)
			coeff = ABCD % coeff
			if (coeff .NE. 0) then
				doit = .TRUE.
				if (excludePosition.GE.0) then
					doit = (ABCD % NS(excludePosition) .EQ. 1) ! coeff for N1
				end if
				if (doit) then
					term = 1.d0
					do i=0, this%order-1
						NSd = ABCD % NS(i)
						if (NSd.EQ.1) then ! * n
							if (i.NE.excludePosition) then ! dn/dn = 1 -> no multiply
								factor = n(i)
								term = term * factor
							end if
						end if
					end do
					ans = ans + coeff*term
				end if
			end if ! coeff<>0
		end do ! numTerms
	end function QTSumABCD_getValueForKronecker
	
	
	