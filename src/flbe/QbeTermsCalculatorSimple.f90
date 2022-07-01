
	subroutine QTC_Simple_destructor( this )
	use omp_lib
	type( QbeTermsCalculatorSimple ) :: this
	
	integer pp, half

!		print*, 'QTC_Simple destructor'

		do pp=0, this % numP - 1
			do half=1, this % numHalves(pp)
				deallocate( this % term(half,pp) % obj)
			end do
		end do
		deallocate( this % term )

	end subroutine QTC_Simple_destructor



	function QTC_Simple_init( sizes, hamilt, F, useKroneckers, qt_version_in ) result(this)

	use omp_lib
	
	type(Geometry), pointer :: sizes
	type(Hamiltonian), pointer :: hamilt
	type(Efactor), pointer :: F
	logical useKroneckers
	integer, optional :: qt_version_in	
	type( QbeTermsCalculatorSimple ), pointer :: this
	
	type(Perturbation), pointer :: p		
	class(QbeTermAsIs), pointer :: qt
	integer pp

	integer qt_version, numHalves, half

		allocate( this )
		
		call this % initFromGHFK( sizes, hamilt, F, useKroneckers )

		qt_version = 3 ! default version
		if (present(qt_version_in)) qt_version = qt_version_in
		
		this % qt_version = qt_version_in
		
		allocate( this % term(1:2, 0:this % numP-1) )

		do pp=0, this % numP - 1
			p => hamilt % perturbations(pp) % obj

			do half=1, this % numHalves(pp)
!				print*, ' pp half=', pp, half
				if (half.EQ.2) p => p % hc
				select case (qt_version)
					case( VERSION_AS_IS_v1 )
						qt => QbeTermAsIs( p, useKroneckers )
					case( VERSION_AS_IS_v2 )
						qt => QbeTermGeneralSum( p, useKroneckers )
					case( VERSION_AS_IS_v3 )
						qt => QbeTermSumOfABCD( p, useKroneckers )
				end select
				this % term(half,pp) % obj => qt
			end do
		end do
!		print*, ' exit'		
	end function QTC_Simple_init
		
! --

	subroutine QTC_Simple_addValue( this, t, n, dn, rta )

	use omp_lib
	
	class( QbeTermsCalculatorSimple ), intent(INOUT) :: this
	type( Occupations ), pointer :: n, dn
	double precision t
	logical rta

	integer pp, ss, i, half, i2, coeff, excludePosition

	type(Subsystem), pointer :: s
	type(Perturbation), pointer :: p
	type(MultidimensionalCycle), pointer :: mdc
	class( QbeTermAsIs ), pointer :: qt

! add terms from perturbations
		do pp=0, this % numP-1
			do ss=0, this % numS-1
				s => this % hamilt % subsystems( ss ) % obj
				if (this % term_ps(pp,ss)) then
					p => this % hamilt % perturbations(pp) % obj
					do half=1, this % numHalves(pp)
						if (half.EQ.2) p => p % hc
						i2 = p % indexOf(s)
						qt => this % term(half,pp) % obj
						call qt % setCurrentSubsystemIndex( i2 )
						excludePosition = -1
						if (rta) excludePosition = p % firstOperatorEntry( i2 )
						coeff = qt % partMultiplier( i2 )
!						print*, ' term: pp=', pp, ' ss=', ss, ' i2=', i2, ' coeff=', coeff
						if (coeff.NE.0) then
!$omp parallel do  private(i,mdc)
							do i=0, this % numT-1
								mdc => this % mdc(i,half,pp) % obj
								call mdc % run( n, qt, s % occupationIndex, &
&									rta, dn, i )
							end do
!$omp end parallel do
						end if
					end do ! halves: 1 or 2 (with H.c.)
 				end if
			end do
		end do	 

	end subroutine QTC_Simple_addValue



















! -- for debug

	subroutine QTC_Simple_show( this, o )
	class( QbeTermsCalculatorSimple ) :: this
	class( Output ), pointer :: o
	
	integer i, j, i2, half, coeff, numNonzero
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	class( QbeTermAsIs ), pointer :: qt
	character*1 a, k
	character*10 kbold, a_brackets

		call o % write( 'QbeTermsCalculator: Simple, ' )
		select case(this % qt_version)
			case( VERSION_AS_IS_v1 )
				call o % write( 'v1 (original formula)' )
			case( VERSION_AS_IS_v2 )
				call o % write( 'v2 (formula with brackets)' )
			case( VERSION_AS_IS_v3 )
				call o % write( 'v3 (sum of ABCD terms)' )
		end select
		do i=0, this % hamilt % numSubsystems-1
			s => this % hamilt % subsystems(i) % obj
			a = s % operatorLetter
			a_brackets = '(' // a // ')'
!			if (subsystem_count.EQ.1) a_brackets=''
			k = 'Y' ! s % indexLetter
			call o % bold( k, kbold )

			do j=0, this % hamilt % numPerturbations-1
				p => this % hamilt % perturbations(j) % obj
				i2 = p % indexof(s)
				if (i2.GE.0) then
					call o % beginEquation()
					call o % switchSuperscriptText( .TRUE. )
					call p % showQbeTermTitle( i2, o )
					call o % formulaItem( '=' )

					numNonzero = 0
					do half=1, this % numHalves(j)
						if (half.GT.1) then
							p => p % hc
							if (numNonzero.GT.0) call o % equationNewLine
							call o % formulaItem('+')
						end if
						qt => this % term(half,j) % obj
						call qt % setCurrentSubsystemIndex( i2 )
						coeff = qt % partMultiplier( i2 )
						call o % writeIntegerExcludingPM1( coeff )
						if ( coeff.NE.0) then
							call qt % show( o )
							numNonzero = numNonzero + 1
						end if
					end do
					call o % endEquation()

!					title => p % getQbeTermFormulaPart( i2, o )
!
!					qt => this % term(j) % obj
!					call qt % setCurrentSubsystemIndex( i2 )
!					body => qt % getFormulaPart( o )
! 
!					eqn => OutputEquation( title, body )
!
!					call o % beginEquation()
!					call eqn % show( o )
!					call o % endEquation()
				end if
			end do
		end do
	
	end subroutine QTC_Simple_show


	

