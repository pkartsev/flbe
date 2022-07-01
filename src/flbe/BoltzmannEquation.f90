module BoltzmannEquation_module


	use Output_module
	use Efactor_module
	use Geometry_module
	use Subsystem_module
	use Perturbation_module
	use QbeTermsCalculator_module
	use Hamiltonian_module
	use Problem_module
	use QbeTermsCalculator_module
	
	implicit none

	public		:: FD_RightPart
	public		:: FD_RightPart_Test
	public		:: BoltzmannEquation
	
	type, abstract :: FD_RightPart
		type( Geometry ), pointer :: sizes
		type( Hamiltonian ), pointer :: hamilt
	contains
		procedure(FD_RP_getRP), deferred :: getRightPart
	end type FD_RightPart

	type, extends(FD_RightPart) :: FD_RightPart_Test
		integer version !  v0/v1/v2 : const/exp/sin
	contains
		procedure :: getRightPart => FD_RP_Test_getRightPart
		final :: FD_RP_Test_destructor
	end type FD_RightPart_Test

	type, extends(FD_RightPart) :: BoltzmannEquation
		type( Problem ), pointer :: prob
		integer numS
		integer version !  as-is (v1/v2/v3) or fast (with FFT)

		class(QbeTermsCalculator), pointer :: qtc
		logical useKroneckers
	contains
		procedure :: getRightPart => BE_getRightPart
		procedure :: getRelaxationTimes => BE_getRelaxationTimes
		procedure :: getCoeffsForRTA => BE_getCoeffsForRTA
		procedure :: show => BE_show
		procedure :: compareTo => BE_compareTo ! equation2, occ
		
		final :: BE_destructor
	end type  BoltzmannEquation	
	
! polymorphic subroutine
	abstract interface
		subroutine FD_RP_getRP(this, t, n, dn)
			import FD_RightPart
			import Occupations
			
			class( FD_RightPart ) :: this
			double precision t
			type( Occupations ), pointer :: n, dn
		end subroutine
	end interface

	interface BoltzmannEquation
		module procedure BE_init
	end interface BoltzmannEquation

	interface FD_RightPart_Test
		module procedure FD_RP_Test_init
	end interface FD_RightPart_Test
	
contains

!----------------------------------------------------------------

	function FD_RP_Test_init( version, prob ) result(this)
	integer version
	type( Problem ), pointer :: prob

	type( FD_RightPart_Test ), pointer :: this

		allocate( this )
		this % version = version
		this % sizes => prob % sizes
		this % hamilt => prob % hamilt

	end function FD_RP_Test_init


	subroutine FD_RP_Test_getRightPart( this, t, n, dn ) ! for testing purposes
	class( FD_RightPart_Test ) :: this
	type( Occupations ), pointer :: n, dn
	double precision t

		select case( this % version )
			case(1) ! n' = 0   ->   solution = const
				dn % values = 0
			case(2) ! n' = n   ->   solution ~ exp()
				dn % values = n % values
			case(3) ! sin' = cos   ->   solution = sin() + const
				dn % values = cos(t)
		end select

	end subroutine FD_RP_Test_getRightPart


	subroutine FD_RP_Test_destructor(this)
	type( FD_RightPart_Test ) :: this
	
!		print*, ' FD_RightPart_destructor'
		nullify ( this % hamilt )
		nullify ( this % sizes )
		
	end subroutine FD_RP_Test_destructor
	

!----------------------------------------------------------------
	
	function BE_init( prob, de_in, version ) result(this)
	type( Problem ), pointer :: prob
	integer, optional :: version
	double precision, optional :: de_in
	
	type( BoltzmannEquation ), pointer :: this
	double precision de
	
		if (present(de_in)) then
			de = de_in
		else
			de = 1.d0
		end if
		allocate( this )
		
		this % prob  => prob
		this % sizes => prob % sizes
		this % hamilt => prob % hamilt
		this % numS = prob % hamilt % numSubsystems

		this % useKroneckers = prob % useKroneckers
				
		if ( present( version )) then
			this % version = version
		else
			this % version = VERSION_FAST ! by default use the fastest version
		end if
		select case( this % version )
			case( VERSION_AS_IS_v1, VERSION_AS_IS_v2, VERSION_AS_IS_v3 )
				this % qtc => QbeTermsCalculatorSimple( &
&					this % sizes, this % hamilt, prob % F, this % useKroneckers, this % version )
			case( VERSION_FAST )
				this % qtc => QbeTermsCalculatorFFT( &
&					this % sizes, this % hamilt, prob % F, this % useKroneckers, de )
		end select
		
	end function BE_init


	subroutine BE_destructor(this)
	type( BoltzmannEquation ) :: this
	
!		print*, ' BE_destructor'
		deallocate( this % qtc )
		nullify ( this % prob )
!		nullify ( this % hamilt )
!		nullify ( this % sizes )
		
	end subroutine BE_destructor
	

! -- main function

	subroutine BE_getRightPart( this, t, n, dn )
	use omp_lib	
	class( BoltzmannEquation ) :: this
	type( Occupations ), pointer :: n, dn
	double precision t

	type( Occupations ), pointer :: n0
	integer, dimension(1:3) :: k
	integer a, ss
	integer volume
	double precision rp, deltaN, tau

	type(Subsystem), pointer :: s

		volume = this % sizes % volume
		n0 => this % prob % equilibriumOccupations
		if (.NOT. this % prob % equilibriumValid) then
			call this % prob % updateEquilibriumOccupations
		end if
		
		! apply lifetime term
		do ss=0, this % numS-1
			s => this % hamilt % subsystems( ss ) % obj
			if (s % lifetimePresent) then ! start from lifetime term
				do a=0, volume-1
					call this % sizes % toVector( a, k )

					deltaN = n % getValue(k, ss) - n0 % getValue(k, ss) ! (n-n0)/tau		
					tau = s % lifetime % getValue( k )
					rp = - deltaN/ tau
					
					call dn % setValue(k, ss, rp )
				end do
			else
				dn % values(:,:,:,ss) = 0.d0 ! or set zero
			end if
				
		end do

		! add perturbation terms
		call this % qtc % addValue( t, n, dn, .FALSE. ) ! rta = .FALSE.

	end subroutine BE_getRightPart



	subroutine BE_getRelaxationTimes( this, times )
	use omp_lib	
	class( BoltzmannEquation ) :: this
	type( Occupations ), pointer :: times ! times % values = calculated times
	double precision t

		call this % getCoeffsForRTA( times )

		! n = n0 + deltaN
		! dn/dt = sum[ (n+1)(...) - n(...) ]
		! d(n0+deltaN)/dt = sum[ (n0+deltaN+1)(...) - (n0+deltaN)(...) ] = 0 + deltaN*(-1/tau)
		! d(deltaN)/dt = deltaN*coeff -> tau = -1/coeff
		times % values(:,:,:,:) = - 1.d0 / times % values(:,:,:,:)

	end subroutine BE_getRelaxationTimes



	subroutine BE_getCoeffsForRTA( this, coeffs )
	use omp_lib	
	class( BoltzmannEquation ) :: this
	type( Occupations ), pointer :: coeffs

	double precision t
	type( Occupations ), pointer :: n0

		n0 => this % prob % equilibriumOccupations
		if (.NOT. this % prob % equilibriumValid) then
			call this % prob % updateEquilibriumOccupations
		end if
		
		t = 0.d0
		coeffs % values(:,:,:,:) = 0.d0 ! set zero before adding

		! use perturbation
		! call with rta = .TRUE.
		call this % qtc % addValue( t, n0, coeffs, .TRUE. ) 

	end subroutine BE_getCoeffsForRTA





! -- for debug

	subroutine BE_show( this, o )
	class( BoltzmannEquation ) :: this
	class( Output ), pointer :: o
	
	integer i, j, i2, cnt
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	character*1 a, k
	character*10 l_id, kbold, a_brackets, acomma

		call o % horizontalLine
		call o % newSection( 'Boltzmann equation' )
if (1.EQ.1) then
		do i=0, this % hamilt % numSubsystems-1
			s => this % hamilt % subsystems(i) % obj
			a = s % operatorLetter
			a_brackets = '(' // a // ')'
			if (subsystem_count.EQ.1) a_brackets=''
			acomma = ',' // a
			if (subsystem_count.EQ.1) acomma=''
			k = s % indexLetter
			call o % bold( k, kbold )

			call o % beginEquation()
			call o % switchSuperscriptText( .TRUE. )
			
			call o % beginFrac()
			call o % formulaItem( 'dn' )
			call o % upperLowerIndex( a_brackets, kbold )
			call o % middleFrac()
			call o % formulaItem( 'dt' )
			call o % endFrac()
			call o % formulaItem( '=' )
			call o % endFormulaItem()

			cnt = 0
			if (s % lifetimePresent) then
				call o % formulaItem( o % minus )
				call o % beginFrac()
				call o % formulaItem( 'n' )
				call o % upperLowerIndex( a_brackets, kbold )
				call o % middleFrac()
				call o % formulaItem( o % tau )
				write(l_id, '(A,I0,A)') '(', s % lifetime % id, ')'
				call o % upperLowerIndex( l_id, kbold )
				call o % endFrac()
				call o % endFormulaItem()
				cnt = cnt + 1
			end if

			do j=0, this % hamilt % numPerturbations-1
				p => this % hamilt % perturbations(j) % obj
				i2 = p % indexof(s)
				if (i2.GE.0) then
					if (cnt.GT.0) call o % formulaItem( '+' )
					call o % formulaItem( 'J' )
					call o % upperLowerIndex( '('//TRIM(p%title)// TRIM(acomma) // ')', kbold )
					cnt = cnt + 1
				end if
			end do

			call o % endEquation()

		end do
	
		call o % write( 'Occupation numbers:' )
		call o % beginEquation()
		do i=0, this % hamilt % numSubsystems-1
			if (i>0) call o % equationNewLine()
			s => this % hamilt % subsystems(i) % obj
			call s % show_N_Operator( o )
		end do
		call o % endEquation()
end if
		call this % qtc % show( o )

	end subroutine BE_show


	




!----------------------------------------------------------------

	double precision function BE_compareTo( this, eq2, t, n )
	class( BoltzmannEquation ) :: this
	class( BoltzmannEquation ), pointer :: eq2
	double precision t
	type( Occupations ), pointer :: n
	
	class( Occupations ), pointer :: dn1, dn2
	
	double precision v
	integer i, j, k
	logical verbose		

		dn1 => Occupations( n )
		dn2 => Occupations( n )
		
		call this % getRightPart( t, n, dn1 )
		call eq2  % getRightPart( t, n, dn2 )
		v = dn1 % compareToRelative( dn2 )

		verbose = (v.GT.1d-8)
!		verbose = .TRUE.
		if (verbose) then
!			print*, ' dn1: sum=', (sum(dn1%values(:,:,:,i)), i=0, n%N4-1)
!			print*, ' dn2: sum=', (sum(dn2%values(:,:,:,i)), i=0, n%N4-1)
!			print*, ' sum(diff2): ', (sum((dn1%values(:,:,:,i)-dn2%values(:,:,:,i))**2), i=0, n%N4-1)
			do i=0, n%N4-1
				write(*, '(A,I0,A)', advance='no') ' (:00',i,')='
				k=0
				do k=0, SIZE(dn1%values(0,:,0,i))-1
				do j=0, SIZE(dn1%values(:,0,0,i))-1
					write(*, '(F24.13)', advance='no') dn1%values(j,k,0,i)
				end do
				print*
				end do
				print*
				write(*, '(A,I0,A)', advance='no') ' (:00',i,')='
				do k=0, SIZE(dn2%values(0,:,0,i))-1
				do j=0, SIZE(dn2%values(:,0,0,i))-1
					write(*, '(F24.13)', advance='no') dn2%values(j,k,0,i)
				end do
				print*
				end do
				print*
			end do
		end if

		deallocate( dn2 )
		deallocate( dn1 )

		BE_compareTo = v
	end function BE_compareTo

end module BoltzmannEquation_module



