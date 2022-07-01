module Efactor_module

	use Output_module
	implicit none

	integer, parameter :: Efactor_LORENTZ = 1
	integer, parameter :: Efactor_GAUSS = 2
	integer, parameter :: Efactor_EXACT = 0	
	double precision, parameter :: PI = 4.d0*atan(1.d0)

	type Efactor
		double precision sigma, emax, step
		integer mode
		integer numPoints
		double precision, allocatable :: table(:)
		
	contains
		procedure :: getValue => Efactor_getValue
		procedure :: show => Efactor_show
		final :: Efactor_destructor
	end type Efactor

! constructors

	interface Efactor
		module procedure Efactor_init
	end interface Efactor
	interface EfactorGauss
		module procedure EfactorGauss_init
	end interface EfactorGauss
	interface EfactorLorentz
		module procedure EfactorLorentz_init
	end interface EfactorLorentz
	interface EfactorExact
		module procedure EfactorExact_init
	end interface EfactorExact


contains

! constructors: body
!  Efactor( Efactor_GAUSS, 0.01 )
! or
!  Efactor( Efactor_GAUSS, 0.01, 10, 40000 ) to use table and linear interpolation

	function EfactorExact_init( sigma ) result(this)
	type(Efactor), pointer :: this
	double precision sigma

		this => Efactor( Efactor_EXACT, sigma )
		
	end function EfactorExact_init
	
	
	subroutine Efactor_destructor( this )
	type(Efactor) :: this

!		print*, 'Efactor: destructor'
		if ( this % mode .NE. Efactor_EXACT ) then
			deallocate( this % table )
		end if

	end subroutine Efactor_destructor	
! --
	
	function EfactorGauss_init( sigma, Emax, numPoints ) result(this)
	type(Efactor), pointer :: this
	double precision sigma
	double precision, optional :: Emax
	integer, optional :: numPoints

		this => Efactor( Efactor_GAUSS, sigma, Emax, numPoints )
		
	end function EfactorGauss_init
	
! --
	
	function EfactorLorentz_init( sigma, Emax, numPoints ) result(this)
	type(Efactor), pointer :: this
	double precision sigma
	double precision, optional :: Emax
	integer, optional :: numPoints

		this => Efactor( Efactor_LORENTZ, sigma, Emax, numPoints )
		
	end function EfactorLorentz_init
	
! --
	
	function Efactor_init( mode, sigma, Emax_in, numPoints_in ) result(this)
	type(Efactor), pointer :: this
	integer mode
	double precision sigma
	double precision, optional :: Emax_in
	integer, optional :: numPoints_in

	double precision step, Emax, arg, norm_factor, v
	integer numPoints, i

		allocate(this)
		
		this % mode = mode
		this % sigma = sigma

		if (present(Emax_in)) then
			Emax = Emax_in
		else
			Emax = 20*sigma
		end if
		if (present(numPoints_in)) then
			numPoints = numPoints_in
		else
			numPoints = 65536
		end if

		step = Emax/numPoints
		this % step = step
		this % Emax = Emax
		this % numPoints = numPoints
		
		if (this % mode.NE.Efactor_EXACT) then
			allocate( this % table(0:numPoints) ) ! including last point

			select case(this % mode)
				case(Efactor_GAUSS)
					norm_factor = 1.d0 / (sqrt(2*PI)*sigma)
				case(Efactor_LORENTZ)
					norm_factor = 1.d0 / (PI*sigma)
			end select

			do i=0, numPoints
				v = 0.d0
				arg = ( step*i/sigma ) **2
				select case(this % mode)			
					case(Efactor_GAUSS)
						if (arg.LT.200) v = exp(-arg/2.d0) ! up to exp(-100)
					case(Efactor_LORENTZ)
 						v= 1.d0 / (1.d0 + arg)
 				end select
				this % table(i) = v * norm_factor 
			end do
		end if
		
	end function Efactor_init

!--

	double precision function Efactor_getValue( this, deltaE )
	class( Efactor ) :: this
	double precision ans, deltaE, arg, frac, v1, absDeltaE
	integer pos, posp1
	
		ans = 0.d0
		absDeltaE = abs(deltaE)
		if (this % mode .EQ. Efactor_EXACT) then
			if (absDeltaE.LT.this%sigma) ans = 1.d0
		else
			! Efactor_GAUSS, Efactor_LORENTZ
			if (absDeltaE.LT.this%Emax) then
				arg = absDeltaE / this % step
				pos = int(arg + 0.5d0)
				ans = this % table( pos )
				posp1 = pos+1
				if (posp1.LT.this%numPoints) then
					frac = arg-pos
					v1 = this%table( posp1 )
					ans = ans + (v1-ans)*frac ! linear interpolation
				end if
			end if
		end if
		Efactor_getValue = ans ! this % step

	end function Efactor_getValue

! -- for debug

	subroutine Efactor_show( this, o )
	class( Efactor ) :: this
	class( Output ), pointer :: o

	character*100 info
	
		write(info, '(A,ES25.15)') TRIM(o%sigma) // '=', this % sigma
		select case(this % mode)			
			case(Efactor_GAUSS)
				call o % writeln ( 'Efactor: Gaussian' )
				call o % beginEquation
				call o % formulaItem( info )
				call o % endEquation
			case(Efactor_LORENTZ)
				call o % writeln ( 'Efactor: Lorentzian' )
				call o % beginEquation
				call o % formulaItem( info )
				call o % endEquation
			case(Efactor_EXACT)
				call o % writeln ( 'Efactor: Exact' )
		end select
		

	end subroutine Efactor_show
	


end module Efactor_module


