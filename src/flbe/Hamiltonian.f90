module Hamiltonian_module

	use Output_module
	use Subsystem_module
	use Perturbation_module

	implicit none

	integer, parameter :: max_num_perturbations = 10

	type Hamiltonian
		type(Subsystem_p), dimension(0:max_num_subsystems-1) :: subsystems
		type(Perturbation_p), dimension(0:max_num_perturbations-1) :: perturbations
		integer numSubsystems
		integer numPerturbations
		logical prepared

	contains
	
		procedure :: addSubsystem => Hamiltonian_addSubsystem
		procedure :: addPerturbation => Hamiltonian_addPerturbation
		procedure :: prepare => Hamiltonian_prepare
		procedure :: dropPreparations => Hamiltonian_dropPreparations

		procedure :: toString => Hamiltonian_toString
		procedure :: show => Hamiltonian_show

		final :: Hamiltonian_destructor
		
	end type Hamiltonian

! constructors

	interface Hamiltonian
		module procedure Hamiltonian_init
	end interface Hamiltonian


contains

! constructors: bodies

	function Hamiltonian_init() result(this)
	type(Hamiltonian), pointer :: this

		allocate( this )

		this % numSubsystems = 0
		this % numPerturbations = 0
		this % prepared = .FALSE.

	end function Hamiltonian_init

	subroutine Hamiltonian_destructor( this )
	type( Hamiltonian ) :: this
	
	integer i

!		print*, ' H.destructor -1-'
		! unregister subsystems for later use in future calculations
		do i=0, this % numSubsystems-1
!			print*, ' H.destructor -2- i=', i
			this % subsystems( i ) % obj % occupationIndex = -1
		end do
!		print*, ' H.destructor -3-'
	end subroutine Hamiltonian_destructor


!-- setters --

	subroutine Hamiltonian_addSubsystem( this, s )
	class(Hamiltonian), intent(INOUT) :: this
	type( Subsystem ), pointer :: s
	integer pos

		if (s % occupationIndex .LT. 0) then ! not registered yet for calculation
			pos = this % numSubsystems
			call s % setOccupationIndex( pos )
			if (pos .GE. max_num_subsystems) then
				stop 'Hamiltonian_addSubsystem : max_num_subsystems exceeded'
			end if
			this % subsystems(pos) % obj => s
			this % numSubsystems = pos + 1
		end if
		call this % dropPreparations

	end subroutine Hamiltonian_addSubsystem

!--

	subroutine Hamiltonian_addPerturbation( this, p )
	class(Hamiltonian), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	integer pos, i
	
		do i=0, p % numSubsystems-1
			call this % addSubsystem( p % subsystems(i) % obj )
		end do

		pos = this % numPerturbations
		if (pos .GE. max_num_perturbations) then
			stop 'Hamiltonian_addPerturbation : max_num_perturbations exceeded'
		end if
		
		this % perturbations(pos) % obj => p
		this % numPerturbations = pos + 1
		call this % dropPreparations
		
	end subroutine Hamiltonian_addPerturbation


	subroutine Hamiltonian_prepare( this )
	class(Hamiltonian), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	integer i

		if ( .NOT. this % prepared ) then
			do i=0, this % numPerturbations-1
				p => this % perturbations(i) % obj
				call p % prepare
			end do
			this % prepared = .TRUE.
		end if

	end subroutine Hamiltonian_prepare
		
	subroutine Hamiltonian_dropPreparations( this )
	class(Hamiltonian), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	integer i

		if ( this % prepared ) then
			do i=0, this % numPerturbations-1
				p => this % perturbations(i) % obj
				call p % dropPreparations
			end do
			this % prepared = .FALSE.
		end if

	end subroutine Hamiltonian_dropPreparations
		
! toString() useful for information and debug

	character*200 function Hamiltonian_toString( this )
	class(Hamiltonian) :: this

	character*100 s1, s2
	
		write(s1, '(I5)') this % numSubsystems
		write(s2, '(I5)') this % numPerturbations

		Hamiltonian_toString = 'Hamiltonian: num.subsystems=' // TRIM(s1) // ' np=' // TRIM(s2)

	end function Hamiltonian_toString

!--

	subroutine Hamiltonian_show(this, o)
    class( Hamiltonian ) :: this
	class( Output ), pointer :: o

    integer i, cnt, lineLength
	character*10 hhat
	character*30 k, kbold
	    
    type(Subsystem), pointer :: s
    type(Perturbation), pointer :: p
	
	call o % write( 'Hamiltonian:' )

	call this % prepare
	
	k = 'k'
	call o % bold( k, kbold )
	call o % hat( 'H', hhat )
	
	call o % write( 'Single-particle energies:' )
    call o % beginEquation()
    do i=0, this%numSubsystems-1
		if (i.GT.0) call o % equationNewLine()
		s => this % subsystems(i) % obj
		if (this%numSubsystems.EQ.1) then
			call s % dispersionLaw % show( '', o )
		else
			call s % dispersionLaw % show( '(' // s % operatorLetter // ')', o )
		end if
	end do    
	call o % endEquation()

	call o % write( 'Perturbations:' )
    do i=0, this%numPerturbations-1
	    call o % beginEquation()
		p => this % perturbations(i) % obj
		call o % formulaItem( hhat )
		call o % upperLowerIndex( '', TRIM(p % title)  )
		call o % formulaItem( '=' )
		call p % show( o )
		call o % endEquation()
	end do
    
    
	end subroutine Hamiltonian_show
!--

end module Hamiltonian_module


