module MultidimensionalCycle_module

	use Geometry_module
	use Efactor_module
	use QbeTerm_module
	use Perturbation_module
	use Occupations_module

	implicit none

	double precision, parameter :: defaultFactorBoundary = 1d-10

	type :: MDC_CachedListLine ! for each (kroIndex, k1, s_index)
		integer num, d
		! list
		integer, allocatable :: k3i4(:,:,:) ! 0=k/1=i4 x 0:d x position
		double precision, allocatable :: coeff(:) ! posision
	contains
		procedure :: init => MDC_CachedListLine_init
		procedure :: takeSpace => MDC_CachedListLine_takeSpace
		procedure :: dropSpace => MDC_CachedListLine_dropSpace
		procedure :: addEntry => MDC_CachedListLine_addEntry
		procedure :: incSize => MDC_CachedListLine_incSize
		procedure :: getEntry => MDC_CachedListLine_getEntry
		procedure :: getMemoryFootprint => MDC_CachedListLine_getMemoryFootprint

		final :: MDC_CachedListLine_destructor
	end type MDC_CachedListLine

	type :: MDC_CachedList
		integer d, kroVolume, volume, numS
		integer state ! 0=none, 1=counting entries, 2=ready
		type(MDC_CachedListLine), allocatable :: line(:,:,:) ! kroIndex, k1, s_index
	contains
		procedure :: setState  => MDC_CachedList_setState
		procedure :: getMemoryFootprint => MDC_CachedList_getMemoryFootprint
		
		final :: MDC_CachedList_destructor
	end type MDC_CachedList

	type MDC_CachedList_p
		class(MDC_CachedList), pointer :: obj
	end type MDC_CachedList_p






	type MultidimensionalCycle
	! fixed external parameters 
		type(Geometry), pointer :: sizes
		type(Perturbation), pointer :: p		
		type(DispersionRelation_p), allocatable :: dr(:)
		type(Efactor), pointer :: F
		double precision :: factorBoundary
		logical rta
				
	! dimensions
		integer d, volume
	! relation to subsystems operators
		integer, allocatable :: s_occupationIndex(:) ! subsystem index 4 for Occupations
		integer, allocatable :: operatorPosition(:) ! operator position in Perturbation
			! i.e.  0=k0, 1=p0, 2=p1, 3=k1 or otherwise
		
	! run parameters (including openmp version)
		type(Occupations), pointer :: occ, docc ! n_from, dn_to
		integer occupationIndex, deltaIndex1
		integer, allocatable :: replacementPosition(:) ! delta_12 -> replaces 2 with 1

	! work variables
		logical running, firstIndexChanged, accumulatorChanged
		integer :: firstIndex(1:3) ! for dn(...)=
		double precision :: accumulator
		integer earliestIndexChanged ! for chain of calculations
		integer numTerms
		
	! chain of calculations: k0+k1-k2 ... , E0+E1-E2... 
		integer, allocatable :: chainForLastIndex3d(:,:)   ! 1:3, -1:d
		double precision, allocatable :: chainForDE(:)
				
	! registers
		integer, allocatable :: digit(:)
		integer, allocatable :: digit3d(:,:)

	! occupation values for these digit : lazy values
		double precision, allocatable :: n(:)
		integer, allocatable :: kForCalc(:)
		integer, allocatable :: n_key(:)
	! energies for these coordinates
		double precision, allocatable :: e(:)
	! delta E for these coordinates
		double precision :: deltaE
		
	! cached list for calculation with kroneckers

	contains
		procedure :: setFactorBoundary => MDC_setFactorBoundary
		procedure :: start => MDC_start
		procedure :: next  => MDC_next
		procedure :: run   => MDC_run
		procedure :: makeLastIndex  => MDC_makeLastIndex
		procedure :: getVfactor => MDC_getVfactor
		procedure :: useQbeTerm => MDC_useQbeTerm
		procedure :: updateNs => MDC_updateNs
		procedure :: updateAccumulator => MDC_updateAccumulator
		
		procedure :: prepareCacheForKroneckers => MDC_prepareCacheForKroneckers
		procedure :: updateCoeffsForKroneckers => MDC_updateCoeffsForKroneckers
		procedure :: runWithKroneckers   => MDC_runWithKroneckers
		
		procedure :: show  => MDC_show

		final :: MDC_destructor
		
	end type MultidimensionalCycle

	type MultidimensionalCycle_p
		class(MultidimensionalCycle), pointer :: obj
	end type MultidimensionalCycle_p

! constructors

	interface MultidimensionalCycle
		module procedure MDC_init
	end interface MultidimensionalCycle
	interface MDC_CachedList
		module procedure MDC_CachedList_init
	end interface MDC_CachedList

contains

! ===== 1. CachedListLine ============================

	subroutine MDC_CachedListLine_init( this, d )
	class( MDC_CachedListLine ), intent(INOUT) :: this
	integer d
	
		this % num = 0
		this % d = d
	
	end subroutine MDC_CachedListLine_init

	subroutine MDC_CachedListLine_destructor( this )
	type( MDC_CachedListLine ) :: this

		call this % dropSpace
	
	end subroutine MDC_CachedListLine_destructor

	subroutine MDC_CachedListLine_dropSpace( this )
	class( MDC_CachedListLine ), intent(INOUT) :: this
	
		if (this % num .GT. 0) then
			deallocate( this % coeff )
			deallocate( this % k3i4 )
		end if
		this % num = 0
	
	end subroutine MDC_CachedListLine_dropSpace

	subroutine MDC_CachedListLine_takeSpace( this )
	class( MDC_CachedListLine ), intent(INOUT) :: this
	
		if (this % num .GT. 0) then
			allocate( this % coeff( 0:this%num-1 ) )
			allocate( this % k3i4(0:1, 0:this%d,  0:this%num-1 )  )
		end if
		this % num = 0 ! and prepare to fill: num is current position
	
	end subroutine MDC_CachedListLine_takeSpace

	subroutine MDC_CachedListLine_incSize( this )
	class( MDC_CachedListLine ), intent(INOUT) :: this
	
		this % num = this % num + 1
	
	end subroutine MDC_CachedListLine_incSize

	subroutine MDC_CachedListLine_addEntry( this, k, i4, coeff )
	class( MDC_CachedListLine ), intent(INOUT) :: this
	integer :: k(0:this%d)
	integer :: i4(0:this%d)
	double precision coeff
	
	integer pos
	
		pos = this % num
		this % num = pos + 1
		this % k3i4(0,:,pos) = k(:)
		this % k3i4(1,:,pos) = i4(:)
		this % coeff(pos) = coeff		
	
	end subroutine MDC_CachedListLine_addEntry

	subroutine MDC_CachedListLine_getEntry( this, pos, k, i4, coeff )
	class( MDC_CachedListLine ), intent(IN) :: this
	integer pos
	integer, intent(OUT) :: k(0:this%d)
	integer, intent(OUT) :: i4(0:this%d)
	double precision, intent(OUT) :: coeff
	
		k(:)  = this % k3i4(0,:,pos)
		i4(:) = this % k3i4(1,:,pos)
		coeff = this % coeff(pos)		
	
	end subroutine MDC_CachedListLine_getEntry


	function MDC_CachedListLine_getMemoryFootprint( this ) result(mem)
	class( MDC_CachedListLine ), intent(INOUT) :: this
	integer*8 mem
	
		! coeff : real*8
		! k, i4 : integer*4 * (d+1)
		mem = this % num * ( 8 + 2*4*(this%d+1) )
	
	end function MDC_CachedListLine_getMemoryFootprint

! ===== 2. CachedList  , collection of CachedListLines



	function MDC_CachedList_init( d, kroVolume, volume, numS ) result(this)
	integer d, kroVolume, volume, numS 
	type( MDC_CachedList ), pointer :: this
	
		allocate(this)
		
		this % d = d ! number of digits in k(:)
		this % kroVolume = kroVolume
		this % volume = volume
		this % numS = numS

		this % state = 0 ! empty
	
	end function MDC_CachedList_init

	subroutine MDC_CachedList_destructor( this )
	type( MDC_CachedList ) this

		if (this % state.GT.0) then
			call this % setState( 0 )
		end if

	end subroutine MDC_CachedList_destructor

	subroutine MDC_CachedList_setState( this, state )
	class( MDC_CachedList ), intent(INOUT) :: this
	integer state
	
	integer ss, k1, kroIndex
	
		select case(state)
			case(0)
				if (this % state .EQ. 2) then
					call this % setState( 1 )
				end if
				if (this % state .GE. 1) then
					deallocate( this % line )
				end if
			case(1)
				if (this % state .EQ.0) then
					allocate( this % line( 0:this % kroVolume, 0:this%volume-1, 0:this % numS-1 ) )
					do ss=0, this%numS-1
					do k1=0, this % volume-1
					do kroIndex=0, this % kroVolume
						call this % line( kroIndex, k1, ss ) % init( this % d ) ! now line is empty
					end do
					end do
					end do
				end if
				if (this % state .EQ. 2) then
					do ss=this % numS-1, 0, -1
					do k1=this % volume-1, 0, -1
					do kroIndex=this % kroVolume, 0, -1
						call this % line( kroIndex, k1, ss) % dropSpace() ! if needed (num .GT. 0)
					end do
					end do
					end do
				end if
			case(2)
				if (this % state .EQ.0) then
					call this % setState( 1 )
				end if
				if (this % state .EQ. 1) then
					do ss=0, this%numS-1
					do k1=0, this % volume-1
					do kroIndex=0, this % kroVolume
						call this % line( kroIndex, k1, ss) % takeSpace() ! if needed (num .GT. 0)
					end do
					end do
					end do
				end if
		end select		
		this % state = state	
	
	end subroutine MDC_CachedList_setState



	function MDC_CachedList_getMemoryFootprint( this ) result(mem)
	class( MDC_CachedList ), intent(INOUT) :: this
	
	integer ss, k1, kroIndex

	integer*8 mem
	
!		state=0/1 ! only main array and empty objects
		mem = 16 * this % kroVolume * this%volume-1 * this % numS
		if (this%state.EQ.2) then ! objects with data
			do ss=0, this%numS-1
			do k1=0, this % volume-1
			do kroIndex=1, this % kroVolume
				mem = mem + this % line( kroIndex, k1, ss) % getMemoryFootprint()
			end do
			end do
			end do
		end if
	end function MDC_CachedList_getMemoryFootprint
	
! ===== 3. MultidimensionalCycle


	function MDC_init( sizes, p, F, factorBoundary ) result(this)

	type(Geometry), pointer :: sizes
	type(Perturbation), pointer :: p		
	type( Efactor ), pointer :: F
	double precision, optional :: factorBoundary
	
	type(MultidimensionalCycle), pointer :: this

	integer d

		allocate( this )

		call p % prepare

		this % F => F
		
		this % sizes => sizes
		this % p => p
		
		d = p % numAllOperators-1 ! last index depends on previous 

		this % d = d
		this % volume = sizes % volume

		this % factorBoundary = defaultFactorBoundary
		if (present(factorBoundary)) then
			this % factorBoundary = factorBoundary
		end if

		allocate( this % replacementPosition(0:d) )
		allocate( this % digit(0:d) )
		allocate( this % digit3d(1:3, 0:d) )
		allocate( this % n_key(0:d) )
		allocate( this % n(0:d) )
		allocate( this % kForCalc(0:d) )
		this % n_key(:) = -1
		allocate( this % e(0:d) )
		allocate( this % operatorPosition(0:d) ) ! 0:k0, 1:p0, 2:p1, 3:k1 or otherwise
		allocate( this % dr(0:d) )
		allocate( this % s_occupationIndex(0:d) )

		allocate( this%chainForLastIndex3d(1:3, -1:d ) )
		allocate( this%chainForDE(-1:d) )
		call this%sizes%zero( this%chainForLastIndex3d(:, -1 ) )
		this%chainForDE(-1) = 0.d0

		this % earliestIndexChanged = 0
		this % accumulatorChanged = .FALSE.

	end function MDC_init



	subroutine MDC_destructor( this )
	type( MultidimensionalCycle ) :: this

!		print*, ' MDC: destructor'

		deallocate( this%chainForDE )
		deallocate( this%chainForLastIndex3d )

		deallocate( this % s_occupationIndex )
		deallocate( this % dr )
		deallocate( this % operatorPosition )
		deallocate( this % e )
		deallocate( this % n_key )
		deallocate( this % digit3d )
		deallocate( this % digit )
		
	end subroutine MDC_destructor


	subroutine MDC_setFactorBoundary( this, fb )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	double precision fb
	
		this % factorBoundary = fb

	end subroutine MDC_setFactorBoundary	

!----------- iterator functions

	function MDC_start( this, occupationIndex, openMpStart1 ) result(not_empty)

	use omp_lib
	
	class( MultidimensionalCycle ), intent(INOUT) :: this
	integer occupationIndex
	integer, optional :: openMpStart1
	
	type(Subsystem), pointer :: s
	integer i, input_s_index, current_occupationIndex, sw, startIndex1, rep
	integer lastIndexCandidate
	logical not_empty

		this % occupationIndex = occupationIndex

! init s_occupationIndex(i) , moving input subsystem occupationIndex to 1st
		input_s_index = -1
		do i=0, this % d
			this % operatorPosition(i) = i

			s =>  this % p % info(i) % subsystem

			current_occupationIndex = s % occupationIndex
			this % dr(i) % obj => s % dispersionLaw
			this % s_occupationIndex(i) = s % occupationIndex

			if ((current_occupationIndex .EQ. occupationIndex ).AND.(input_s_index.LT.0)) then 
				! first time to encounter this subsystem
				input_s_index = i
			end if
		end do
		
		not_empty = .FALSE.
!		if (input_s_index.LT.0) stop 'mdc_start: occupationIndex not found'
		if (input_s_index.GE.0) then
			not_empty = .TRUE.

			lastIndexCandidate = 1
			! current index not 0 -> move to 0
			if (input_s_index.GT.0) then
				sw = this % s_occupationIndex(0)
				this % s_occupationIndex(0) = this % s_occupationIndex(input_s_index)
				this % s_occupationIndex(input_s_index) = sw

				sw = this % operatorPosition(0)
				this % operatorPosition(0) = this % operatorPosition(input_s_index)
				this % operatorPosition(input_s_index) = sw

				lastIndexCandidate = 0
			end if

			startIndex1 = this % volume
			this % deltaIndex1 = 1
			if (present(openMpStart1)) then
				startIndex1 = this % volume - openMpStart1	
				this % deltaIndex1 = omp_get_num_threads()
			end if

			! start calculation
			this % digit(:) = this % volume - 1 ! so that first next() gives zeros	
			this % digit3D(:,:) = 777 ! bad value to be replaced
			this % digit(0) = startIndex1 - 1 ! start point using openmp
			this % running = .FALSE. ! stop when zeros appear second time
			this % firstIndexChanged = .TRUE.
			this % accumulator = -666.d0 ! for debug
		
			this % n_key(:) = -1
		end if

	end function MDC_start

!--

	logical function MDC_next( this )
	class( MultidimensionalCycle ), intent(INOUT) :: this

	logical ans, carry
	integer i, j, rep, delta
	
		ans = .TRUE.
		
		this % firstIndexChanged = .FALSE.
		carry = .TRUE.
		i = this%d - 1
		do while ( (i.GE.0 ).AND. carry )
			if (carry) then
				j = this % operatorPosition(i)

				delta = 1
				if (i.EQ.0) delta = this % deltaIndex1
				
				if (this%replacementPosition(j).LT.0) then
					this % digit(i) = this % digit(i) + delta
					carry = (this % digit(i) .GE. this % volume)
					if (carry) this % digit(i) = this % digit(i) - this % volume
					call this % sizes % toVector( this % digit(i), this % digit3d(:,i) )
					this % e(i) = this % dr( j ) % obj % getEnergy( this % digit3d(:,i), this%sizes )
				end if
				if (i.EQ.0) this % firstIndexChanged = .TRUE.
			end if
			i = i - 1
		end do
		this % earliestIndexChanged = max(i,0)
		ans = .NOT.(carry.AND.this%running) ! first time zeros are not the end

		! use replacement positions! excluding last index (not calculated yet)
		do i=1, this%d-1
			j = this % operatorPosition(i)
			rep = this%replacementPosition(j)
			if (rep.GE.0) then
				this % digit(this%d) = 666
				this % digit3d(:,this%d) = 666
				this % digit(i) = this % digit(rep)
				this % digit3d(:,i) = this % digit3d(:,rep) ! delta_12 : copy 1 to 2
				this % e(i) = this % dr( j ) % obj % getEnergy( this % digit3d(:,i), this%sizes )
			end if
		end do
		if (this % firstIndexChanged) then
			if (this%running .AND. this % accumulatorChanged) then
				call this % docc % incValue(this%firstIndex(:), &
&					this % s_occupationIndex(0), this % accumulator)
				this % accumulator = 0.d0
				this % accumulatorChanged = .FALSE.
			end if
			if (ans) then
				this % firstIndex(:) = this % digit3d(:,0)
				this % accumulator = 0.d0
				this % accumulatorChanged = .FALSE.
			end if
		end if

		this % running = ans ! while ans = .TRUE.
				
		MDC_next = ans
		
	end function MDC_next



! calculate last index (k4=k1+k2-k3 etc.)
! and delta E
	logical function MDC_makeLastIndex(this)
	class( MultidimensionalCycle ), intent(INOUT) :: this
	logical good
	
	integer i, j, istart, rep
	logical minus

			istart = this % earliestIndexChanged

			! update chain of calculations
			do i=istart, this%d-1
				j = this % operatorPosition(i)
				minus = (.NOT.this % p % info(j) % dagger)

				if (minus) then
					call this % sizes % minus1(this % chainForLastIndex3D(:,i-1), &
&						this % digit3d(:,i), this % chainForLastIndex3D(:,i) )
					this % chainForDE(i) = this % chainForDE(i-1) - this % e(i)
				else
					call this % sizes % plus1(this % chainForLastIndex3D(:,i-1), &
&						this % digit3d(:,i), this % chainForLastIndex3D(:,i) )
					this % chainForDE(i) = this % chainForDE(i-1) + this % e(i)
				end if
			end do

			! last index almost ready!
			! invert index for a(+)
			i = this%d
			j = this % operatorPosition(i)
			this % chainForLastIndex3D(:,i) = this % chainForLastIndex3D(:,i-1)
			minus = this % p % info( j ) % dagger
			if (minus) then
				call this % sizes % minus( this % chainForLastIndex3D(:,i) )
			end if


			this % digit3d(:, i) = this % chainForLastIndex3D(:,i)
			this % digit( i ) = this % sizes % fromVector( this % digit3d(:, i) )
			this % e(i) = this % dr( j ) % obj % getEnergy( this % digit3d(:,i), this%sizes )
			if (.NOT.minus) then
				this % chainForDE(i) = this % chainForDE(i-1) - this % e(i)
			else
				this % chainForDE(i) = this % chainForDE(i-1) + this % e(i)
			end if

			this % deltaE = this % chainForDE(i)

! check if delta_34 is OK (if last digit is restricted by Kronecker) ! TODO: need to research more (mdc+abcd+cache)
			good = .TRUE.
			rep = this%replacementPosition(j) ! TODO: check (i) vs (j) when delta_34
			if (rep.GE.0) then
				good = (this%digit(i) .EQ. this%digit(rep) )
			end if
			
			MDC_makeLastIndex = good

	end function MDC_makeLastIndex

!--------------------

	subroutine MDC_run( this, occ, qt, occupationIndex, rta, docc, openMpStart1 )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermAsIs), pointer :: qt
	type(Occupations), pointer :: occ, docc
	integer occupationIndex
	logical rta
	integer, optional :: openMpStart1
	
		this % occ => occ		! work data for this run
		this % docc => docc
		this % rta = rta
		
		this % replacementPosition(:) = -1 ! nothing to replace

		if (this % start( occupationIndex, openMpStart1 )) then ! start iterator, prepare work data

			this % numTerms = 0
			do while (this % next() ) 
				call this % updateAccumulator( qt )
			end do
		end if
		this % rta = .FALSE.
		
	end subroutine MDC_run


	subroutine MDC_prepareCacheForKroneckers( this, qt, &
&		occupationIndex, rta, cachedList, iter, openMpStart1 )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermSumOfABCD), pointer :: qt
	class(MDC_CachedList), pointer :: cachedList
	
	integer occupationIndex, iter ! iter: 1=count, 2=fill
	logical rta
	integer, optional :: openMpStart1

	integer kroIndex, m, minKroIndex, num, maxKroIndex
	
		this % rta = rta
		m = qt % partMultiplier( occupationIndex )

		minKroIndex = 1
		if (rta) minKroIndex = 0
		maxKroIndex = this % p % numDeltaCombinations - 1
		do kroIndex=minKroIndex, maxKroIndex
			this % replacementPosition(:) = -1
			if (kroIndex.GT.0) this % replacementPosition(:) = qt % replacementPosition(:, kroIndex)
			if (this % start( occupationIndex, openMpStart1 )) then
				do while (this % next() ) 
					call this % updateCoeffsForKroneckers( qt, &
&						kroIndex, cachedList, iter, m )
				end do ! iterator
			end if
		end do ! kroIndex
		this % rta = .FALSE.

	end subroutine MDC_prepareCacheForKroneckers


	subroutine MDC_runWithKroneckers( this, occ, qt, occupationIndex, &
&		rta, docc, cachedList, openMpStart1 )

	use omp_lib
	
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermSumOfABCD), pointer :: qt
	
	type(Occupations), pointer :: occ, docc
	integer occupationIndex
	logical rta
	class(MDC_CachedList), pointer :: cachedList
	integer, optional :: openMpStart1

	integer kroIndex, maxKroIndex, k1, pos, num, i, minKroIndex
	double precision accumulator, coeff, n_factor
	class(QbeTermAsIs), pointer :: qt_asis ! to change type of pointer
	integer k1start, deltak1
	type(Subsystem), pointer :: s
	
		this % occ => occ		! work data for this run
		this % docc => docc
		this % rta = rta

		qt % currentSubsystemIndex = occupationIndex
		qt_asis => qt % QBeTermAsIs
		maxKroIndex = qt % kroVolume - 1
		do k1=0, this%d
			this % operatorPosition(k1) = k1
		end do

		minKroIndex = 1
		if (rta) minKroIndex = 0
	
		k1start = 0
		deltak1 = 1
		if (present(openMpStart1)) then
			k1start = openMpStart1
			deltak1 = omp_get_num_threads()
		end if
		do k1=k1start, this % volume-1, deltak1
			accumulator = 0.d0
			do kroIndex=minKroIndex, maxKroIndex

				num = cachedList % line(kroIndex,k1,occupationIndex) % num
				do pos=0, num-1
					call cachedList % line(kroIndex,k1,occupationIndex) % getEntry( pos, &
&						this % digit( : ), &
&						this % s_occupationIndex(:), &
&						coeff )
					do i=0, this%d
						call this % sizes % toVector( this % digit(i), this % digit3d(:,i) )
					end do
					call this % updateNs( .TRUE. )
					if (rta) then
						n_factor = qt % addValueForRTAForKronecker( this%n, this%kForCalc, kroIndex, occupationIndex )
					else
						n_factor = qt % getValueForKronecker( this%n, kroIndex, -1 )
					end if
					accumulator = accumulator + n_factor * coeff
				end do
			end do
			call this % sizes % toVector( k1, this % digit3d(:,0) )
			call docc % incValue( this % digit3d(:,0), occupationIndex, accumulator )
		end do

	end subroutine MDC_runWithKroneckers


	subroutine MDC_updateCoeffsForKroneckers( this, qt, &
&		kroIndex, cachedList, iter, partMultiplier )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermSumOfABCD), pointer :: qt
	integer kroIndex
	class(MDC_CachedList), pointer :: cachedList
	integer iter, partMultiplier
		
	double precision F, V_factor, coeff
	class(QbeTermAsIs), pointer :: qt_asis ! to change type of pointer

	integer b, c

		if ( this % makeLastIndex() ) then

			if (qt % calculateKroneckerFactor( this%digit( this%operatorPosition(:)), kroIndex ).NE.0) then
			qt_asis => qt % QBeTermAsIs

				F = this % F % getValue( this%deltaE )
				V_factor = this % getVfactor( qt_asis )
				coeff = F * (V_factor**2)
				if ( coeff .GT. this % factorBoundary ) then
					coeff = coeff * partMultiplier
					b = this % digit(0)
					c = this % occupationIndex
					if (iter.EQ.1) then
						call cachedList % line( kroIndex, b, c ) % incSize
					else
						call cachedList % line( kroIndex, b, c ) % addEntry( &
&							this % digit( this % operatorPosition(:) ), &
&							this % s_occupationIndex( this % operatorPosition(:) ), &
&							coeff )
					end if ! iter=1/2
				end if ! coeff .GT. factorBoundary
			else
				print*, ' ABCD: invalid kroFactor: kroIndex=', kroIndex
				print*, '   k=', this%digit( this%operatorPosition(:))
				print*, '   p%deltaPos=', this%p%deltaPositions
				stop
			end if ! valid kroneckerFactor
		end if ! makeLastIndex

	end subroutine MDC_updateCoeffsForKroneckers

!--
	subroutine MDC_updateAccumulator( this, qt )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermAsIs), pointer :: qt
	
	double precision F, n_factor, V_factor
	
		if ( this % makeLastIndex() ) then
			F = this % F % getValue( this%deltaE )
			if ( F .GT. this % factorBoundary ) then
				n_factor = this % useQbeTerm( qt )
				V_factor = this % getVfactor( qt )
				this % accumulator = this % accumulator + F * n_factor * (V_factor ** 2)
				this % accumulatorChanged  = .TRUE.
				this % numTerms = this % numTerms + 1
			end if
		end if
	end subroutine MDC_updateAccumulator



	double precision function MDC_getVfactor( this, qt )
	class( MultidimensionalCycle ) :: this
	class(QbeTermAsIs), pointer :: qt
	
	integer q(1:3), posForQ(0:1), s_index, j

		if (this%p%mode .EQ. Q_DEPENDENT) then
			s_index = qt % currentSubsystemIndex

			posForQ(:) = this % p % operatorPositionForQ(:,0)
			if ( (s_index .EQ. posForQ(0)) .OR. (s_index .EQ. posForQ(1))) then
				posForQ(:) = this % p % operatorPositionForQ(:,1)
			end if

			q(:) = 0
			if (posForQ(0).GE.0) then
				j = this % operatorPosition( posForQ(0) )
				call this % sizes % add1( q, this%digit3d( :, j ) )
			end if
			if (posForQ(1).GE.0) then
				j = this % operatorPosition( posForQ(1) )
				call this % sizes % sub1( q, this%digit3d( :, j ) )
			end if
			MDC_getVfactor = this % p % amplitudeTable % getValue( q )
		else
			MDC_getVfactor = this % p % amplitudeValue
		end if
		
	end function MDC_getVfactor
	
!--

	subroutine MDC_updateNs( this, updateAnyways )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	logical, optional :: updateAnyways
	
	integer i, i4, opPos
	double precision n

		do i=0, this%d
			if ((present(updateAnyways)) .OR. (this % n_key(i) .NE. this % digit(i) )) then
				i4 = this % s_occupationIndex(i)
				n = this % occ % getValue( this % digit3d(:,i), i4 )
				opPos = this % operatorPosition(i)
				this % n( opPos ) = n
				this % kForCalc( opPos ) = this % digit( i )
				this % n_key(i) = this % digit(i)
			end if
		end do

	end subroutine MDC_updateNs
	
	
	
	double precision function MDC_useQbeTerm( this, qt )
	class( MultidimensionalCycle ), intent(INOUT) :: this
	class(QbeTermAsIs), pointer :: qt

		call this % updateNs()
		if (this % rta) then
			MDC_useQbeTerm = qt % addValueForRTA( this%n, this%kForCalc, this % occupationIndex )
		else
			MDC_useQbeTerm = qt % getValue( this%n, this%kForCalc, -1 )
		end if
	end function MDC_useQbeTerm
	

!---

	subroutine MDC_show( this )
	class( MultidimensionalCycle ) :: this
	integer i, j
	
		write(*, '(A)', advance='no') ' k[]=('
		do i=0, this%d
			j = this % operatorPosition(i)
			write(*, '(A,I1,A)', advance='no') ' ', this%s_occupationIndex(j),':'
			call Geometry_print( this%digit3d(:,j) )
		end do
		print*, '), de=', this%deltaE
		
	end subroutine MDC_show

end module MultidimensionalCycle_module


