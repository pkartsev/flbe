	subroutine QTC_FFT_destructor( this )
	type( QbeTermsCalculatorFFT ) :: this

	integer pp, half, rtaIndex
	
!		print*, 'QBE_FFT destructor'

		if (this%prepared) then ! more arrays to release

			call fftw_destroy_plan( this % plan3dManyBack )
			call fftw_destroy_plan( this % plan3dMany )

			call fftw_destroy_plan( this % plan4dBackC2R )
			call fftw_destroy_plan( this % plan4dR2C )

			call fftw_destroy_plan( this % plan1dBack )
			call fftw_destroy_plan( this % plan1d )

			deallocate( this % calc_k_many )
			deallocate( this % calc_r_many )
			deallocate( this % calc_keR )
			deallocate( this % calc_rgR )
			deallocate( this % n_rgR )
			deallocate( this % n_keR )
			deallocate( this % S_rgR )

			deallocate( this % emin )

			deallocate( this % lineshape_g )
			deallocate( this % lineshape_e )
			
			this % prepared = .FALSE.

!			call fftw_cleanup()
!			call fftw_cleanup_threads()
		end if
		
		do pp=0, this % numP - 1
			do half=1,2
				do rtaIndex=0,1
					if (associated( this % mdc_cache(half,pp,rtaIndex) % obj) ) then
						deallocate( this % mdc_cache(half,pp,rtaIndex) % obj)
					end if
				end do
			end do
		end do
		
		deallocate( this % mdc_cache )

!		print*, 'QTC_Simple destructor'

		deallocate( this % term )
		deallocate( this % minusR )

	end subroutine QTC_FFT_destructor


	function QTC_FFT_init( sizes, hamilt, F, useKroneckers, de ) result(this)

	use omp_lib
	use fftw3
	
	type( QbeTermsCalculatorFFT ), pointer :: this
	type(Geometry), pointer :: sizes
	type(Hamiltonian), pointer :: hamilt
	type(Efactor), pointer :: F
	logical useKroneckers
	double precision :: de

		allocate( this )
		
		call this % initFromGHFK( sizes, hamilt, F, useKroneckers )
		call this % init2( de )
		
	end function QTC_FFT_init


	subroutine QTC_FFT_init2( this, de )
	
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this
	double precision de ! energy axis discretization
	
	type( OperatorInfo ), pointer :: info
	type( Perturbation ), pointer :: p		
	type( Subsystem ), pointer :: s
	type( DispersionRelation ), pointer :: dr
	class( QbeTermFFT ), pointer :: qt
	integer pp, ss, pos, LX, LY, LZ, i
	logical as_current, current_dagger
	
	double precision emax, emin, sumDE
	integer numE, half
	double precision sumDE_C, sumDE_A ! creation/annihilation
		

		allocate( this % minusR(0:max_num_operators-1, 1:2, 0 : this % numS - 1, 0 : this % numP - 1) )
		allocate( this % term(1:2, 0 : this%numP - 1) )
		allocate( this % mdc_cache( 1:2, 0 : this%numP - 1, 0:1) ) ! half, p, rta
		
		
		sumDE = 0.d0

		this % maxk1multiplier = 1
		do pp=0, this % numP - 1
			nullify( this % mdc_cache( 1,pp,0 ) % obj ) ! so that can we check if associated()
			nullify( this % mdc_cache( 2,pp,0 ) % obj )	! in destructor
			nullify( this % mdc_cache( 1,pp,1 ) % obj )
			nullify( this % mdc_cache( 2,pp,1 ) % obj )
			p => this % hamilt % perturbations(pp) % obj
			call p % prepare
			this % maxk1multiplier = max( this % maxk1multiplier, p % maxNumEqualMomenta )

			sumDE_C = 0.d0 ! actually maxDE
			sumDE_A = 0.d0
			do ss=0, p%numAllOperators-1
				info => p % info(ss)
				s => info % subsystem
				dr => s % dispersionLaw
				emax = dr % getMaxEnergy( this % sizes )
				emin = dr % getMinEnergy( this % sizes )
!				print*, ' ss=', ss, ' emin/max=', emin, emax
				if (info % dagger) then
					sumDE_C = max( sumDE_C, (de + emax - emin) * p % numOperators( CREATION, ss ))
				else
					sumDE_A = max( sumDE_A, (de + emax - emin) * p % numOperators( ANNIHILATION, ss ))
				end if
			end do
			sumDE_C = sumDE_C * p % numSubsystems
			sumDE_A = sumDE_A * p % numSubsystems
			if (sumDE_C.GT.sumDE) sumDE = sumDE_C
			if (sumDE_A.GT.sumDE) sumDE = sumDE_A

			do half=1, this%numHalves(pp)
				if (half.EQ.2) p => p % hc
				qt => QbeTermFFT( p, this % useKroneckers )
				this % term(half,pp) % obj => qt
				do ss=0, this % numS - 1 
					s => this % hamilt % subsystems( ss ) % obj
					pos = p % indexOf( s )
					if (this % term_ps(pp,ss)) then
					    current_dagger = qt % partDagger( pos )
						do i=0, qt % order - 1 ! operators in this perturbation
							as_current = qt % dagger(i).EQV.current_dagger
							this % minusR(i,half,ss,pp) = as_current
						end do
					end if ! term present for p,s
				end do ! ss


			end do ! half
		end do ! pp
!		print*, ' sumDE=', sumDE

		! FFT specific stuff
		
		! max.k = L/2 -> sum(k^2) ; 0...Emax -> twice and little more
!		LX = this % sizes % L(1)
!		LY = this % sizes % L(2)
!		LZ = this % sizes % L(3)
!		sumEmax = 2*( (LX*LX + LY*LY + LZ*LZ)/4) 
!		
!		print*, 'old numE=', numE		

!		numE = 4 + int(sumDE/de + 0.5d0)
		numE = 0 + int(sumDE/de + 0.5d0)
!		print*, 'numE=', numE	, ' de=', de
!stop
		if (numE.GT.1000000000) then
			print*, 'numE=', numE, ' de=', de, ' sumDE=', sumDE
			stop 'probably bad DE'
		end if
		
		this % de = de
		this % numE = numE
		

		allocate( this % emin(0 : this%numS - 1) )
		do ss=0, this%numS-1
			s => this % hamilt % subsystems(ss) % obj
			dr => s % dispersionLaw
			emin = dr % getMinEnergy( this % sizes )
			this % emin( ss ) = emin
		end do

		this % prepared = .FALSE.
		this % cachePrepared(:) = .FALSE.

		call this % latePreparation ! can be omitted to speed up initialization; is being called again later
		
	end subroutine QTC_FFT_init2


	subroutine QTC_FFT_latePreparation( this )
	implicit none
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

		if (.NOT. this%prepared) then
			call this % prepareWorkArrays
			call this % prepareFftPlans
			call this % prepareLineshape
			call this % make_S_rg
			call this % prepareMoreArrays

			this % prepared = .TRUE.
		end if
		
	end subroutine QTC_FFT_latePreparation
		
	subroutine QTC_FFT_prepareWorkArrays( this )
	implicit none
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

	integer LX, LY, LZ, numE, MB
	integer*8 mem
	
		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

		allocate( this % lineshape_e( 0:numE-1 ) ) 
		allocate( this % lineshape_g( 0:numE-1 ) ) 

		this % volume3D = LX*LY*LZ
		this % volume4D = LX*LY*LZ*numE

		mem = this % volume4D * (1 + 2 + 2 * this % numS)  ! (2*NS+3+1) of 4D arrays of real*8
		! complex*16 volume/2: *_rg = 3 arrays (1+NS+NS), 
		! real*8 volume*1: *_ke = 2 arrays (+ 1 temp)
		MB = int( (8.d0*mem)/1048576 + 0.5d0)
		if (MB.GT.8000) then
			write(*, '(A,I0,A)') 'memory size for FFT: ', MB, 'MB'
		end if

		allocate( this % calc_rgR( 0:LX/2, 0:LY-1, 0:LZ-1, 0:numE-1 ) )
		allocate( this % S_rgR( 0:LX/2, 0:LY-1, 0:LZ-1, 0:numE-1, 0:this % numS - 1 ) ) 

		allocate( this % n_keR( 0:LX-1, 0:LY-1, 0:LZ-1, 0:numE-1 ) ) 
		allocate( this % calc_keR( 0:LX-1, 0:LY-1, 0:LZ-1, 0:numE-1 ) )

		allocate( this % calc_r_many( 0:LX-1, 0:LY-1, 0:LZ-1, 0:numE-1 ) )
		allocate( this % calc_k_many( 0:LX-1, 0:LY-1, 0:LZ-1, 0:numE-1 ) )
		
	end subroutine QTC_FFT_prepareWorkArrays



	subroutine QTC_FFT_prepareMoreArrays( this )
	implicit none
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

	integer LX, LY, LZ, numE

		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

		! many S : different E(k) for each subsystem 
		allocate( this % n_rgR( 0:LX/2, 0:LY-1, 0:LZ-1, 0:numE-1,  0:this % numS - 1 ) )

	end subroutine QTC_FFT_prepareMoreArrays


	subroutine QTC_FFT_prepareFftPlans( this )
	use omp_lib
	implicit none
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

	integer dims(0:3) ! 4D
	integer dims3d(0:2), inembed(0:2), onembed(0:2), howmany, idist, odist, istride, ostride ! 3D x numE
	integer plannerFlags
	logical oldWisdom
	
!	integer void
!		void = fftw_init_threads()
!		if (void.EQ.0) stop 'fail'
		call fftw_plan_with_nthreads( this % numT )
		dims(3) = this % sizes % L(1)
		dims(2) = this % sizes % L(2)
		dims(1) = this % sizes % L(3)
		dims(0) = this % numE
		plannerFlags = FFTW_ESTIMATE
!		plannerFlags = FFTW_MEASURE
!		plannerFlags = FFTW_PATIENT + FFTW_DESTROY_INPUT
	
		oldWisdom = QTC_FFT_loadWisdom( dims, this % numT )
		if (oldWisdom) then
			plannerFlags = FFTW_ESTIMATE
		end if

		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 1       ' // CHAR(13)
		this % plan1d = fftw_plan_dft_1d( this%numE, this%lineshape_e, this%lineshape_g, FFTW_FORWARD, plannerFlags)
		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 2       ' // CHAR(13)
		this % plan1dBack = fftw_plan_dft_1d( this%numE, this%lineshape_g, this%lineshape_e, FFTW_BACKWARD, plannerFlags)

		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 3       ' // CHAR(13)
		this % plan4dR2C = fftw_plan_dft_r2c(4,dims, this % calc_keR, this % calc_rgR, plannerFlags)
		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 4       ' // CHAR(13)
		this % plan4dBackC2R = fftw_plan_dft_c2r(4,dims, this % calc_rgR, this % calc_keR, plannerFlags)

		! Advanced interface: numE 3D arrays
		dims3d(2) = this % sizes % L(1)
		dims3d(1) = this % sizes % L(2)
		dims3d(0) = this % sizes % L(3)
		idist = dims3d(0)*dims3d(1)*dims3d(2)
		odist = idist
		howmany = this % numE 
		inembed(:) = dims3d(:)
		onembed(:) = dims3d(:)
		istride = 1
		ostride = 1
		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 5       ' // CHAR(13)
		this % plan3dMany = &
&			fftw_plan_many_dft(3, dims3d, howmany, &
&			this % calc_k_many, inembed, istride, idist, &
&			this % calc_r_many, onembed, ostride, odist, &
&			FFTW_FORWARD, plannerFlags)
		write(*, '(A)', advance='no') CHAR(13) // 'FFT plan 6       ' // CHAR(13)
		this % plan3dManyBack = &
&			fftw_plan_many_dft(3, dims3d, howmany, &
&			this % calc_r_many, onembed, ostride, odist, &
&			this % calc_k_many, inembed, istride, idist, &
&			FFTW_BACKWARD, plannerFlags)
		write(*, '(A)', advance='no') CHAR(13) // 'FFT plans ready               ' // CHAR(13)

		if (.NOT.oldWisdom) call QTC_FFT_saveWisdom( dims, this % numT )
	end subroutine QTC_FFT_prepareFftPlans


	
	subroutine QTC_FFT_prepareLineshape( this )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

	integer eindex, numE
	double precision e, v, sum_lineshape

		numE = this % numE

		this % lineshape_e(:) = 0.d0
		do eindex=0, numE/4-1
			e = eindex * this % de
			v = this % F % getValue( e )
			this % lineshape_e(eindex) = v
			if (eindex.GT.0) this % lineshape_e(numE-eindex) = v
		end do
		sum_lineshape = sum( this % lineshape_e(:) )
!		print*, ' sum_lineshape=', sum_lineshape
!		this % lineshape_e(:) = this % lineshape_e(:) /sum_lineshape

		call fftw_execute_dft(this % plan1d, this % lineshape_e, this % lineshape_g)
		
	end subroutine QTC_FFT_prepareLineshape
	



	subroutine QTC_FFT_make_S_rg( this )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

	integer LX, LY, LZ, numE, kx, ky, kz, ss, k(1:3), eindex
	double precision e, emin
	type( Subsystem ), pointer :: s

		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

		allocate( this % S_keR( 0:LX-1, 0:LY-1, 0:LZ-1, 0:numE-1 ) ) ! one-use array

		do ss=0, this % numS - 1
			s => this % hamilt % subsystems( ss ) % obj
			emin = this % emin(ss)

			this % S_keR(:,:,:,:) = 0.d0
			do kz=0, LZ-1
				k(3) = kz
				do ky=0, LY-1
					k(2) = ky
					do kx=0, LX-1
						k(1) = kx
						e = s % dispersionLaw % getEnergy( k, this % sizes ) - emin
						eindex = int( 0.5d0 + e / this % de )
!print*, ' ss=', ss, ' eindex=', eindex, ' for e=', e, ' de=', this%de, ' emin=', emin
						if (eindex.LT.numE/2) then ! for zero-padding in convolution
							this % S_keR(kx,ky,kz,eindex) = 1.d0
						end if
					end do
				end do
			end do

			call fftw_execute_dft_r2c(this % plan4dR2C, this % S_keR(:,:,:,:), this % S_rgR(:,:,:,:,ss) )

		end do

		deallocate( this % S_keR ) ! not needed anymore
		
	end subroutine QTC_FFT_make_S_rg


	subroutine QTC_FFT_make_N_rg( this, n )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this
	type( Occupations ), pointer :: n

	integer LX, LY, LZ, numE, kx, ky, kz, k(1:3), eindex
	double precision e, emin
	type( Subsystem ), pointer :: s
	integer ss

		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

		do ss=0, this % numS - 1
			s => this % hamilt % subsystems( ss ) % obj
			emin = this % emin(ss)

			this % N_keR(:,:,:,:) = 0.d0

!$omp parallel do private(kx,ky,kz,k,e,eindex)
			do kz=0, LZ-1
				k(3) = kz
				do ky=0, LY-1
					k(2) = ky
					do kx=0, LX-1
						k(1) = MOD(kx, LX)
						e = s % dispersionLaw % getEnergy( k, this % sizes ) - emin
						eindex = int( 0.5d0 + e / this % de )
						if (eindex.LT.numE/2) then ! for zero-padding in convolution
							this % N_keR(kx,ky,kz,eindex) = n % getValue(k,s%occupationIndex)
						end if
					end do
				end do
			end do
!$omp end parallel do

			call fftw_execute_dft_r2c(this % plan4dR2C, this % N_keR, this % N_rgR(:,:,:,:,ss) )

		end do
		
	end subroutine QTC_FFT_make_N_rg


	logical function QTC_FFT_loadWisdom( dims, nt )

	use fftw3
	implicit none
	integer dims(0:3)
	integer nt

	integer void
	character*80 filename

		write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A)') './fftw_', &
&			dims(0),'_',dims(1),'_',dims(2),'_',dims(3),'_',nt,'.dat'
		void = fftw_import_wisdom_from_filename( TRIM(filename) // C_NULL_CHAR )
!		print*, ' load from filename=[' // TRIM(filename) // ']: ', void
		QTC_FFT_loadWisdom = (void.NE.0)
	
	end function QTC_FFT_loadWisdom

!---

	subroutine QTC_FFT_saveWisdom( dims, nt )

	use fftw3
	implicit none
	integer dims(0:3)
	integer nt

	integer void
	character*80 filename

		write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A)') './fftw_', &
&			dims(0),'_',dims(1),'_',dims(2),'_',dims(3),'_',nt,'.dat'
		void = fftw_export_wisdom_to_filename( TRIM(filename) // C_NULL_CHAR )
!		print*, '   save to filename=[' // TRIM(filename) // ']: ', void
		if (void.EQ.0) stop 'save wisdom fail'
		
	end subroutine QTC_FFT_saveWisdom


!========================= calculation ==============================================

	subroutine QTC_FFT_addValue( this, t, n, dn, rta )

	use omp_lib
	
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this
	type( Occupations ), pointer :: n, dn
	double precision t
	logical rta
	
	integer ss, k1multiplier

!		print*, ' addValue: rta=', rta

		call this % latePreparation ! create arrays if needed
		call this % prepareCacheForKroneckers( rta ) ! prepare cache if needed
		
		call this % make_N_rg( n )
		do ss=0, this % numS - 1

!			print*, ' ss=', ss
			call this % calculateTermForRightPart( n, dn, ss, .TRUE.,  1, rta )
			if (.NOT.rta) then
				call this % calculateTermForRightPart( n, dn, ss, .FALSE., 1, rta )
				if (this % useKroneckers) then
					do k1multiplier=2, this % maxk1multiplier
						call this % calculateTermForRightPart( n, dn, ss, .TRUE.,  k1multiplier, rta )
						call this % calculateTermForRightPart( n, dn, ss, .FALSE., k1multiplier, rta )
					end do
				end if
			end if
		end do
		call this % simpleCalculationForKroneckers( n, dn, rta )
		
	end subroutine QTC_FFT_addValue

! -- straight approach : plain calculation of original sum -----------------------------------
!  case 2.1 : 1 or more Kroneckers -> less terms

	subroutine QTC_FFT_simpleCalculationForKroneckers( this, n, dn, rta )
	class( QbeTermsCalculatorFFT ) :: this
	type( Occupations ), pointer :: n, dn
	logical rta
	
	integer ss, pp
	
	type(Subsystem), pointer :: s
	type( Perturbation ), pointer :: p

	class( QbeTermFFT ), pointer :: qt
	class( QbeTermSumOfABCD ), pointer :: qt_abcd
	class( MDC_CachedList ), pointer :: cache
	integer i, half, i2, rtaIndex
	type(MultidimensionalCycle), pointer :: mdc

		rtaIndex = 0
		if (rta) rtaIndex = 1

		do ss=0, this % numS - 1
			s => this % hamilt % subsystems( ss ) % obj
			do pp=0, this % numP - 1
				if (this % term_ps(pp,ss)) then
					p => this % hamilt % perturbations(pp) % obj
					if (rta .OR. (p % mode .EQ. Q_DEPENDENT)) then
						do half=1, this%numHalves(pp)
							if (half.EQ.2) p => p % hc
							i2 = p % indexOf(s)
							qt => this % term(half,pp) % obj

							if ( qt % partMultiplier( i2 ).NE.0 ) then
								qt_abcd => qt % QbeTermSumOfABCD
	
								cache => this % mdc_cache(half,pp,rtaIndex) % obj
!print*, ' cache(', pp,half,rtaIndex, ') : associated=', associated(cache)
								call qt % setCurrentSubsystemIndex( i2 )
!$omp parallel do private(i,mdc)
								do i=0, this % numT-1
								mdc => this % mdc(i,half, pp) % obj
									call mdc % runWithKroneckers( n, qt_abcd, &
&										i2, rta, dn, cache, i )
!										s % occupationIndex, rta, dn, cache, i )
								end do
!$omp end parallel do
							end if ! coeff.NE.0
						end do ! half
					end if ! rta or V(q)
				end if ! term_ps		
			end do ! pp
		end do ! ss

	end subroutine QTC_FFT_simpleCalculationForKroneckers





	subroutine QTC_FFT_prepareCacheForKroneckers( this, rta )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this
	logical rta
	
	integer pp, iter

	integer ss
	type(Subsystem), pointer :: s
	type( Perturbation ), pointer :: p

	class( QbeTermFFT ), pointer :: qt
	class( QbeTermSumOfABCD ), pointer :: qt_abcd
	class( MDC_CachedList ), pointer :: cache
	integer i, half, i2, rtaIndex
	type(MultidimensionalCycle), pointer :: mdc

	integer*8 mem
	integer MB

		rtaIndex = 0
		if (rta) rtaIndex = 1
!print*, ' PREPARED CACHE? ', this % cachePrepared(rtaIndex)
		if (.NOT. this % cachePrepared(rtaIndex)) then
		
			mem = 0
			do pp=0, this % numP - 1
				p => this % hamilt % perturbations(pp) % obj
				if (rta .OR. (p % mode .EQ. Q_DEPENDENT)) then
					do half=1, this%numHalves(pp)
!						print*, ' prepareCacheForKroneckers: pp=', pp, ' rta=', rta, ' half=', half
						if (half.EQ.2) p => p % hc

						do iter=1,2
							if (iter.EQ.1) then
								cache => MDC_CachedList( &
&	p%numAllOperators-1, p%numDeltaCombinations, &
&	this%sizes%volume, p%numSubsystems )
								this % mdc_cache(half,pp,rtaIndex) % obj => cache 
							end if
							call cache % setState( iter )

							do ss=0, this % numS - 1
!					do ss=this % numS - 1, 0, -1
								s => this % hamilt % subsystems( ss ) % obj
								if (this % term_ps(pp,ss)) then
									i2 = p % indexOf(s)
									qt => this % term(half,pp) % obj

									if ( qt % partMultiplier( i2 ).NE.0 ) then
										qt_abcd => qt % QbeTermSumOfABCD
										call qt % setCurrentSubsystemIndex( i2 )
!$omp parallel do private(i,mdc)
						do i=0, this % numT-1
							mdc => this % mdc(i,half, pp) % obj
							call mdc % prepareCacheForKroneckers( qt_abcd, &
&								i2, rta, cache, iter, i )
!								s % occupationIndex, rta, cache, iter, i )
						end do
!$omp end parallel do
										if (iter.EQ.2) mem = mem + cache % getMemoryFootprint()
									end if ! coeff.NE.0
								end if ! term_ps		
							end do ! ss
						end do ! iter
					end do ! half
								
!					print*, 'DONE prepareCacheForKroneckers: pp=', pp, ' rta=', rta
				end if ! rta .OR. V(q)
			end do ! pp
			this % cachePrepared(rtaIndex) = .TRUE.
			if (mem.GT.0) then
				MB = int( (8.d0*mem)/1048576 + 0.5d0)
				write(*, '(A,I0,A)', advance='no') 'memory size for cache: ', MB, 'MB      '//CHAR(13)
			end if
		end if
	end subroutine QTC_FFT_prepareCacheForKroneckers


! -- FFT approach -------------------------------------------------------------------------------------
!  case 1.1 :  U=const + all terms including Kroneckers
!  case 1.2 :  U=U(q)  + first term without Kroneckers


	subroutine QTC_FFT_calculateTermForRightPart( this, n, dn, ss, with_Nk, k1multiplier, rta )
	class( QbeTermsCalculatorFFT ) :: this
	type( Occupations ), pointer :: n, dn
	integer ss, k1multiplier
	logical with_Nk, rta, enableNk

		if (this % make_calc_rg( ss, with_Nk, k1multiplier, n, rta )) then
!			print*, ' nonzero'   ! i.e. if number of terms .GT. 0

			call this % make_calc_ke

			enableNk = (.NOT.rta)
			call this % updateDN( ss, with_Nk.AND.enableNk, &
&				k1multiplier, n, dn )
		end if

	end subroutine QTC_FFT_calculateTermForRightPart



	subroutine QTC_FFT_constant_interaction_XYZ( this, pp,ss,externalIndex, &
&		posForQ, p, qt, ABCD, g, kroIndex, partMultiplier, half )
	class( QbeTermsCalculatorFFT ) :: this
	type( Perturbation ), pointer :: p
	class( QbeTermFFT ), pointer :: qt
	type(ABCD_info), pointer :: ABCD

	integer pp,ss,externalIndex,posForQ(0:1),g, kroIndex, partMultiplier, half
	
	integer LX, LY, LZ, numE, multiplier
	integer x,y,z,mx,my,mz,i,NorS, mg, sss
	complex*16 result, term, factor
	logical secondHalf
	logical dbg
	
		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE
		
		do z=0, LZ-1
		do y=0, LY-1
		do x=0, LX/2
			dbg = ((x.EQ.0).AND.(y.EQ.0).AND.(z.EQ.0).AND.(g.EQ.0))
			dbg = .FALSE.
!			if (dbg) print*, 'term factors:  coeff=', ABCD%coeff , ' partMultiplier=', partMultiplier, &
!&				' posM=', qt % positionMultiplier(0:3, kroIndex)
		
			result = 0.d0

			if (p%mode .EQ. Q_DEPENDENT) then
				stop 'not expected to have U(q)'
				term = this % calc_r_many(x,y,z,g)
			else                           
				term = 1.d0
			end if

			do i=0, qt % order - 1 ! operators
				if ( i .NE. externalIndex ) then
					NorS = ABCD % NS(i) ! 1=N or 0=S 
	
					sss = p % info(i) % subsystem % occupationIndex
					if (sss.GE.this%numS) stop ' sss GT numS'

					multiplier = qt % positionMultiplier(i, kroIndex)
					if (multiplier.NE.0) then
						mx = MOD(x*multiplier, LX)
						my = MOD(y*multiplier, LY)
						mz = MOD(z*multiplier, LZ)
						mg = MOD(g*multiplier, numE)
						secondHalf = ( mx.GT.(LX/2) )
						if (secondHalf) then
							mx = MOD(LX-mx, LX)
							my = MOD(LY-my, LY)
							mz = MOD(LZ-mz, LZ)
							mg = MOD(numE-mg, numE)
						end if
						factor = -1.d3 ! debug
						if (NorS.EQ.0) then
							if (dbg) print*, 'i=', i, ' factor=S'
							factor = this % S_rgR(  mx, my, mz, mg, sss )
						else
							if (dbg) print*, 'i=', i, ' factor=N'
							factor = this % N_rgR(  mx, my, mz, mg, sss )
						end if
						if (this % minusR(i,half,ss,pp) .NEQV. secondHalf) then !  -R => conjg(R) 
							factor = conjg( factor )
						end if
						term = term * factor
					end if

				end if ! not k0
			end do ! operators

			result = term * ABCD%coeff * partMultiplier * this % lineshape_g(g)
			result = result * ((p % amplitudeValue)**2 )
!			if (part.EQ.1) result=-result

			this % calc_rgR(x,y,z,g) = this % calc_rgR(x,y,z,g) + result
			
		end do
		end do
		end do

	end subroutine QTC_FFT_constant_interaction_XYZ




	function QTC_FFT_make_calc_rg( this, ss, with_Nk, k1multiplier, n, rta ) result(nonzero)
	implicit none
	class( QbeTermsCalculatorFFT ) :: this
	integer ss
	logical with_Nk
	integer k1multiplier
	type( Occupations ), pointer :: n
	logical rta

	integer LX, LY, LZ, numE
	integer x,y,z,g, mx, my, mz, mg, pp, i, t, sss, kroIndex, maxKroIndex
	integer coeff, s_index, NorS, externalIndex
	logical minusRforBC(0:1)
	type( Perturbation ), pointer :: p
	class( QbeTermFFT ), pointer :: qt
	type(ABCD_info), pointer :: ABCD
	type(Subsystem), pointer :: s

	complex*16 term, result, factor
	integer versionForQ, posForQ(0:1), partMultiplier, half
	logical nonzero, rtaUseDagger, sameDagger, good
	
		nonzero = .FALSE.
		s => this % hamilt % subsystems( ss ) % obj

		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE
!		print*, ' LX LY LZ numE numS numP=', LX, LY, LZ, numE, this % numS, this % numP

		this % calc_rgR(:,:,:,:) = 0.d0

		do pp=0, this % numP - 1
			if (this % term_ps(pp,ss)) then
			  p => this % hamilt % perturbations(pp) % obj
			  do half=1, this%numHalves(pp)
				if (half.EQ.2) p => p % hc
				qt => this % term(half,pp) % obj
				s_index = p % indexOf( s )
				externalIndex = p % firstOperatorEntry( s_index )
				
				if (p%mode .EQ. Q_DEPENDENT) then
					versionForQ = 0
					if ( (s_index.EQ.p%operatorPositionForQ(0,0)) .OR. &
&							(s_index.EQ.p%operatorPositionForQ(1,0))) then
						versionForQ = 1
					end if
					posForQ(:) = p % operatorPositionForQ(:,versionForQ)
!					print*, ' posForQ=', posForQ
				else
					posForQ(:) = -1
				end if			
				maxKroIndex = 0
				if (this % useKroneckers) then
					if ((p%mode .EQ. Q_INDEPENDENT).AND.(.NOT.rta))  then ! only V(0), dn/dt
						maxKroIndex = qt % kroVolume-1
					end if
				end if
				do kroIndex = 0, maxKroIndex
				  if (k1multiplier .EQ. qt % positionMultiplier(externalIndex, kroIndex) ) then
						do t=0, qt % numTerms(kroIndex, s_index)-1
							ABCD => qt % info(t, kroIndex, s_index)
							coeff = ABCD % coeff
							if (coeff .NE. 0) then
								partMultiplier = qt % partMultiplier( s_index) 
								NorS = ABCD % NS(externalIndex)
								if ( (NorS .EQ. 1) .EQV. with_Nk) then

									nonzero = .TRUE.
									if (p%mode.EQ.Q_INDEPENDENT) then
!$omp parallel do private(x,y,z,sss,result,term,i,NorS,factor,mx,my,mz,g,mg)
										do g=0, numE-1
											call this % constant_interaction_XYZ( &
&	pp, ss, externalIndex, posForQ, &
&	p, qt, ABCD, g, kroIndex, partMultiplier, half )
										end do
!$omp end parallel do
									else
					
	minusRforBC(:) = .FALSE.
	if (posForQ(0).GE.0) minusRforBC(0) = this % minusR( posForQ(0),half,ss,pp )
	if (posForQ(1).GE.0) minusRforBC(1) = this % minusR( posForQ(1),half,ss,pp )
!$omp parallel do private(g)
										do g=0, numE-1	
		call this % make_BC_for_Q( p, qt, ABCD, g, posForQ, minusRforBC, kroIndex ) ! calc_r
										end do ! g
!$omp end parallel do
	call this % make_calc_k_3D_many

!$omp parallel do private(g)
										do g=0, numE-1	
			this % calc_k_many(:,:,:,g) = this % calc_k_many(:,:,:,g) &
&				 * ((p % amplitudeTable % values(:,:,:))**2 )
										end do ! g
!$omp end parallel do

	call this % make_calc_r_3D_many

!$omp parallel do private(z,g)
										do g=0, numE-1	
											do z=0, LZ-1
												call this % general_interaction_XY( &
&		pp, ss, externalIndex, posForQ, &
&		p, qt, ABCD, g, z, kroIndex, partMultiplier, half )
											end do
										end do ! g
!$omp end parallel do
									end if ! Q_DEPENDENT/Q_INDEPENDENT
								end if ! with_Nk
							end if ! coeff <> 0
						end do ! terms in ABCD
				  end if ! expected multiplier for k1
				end do ! kroIndex
			  end do ! half
			end if ! pp contains ss
		end do ! perturbations

!		print*, ' calc_rg: nonzero=', nonzero
!		print*, ' calc_rg=', sum(this%calc_rg(:,:,:,:)**2)

	end function QTC_FFT_make_calc_rg


	subroutine QTC_FFT_general_interaction_XY( this, &
&		pp,ss,externalIndex, posForQ, &
&		p, qt, ABCD, g, z, kroIndex, partMultiplier, half )
	class( QbeTermsCalculatorFFT ) :: this
	type( Perturbation ), pointer :: p
	class( QbeTermFFT ), pointer :: qt
	type(ABCD_info), pointer :: ABCD

	integer pp,ss,externalIndex,posForQ(0:1),g,z, kroIndex, partMultiplier, half
	
	integer LX, LY, LZ, numE, kMultiplier
	integer x,y,mx,my,i,NorS, mz,mg, sss
	complex*16 result, term, factor
	logical secondHalf, dbg
	
!		print*, 'XY not ready'
		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

		do y=0, LY-1
		do x=0, LX/2
			dbg = ((x.EQ.0).AND.(y.EQ.0).AND.(z.EQ.0).AND.(g.EQ.0))
			dbg = .FALSE.
			if (dbg) print*, 'term factors:  coeff=', ABCD%coeff , ' partMultiplier=', partMultiplier, &
&				' posM=', qt % positionMultiplier(0:3, kroIndex)
		
			result = 0.d0

			if (p%mode .EQ. Q_DEPENDENT) then
				term = this % calc_r_many(x,y,z,g)
			else                           
				term = 1.d0
				stop 'unexpected: XY V_0'
			end if

			do i=0, qt % order - 1 ! operators
				if ( (i .NE. externalIndex).AND.(i.NE.posForQ(0)).AND.(i.NE.posForQ(1)) ) then
					dbg = .FALSE.
					kMultiplier = qt % positionMultiplier(i, kroIndex)
					if (kMultiplier.NE.0) then
						mx = MOD(x*kMultiplier, LX)
						my = MOD(y*kMultiplier, LY)
						mz = MOD(z*kMultiplier, LZ)
						mg = MOD(g*kMultiplier, numE)
						secondHalf = ( mx.GT.(LX/2) )
						if (secondHalf) then
							mx = MOD(LX-mx, LX)
							my = MOD(LY-my, LY)
							mz = MOD(LZ-mz, LZ)
							mg = MOD(numE-mg, numE)
						end if

						NorS = ABCD % NS(i) ! 1=N or 0=S
	
						sss = p % info(i) % subsystem % occupationIndex
						if (sss.GE.this%numS) stop ' sss GT numS'

						factor = 1
						if (NorS.EQ.0) then
if (dbg) print*, 'i=', i, ' factor=S'
							factor = this % S_rgR(  mx, my, mz, mg, sss )
						else
if (dbg) print*, 'i=', i, ' factor=N'
							factor = this % N_rgR(  mx, my, mz, mg, sss )
						end if
						if ( this % minusR(i,half,ss,pp) .NEQV. secondHalf ) then !  -R => conjg(R) 
							factor = conjg( factor )
						end if
						
						term = term * factor
					end if

				end if ! not k0
			end do ! operators

			result = term * ABCD%coeff * partMultiplier * this % lineshape_g(g)
			this % calc_rgR(x,y,z,g) = this % calc_rgR(x,y,z,g) + result
			
		end do
		end do

	end subroutine QTC_FFT_general_interaction_XY
	


	subroutine QTC_FFT_make_BC_for_Q(this, p, qt, ABCD, g, pos, minusR, kroIndex )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this
	type( Perturbation ), pointer :: p
	class( QbeTermFFT ), pointer :: qt
	type(ABCD_info), pointer :: ABCD
	integer g, pos(0:1)
	logical minusR(0:1)
	integer kroIndex

	integer x,y,z, LX,LY,LZ, s_occupationIndex, NorS, i, mx,my,mz, mg, multiplier, numE
	complex*16 v, factor
	logical secondHalf, dbg

		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE
!		print*, ' ...forQ: multipliers = ', qt % positionMultiplier( pos(:) )
!		print*, ' ...forQ: minusR = ', minusR(:)

!$omp parallel do private(x,y,z,v,i,s_occupationIndex,multiplier,secondHalf,mx,my,mz,mg,NorS,factor,dbg)
		do z=0, LZ-1
		do y=0, LY-1
		do x=0, LX-1
			dbg = ((x.EQ.0).AND.(y.EQ.0).AND.(z.EQ.0).AND.(g.EQ.0))
			dbg = .FALSE.
			if (dbg) print*, 'term factors:  coeff=', ABCD%coeff , &
&				' posM=', qt % positionMultiplier(0:3, kroIndex)

			v = 1.d0
			do i=0,1
				if (pos(i).GE.0) then
					s_occupationIndex = p % info( pos(i) ) % subsystem % occupationIndex

					multiplier = qt % positionMultiplier( pos(i), kroIndex )
					if (multiplier.NE.0) then
				
						mx = MOD(x*multiplier, LX)
						my = MOD(y*multiplier, LY)
						mz = MOD(z*multiplier, LZ)
						mg = MOD(g*multiplier, numE)
						secondHalf = ( mx.GT.(LX/2) )
						if (secondHalf) then
							mx = MOD(LX-mx, LX)
							my = MOD(LY-my, LY)
							mz = MOD(LZ-mz, LZ)
							mg = MOD(numE-mg, numE)
						end if

						factor = 1
						NorS = ABCD % NS( pos(i) )
						if (NorS.EQ.0) then
							if (dbg) print*, 'i=', pos(i), ' factor=S'
							factor = this % S_rgR(  mx, my, mz, mg, s_occupationIndex )
						else
							if (dbg) print*, 'i=', pos(i), ' factor=N'
							factor = this % N_rgR(  mx, my, mz, mg, s_occupationIndex )
						end if
						if (minusR(i) .NEQV. secondHalf) then ! 1=> -R => conjg(R)
							factor = conjg( factor )
						end if
						v = v*factor
					end if
				end if
			end do

			this%calc_r_many(x,y,z,g) = v

		end do
		end do
		end do
!$omp end parallel do

	end subroutine QTC_FFT_make_BC_for_Q










	subroutine QTC_FFT_make_calc_ke( this )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

		call fftw_execute_dft_c2r(this % plan4dBackC2R, this % calc_rgR, this % calc_keR )

		this % calc_keR(:,:,:,:) = this % calc_keR(:,:,:,:) / this % volume4D ! 1/V for backwards transform

	end subroutine QTC_FFT_make_calc_ke








	subroutine QTC_FFT_updateDN( this, ss, with_Nk, k1multiplier, n, dn )
	class( QbeTermsCalculatorFFT ) :: this
	integer ss
	logical with_Nk
	integer k1multiplier
	type( Occupations ), pointer :: n, dn

	integer LX, LY, LZ, numE
	integer kx, ky, kz, eindex, k(1:3), ind, me,mx,my,mz
	double precision e, rp, emin
	type(Subsystem), pointer :: s

		s => this % hamilt % subsystems( ss ) % obj
		emin = this % emin( ss )
		ind = s % occupationIndex
!		print*, ' updateDN: ss=', ss, ' ind=', ind
		LX = this % sizes % L(1)
		LY = this % sizes % L(2)
		LZ = this % sizes % L(3)
		numE = this % numE

!		print*,' updateDN: multiplier=', k1multiplier, ' with_Nk=', with_Nk
!$omp parallel do private(kx,ky,kz,k,e,eindex,rp,mx,my,mz,me)
		do kz=0, LZ-1
			k(3) = kz
			mz = MOD(kz*k1multiplier, LZ)
		do ky=0, LY-1
			k(2) = ky
			my = MOD(ky*k1multiplier, LY)
		do kx=0, LX-1
			k(1) = kx
			mx = MOD(kx*k1multiplier, LX)

			e = s % dispersionLaw % getEnergy( k, this % sizes ) - emin
			eindex = int( 0.5d0 + e / this % de )
			if (eindex.LT.numE/2) then ! for zero-padding in convolution
				me = MOD(eindex*k1multiplier, numE)
				rp = this % calc_keR(mx,my,mz,me) 
				if (with_Nk) then
					rp = rp * n % values(kx,ky,kz,ind)
				end if
!				dn % values(kx, ky, kz, ind ) = dn % values(kx, ky, kz, ind )  +  rp
				call dn % incValue( k, ind,    rp )
!				if ((kx.EQ.0).AND.(ky.EQ.0).AND.(kz.EQ.0)) print*, ' dN_000=', rp
			end if
			
		end do
		end do
		end do
!$end parallel do
		
	end subroutine QTC_FFT_updateDN




! --------------

	subroutine QTC_FFT_make_calc_k_3D_many( this )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

		call fftw_execute_dft(this % plan3dManyBack, this % calc_r_many, this % calc_k_many )
			
		this % calc_k_many(:,:,:,:) = this % calc_k_many(:,:,:,:) / this % volume3D ! 1/V for backwards transform
		
	end subroutine QTC_FFT_make_calc_k_3D_many

	subroutine QTC_FFT_make_calc_r_3D_many( this )
	class( QbeTermsCalculatorFFT ), intent(INOUT) :: this

		call fftw_execute_dft(this % plan3dMany, this % calc_k_many, this % calc_r_many )
			
	end subroutine QTC_FFT_make_calc_r_3D_many




! -- for debug

	subroutine QTC_FFT_show( this, o )
	class( QbeTermsCalculatorFFT ) :: this
	class( Output ), pointer :: o
	
	integer i, j, i2, half
	type( Subsystem ), pointer :: s
	type( Perturbation ), pointer :: p
	class( QbeTermFFT ), pointer :: qt
	character*1 a, k
	character*10 kbold, a_brackets

		call o % write( 'QbeTermsCalculator: fast with FFT' )

		do i=0, this % hamilt % numSubsystems-1
			s => this % hamilt % subsystems(i) % obj
			a = s % operatorLetter
			a_brackets = '(' // a // ')'
			k = s % indexLetter
			call o % bold( k, kbold )

			do j=0, this % hamilt % numPerturbations-1
				p => this % hamilt % perturbations(j) % obj
				i2 = p % indexof(s)
				if (i2.GE.0) then
					call o % beginEquation()
					call o % switchSuperscriptText( .TRUE. )
					call p % showQbeTermTitle( i2, o )

					do half=1, this%numHalves(j)
						qt => this % term(half,j) % obj

						call qt % setCurrentSubsystemIndex( i2 )
						call qt % show( o )
					end do
					call o % endEquation()
				end if
			end do
		end do
	
	end subroutine QTC_FFT_show


