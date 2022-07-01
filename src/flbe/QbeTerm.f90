module QbeTerm_module

	use Perturbation_module
	
	implicit none

	type, abstract :: QbeTerm
		character*50 title
		type( Perturbation ), pointer :: p
		logical useKroneckers
		integer kroVolume
	contains
		procedure( QbeTermShow ), deferred :: show
	end type QbeTerm

	type, extends( QbeTerm ) :: QbeTermAsIs
		logical, allocatable :: subsystemHasDaggers(:)
		integer, allocatable :: subsystemIndex(:)

		integer order

		logical, allocatable :: bose(:) ! operator(i) bose from Perturbation % info(:)
		logical, allocatable :: dagger(:) ! operator(i) dagger from Perturbation % info(:)
!		integer, allocatable :: numParts(:) ! 1 or 2 (with H.c.)
		integer, allocatable :: partMultiplier(:) ! part x s_index ; 1 1 for a+ a, 2 2 for a+a+ aa, 2 0 for a+a+ etc.
		logical, allocatable :: partDagger(:) ! a+a -> T,F ; aa -> F,x ; a+a+ -> T,x
							
		integer currentSubsystemIndex ! current value for show() and getValue()
	contains
		procedure :: setPerturbation => QTAI_setPerturbation
		procedure :: getValue => QTAI_getValue
		procedure :: addValueForRTA => QTAI_addValueForRTA
		procedure :: setCurrentSubsystemIndex => QTAI_setCurrentSubsystemIndex

		procedure :: show => QTAI_show
		procedure :: showHead => QTAI_showHead
		procedure :: showTail => QTAI_showTail
		procedure :: showSingleTerm => QTAI_showSingleTerm
		procedure :: showVFDeltaKDeltaE => QTAI_showVFDeltaKDeltaE

		final :: QTAI_destructor
		
	end type QbeTermAsIs

	type, extends( QbeTermAsIs ) :: QbeTermGeneralSum
		! 0:2, 0:L-1, 0:1
		! where  0:2=coeffs for 1, n  and delta
		! 	0:L-1 is operator
		! 	0:1 is half of formula: nnn-nnn
		integer, allocatable :: A_b(:,:,:)
		! 0:L-1, 0:1
		integer, allocatable :: A_b_desc(:,:)	! description index for factor
		character*30, pointer :: desc(:) ! 0='n', 1='n+1', 2='1-n'
	contains
		procedure :: getValue => QTGS_getValue

		procedure :: setTwoCoeffs => QTGS_setTwoCoeffs ! (1 + n) / (n-1) / n / 1 etc
		procedure :: showSingleTerm => QTGS_showSingleTerm

		final :: QTGS_destructor
		
	end type QbeTermGeneralSum

	type QbeTermAsIs_p
		class( QbeTermAsIs ), pointer :: obj
	end type QbeTermAsIs_p

	type ABCD_info
		integer coeff
		integer*1, dimension(0:max_num_operators-1) :: NS ! S=0, N=1, kronecker=2+
		integer addr ! its own index in QbeTerm % info(), for fast access
	end type ABCD_info
	

	type, extends( QbeTermAsIs ) :: QbeTermSumOfABCD
		integer maxNumTerms

		integer, allocatable :: numTerms(:,:)	! (kroneckerIndex,subsystemIndex)

		type(ABCD_info), allocatable :: info(:,:,:)	! (operatorIndex,kroneckerIndex,subsystemIndex)

!		integer kroIndex
		logical, allocatable :: skipIndex(:,:) ! i, kroIndex
		integer, allocatable :: positionMultiplier(:,:)
		integer, allocatable :: replacementPosition(:,:) ! delta_12 -> [2]=1
	contains
		procedure :: initFromPerturbation => QTSumABCD_initFromPerturbation
		procedure :: initFromGeneralSum => QTSumABCD_initFromGeneralSum
		procedure :: simplify => QTSumABCD_simplify

		procedure :: fillKroneckerInfo => QTSumABCD_fillKroneckerInfo

		procedure :: show => QTSumABCD_show
		procedure :: showHeadForKronecker => QTSumABCD_showHeadForKronecker
		procedure :: showSingleTermForKronecker => QTSumABCD_showSingleTermForKronecker

		procedure :: showVdKdEforKroneckers => QTSumABCD_showVdKdEforKroneckers

		procedure :: writeKroneckerFactor => QTSumABCD_writeKroneckerFactor
		procedure :: calculateKroneckerFactor => QTSumABCD_calculateKroneckerFactor
		
		procedure :: getValue => QTSumABCD_getValue
		procedure :: getValueForKronecker => QTSumABCD_getValueForKronecker
		procedure :: addValueForRTAForKronecker => QTSumABCD_addValueForRTAForKronecker

		final :: QtSumABCD_destructor		
		
	end type QbeTermSumOfABCD

	type QbeTermABCD_p
		class( QbeTermSumOfABCD ), pointer :: obj
	end type QbeTermABCD_p

	type, extends( QbeTermSumOfABCD ) :: QbeTermFFT
	contains
		procedure :: show => QbeTermFFT_show
		procedure :: showTermsForKronecker => QbeTermFFT_showTermsForKronecker
		procedure :: countTerms => QbeTermFFT_countTerms
		procedure :: showHeadWithMultiplier => QbeTermFFT_showHeadWithMultiplier
		procedure :: showTail => QbeTermFFT_showTail
		procedure :: showVFDeltaKDeltaE => QbeTermFFT_showVFDeltaKDeltaE

	end type QbeTermFFT

	type QbeTermFFT_p
		class( QbeTermFFT ), pointer :: obj
	end type QbeTermFFT_p

! constructors

	interface QbeTermAsIs
		module procedure QTAI_init
	end interface QbeTermAsIs
	interface QbeTermGeneralSum
		module procedure QTGS_init
	end interface QbeTermGeneralSum
	interface QbeTermSumOfABCD
		module procedure QTSumABCD_fromPerturbation
	end interface QbeTermSumOfABCD
	interface QbeTermFFT
		module procedure QbeTermFFT_fromPerturbation ! from perturbation -> use ABCD
	end interface QbeTermFFT

! polymorphic subroutine
	abstract interface
		subroutine QbeTermShow( this, o )
			import QbeTerm
			import Output
			class(QbeTerm) :: this
			class(Output), pointer :: o
		end subroutine
	end interface

contains

!---------------------------------------------------------

	include 'QbeTermAsIs.f90'
	include 'QbeTermGeneralSum.f90'
	include 'QbeTermABCD.f90'
	include 'QbeTermFFT.f90'

end module QbeTerm_module


