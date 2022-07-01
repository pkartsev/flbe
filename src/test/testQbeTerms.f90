subroutine testQbeTerms

use flbe
implicit none

type(Subsystem), pointer :: sf,sb,sf2,sb2
class(DispersionRelation), pointer :: dr
class(Hamiltonian), pointer :: hamilt
type(Perturbation), pointer :: v_ee, v_eph, v_bb, v_21, v_12, p, v_20, v_31, v_13
class(QbeTermAsIs), pointer :: qtai, qtgs, qtabcd

class(Output), pointer :: oH, oT

integer, parameter :: numTry = 1
integer version, t, i, numOp, numSubsystems, ss
real r
double precision, allocatable :: n(:)
integer, allocatable :: k(:)
integer kro, excludePosition
logical useKroneckers, good
double precision theor, theor2, v0,v1,v2,v3, n0,n1,n2,n3

	call printTwice_openFile('testQbeTerms.log', ' [ test QbeTerms...')

!-----------------------------------

	excludePosition = -1

	dr => DispersionRelationConst( 0.d0 )

	sf => SubsystemFermi( dr )
	sf2 => SubsystemFermi( dr )
	sb => SubsystemBose( dr )
	sb2 => SubsystemBose( dr )

	v_ee => InteractionPairwise( sf, sf ) ! spin-up x spin-down
	v_bb => InteractionPairwise( sb ) ! same Bose system
	v_20 => Perturbation( '2-0') ! general case
	call v_20 % addSubsystem( sb, 2, 0)
	v_21 => Perturbation( '2-1') ! general case
	call v_21 % addSubsystem( sb, 2, 1)
	v_12 => Perturbation( '1-2') ! general case
	call v_12 % addSubsystem( sb, 1, 2)
	v_31 => Perturbation( '3-1') ! general case
	call v_31 % addSubsystem( sb, 3, 1)
	v_13 => Perturbation( '1-3') ! general case
	call v_13 % addSubsystem( sb, 1, 3)

	numSubsystems = 1
	ss = 0

  do kro=0, 1
  	useKroneckers = kro.GT.0
	call printTwiceSL( ' useKroneckers: ', useKroneckers)
!	write(*, *) ' useKroneckers: ', useKroneckers
	
	do version=0, 6
		select case(version)
			case(0)
				p => v_ee
				call printTwice('  * V_ee')
			case(1)
				p => v_bb
				call printTwice('  * V_bb')
			case(2)
				p => v_20
				call printTwice('  * V_20')
			case(3)
				p => v_21
				call printTwice('  * V_21')
			case(4)
				p => v_12
				call printTwice('  * V_12')
			case(5)
				p => v_31
				call printTwice('  * V_31')
			case(6)
				p => v_13
		end select
		call p % setAmplitudeConst( 1.d0 )
		hamilt => Hamiltonian()
		call hamilt % addPerturbation( p )
		
		call p % prepare
		qtai   => QbeTermAsIs( p, useKroneckers )
		qtgs   => QbeTermGeneralSum( p, useKroneckers )
		qtabcd => QbeTermSumOfABCD( p, useKroneckers )

		numOp = p % numAllOperators
		allocate( k(0:numOp-1) )
		allocate( n(0:numOp-1) )
		do t=0, numTry-1
			write(debugFileNumber, *) '  try # ', t
			if (verbose) write(*, '(A)', advance='no') '  n(:)='
			if (useKroneckers) then
				k(:) = 666
				call RANDOM_NUMBER(r)
				n(:) = (1+int(r*7))/8.d0 ! 0.125, 0.25, ... : exact values, easy to check
				call printTwiceSFS( '', n(0), '')
			else
				do i=0, numOp-1
					k(i) = i ! all different

					call RANDOM_NUMBER(r)
					n(i) = (1+int(r*7))/8.d0 ! 0.125, 0.25, ... : exact values, easy to check
					call printTwiceSFS( ' ', n(i), ' ')
				end do
			end if
			call printTwice( ' ' )
			
			theor = 18
			theor2 = 19

			n0 = n(0)
			n1 = n(1)
			select case(version)
				case(0)
					n2 = n(2)
					n3 = n(3)
					theor  = 1*( (1-n0) * (1-n1) * n2 * n3 - n0 * n1 * (1-n2) * (1-n3) )
				case(1)
					n2 = n(2)
					n3 = n(3)
					if (useKroneckers) then
						theor  = 2*( (1+n0) * (2+n1) * n2 * (n3-1) - n0 * (n1-1) * (1+n2) * (2+n3) )
					else
						theor  = 2*( (1+n0) * (1+n1) * n2 * n3 - n0 * n1 * (1+n2) * (1+n3) )
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F8.3)') &
&							'  theor=2*[', 1+n0, ' *', 1+n1, ' *', n2, ' *', n3
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F8.3,A,F10.6)') &
&							' -', n0, ' *', n1, ' *', 1+n2, ' *', 1+n3, ' ]=', theor
					end if
				case(2)
					if (useKroneckers) then
						theor  = 2*( (1+n0) * (2+n1) - n0 * (n1-1) )
						theor2 = 666.d0
					else
						theor  = 2*( (1+n0) * (1+n1) - n0 * n1 )
						theor2 = 666.d0
					end if
				case(3)
					n2 = n(2)
					if (useKroneckers) then
						theor  = 2*( (1+n0) * (2+n1) * n2   -   n0 * (n1-1) * (1+n2) )
						theor2 = 1*( (1+n0) * (n1-1) * n2   -   n0 * (2+n1) * (1+n2) )
					else
						theor  = 2*( (1+n0) * (1+n1) * n2   -   n0 * n1 * (1+n2) )
						theor2 = 1*( (1+n0) * n1 * n2   -   n0 * (1+n1) * (1+n2) )
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3)') &
&							'  theor=2*[', 1+n0, ' *', 1+n1, ' *', n2
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F10.6)') &
&							' -', n0, ' *', n1, ' *', 1+n2, ' ]=', theor
					end if
				case(4)
					n2 = n(2)
					if (useKroneckers) then
						theor  = 1*( (1+n0) * n1 * (n2-1)   -   n0 * (1+n1) * (2+n2) )
						theor2 = 2*( (2+n0) * (1+n1) * n2   -   (n0-1) * n1 * (1+n2) )
					else
						theor  = 1*( (1+n0) * n1 * n2   -   n0 * (1+n1) * (1+n2) )
						theor2 = 2*( (1+n0) * (1+n1) * n2   -   n0 * n1 * (1+n2) )
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3)') &
&							'  theor=1*[', 1+n0, ' *', n1, ' *', n2
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F10.6)') &
&							' -', n0, '*', 1+n1, '*', 1+n2, ' ]=', theor
					end if
				case(5)
					n2 = n(2)
					n3 = n(3)
					if (useKroneckers) then
						theor  = 3*( (1+n0) * (2+n1) * (3+n2) * n3   -   n0 * (n1-1) * (n2-2) * (1+n3) )
						theor2 = 1*( (1+n0) * (n1-2) * (n2-1) * n3   -   n0*(3+n1) * (2+n2) * (1+n3) )
					else
						theor  = 3*( (1+n0) * (1+n1) * (1+n2) * n3   -   n0 * n1 * n2 * (1+n3) )
						theor2 = 1*( (1+n0) * n1 * n2 * n3   -   n0*(1+n1) * (1+n2) * (1+n3) )
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F8.3)') &
&							'  theor=3*[', 1+n0, ' *', 1+n1, ' *', 1+n2, ' *', n3
						write(debugFileNumber, '(A,F8.3,A,F8.3,A,F8.3,A,F8.3,A,F10.6)') &
&							' -', n0, '*', n1, '*', n2, '*', 1+n3, ' ]=', theor
					end if
				case(6)
					n2 = n(2)
					n3 = n(3)
					if (useKroneckers) then
						theor2 = 3*( (1+n0) * (2+n1) * (3+n2) * n3   -   n0 * (n1-1) * (n2-2) * (1+n3) )
						theor  = 1*( (1+n0) * (n1-2) * (n2-1) * n3   -   n0*(3+n1) * (2+n2) * (1+n3) )
					else
						theor2 = 3*( (1+n0) * (1+n1) * (1+n2) * n3   -   n0 * n1 * n2 * (1+n3) )
						theor  = 1*( (1+n0) * n1 * n2 * n3   -   n0*(1+n1) * (1+n2) * (1+n3) )
					end if
			end select
			
			call qtai   % setCurrentSubsystemIndex( ss )
			call qtgs   % setCurrentSubsystemIndex( ss )
			call qtabcd % setCurrentSubsystemIndex( ss )
			
			v0 = theor
			v1 = qtai   % getValue(n, k, excludePosition)
			v2 = qtgs   % getValue(n, k, excludePosition)
			v3 = qtabcd % getValue(n, k, excludePosition)
!			print*, ' abs(v0-v1,v1-v2)=', abs(v0-v1), abs(v1-v2)
			write(debugFileNumber, '(A,4F12.6)' ) 'values: v0/v1/v2/v3=', v0, v1, v2, v3
			write(debugFileNumber, *) 'difference: abs(dv)=', abs(v0-v1), abs(v1-v2), abs(v2-v3)
			if (verbose) write(*, '(A,4F12.6)' ) '   values: v0/v1/v2/v3=', v0, v1, v2, v3
			if (verbose) write(*, *) '   difference: abs(dv)=', abs(v0-v1), abs(v1-v2), abs(v2-v3)
			
			good = ((abs(v0-v1).LT.1d-14).AND.(abs(v1-v2).LT.1d-14).AND.(abs(v2-v3).LT.1d-14))
			if (.NOT.good) then

				oH => OutputHTML( 'testQbeTerms.html', debugFileNumber+1 )	
				call hamilt % show( oH )

				call oH%beginEquation
				call qtai % show( oH )
				call oH%endEquation

				call oH%beginEquation
				call qtgs % show( oH )
				call oH%endEquation

				call oH%beginEquation
				call qtabcd % show( oH )
				call oH%endEquation
				deallocate( oH )

				oT => OutputTeX( 'testQbeTerms.tex', debugFileNumber+2 )
				call hamilt % show( oT )

				call oT%beginEquation
				call qtai % show( oT )
				call oT%endEquation

				call oT%beginEquation
				call qtgs % show( oT )
				call oT%endEquation

				call oT%beginEquation
				call qtabcd % show( oT )
				call oT%endEquation
				deallocate( oT )

				stop 'see details in testQbeTerms.log/html/tex'
			end if
			deallocate( hamilt )
		end do
		deallocate( n )
		deallocate( k )
		deallocate( qtabcd )
		deallocate( qtgs )
		deallocate( qtai )
	end do
  end do
	deallocate( v_13 )
	deallocate( v_31 )
	deallocate( v_12 )
	deallocate( v_21 )
	deallocate( v_bb )
	deallocate( v_ee )
	deallocate( sf2 )
	deallocate( sb2 )
	deallocate( sf )
	deallocate( sb )
	deallocate( dr )

!	print*, '---- ALL OK ----'

	call printTwice_closeFile( ' OK ]' )
	
end subroutine testQbeTerms

