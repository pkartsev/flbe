subroutine testDispersionRelation

use flbe
implicit none

type(Geometry), pointer :: geom
type(DispersionRelation), pointer :: dr0, dr1, dr2, drTable

integer kx,ky,kz, L, k(1:3), ksym(1:3)
double precision e, coeff0, coeff1, coeff2, e0
real r
logical bad

	write(*, '(A)', advance='no') ' [ test DispersionRelation...'
	L=8
	geom => Geometry(L,L,L)

	call RANDOM_NUMBER(r)
	coeff0 = r + 1.d0
	call RANDOM_NUMBER(r)
	coeff1 = r + 2.d0
	call RANDOM_NUMBER(r)
	coeff2 = r + 3.d0
!	coeff1 = 1.2345678d1
!	coeff2 = 7.654321d-1
!	coeff0 = 2.44244d0

	dr0 => DispersionRelationConst( coeff0 )
	dr1 => DispersionRelationLinear(  )
	dr2 => DispersionRelationParabolic(  )
	call dr1 % setCoeff( coeff1 )
	call dr2 % setCoeff( coeff2 )

	bad = .FALSE.
	
	do kz=0, L-1
	do ky=0, L-1
	do kx=0, L-1
		k = (/ kx, ky, kz /)
		call geom % symmetrize( k, ksym )
		
		e = dr2 % getEnergy( k, geom )
		e0 = coeff2 * sum( ksym(:)**2 )
		if (abs(e-e0).GT.1d-5) bad = .TRUE.
		
		e = dr1 % getEnergy( k, geom )
		e0 = coeff1 * sqrt( 0.d0 + sum( ksym(:)**2 ) )
		if (abs(e-e0).GT.1d-5) bad = .TRUE.
		
		e = dr0 % getEnergy( k, geom )
		e0 = coeff0
		if (abs(e-e0).GT.1d-5) bad = .TRUE.
		
	end do
	end do
	end do

	deallocate(dr2)
	deallocate(dr1)
	deallocate(dr0)
	deallocate(geom)

	if (bad) stop 'FAILED'
	print*, ' OK ]'

end subroutine testDispersionRelation
