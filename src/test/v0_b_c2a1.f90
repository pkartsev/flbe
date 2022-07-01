	subroutine v0_b_c2a1(L, n_k, dn_dt, V0, rta, useKroneckers)
	integer L
	logical rta, useKroneckers

    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0) :: n_k
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0), intent(OUT) :: dn_dt

	integer L2
    integer mx,my,mz
    integer m1x,m1y,m1z
    integer m2x,m2y,m2z, m12x,m12y,m12z, m23x,m23y,m23z
    integer m3x,m3y,m3z
    integer factor, d12,d23,d13
    real*8 e, term, n1, n2, n3, e1, e2, e3, occ, arg, v, V0, V0_sqr

	V0_sqr = V0**2
	
	L2=L/2
	dn_dt(:,:,:,:) = 0.d0

	!$omp parallel do private(m1x,m1y,m1z,n1,e1, v, m3x,m3y,m3z,n3,e3, &
	!$omp  &  m2x,m2y,m2z, e2, n2, d12,d23,d13,term )
	do m1z = -L2, L2-1
    do m1y = -L2, L2-1
	do m1x = -L2, L2-1
		n1 = n_k(m1x,m1y,m1z, 0)
		e1 = e_k(m1x,m1y,m1z)

		v = 0.d0
		do m3z = -L2, L2-1
		do m3y = -L2, L2-1
		do m3x = -L2, L2-1
			n3 = n_k(m3x,m3y,m3z, 0)
			e3 = e_k(m3x,m3y,m3z)

			d13 = 0
			if ( (m1x.EQ.m3x) .AND. (m1y.EQ.m3y) .AND. (m1z.EQ.m3z) ) d13 = 1
			
			!  a+a+ a : k1+k2 = k3 -> k2 = k3-k1
			m2z = k_sub( m3z, m1z, L)
			m2y = k_sub( m3y, m1y, L)
			m2x = k_sub( m3x, m1x, L)
			e2 = e_k(m2x,m2y,m2z)
			if (e1+e2 .EQ. e3) then
				n2 = n_k(m2x,m2y,m2z, 0)

				d12 = 0
				if ( (m1x.EQ.m2x) .AND. (m1y.EQ.m2y) .AND. (m1z.EQ.m2z) ) d12 = 1

				d23 = 0
				if ( (m2x.EQ.m3x) .AND. (m2y.EQ.m3y) .AND. (m2z.EQ.m3z) ) d23 = 1

				if (rta) then
					if (useKroneckers) then
						term = 2*( 1*(n2+1+d12)*n3 - 1*(n2-d12)*(n3+1) )
						if (d12.EQ.1) term = term + 2*( (n1+1)*1*n3 - n1*1*(n3+1) )
						if (d13.EQ.1) term = term + 2*( (n1+1)*(n2+1+d12)*1 - n1*(n2-d12)*1 )
					else
						term = 2*( 1*(n2+1)*n3 - 1*n2*(n3+1) )
						if (d12.EQ.1) term = term + 2*( (n1+1)*1*n3 - n1*1*(n3+1) )
						if (d13.EQ.1) term = term + 2*( (n1+1)*(n2+1)*1 - n1*n2*1 )
					end if
				else
					if (useKroneckers) then
						term = 2*( (n1+1)*(n2+1+d12)*n3 - n1*(n2-d12)*(n3+1) )
					else
						term = 2*( (n1+1)*(n2+1)*n3 - n1*n2*(n3+1) )
					end if
				end if
				v = v + term * V0_sqr
			end if
							
			! (H.c.) : a+ a a : k1 = k2+k3 -> k2 = k1-k3
			m2z = k_sub( m1z, m3z, L)
			m2y = k_sub( m1y, m3y, L)
			m2x = k_sub( m1x, m3x, L)
			e2 = e_k(m2x,m2y,m2z)
			if (e1 .EQ. e2+e3) then
				n2 = n_k(m2x,m2y,m2z, 0)

				d12 = 0
				if ( (m1x.EQ.m2x) .AND. (m1y.EQ.m2y) .AND. (m1z.EQ.m2z) ) d12 = 1

				d23 = 0
				if ( (m2x.EQ.m3x) .AND. (m2y.EQ.m3y) .AND. (m2z.EQ.m3z) ) d23 = 1

				if (rta) then
					if (useKroneckers) then
						term = ( 1*n2*(n3-d23) - 1*(n2+1)*(n3+1+d23) )
						if (d12.EQ.1) term = term + ( (n1+1)*1*(n3-d23) - n1*1*(n3+1+d23) )
						if (d13.EQ.1) term = term + ( (n1+1)*n2*1 - n1*(n2+1)*1 )
					else
						term = ( 1*n2*n3 - 1*(n2+1)*(n3+1) )
						if (d12.EQ.1) term = term + ( (n1+1)*1*n3 - n1*1*(n3+1) )
						if (d13.EQ.1) term = term + ( (n1+1)*n2*1 - n1*(n2+1)*1 )
					end if
				else
					if (useKroneckers) then
						term = ( (n1+1)*n2*(n3-d23) - n1*(n2+1)*(n3+1+d23) )
					else
						term = ( (n1+1)*n2*n3 - n1*(n2+1)*(n3+1) )
					end if
				end if
				v = v + term * V0_sqr
			end if
		end do
		end do
		end do
		dn_dt( m1x, m1y, m1z, 0) = v
	end do
	end do
    end do
	!$omp end parallel do

	end subroutine v0_b_c2a1

