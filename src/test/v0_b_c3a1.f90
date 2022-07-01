	subroutine v0_b_c3a1(L, n_k, dn_dt, V0, rta, useKroneckers)
	integer L
	logical rta, useKroneckers

    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0) :: n_k
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0), intent(OUT) :: dn_dt

	integer L2
    integer mx,my,mz
    integer m1x,m1y,m1z
    integer m2x,m2y,m2z, m12x,m12y,m12z, m23x,m23y,m23z
    integer m3x,m3y,m3z
    integer m4x,m4y,m4z
    integer factor, d12,d23,d13, d14,d24,d34
    real*8 e, term, n1, n2, n3, n4, e1, e2, e3, e4, occ, arg, v, V0, V0_sqr

	V0_sqr = V0**2
	
	L2=L/2
	dn_dt(:,:,:,:) = 0.d0

	!$omp parallel do private(m1x,m1y,m1z,n1,e1, v, m2x,m2y,m2z,n2,e2, d12, &
	!$omp  &  m3x,m3y,m3z,n3,e3, d13,d23, m23x,m23y,m23z, m4x,m4y,m4z, e4,n4, d14,d24,d34,term )
	do m1z = -L2, L2-1
    do m1y = -L2, L2-1
	do m1x = -L2, L2-1
		n1 = n_k(m1x,m1y,m1z, 0)
		e1 = e_k(m1x,m1y,m1z)

		v = 0.d0

		do m2z = -L2, L2-1
		do m2y = -L2, L2-1
		do m2x = -L2, L2-1
			n2 = n_k(m2x,m2y,m2z, 0)
			e2 = e_k(m2x,m2y,m2z)

			d12 = 0
			if ( (m1x.EQ.m2x) .AND. (m1y.EQ.m2y) .AND. (m1z.EQ.m2z) ) d12 = 1

			do m3z = -L2, L2-1
			do m3y = -L2, L2-1
			do m3x = -L2, L2-1
				n3 = n_k(m3x,m3y,m3z, 0)
				e3 = e_k(m3x,m3y,m3z)

				d13 = 0
				if ( (m1x.EQ.m3x) .AND. (m1y.EQ.m3y) .AND. (m1z.EQ.m3z) ) d13 = 1

				d23 = 0
				if ( (m2x.EQ.m3x) .AND. (m2y.EQ.m3y) .AND. (m2z.EQ.m3z) ) d23 = 1
			
				m23z = k_add( m2z, m3z, L)
				m23y = k_add( m2y, m3y, L)
				m23x = k_add( m2x, m3x, L)

				!  a+a+a+ a : k1+k2+k3 = k4 -> k4 = k1+k2+k3
				m4z = k_add( m1z, m23z, L)
				m4y = k_add( m1y, m23y, L)
				m4x = k_add( m1x, m23x, L)
				e4 = e_k(m4x,m4y,m4z)
				if (e1+e2+e3 .EQ. e4) then

					n4 = n_k(m4x,m4y,m4z, 0)

					d14 = 0
					if ( (m1x.EQ.m4x) .AND. (m1y.EQ.m4y) .AND. (m1z.EQ.m4z) ) d14 = 1

					if (rta) then
						if (useKroneckers) then
							term = 3*( 1*(n2+1+d12)*(n3+1+d13+d23)*n4 - 1*(n2-d12)*(n3-d13-d23)*(n4+1) )
							if (d12.EQ.1) term = term + &
&								3*( (n1+1)*1*(n3+1+d13+d23)*n4 - n1*1*(n3-d13-d23)*(n4+1) )
							if (d13.EQ.1) term = term + &
&								3*( (n1+1)*(n2+1+d12)*1*n4 - n1*(n2-d12)*1*(n4+1) )
							if (d14.EQ.1) term = term + &
&								3*( (n1+1)*(n2+1+d12)*(n3+1+d13+d23)*1 - n1*(n2-d12)*(n3-d13-d23)*1 )
						else
							term = 3*( 1*(n2+1)*(n3+1)*n4 - 1*n2*n3*(n4+1) )
							if (d12.EQ.1) term = term + &
&								3*( (n1+1)*1*(n3+1)*n4 - n1*1*n3*(n4+1) )
							if (d13.EQ.1) term = term + &
&								3*( (n1+1)*(n2+1)*1*n4 - n1*n2*1*(n4+1) )
							if (d14.EQ.1) term = term + &
&								3*( (n1+1)*(n2+1)*(n3+1)*1 - n1*n2*n3*1 )
						end if
					else
						if (useKroneckers) then
							term = 3*( (n1+1)*(n2+1+d12)*(n3+1+d13+d23)*n4 - n1*(n2-d12)*(n3-d13-d23)*(n4+1) )
						else
							term = 3*( (n1+1)*(n2+1)*(n3+1)*n4 - n1*n2*n3*(n4+1) )
						end if
					end if
					v = v + term * V0_sqr
				end if
							
				! (H.c.) : a+ a a a : k1 = k2+k3+k4 -> k4 = k1-(k2+k3)
				m4z = k_sub( m1z, m23z, L)
				m4y = k_sub( m1y, m23y, L)
				m4x = k_sub( m1x, m23x, L)
				e4 = e_k(m4x,m4y,m4z)
				if (e1 .EQ. e2+e3+e4) then
!					n4 = 0.d0
					n4 = n_k(m4x,m4y,m4z, 0)
!					n4 = n_k(m2x,m2y,m2z, 0)

					d14 = 0
					if ( (m1x.EQ.m4x) .AND. (m1y.EQ.m4y) .AND. (m1z.EQ.m4z) ) d14 = 1

					d24 = 0
					if ( (m2x.EQ.m4x) .AND. (m2y.EQ.m4y) .AND. (m2z.EQ.m4z) ) d24 = 1

					d34 = 0
					if ( (m3x.EQ.m4x) .AND. (m3y.EQ.m4y) .AND. (m3z.EQ.m4z) ) d34 = 1

					if (rta) then
						if (useKroneckers) then
							term = ( 1*n2*(n3-d23)*(n4-d24-d34) - 1*(n2+1)*(n3+1+d23)*(n4+1+d24+d34) )
							if (d12.EQ.1) term = term + &
&								( (n1+1)*1*(n3-d23)*(n4-d24-d34) - n1*1*(n3+1+d23)*(n4+1+d24+d34) )
							if (d13.EQ.1) term = term + &
&								( (n1+1)*n2*1*(n4-d24-d34) - n1*(n2+1)*1*(n4+1+d24+d34) )
							if (d14.EQ.1) term = term + &
&								( (n1+1)*n2*(n3-d23)*1 - n1*(n2+1)*(n3+1+d23)*1 )
						else
							term = ( 1*n2*n3*n4 - 1*(n2+1)*(n3+1)*(n4+1) )
							if (d12.EQ.1) term = term + &
&								( (n1+1)*1*n3*n4 - n1*1*(n3+1)*(n4+1) )
							if (d13.EQ.1) term = term + &
&								( (n1+1)*n2*1*n4 - n1*(n2+1)*1*(n4+1) )
							if (d14.EQ.1) term = term + &
&								( (n1+1)*n2*n3*1 - n1*(n2+1)*(n3+1)*1 )
						end if
					else
						if (useKroneckers) then
							term = ( (n1+1)*n2*(n3-d23)*(n4-d24-d34) - n1*(n2+1)*(n3+1+d23)*(n4+1+d24+d34) )
						else
							term = ( (n1+1)*n2*n3*n4 - n1*(n2+1)*(n3+1)*(n4+1) )
						end if
					end if
					v = v + term * V0_sqr
				end if
			end do
			end do
			end do ! m3 (x,y,z)
		end do
		end do
		end do ! m2 (x,y,z)
		dn_dt( m1x, m1y, m1z, 0) = v
	end do
	end do
    end do ! m1 (x,y,z)
	!$omp end parallel do

	end subroutine v0_b_c3a1

