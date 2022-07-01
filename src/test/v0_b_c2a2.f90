	subroutine v0_b_c2a2(L, n_k, dn_dt, V_q, V_q_dependent, rta, useKroneckers)
	integer L
	logical rta, useKroneckers, V_q_dependent

    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0) :: n_k
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:0), intent(OUT) :: dn_dt
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1) :: V_q

	integer L2
    integer mx,my,mz
    integer m1x,m1y,m1z
    integer m2x,m2y,m2z, m12x,m12y,m12z
    integer m3x,m3y,m3z, qx,qy,qz
    integer m4x,m4y,m4z, m1xy, factor, d12, d34, d13, d14
    real*8 e, term, n1, n2, n3, n4, e1, e2, e3, e4, occ, arg, v, V_sqr

	L2=L/2
	dn_dt(:,:,:,:) = 0.d0

	!$omp parallel do private(m1x,m1y,m1z,n1,e1, v, m2x,m2y,m2z,m12x,m12y,m12z,d12,n2,e2, &
	!$omp  & m3x,m3y,m3z, m4x,m4y,m4z, d34,n3,n4,e3,e4, term, d13,d14, qx,qy,qz, V_sqr )
	do m1z = -L2, L2-1
    do m1y = -L2, L2-1
	do m1x = -L2, L2-1
		n1 = n_k(m1x,m1y,m1z, 0)
		e1 = e_k(m1x,m1y,m1z)

		v = 0.d0
        do m2z = -L2, L2-1
			m12z = k_add( m1z, m2z, L)
	    do m2y = -L2, L2-1
			m12y = k_add( m1y, m2y, L)
	   	do m2x = -L2, L2-1
			m12x = k_add( m1x, m2x, L)
			
			d12 = 0
			if ( (m1x.EQ.m2x) .AND. (m1y.EQ.m2y) .AND. (m1z.EQ.m2z) ) d12 = 1

			n2 = n_k(m2x,m2y,m2z,0)
			e2 = e_k(m2x,m2y,m2z)

			do m3z = -L2, L2-1
				m4z = k_sub( m12z, m3z, L)
			do m3y = -L2, L2-1
				m4y = k_sub( m12y, m3y, L)
			do m3x = -L2, L2-1
				m4x = k_sub( m12x, m3x, L)

				d13 = 0
				if ( (m1x.EQ.m3x) .AND. (m1y.EQ.m3y) .AND. (m1z.EQ.m3z) ) d13 = 1

				d14 = 0
				if ( (m1x.EQ.m4x) .AND. (m1y.EQ.m4y) .AND. (m1z.EQ.m4z) ) d14 = 1

				d34 = 0
				if ( (m3x.EQ.m4x) .AND. (m3y.EQ.m4y) .AND. (m3z.EQ.m4z) ) d34 = 1

				n3 = n_k(m3x,m3y,m3z,0)
                e3 = e_k(m3x,m3y,m3z)
                e4 = e_k(m4x,m4y,m4z)
				if (e1+e2 .EQ. e3+e4) then

					n4 = n_k(m4x,m4y,m4z,0)
					term = 0.d0
					if (.NOT.rta) then ! dn/dt
						if (useKroneckers) then
							term = 2*( (n1+1)*(n2+1+d12)*n3*(n4-d34) - n1*(n2-d12)*(n3+1)*(n4+1+d34) ) 
						else
							term = 2*( (n1+1)*(n2+1)*n3*n4 - n1*n2*(n3+1)*(n4+1) )
						end if
					else ! RTA
						if (useKroneckers) then
							term = 2*( 1*(n2+1+d12)*n3*(n4-d34) - 1*(n2-d12)*(n3+1)*(n4+1+d34) ) 
! test
							if (d12.EQ.1) term = term + 2*( (n1+1)*1*n3*(n4-d34) - n1*1*(n3+1)*(n4+1+d34) ) 
							if (d13.EQ.1) term = term + 2*( (n1+1)*(n2+1+d12)*1*(n4-d34) - n1*(n2-d12)*1*(n4+1+d34) ) 
							if (d14.EQ.1) term = term + 2*( (n1+1)*(n2+1+d12)*n3*1 - n1*(n2-d12)*(n3+1)*1 ) 
						else
							term = 2*( 1*(n2+1)*n3*n4 - 1*n2*(n3+1)*(n4+1) ) 
! test
							if (d12.EQ.1) term = term + 2*( (n1+1)*1*n3*n4 - n1*1*(n3+1)*(n4+1) ) 
							if (d13.EQ.1) term = term + 2*( (n1+1)*(n2+1)*1*n4 - n1*(n2)*1*(n4+1) ) 
							if (d14.EQ.1) term = term + 2*( (n1+1)*(n2+1)*n3*1 - n1*n2*(n3+1)*1 ) 
						end if
					end if
					if (V_q_dependent) then
						qx = k_sub( m1x, m4x, L )
						qy = k_sub( m1y, m4y, L )
						qz = k_sub( m1z, m4z, L )
						V_sqr = V_q( qx,qy,qz )**2
					else
						V_sqr = V_q( 0,0,0 )**2
					end if
					v = v + term * V_sqr
				end if
			end do
			end do
			end do
		end do
		end do
		end do
		dn_dt( m1x, m1y, m1z, 0) = v
	end do
	end do
    end do
	!$omp end parallel do

	end subroutine v0_b_c2a2

