	subroutine v0_ff2_c1a1c1a1(L, n_k, dn_dt, V_q, V_q_dependent, rta)
	integer L
	logical rta, V_q_dependent

    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:1) :: n_k
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:1), intent(OUT) :: dn_dt
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1) :: V_q

	integer L2
    integer mx,my,mz
    integer m1x,m1y,m1z
    integer m2x,m2y,m2z, m12x,m12y,m12z
    integer m3x,m3y,m3z
    integer m4x,m4y,m4z, m1xy, factor, s, qx,qy,qz
    real*8 e, term, n1, n2, n3, n4, e1, e2, e3, e4, occ, arg, v, V_sqr

	L2=L/2
	dn_dt(:,:,:,:) = 0.d0

	!$omp parallel do private(m1x,m1y,m1z,n1,e1, v, m2x,m2y,m2z,m12x,m12y,m12z,n2,e2, &
	!$omp  & m3x,m3y,m3z, m4x,m4y,m4z, n3,n4,e3,e4, term, s, qx,qy,qz, V_sqr )
	do m1z = -L2, L2-1
    do m1y = -L2, L2-1
	do m1x = -L2, L2-1
	do s=0,1
		n1 = n_k(m1x,m1y,m1z, s)
		e1 = e_k_multi(m1x,m1y,m1z,s)
!		if (s.EQ.1) e1 = 0.d0
		v = 0.d0
        do m2z = -L2, L2-1
			m12z = k_add( m1z, m2z, L)
	    do m2y = -L2, L2-1
			m12y = k_add( m1y, m2y, L)
	   	do m2x = -L2, L2-1
			m12x = k_add( m1x, m2x, L)
			n2 = n_k(m2x,m2y,m2z,1-s)
			e2 = e_k_multi(m2x,m2y,m2z, 1-s)
!			if ((1-s).EQ.1) e2 = 0.d0
			do m3z = -L2, L2-1
				m4z = k_sub( m12z, m3z, L)
			do m3y = -L2, L2-1
				m4y = k_sub( m12y, m3y, L)
			do m3x = -L2, L2-1
				m4x = k_sub( m12x, m3x, L)

				n3 = n_k(m3x,m3y,m3z,1-s)
                e3 = e_k_multi(m3x,m3y,m3z, 1-s)
!				if ((1-s).EQ.1) e3 = 0.d0
                e4 = e_k_multi(m4x,m4y,m4z, s)
!				if (s.EQ.1) e4 = 0.d0
				if (e1+e2 .EQ. e3+e4) then

					n4 = n_k(m4x,m4y,m4z,s)
					if (rta) then
						term = (-1)*(1-n2)*n3*n4 - 1*n2*(1-n3)*(1-n4)
						if ( (m1x.EQ.m4x) .AND. (m1y.EQ.m4y) .AND. (m1z.EQ.m4z) ) then
							term = term + (1-n1)*(1-n2)*n3*1 - n1*n2*(1-n3)*(-1)
						end if
					else
						term = (1-n1)*(1-n2)*n3*n4 - n1*n2*(1-n3)*(1-n4)
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
		dn_dt( m1x, m1y, m1z, s) = v
	end do ! s : subsystem(0/1)
	end do
	end do
    end do
	!$omp end parallel do

	end subroutine v0_ff2_c1a1c1a1

