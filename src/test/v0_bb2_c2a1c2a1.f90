	subroutine v0_bb2_c2a1c2a1(L, n_k, dn_dt, V0, rta, useKroneckers)
	use omp_lib
	integer L
	logical rta, useKroneckers, V_q_dependent

    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:1) :: n_k
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1, 0:1), intent(OUT) :: dn_dt

	integer L2, Lcube, m1xyzs
    integer m1x,m1y,m1z, s
    integer m3x,m3y,m3z
    integer m2x,m2y,m2z, m12x,m12y,m12z, d12
    integer m4x,m4y,m4z, m34x,m34y,m34z, d34
    integer m5x,m5y,m5z, m125x,m125y,m125z
    integer m6x,m6y,m6z, d56, d16, d15, id
    real*8 e, term, n1, n2, n3, n4, n5, n6, e1, e2, e3, e4, e5, e6, v, V0, V0_sqr

	V0_sqr = V0**2

	L2=L/2
	dn_dt(:,:,:,:) = 0.d0

	! two terms: 
	!	AA BB B A : k1+k2+p3+p4 = p5+k6  -> k6 = (k1+k2) + (p3+p4) - p5
	!	A B BB AA : k1+p2 = p3+p4+k5+k6  -> k6 = (k1+p2) - (p3+p4) - k5
	!
	! we see k1, p3, p4 in both terms -> n1,n3,n4 read early
	!
	
	Lcube = L*L*L	
	!$omp parallel do private(m1x,m1y,m1z, s, e1,n1, v, m2x,m2y,m2z, m12x,m12y,m12z, e2,d12, &
	!$omp  & m3x,m3y,m3z,e3,n3, m4x,m4y,m4z, m34x,m34y,m34z, e4,n4,d34, d15, id, &
	!$omp  & m5x,m5y,m5z, m125x,m125y,m125z, e5, m6x,m6y,m6z, e6, n2,n5,n6, term, d56, d16, m1xyzs )
	do m1xyzs=0, L*L*L*2-1
		m1x = MOD(m1xyzs, L) - L2
		m1y = MOD(m1xyzs/L, L) - L2
		m1z = MOD(m1xyzs/(L*L), L) - L2
		s = m1xyzs/Lcube

		id = omp_get_thread_num()
		if (id.EQ.0) then
			write(*, '(I0,A,I0,A)', advance='no') m1xyzs, '/', Lcube*2/omp_get_max_threads(), CHAR(13)
		end if
!	do m1z = -L2, L2-1
!	do m1y = -L2, L2-1
!	do m1x = -L2, L2-1
!	do s=0,1
		e1 = e_k(m1x,m1y,m1z)
		n1 = n_k(m1x,m1y,m1z, s)

		v = 0.d0
        do m2z = -L2, L2-1
			m12z = k_add( m1z, m2z, L)
	    do m2y = -L2, L2-1
			m12y = k_add( m1y, m2y, L)
	   	do m2x = -L2, L2-1
			m12x = k_add( m1x, m2x, L)

			e2 = e_k(m2x,m2y,m2z)

			d12 = 0
			if ( (m1x.EQ.m2x) .AND. (m1y.EQ.m2y) .AND. (m1z.EQ.m2z) ) d12 = 1

			do m3z = -L2, L2-1
			do m3y = -L2, L2-1
			do m3x = -L2, L2-1
				e3 = e_k(m3x,m3y,m3z)
				n3 = n_k(m3x,m3y,m3z,1-s)

				do m4z = -L2, L2-1
					m34z = k_add( m3z, m4z, L)
				do m4y = -L2, L2-1
					m34y = k_add( m3y, m4y, L)
				do m4x = -L2, L2-1
					m34x = k_add( m3x, m4x, L)

	                e4 = e_k(m4x,m4y,m4z)
					n4 = n_k(m4x,m4y,m4z,1-s)

					d34 = 0
					if ( (m3x.EQ.m4x) .AND. (m3y.EQ.m4y) .AND. (m3z.EQ.m4z) ) d34 = 1
	                
					do m5z = -L2, L2-1
						m125z = k_sub( m12z, m5z, L)
					do m5y = -L2, L2-1
						m125y = k_sub( m12y, m5y, L)
					do m5x = -L2, L2-1
						m125x = k_sub( m12x, m5x, L)
		                e5 = e_k(m5x,m5y,m5z)
	                
	! 1st term
	!	AA BB B A : k1+k2+p3+p4 = p5+k6  -> k6 = (k1+k2) + (p3+p4) - p5
						m6z = k_add( m125z, m34z, L)
						m6y = k_add( m125y, m34y, L)
						m6x = k_add( m125x, m34x, L)
						e6 = e_k(m6x,m6y,m6z)
	                
						if (e1+e2+e3+e4 .EQ. e5+e6) then
	
							n2 = n_k(m2x,m2y,m2z,s)
							n5 = n_k(m5x,m5y,m5z,1-s)
							n6 = n_k(m6x,m6y,m6z,s)

							d16 = 0
							if ((m1x.EQ.m6x) .AND. (m1y.EQ.m6y) .AND. (m1z.EQ.m6z)) d16 = 1
							
							if (.NOT.rta) then
								if (useKroneckers) then
									term = 2 * ( &
&										(n1+1)*(n2+1+d12)*(n3+1)*(n4+1+d34)*n5*n6 &
&										 - n1*(n2-d12)*n3*(n4-d34)*(n5+1)*(n6+1) )
								else
									term = 2 * ( &
&										(n1+1)*(n2+1)*(n3+1)*(n4+1)*n5*n6 &
&										 - n1*n2*n3*n4*(n5+1)*(n6+1) )
								end if
							else
								if (useKroneckers) then
									term = 2 * ( &
&										1*(n2+1+d12)*(n3+1)*(n4+1+d34)*n5*n6 &
&										 - 1*(n2-d12)*n3*(n4-d34)*(n5+1)*(n6+1) )
									if (d12.EQ.1) term = term + 2 * ( &
&										(n1+1)*1*(n3+1)*(n4+1+d34)*n5*n6 &
&										 - n1*1*n3*(n4-d34)*(n5+1)*(n6+1) )
									if (d16.EQ.1) term = term + 2 * ( &
&										(n1+1)*(n2+1+d12)*(n3+1)*(n4+1+d34)*n5*1 &
&										 - n1*(n2-d12)*n3*(n4-d34)*(n5+1)*1 )
								else
									term = 2 * ( &
&										1*(n2+1)*(n3+1)*(n4+1)*n5*n6 &
&										 - 1*n2*n3*n4*(n5+1)*(n6+1) )
									if (d12.EQ.1) term = term + 2 * ( &
&										(n1+1)*1*(n3+1)*(n4+1)*n5*n6 &
&										 - n1*1*n3*n4*(n5+1)*(n6+1) )
									if (d16.EQ.1) term = term + 2 * ( &
&										(n1+1)*(n2+1)*(n3+1)*(n4+1)*n5*1 &
&										 - n1*n2*n3*n4*(n5+1)*1 )
								end if
							end if
							
							v = v + term * V0_sqr
						end if

	! 2nd term
	!	A B BB AA : k1+p2 = p3+p4+k5+k6  -> k6 = (k1+p2) - (p3+p4) - k5
						m6z = k_sub( m125z, m34z, L)
						m6y = k_sub( m125y, m34y, L)
						m6x = k_sub( m125x, m34x, L)
						e6 = e_k(m6x,m6y,m6z)
	                
						if (e1+e2 .EQ. e3+e4+e5+e6) then
	
							d56 = 0
							if ( (m5x.EQ.m6x) .AND. (m5y.EQ.m6y) .AND. (m5z.EQ.m6z) ) d56 = 1
							d15 = 0
							if ( (m1x.EQ.m5x) .AND. (m1y.EQ.m5y) .AND. (m1z.EQ.m5z) ) d15 = 1
							d16 = 0
							if ( (m1x.EQ.m6x) .AND. (m1y.EQ.m6y) .AND. (m1z.EQ.m6z) ) d16 = 1
							
							n2 = n_k(m2x,m2y,m2z,1-s)
							n5 = n_k(m5x,m5y,m5z,s)
							n6 = n_k(m6x,m6y,m6z,s)
							if (.NOT.rta) then
								if (useKroneckers) then
									term = ( &
&										(n1+1)*(n2+1)*n3*(n4-d34)*n5*(n6-d56) &
&										 - n1*n2*(n3+1)*(n4+1+d34)*(n5+1)*(n6+1+d56) )
								else
									term = ( &
&										(n1+1)*(n2+1)*n3*n4*n5*n6 &
&										 - n1*n2*(n3+1)*(n4+1)*(n5+1)*(n6+1) )
								end if
							else
								if (useKroneckers) then
									term = ( &
&										1*(n2+1)*n3*(n4-d34)*n5*(n6-d56) &
&										 - 1*n2*(n3+1)*(n4+1+d34)*(n5+1)*(n6+1+d56) )
									if (d15.EQ.1) term = term + ( &
&										(n1+1)*(n2+1)*n3*(n4-d34)*1*(n6-d56) &
&										 - n1*n2*(n3+1)*(n4+1+d34)*1*(n6+1+d56) )
									if (d16.EQ.1) term = term + ( &
&										(n1+1)*(n2+1)*n3*(n4-d34)*n5*1 &
&										 - n1*n2*(n3+1)*(n4+1+d34)*(n5+1)*1 )
								else
									term = ( &
&										1*(n2+1)*n3*n4*n5*n6 &
&										 - 1*n2*(n3+1)*(n4+1)*(n5+1)*(n6+1) )
									if (d15.EQ.1) term = term + ( &
&										(n1+1)*(n2+1)*n3*n4*1*n6 &
&										 - n1*n2*(n3+1)*(n4+1)*1*(n6+1) )
									if (d16.EQ.1) term = term + ( &
&										(n1+1)*(n2+1)*n3*n4*n5*1 &
&										 - n1*n2*(n3+1)*(n4+1)*(n5+1)*1 )
								end if
							end if
							v = v + term * V0_sqr
						end if
					end do 
					end do
					end do ! m5
				end do
				end do
				end do ! m4
			end do
			end do
			end do ! m2
		end do
		end do
		end do ! m3
		dn_dt( m1x, m1y, m1z, s) = v
!	end do ! s : subsystem(0/1)
!	end do
!	end do
    end do ! m1
	!$omp end parallel do

	end subroutine v0_bb2_c2a1c2a1

