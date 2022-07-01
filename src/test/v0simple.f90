module v0simple

    implicit none

contains

	function k_add( k1, k2, L ) result(r)
	integer k1, k2, L, r

		r = k1 + k2
	    if (r .GE. L/2)  r = r - L
		if (r .LT. -L/2) r = r + L

	end function k_add

	function k_sub( k1, k2, L ) result(r)
	integer k1, k2, L, r

		r = k1 - k2
	    if (r .GE. L/2)  r = r - L
		if (r .LT. -L/2) r = r + L

	end function k_sub

	function e_k( x,y,z ) result(e)
	integer x,y,z, e

		e = x*x + y*y + z*z
	
	end function e_k

	function e_k_multi( x,y,z, s ) result(e)
	integer x,y,z, s, e

		e = (x*x + y*y + z*z)*(1+s)
	
	end function e_k_multi

!-------------- service functions ------------------------

	subroutine symm2asymm( L,  symm, asymm )
	integer L
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1) :: symm
    real*8, dimension(0:L-1, 0:L-1, 0:L-1), intent(OUT) :: asymm
	
	integer x,y,z, sx,sy,sz, L2
	
		L2 = L/2
		do sz=-L2, L2-1
			z = MOD(sz + L, L)
			do sy=-L2, L2-1
				y = MOD(sy + L, L)
				do sx=-L2, L2-1
					x = MOD(sx + L, L)
					
					asymm(x,y,z) = symm(sx,sy,sz)
				end do
			end do
		end do

	end subroutine symm2asymm

	subroutine asymm2symm( L,  asymm, symm )
	integer L
    real*8, dimension(0:L-1, 0:L-1, 0:L-1) :: asymm
    real*8, dimension(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1), intent(OUT) :: symm
	
	integer x,y,z, sx,sy,sz, L2
	
		L2 = L/2
		do sz=-L2, L2-1
			z = MOD(sz + L, L)
			do sy=-L2, L2-1
				y = MOD(sy + L, L)
				do sx=-L2, L2-1
					x = MOD(sx + L, L)
					
					symm(sx,sy,sz) = asymm(x,y,z)
				end do
			end do
		end do

	end subroutine asymm2symm

!---------------- functions for various interactions --------------------

	include 'v0_b_c2a2.f90'
	include 'v0_b_c2a1.f90'
	include 'v0_b_c3a1.f90'
	include 'v0_f_c2a2.f90'
	include 'v0_bb2_c1a1c1a1.f90'
	include 'v0_bb2_c2a1c2a1.f90'
	include 'v0_ff2_c1a1c1a1.f90'


end module v0simple
