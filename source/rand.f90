!==============================================================
!
double precision function rand(ix)
!     function rand(ix)
!
!==============================================================
!     portable random number generator using the recursion
!     ix=ix*a mod p
!     implicit real*8(a-h,o-z)
!
  implicit none(type, external)
  integer a, p, ix, b15, b16, xhi, xalo, leftlo, fhi, k
!
!     7**5, 2**15, 2**16, 2**31-1
  data a/16807/, b15/32768/, b16/65536/, p/2147483647/
!
!     get 15 high order bits of ix
  xhi = ix/b16
!     get low order bits into ix and form low order product
  xalo = (ix - xhi*b16)*a
!     get 15 high order bits of lower order product
  leftlo = xalo/b16
!     form the 31 highest bits of full product
  fhi = xhi*a + leftlo
!     get overflow past 31st bit of full product
  k = fhi/b15
!     assemble all the parts and presubtract p
!     the parentheses are all essential
  ix = (((xalo - leftlo*b16) - p) + (fhi - k*b15)*b16) + k
!     add p back in if necessary
  if (ix .lt. 0) ix = ix + p
!     multiply by 1/(2**31-1)
  rand = dble(ix)/dble(2.**31 - 1)
  return
end function rand
