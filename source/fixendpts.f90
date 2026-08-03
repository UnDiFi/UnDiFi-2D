subroutine fixendpts(xy, n1, n2, ndim)

!     this routine does the following:
!     1) identifies the bndry segments cut by the shock lines
!     2) modifies the bndry structured (i.e. the pointer
!        ibndfac)
!        2a) to include the shock segments
!        2b) to account for the bndry segments of the background mesh
!            cut by the shock

  use mod_kinds, only: wp, i4
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) ndim, n1, n2

!     .. array arguments ..
  real(wp) xy(ndim, *)

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
  real(wp) x1, y1, x2, y2, dx, dy, x4, y4, s, t, help
  integer(i4) i, ifail

!     check the intersections of the shock lines
  x1 = xy(1, n1)
  y1 = xy(2, n1)
  x2 = xy(1, n2)
  y2 = xy(2, n2)
  dx = x2 - x1
  dy = y2 - y1
  xy(1, n1) = 0.d0
  xy(2, n1) = y1 - dy/dx*x1

  write (6, *) 'shock point ', n1, ' was at', x1, y1,&
  &' is now at', xy(1, n1), xy(2, n1)
!     call x04eaf('general',' ',3,nbfac,ibndfac,3,
!    +            'bndry pointer',ifail)
!
!
  return
end subroutine fixendpts
