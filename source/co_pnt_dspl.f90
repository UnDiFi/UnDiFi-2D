! Compute displaced positions of the shock/discontinuity points

subroutine co_pnt_dspl(&                   ! not used
&xysh,&
&xyshu,&
&xyshd,&
&zroeshu,&
&nodcodsh,&
&vshnor,&
&nshocks,&
&nshockpoints,&
&nshockedges,&       ! not used
&typesh,&
&nspecpoints,&
&typespecpoints,&
&shinspps,&
&ispclr)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_unpack
  implicit none(type, external)
  include 'paramt.h'

!     .. scalar arguments ..
  integer(i4) nspecpoints, shinspps(2, 5, *), ispclr(5, *)
  integer(i4) nshocks, nshockedges(*), nshockpoints(*)
  character*1 typesh(*)
  character*5 typespecpoints(*)

!     .. array arguments ..
  real(wp) xysh(ndim, npshmax, *),&
  &xyshu(ndim, npshmax, *),&
  &xyshd(ndim, npshmax, *),&
  &vshnor(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *)

  integer(i4) nodcodsh(npshmax, *)

!     .. local scalars ..
  real(wp) dx, dy
  integer(i4) i, isppnts, ish
  class(special_point_t), allocatable :: sp

  do ish = 1, nshmax
    do i = 1, npshmax
      nodcodsh(i, ish) = -99
    end do
  end do

!     add a second layer of shock nodes
!     place extra shock nodes on the shock normal
  do ish = 1, nshocks
    do i = 1, nshockpoints(ish)
      dx = vshnor(1, i, ish)
      dy = vshnor(2, i, ish)
      xyshu(1, i, ish) = xysh(1, i, ish) + 0.5d0*eps*dx
      xyshu(2, i, ish) = xysh(2, i, ish) + 0.5d0*eps*dy
      xyshd(1, i, ish) = xysh(1, i, ish) - 0.5d0*eps*dx
      xyshd(2, i, ish) = xysh(2, i, ish) - 0.5d0*eps*dy
      nodcodsh(i, ish) = 10
    end do
  end do

!     correct the displacement in the boundary special points
!     dx/dy are intentionally NOT reset per point below -- they carry
!     forward across dispatch calls exactly as the pre-refactor
!     subroutine's own shared locals did (see mod_special_point.f90's
!     sp_displace_if for why regular_reflection_t's RRX path depends on
!     reading this exact stale value).
  do isppnts = 1, nspecpoints
    call make_special_point(typespecpoints(isppnts), sp)
    call sp_unpack(sp, shinspps(:, :, isppnts))
    call sp%displace(xysh(:, :, 1:nshmax), xyshu(:, :, 1:nshmax),&
    &xyshd(:, :, 1:nshmax), vshnor(:, :, 1:nshmax), zroeshu(:, :, 1:nshmax),&
    &nshockpoints(1:nshmax), dx, dy)
  end do

  return
end subroutine co_pnt_dspl

!     find interpolation point
subroutine co_intr_pnt(xi, yi, xc, yc, xs, ys)
  use mod_kinds, only: wp, i4

  integer(i4) nn
  parameter(nn=2)
  real(wp) a(nn, nn), b(nn), x(nn)
  real(wp) xi, yi
  real(wp) xc(2), yc(2), xs(2), ys(2)

  rdshp = -1.0

  a(1, 1) = (ys(2) - ys(1))
  a(1, 2) = (xs(1) - xs(2))
  b(1) = xs(2)*(ys(2) - ys(1)) + ys(2)*(xs(1) - xs(2))
  a(2, 1) = (yc(2) - yc(1))
  a(2, 2) = (xc(1) - xc(2))
  b(2) = xc(2)*(yc(2) - yc(1)) + yc(2)*(xc(1) - xc(2))

  call solg(nn, nn, a, b, x)

  xi = x(1)
  yi = x(2)

  return
end subroutine co_intr_pnt
