! Compute an unsteady regular reflection point

module co_urr_ctx_m
! Phase 3.2 increment 1 (ROADMAP.md #14): closure payload for futp1's
! residual evaluation, replacing the (a, b, tauwx, tauwy) argument list
! co_urr used to pass on every call into futp1, now that both go
! through mod_newton_solve's generic residual_if(i, y, ctx) interface.
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: urr_ctx_t

  type :: urr_ctx_t
    real(wp) :: a(10, 15) = 0.0_wp
    real(wp) :: b(10) = 0.0_wp
    real(wp) :: tauwx = 0.0_wp, tauwy = 0.0_wp
  end type urr_ctx_t

end module co_urr_ctx_m

subroutine co_urr(y, tauwx, tauwy, yn1)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_newton_solve, only: newton_solve, residual_if
  use co_urr_ctx_m, only: urr_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: futp1
  include 'paramt.h'

  integer(i4) i, j, nn
  real(wp) gam, delta
  real(wp) tauwx, tauwy
  real(wp) y, yn1
  dimension y(15), yn1(15)
  type(urr_ctx_t) :: ctx

!     assign constants
  gam = ga
  delta = (gam - 1.0)/2.

  nn = 15

! build the matrix a and the vecotr b of the additional conditions
! (which varies according to the case)

!     compute oblique shock given the upsteam state and deviation angle
  do i = 1, nn - (4 + 1)
    do j = 1, nn
      ctx%a(i, j) = 0.d+0
    end do
    ctx%b(i) = 0.d+0
  end do

!     state 1 is known
  ctx%a(1, 1) = 1.0d+0
  ctx%a(2, 2) = 1.0d+0
  ctx%a(3, 3) = 1.0d+0
  ctx%a(4, 4) = 1.0d+0
  ctx%b(1) = y(1)
  ctx%b(2) = y(2)
  ctx%b(3) = y(3)
  ctx%b(4) = y(4)

!     state 2 is known
  ctx%a(5, 5) = 1.0d+0
  ctx%a(6, 6) = 1.0d+0
  ctx%a(7, 7) = 1.0d+0
  ctx%a(8, 8) = 1.0d+0
  ctx%b(5) = y(5)
  ctx%b(6) = y(6)
  ctx%b(7) = y(7)
  ctx%b(8) = y(8)

!     value of sx12 is known
  ctx%a(9, 13) = 1.0d+0
  ctx%b(9) = y(13)

!     velocity of the reflection point
  ctx%a(10, 15) = 1.0d+0
  ctx%b(10) = y(15)

  ctx%tauwx = tauwx
  ctx%tauwy = tauwy

! compute downstream state and shock velocity with the Newton-Raphson method
  call newton_solve(nn, y, futp1, ctx, 1.0_wp, 1.0e-6_wp, 0.01_wp, yn1)

end subroutine co_urr

real(wp) function futp1(i, y, ctx) result(r)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use co_urr_ctx_m, only: urr_ctx_t
  implicit none(type, external)
  include 'paramt.h'

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  integer(i4) ii, j
  real(wp) a(10, 15), b(10), tauwx, tauwy, wrr, wn
  real(wp) ro2, ro3, p2, p3, u2, u3, gam, delta
  real(wp) v2, v3, un2, un3, ut2, ut3
  real(wp) nwx, nwy
  real(wp) taux23, tauy23, nx23, ny23, sx23

  select type (ctx)
  type is (urr_ctx_t)
    a = ctx%a
    b = ctx%b
    tauwx = ctx%tauwx
    tauwy = ctx%tauwy
  end select

!     assign constants and variables
  gam = ga
  delta = (gam - 1.0)/2.

  ro2 = y(5)
  p2 = y(6)
  u2 = y(7)
  v2 = y(8)
  ro3 = y(9)
  p3 = y(10)
  u3 = y(11)
  v3 = y(12)
  sx23 = y(14)
  wrr = y(15)

!     compute normal and tangent versor
  taux23 = cos(sx23)
  tauy23 = sin(sx23)
  nx23 = -tauy23
  ny23 = taux23

!     compute the normal and tangential velocity components
  un2 = u2*nx23 + v2*ny23
  un3 = u3*nx23 + v3*ny23
  ut2 = u2*taux23 + v2*tauy23
  ut3 = u3*taux23 + v3*tauy23

!     compute the normal velocity component of the reflection point
  wn = wrr*(tauwx*nx23 + tauwy*ny23)

  r = 0.d0
  if (i .eq. 1) then
    r = ro2*un2 - ro3*un3 - wn*(ro2 - ro3)
  elseif (i .eq. 2) then
    r = p2 + ro2*(un2 - wn)**2 - p3 - ro3*(un3 - wn)**2
  elseif (i .eq. 3) then
    r = ut2 - ut3
  elseif (i .eq. 4) then
    r = gam/(gam - 1.0)*p2/ro2 + 0.5*(un2 - wn)**2&
    &- gam/(gam - 1.0)*p3/ro3 - 0.5*(un3 - wn)**2

  elseif (i .eq. 5) then
    nwx = -tauwy
    nwy = tauwx
    r = u3*nwx + v3*nwy

  elseif (i .ge. 6) then
    r = 0.d0
    ii = i - (1*4 + 1)
    do j = 1, 15
      r = r + a(ii, j)*y(j)
    end do
    r = r - b(ii)
  end if

end function futp1
