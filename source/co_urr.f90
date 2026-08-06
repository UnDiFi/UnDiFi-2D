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
  use mod_newton_solve, only: newton_solve, residual_if, jacobian_if
  use co_urr_ctx_m, only: urr_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: futp1
  procedure(jacobian_if) :: jfutp1
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
  call newton_solve(nn, y, futp1, ctx, 1.0_wp, 1.0e-6_wp, 0.01_wp, yn1, jac=jfutp1)

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

! ************************************
! Phase 3.5 increment 3 (ROADMAP.md #14): analytic Jacobian of futp1.
! Rows 1-5 are the genuinely nonlinear equations (mass/momentum/
! tangential-velocity/energy across the reflected shock, plus the
! reflection-point tangency condition); rows 6-15 are the linear
! "known-value" constraint rows built into ctx%a/ctx%b before
! newton_solve is ever called, so their Jacobian is exactly ctx%a's row
! -- no derivation needed there. un2/un3/ut2/ut3/wn are built from the
! same rotated (normal, tangent) basis (nx,ny)/(tx,ty) parameterized by
! sx23=y(14); differentiating that rotation gives the identities
! d(un)/d(sx23) = -ut and d(ut)/d(sx23) = un reused throughout below.
! Validated the same way as co_shock's jf / co_dc's jfdc (see
! mod_newton_solve.f90) before being trusted.
subroutine jfutp1(y, ctx, g)
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use co_urr_ctx_m, only: urr_ctx_t
  implicit none(type, external)
  include 'paramt.h'

  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx
  real(wp), intent(out) :: g(:, :)

  real(wp) :: a(10, 15), tauwx, tauwy
  real(wp) :: ro2, p2, u2, v2, ro3, p3, u3, v3, sx23, wrr
  real(wp) :: gam, c
  real(wp) :: nx, ny, tx, ty, un2, un3, ut2, ut3, wn, d2, d3, kk, kt
  integer(i4) :: i, j

  select type (ctx)
  type is (urr_ctx_t)
    a = ctx%a
    tauwx = ctx%tauwx
    tauwy = ctx%tauwy
  end select

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

  tx = cos(sx23)
  ty = sin(sx23)
  nx = -ty
  ny = tx

  un2 = u2*nx + v2*ny
  un3 = u3*nx + v3*ny
  ut2 = u2*tx + v2*ty
  ut3 = u3*tx + v3*ty
  kk = tauwx*nx + tauwy*ny ! wn = wrr*kk
  kt = tauwx*tx + tauwy*ty ! d(wn)/d(sx23) = -wrr*kt
  wn = wrr*kk
  d2 = un2 - wn
  d3 = un3 - wn

  gam = ga
  c = gam/(gam - 1.0_wp)

  g = 0.0_wp

! r1 = ro2*un2 - ro3*un3 - wn*(ro2-ro3)
  g(1, 5) = un2 - wn
  g(1, 9) = wn - un3
  g(1, 7) = ro2*nx
  g(1, 8) = ro2*ny
  g(1, 11) = -ro3*nx
  g(1, 12) = -ro3*ny
  g(1, 14) = -ro2*ut2 + ro3*ut3 + wrr*kt*(ro2 - ro3)
  g(1, 15) = -kk*(ro2 - ro3)

! r2 = p2 + ro2*(un2-wn)**2 - p3 - ro3*(un3-wn)**2
  g(2, 6) = 1.0_wp
  g(2, 10) = -1.0_wp
  g(2, 5) = d2**2
  g(2, 9) = -d3**2
  g(2, 7) = 2.0_wp*ro2*d2*nx
  g(2, 8) = 2.0_wp*ro2*d2*ny
  g(2, 11) = -2.0_wp*ro3*d3*nx
  g(2, 12) = -2.0_wp*ro3*d3*ny
  g(2, 14) = 2.0_wp*ro2*d2*(-ut2 + wrr*kt) - 2.0_wp*ro3*d3*(-ut3 + wrr*kt)
  g(2, 15) = -2.0_wp*kk*(ro2*d2 - ro3*d3)

! r3 = ut2 - ut3
  g(3, 7) = tx
  g(3, 8) = ty
  g(3, 11) = -tx
  g(3, 12) = -ty
  g(3, 14) = un2 - un3

! r4 = c*p2/ro2 + 0.5*(un2-wn)**2 - c*p3/ro3 - 0.5*(un3-wn)**2
  g(4, 5) = -c*p2/ro2**2
  g(4, 6) = c/ro2
  g(4, 9) = c*p3/ro3**2
  g(4, 10) = -c/ro3
  g(4, 7) = d2*nx
  g(4, 8) = d2*ny
  g(4, 11) = -d3*nx
  g(4, 12) = -d3*ny
  g(4, 14) = d2*(-ut2 + wrr*kt) - d3*(-ut3 + wrr*kt)
  g(4, 15) = -kk*(d2 - d3)

! r5 = u3*nwx + v3*nwy, nwx=-tauwy, nwy=tauwx (constants)
  g(5, 11) = -tauwy
  g(5, 12) = tauwx

! rows 6-15: linear "known-value" rows, Jacobian = ctx%a's row directly
  do i = 6, 15
    do j = 1, 15
      g(i, j) = a(i - 5, j)
    end do
  end do

end subroutine jfutp1
