! Compute an unsteady quadruple point (QP)

module co_uqp_ctx_m
! Phase 3.2 increment 1 (ROADMAP.md #14): closure payload for futp2's
! residual evaluation, replacing the (a, b, wqpx, wqpy) argument list
! co_uqp used to pass on every call into futp2, now that both go
! through mod_newton_solve's generic residual_if(i, y, ctx) interface.
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: uqp_ctx_t

  type :: uqp_ctx_t
    real(wp) :: a(14, 24) = 0.0_wp
    real(wp) :: b(14) = 0.0_wp
    real(wp) :: wqpx = 0.0_wp, wqpy = 0.0_wp
  end type uqp_ctx_t

end module co_uqp_ctx_m

subroutine co_uqp(y, wqpx, wqpy, yn1)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_newton_solve, only: newton_solve, residual_if, jacobian_if
  use co_uqp_ctx_m, only: uqp_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: futp2
  procedure(jacobian_if) :: jfutp2
  include 'paramt.h'

  integer(i4) i, j, nn
  real(wp) gam, delta
  real(wp) wqpx, wqpy
  real(wp) y, yn1
  dimension y(24), yn1(24)
  type(uqp_ctx_t) :: ctx

  open (8, file='log/co_uqp.log')

!     assign costans
  gam = ga
  delta = (gam - 1.0)/2.

  nn = 24

! build matrix a and vector b of the additional conditions
! which varies depending on the case

! calculation of incident shock 1, incident shock 2 being the upstream
! state and deviation known
  do i = 1, (3*4 + 2)
    do j = 1, nn
      ctx%a(i, j) = 0.d+0
    end do
    ctx%b(i) = 0.d+0
  end do

! state 1 known
  ctx%a(1, 1) = 1.0d+0
  ctx%a(2, 2) = 1.0d+0
  ctx%a(3, 3) = 1.0d+0
  ctx%a(4, 4) = 1.0d+0
  ctx%b(1) = y(1)
  ctx%b(2) = y(2)
  ctx%b(3) = y(3)
  ctx%b(4) = y(4)

! state 2 known
  ctx%a(5, 5) = 1.0d+0
  ctx%a(6, 6) = 1.0d+0
  ctx%a(7, 7) = 1.0d+0
  ctx%a(8, 8) = 1.0d+0
  ctx%b(5) = y(5)
  ctx%b(6) = y(6)
  ctx%b(7) = y(7)
  ctx%b(8) = y(8)

! state 3 known
  ctx%a(9, 9) = 1.0d+0
  ctx%a(10, 10) = 1.0d+0
  ctx%a(11, 11) = 1.0d+0
  ctx%a(12, 12) = 1.0d+0
  ctx%b(9) = y(9)
  ctx%b(10) = y(10)
  ctx%b(11) = y(11)
  ctx%b(12) = y(12)

! value sx12 known
  ctx%a(13, 21) = 1.0d+0
  ctx%b(13) = y(21)

! value sx13 known
  ctx%a(14, 23) = 1.0d+0
  ctx%b(14) = y(23)

  ctx%wqpx = wqpx
  ctx%wqpy = wqpy

! calcuate the downstream state and shock velocity
! with the newton-raphson method
  call newton_solve(nn, y, futp2, ctx, 0.5_wp, 1.0e-11_wp, 0.001_wp, yn1,&
  &log_unit=8_i4, jac=jfutp2)

  write (8, *)
  do i = 1, nn
    write (8, *) i, y(i), yn1(i)
  end do

  close (8)

end subroutine co_uqp

real(wp) function futp2(i, y, ctx) result(r)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use co_uqp_ctx_m, only: uqp_ctx_t
  implicit none(type, external)
  include 'paramt.h'

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  integer(i4) ii, j
  real(wp) a(14, 24), b(14), wqpx, wqpy, wn
  real(wp) ro1, ro2, p1, p2, u1, u2, gam, delta
  real(wp) v1, v2, un1, un2, ut1, ut2
  real(wp) taux12, tauy12, nx12, ny12, sx12

  select type (ctx)
  type is (uqp_ctx_t)
    a = ctx%a
    b = ctx%b
    wqpx = ctx%wqpx
    wqpy = ctx%wqpy
  end select

! assign constants and variables
  gam = ga
  delta = (gam - 1.0)/2.

  if (i .ge. 1 .and. i .le. 4) then
    ro1 = y(5)
    p1 = y(6)
    u1 = y(7)
    v1 = y(8)
    ro2 = y(13)
    p2 = y(14)
    u2 = y(15)
    v2 = y(16)
    sx12 = y(22)

! calculation of the tangent and normal versor to shock 24 (reflected shock 1)
    taux12 = cos(sx12)
    tauy12 = sin(sx12)
    nx12 = -tauy12
    ny12 = taux12

! calculation of the normal and tangential velocity components
    un1 = u1*nx12 + v1*ny12
    un2 = u2*nx12 + v2*ny12
    ut1 = u1*taux12 + v1*tauy12
    ut2 = u2*taux12 + v2*tauy12

! calculation of the normal velocity component to shock 24 in the quadruple point
    wn = wqpx*nx12 + wqpy*ny12

    r = 0.d0
    if (i .eq. 1) then
      r = ro1*un1 - ro2*un2 - wn*(ro1 - ro2)
    elseif (i .eq. 2) then
      r = p1 + ro1*(un1 - wn)**2 - p2 - ro2*(un2 - wn)**2
    elseif (i .eq. 3) then
      r = ut1 - ut2
    elseif (i .eq. 4) then
      r = gam/(gam - 1.0)*p1/ro1 + 0.5*(un1 - wn)**2&
      &- gam/(gam - 1.0)*p2/ro2 - 0.5*(un2 - wn)**2

    end if
  elseif (i .ge. 5 .and. i .le. 8) then
    ro1 = y(9)
    p1 = y(10)
    u1 = y(11)
    v1 = y(12)
    ro2 = y(17)
    p2 = y(18)
    u2 = y(19)
    v2 = y(20)
    sx12 = y(24)

! calculation of the tangent and normal versor to shock 35 (reflected shock 2)
    taux12 = cos(sx12)
    tauy12 = sin(sx12)
    nx12 = -tauy12
    ny12 = taux12

! calculation of the normal velocity component to shock 35 in the quadruple point
    wn = wqpx*nx12 + wqpy*ny12

! calculation of the normal and tangential velocity components
    un1 = u1*nx12 + v1*ny12
    un2 = u2*nx12 + v2*ny12
    ut1 = u1*taux12 + v1*tauy12
    ut2 = u2*taux12 + v2*tauy12

    r = 0.d0
    if (i .eq. 5) then
      r = ro1*un1 - ro2*un2 - wn*(ro1 - ro2)
    elseif (i .eq. 6) then
      r = p1 + ro1*(un1 - wn)**2 - p2 - ro2*(un2 - wn)**2

    elseif (i .eq. 7) then
      r = ut1 - ut2
    elseif (i .eq. 8) then
      r = gam/(gam - 1.0)*p1/ro1 + 0.5*(un1 - wn)**2&
      &- gam/(gam - 1.0)*p2/ro2 - 0.5*(un2 - wn)**2
    end if
  elseif (i .eq. 9) then
    ro1 = y(13)
    p1 = y(14)
    u1 = y(15)
    v1 = y(16)
    ro2 = y(17)
    p2 = y(18)
    u2 = y(19)
    v2 = y(20)
    r = p1 - p2
  elseif (i .eq. 10) then
    ro1 = y(13)
    p1 = y(14)
    u1 = y(15)
    v1 = y(16)
    ro2 = y(17)
    p2 = y(18)
    u2 = y(19)
    v2 = y(20)
    r = (u1*u2 + v1*v2)**2 - (u1*u1 + v1*v1)*(u2*u2 + v2*v2)
  elseif (i .ge. 11) then
    r = 0.d0
    ii = i - (2*4 + 1 + 1)
    do j = 1, 24
      r = r + a(ii, j)*y(j)
    end do
    r = r - b(ii)
  end if

end function futp2

! ************************************
! Phase 3.5 increment 4 (ROADMAP.md #14): analytic Jacobian of futp2.
! Equations 1-4 and 5-8 are two independent copies of the same R-H
! block used in co_urr/co_shock (mass/momentum/tangential-velocity/
! energy), each parameterized by its own shock angle (sx=y(22) for
! block A/shock 24, sx=y(24) for block B/shock 35); unlike co_urr's
! wn=wrr*(tauwx*nx+tauwy*ny), here wn=wqpx*nx+wqpy*ny directly (wqpx/
! wqpy are ctx constants, no separate unknown scaling it), so d(wn)/dsx
! = -(wqpx*tx+wqpy*ty) with no extra factor. Equation 9 is a plain
! pressure-equality row (linear). Equation 10 is a parallel-velocity
! condition between the two blocks' downstream states, differentiated
! directly (not via its (a.b)^2-|a|^2|b|^2 = -(a x b)^2 identity, to
! keep the derivation mechanical and easy to re-check against futp2's
! own source). Rows 11-24 are the linear known-value rows, Jacobian =
! ctx%a's row directly. Validated the same way as the other analytic
! Jacobians (see mod_newton_solve.f90) before being trusted.
subroutine jfutp2(y, ctx, g)
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use co_uqp_ctx_m, only: uqp_ctx_t
  implicit none(type, external)
  include 'paramt.h'

  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx
  real(wp), intent(out) :: g(:, :)

  real(wp) :: a(14, 24), wqpx, wqpy
  real(wp) :: ro1, p1, u1, v1, ro2, p2, u2, v2, sx
  real(wp) :: gam, c
  real(wp) :: nx, ny, tx, ty, un1, un2, ut1, ut2, wn, d1, d2, kk, kt
  real(wp) :: uA, vA, uB, vB, aa, bb, cc
  integer(i4) :: i, j

  select type (ctx)
  type is (uqp_ctx_t)
    a = ctx%a
    wqpx = ctx%wqpx
    wqpy = ctx%wqpy
  end select

  gam = ga
  c = gam/(gam - 1.0_wp)

  g = 0.0_wp

! ---- block A: eq 1-4, shock 24 ----
  ro1 = y(5); p1 = y(6); u1 = y(7); v1 = y(8)
  ro2 = y(13); p2 = y(14); u2 = y(15); v2 = y(16)
  sx = y(22)

  tx = cos(sx); ty = sin(sx)
  nx = -ty; ny = tx

  un1 = u1*nx + v1*ny
  un2 = u2*nx + v2*ny
  ut1 = u1*tx + v1*ty
  ut2 = u2*tx + v2*ty
  kk = wqpx*nx + wqpy*ny ! wn = kk
  kt = wqpx*tx + wqpy*ty ! d(wn)/d(sx) = -kt
  wn = kk
  d1 = un1 - wn
  d2 = un2 - wn

  g(1, 5) = un1 - wn
  g(1, 13) = wn - un2
  g(1, 7) = ro1*nx
  g(1, 8) = ro1*ny
  g(1, 15) = -ro2*nx
  g(1, 16) = -ro2*ny
  g(1, 22) = -ro1*ut1 + ro2*ut2 + kt*(ro1 - ro2)

  g(2, 6) = 1.0_wp
  g(2, 14) = -1.0_wp
  g(2, 5) = d1**2
  g(2, 13) = -d2**2
  g(2, 7) = 2.0_wp*ro1*d1*nx
  g(2, 8) = 2.0_wp*ro1*d1*ny
  g(2, 15) = -2.0_wp*ro2*d2*nx
  g(2, 16) = -2.0_wp*ro2*d2*ny
  g(2, 22) = 2.0_wp*ro1*d1*(-ut1 + kt) - 2.0_wp*ro2*d2*(-ut2 + kt)

  g(3, 7) = tx
  g(3, 8) = ty
  g(3, 15) = -tx
  g(3, 16) = -ty
  g(3, 22) = un1 - un2

  g(4, 5) = -c*p1/ro1**2
  g(4, 6) = c/ro1
  g(4, 13) = c*p2/ro2**2
  g(4, 14) = -c/ro2
  g(4, 7) = d1*nx
  g(4, 8) = d1*ny
  g(4, 15) = -d2*nx
  g(4, 16) = -d2*ny
  g(4, 22) = d1*(-ut1 + kt) - d2*(-ut2 + kt)

! ---- block B: eq 5-8, shock 35 ----
  ro1 = y(9); p1 = y(10); u1 = y(11); v1 = y(12)
  ro2 = y(17); p2 = y(18); u2 = y(19); v2 = y(20)
  sx = y(24)

  tx = cos(sx); ty = sin(sx)
  nx = -ty; ny = tx

  un1 = u1*nx + v1*ny
  un2 = u2*nx + v2*ny
  ut1 = u1*tx + v1*ty
  ut2 = u2*tx + v2*ty
  kk = wqpx*nx + wqpy*ny
  kt = wqpx*tx + wqpy*ty
  wn = kk
  d1 = un1 - wn
  d2 = un2 - wn

  g(5, 9) = un1 - wn
  g(5, 17) = wn - un2
  g(5, 11) = ro1*nx
  g(5, 12) = ro1*ny
  g(5, 19) = -ro2*nx
  g(5, 20) = -ro2*ny
  g(5, 24) = -ro1*ut1 + ro2*ut2 + kt*(ro1 - ro2)

  g(6, 10) = 1.0_wp
  g(6, 18) = -1.0_wp
  g(6, 9) = d1**2
  g(6, 17) = -d2**2
  g(6, 11) = 2.0_wp*ro1*d1*nx
  g(6, 12) = 2.0_wp*ro1*d1*ny
  g(6, 19) = -2.0_wp*ro2*d2*nx
  g(6, 20) = -2.0_wp*ro2*d2*ny
  g(6, 24) = 2.0_wp*ro1*d1*(-ut1 + kt) - 2.0_wp*ro2*d2*(-ut2 + kt)

  g(7, 11) = tx
  g(7, 12) = ty
  g(7, 19) = -tx
  g(7, 20) = -ty
  g(7, 24) = un1 - un2

  g(8, 9) = -c*p1/ro1**2
  g(8, 10) = c/ro1
  g(8, 17) = c*p2/ro2**2
  g(8, 18) = -c/ro2
  g(8, 11) = d1*nx
  g(8, 12) = d1*ny
  g(8, 19) = -d2*nx
  g(8, 20) = -d2*ny
  g(8, 24) = d1*(-ut1 + kt) - d2*(-ut2 + kt)

! r9 = p1 - p2 (block A's downstream p vs block B's downstream p)
  g(9, 14) = 1.0_wp
  g(9, 18) = -1.0_wp

! r10 = (u1*u2+v1*v2)**2 - (u1**2+v1**2)*(u2**2+v2**2)
  uA = y(15); vA = y(16); uB = y(19); vB = y(20)
  aa = uA*uB + vA*vB
  bb = uA**2 + vA**2
  cc = uB**2 + vB**2
  g(10, 15) = 2.0_wp*(aa*uB - uA*cc)
  g(10, 16) = 2.0_wp*(aa*vB - vA*cc)
  g(10, 19) = 2.0_wp*(aa*uA - bb*uB)
  g(10, 20) = 2.0_wp*(aa*vA - bb*vB)

! rows 11-24: linear known-value rows, Jacobian = ctx%a's row directly
  do i = 11, 24
    do j = 1, 24
      g(i, j) = a(i - 10, j)
    end do
  end do

end subroutine jfutp2
