! Compute a material interface point (Issue #23/E3: 'interface' kind --
! a contact-like discontinuity where the two sides carry DIFFERENT
! specific-heat ratios, e.g. two different gases meeting across a
! shock-bubble interface).
!
! This is co_dc (co_shock.f90) generalized from one shared `gam` to two
! independent `gam1`/`gam2`: same 7-equation structure (Riemann
! invariant + isentrope on each side, pressure/normal-velocity/
! interface-speed matching across the two), just with each side's own
! gamma governing its own two equations. Deliberately a SEPARATE kernel
! rather than adding gam1/gam2 arguments to co_dc itself -- the
! ubiquitous 'D' call sites (co_state_dps.f90 etc.) all read a single
! shared GA from COMMON/PARAMT/ with no per-shock gamma storage at all;
! wiring per-shock gamma through the real dispatch chain is E1/#21
! (multi-species) territory ("mod_gas: type gas_t with gamma, R...",
! not built yet), not this kernel's own concern. This file is the
! solver kernel only, verified standalone (tests/unit/test_interface.f90)
! -- NOT dispatched from co_state_dps.f90/co_norm.f90's typesh checks,
! unlike 'L' (slip). Wiring it in is exactly the "Assess per solver"/
! gas-registry groundwork #21 already calls out as its own prerequisite.
module co_interface_ctx_m
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: interface_ctx_t

  type :: interface_ctx_t
    real(wp) :: gam1 = 0.0_wp, gam2 = 0.0_wp
    real(wp) :: delta1 = 0.0_wp, delta2 = 0.0_wp
    real(wp) :: r1 = 0.0_wp, r2 = 0.0_wp, s1 = 0.0_wp, s2 = 0.0_wp
  end type interface_ctx_t

end module co_interface_ctx_m

subroutine co_interface(x1, x2, wintf, gam1, gam2)

!     x1(1) and x2(1) upstream and downstream density
!     x1(2) and x2(2) upstream and downstream pressure
!     x1(3) and x2(3) upstream and downstream normal velocity
!     gam1, gam2  specific heat ratio on side 1 (x1) and side 2 (x2)
!     wintf  interface velocity

  use mod_kinds, only: wp, i4
  use mod_newton_solve, only: newton_solve, residual_if, jacobian_if
  use co_interface_ctx_m, only: interface_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: fintf
  procedure(jacobian_if) :: jfintf

  real(wp) x1, x2, wintf, gam1, gam2
  dimension x1(4), x2(4)

  real(wp) ro1, ro2, p1, p2, u1, u2
  real(wp) y0(7), yn1(7)
  type(interface_ctx_t) :: ctx

  ro1 = x1(1)
  p1 = x1(2)
  u1 = x1(3)
  ro2 = x2(1)
  p2 = x2(2)
  u2 = x2(3)

! compute invariants, one gamma per side
  ctx%gam1 = gam1
  ctx%gam2 = gam2
  ctx%delta1 = (gam1 - 1.0_wp)/2.0_wp
  ctx%delta2 = (gam2 - 1.0_wp)/2.0_wp
  ctx%r1 = sqrt(gam1*p1/ro1) + ctx%delta1*u1
  ctx%r2 = sqrt(gam2*p2/ro2) - ctx%delta2*u2
  ctx%s1 = p1/ro1**gam1
  ctx%s2 = p2/ro2**gam2

! initialize the vector of unknowns
  y0(1) = ro1
  y0(2) = p1
  y0(3) = u1
  y0(4) = ro2
  y0(5) = p2
  y0(6) = u2
  y0(7) = wintf

  call newton_solve(7_i4, y0, fintf, ctx, 1.0_wp, 1.0e-10_wp, 0.01_wp, yn1, jac=jfintf)

  wintf = yn1(7)
  x1(1) = yn1(1)
  x1(2) = yn1(2)
  x1(3) = yn1(3)
  x2(1) = yn1(4)
  x2(2) = yn1(5)
  x2(3) = yn1(6)

  return
end subroutine co_interface

real(wp) function fintf(i, y, ctx) result(r)
  use mod_kinds, only: wp, i4
  use co_interface_ctx_m, only: interface_ctx_t
  implicit none(type, external)

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  real(wp) ro1, ro2, p1, p2, u1, u2, w
  real(wp) gam1, gam2, delta1, delta2, R1, R2, S1, S2

  select type (ctx)
  type is (interface_ctx_t)
    gam1 = ctx%gam1; gam2 = ctx%gam2
    delta1 = ctx%delta1; delta2 = ctx%delta2
    R1 = ctx%r1; R2 = ctx%r2
    S1 = ctx%s1; S2 = ctx%s2
  end select

  ro1 = y(1)
  p1 = y(2)
  u1 = y(3)
  ro2 = y(4)
  p2 = y(5)
  u2 = y(6)
  w = y(7)

  r = 0.0_wp
  if (i .eq. 1) then
    r = sqrt(gam1*p1/ro1) + delta1*u1 - R1
  elseif (i .eq. 2) then
    r = p1/ro1**gam1 - S1
  elseif (i .eq. 3) then
    r = sqrt(gam2*p2/ro2) - delta2*u2 - R2
  elseif (i .eq. 4) then
    r = p2/ro2**gam2 - S2
  elseif (i .eq. 5) then
    r = p1 - p2
  elseif (i .eq. 6) then
    r = u1 - u2
  elseif (i .eq. 7) then
    r = w - u1
  end if

  return
end function fintf

! Analytic Jacobian of fintf -- identical structure to co_dc's jfdc
! (co_shock.f90), with gam1 governing rows 1-2/columns 1-3 (side 1) and
! gam2 governing rows 3-4/columns 4-6 (side 2) instead of one shared gam.
subroutine jfintf(y, ctx, g)
  use mod_kinds, only: wp, i4
  use co_interface_ctx_m, only: interface_ctx_t
  implicit none(type, external)

  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx
  real(wp), intent(out) :: g(:, :)

  real(wp) ro1, p1, ro2, p2, gam1, gam2, delta1, delta2, s1, s3, t2, t4

  select type (ctx)
  type is (interface_ctx_t)
    gam1 = ctx%gam1; gam2 = ctx%gam2
    delta1 = ctx%delta1; delta2 = ctx%delta2
  end select

  ro1 = y(1)
  p1 = y(2)
  ro2 = y(4)
  p2 = y(5)

  s1 = sqrt(gam1*p1/ro1)
  s3 = sqrt(gam2*p2/ro2)
  t2 = p1/ro1**gam1
  t4 = p2/ro2**gam2

  g = 0.0_wp

! r1 = sqrt(gam1*p1/ro1) + delta1*u1 - R1
  g(1, 1) = -0.5_wp*s1/ro1
  g(1, 2) = 0.5_wp*s1/p1
  g(1, 3) = delta1

! r2 = p1/ro1**gam1 - S1
  g(2, 1) = -gam1*t2/ro1
  g(2, 2) = t2/p1

! r3 = sqrt(gam2*p2/ro2) - delta2*u2 - R2
  g(3, 4) = -0.5_wp*s3/ro2
  g(3, 5) = 0.5_wp*s3/p2
  g(3, 6) = -delta2

! r4 = p2/ro2**gam2 - S2
  g(4, 4) = -gam2*t4/ro2
  g(4, 5) = t4/p2

! r5 = p1 - p2
  g(5, 2) = 1.0_wp
  g(5, 5) = -1.0_wp

! r6 = u1 - u2
  g(6, 3) = 1.0_wp
  g(6, 6) = -1.0_wp

! r7 = w - u1
  g(7, 3) = -1.0_wp
  g(7, 7) = 1.0_wp

end subroutine jfintf
