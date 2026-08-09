! Compute a shock crossing a material interface (Issue #21/E1: the R-H
! residual generalized to different gamma on the two sides, needed for
! a shock that is simultaneously a gas-composition change -- e.g. a
! shock transmitting across a shock-bubble interface).
!
! co_shock (co_shock.f90) with one shared `gam` split into `gamm`
! (upstream, x2's gas) and `gamv` (downstream, x1's gas). The mass (r1)
! and momentum (r2) R-H equations are gamma-independent, unchanged; only
! the energy equation (r3) actually mixes the two sides' thermodynamics
! and needs separate c=gam/(gam-1) coefficients. The Riemann-invariant
! closure (r4) already only ever involved the downstream side's own gas
! (gamv), so it is unchanged in form, just explicit that it's gamv, not
! a shared gam.
!
! Deliberately a separate kernel, not a co_shock signature change --
! same reasoning as co_interface.f90 (issue #23): every real call site
! of co_shock (co_state_dps.f90 etc.) reads one shared GA from
! COMMON/PARAMT/, with no per-shock gas storage. Wiring a real typesh
! dispatch to this kernel needs mod_gas/gas_registry actually threaded
! through those call sites -- this file is the solver kernel only,
! verified standalone (tests/unit/test_co_shock_2gas.f90).
!
! Not every (M1, gamm, gamv) combination has a real solution here --
! confirmed empirically (continuation-method scipy solve, see this
! commit's own session notes) while deriving the manufactured test
! case: for a fixed upstream Mach, sufficiently large |gamv-gamm| loses
! convergence (the 3-equation mass/momentum/energy system's real root
! branch appears to fold/vanish), not a solver bug -- pick modest gamma
! differences and/or weaker shocks when using this kernel until that
! solvability boundary is characterized properly (not done here, out of
! scope for this increment).
module co_shock_2gas_ctx_m
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: shock_2gas_ctx_t

  type :: shock_2gas_ctx_t
    real(wp) :: rom = 0.0_wp, pm = 0.0_wp, um = 0.0_wp
    real(wp) :: gamm = 0.0_wp, gamv = 0.0_wp, delta_v = 0.0_wp, r2 = 0.0_wp
  end type shock_2gas_ctx_t

end module co_shock_2gas_ctx_m

subroutine co_shock_2gas(x1, x2, wshk, R14, gamm, gamv)

!     x1(1) and x2(1) downstream and upstream density
!     x1(2) and x2(2) downstream and upstream pressure
!     x1(3) and x2(3) downstream and upstream normal velocity
!     wshk shock velocity
!     gamm, gamv  specific heat ratio upstream (x2's gas), downstream (x1's gas)

  use mod_kinds, only: wp, i4
  use mod_newton_solve, only: newton_solve, residual_if, jacobian_if
  use co_shock_2gas_ctx_m, only: shock_2gas_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: f2gas
  procedure(jacobian_if) :: jf2gas

  real(wp) x1, x2, wshk, R14, gamm, gamv
  dimension x1(4), x2(4)
  real(wp) rov, rom, pv, pm, uv, um
  real(wp) y0(4), yn1(4)
  type(shock_2gas_ctx_t) :: ctx

  rov = x1(1)
  pv = x1(2)
  uv = x1(3)

  rom = x2(1)
  pm = x2(2)
  um = x2(3)

  ctx%r2 = R14
  ctx%rom = rom
  ctx%pm = pm
  ctx%um = um
  ctx%gamm = gamm
  ctx%gamv = gamv
  ctx%delta_v = (gamv - 1.0_wp)/2.0_wp

  y0(1) = rov
  y0(2) = pv
  y0(3) = uv
  y0(4) = -0.001_wp

  call newton_solve(4_i4, y0, f2gas, ctx, 0.2_wp, 1.0e-7_wp, 0.001_wp, yn1, jac=jf2gas)

  wshk = yn1(4)
  x1(1) = yn1(1)
  x1(2) = yn1(2)
  x1(3) = yn1(3)

  return
end subroutine co_shock_2gas

real(wp) function f2gas(i, y, ctx) result(r)
  use mod_kinds, only: wp, i4
  use co_shock_2gas_ctx_m, only: shock_2gas_ctx_t
  implicit none(type, external)

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  real(wp) rov, pv, uv, w, rom, pm, um, gamm, gamv, delta_v, R2, cm, cv

  select type (ctx)
  type is (shock_2gas_ctx_t)
    rom = ctx%rom
    pm = ctx%pm
    um = ctx%um
    gamm = ctx%gamm
    gamv = ctx%gamv
    delta_v = ctx%delta_v
    R2 = ctx%r2
  end select

  rov = y(1)
  pv = y(2)
  uv = y(3)
  w = y(4)

  cm = gamm/(gamm - 1.0_wp)
  cv = gamv/(gamv - 1.0_wp)

  r = 0.0_wp
  if (i .eq. 1) then
    r = rov*(uv - w) - rom*(um - w)
  elseif (i .eq. 2) then
    r = pv + rov*(uv - w)**2 - pm - rom*(um - w)**2
  elseif (i .eq. 3) then
    r = cv*pv/rov + 0.5_wp*(uv - w)**2 - cm*pm/rom - 0.5_wp*(um - w)**2
  elseif (i .eq. 4) then
    r = sqrt(gamv*pv/rov) + delta_v*uv - R2
  end if
  return
end function f2gas

! Analytic Jacobian of f2gas. Rows 1-2 (mass/momentum) and row 4
! (Riemann invariant) are structurally identical to co_shock.f90's own
! jf, just with gamv explicit where jf used a shared gam. Row 3 (energy)
! only picks up cv in its Jacobian entries -- cm's term (-cm*pm/rom)
! is a function of ctx-fixed rom/pm only, zero derivative w.r.t. the
! unknowns (rov,pv,uv,w).
subroutine jf2gas(y, ctx, g)
  use mod_kinds, only: wp, i4
  use co_shock_2gas_ctx_m, only: shock_2gas_ctx_t
  implicit none(type, external)

  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx
  real(wp), intent(out) :: g(:, :)

  real(wp) rov, pv, uv, w, rom, um, gamv, delta_v, cv, s

  select type (ctx)
  type is (shock_2gas_ctx_t)
    rom = ctx%rom
    um = ctx%um
    gamv = ctx%gamv
    delta_v = ctx%delta_v
  end select

  rov = y(1)
  pv = y(2)
  uv = y(3)
  w = y(4)

  cv = gamv/(gamv - 1.0_wp)
  s = sqrt(gamv*pv/rov)

! r1 = rov*(uv-w) - rom*(um-w)
  g(1, 1) = uv - w
  g(1, 2) = 0.0_wp
  g(1, 3) = rov
  g(1, 4) = -rov + rom

! r2 = pv + rov*(uv-w)**2 - pm - rom*(um-w)**2
  g(2, 1) = (uv - w)**2
  g(2, 2) = 1.0_wp
  g(2, 3) = 2.0_wp*rov*(uv - w)
  g(2, 4) = -2.0_wp*rov*(uv - w) + 2.0_wp*rom*(um - w)

! r3 = cv*pv/rov + 0.5*(uv-w)**2 - cm*pm/rom - 0.5*(um-w)**2
! (cm's term -cm*pm/rom is ctx-fixed, zero derivative)
  g(3, 1) = -cv*pv/rov**2
  g(3, 2) = cv/rov
  g(3, 3) = uv - w
  g(3, 4) = -(uv - w) + (um - w)

! r4 = sqrt(gamv*pv/rov) + delta_v*uv - R2
  g(4, 1) = -0.5_wp*s/rov
  g(4, 2) = 0.5_wp*s/pv
  g(4, 3) = delta_v
  g(4, 4) = 0.0_wp

end subroutine jf2gas
