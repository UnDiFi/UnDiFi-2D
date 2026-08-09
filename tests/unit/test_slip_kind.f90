module test_slip_kind
! Issue #23 (E3, ROADMAP.md Part III): the new 'L' (slip line)
! discontinuity kind added to co_state_dps.f90/co_norm.f90, dispatched
! through the exact same co_dc/dcp_dpndnc calls as 'D' (contact) --
! see co_state_dps.f90's own header comment on that dispatch for why
! this is a semantic/architectural distinction (letting future code,
! e.g. DMR's second slipstream, identify a dynamically-generated
! slipstream apart from a generic material contact) rather than new
! numerics: co_dc's own 7-equation system never constrains tangential
! velocity on either side, so a tangential jump was already
! representable under 'D'.
!
! Two things checked here:
! 1. Regression safety: 'L' produces BIT-IDENTICAL output to 'D' for
!    the same input, at co_state_dps's own level (not just "the same
!    function gets called" by inspection -- run both, diff the result).
! 2. The physics the issue text actually asks for ("zero pressure jump,
!    tangential velocity jump"): feeding co_state_dps a state with
!    matched pressure and a genuine tangential-velocity difference
!    across the discontinuity (built from the same theta-beta-M
!    triple-point construction as test_three_shock_theory.f90 --
!    state3/state4 there already satisfy p3=p4 by construction, with
!    zero normal and nonzero, DIFFERENT tangential velocity once the
!    slip line's own tangent is aligned with phi_common) converges to a
!    result that still shows a tangential jump, not one erased by the
!    solve -- 'L' is not silently behaving like a shock's tangential-
!    velocity-continuity rule (co_state_dps.f90's "impose equality of
!    tangential components" line, which is 'S'-only).
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_slip_kind

contains

  subroutine collect_slip_kind(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("slip_kind_matches_contact_kind", test_regression),&
    &new_unittest("slip_kind_preserves_tangential_jump", test_physics)&
    &]
  end subroutine collect_slip_kind

! Runs co_state_dps on a single one-point discontinuity with the given
! typesh code and (rho,p,u,v) states (global frame) on each side, plus
! a chosen unit normal (nx,ny). Returns the resulting downstream/
! upstream (rho,p,u,v).
  subroutine run_one_point(kind_code, nx, ny, rho_d, p_d, u_d, v_d,&
  &rho_u, p_u, u_u, v_u, out_rho_d, out_p_d, out_u_d, out_v_d,&
  &out_rho_u, out_p_u, out_u_u, out_v_u)
    character(len=1), intent(in) :: kind_code
    real(wp), intent(in) :: nx, ny, rho_d, p_d, u_d, v_d, rho_u, p_u, u_u, v_u
    real(wp), intent(out) :: out_rho_d, out_p_d, out_u_d, out_v_d
    real(wp), intent(out) :: out_rho_u, out_p_u, out_u_u, out_v_u

    include 'paramt.h'
    external :: co_state_dps

    real(wp), allocatable :: xysh(:, :, :), zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), allocatable :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), allocatable :: vshnor(:, :, :), wsh(:, :, :)
    integer(i4) :: nshockpoints(nshmax), nshockedges(nshmax)
    character(len=1) :: typesh(nshmax)
    real(wp) :: gam, hh, z1, z2, z3, z4
    real(wp) :: out_z1d, out_z2d, out_z3d, out_z4d
    real(wp) :: out_z1u, out_z2u, out_z3u, out_z4u

    allocate (xysh(ndim, npshmax, nshmax), zroeshu(ndof, npshmax, nshmax),&
    &zroeshd(ndof, npshmax, nshmax), zroeshuold(ndof, npshmax, nshmax),&
    &zroeshdold(ndof, npshmax, nshmax), vshnor(ndim, npshmax, nshmax),&
    &wsh(ndim, npshmax, nshmax))
    xysh = 0.0_wp; zroeshu = 0.0_wp; zroeshd = 0.0_wp
    zroeshuold = 0.0_wp; zroeshdold = 0.0_wp; vshnor = 0.0_wp; wsh = 0.0_wp

    gam = 1.4_wp
    ga = gam
    gm1 = gam - 1.0_wp

    vshnor(1, 1, 1) = nx; vshnor(2, 1, 1) = ny

    ! downstream (zroeshd), global (rho,p,u,v) -> Roe sqrt-variables
    z1 = sqrt(rho_d)
    hh = gam/(gam - 1.0_wp)*p_d/rho_d + 0.5_wp*(u_d**2 + v_d**2)
    z2 = z1*hh; z3 = z1*u_d; z4 = z1*v_d
    zroeshd(:, 1, 1) = [z1, z2, z3, z4]

    ! upstream (zroeshu)
    z1 = sqrt(rho_u)
    hh = gam/(gam - 1.0_wp)*p_u/rho_u + 0.5_wp*(u_u**2 + v_u**2)
    z2 = z1*hh; z3 = z1*u_u; z4 = z1*v_u
    zroeshu(:, 1, 1) = [z1, z2, z3, z4]

    nshockpoints = 0; nshockpoints(1) = 1
    nshockedges = 0; nshockedges(1) = 0
    typesh = ' '; typesh(1) = kind_code

    call co_state_dps(xysh, zroeshu, zroeshd, zroeshuold, zroeshdold,&
    &vshnor, wsh, 1_i4, nshockpoints, nshockedges, typesh, 1_i4)

    out_z1d = zroeshd(1, 1, 1); out_z2d = zroeshd(2, 1, 1)
    out_z3d = zroeshd(3, 1, 1); out_z4d = zroeshd(4, 1, 1)
    out_rho_d = out_z1d**2
    out_u_d = out_z3d/out_z1d; out_v_d = out_z4d/out_z1d
    out_p_d = (gam - 1.0_wp)/gam*(out_z1d*out_z2d - 0.5_wp*(out_z3d**2 + out_z4d**2))

    out_z1u = zroeshu(1, 1, 1); out_z2u = zroeshu(2, 1, 1)
    out_z3u = zroeshu(3, 1, 1); out_z4u = zroeshu(4, 1, 1)
    out_rho_u = out_z1u**2
    out_u_u = out_z3u/out_z1u; out_v_u = out_z4u/out_z1u
    out_p_u = (gam - 1.0_wp)/gam*(out_z1u*out_z2u - 0.5_wp*(out_z3u**2 + out_z4u**2))
  end subroutine run_one_point

  subroutine test_regression(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: rho_d1, p_d1, u_d1, v_d1, rho_u1, p_u1, u_u1, v_u1
    real(wp) :: rho_d2, p_d2, u_d2, v_d2, rho_u2, p_u2, u_u2, v_u2

    ! arbitrary non-degenerate two-sided state, normal along (0,1)
    call run_one_point('D', 0.0_wp, 1.0_wp,&
    &1.3_wp, 1.2_wp, 0.10_wp, 0.30_wp, 1.0_wp, 1.0_wp, 0.10_wp, -0.05_wp,&
    &rho_d1, p_d1, u_d1, v_d1, rho_u1, p_u1, u_u1, v_u1)

    call run_one_point('L', 0.0_wp, 1.0_wp,&
    &1.3_wp, 1.2_wp, 0.10_wp, 0.30_wp, 1.0_wp, 1.0_wp, 0.10_wp, -0.05_wp,&
    &rho_d2, p_d2, u_d2, v_d2, rho_u2, p_u2, u_u2, v_u2)

    call check(error, rho_d2, rho_d1, thr=1.0e-14_wp, message="rho_d differs between D and L")
    if (allocated(error)) return
    call check(error, p_d2, p_d1, thr=1.0e-14_wp, message="p_d differs between D and L")
    if (allocated(error)) return
    call check(error, u_d2, u_d1, thr=1.0e-14_wp, message="u_d differs between D and L")
    if (allocated(error)) return
    call check(error, v_d2, v_d1, thr=1.0e-14_wp, message="v_d differs between D and L")
    if (allocated(error)) return
    call check(error, rho_u2, rho_u1, thr=1.0e-14_wp, message="rho_u differs between D and L")
    if (allocated(error)) return
    call check(error, p_u2, p_u1, thr=1.0e-14_wp, message="p_u differs between D and L")
    if (allocated(error)) return
    call check(error, u_u2, u_u1, thr=1.0e-14_wp, message="u_u differs between D and L")
    if (allocated(error)) return
    call check(error, v_u2, v_u1, thr=1.0e-14_wp, message="v_u differs between D and L")
    if (allocated(error)) return
  end subroutine test_regression

! State built the same way as test_three_shock_theory.f90's shock legs:
! two states sharing phi_common=0 (both velocities purely along +x), so
! with the slip line's normal along +y, un=0 on both sides (trivially,
! since v=0) and ut=u (nonzero, DIFFERENT on each side) -- exactly
! "zero pressure jump [by construction], tangential velocity jump".
!
! Unlike co_shock's R14 or co_utp/co_uqp's shared angle, co_dc derives
! its own Riemann invariants (R1/S1/R2/S2) directly from whatever state
! is passed in as the Newton seed -- there is no separate externally-
! supplied target to perturb away from and converge back to (any input
! satisfying p1=p2/un1=un2 is already an exact fixed point of co_dc's
! own equations by construction). So this test isn't a "recovers from a
! perturbed guess" check like the others -- it verifies co_dc/
! co_state_dps correctly PRESERVES a valid slip-line configuration
! through the full solve-and-rotate-back-to-global-frame pipeline,
! without corrupting or force-equalizing the tangential component
! (which co_dc never even looks at) along the way.
  subroutine test_physics(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp), parameter :: rho3 = 2.0847227117573404_wp, p_common = 3.77125746308266_wp
    real(wp), parameter :: u3 = 2.7775971951838696_wp
    real(wp), parameter :: rho4 = 2.4180659314079604_wp
    real(wp), parameter :: u4 = 2.946638574487431_wp

    real(wp) :: out_rho_d, out_p_d, out_u_d, out_v_d
    real(wp) :: out_rho_u, out_p_u, out_u_u, out_v_u

    call run_one_point('L', 0.0_wp, 1.0_wp,&
    &rho3, p_common, u3, 0.0_wp, rho4, p_common, u4, 0.0_wp,&
    &out_rho_d, out_p_d, out_u_d, out_v_d, out_rho_u, out_p_u, out_u_u, out_v_u)

    ! pressure stays matched (co_dc's r5, p1=p2)
    call check(error, out_p_d, out_p_u, thr=1.0e-6_wp,&
    &message="slip line did not preserve pressure matching")
    if (allocated(error)) return
    ! tangential velocity jump preserved, not erased (u3 /= u4 driven
    ! through unchanged -- co_dc never touches x1(4)/x2(4), and 'L'
    ! does not get the 'S'-only tangential-equality line)
    call check(error, out_u_d, u3, thr=1.0e-10_wp,&
    &message="downstream tangential velocity was altered")
    if (allocated(error)) return
    call check(error, out_u_u, u4, thr=1.0e-10_wp,&
    &message="upstream tangential velocity was altered")
    if (allocated(error)) return
    call check(error, abs(out_u_d - out_u_u) > 0.1_wp,&
    &message="tangential jump was erased -- 'L' behaved like a no-slip contact")
    if (allocated(error)) return
  end subroutine test_physics

end module test_slip_kind
