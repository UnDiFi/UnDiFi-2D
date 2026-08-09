module test_co_uqp
! Phase 7.5 increment 5 (ROADMAP.md #20): co_uqp -- the quadruple-point
! (shock-shock crossing, opposite family) special point's Newton solve --
! against a manufactured two-shock solution.
!
! co_uqp's 24-unknown y vector: y(1:4) an unused passthrough slot (like
! co_urr's state1); y(5:8)/y(9:12) the two shocks' upstream states
! (block A / shock 24, block B / shock 35 -- fixed inputs, held via
! constraint rows exactly like co_urr's state2); y(13:16)/y(17:20) the
! two shocks' downstream states (the real per-block unknowns); y(21)/
! y(23) unused passthrough slots (like y(1:4)); y(22)/y(24) the two
! shocks' orientations sx (the other real unknowns). wqpx/wqpy (the
! quadruple point's own velocity) are plain subroutine arguments here,
! not part of y -- unlike co_urr's wrr, there's no y-component to pin.
!
! Each block (eqs 1-4 / 5-8) is algebraically the SAME 4-equation R-H
! block as co_urr (mass/momentum/tangential-velocity/energy across a
! shock, wn=wqpx*nx+wqpy*ny here rather than wrr*(...)), reused via the
! same M1-based normal-shock construction as co_shock/co_urr, again with
! wn=0 (wqpx=wqpy=0 here, i.e. a stationary quadruple point) to keep the
! construction closed-form. What's new is the two CROSS-block equations
! coupling the blocks' downstream states across the quadruple point's
! central slip line: eq9 (equal downstream pressure) and eq10 (parallel
! downstream velocity, via the (a.b)**2=|a|**2|b|**2 collinearity
! identity). Satisfied by construction rather than solved: both shocks
! share the same orientation sx, and each block's tangential velocity is
! set to ut=k*un2 for the SAME k -- so both downstream velocity vectors
! are un2*[(nx,ny)+k*(taux,tauy)], i.e. exactly parallel (different
! magnitude allowed, since un2 differs per block) -- and block B's base
! pressure pmB is solved in closed form (pmB=p2A/p_ratio(M1B)) so its
! downstream pressure comes out exactly equal to block A's. Verified
! numerically (all 10 residuals ~1e-15) before writing the Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_uqp

contains

  subroutine collect_co_uqp(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("quadruple_point_asymmetric", test_asymmetric)&
    &]
  end subroutine collect_co_uqp

! Computes one shock block's (upstream, downstream) states given a base
! (rom,pm), Mach number M1 and tangential-velocity proportionality k
! (ut=k*un2, applied identically upstream and downstream so tangential
! continuity, eq3/eq7, holds by construction), all in the shock-normal
! frame (wn=0).
  subroutine shock_block(gam, rom, pm, M1, k, un1, ut1, ro2, p2, un2, ut2)
    real(wp), intent(in) :: gam, rom, pm, M1, k
    real(wp), intent(out) :: un1, ut1, ro2, p2, un2, ut2
    real(wp) :: am, rho_ratio, p_ratio

    am = sqrt(gam*pm/rom)
    un1 = M1*am
    rho_ratio = (gam + 1.0_wp)*M1**2/((gam - 1.0_wp)*M1**2 + 2.0_wp)
    p_ratio = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(M1**2 - 1.0_wp)
    ro2 = rom*rho_ratio
    p2 = pm*p_ratio
    un2 = un1/rho_ratio
    ut2 = k*un2
    ut1 = ut2
  end subroutine shock_block

  subroutine test_asymmetric(error)
    type(error_type), allocatable, intent(out) :: error

    include 'paramt.h'
    external :: co_uqp

    real(wp), parameter :: gam = 1.4_wp, sx_true = 0.8_wp, kk = 0.3_wp
    real(wp), parameter :: romA = 1.0_wp, pmA = 1.0_wp, M1A = 2.5_wp
    real(wp), parameter :: romB = 1.2_wp, M1B = 1.8_wp

    real(wp) :: un1A, ut1A, ro2A, p2A, un2A, ut2A
    real(wp) :: pmB, un1B, ut1B, ro2B, p2B, un2B, ut2B
    real(wp) :: p_ratioB, nx, ny, tx, ty
    real(wp) :: u1A, v1A, u2A, v2A, u1B, v1B, u2B, v2B
    real(wp) :: y(24), yn1(24)

    ga = gam
    gm1 = gam - 1.0_wp

    call shock_block(gam, romA, pmA, M1A, kk, un1A, ut1A, ro2A, p2A, un2A, ut2A)

    ! solve pmB in closed form so block B's downstream pressure matches
    ! block A's exactly (eq9), then build block B the same way
    p_ratioB = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(M1B**2 - 1.0_wp)
    pmB = p2A/p_ratioB
    call shock_block(gam, romB, pmB, M1B, kk, un1B, ut1B, ro2B, p2B, un2B, ut2B)

    nx = -sin(sx_true); ny = cos(sx_true)
    tx = cos(sx_true); ty = sin(sx_true)

    u1A = un1A*nx + ut1A*tx; v1A = un1A*ny + ut1A*ty
    u2A = un2A*nx + ut2A*tx; v2A = un2A*ny + ut2A*ty
    u1B = un1B*nx + ut1B*tx; v1B = un1B*ny + ut1B*ty
    u2B = un2B*nx + ut2B*tx; v2B = un2B*ny + ut2B*ty

    y(1:4) = [romA, pmA, un1A, 0.0_wp]  ! unused passthrough
    y(5:8) = [romA, pmA, u1A, v1A]      ! block A upstream: exact fixed input
    y(9:12) = [romB, pmB, u1B, v1B]     ! block B upstream: exact fixed input
    ! perturbed initial guesses for the real unknowns
    y(13:16) = [1.05_wp*ro2A, 0.95_wp*p2A, 1.05_wp*u2A, 1.05_wp*v2A]
    y(17:20) = [0.95_wp*ro2B, 1.05_wp*p2B, 0.95_wp*u2B, 0.95_wp*v2B]
    y(21) = 0.0_wp                      ! unused passthrough
    y(22) = sx_true - 0.03_wp           ! perturbed shock-24-angle guess
    y(23) = 0.0_wp                      ! unused passthrough
    y(24) = sx_true + 0.02_wp           ! perturbed shock-35-angle guess

    call co_uqp(y, 0.0_wp, 0.0_wp, yn1)  ! stationary quadruple point

    call check(error, yn1(13), ro2A, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(14), p2A, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(15), u2A, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(16), v2A, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(17), ro2B, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(18), p2B, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(19), u2B, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(20), v2B, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(22), sx_true, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(24), sx_true, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine test_asymmetric

end module test_co_uqp
