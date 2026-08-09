module test_co_utp
! Phase 7.5 increment 6 (ROADMAP.md #20): co_utp -- the triple-point
! special point's Newton solve (incident shock 12 / reflected shock 23 /
! Mach stem 14, joined by a slip line between states 3 and 4) -- against
! a manufactured two-shock-plus-invariant solution.
!
! co_utp's 20-unknown y vector: y(1:4)/y(5:8) states 1/2 (upstream of the
! Mach stem / reflected shock respectively -- fixed inputs, held via
! constraint rows exactly like co_urr/co_uqp's upstream states); y(9:12)/
! y(13:16) states 3/4 (downstream of shock 23 / shock 14 -- the real
! per-block unknowns); y(17)=sx12 (incident-shock angle -- despite
! appearing in both blocks' wn formula, it is ALSO pinned via a
! constraint row, i.e. an input like state1/state2, not solved here);
! y(18)/y(19) the two solved shocks' own orientations sx23/sx14; y(20)=
! wsh (triple-point velocity along shock12's tangent -- the third real
! shared unknown, appearing in both blocks' wn=wsh*(taux12*n+tauy12*n)).
!
! Each block (eqs 1-4 shock23, 5-8 shock14) is the same R-H block as
! co_urr/co_uqp; wn=0 is achieved here by manufacturing wsh_true=0
! (matching every other kernel in this suite's near-zero-seed pattern),
! which makes wn zero for BOTH blocks regardless of sx12/sx23/sx14 --
! avoiding having to reason about wsh's genuinely shared, nonzero-case
! coupling. What's new relative to co_uqp: eq9, a co_shock-style Riemann-
! invariant closure on state4 alone (using a FIXED direction (dxr14,
! dyr14), not the shock's own dynamic normal -- confirmed dead code in
! futp's own nx14/ny14 right before eq9, see jfutp's header comment),
! and eq11 uses a direct (unsquared) cross product for the state3/state4
! parallel-velocity condition rather than co_uqp's squared collinearity
! form. Constructed with the SAME simplifications as co_uqp (shared shock
! angle sx23=sx14=sx_true, tangential velocity ut=k*un2 with a shared k,
! block2's base pressure solved in closed form so p4 matches p3 exactly)
! plus r14 computed directly from the manufactured state4. Verified
! numerically (all 11 residuals ~1e-15) before writing the Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_utp

contains

  subroutine collect_co_utp(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("triple_point_asymmetric", test_asymmetric)&
    &]
  end subroutine collect_co_utp

! Same helper as test_co_uqp.f90's shock_block: one shock's (upstream,
! downstream) states given a base (rom,pm), Mach number M1 and tangential
! proportionality k (ut=k*un2, applied both sides so tangential
! continuity holds by construction), in the shock-normal frame (wn=0).
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
    external :: co_utp

    real(wp), parameter :: gam = 1.4_wp, sx_true = 0.8_wp, kk = 0.3_wp
    real(wp), parameter :: dxr14 = 1.0_wp, dyr14 = 0.0_wp
    real(wp), parameter :: rom2 = 1.0_wp, pm2 = 1.0_wp, M1_23 = 2.5_wp
    real(wp), parameter :: rom1 = 1.3_wp, M1_14 = 1.9_wp

    real(wp) :: un2_23, ut2_23, ro3, p3, un3, ut3
    real(wp) :: pm1, un1_14, ut1_14, ro4, p4, un4, ut4
    real(wp) :: p_ratio14, nx, ny, tx, ty
    real(wp) :: u1, v1, u2, v2, u3, v3, u4, v4
    real(wp) :: r14, un4_proj
    real(wp) :: y(20), yn1(20)
    logical :: ifail

    ga = gam
    gm1 = gam - 1.0_wp

    ! block 1: shock23, state2 (fixed input) -> state3 (unknown)
    call shock_block(gam, rom2, pm2, M1_23, kk, un2_23, ut2_23, ro3, p3, un3, ut3)

    ! block 2: shock14, state1 (fixed input) -> state4 (unknown); pm1 solved
    ! in closed form so p4 comes out exactly equal to p3 (eq10)
    p_ratio14 = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(M1_14**2 - 1.0_wp)
    pm1 = p3/p_ratio14
    call shock_block(gam, rom1, pm1, M1_14, kk, un1_14, ut1_14, ro4, p4, un4, ut4)

    nx = -sin(sx_true); ny = cos(sx_true)
    tx = cos(sx_true); ty = sin(sx_true)

    u2 = un2_23*nx + ut2_23*tx; v2 = un2_23*ny + ut2_23*ty
    u3 = un3*nx + ut3*tx; v3 = un3*ny + ut3*ty
    u1 = un1_14*nx + ut1_14*tx; v1 = un1_14*ny + ut1_14*ty
    u4 = un4*nx + ut4*tx; v4 = un4*ny + ut4*ty

    un4_proj = u4*dxr14 + v4*dyr14
    r14 = sqrt(gam*p4/ro4) + (gam - 1.0_wp)/2.0_wp*un4_proj

    y(1:4) = [rom1, pm1, u1, v1]        ! state1: exact fixed input
    y(5:8) = [rom2, pm2, u2, v2]        ! state2: exact fixed input
    ! perturbed initial guesses for the real unknowns
    y(9:12) = [1.05_wp*ro3, 0.95_wp*p3, 1.05_wp*u3, 1.05_wp*v3]
    y(13:16) = [0.95_wp*ro4, 1.05_wp*p4, 0.95_wp*u4, 0.95_wp*v4]
    y(17) = 0.4_wp                      ! sx12: exact fixed input (dead-code
                                         ! angle here since wsh_true=0 makes
                                         ! wn independent of it -- any value
                                         ! exercises the pinning row honestly)
    y(18) = sx_true - 0.03_wp           ! perturbed shock-23-angle guess
    y(19) = sx_true + 0.02_wp           ! perturbed shock-14-angle guess
    y(20) = 0.01_wp                     ! perturbed wsh guess (true value 0)

    call co_utp(y, r14, dxr14, dyr14, 0.0_wp, 0.0_wp, yn1, .false., ifail)

    call check(error, .not. ifail, message="co_utp reported ifail")
    if (allocated(error)) return
    call check(error, yn1(9), ro3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(10), p3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(11), u3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(12), v3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(13), ro4, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(14), p4, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(15), u4, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(16), v4, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(18), sx_true, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(19), sx_true, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(20), 0.0_wp, thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine test_asymmetric

end module test_co_utp
