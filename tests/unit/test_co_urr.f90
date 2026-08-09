module test_co_urr
! Phase 7.5 increment 4 (ROADMAP.md #20): co_urr -- the regular-reflection
! special point's Newton solve (sp%solve_state for regular_reflection_t,
! dispatching through fx_state_dps.f90 in production) -- against a
! manufactured reflected-shock solution.
!
! co_urr's 15-unknown y vector is: y(1:4) an unrelated, physically-unused
! state1 slot (rows 1-4 just pin it to whatever's passed in -- dead
! weight in this Newton system, not exercised here beyond needing SOME
! value); y(5:8) state2 (behind the incident shock, held fixed at
! whatever's passed in via rows 5-8 -- a genuine INPUT, not a seed);
! y(9:12) state3 (behind the reflected shock, the real unknown);
! y(13)=sx12 (unused in the physics, pinned like state1); y(14)=sx23
! (reflected-shock orientation, the other real unknown); y(15)=wrr
! (reflection-point velocity along the wall, pinned like state2 -- an
! input, here fixed at 0 to match [[co_shock's own near-zero-seed
! regime]] -- a stationary reflection point in this frame).
!
! Manufactured solution: state2<->state3 must satisfy the same 3
! algebraic mass/momentum/energy equations as co_shock's R-H system (in
! the shock-normal frame, with wn=0 exactly since wrr=0), which is
! algebraically symmetric under swapping which side is called "upstream"
! -- so this test reuses co_shock's own M1-based normal-shock construction
! verbatim, just relabeled (state2 plays co_shock's "upstream/m" role,
! state3 its "downstream/v" role). The 4th equation (tangential velocity
! continuity, ut2=ut3) is trivially satisfied by construction: state2 and
! state3 are built by rotating the SAME tangential component ut back into
! (x,y) with un2/un3 from the normal-shock relations. The 5th equation
! (state3's velocity exactly tangent to the wall) is solved in closed
! form for ut, not iteratively: with the wall fixed at (1,0), tangency
! needs v3=un3*ny23+ut*tauy23=0, linear in the one remaining free
! parameter ut, so ut = -un3*ny23/tauy23 makes it exact rather than
! approximate. Verified numerically (residuals ~1e-15) before writing
! the Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_urr

contains

  subroutine collect_co_urr(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("regular_reflection_moderate", test_moderate),&
    &new_unittest("regular_reflection_strong", test_strong)&
    &]
  end subroutine collect_co_urr

  subroutine check_regular_reflection(error, gam, rom, pm, M1, sx23_true)
    type(error_type), allocatable, intent(out) :: error
    real(wp), intent(in) :: gam, rom, pm, M1, sx23_true

    include 'paramt.h'
    external :: co_urr

    real(wp) :: am, un2, rho_ratio, p_ratio, ro3, p3, un3
    real(wp) :: nx23, ny23, taux23, tauy23, ut
    real(wp) :: u2, v2, u3, v3
    real(wp) :: y(15), yn1(15)

    ga = gam
    gm1 = gam - 1.0_wp

    am = sqrt(gam*pm/rom)
    un2 = M1*am

    rho_ratio = (gam + 1.0_wp)*M1**2/((gam - 1.0_wp)*M1**2 + 2.0_wp)
    p_ratio = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(M1**2 - 1.0_wp)
    ro3 = rom*rho_ratio
    p3 = pm*p_ratio
    un3 = un2/rho_ratio

    nx23 = -sin(sx23_true); ny23 = cos(sx23_true)
    taux23 = cos(sx23_true); tauy23 = sin(sx23_true)

    ut = -un3*ny23/tauy23   ! makes v3 (below) exactly 0: wall tangency

    u3 = un3*nx23 + ut*taux23
    v3 = un3*ny23 + ut*tauy23
    u2 = un2*nx23 + ut*taux23
    v2 = un2*ny23 + ut*tauy23

    y(1:4) = [rom, pm, un2, 0.0_wp]  ! state1: physically unused, any value
    y(5:8) = [rom, pm, u2, v2]       ! state2: exact -- this is a fixed input
    ! perturbed initial guess for the real unknowns (state3, sx23)
    y(9:12) = [1.05_wp*ro3, 0.95_wp*p3, 1.05_wp*u3, u3*0.05_wp - v3*0.05_wp]
    y(13) = 0.0_wp                   ! sx12: physically unused
    y(14) = sx23_true + 0.03_wp      ! perturbed shock-angle guess
    y(15) = 0.0_wp                   ! wrr: exact fixed input (stationary)

    call co_urr(y, 1.0_wp, 0.0_wp, yn1)  ! wall tangent (tauwx,tauwy)=(1,0)

    call check(error, yn1(9), ro3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(10), p3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(11), u3, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(12), v3, thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, yn1(14), sx23_true, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine check_regular_reflection

  subroutine test_moderate(error)
    type(error_type), allocatable, intent(out) :: error
    call check_regular_reflection(error, 1.4_wp, 1.0_wp, 1.0_wp, 2.5_wp, 1.0_wp)
  end subroutine test_moderate

  subroutine test_strong(error)
    type(error_type), allocatable, intent(out) :: error
    call check_regular_reflection(error, 1.4_wp, 1.0_wp, 1.0_wp, 8.0_wp, 1.3_wp)
  end subroutine test_strong

end module test_co_urr
