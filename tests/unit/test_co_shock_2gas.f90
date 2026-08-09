module test_co_shock_2gas
! Issue #21 (E1): co_shock_2gas -- co_shock's R-H residual generalized
! to different gamma on the two sides, against a manufactured two-gas
! normal-shock solution.
!
! Unlike co_shock's own test (test_co_shock.f90, issue #20), the
! manufactured solution here can't come from the classical single-gamma
! rho_ratio/p_ratio closed-form formulas -- those were derived assuming
! one shared gamma in the energy equation. With gamm (upstream) /=
! gamv (downstream), the 3-equation mass/momentum/energy system
! (mass and momentum are still gamma-independent; only energy mixes the
! two) has no simple closed form and was solved numerically (scipy
! fsolve, continuation in gamv from the known single-gamma answer at
! gamv=gamm up to the target, to stay on the correct solution branch --
! a cold-start solve from a generic guess repeatedly diverged or lost
! convergence, see below) before writing the Fortran assertions.
!
! Real finding worth keeping in mind for anyone extending this kernel:
! not every (M1, gamm, gamv) combination has a nearby real solution.
! Continuation from gamv=1.4 (M1~2.1 case) toward gamv=1.667 lost
! convergence partway (residual growing past ~1e-3 around gamv~1.65),
! consistent with the solution branch folding or the Jacobian
! approaching singularity rather than a numerical-method failure --
! this test deliberately uses a more modest gamma gap (1.4 -> 1.5) and
! a weaker shock (M1~1.7) where continuation stayed converged to full
! double precision at every step.
  use mod_kinds, only: wp, i4
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_shock_2gas

contains

  subroutine collect_co_shock_2gas(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("normal_shock_two_gases", test_2gas)&
    &]
  end subroutine collect_co_shock_2gas

  subroutine test_2gas(error)
    type(error_type), allocatable, intent(out) :: error
    external :: co_shock_2gas

    ! manufactured solution (scipy continuation solve, see module header)
    real(wp), parameter :: gamm = 1.4_wp, gamv = 1.5_wp
    real(wp), parameter :: rom = 1.0_wp, pm = 1.0_wp, um = 2.0_wp
    real(wp), parameter :: rov_exact = 1.5669152706817999_wp
    real(wp), parameter :: pv_exact = 2.4472135954999596_wp
    real(wp), parameter :: uv_exact = 1.2763932022500204_wp
    real(wp), parameter :: delta_v = (gamv - 1.0_wp)/2.0_wp

    real(wp) :: R14, x1(4), x2(4), wshk

    R14 = sqrt(gamv*pv_exact/rov_exact) + delta_v*uv_exact

    x2 = [rom, pm, um, 0.0_wp]
    ! perturbed initial guess for the real unknowns
    x1 = [1.05_wp*rov_exact, 0.95_wp*pv_exact, 1.05_wp*uv_exact, 0.0_wp]

    call co_shock_2gas(x1, x2, wshk, R14, gamm, gamv)

    call check(error, x1(1), rov_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, x1(2), pv_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, x1(3), uv_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine test_2gas

end module test_co_shock_2gas
