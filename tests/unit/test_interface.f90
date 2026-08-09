module test_interface
! Issue #23 (E3): co_interface -- the material-interface jump-relation
! kernel (co_dc generalized to two independent specific-heat ratios,
! gam1/gam2, instead of one shared gam).
!
! Scope note (see co_interface.f90's own header comment): this tests the
! solver KERNEL only. It is NOT dispatched from co_state_dps.f90/
! co_norm.f90's typesh checks the way 'L' (slip, this session's earlier
! increment) is -- doing that needs per-shock gamma storage that
! doesn't exist yet (E1/#21's own "mod_gas"/gas-registry groundwork,
! not built). Wiring a real typesh code (e.g. 'I') to this kernel is
! deferred to that prerequisite landing.
!
! Like co_dc (see test_slip_kind.f90's header comment on the same
! issue), co_interface derives its own Riemann invariants directly from
! whatever state is passed in as the Newton seed -- there is no
! externally-supplied target to perturb away from and converge back to.
! So this suite checks two different things instead: (1) the residual
! EQUATIONS themselves, called directly (bypassing Newton) at a hand-
! constructed point that satisfies all 7 by construction with gam1 /=
! gam2 -- catches a wrong-gamma-on-the-wrong-side transcription bug,
! which a self-referential-invariant kernel's own Newton convergence
! can't; and (2) that swapping gam1/gam2 between two otherwise-identical
! calls changes the result -- proving both gammas are genuinely
! independent inputs, not one silently ignored or the two accidentally
! aliased to the same value.
  use mod_kinds, only: wp, i4
  use mod_newton_solve, only: residual_if
  use co_interface_ctx_m, only: interface_ctx_t
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_interface

contains

  subroutine collect_interface(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("residual_equations_use_correct_side_gamma", test_residuals),&
    &new_unittest("solve_preserves_a_valid_interface_state", test_solve),&
    &new_unittest("gammas_are_independent_not_aliased", test_gamma_independence)&
    &]
  end subroutine collect_interface

  subroutine test_residuals(error)
    type(error_type), allocatable, intent(out) :: error
    procedure(residual_if) :: fintf

    real(wp), parameter :: gam1 = 1.4_wp, gam2 = 1.667_wp
    real(wp), parameter :: ro1 = 1.0_wp, p1 = 1.0_wp, u1 = 0.3_wp
    real(wp), parameter :: ro2 = 1.2_wp, p2 = 1.0_wp, u2 = 0.3_wp
    type(interface_ctx_t) :: ctx
    real(wp) :: y(7)
    integer(i4) :: i

    ctx%gam1 = gam1; ctx%gam2 = gam2
    ctx%delta1 = (gam1 - 1.0_wp)/2.0_wp; ctx%delta2 = (gam2 - 1.0_wp)/2.0_wp
    ctx%r1 = sqrt(gam1*p1/ro1) + ctx%delta1*u1
    ctx%s1 = p1/ro1**gam1
    ctx%r2 = sqrt(gam2*p2/ro2) - ctx%delta2*u2
    ctx%s2 = p2/ro2**gam2

    y = [ro1, p1, u1, ro2, p2, u2, u1]

    do i = 1, 7
      call check(error, fintf(i, y, ctx), 0.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_residuals

  subroutine test_solve(error)
    type(error_type), allocatable, intent(out) :: error
    external :: co_interface

    real(wp), parameter :: gam1 = 1.4_wp, gam2 = 1.667_wp
    real(wp) :: x1(4), x2(4), wintf

    x1 = [1.0_wp, 1.0_wp, 0.3_wp, 0.0_wp]
    x2 = [1.2_wp, 1.0_wp, 0.3_wp, 0.0_wp]
    wintf = 0.3_wp

    call co_interface(x1, x2, wintf, gam1, gam2)

    call check(error, x1(1), 1.0_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, x1(2), 1.0_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, x1(3), 0.3_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, x2(1), 1.2_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, x2(2), 1.0_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, x2(3), 0.3_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, wintf, 0.3_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
  end subroutine test_solve

! A state that is NOT already cross-matched (u1 /= u2, unlike
! test_solve's already-self-consistent input) is solved twice, once as
! (gam1=1.4, gam2=1.667) and once with the two swapped. With a
! self-referential-invariant kernel like this one, an already-matched
! input is a trivial fixed point regardless of which gamma governs
! which side (see test_solve) -- that would NOT distinguish the two
! calls at all. Starting from a genuine mismatch forces Newton to move
! each side along its OWN Riemann-invariant curve (a different curve
! depending on that side's gamma) until they meet, so the converged
! common state actually depends on which gamma was assigned to which
! side. Verified independently via a from-scratch Python Newton solve
! before writing this assertion (~3% difference in the converged
! density/velocity between the two gamma assignments) -- this is a
! differential check (the two results must differ), not a comparison
! against a precomputed target.
  subroutine test_gamma_independence(error)
    type(error_type), allocatable, intent(out) :: error
    external :: co_interface

    real(wp) :: x1a(4), x2a(4), wa
    real(wp) :: x1b(4), x2b(4), wb

    x1a = [1.0_wp, 1.0_wp, 0.5_wp, 0.0_wp]
    x2a = [1.2_wp, 1.0_wp, 0.1_wp, 0.0_wp]
    wa = 0.3_wp
    call co_interface(x1a, x2a, wa, 1.4_wp, 1.667_wp)

    x1b = [1.0_wp, 1.0_wp, 0.5_wp, 0.0_wp]
    x2b = [1.2_wp, 1.0_wp, 0.1_wp, 0.0_wp]
    wb = 0.3_wp
    call co_interface(x1b, x2b, wb, 1.667_wp, 1.4_wp)

    ! both still solve their own equations to a valid, cross-matched root
    call check(error, x1a(2), x2a(2), thr=1.0e-8_wp, message="run A: p1/=p2")
    if (allocated(error)) return
    call check(error, x1a(3), x2a(3), thr=1.0e-8_wp, message="run A: u1/=u2")
    if (allocated(error)) return
    call check(error, x1b(2), x2b(2), thr=1.0e-8_wp, message="run B: p1/=p2")
    if (allocated(error)) return
    call check(error, x1b(3), x2b(3), thr=1.0e-8_wp, message="run B: u1/=u2")
    if (allocated(error)) return
    ! but swapping which gamma governs which side gives a genuinely
    ! different converged state -- proves both are real, independent
    ! inputs, not one silently ignored or the two aliased together
    call check(error, abs(x1a(1) - x1b(1)) > 0.01_wp,&
    &message="swapping gam1/gam2 did not change the converged density")
    if (allocated(error)) return
    call check(error, abs(x1a(3) - x1b(3)) > 0.005_wp,&
    &message="swapping gam1/gam2 did not change the converged velocity")
    if (allocated(error)) return
  end subroutine test_gamma_independence

end module test_interface
