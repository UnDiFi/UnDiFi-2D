module test_newton_solve
! Phase 7.5 increment 1 (ROADMAP.md #20): coverage for the shared
! damped-Newton engine mod_newton_solve.f90 extracted in Phase 3.2
! (issue #14) that co_utp/co_uqp/co_urr/co_shock/co_dc all now call --
! this is the loop mechanics every later physics-kernel test (co_shock
! etc.) implicitly depends on, so it gets its own coverage first,
! against a synthetic system with a hand-checkable analytic root rather
! than any real R-H physics.
  use mod_kinds, only: wp, i4
  use mod_newton_solve, only: newton_solve, residual_if, jacobian_if
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_newton_solve

contains

  subroutine collect_newton_solve(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("circle_line_analytic_jacobian", test_circle_line_analytic),&
    &new_unittest("circle_line_fd_jacobian", test_circle_line_fd)&
    &]
  end subroutine collect_newton_solve

! r1 = y1**2 + y2**2 - 25  (circle, radius 5, centered at origin)
! r2 = y1 - y2 - 1         (line)
! Two real intersections, (4,3) and (-3,-4); seeded near (4,3) below so
! both the analytic- and FD-Jacobian paths converge to the same one.
  real(wp) function circle_line_residual(i, y, ctx) result(r)
    integer(i4), intent(in) :: i
    real(wp), intent(in) :: y(:)
    class(*), intent(in) :: ctx

    r = 0.0_wp
    if (i == 1) then
      r = y(1)**2 + y(2)**2 - 25.0_wp
    else if (i == 2) then
      r = y(1) - y(2) - 1.0_wp
    end if
  end function circle_line_residual

  subroutine circle_line_jacobian(y, ctx, g)
    real(wp), intent(in) :: y(:)
    class(*), intent(in) :: ctx
    real(wp), intent(out) :: g(:, :)

    g(1, 1) = 2.0_wp*y(1)
    g(1, 2) = 2.0_wp*y(2)
    g(2, 1) = 1.0_wp
    g(2, 2) = -1.0_wp
  end subroutine circle_line_jacobian

  subroutine test_circle_line_analytic(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: y0(2), y(2)
    integer :: ctx

    ctx = 0
    y0 = [3.0_wp, 2.0_wp]
    call newton_solve(2_i4, y0, circle_line_residual, ctx, 1.0_wp, 1.0e-10_wp,&
    &1.0e-6_wp, y, jac=circle_line_jacobian)

    call check(error, y(1), 4.0_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, y(2), 3.0_wp, thr=1.0e-8_wp)
    if (allocated(error)) return
  end subroutine test_circle_line_analytic

! No jac= argument -- exercises newton_solve's central-finite-difference
! Jacobian path (fd_jacobian) instead, the one every solver used before
! Phase 3.5 added analytic Jacobians and still the only path for
! co_utp/co_uqp/co_urr today.
  subroutine test_circle_line_fd(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: y0(2), y(2)
    integer :: ctx

    ctx = 0
    y0 = [3.0_wp, 2.0_wp]
    call newton_solve(2_i4, y0, circle_line_residual, ctx, 1.0_wp, 1.0e-10_wp,&
    &1.0e-6_wp, y)

    call check(error, y(1), 4.0_wp, thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, y(2), 3.0_wp, thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine test_circle_line_fd

end module test_newton_solve
