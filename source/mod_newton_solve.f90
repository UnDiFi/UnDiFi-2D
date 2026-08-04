module mod_newton_solve
! Phase 3.2 increment 1 (ROADMAP.md #14): the one damped-Newton +
! central-finite-difference-Jacobian + solg (fnd_phps.f90) linear-solve
! loop that co_utp, co_uqp, co_urr, co_shock, and co_dc each hand-
! duplicated, with different damping factors, convergence tolerances,
! and finite-difference perturbation scales. This module extracts the
! loop mechanics only -- each caller still builds its own residual
! function and (for co_utp/co_uqp/co_urr) its own linear "known value"
! constraint rows; those are per-point-type physics/setup, not loop
! duplication, and stay in each caller's own file.
!
! Per-solver constants, verified directly from the pre-refactor source
! (this is a mechanical extraction, not a behavior change):
!   co_utp:   n=20, damping=0.5, tol=1e-7,  fd_eps_rel=0.001, and the
!             only one of the five with in-loop divergence detection
!             (residual growing between iterations -> ifail=.true.)
!   co_uqp:   n=24, damping=0.5, tol=1e-11, fd_eps_rel=0.001
!   co_urr:   n=15, damping=1.0 (no damping), tol=1e-6,  fd_eps_rel=0.01
!   co_shock: n=4,  damping=0.2, tol=1e-7,  fd_eps_rel=0.001
!   co_dc:    n=7,  damping=1.0 (no damping), tol=1e-10, fd_eps_rel=0.01
! The absolute perturbation floor (1e-7) is common to all five; the
! pre-refactor floor comparison was ".lt." in co_utp/co_uqp/co_urr and
! ".le." in co_shock/co_dc -- an immaterial difference (only matters at
! exact floating-point equality to 1e-7, never hit by a computed
! abs(y)*rel_eps in practice), standardized here to ".le.".

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: newton_solve, residual_if

  real(wp), parameter :: fd_eps_floor = 1.0e-7_wp

  abstract interface
    function residual_if(i, y, ctx) result(r)
      import :: wp, i4
      integer(i4), intent(in) :: i
      real(wp), intent(in) :: y(:)
      class(*), intent(in) :: ctx
      real(wp) :: r
    end function residual_if
  end interface

contains

  subroutine newton_solve(n, y0, resid, ctx, damping, tol, fd_eps_rel, y,&
  &ifail, log_unit, log_iter, detect_divergence)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: y0(n), damping, tol, fd_eps_rel
    procedure(residual_if) :: resid
    class(*), intent(in) :: ctx
    real(wp), intent(out) :: y(n)
    logical, intent(out), optional :: ifail
    integer(i4), intent(in), optional :: log_unit
    logical, intent(in), optional :: log_iter, detect_divergence
    external solg

    real(wp) :: yn(n), yn1(n), g(n, n), bb(n), dyn(n)
    real(wp) :: dyn1, dum1, dum2, dum, dumold
    integer(i4) :: i, j, k, icont
    logical :: want_divergence, want_log_iter

    want_divergence = .false.
    if (present(detect_divergence)) want_divergence = detect_divergence
    want_log_iter = .false.
    if (present(log_iter)) want_log_iter = log_iter
    if (present(ifail)) ifail = .false.

    yn1 = y0
    icont = 0

    do
      do i = 1, n
        yn(i) = yn1(i)
        bb(i) = resid(i, yn1, ctx)
      end do

! jacobian: central finite differences, relative perturbation fd_eps_rel
! of the current iterate, floored at fd_eps_floor absolute
      do i = 1, n
        do j = 1, n
          do k = 1, n
            yn1(k) = yn(k)
          end do
          dyn1 = abs(yn1(j))*fd_eps_rel
          if (dyn1 .le. fd_eps_floor) dyn1 = fd_eps_floor
          yn1(j) = yn(j) + dyn1
          dum2 = resid(i, yn1, ctx)
          yn1(j) = yn(j) - dyn1
          dum1 = resid(i, yn1, ctx)
          g(i, j) = (dum2 - dum1)/(2.0_wp*dyn1)
        end do
      end do

      call solg(n, n, g, bb, dyn)
      do i = 1, n
        yn1(i) = yn(i) - damping*dyn(i)
      end do

! residual: L1 norm of the step (not of bb) -- matches all five
! pre-refactor solvers' convergence check
      dum = 0.0_wp
      do i = 1, n
        dum = dum + abs(yn1(i) - yn(i))
      end do
      icont = icont + 1

      if (present(log_unit)) then
        if (want_log_iter) then
          write (log_unit, *) 'conv--->', dum, icont
        else
          write (log_unit, *) 'conv--->', dum
        end if
      end if

      if (want_divergence) then
! co_utp only: bail out (ifail=.true.) if the step norm grows between
! iterations, instead of looping forever on a diverging sequence
        if (icont .eq. 1) then
          dumold = dum
          cycle
        end if
        if (dum .gt. dumold) then
          if (present(ifail)) ifail = .true.
          y = yn1
          return
        end if
        if (dum .gt. tol) then
          dumold = dum
          cycle
        end if
        exit
      else
        if (dum .le. tol) exit
      end if
    end do

    y = yn1
  end subroutine newton_solve

end module mod_newton_solve
