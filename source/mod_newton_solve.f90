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
!
! Phase 3.5 increment 1 (ROADMAP.md #14): optional analytic Jacobian
! (`jac`), replacing the FD double-loop below when supplied. Absent (the
! default -- still true for co_dc/co_urr/co_uqp/co_utp until their own
! increments add a `jac`), behavior is byte-for-byte identical to before.
! `verify_jac` is a permanent, opt-in runtime cross-check: when a `jac`
! is supplied, it also computes the FD Jacobian purely for comparison and
! warns on disagreement, while still stepping with the analytic one --
! meant to validate a new analytic Jacobian's derivation against real
! solver traffic (many distinct iterates, not just a few hand-picked
! sample points) before trusting it, then left available as a standing
! safety net for any later analytic Jacobian added the same way.

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: newton_solve, residual_if, jacobian_if

  real(wp), parameter :: fd_eps_floor = 1.0e-7_wp
! jac_check_* is deliberately loose, not a numerics tolerance: the FD
! comparator itself carries roundoff noise (tiny fd_eps_floor-scale
! perturbations divided back out) that a genuinely correct analytic
! Jacobian will still disagree with at the ~1e-4 relative level on
! entries with near-zero arguments -- this only needs to catch a real
! derivation mistake (wrong sign, wrong term, swapped index), which
! shows up as order-1 relative disagreement, not roundoff-level noise.
  real(wp), parameter :: jac_check_rtol = 1.0e-4_wp, jac_check_atol = 1.0e-7_wp

  abstract interface
    function residual_if(i, y, ctx) result(r)
      import :: wp, i4
      integer(i4), intent(in) :: i
      real(wp), intent(in) :: y(:)
      class(*), intent(in) :: ctx
      real(wp) :: r
    end function residual_if

    subroutine jacobian_if(y, ctx, g)
      import :: wp
      real(wp), intent(in) :: y(:)
      class(*), intent(in) :: ctx
      real(wp), intent(out) :: g(:, :)
    end subroutine jacobian_if
  end interface

contains

  subroutine newton_solve(n, y0, resid, ctx, damping, tol, fd_eps_rel, y,&
  &ifail, log_unit, log_iter, detect_divergence, maxiter, jac, verify_jac)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: y0(n), damping, tol, fd_eps_rel
    procedure(residual_if) :: resid
    class(*), intent(in) :: ctx
    real(wp), intent(out) :: y(n)
    logical, intent(out), optional :: ifail
    integer(i4), intent(in), optional :: log_unit
    logical, intent(in), optional :: log_iter, detect_divergence
    integer(i4), intent(in), optional :: maxiter
    procedure(jacobian_if), optional :: jac
    logical, intent(in), optional :: verify_jac
    external solg

    real(wp) :: yn(n), yn1(n), g(n, n), gfd(n, n), bb(n), dyn(n)
    real(wp) :: dum, dumold
    integer(i4) :: i, j, icont, maxiter_eff
    logical :: want_divergence, want_log_iter, want_verify_jac

    want_divergence = .false.
    if (present(detect_divergence)) want_divergence = detect_divergence
    want_log_iter = .false.
    if (present(log_iter)) want_log_iter = log_iter
    if (present(ifail)) ifail = .false.
    maxiter_eff = 500_i4
    if (present(maxiter)) maxiter_eff = maxiter
    want_verify_jac = .false.
    if (present(verify_jac)) want_verify_jac = verify_jac

    yn1 = y0
    icont = 0

    do
      do i = 1, n
        yn(i) = yn1(i)
        bb(i) = resid(i, yn1, ctx)
      end do

      if (present(jac)) then
        call jac(yn, ctx, g)
        if (want_verify_jac) then
          call fd_jacobian(n, yn, resid, ctx, fd_eps_rel, gfd)
          do i = 1, n
            do j = 1, n
              if (abs(g(i, j) - gfd(i, j)) .gt.&
              &jac_check_atol + jac_check_rtol*abs(gfd(i, j))) then
                if (present(log_unit)) then
                  write (log_unit, *) 'JACOBIAN MISMATCH', i, j, g(i, j), gfd(i, j)
                end if
              end if
            end do
          end do
        end if
      else
! jacobian: central finite differences, relative perturbation fd_eps_rel
! of the current iterate, floored at fd_eps_floor absolute
        call fd_jacobian(n, yn, resid, ctx, fd_eps_rel, g)
      end if

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

! Safety net: none of the five pre-refactor solvers had one (the
! do/exit-only-on-convergence loop was extracted verbatim), and real
! end-to-end regression testing -- only possible once the test harness
! itself was fixed to actually run the CMake build (see this session's
! earlier commits) -- found two independent, unrelated inputs where
! that assumption was false: co_shock hanging on NACA0012_M080_A0
! (separately root-caused and fixed: an uninitialized Newton seed) and
! co_dc hanging on SSInteractions2-2 (root cause not yet isolated).
! Bail out with ifail after maxiter_eff iterations rather than loop
! forever; the four callers that don't check ifail (only co_utp does,
! for its own divergence-triggered retry) simply proceed with
! whatever iterate was reached, no worse than hanging and no different
! from what they already do the instant tol is satisfied.
      if (icont .ge. maxiter_eff) then
        if (present(ifail)) ifail = .true.
        exit
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

! Central-finite-difference Jacobian at the base point yn, relative
! perturbation fd_eps_rel of the current iterate, floored at
! fd_eps_floor absolute. Factored out of newton_solve's main loop so it
! can be called either as the sole Jacobian source (no `jac` supplied)
! or, when `verify_jac` is set, purely as a cross-check against an
! analytic `jac` -- same computation either way, not a behavior change
! from the pre-3.5 inline double loop.
  subroutine fd_jacobian(n, yn, resid, ctx, fd_eps_rel, g)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: yn(n), fd_eps_rel
    procedure(residual_if) :: resid
    class(*), intent(in) :: ctx
    real(wp), intent(out) :: g(n, n)

    real(wp) :: yp(n)
    real(wp) :: dyn1, dum1, dum2
    integer(i4) :: i, j

    do i = 1, n
      do j = 1, n
        yp = yn
        dyn1 = abs(yn(j))*fd_eps_rel
        if (dyn1 .le. fd_eps_floor) dyn1 = fd_eps_floor
        yp(j) = yn(j) + dyn1
        dum2 = resid(i, yp, ctx)
        yp(j) = yn(j) - dyn1
        dum1 = resid(i, yp, ctx)
        g(i, j) = (dum2 - dum1)/(2.0_wp*dyn1)
      end do
    end do
  end subroutine fd_jacobian

end module mod_newton_solve
