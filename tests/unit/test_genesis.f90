module test_genesis
! ROADMAP.md item 192 / issue #185: mod_genesis -- shock genesis from the
! singularity-formation and shock-development theory, replacing the tuned
! gradient sensor a capture-to-fit handover otherwise needs.
!
! Everything here is checked against a datum whose preshock is known in closed
! form.  For lambda_0(x) = -sin(x) evolving under Burgers -- and the genuinely
! nonlinear characteristic speed lambda = u + a obeys Burgers exactly in a
! simple wave -- the whole self-similar structure is explicit:
!
!     A = -lambda_0'(0)      = 1            t* = 1/A          = 1
!     B =  lambda_0'''(0)/6  = 1/6          K  = A^2/sqrt(B)  = sqrt(6)
!     E =  lambda_0''''(0)/24 = 0           (odd datum: front sits exactly at 0)
!     min_x lambda_x(t) = -1/(1-t)          so mu(t) = 1-t: linear, slope -1
!     |[lambda]|(t*+dt) = 2 K sqrt(dt)      [s]/c_v ~ dt^(3/2)
!
! so these are manufactured-solution checks against exact values, not
! regression baselines.  The exact solution at t < t* comes from inverting the
! characteristic map by Newton, which is monotone there, so no discretisation of
! the PDE enters and what is being tested is the estimator itself.
!
! The two t* estimators use disjoint information -- the spatial structure at one
! instant versus the time history of one scalar -- so test_agreement is a real
! falsification test rather than a consistency check.
  use mod_kinds, only: wp, i4
  use mod_genesis
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_genesis

  integer(i4), parameter :: NX = 801_i4
  real(wp), parameter :: A_EX = 1.0_wp, B_EX = 1.0_wp/6.0_wp, TS_EX = 1.0_wp
  real(wp), parameter :: K_EX = 2.449489742783178_wp     ! sqrt(6)

contains

  subroutine collect_genesis(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("genesis_selfsim_profile", test_selfsim),&
    &new_unittest("genesis_selfsim_cusp", test_cusp),&
    &new_unittest("genesis_profile_fit", test_profile_fit),&
    &new_unittest("genesis_profile_fit_evolved", test_profile_evolved),&
    &new_unittest("genesis_residual_is_an_error_bar", test_residual_bar),&
    &new_unittest("genesis_history_fit", test_history_fit),&
    &new_unittest("genesis_estimators_agree", test_agreement),&
    &new_unittest("genesis_emit_exponents", test_emit),&
    &new_unittest("genesis_refuses_uncompressed", test_refusal),&
    &new_unittest("genesis_handover_is_consistent", test_handover)&
    &]
  end subroutine collect_genesis

! --- helpers ---------------------------------------------------------------

! exact Burgers solution for lambda_0 = -sin(x): invert x = x0 + t*lambda_0(x0)
  real(wp) function burgers_exact(xt, t) result(lm)
    real(wp), intent(in) :: xt, t
    real(wp) :: x0, g, gp
    integer(i4) :: it
    x0 = xt
    do it = 1, 100
      g = x0 - t*sin(x0) - xt
      gp = 1.0_wp - t*cos(x0)
      if (abs(gp) < 1.0e-300_wp) exit
      x0 = x0 - g/gp
      if (abs(g) < 1.0e-15_wp) exit
    end do
    lm = -sin(x0)
  end function burgers_exact

  subroutine sample(t, x, lam)
    real(wp), intent(in) :: t
    real(wp), intent(out) :: x(NX), lam(NX)
    integer(i4) :: i
    do i = 1, NX
      x(i) = -1.0_wp + 2.0_wp*real(i - 1, wp)/real(NX - 1, wp)
      lam(i) = burgers_exact(x(i), t)
    end do
  end subroutine sample

! --- the attractor ---------------------------------------------------------

! genesis_selfsim must invert y = -U^3 - U exactly, for every y: the depressed
! cubic U^3 + U + y has positive discriminant everywhere, so there is one real
! root and no branch to choose.
  subroutine test_selfsim(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4) :: k
    real(wp) :: y, u
    do k = -40, 40
      y = 0.25_wp*real(k, wp)
      u = genesis_selfsim(y)
      call check(error, y + u**3 + u, 0.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_selfsim

! far field U ~ -y^(1/3): the C^(1/3) preshock cusp, which is what makes the
! newborn front's profile Hoelder-1/3 rather than smooth
  subroutine test_cusp(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: y, u
    y = 1.0e6_wp
    u = genesis_selfsim(-y)
    call check(error, u/y**(1.0_wp/3.0_wp), 1.0_wp, rel=.true., thr=1.0e-4_wp)
  end subroutine test_cusp

! --- the spatial estimator -------------------------------------------------

  subroutine test_profile_fit(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(NX), lam(NX)
    type(genesis_fit_t) :: fit

    call sample(0.0_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fit)

    call check(error, fit%ok)
    if (allocated(error)) return
    call check(error, fit%amod, A_EX, rel=.true., thr=1.0e-8_wp)
    if (allocated(error)) return
    call check(error, fit%bmod, B_EX, rel=.true., thr=1.0e-4_wp)
    if (allocated(error)) return
    call check(error, fit%kmod, K_EX, rel=.true., thr=1.0e-4_wp)
    if (allocated(error)) return
    call check(error, fit%tstar, TS_EX, rel=.true., thr=1.0e-8_wp)
    if (allocated(error)) return
!   the datum is odd about x = 0, so E must come out at zero and the newborn
!   front must sit exactly on the cusp
    call check(error, fit%emod, 0.0_wp, thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, fit%x0, 0.0_wp, thr=1.0e-8_wp)
  end subroutine test_profile_fit

! t* is a PREDICTION: recovered from the field well before the singularity, and
! it must not drift as the field steepens towards it
  subroutine test_profile_evolved(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(NX), lam(NX), t
    type(genesis_fit_t) :: fit
    integer(i4) :: k

    do k = 1, 3
      t = 1.0_wp - 10.0_wp**(-real(k, wp)/3.0_wp)     ! 0.536, 0.785, 0.900
      call sample(t, x, lam)
      call genesis_fit_profile(x, lam, NX, t, fit)
      call check(error, fit%ok)
      if (allocated(error)) return
      call check(error, fit%tstar, TS_EX, rel=.true., thr=1.0e-3_wp)
      if (allocated(error)) return
    end do
  end subroutine test_profile_evolved

! The design claim that matters: the fit RESIDUAL is a usable error bar on t*,
! not merely a diagnostic.  Once the least-squares window no longer fits inside
! the self-similar region -- width K|t-t*|^(3/2) -- the estimate degrades, and
! the estimator has to REFUSE rather than mishandle.  Checked at t = 0.9, where
! the fit is still good, and t = 0.99, where it is not.
  subroutine test_residual_bar(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(NX), lam(NX), err
    type(genesis_fit_t) :: fit

    call sample(0.9_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.9_wp, fit)
    call check(error, fit%ok)
    if (allocated(error)) return
    err = abs(fit%tstar - TS_EX)/TS_EX
!   the residual must bound the actual error to within an order of magnitude --
!   that is what makes it an error bar rather than a number
    call check(error, err < 10.0_wp*fit%resid .and. fit%resid < 100.0_wp*err)
    if (allocated(error)) return

!   window far outside the self-similar region: must be refused, not returned
    call sample(0.99_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.99_wp, fit)
    call check(error, .not. fit%ok)
  end subroutine test_residual_bar

! --- the temporal estimator ------------------------------------------------

! mu(t) = -1/min_x lambda_x decays LINEARLY to zero at t*.  The slope is a
! prediction, not a fitted parameter: it must come out at -1, and its departure
! from -1 is an independent diagnostic of the fit window.
  subroutine test_history_fit(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4), parameter :: NT = 12_i4
    real(wp) :: t(NT), lxmin(NT)
    type(genesis_fit_t) :: fit
    integer(i4) :: k

    do k = 1, NT
      t(k) = 0.05_wp*real(k, wp)
      lxmin(k) = -1.0_wp/(1.0_wp - t(k))
    end do
    call genesis_fit_history(t, lxmin, NT, fit)

    call check(error, fit%ok)
    if (allocated(error)) return
    call check(error, fit%tstar, TS_EX, rel=.true., thr=1.0e-10_wp)
    if (allocated(error)) return
    call check(error, fit%slope, -1.0_wp, thr=1.0e-10_wp)
  end subroutine test_history_fit

! --- the falsification test ------------------------------------------------

! The two estimators share no information, so agreement within their own
! residuals is evidence and disagreement beyond them is a refusal.  Both halves
! are checked: a genuine preshock must agree, and a deliberately wrong history
! must not.
  subroutine test_agreement(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4), parameter :: NT = 12_i4
    real(wp) :: x(NX), lam(NX), t(NT), lxmin(NT)
    type(genesis_fit_t) :: fp, fh
    integer(i4) :: k

    call sample(0.0_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fp)
    do k = 1, NT
      t(k) = 0.05_wp*real(k, wp)
      lxmin(k) = -1.0_wp/(1.0_wp - t(k))
    end do
    call genesis_fit_history(t, lxmin, NT, fh)
    call check(error, genesis_fit_agree(fp, fh))
    if (allocated(error)) return

!   a history whose singularity is at t* = 2 instead of 1 must be rejected,
!   even though each fit is individually perfect
    do k = 1, NT
      lxmin(k) = -1.0_wp/(2.0_wp - t(k))
    end do
    call genesis_fit_history(t, lxmin, NT, fh)
    call check(error, .not. genesis_fit_agree(fp, fh))
  end subroutine test_agreement

! --- the operator ----------------------------------------------------------

! The exponents are the point: strength ~ dt^(1/2) with constant 2K, entropy
! ~ dt^(3/2) with the Landau-Lifshitz constant (gamma^2-1)/(12 gamma^2).  Both
! are checked as CONSTANTS over four decades, which tests the exponent, not just
! the value at one point.
  subroutine test_emit(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp), parameter :: GAM = 1.4_wp, RHO0 = 1.0_wp, P0 = 1.0_wp
    real(wp) :: x(NX), lam(NX), dt, c1, c1ref, c3, c3ref
    type(genesis_fit_t) :: fit
    type(genesis_data_t) :: gd
    integer(i4) :: k

    call sample(0.0_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fit)
    call check(error, fit%ok)
    if (allocated(error)) return

    c1ref = 0.0_wp
    c3ref = 0.0_wp
    do k = 1, 5
      dt = 10.0_wp**(-real(k, wp))
      call genesis_emit(fit, GAM, RHO0, P0, dt, gd)
      call check(error, gd%ok)
      if (allocated(error)) return

!     |[lambda]| / (2 K sqrt(dt))  must be 1, independent of dt
      c1 = gd%dlam/(2.0_wp*K_EX*sqrt(dt))
      call check(error, c1, 1.0_wp, rel=.true., thr=1.0e-4_wp)
      if (allocated(error)) return

!     [s]/c_v / dt^(3/2) must be a constant, independent of dt
      c3 = gd%ds_cv/dt**1.5_wp
      if (k == 1) then
        c3ref = c3
      else
        call check(error, c3, c3ref, rel=.true., thr=1.0e-8_wp)
        if (allocated(error)) return
      end if

!     odd datum: the newborn front sits on the cusp, to the order carried
      call check(error, gd%xs, 0.0_wp, thr=1.0e-8_wp)
      if (allocated(error)) return
!     and it moves at the local characteristic speed
      call check(error, gd%w, fit%lam0, rel=.true., thr=1.0e-12_wp)
      if (allocated(error)) return
      if (c1ref == 0.0_wp) c1ref = c1
    end do
  end subroutine test_emit

! A field with no compression anywhere has no shock forming, and saying so is
! the correct answer.  Likewise a monotone-increasing lambda: the estimator must
! not manufacture a birth event out of an expansion wave.
  subroutine test_refusal(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(NX), lam(NX)
    type(genesis_fit_t) :: fit
    integer(i4) :: i

    do i = 1, NX
      x(i) = -1.0_wp + 2.0_wp*real(i - 1, wp)/real(NX - 1, wp)
      lam(i) = 3.0_wp + 0.5_wp*x(i)            ! expansion: lambda' > 0 everywhere
    end do
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fit)
    call check(error, .not. fit%ok)
    if (allocated(error)) return

!   uniform field: nothing at all
    lam = 3.0_wp
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fit)
    call check(error, .not. fit%ok)
  end subroutine test_refusal

! The handover time is where the expansion's own truncation error, relative
! O(A dt), meets the floor co_shock would leave behind, c_tau/(c_p sqrt(dt)).
! Consistency check: at dt*, those two must be equal.  That is the whole content
! of the formula, and it is checked rather than asserted.
  subroutine test_handover(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp), parameter :: GAM = 1.4_wp, RHO0 = 1.0_wp, P0 = 1.0_wp
    real(wp) :: x(NX), lam(NX), dts, a0, cp, e_expansion, e_solver
    type(genesis_fit_t) :: fit

    call sample(0.0_wp, x, lam)
    call genesis_fit_profile(x, lam, NX, 0.0_wp, fit)
    dts = genesis_handover_time(fit, GAM, RHO0, P0)
    call check(error, dts > 0.0_wp)
    if (allocated(error)) return

    a0 = sqrt(GAM*P0/RHO0)
    cp = 4.0_wp*RHO0*a0*fit%kmod/(GAM + 1.0_wp)
    e_expansion = fit%amod*dts                       ! grows with dt
    e_solver = genesis_ctau_ref/(cp*sqrt(dts))       ! shrinks with dt
    call check(error, e_expansion, e_solver, rel=.true., thr=1.0e-10_wp)
  end subroutine test_handover

end module test_genesis
