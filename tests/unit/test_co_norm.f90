module test_co_norm
! Phase 7.5 increment 3 (ROADMAP.md #20): co_norm -- the upwind-biased
! finite-difference unit-normal-vector calculation for shock/discontinuity
! points -- against analytic curves (straight line, circle, parabola), all
! at non-uniform point spacing. The line and circle cases turn out to be
! EXACT identities for this formula (see test_line/circle_tangent_error
! below for why); the parabola case is the one that actually exercises
! discretization error, checked for the claimed 2nd-order convergence.
!
! co_norm's per-point normal comes from an upwind-biased choice between a
! plain two-sided weighted-central-difference tangent estimate and two
! different one-sided estimates, the choice driven per-side by
! shp_dpndnc/dcp_dpndnc's local domain-of-dependence test on the flow
! state at the neighboring point. Every state in this suite is set
! deliberately QUIESCENT (u=v=0 everywhere): with uu=vv=0, shp_dpndnc's
! `dum` reduces to aa**2*(xx**2+yy**2), positive for any two distinct
! points, and its quadratic's two roots reduce to dt1=-sqrt(dum)/aa**2
! (always negative) and dt2=+sqrt(dum)/aa**2 (always positive) -- so
! shp_dpndnc always returns 1 (dependent) for both neighbors of every
! interior point, meaning `depim1*depip1` is always 1 and the one-sided
! branches in co_norm.f90 never execute. This isolates exactly the plain
! FD tangent estimator the issue asks to check, without also having to
! reverse-engineer the upwind-selection logic to hit it.
  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax, nshmax, naddholesmax, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_norm

contains

  subroutine collect_co_norm(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("straight_line_exact", test_line),&
    &new_unittest("circle_nonuniform_exact", test_circle_exact),&
    &new_unittest("parabola_nonuniform_convergence", test_parabola_convergence)&
    &]
  end subroutine collect_co_norm

! Runs co_norm on a single shock curve of n points (all other production
! arguments set to their "nothing else going on" values: one shock,
! nspecpoints=0 so sp%correct_normal is never dispatched) and returns the
! computed unit normals. Every point is given the same quiescent state
! (rho=1, p=1, u=v=0) in co_norm's expected Roe sqrt-variable layout
! (zroe(1)=sqrt(rho), zroe(2) chosen so p=gm1/ga*(zroe(1)*zroe(2)-half of
! the -- zero here -- momentum term), zroe(3)=zroe(4)=0).
  subroutine run_co_norm_on_curve(x, y, n, vshnor)
    real(wp), intent(in) :: x(:), y(:)
    integer(i4), intent(in) :: n
    real(wp), intent(out) :: vshnor(2, n)

    include 'paramt.h'
    external :: co_norm

    real(wp), allocatable :: xysh(:, :, :), zroesh(:, :, :), zroeshu(:, :, :)
    real(wp), allocatable :: vshnor_full(:, :, :)
    integer(i4) :: nshockpoints(nshmax)
    character(len=1) :: typesh(nshmax)
    character(len=5) :: typespecpoints(1)
    integer(i4) :: shinspps(2, 5, 1), ispclr(5, 1)
    integer(i4) :: ia(1), ja(1), iclr(1)
    real(wp) :: corg(ndim, 1)
    integer(i4) :: i

    allocate (xysh(ndim, npshmax, nshmax), zroesh(ndof, npshmax, nshmax),&
    &zroeshu(ndof, npshmax, nshmax), vshnor_full(ndim, npshmax, nshmax))

    ga = 1.4_wp
    gm1 = 0.4_wp

    xysh = 0.0_wp
    zroesh = 0.0_wp
    zroeshu = 0.0_wp
    vshnor_full = 0.0_wp

    do i = 1, n
      xysh(1, i, 1) = x(i)
      xysh(2, i, 1) = y(i)
      zroesh(1, i, 1) = 1.0_wp    ! sqrt(rho), rho=1
      zroesh(2, i, 1) = ga/gm1    ! => p = gm1/ga*(1*ga/gm1 - 0) = 1
      zroesh(3, i, 1) = 0.0_wp    ! sqrt(rho)*u, u=0
      zroesh(4, i, 1) = 0.0_wp    ! sqrt(rho)*v, v=0
    end do

    nshockpoints = 0
    nshockpoints(1) = n
    typesh = ' '
    typesh(1) = 'S'

    call co_norm(xysh, zroeshu, zroesh, vshnor_full, 1_i4, nshockpoints,&
    &typesh, 0_i4, typespecpoints, shinspps, ispclr, ia, ja, iclr, 1_i4, corg)

    vshnor(:, 1:n) = vshnor_full(:, 1:n, 1)
  end subroutine run_co_norm_on_curve

! A straight line's secant direction is exactly the same on both sides of
! any interior point, for ANY spacing (uniform or not) -- so the weighted
! formula's output is exact, not just asymptotically accurate. Uses
! deliberately non-uniform spacing (s_i = i + 0.15*sin(i)) so this isn't
! secretly also a uniform-spacing special case.
  subroutine test_line(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4), parameter :: n = 11
    real(wp) :: s(n), x(n), y(n), vshnor(2, n)
    real(wp) :: alpha, cosa, sina
    integer(i4) :: i

    alpha = 0.7_wp
    cosa = cos(alpha)
    sina = sin(alpha)
    do i = 1, n
      s(i) = real(i, wp) + 0.15_wp*sin(real(i, wp))
      x(i) = s(i)*cosa
      y(i) = s(i)*sina
    end do

    call run_co_norm_on_curve(x, y, n, vshnor)

    ! s is strictly increasing (0.15 perturbation can't reverse unit
    ! steps), so every secant points in the +(cosa,sina) direction and the
    ! resulting normal is unambiguously (sina, -cosa) at every interior
    ! point, to machine precision.
    do i = 4, n - 3
      call check(error, vshnor(1, i), sina, thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, vshnor(2, i), -cosa, thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_line

! Builds n points on a circle of radius r at angles theta0 + sum of
! non-uniform increments (alternating short/long steps, ratio 13:7) and
! returns the FD-vs-analytic angular error of the unit tangent at the
! fixed interior index i0, for two very different step sizes h.
!
! Numerically verified (see this increment's commit message / session
! notes) rather than assumed: the chord-length**2-weighted formula turns
! out to be EXACT for any 3 points on a common circle, independent of how
! asymmetric the angular spacing is -- not merely 2nd-order accurate. This
! is a genuine circle-specific geometric identity (the same tangent-chord
! angle relationship behind the inscribed-angle theorem), not a property
! this formula has for a general smooth curve. A first version of this
! test tried to measure an O(h**2) convergence RATE here, as for the
! parabola below; the "errors" it measured were pure floating-point
! roundoff (~1e-15, not decaying in any h-dependent way, confirmed with a
! standalone check across step sizes spanning 0.001 to 1.0 rad), i.e. it
! was measuring noise, not discretization error. Checking exactness
! directly, at two unrelated step sizes, is both the more robust test and
! the true property.
  subroutine circle_tangent_error(h, err)
    real(wp), intent(in) :: h
    real(wp), intent(out) :: err

    integer(i4), parameter :: n = 21, i0 = 11
    real(wp) :: r, theta0, theta(n), x(n), y(n), vshnor(2, n)
    real(wp) :: tangent_exact(2), tangent_fd(2), dtheta
    integer(i4) :: i

    r = 2.0_wp
    theta0 = 0.3_wp
    theta(i0) = theta0
    do i = i0 + 1, n
      dtheta = merge(1.3_wp*h, 0.7_wp*h, mod(i, 2_i4) == 0_i4)
      theta(i) = theta(i - 1) + dtheta
    end do
    do i = i0 - 1, 1, -1
      dtheta = merge(1.3_wp*h, 0.7_wp*h, mod(i, 2_i4) == 0_i4)
      theta(i) = theta(i + 1) - dtheta
    end do
    do i = 1, n
      x(i) = r*cos(theta(i))
      y(i) = r*sin(theta(i))
    end do

    call run_co_norm_on_curve(x, y, n, vshnor)

    ! analytic unit tangent at theta0 (counterclockwise, increasing i);
    ! co_norm's normal is the tangent rotated by (tauy,-taux), so the
    ! matching analytic normal is (sin theta0, -cos theta0) -- compare
    ! angularly via the cross product (small-angle sine) rather than
    ! componentwise, since only the FD tangent's DIRECTION is meaningful.
    tangent_exact = [-sin(theta0), cos(theta0)]
    tangent_fd = [-vshnor(2, i0), vshnor(1, i0)]  ! undo the -90 deg rotation
    err = abs(tangent_exact(1)*tangent_fd(2) - tangent_exact(2)*tangent_fd(1))
  end subroutine circle_tangent_error

  subroutine test_circle_exact(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: err_coarse, err_fine

    call circle_tangent_error(0.5_wp, err_coarse)
    call circle_tangent_error(0.02_wp, err_fine)

    call check(error, err_coarse, 0.0_wp, thr=1.0e-10_wp,&
    &message="FD tangent not exact on a circle at coarse non-uniform spacing")
    if (allocated(error)) return
    call check(error, err_fine, 0.0_wp, thr=1.0e-10_wp,&
    &message="FD tangent not exact on a circle at fine non-uniform spacing")
    if (allocated(error)) return
  end subroutine test_circle_exact

! Same idea as circle_tangent_error but for a parabola y=a*x**2, sampled
! at non-uniformly spaced x (alternating short/long steps), checking the
! FD tangent against the analytic derivative dy/dx=2*a*x at a fixed
! interior point.
  subroutine parabola_tangent_error(h, err)
    real(wp), intent(in) :: h
    real(wp), intent(out) :: err

    integer(i4), parameter :: n = 21, i0 = 11
    real(wp), parameter :: a_coeff = 0.5_wp
    real(wp) :: x(n), y(n), vshnor(2, n), dx
    real(wp) :: tangent_exact(2), tangent_fd(2), slope, tlen
    integer(i4) :: i

    x(i0) = 0.4_wp
    do i = i0 + 1, n
      dx = merge(1.3_wp*h, 0.7_wp*h, mod(i, 2_i4) == 0_i4)
      x(i) = x(i - 1) + dx
    end do
    do i = i0 - 1, 1, -1
      dx = merge(1.3_wp*h, 0.7_wp*h, mod(i, 2_i4) == 0_i4)
      x(i) = x(i + 1) - dx
    end do
    do i = 1, n
      y(i) = a_coeff*x(i)**2
    end do

    call run_co_norm_on_curve(x, y, n, vshnor)

    slope = 2.0_wp*a_coeff*x(i0)
    tlen = sqrt(1.0_wp + slope**2)
    tangent_exact = [1.0_wp/tlen, slope/tlen]
    tangent_fd = [-vshnor(2, i0), vshnor(1, i0)]
    err = abs(tangent_exact(1)*tangent_fd(2) - tangent_exact(2)*tangent_fd(1))
  end subroutine parabola_tangent_error

  subroutine test_parabola_convergence(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: err_coarse, err_fine, ratio

    call parabola_tangent_error(0.08_wp, err_coarse)
    call parabola_tangent_error(0.04_wp, err_fine)

    ratio = err_coarse/err_fine
    call check(error, err_fine < err_coarse,&
    &message="finer spacing did not reduce the tangent error")
    if (allocated(error)) return
    call check(error, ratio > 2.5_wp .and. ratio < 6.0_wp,&
    &message="convergence ratio not consistent with 2nd-order accuracy")
    if (allocated(error)) return
  end subroutine test_parabola_convergence

end module test_co_norm
