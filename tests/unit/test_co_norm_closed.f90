module test_co_norm_closed
! Issue #22 (E2): co_norm_closed -- cyclic-indexed normal computation
! for a closed discontinuity curve, on a full circle (the canonical
! closed-curve case: a circular contact discontinuity, matching issue
! #22's own suggested test "a closed circular contact discontinuity
! advected in a uniform flow must retain its shape and produce no
! spurious signal" -- this increment checks the NORMAL computation
! itself, the first building block that test would depend on, not the
! full advection scenario).
!
! Same "quiescent state forces the central-difference branch" technique
! as test_co_norm.f90 (issue #20): with u=v=0 everywhere,
! shp_dpndnc/dcp_dpndnc's dependency test provably returns 1 (dependent)
! for any two distinct points (see that file's header for the algebra),
! so every point here takes the plain weighted-central-tangent branch,
! isolating the cyclic INDEX arithmetic itself as the thing under test
! -- the whole point of this increment -- rather than also having to
! reconstruct the upwind-selection logic on a closed loop.
  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, naddholesmax, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_norm_closed

  real(wp), parameter :: pi = 3.14159265358979323846_wp

contains

  subroutine collect_co_norm_closed(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("circle_normals_point_radially_outward", test_circle),&
    &new_unittest("wraps_at_the_seam_no_endpoint_artifact", test_seam)&
    &]
  end subroutine collect_co_norm_closed

! Builds a full circle of n points (NOT closed by repeating point 1 as
! point n+1 -- co_norm_closed's own cyclic indexing handles the wrap,
! so point n's forward neighbor is point 1 directly) with a quiescent
! (rho=1, p=1, u=v=0) state, calls co_norm_closed, returns the normals.
  subroutine run_on_circle(n, r, theta0, dtheta_pattern, vshnor)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: r, theta0
    real(wp), intent(in) :: dtheta_pattern(n) ! per-step angular increment
    real(wp), intent(out) :: vshnor(2, n)

    include 'paramt.h'
    external :: co_norm_closed

    real(wp) :: xysh(ndim, n), zroesh(ndof, n), theta(n)
    integer(i4) :: i

    ga = 1.4_wp
    gm1 = 0.4_wp

    theta(1) = theta0
    do i = 2, n
      theta(i) = theta(i - 1) + dtheta_pattern(i)
    end do

    do i = 1, n
      xysh(1, i) = r*cos(theta(i))
      xysh(2, i) = r*sin(theta(i))
      zroesh(1, i) = 1.0_wp    ! sqrt(rho), rho=1
      zroesh(2, i) = ga/gm1    ! => p = gm1/ga*(1*ga/gm1 - 0) = 1
      zroesh(3, i) = 0.0_wp
      zroesh(4, i) = 0.0_wp
    end do

    call co_norm_closed(xysh, zroesh, vshnor, n, 'S')
  end subroutine run_on_circle

! Uniform spacing: on a circle, symmetric equal-angle neighbors make the
! weighted-central formula an EXACT geometric identity for the tangent
! (same reasoning as test_co_norm.f90's own circle case), so normals
! should point exactly radially outward at every one of the n points,
! all the way around including the seam between point n and point 1.
  subroutine test_circle(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4), parameter :: n = 40
    real(wp) :: dtheta(n), vshnor(2, n), theta_i, expect_nx, expect_ny
    integer(i4) :: i

    dtheta = 2.0_wp*pi/real(n, wp)
    call run_on_circle(n, 2.0_wp, 0.0_wp, dtheta, vshnor)

    do i = 1, n
      theta_i = real(i - 1, wp)*2.0_wp*pi/real(n, wp)
      ! co_norm's own convention: vshnor=(tauy,-taux); for a CCW circle
      ! (increasing i = increasing theta) the tangent is
      ! (-sin,cos), so the normal is (cos,sin) -- radially outward.
      expect_nx = cos(theta_i)
      expect_ny = sin(theta_i)
      call check(error, vshnor(1, i), expect_nx, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, vshnor(2, i), expect_ny, thr=1.0e-10_wp)
      if (allocated(error)) return
    end do
  end subroutine test_circle

! Non-uniform spacing (alternating short/long steps, same pattern as
! test_co_norm.f90's own non-uniform circle case) -- checks specifically
! AT and AROUND the seam (points n, 1, 2) where an open-curve
! implementation would have hit its i==1/i==nshockpoints special-casing.
! Reuses the SAME "chord-length**2-weighted formula is exact for any 3
! co-circular points" identity already established in test_co_norm.f90
! -- the seam is not special here BECAUSE it is not special in the
! implementation (no endpoint branch exists to treat it differently).
  subroutine test_seam(error)
    type(error_type), allocatable, intent(out) :: error
    integer(i4), parameter :: n = 21
    real(wp) :: dtheta(n), vshnor(2, n), theta(n), expect_nx, expect_ny
    integer(i4) :: i, k, seam_idx(3)

    do i = 1, n
      dtheta(i) = merge(1.3_wp, 0.7_wp, mod(i, 2_i4) == 0_i4)*(2.0_wp*pi/real(n, wp))
    end do
    call run_on_circle(n, 1.5_wp, 0.4_wp, dtheta, vshnor)

    theta(1) = 0.4_wp
    do i = 2, n
      theta(i) = theta(i - 1) + dtheta(i)
    end do

    ! points n, 1, 2 straddle the seam where an open-curve implementation
    ! would have hit its i==1/i==nshockpoints special-casing
    seam_idx = [n, 1, 2]
    do k = 1, 3
      i = seam_idx(k)
      expect_nx = cos(theta(i)); expect_ny = sin(theta(i))
      call check(error, vshnor(1, i), expect_nx, thr=1.0e-10_wp)
      if (allocated(error)) return
      call check(error, vshnor(2, i), expect_ny, thr=1.0e-10_wp)
      if (allocated(error)) return
    end do
  end subroutine test_seam

end module test_co_norm_closed
