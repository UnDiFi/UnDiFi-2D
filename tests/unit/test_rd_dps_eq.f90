module test_rd_dps_eq
! Phase 7.5 increment 8 (ROADMAP.md #20): rd_dps_eq -- full uniform
! re-parametrization of a shock's points by arc length -- idempotency
! and spacing-bound checks.
!
! Unlike rd_dps (at most one insertion/removal per call), rd_dps_eq
! rebuilds the ENTIRE point set at uniform arc-length spacing dxcell in
! one call: new point count = floor(total_length/dxcell)+1, new spacing
! = total_length/(new_count-1), every interior point linearly
! interpolated (position AND state) between the two bracketing original
! points by arc length. Both tests here use points laid out along the
! x-axis (y=0), so arc length equals x-coordinate exactly -- after
! redistribution the new x-coordinates must equal the new uniform
! abscissa values exactly, a very direct, self-checking property that
! doesn't depend on trusting a separate hand (or Python) reference
! computation for the geometry, only for the interpolated state values.
! The exact bracket-search index arithmetic (which original edge each
! new point's arc-length falls into, including the boundary case where a
! new abscissa lands EXACTLY on an old point) was verified against a
! standalone Python simulation before writing the Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_rd_dps_eq

contains

  subroutine collect_rd_dps_eq(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("uniform_input_is_idempotent", test_idempotent),&
    &new_unittest("nonuniform_input_becomes_uniform", test_redistribute)&
    &]
  end subroutine collect_rd_dps_eq

! Same run/tracer convention as test_rd_dps.f90's run_rd_dps: points
! along the x-axis, zroeshd dof 1 = 10*original point index.
  subroutine run_rd_dps_eq(x, n, dxcell_in, xysh_out, tracer_out, nout)
    real(wp), intent(in) :: x(:)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: dxcell_in
    real(wp), intent(out) :: xysh_out(:), tracer_out(:)
    integer(i4), intent(out) :: nout

    include 'paramt.h'
    external :: rd_dps_eq

    real(wp), allocatable :: xysh(:, :, :), zroeshu(:, :, :), zroeshd(:, :, :)
    integer(i4) :: nshockpoints(nshmax), nshockedges(nshmax)
    integer(i4) :: i

    allocate (xysh(ndim, npshmax, nshmax), zroeshu(ndof, npshmax, nshmax),&
    &zroeshd(ndof, npshmax, nshmax))
    xysh = 0.0_wp; zroeshu = 0.0_wp; zroeshd = 0.0_wp

    dxcell = dxcell_in

    do i = 1, n
      xysh(1, i, 1) = x(i)
      xysh(2, i, 1) = 0.0_wp
      zroeshd(1, i, 1) = 10.0_wp*real(i, wp)
    end do

    nshockpoints = 0; nshockpoints(1) = n
    nshockedges = 0; nshockedges(1) = n - 1

    call rd_dps_eq(xysh, zroeshu, zroeshd, 1_i4, nshockpoints, nshockedges)

    nout = nshockpoints(1)
    do i = 1, nout
      xysh_out(i) = xysh(1, i, 1)
      tracer_out(i) = zroeshd(1, i, 1)
    end do
  end subroutine run_rd_dps_eq

! 6 points already uniform at exactly dxcell: total length (n-1)*dxcell
! is an exact multiple of dxcell, so the new point count/spacing/abscissa
! all come out identical to the input -- a genuine no-op.
  subroutine test_idempotent(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(6), xysh_out(6), tracer_out(6)
    integer(i4) :: nout, i

    x = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
    call run_rd_dps_eq(x, 6_i4, 1.0_wp, xysh_out, tracer_out, nout)

    call check(error, nout, 6_i4)
    if (allocated(error)) return
    do i = 1, 6
      call check(error, xysh_out(i), x(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, tracer_out(i), 10.0_wp*real(i, wp), thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_idempotent

! 4 points, non-uniformly spaced (edges 1.0, 0.5, 2.5; total length 4.0)
! -> 5 new points at exactly dx=1.0 apart. Since points lie on the
! x-axis, the redistributed x-coordinates must be exactly [0,1,2,3,4] --
! the actual "uniform spacing restored" property this increment checks.
  subroutine test_redistribute(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(4), xysh_out(5), tracer_out(5)
    real(wp) :: expect_x(5), expect_tracer(5)
    integer(i4) :: nout, i

    x = [0.0_wp, 1.0_wp, 1.5_wp, 4.0_wp]
    call run_rd_dps_eq(x, 4_i4, 1.0_wp, xysh_out, tracer_out, nout)

    call check(error, nout, 5_i4, message="new point count")
    if (allocated(error)) return

    expect_x = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
    expect_tracer = [10.0_wp, 20.0_wp, 32.0_wp, 36.0_wp, 40.0_wp]
    do i = 1, 5
      call check(error, xysh_out(i), expect_x(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, tracer_out(i), expect_tracer(i), thr=1.0e-10_wp)
      if (allocated(error)) return
    end do

    ! spacing bound: every new edge is exactly dx=1.0 (uniform)
    do i = 2, 5
      call check(error, xysh_out(i) - xysh_out(i - 1), 1.0_wp, thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_redistribute

end module test_rd_dps_eq
