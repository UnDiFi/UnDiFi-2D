module test_rd_dps
! Phase 7.5 increment 7 (ROADMAP.md #20): rd_dps -- shock/discontinuity
! point redistribution (remove a too-short edge's point, insert a point
! at a too-long edge's midpoint) -- idempotency and spacing-bound checks.
!
! Per call, rd_dps finds AT MOST one too-short edge (length < 0.5*dxcell)
! and AT MOST one too-long edge (length > 1.5*dxcell) and fixes each
! independently -- not a full sweep to convergence, so "idempotency"
! here means: a curve with every edge already within [0.5,1.5]*dxcell is
! left completely unmodified (neither search fires). The insertion and
! removal index arithmetic is genuinely subtle (the insertion step reuses
! an array slot the shift loop *just* overwrote to compute the midpoint,
! rather than saving the old value first) -- verified against a standalone
! Python simulation of the exact same index arithmetic before writing the
! Fortran below, not just reasoned about on paper.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_rd_dps

contains

  subroutine collect_rd_dps(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("well_spaced_is_idempotent", test_idempotent),&
    &new_unittest("short_edge_removes_a_point", test_removal),&
    &new_unittest("long_edge_inserts_a_midpoint", test_insertion)&
    &]
  end subroutine collect_rd_dps

! Builds a single shock of n points along the x-axis at the given
! x-coordinates (y=0), with a two-component "tracer" in zroeshu/zroeshd
! (dof 1 = 100+original point index, dof 2 = 200+original point index,
! dof 3/4 left at 0) so post-call values can be traced back to exactly
! which original point (or average of which two) they came from.
  subroutine run_rd_dps(x, n, dxcell_in, xysh_out, tracer1_out, tracer2_out,&
  &nshockpoints_out, nshockedges_out)
    real(wp), intent(in) :: x(:)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: dxcell_in
    real(wp), intent(out) :: xysh_out(:)
    real(wp), intent(out) :: tracer1_out(:), tracer2_out(:)
    integer(i4), intent(out) :: nshockpoints_out, nshockedges_out

    include 'paramt.h'
    external :: rd_dps

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
      zroeshu(1, i, 1) = 100.0_wp + real(i, wp)
      zroeshu(2, i, 1) = 200.0_wp + real(i, wp)
      zroeshd(1, i, 1) = 100.0_wp + real(i, wp)
      zroeshd(2, i, 1) = 200.0_wp + real(i, wp)
    end do

    nshockpoints = 0; nshockpoints(1) = n
    nshockedges = 0; nshockedges(1) = n - 1

    call rd_dps(xysh, zroeshu, zroeshd, 1_i4, nshockpoints, nshockedges)

    nshockpoints_out = nshockpoints(1)
    nshockedges_out = nshockedges(1)
    do i = 1, nshockpoints_out
      xysh_out(i) = xysh(1, i, 1)
      tracer1_out(i) = zroeshd(1, i, 1)
      tracer2_out(i) = zroeshd(2, i, 1)
    end do
  end subroutine run_rd_dps

! Six points, every edge exactly dxcell (ratio 1.0, inside [0.5,1.5] both
! ways) -- neither the removal nor the insertion search should fire.
  subroutine test_idempotent(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(6), xysh_out(6), t1(6), t2(6)
    integer(i4) :: np, ne, i

    x = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
    call run_rd_dps(x, 6_i4, 1.0_wp, xysh_out, t1, t2, np, ne)

    call check(error, np, 6_i4)
    if (allocated(error)) return
    call check(error, ne, 5_i4)
    if (allocated(error)) return
    do i = 1, 6
      call check(error, xysh_out(i), x(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, t1(i), 100.0_wp + real(i, wp), thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_idempotent

! Edge 3 (points 3-4) is 0.2*dxcell, well under the 0.5 threshold; all
! other edges are exactly dxcell. Point 3 is the one that disappears
! (points 4,5,6 shift down to 3,4,5) -- verified against the standalone
! index-arithmetic simulation, not derived from the Fortran itself.
  subroutine test_removal(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(6), xysh_out(5), t1(5), t2(5)
    real(wp) :: expect_x(5), expect_t1(5), expect_t2(5)
    integer(i4) :: np, ne, i

    x = [0.0_wp, 1.0_wp, 2.0_wp, 2.2_wp, 3.2_wp, 4.2_wp]
    call run_rd_dps(x, 6_i4, 1.0_wp, xysh_out, t1, t2, np, ne)

    call check(error, np, 5_i4, message="point count after removal")
    if (allocated(error)) return
    call check(error, ne, 4_i4, message="edge count after removal")
    if (allocated(error)) return

    expect_x = [0.0_wp, 1.0_wp, 2.2_wp, 3.2_wp, 4.2_wp]
    expect_t1 = [101.0_wp, 102.0_wp, 104.0_wp, 105.0_wp, 106.0_wp]
    expect_t2 = [201.0_wp, 202.0_wp, 204.0_wp, 205.0_wp, 206.0_wp]
    do i = 1, 5
      call check(error, xysh_out(i), expect_x(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, t1(i), expect_t1(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, t2(i), expect_t2(i), thr=1.0e-12_wp)
      if (allocated(error)) return
    end do
  end subroutine test_removal

! Edge 3 (points 3-4) is 2.5*dxcell, well over the 1.5 threshold; all
! other edges are exactly dxcell. A new point is inserted between points
! 3 and 4 at their exact midpoint (position AND state, i.e. the
! interpolated tracer), and its two new edges (1.25*dxcell each) land
! back inside the "safe" [0.5,1.5] spacing band -- the actual "spacing
! bound restored" property this increment is meant to check.
  subroutine test_insertion(error)
    type(error_type), allocatable, intent(out) :: error
    real(wp) :: x(6), xysh_out(7), t1(7), t2(7)
    real(wp) :: expect_x(7), expect_t1(7), expect_t2(7)
    integer(i4) :: np, ne, i

    x = [0.0_wp, 1.0_wp, 2.0_wp, 4.5_wp, 5.5_wp, 6.5_wp]
    call run_rd_dps(x, 6_i4, 1.0_wp, xysh_out, t1, t2, np, ne)

    call check(error, np, 7_i4, message="point count after insertion")
    if (allocated(error)) return
    call check(error, ne, 6_i4, message="edge count after insertion")
    if (allocated(error)) return

    expect_x = [0.0_wp, 1.0_wp, 2.0_wp, 3.25_wp, 4.5_wp, 5.5_wp, 6.5_wp]
    expect_t1 = [101.0_wp, 102.0_wp, 103.0_wp, 103.5_wp, 104.0_wp, 105.0_wp, 106.0_wp]
    expect_t2 = [201.0_wp, 202.0_wp, 203.0_wp, 203.5_wp, 204.0_wp, 205.0_wp, 206.0_wp]
    do i = 1, 7
      call check(error, xysh_out(i), expect_x(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, t1(i), expect_t1(i), thr=1.0e-12_wp)
      if (allocated(error)) return
      call check(error, t2(i), expect_t2(i), thr=1.0e-12_wp)
      if (allocated(error)) return
    end do

    ! spacing bound: the two edges either side of the new point (index 4)
    ! are each 1.25*dxcell, back inside [0.5,1.5]
    call check(error, abs(xysh_out(4) - xysh_out(3)), 1.25_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
    call check(error, abs(xysh_out(5) - xysh_out(4)), 1.25_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
  end subroutine test_insertion

end module test_rd_dps
