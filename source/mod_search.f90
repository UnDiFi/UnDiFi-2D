module mod_search
! Phase 4.1 (ROADMAP.md #16): a uniform bin grid replacing the brute-force
! O(N_elem) scans in fnd_phps.f90 (shock-segment/cell crossing) and
! interp.f90's finder (phantom-point/cell containment) -- finding F5.
! Not a k-d-tree: mesh cell sizes are roughly uniform within one
! DXCELL-controlled refinement level, so a uniform grid is a simpler,
! easier-to-verify fit; a k-d-tree could replace this behind the same
! query_point/query_segment interface later if a non-uniform case needs
! it. This module only NARROWS the candidate element set -- every
! existing geometric decision (ishel1/ishel2/rdshp in fnd_phps.f90, the
! barycentric test in interp.f90) is untouched by design, so results are
! unchanged as long as query_point/query_segment never omit a true hit
! (false positives are just extra, cheap re-checks of the same geometry
! test that used to run on every element anyway).
!
! Usage: build_bin_grid once per caller invocation (the mesh changes
! every outer iteration, so there's no cross-call caching), then
! query_point/query_segment as many times as needed against that one
! grid. bin_grid_t's allocatable components are freed automatically when
! a local `type(bin_grid_t)` variable goes out of scope -- callers don't
! need to explicitly tear it down.

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private

  public :: bin_grid_t, build_bin_grid, query_point, query_segment

  type :: bin_grid_t
    integer(i4) :: nx = 0, ny = 0
    real(wp) :: xmin = 0.0_wp, ymin = 0.0_wp
    real(wp) :: hx = 1.0_wp, hy = 1.0_wp
    ! CSR layout: bin b's elements are bin_elems(bin_start(b):bin_start(b+1)-1)
    integer(i4), allocatable :: bin_start(:)
    integer(i4), allocatable :: bin_elems(:)
  end type bin_grid_t

contains

  ! icelnod(nvt,nelem)/xy(2,*): same arrays fnd_phps.f90/interp.f90 already
  ! use for their own scans -- the bounding box is computed purely from
  ! the nodes each element actually references, so no separate total-
  ! point-count argument is needed (interp.f90 in particular has no clean
  ! single "total points in this mesh" value to hand in).
  subroutine build_bin_grid(grid, icelnod, nvt, xy, nelem)
    type(bin_grid_t), intent(out) :: grid
    integer(i4), intent(in) :: nvt, nelem
    integer(i4), intent(in) :: icelnod(nvt, nelem)
    real(wp), intent(in) :: xy(2, *)

    integer(i4) :: ielem, ix0, ix1, iy0, iy1, ix, iy, ibin, nbins
    real(wp) :: xmin, xmax, ymin, ymax, exlo, exhi, eylo, eyhi
    integer(i4), allocatable :: counts(:), cursor(:)

    xmin = huge(1.0_wp); xmax = -huge(1.0_wp)
    ymin = huge(1.0_wp); ymax = -huge(1.0_wp)
    do ielem = 1, nelem
      call element_bbox(icelnod, nvt, xy, ielem, exlo, exhi, eylo, eyhi)
      xmin = min(xmin, exlo); xmax = max(xmax, exhi)
      ymin = min(ymin, eylo); ymax = max(ymax, eyhi)
    end do

    grid%nx = max(1, ceiling(sqrt(real(nelem, wp))))
    grid%ny = grid%nx
    grid%xmin = xmin
    grid%ymin = ymin
    grid%hx = merge((xmax - xmin)/real(grid%nx, wp), 1.0_wp, xmax > xmin)
    grid%hy = merge((ymax - ymin)/real(grid%ny, wp), 1.0_wp, ymax > ymin)

    nbins = grid%nx*grid%ny
    allocate (counts(nbins), source=0_i4)

    ! pass 1: count each bin's element membership (an element spanning
    ! multiple bins is counted once per bin it overlaps)
    do ielem = 1, nelem
      call element_bbox(icelnod, nvt, xy, ielem, exlo, exhi, eylo, eyhi)
      ix0 = bin_coord(exlo, grid%xmin, grid%hx, grid%nx)
      ix1 = bin_coord(exhi, grid%xmin, grid%hx, grid%nx)
      iy0 = bin_coord(eylo, grid%ymin, grid%hy, grid%ny)
      iy1 = bin_coord(eyhi, grid%ymin, grid%hy, grid%ny)
      do iy = iy0, iy1
        do ix = ix0, ix1
          ibin = (iy - 1)*grid%nx + ix
          counts(ibin) = counts(ibin) + 1
        end do
      end do
    end do

    allocate (grid%bin_start(nbins + 1))
    grid%bin_start(1) = 1
    do ibin = 1, nbins
      grid%bin_start(ibin + 1) = grid%bin_start(ibin) + counts(ibin)
    end do

    allocate (grid%bin_elems(grid%bin_start(nbins + 1) - 1))
    allocate (cursor(nbins))
    cursor = grid%bin_start(1:nbins)

    ! pass 2: scatter each element's index into every bin it overlaps
    do ielem = 1, nelem
      call element_bbox(icelnod, nvt, xy, ielem, exlo, exhi, eylo, eyhi)
      ix0 = bin_coord(exlo, grid%xmin, grid%hx, grid%nx)
      ix1 = bin_coord(exhi, grid%xmin, grid%hx, grid%nx)
      iy0 = bin_coord(eylo, grid%ymin, grid%hy, grid%ny)
      iy1 = bin_coord(eyhi, grid%ymin, grid%hy, grid%ny)
      do iy = iy0, iy1
        do ix = ix0, ix1
          ibin = (iy - 1)*grid%nx + ix
          grid%bin_elems(cursor(ibin)) = ielem
          cursor(ibin) = cursor(ibin) + 1
        end do
      end do
    end do
  end subroutine build_bin_grid

  ! Candidates for "which element contains (x,y)": the point's own home
  ! bin's list is sufficient -- if an element contains the point, that
  ! element's bounding box necessarily overlaps the point's bin, so
  ! build_bin_grid already put it there. No neighbor-bin fallback needed.
  subroutine query_point(grid, x, y, candidates)
    type(bin_grid_t), intent(in) :: grid
    real(wp), intent(in) :: x, y
    integer(i4), allocatable, intent(out) :: candidates(:)

    integer(i4) :: ix, iy, ibin

    ix = bin_coord(x, grid%xmin, grid%hx, grid%nx)
    iy = bin_coord(y, grid%ymin, grid%hy, grid%ny)
    ibin = (iy - 1)*grid%nx + ix
    candidates = grid%bin_elems(grid%bin_start(ibin):grid%bin_start(ibin + 1) - 1)
  end subroutine query_point

  ! Candidates for "which elements might be crossed by segment
  ! (x1,y1)-(x2,y2)": the union of every bin the segment's bounding box
  ! overlaps. May contain duplicates when an element spans more than one
  ! of those bins -- harmless, since the geometric tests applied to each
  ! candidate afterward (ishel1/ishel2/rdshp) are idempotent.
  subroutine query_segment(grid, x1, y1, x2, y2, candidates)
    type(bin_grid_t), intent(in) :: grid
    real(wp), intent(in) :: x1, y1, x2, y2
    integer(i4), allocatable, intent(out) :: candidates(:)

    integer(i4) :: ix0, ix1, iy0, iy1, ix, iy, ibin, n, pos

    ix0 = bin_coord(min(x1, x2), grid%xmin, grid%hx, grid%nx)
    ix1 = bin_coord(max(x1, x2), grid%xmin, grid%hx, grid%nx)
    iy0 = bin_coord(min(y1, y2), grid%ymin, grid%hy, grid%ny)
    iy1 = bin_coord(max(y1, y2), grid%ymin, grid%hy, grid%ny)

    n = 0
    do iy = iy0, iy1
      do ix = ix0, ix1
        ibin = (iy - 1)*grid%nx + ix
        n = n + grid%bin_start(ibin + 1) - grid%bin_start(ibin)
      end do
    end do

    allocate (candidates(n))
    pos = 0
    do iy = iy0, iy1
      do ix = ix0, ix1
        ibin = (iy - 1)*grid%nx + ix
        n = grid%bin_start(ibin + 1) - grid%bin_start(ibin)
        if (n > 0) then
          candidates(pos + 1:pos + n) = grid%bin_elems(grid%bin_start(ibin):grid%bin_start(ibin + 1) - 1)
          pos = pos + n
        end if
      end do
    end do
  end subroutine query_segment

  subroutine element_bbox(icelnod, nvt, xy, ielem, xlo, xhi, ylo, yhi)
    integer(i4), intent(in) :: nvt, ielem
    integer(i4), intent(in) :: icelnod(nvt, *)
    real(wp), intent(in) :: xy(2, *)
    real(wp), intent(out) :: xlo, xhi, ylo, yhi

    integer(i4) :: iv, ipoin
    real(wp) :: x, y

    xlo = huge(1.0_wp); xhi = -huge(1.0_wp)
    ylo = huge(1.0_wp); yhi = -huge(1.0_wp)
    do iv = 1, nvt
      ipoin = icelnod(iv, ielem)
      x = xy(1, ipoin); y = xy(2, ipoin)
      xlo = min(xlo, x); xhi = max(xhi, x)
      ylo = min(ylo, y); yhi = max(yhi, y)
    end do
  end subroutine element_bbox

  ! Maps a coordinate to a bin index in [1,n], clamping so that a
  ! coordinate exactly at the bounding box's max edge (v=vmin+n*h) lands
  ! in bin n rather than one past the end.
  integer(i4) function bin_coord(v, vmin, h, n) result(idx)
    real(wp), intent(in) :: v, vmin, h
    integer(i4), intent(in) :: n

    idx = int((v - vmin)/h) + 1
    if (idx < 1) idx = 1
    if (idx > n) idx = n
  end function bin_coord

end module mod_search
