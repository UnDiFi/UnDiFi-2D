module test_interp
! Phase 7.5 increment 9 (ROADMAP.md #20): interp -- background-mesh
! phantom-node interpolation -- accuracy/conservation on a manufactured
! field.
!
! interp.f90's own per-phantom-node loop calls finder_search (the
! candidate-list-restricted twin of finder, both defined in interp.f90),
! which locates the enclosing triangle via barycentric coordinates and
! interpolates zroe linearly over that triangle -- this is the actual
! interpolation kernel; interp itself is mostly orchestration (bin-grid
! spatial search, boundary-edge fallback, OpenMP loop structure) around
! it, needing a full background mesh + boundary-face array to exercise
! end-to-end. Tested here directly against finder_search on a small
! 2-triangle manufactured mesh (unit square split along its diagonal),
! with two manufactured affine fields (z=a+b*x+c*y per dof, different
! coefficients per dof) at the 4 nodes.
!
! Barycentric interpolation is exact for any affine function -- not just
! asymptotically accurate -- so this checks exact reproduction at
! interior points, not a convergence rate. "Conservation" is checked as
! the affine field's degenerate case: a spatially uniform field is
! conserved exactly at every interior point. A point outside every
! triangle is checked to correctly fail the search (info /= 0) rather
! than silently returning a stale/garbage interpolant.
  use mod_kinds, only: wp, i4
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_interp

contains

  subroutine collect_interp(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("affine_field_exact_in_each_triangle", test_affine),&
    &new_unittest("uniform_field_conserved", test_uniform),&
    &new_unittest("point_outside_mesh_fails_search", test_outside)&
    &]
  end subroutine collect_interp

! Unit square (0,0)-(1,0)-(1,1)-(0,1), split along the (0,0)-(1,1)
! diagonal into T1=(1,2,3) (lower-right) and T2=(1,3,4) (upper-left),
! both nodes listed counterclockwise (AREA, lib/libmylib/area.f, is
! positive for CCW -- required for finder_search's barycentric ratios to
! land in [0,1] for an interior point).
  subroutine build_mesh(coor, icelnod, zroe)
    real(wp), intent(out) :: coor(2, 4)
    integer(i4), intent(out) :: icelnod(3, 2)
    real(wp), intent(out) :: zroe(2, 4)
    integer(i4) :: i
    real(wp) :: x, y

    coor(:, 1) = [0.0_wp, 0.0_wp]
    coor(:, 2) = [1.0_wp, 0.0_wp]
    coor(:, 3) = [1.0_wp, 1.0_wp]
    coor(:, 4) = [0.0_wp, 1.0_wp]

    icelnod(:, 1) = [1, 2, 3]
    icelnod(:, 2) = [1, 3, 4]

    ! dof 1: z = 2.0 + 3.0*x - 1.5*y ; dof 2: z = 1.0 - 0.5*x + 2.0*y
    do i = 1, 4
      x = coor(1, i); y = coor(2, i)
      zroe(1, i) = 2.0_wp + 3.0_wp*x - 1.5_wp*y
      zroe(2, i) = 1.0_wp - 0.5_wp*x + 2.0_wp*y
    end do
  end subroutine build_mesh

  subroutine test_affine(error)
    type(error_type), allocatable, intent(out) :: error

    external :: finder_search
    real(wp) :: coor(2, 4), zroe(2, 4), zout(2), xyin(2)
    integer(i4) :: icelnod(3, 2), candidates(2), ielem, info

    call build_mesh(coor, icelnod, zroe)
    candidates = [1, 2]

    ! P=(0.3,0.2): below the diagonal (y<x) -> triangle T1
    xyin = [0.3_wp, 0.2_wp]
    call finder_search(icelnod, coor, 2_i4, zroe, 2_i4, xyin, zout,&
    &candidates, 2_i4, ielem, info)
    call check(error, info, 0_i4, message="P search failed")
    if (allocated(error)) return
    call check(error, ielem, 1_i4, message="P not found in T1")
    if (allocated(error)) return
    call check(error, zout(1), 2.0_wp + 3.0_wp*0.3_wp - 1.5_wp*0.2_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
    call check(error, zout(2), 1.0_wp - 0.5_wp*0.3_wp + 2.0_wp*0.2_wp, thr=1.0e-12_wp)
    if (allocated(error)) return

    ! Q=(0.2,0.7): above the diagonal (y>x) -> triangle T2
    xyin = [0.2_wp, 0.7_wp]
    call finder_search(icelnod, coor, 2_i4, zroe, 2_i4, xyin, zout,&
    &candidates, 2_i4, ielem, info)
    call check(error, info, 0_i4, message="Q search failed")
    if (allocated(error)) return
    call check(error, ielem, 2_i4, message="Q not found in T2")
    if (allocated(error)) return
    call check(error, zout(1), 2.0_wp + 3.0_wp*0.2_wp - 1.5_wp*0.7_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
    call check(error, zout(2), 1.0_wp - 0.5_wp*0.2_wp + 2.0_wp*0.7_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
  end subroutine test_affine

  subroutine test_uniform(error)
    type(error_type), allocatable, intent(out) :: error

    external :: finder_search
    real(wp) :: coor(2, 4), zroe_affine(2, 4), zroe(2, 4), zout(2), xyin(2)
    integer(i4) :: icelnod(3, 2), candidates(2), ielem, info

    call build_mesh(coor, icelnod, zroe_affine)
    zroe = 7.0_wp   ! spatially uniform field, both dof
    candidates = [1, 2]

    xyin = [0.45_wp, 0.55_wp]  ! interior to T2, arbitrary
    call finder_search(icelnod, coor, 2_i4, zroe, 2_i4, xyin, zout,&
    &candidates, 2_i4, ielem, info)
    call check(error, info, 0_i4)
    if (allocated(error)) return
    call check(error, zout(1), 7.0_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
    call check(error, zout(2), 7.0_wp, thr=1.0e-12_wp)
    if (allocated(error)) return
  end subroutine test_uniform

  subroutine test_outside(error)
    type(error_type), allocatable, intent(out) :: error

    external :: finder_search
    real(wp) :: coor(2, 4), zroe(2, 4), zout(2), xyin(2)
    integer(i4) :: icelnod(3, 2), candidates(2), ielem, info

    call build_mesh(coor, icelnod, zroe)
    candidates = [1, 2]

    xyin = [5.0_wp, 5.0_wp]  ! well outside the unit square
    call finder_search(icelnod, coor, 2_i4, zroe, 2_i4, xyin, zout,&
    &candidates, 2_i4, ielem, info)
    call check(error, info, 1_i4, message="search should fail outside the mesh")
    if (allocated(error)) return
  end subroutine test_outside

end module test_interp
