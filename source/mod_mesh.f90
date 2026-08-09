module mod_mesh
! Phase 2.2 (ROADMAP.md #13): the allocatable replacement for the
! dstak/istak pointer-offset arrays readmesh.f90/main.f90 use today
! (lcorg, lzroe, lcelnod, lcelcel, lbndfac, lnodcod, lnodptr, ledgptr,
! lia, lja, liclr, lpmap). Not yet wired into readmesh.f90/main.f90 --
! that conversion is later Phase 2 sub-items.

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: mesh_t

  type :: mesh_t
    integer(i4)              :: npoin = 0, nelem = 0, nbfac = 0, nbpoin = 0
    integer(i4)              :: nedge = 0, nhole = 0, nvt = 3, nclr = 0, npnod = 0
    real(wp),    allocatable :: xy(:, :)     ! (ndim, npoin)
    real(wp),    allocatable :: zroe(:, :)   ! (ndof, npoin)
    integer(i4), allocatable :: celnod(:, :) ! (nvt, nelem)
    integer(i4), allocatable :: celcel(:, :) ! (nvt, nelem)
    integer(i4), allocatable :: bndfac(:, :) ! (3, nbfac)
    integer(i4), allocatable :: nodcod(:)    ! (npoin)
    integer(i4), allocatable :: nodptr(:, :) ! (nbpoin,3), from setbndrynodeptr
    integer(i4), allocatable :: edgptr(:, :) ! (3, nedge)
    integer(i4), allocatable :: ia(:), ja(:), iclr(:) ! boundary-patch CSR, from setbndrynodeptr
    integer(i4), allocatable :: pmap(:)      ! (npoin), periodic-node map, see readpmap.f90
  contains
    procedure :: alloc => mesh_alloc
    procedure :: free  => mesh_free
    final     :: mesh_finalize
  end type mesh_t

  ! No copy_from: plain intrinsic assignment (bak = bkg) already deep-copies
  ! every allocatable component -- a wrapper would need a non-polymorphic
  ! passed-object dummy, which a type-bound procedure cannot have.

contains

  subroutine mesh_alloc(this, ndim, ndof, npoin, nelem, nbfac, nvt, nedge)
    class(mesh_t), intent(inout) :: this
    integer(i4),   intent(in)    :: ndim, ndof, npoin, nelem, nbfac, nvt, nedge

    call this%free()

    this%npoin = npoin
    this%nelem = nelem
    this%nbfac = nbfac
    this%nvt   = nvt
    this%nedge = nedge

    allocate (this%xy(ndim, npoin),    source=0.0_wp)
    allocate (this%zroe(ndof, npoin),  source=0.0_wp)
    allocate (this%celnod(nvt, nelem), source=0_i4)
    allocate (this%celcel(nvt, nelem), source=0_i4)
    allocate (this%bndfac(3, nbfac),   source=0_i4)
    allocate (this%nodcod(npoin),      source=0_i4)
    allocate (this%edgptr(3, nedge),   source=0_i4)
  end subroutine mesh_alloc

  subroutine mesh_free(this)
    class(mesh_t), intent(inout) :: this

    this%npoin = 0; this%nelem = 0; this%nbfac = 0; this%nbpoin = 0
    this%nedge = 0; this%nhole = 0; this%nvt = 3; this%nclr = 0; this%npnod = 0

    if (allocated(this%xy))     deallocate (this%xy)
    if (allocated(this%zroe))   deallocate (this%zroe)
    if (allocated(this%celnod)) deallocate (this%celnod)
    if (allocated(this%celcel)) deallocate (this%celcel)
    if (allocated(this%bndfac)) deallocate (this%bndfac)
    if (allocated(this%nodcod)) deallocate (this%nodcod)
    if (allocated(this%nodptr)) deallocate (this%nodptr)
    if (allocated(this%edgptr)) deallocate (this%edgptr)
    if (allocated(this%ia))     deallocate (this%ia)
    if (allocated(this%ja))     deallocate (this%ja)
    if (allocated(this%iclr))   deallocate (this%iclr)
    if (allocated(this%pmap))   deallocate (this%pmap)
  end subroutine mesh_free

  subroutine mesh_finalize(this)
    type(mesh_t), intent(inout) :: this

    call this%free()
  end subroutine mesh_finalize

end module mod_mesh
