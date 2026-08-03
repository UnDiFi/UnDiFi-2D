!     after readmesh():
!
!     call readpmap(mesh)
!
!     for a given meshpoint i
!     j = mesh%pmap(i) /= 0 is the periodic meshpoint corresponding to i
!     j = mesh%pmap(i) == 0 means that i is not among the periodic meshpoints
!     mesh%npnod is the nof periodic gridpoints in only one of the two sets

subroutine readpmap(mesh)
  use mod_kinds, only: wp, i4
  use mod_mesh, only: mesh_t
  implicit none(type, external)

  type(mesh_t), intent(inout) :: mesh

  logical lflag
  integer(i4) i, j, n

!     read the file with the index of all nodes and in case
!     the corresponding node on the periodic boundary
  inquire (file="pnodes0.dat", exist=lflag)
  if (.not. lflag) then
!     if the file is absent, build a table with all zero (i.e. with any periodic node)
    allocate (mesh%pmap(mesh%npoin), source=0_i4)
    return
  end if
  open (13, file="pnodes0.dat")
  read (13, *) n, mesh%npnod
  if (n .ne. mesh%npoin) then
    write (6, *) 'the nof meshpoints in the dataset and in pnodes0.d&
    &at do not match'
    error stop 13
  end if
  allocate (mesh%pmap(n), source=0_i4)
  do i = 1, n
    read (13, *) j, mesh%pmap(i)
  end do
  close (13)
  return
end subroutine readpmap
