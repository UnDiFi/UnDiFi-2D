subroutine readmesh(mesh, fname, fndbnds)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, neshmax, npshmax, nshmax, zero
  use mod_mesh, only: mesh_t
  implicit none(type, external)
  external check, rtri, setbndrynodeptr

!     this subroutine reads a mesh written by the triangle code, straight
!     into mesh's allocatable arrays (ROADMAP.md Phase 2.3/2.4, issue #13).
!
!     there are various triangle file that can be read:
!     node neigh edge poly
!     their format is explained in triangle/triangle.ps
!
!     ndim  is the space dimension = 2
!     nvt = ndim+1 is the number of vertices
!
!     nshmax  : max number of shock
!     npshmax : max number of shock points for each shock.
!     neshmax : max number of shock element for each shocks

  type(mesh_t), intent(inout) :: mesh
  character*(*) fname
  logical fndbnds

  integer(i4) k, nvt, npad, nbfacpad

  integer(i4) lenstr
  external lenstr

!     placeholders for the mode-0 rtri calls, which never touch their
!     array dummies (only the scalar counts are learned in that mode)
  real(wp) dummyr(1)
  integer(i4) dummyi(1)

!     open log file
  open (8, file='log/readmesh.log')

  call mesh%free()
  nvt = 3

!     ********** node file ***************
!     mode 0: just learn npoin
  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, dummyi&
  &, dummyi, nvt, mesh%nelem, dummyr, dummyr,&
  &dummyi, mesh%npoin, "node", 0, fname)

!     allocate nodal and shock coordinates  xy(1:ndim,1:npoin+nshmax*npshmax)
  npad = mesh%npoin + 2*nshmax*npshmax
  allocate (mesh%xy(ndim, npad), source=zero)
  allocate (mesh%nodcod(npad), source=0_i4)
  allocate (mesh%zroe(ndof, npad), source=zero)

  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, dummyi&
  &, dummyi, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "node", 1, fname)

  write (8, *) 'rtri node  --> ok'

!     ********** ele file ***************
  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, dummyi&
  &, dummyi, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "ele", 0, fname)

!     allocate cell to node pointer, then read celnod(1:nvt,1:nelem)
!     gives the nvt vertices of a given cell
  allocate (mesh%celnod(nvt, mesh%nelem), source=0_i4)

  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, dummyi, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "ele", 1, fname)

  write (8, *) 'rtri ele  --> ok'

!     ********** neigh file ***************
  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, dummyi, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "neigh", 0, fname)

!     allocate cell to cell pointer, then read celcel(1:nvt,1:nelem)
!     gives the nvt neighbouring triangles of a given cell
!     celcel(i,ielem) is the triangle sharing the edge opposite
!     the i-th vertex of triangle ielem
  allocate (mesh%celcel(nvt, mesh%nelem), source=0_i4)

  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, mesh%celcel, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "neigh", 1, fname)

  write (8, *) 'rtri neigh --> ok'

!     ********** edge file ***************
  call rtri(dummyi, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, mesh%celcel, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "edge", 0, fname)

!     allocate edge to vertex pointer, then read
  allocate (mesh%edgptr(3, mesh%nedge), source=0_i4)

  call rtri(mesh%edgptr, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, mesh%celcel, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "edge", 1, fname)

  write (8, *) 'rtri edge  --> ok'

!     ********** poly file ***************
  call rtri(mesh%edgptr, mesh%nedge, dummyi, mesh%nbfac, mesh%celnod&
  &, mesh%celcel, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "poly", 0, fname)

!     boundary structure; bndfac(1:3,1:nbfac+nshmax*neshmax)
!     bndfac(1,ibfac) and bndfac(2,ibfac) are the endpoints
!     of bndry segment ibfac; bndfac(3,ibfac) is the colour
!     bndfac stores not only boundary edges, but also shock edges
!     shock edges are coloured 10
  nbfacpad = mesh%nbfac + 2*nshmax*neshmax ! leave room for duplicated nodes
  allocate (mesh%bndfac(3, nbfacpad), source=0_i4)

  call rtri(mesh%edgptr, mesh%nedge, mesh%bndfac, mesh%nbfac, mesh%celnod&
  &, mesh%celcel, nvt, mesh%nelem, mesh%xy, mesh%zroe,&
  &mesh%nodcod, mesh%npoin, "poly", 1, fname)

  mesh%nvt = nvt
  k = lenstr(fname)
  write (8, *) 'done reading ', fname(1:k), ' triangle files'

!     nodptr,ia,ja,iclr are computed by setbndrynodeptr; the number of
!     colours within the mesh is nclr
  if (fndbnds) then
    call setbndrynodeptr(mesh%bndfac, mesh%nodcod, mesh%nbfac, mesh%npoin,&
    &mesh%nbpoin, mesh%nodptr, mesh%ia, mesh%ja, mesh%iclr, mesh%nclr)

!     the following subroutine writes the sequence of gridpoints of each
!     boundary (patch) separately within a different file bndry??.dat
!     it also illustrates how to use the various pointers
    call check(mesh%ia, mesh%ja, mesh%iclr, mesh%nclr, mesh%xy, ndim)
  end if

  return
end subroutine readmesh
