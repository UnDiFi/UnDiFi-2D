subroutine readmesh(lbndfac, lcelcel, lcelnod, lcorg, ledgptr,&
&lnodcod, lnodptr, lzroe, nbfac, npoin, nelem, nhole, nbpoin,&
&nvt, nedge, fname, lia, lja, liclr, nclr, fndbnds)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, neshmax, npshmax, nshmax, zero
  implicit none(type, external)
  external check, rtri, setbndrynodeptr

!     this subroutine reads a mesh written by the triangle code.
!     the mesh itself is built locally into plain allocatable arrays
!     (no istak/dstak/istkgt); the result is then copied into the
!     shared dstak arena so that the pointer-offset outputs below stay
!     valid for callers that have not yet been converted to mod_mesh
!     (ROADMAP.md Phase 2.2/2.3, issue #13) -- this bridging copy is
!     scaffolding, removed once main.f90 itself holds a type(mesh_t).
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

!     .. pointer-offset dummy arguments (bridged into dstak at the end) ..
  integer(i4) lbndfac, lcelcel, lcelnod, lcorg, ledgptr, lnodcod, lnodptr,&
  &lzroe, lia, lja, liclr, nclr
  integer(i4) nbfac, npoin, nelem, nhole, nbpoin, nvt, nedge
  integer(i4) k
  logical fndbnds

  integer(i4) lenstr
  external lenstr

  character*(*) fname

!     .. arrays in common (bridge target only) ..
  real(wp) dstak(1)
  integer(i4) istak(1)
  common/cstak/dstak
  equivalence(dstak(1), istak(1))
  integer(i4) istkgt
  external istkgt

!     .. local mesh arrays (the real work happens on these) ..
  real(wp), allocatable :: corg(:, :), zroe(:, :)
  integer(i4), allocatable :: nodcod(:), celnod(:, :), celcel(:, :),&
  &edgptr(:, :), bndfac(:, :)
  integer(i4), allocatable :: nodptr(:, :), ia(:), ja(:), iclr(:)
  integer(i4) npad, nbfacpad

!     open log file
  open (8, file='log/readmesh.log')

!     ********** node file ***************
!     mode 0: just learn npoin (rtri never touches the array dummies here)
  call rtri(istak(1), nedge, istak(1), nbfac, istak(1)&
  &, istak(1), nvt, nelem, dstak(1), dstak(1),&
  &istak(1), npoin, "node", 0, fname)

!     allocate nodal and shock coordinates  corg(1:ndim,1:npoin+nshmax*npshmax)
  npad = npoin + 2*nshmax*npshmax
  allocate (corg(ndim, npad), source=zero)
  allocate (nodcod(npad), source=0_i4)
  allocate (zroe(ndof, npad), source=zero)

  call rtri(istak(1), nedge, istak(1), nbfac, istak(1)&
  &, istak(1), nvt, nelem, corg, zroe,&
  &nodcod, npoin, "node", 1, fname)

  write (8, *) 'rtri node  --> ok'

!     ********** ele file ***************
  call rtri(istak(1), nedge, istak(1), nbfac, istak(1)&
  &, istak(1), nvt, nelem, corg, zroe,&
  &nodcod, npoin, "ele", 0, fname)

!     allocate cell to node pointer, then read icelnod(1:nvt,1:nelem)
!     gives the nvt vertices of a given cell
  allocate (celnod(nvt, nelem), source=0_i4)

  call rtri(istak(1), nedge, istak(1), nbfac, celnod&
  &, istak(1), nvt, nelem, corg, zroe,&
  &nodcod, npoin, "ele", 1, fname)

  write (8, *) 'rtri ele  --> ok'

!     ********** neigh file ***************
  call rtri(istak(1), nedge, istak(1), nbfac, celnod&
  &, istak(1), nvt, nelem, corg, zroe,&
  &nodcod, npoin, "neigh", 0, fname)

!     allocate cell to cell pointer, then read icelcel(1:nvt,1:nelem)
!     gives the nvt neighbouring triangles of a given cell
!     icelcel(i,ielem) is the triangle sharing the edge opposite
!     the i-th vertex of triangle ielem
  allocate (celcel(nvt, nelem), source=0_i4)

  call rtri(istak(1), nedge, istak(1), nbfac, celnod&
  &, celcel, nvt, nelem, corg, zroe,&
  &nodcod, npoin, "neigh", 1, fname)

  write (8, *) 'rtri neigh --> ok'

!     ********** edge file ***************
  call rtri(istak(1), nedge, istak(1), nbfac, celnod&
  &, celcel, nvt, nelem, corg, zroe,&
  &nodcod, npoin, "edge", 0, fname)

!     allocate edge to vertex pointer, then read
  allocate (edgptr(3, nedge), source=0_i4)

  call rtri(edgptr, nedge, istak(1), nbfac, celnod&
  &, celcel, nvt, nelem, corg, zroe,&
  &nodcod, npoin, "edge", 1, fname)

  write (8, *) 'rtri edge  --> ok'

!     ********** poly file ***************
  call rtri(edgptr, nedge, istak(1), nbfac, celnod&
  &, celcel, nvt, nelem, corg, zroe,&
  &nodcod, npoin, "poly", 0, fname)

!     boundary structure; ibndfac(1:3,1:nbfac+nshmax*neshmax)
!     ibndfac(1,ibfac) and ibndfac(2,ibfac) are the endpoints
!     of bndry segment ibfac; ibndfac(3,ibfac) is the colour
!     ibndfac stores not only boundary edges, but also shock edges
!     shock edges are coloured 10
  nbfacpad = nbfac + 2*nshmax*neshmax ! leave room for duplicated nodes
  allocate (bndfac(3, nbfacpad), source=0_i4)

  call rtri(edgptr, nedge, bndfac, nbfac, celnod&
  &, celcel, nvt, nelem, corg, zroe,&
  &nodcod, npoin, "poly", 1, fname)

  k = lenstr(fname)
  write (8, *) 'done reading ', fname(1:k), ' triangle files'

!     nodptr,ia,ja,iclr are computed by setbndrynodeptr; the number of
!     colours within the mesh is nclr
  if (fndbnds) then
    call setbndrynodeptr(bndfac, nodcod, nbfac, npoin,&
    &nbpoin, nodptr, ia, ja, iclr, nclr)

!     the following subroutine writes the sequence of gridpoints of each
!     boundary (patch) separately within a different file bndry??.dat
!     it also illustrates how to use the various pointers
    call check(ia, ja, iclr, nclr, corg, ndim)
  end if

!     .. bridge: copy the locally-built mesh into the shared dstak arena
!        so that pointer-offset outputs stay valid for unconverted callers ..
  lcorg = istkgt(ndim*npad, 4)
  dstak(lcorg:lcorg + ndim*npad - 1) = reshape(corg, [ndim*npad])

  lzroe = istkgt(ndof*npad, 4)
  dstak(lzroe:lzroe + ndof*npad - 1) = reshape(zroe, [ndof*npad])

  lnodcod = istkgt(npad, 2)
  istak(lnodcod:lnodcod + npad - 1) = nodcod

  lcelnod = istkgt(nvt*nelem, 2)
  istak(lcelnod:lcelnod + nvt*nelem - 1) = reshape(celnod, [nvt*nelem])

  lcelcel = istkgt(nvt*nelem, 2)
  istak(lcelcel:lcelcel + nvt*nelem - 1) = reshape(celcel, [nvt*nelem])

  ledgptr = istkgt(3*nedge, 2)
  istak(ledgptr:ledgptr + 3*nedge - 1) = reshape(edgptr, [3*nedge])

  lbndfac = istkgt(3*nbfacpad, 2)
  istak(lbndfac:lbndfac + 3*nbfacpad - 1) = reshape(bndfac, [3*nbfacpad])

  if (fndbnds) then
    lnodptr = istkgt(3*nbpoin, 2)
    istak(lnodptr:lnodptr + 3*nbpoin - 1) = reshape(nodptr, [3*nbpoin])

    lia = istkgt(nclr + 1, 2)
    istak(lia:lia + nclr) = ia

    liclr = istkgt(nclr, 2)
    istak(liclr:liclr + nclr - 1) = iclr

    lja = istkgt(ia(nclr + 1) - ia(1), 2)
    istak(lja:lja + (ia(nclr + 1) - ia(1)) - 1) = ja
  end if

  return
end subroutine readmesh
