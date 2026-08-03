!> @bndfac[in] bndry face pointer read from the poly file, bndfac(3,nbfac)
!> @nodcod[in] nodal pointer read from the node or poly file, nodcod(npoin)
!> @nbfac[in] nof bndry faces
!> @npoin[in] nof vertices
!> @nbpoin[out] nof boundary vertices
!> @nodptr[out] pointer for boundary vertices, nodptr(nbpoin,3)
!> @ia[out] ia(1:nclr+1)
!> @ja[out] ja(1:nnzr)
!> @iclr[out] iclr(1:nclr) colour of the nclr patches
subroutine setbndrynodeptr(bndfac, nodcod, nbfac, npoin,&
&nbpoin, nodptr, ia, ja, iclr, nclr)
!
  use mod_kinds, only: wp, i4
  implicit none(type, external)
!
!     ndim  is the space dimension =2
!     nvt = ndim+1 is the number of vertices
!
!     npshmax     : max number of shock points.
!     neshmax     : max number of shock elements
!
!     .. scalar arguments ..
  integer(i4), intent(in) :: nbfac, npoin
  integer(i4), intent(out) :: nbpoin, nclr

!     .. array arguments ..
  integer(i4), intent(in) :: bndfac(3, nbfac), nodcod(npoin)
  integer(i4), allocatable, intent(out) :: nodptr(:, :), ia(:), ja(:), iclr(:)

!     .. local scalars ..
  integer(i4) ipoin, k, nnzr

!     .. local arrays ..
  integer(i4), parameter :: maxpatches = 50
  integer(i4) ic(0:maxpatches) ! number of bndry gridpoints within patch coloured i, 0 <=i<= maxpatches
  logical closed(0:maxpatches) ! .true. if the boundary patch is closed (i.e. a profile)
  integer(i4), allocatable :: iwork(:)
!     ..
!     .. external subroutines ..
  external findcolours, myroutine, setbndrynodelist
!
!     count the nof boundary vertices
!
  nbpoin = 0
  do ipoin = 1, npoin
    if (nodcod(ipoin) .gt. 0) nbpoin = nbpoin + 1
  end do
  write (6, *) 'setbndrynodeptr: found ', nbpoin, ' boundary points'
  allocate (nodptr(nbpoin, 3), source=0_i4)
  allocate (iwork(2*nbpoin), source=0_i4)
  call myroutine(nodcod, bndfac, nbfac, nodptr, npoin, nbpoin, iwork)
  deallocate (iwork) ! workspace
  call findcolours(bndfac, nodptr, ic, closed, maxpatches, nbfac, nbpoin, nclr)
  write (6, *) 'setbndrynodeptr: found ', nclr, ' colours or bndry patches'
!
  allocate (ia(nclr + 1), source=0_i4) ! pointer in csr format
  allocate (iclr(nclr), source=0_i4) ! colour of the k-th boundary patch
!
!     we use two pointers ia(1:nclr+1) and ja(1:nnzr) to address the bndry gridpoints
!     gridpoints belonging to boundary (or patch) i, where 1<=i<=nclr
!     are stored in ja(j), j=jbgn,jend where
!     jbgn = ia(i), jend = ia(i+1)-1
!     the total nof entries in ja is nnzr = ia(nclr+1)-ia(1);
!     note that ia(1) = 1
!
!     therefore, in order to address all gridpoints of patch coloured 3, we do:
!
!     jbgn = ia(3)
!     jend = ia(4) -1
!     do j = jbgn,jend
!        ipoin = ja(j) ! global node number
!     enddo
!
!     here we first setup the ia array
!
  ia(1) = 1
  k = 0
  do ipoin = 0, maxpatches
    if (ic(ipoin) .gt. 0) then ! skip empty colours
      k = k + 1
      if (k .gt. nclr) then
        write (6, *) 'setbndrynodeptr: found too many colours; check&
        & nclr !'
        error stop 1
      end if
      ia(k + 1) = ia(k) + ic(ipoin)
      iclr(k) = ipoin ! colour of the k-th boundary patch
    end if
  end do
  nnzr = ia(nclr + 1) - ia(1)
  write (6, *) 'setbndrynodeptr: ', nnzr, ' entries expected in ja'
!
!     allocate ja
!
  allocate (ja(nnzr), source=0_i4)
!
  write (6, *) 'setbndrynodeptr: finished initializing'
  write (6, *) 'setbndrynodeptr: now calling setbndrynodelist'
  write (6, *)
!
!     here we fill the ja array
!
  call setbndrynodelist(ia, ja, iclr, nclr, closed, bndfac, nbfac, nodptr, nbpoin)
end subroutine setbndrynodeptr
!
subroutine myroutine(nodcode, ibndptr, nbfac, inodptr, npoin, nbpoin,&
&iwork)
!
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  external binsrc, isortrx
!
  integer(i4) nbfac, npoin, nbpoin
  integer(i4) nodcode(*)
  integer(i4) ibndptr(3, nbfac), inodptr(nbpoin, 3), iwork(2*nbpoin)
  integer(i4) ipoin, ipos, last, ifail, j, k, iface
  logical verbose
  parameter(verbose=.false.)
!     parameter(verbose=.true.)
!
!     inodptr is a nodal pointer for boundary nodes
!     inodptr(1,*) addresses the global nodenumber
!     inodptr(2,*) addresses one of the two edges it belongs to
!     inodptr(3,*) addresses the other edge it belongs to
!
!     ibndptr is a nodal pointer for boundary edges
!     ibndptr(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
!     ibndptr(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
!     ibndptr(3,*) is the colour of the boundary face
!
!     fill the first entry of the pointer with the node number
!
  last = 0
  do ipoin = 1, npoin
    if (nodcode(ipoin) .gt. 0) then
      last = last + 1
      iwork(last) = ipoin
    end if
  end do
  if (last .eq. nbpoin) then
    write (8, *) 'found ', nbpoin, ' boundary points'
  else
    error stop 'last .ne. nbpoin'
  end if
!
!     iwork(1:nbpoin) stores the nbpoin node numbers
!
  call isortrx(nbpoin, iwork, iwork(nbpoin + 1))
  do ipoin = 1, nbpoin
    inodptr(ipoin, 1) = iwork(iwork(nbpoin + ipoin))
  end do
!
  ifail = 0
  outer: do iface = 1, nbfac
    inner: do j = 1, 2
      ipoin = ibndptr(j, iface)
      call binsrc(ipoin, inodptr(1, 1), nbpoin, ipos, last)
      if (ipos .eq. 0) then
        write (6, *) 'subr. myroutine: entry not found for ',&
        &ipoin
        error stop 1
      end if
      do k = 2, 3
        if (inodptr(ipos, k) .eq. 0) then
          inodptr(ipos, k) = iface
          cycle inner
        end if
      end do
!
!    the node seems to belong to more than 2 boundary faces
!
      ifail = ipoin
      write (6, *) 'subr. myroutine: the node seems to belong to more&
      &    than 2 boundary faces'
      write (6, *) 'subr. myroutine: face no. ', iface, (ibndptr(k, iface),&
      &k=1, 3)
      write (6, *) 'subr. myroutine: node no. ', ipoin, (inodptr(ipos, k), k&
      &=1, 3)
      exit outer
    end do inner
  end do outer
!
  if (ifail .ne. 0 .or. verbose) then
    if (ifail .ne. 0) error stop 'unrecoverable error in checkbndrypntr'
  end if
!
end subroutine myroutine

subroutine findcolours(ibndptr, inodptr, ic, closed, maxpatches,&
&nbfac, nbpoin, nc)
!
  use mod_kinds, only: wp, i4
  implicit none(type, external)
!
  integer(i4) nbfac, nbpoin, nc, maxpatches
  integer(i4) ibndptr(3, nbfac), inodptr(nbpoin, 3), ic(0:*)
  logical closed(0:*)
  integer(i4) ipoin, ipos, i, j, ibc, ie
  integer(i4) iflg(2)

  ! 13/06/2020 -- bugfix by prof. bonfiglioli
  ! do ibc = 1, maxpatches
  ! to avoid nasty initialization of ibc(0)
  do ibc = 0, maxpatches
    closed(ibc) = .true. ! true if the boundary is closed (e.g. airfoil)
    ic(ibc) = 0
  end do
!
!     here we count the nof boundary patches (nc) i.e. boundary surfaces
!     with different colours
!
  nc = 0
  do ipos = 1, nbfac
    ibc = ibndptr(3, ipos)
    if ((ibc .lt. 0) .or. (ibc .gt. maxpatches)) then
      write (6, *) 'subr. findcolours: boundary colour ', ibc, ' is&
      &          outside the range ', 0, maxpatches
      error stop 2
    end if
    ic(ibc) = ic(ibc) + 1
  end do
  do ibc = 1, maxpatches
    if (ic(ibc) .ne. 0) then
      nc = nc + 1
      write (6, *) 'subr. findcolours: boundary coloured ', ibc, ' ha&
      &s ', ic(ibc), ' edges'
    end if
  end do
  write (6, *)
  write (6, *) 'subr. findcolours: there are ', nc, ' different boundary&
  & patches'
  write (6, *)
!
!     check whether the patches are open or closed:
!     whenever a boundary gridpoint belongs to two edges coloured differently,
!     those two colours belong to an open boundary
!     otherwise the boundary is closed
!
  do j = 1, nbpoin ! loop over boundary points
    ipoin = inodptr(j, 1) ! the global nodenumber of the j-th boundary point
    do i = 2, 3 ! loop over the edges that share boundary point j
      ie = inodptr(j, i)
      iflg(i - 1) = ibndptr(3, ie) ! this is the colour
    end do
    if (iflg(1) .ne. iflg(2)) then ! identify the boundary gridpoints that belong to bndry edges of different colours
      write (6, *) 'subr. findcolours: gridpoint ', ipoin, ' belongs to&
      & patches ', iflg(1), ' and ', iflg(2)
      closed(iflg(1)) = .false.
      closed(iflg(2)) = .false.
    end if
  end do ! end loop over bndry gridpoints
!
  do ibc = 1, maxpatches
    if (ic(ibc) .ne. 0) then
      write (6, *) 'subr. findcolours: bndry patch ', ibc, ' has ',&
      &ic(ibc), ' bndry edges; closed is ', closed(ibc)
!
!     if the boundary patch is closed, the nof boundary points equals the nof bndry edges
!     otherwise it equals the nof bndry edges+1
!     here we reset ic which is no longer the nof bndry edges, but the nof boundary points
!     with color ibc
!
      if (closed(ibc)) then
!                 ic(ibc) = ic(ibc)
      else
        ic(ibc) = ic(ibc) + 1 ! ic will return the nof bndry points lying on the bndry coloured ibc
      end if
    end if  !
  end do
end subroutine findcolours
!
subroutine setbndrynodelist(ia, ja, iclr, nclr, closed, ibndptr, nbfac,&
&inodptr, nbpoin)
!
!     this routine finds the list of bndry gridpoints belonging to bndry iclr
!
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  external nearby
!
  integer(i4) nbfac, nbpoin, nclr
  integer(i4) ia(*), ja(*), iclr(nclr)
  integer(i4) ibndptr(3, nbfac), inodptr(nbpoin, 3)
  integer(i4) ipoin, i, j, last, now, iend, ibgn, istart, ipatch, ibc, nc, ie
  integer(i4) iflg(2), neighb(2)
  logical closed(0:*)
  logical verbose
  logical found
  logical advanced
!     parameter(verbose=.true.)
  parameter(verbose=.false.)
!
!     input
!     inodptr is a nodal pointer for boundary nodes
!     inodptr(1,*) addresses the global nodenumber
!     inodptr(2,*) addresses one of the two edges it belongs to
!     inodptr(3,*) addresses the other edge it belongs to
!
!     input
!     ibndptr is a nodal pointer for boundary edges
!     ibndptr(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
!     ibndptr(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
!     ibndptr(3,*) is the colour of the boundary face
!
  do ipatch = 1, nclr  ! loo over all patches
    ibc = iclr(ipatch)
    write (6, *) 'setbndrynodelist: patch ', ipatch, ' has colour ', ibc,&
    &' expecting ', ia(ipatch + 1) - ia(ipatch), ' vertices'
    if (closed(ibc)) then
      write (6, *) 'setbndrynodelist: patch ', ipatch, ' is closed ',&
      &closed(ibc)
!
!     any bndry gridpoint that belongs to an adge coloured ibc is all right
!
      found = .false.
      search_closed: do j = 1, nbpoin ! loop over bndry gridpoints
        ipoin = inodptr(j, 1) ! the global nodenumber of the j-th boundary point
        do i = 2, 3 ! loop over the edges that share boundary point ipoin
          ie = inodptr(j, i)
          if (ibndptr(3, ie) .eq. ibc) then
            istart = ipoin
            found = .true.
            exit search_closed
          end if
        end do
      end do search_closed
      if (.not. found) then
        write (6, *) 'cannot find node; un-recoverable error in setbndr&
        &ynodelist (1)'
      end if
    else ! the boundary is open
      write (6, *) 'setbndrynodelist: patch ', ipatch, ' is open ',&
      &closed(ibc)
!
!     identify one of the two endpoints lying on patch coloured ibc:
!     for an open boundary, this is the endpoint that is shared btw two
!     bndry faces of different colour
!
      found = .false.
      search_open: do j = 1, nbpoin ! loop over bndry gridpoints
        ipoin = inodptr(j, 1) ! the global nodenumber of the j-th boundary point
        do i = 2, 3 ! loop over the edges that share boundary point ipoin
          ie = inodptr(j, i)
          iflg(i - 1) = ibndptr(3, ie)
        end do
        if (iflg(1) .ne. iflg(2)) then ! we have found the boundary gridpoints that belong to bndry edges of different colours
          write (6, *) 'setbndrynodelist: gridpoint ', ipoin,&
          &' belongs to patch ', iflg(1), ' and ', iflg(2)
          if ((iflg(1) .eq. ibc) .or. (iflg(2) .eq. ibc)) then
            istart = ipoin ! set the starting point
            found = .true.
            exit search_open
          end if
        end if
      end do search_open ! end loop over bndry gridpoints
!
      if (.not. found) then
        write (6, *) 'cannot find node; un-recoverable error in setbndr&
        &ynodelist (2)'
        error stop 2
      end if
    end if ! test on closed(*)
    ibgn = ia(ipatch)
    iend = ia(ipatch + 1) - 1
    ja(ibgn) = istart
    now = istart
    last = -1
    if (verbose) write (6, 200) ipatch, now, last, (neighb(j), j=1, 2)
    nc = iend - ibgn + 1  ! number of vertices expected on current patch
    traverse: do
      call nearby(now, ibndptr, inodptr, nbpoin, neighb, ibc)
      if (verbose) write (6, 200) ipatch, now, last, (neighb(j), j=1, 2)
!
!     neighb(1:2) gives the two vertices that surrount gridpoint now
!     and have the same bndry colour
!
      advanced = .false.
      do j = 1, 2
        if ((neighb(j) .ne. 0) .and. (neighb(j) .ne. last)) then
          last = now
          now = neighb(j)
          if (now .eq. istart) exit traverse ! this should only occur when the patch is closed
          ibgn = ibgn + 1
          ja(ibgn) = now
          advanced = .true.
          exit
        end if
      end do
      if (.not. advanced) exit traverse
    end do traverse
!
    if (ibgn .ne. iend) then
      write (6, *) 'is ', ibgn, ' = ', iend, ' ?????'
      write (6, *) (ja(j), j=ibgn, iend)
      error stop 4
    end if
    ibgn = ia(ipatch)
    write (6, *) 'subr. setbndrynodelist: found ', iend - ibgn + 1,&
    &' vertices in patch ', ipatch
    write (6, *) 'subr. setbndrynodelist: vertices are: ',&
    &(ja(j), j=ibgn, iend)
    write (6, *)
  end do ! end the outermost loop over patches
200 format(1x, 'patch ', i2, ' curr, prev, vertices and neighb are ', 4(i3, 1x))
end subroutine setbndrynodelist
!
subroutine nearby(inode, ibndptr, inodptr, nbpoin, neighb, ibc)
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  external binsrc
  integer(i4) inode, nbpoin, ibc ! input
!     inode is a global nodenumber
  integer(i4) neighb(*) ! output
  integer(i4) ibndptr(3, *), inodptr(nbpoin, 3) ! input
!
!     inodptr is a nodal pointer for boundary nodes
!     inodptr(1,*) addresses the global nodenumber
!     inodptr(2,*) addresses one of the two edges it belongs to
!     inodptr(3,*) addresses the other edge it belongs to
  integer(i4) ipos, n1, jbc, ie, i, j, last, k
  logical verbose
  parameter(verbose=.false.)
!     parameter(verbose=.true.)
!
!     we need to find one of the two bndry vertices neighbouring inode
!     that shares the same bndry code ibc; if the boundary is closed
!     there are two candidates, if it is open there must be only one
!
  neighb(1) = 0
  neighb(2) = 0
!     find the location ipos in inodptr where inode is stored
  call binsrc(inode, inodptr(1, 1), nbpoin, ipos, last)
  if (ipos .eq. 0) then
    write (6, *) 'subr. nearby(): entry not found for ',&
    &inode
    error stop 1
  else
    if (verbose)&
    &write (6, *) 'node ', inode, ' found in entry ', ipos,&
    &inodptr(ipos, 1)
  end if
!
!                 inode
!     +-------------+-------------+
!            ^             ^
!            |             |
!     inodptr(ipos,2) inodptr(ipos,3)
!
  last = 0
  do i = 2, 3 ! loop over the edges that share boundary point inode
    ie = inodptr(ipos, i)
    jbc = ibndptr(3, ie) ! colour of the neighbouring boundary face
    if (verbose)&
    &write (6, *) 'node ', (inodptr(ipos, k), k=1, 3), ' edge', i - 1, ' is ', ie&
    &, ' coloured ', jbc, ' with verts ', (ibndptr(j, ie), j=1, 2)
!
!        pick up the edge that has colour ibc
!
    if (jbc .eq. ibc) then ! the neighbouring face has the same colour
      do j = 1, 2  ! loop over the two vertices of the boundary face
        n1 = ibndptr(j, ie)
        if ((n1 .ne. inode)) then
          last = last + 1
          if (last .gt. 2) then ! a node on a boundary must not have more than two neighbours
            error stop 'there is smthg very wrong'
          end if
          neighb(last) = n1
        end if
      end do ! end loop over the two vertices of the neighbouring edges
    else
      if (verbose)&
      &write (6, *) 'skipping face ', ie, ' has colour ', jbc,&
      &'rather than ', ibc
!                    neighb(j) = 0 ! j might be uninitialized
    end if
  end do ! end loop over the two edges that meet at a bndry gridpoint
  if (verbose)&
  &write (6, *) 'node ', inode, ' has neighb ', (neighb(k), k=1, 2)
end subroutine nearby
!
subroutine check(ia, ja, iclr, nclr, corg, ndim)
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  integer(i4) ndim, nclr
  integer(i4) ia(*), ja(*), iclr(nclr)
  real(wp) corg(ndim, *)
  integer(i4) i, j, k, jbgn, jend, l, nnzr
  character*24 fname
  fname = "bndry00.dat"
  nnzr = ia(nclr + 1) - ia(1)
  write (6, *) 'nnzr = ', nnzr
  do i = 1, nclr
    write (fname(6:7), fmt="(I2.2)") i
    open (10, file=fname)
    write (10, *) '# patch ', i, ' has colour ', iclr(i)
    jbgn = ia(i)
    jend = ia(i + 1) - 1
    write (6, *) 'jbgn, jend = ', jbgn, jend
    do j = jbgn, jend
      k = ja(j)
      write (10, *) (corg(l, k), l=1, ndim)
    end do
    close (10)
  end do
end subroutine check
