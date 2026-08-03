!> @LBNDFAC[in] bndry face pointer read from the poly file
!> @LNODCOD[in] nodal pointe read from the node or poly file
!> @NBFAC[in] nof bndry faces
!> @NPOIN[in] nof vertices
!> @NBPOIN[in] nof boundary vertices
!> @LNODPTR[out] pointer for boundary vertices
!> @LIAO[out] IA(1:NCLR+1)
!> @LJAO[out] JA(1:NNZR)
!> @LICLR[out] ICLR(1:NCLR) colour of the NCLR patches
subroutine SetBndryNodePtr(LBNDFAC, LNODCOD, NBFAC, NPOIN,&
&NBPOIN, LNODPTR, LIAO, LJAO, LICLR, NCLR)
!
  implicit none(type, external)
!
!     $Id: setbndrynodeptr.f,v 1.4 2018/08/06 09:18:23 abonfi Exp abonfi $
!
!     NDIM  is the space dimension =2
!     NVT = NDIM+1 is the number of vertices
!
!     NPSHMAX     : max number of shock points.
!     NESHMAX     : max number of shock elements
!
!     .. Parameters ..

!     .. Local Scalars ..
  integer LBNDFAC, LNODCOD
  integer LNODPTR, LWKSP
  integer NBFAC, NPOIN, IPOIN, NBPOIN
  integer IFAIL, I, K, NCLR, NNZR
  integer LIAO, LJAO, LICLR

!     .. Local Arrays ..

  integer MAXPATCHES
  parameter(MAXPATCHES=50)
  integer IC(0:MAXPATCHES) ! number of bndry gridpoints within patch coloured i, 0 <=i<= MAXPATCHES
  logical CLOSED(0:MAXPATCHES) ! .TRUE. if the boundary patch is closed (i.e. a profile)
!     ..
!     .. External Subroutines ..
  external DINIT, IINIT, ISTKIN, ISTKRL
  external FINDCOLOURS, MYROUTINE, SetBndryNodeList
!     ..
!
!     CHARACTER*(*) FNAME
!
!     ..
!     .. Arrays in Common ..
  double precision DSTAK(1)
  integer ISTAK(1)
  common/CSTAK/DSTAK
!     ..
!     .. Equivalences ..
  equivalence(DSTAK(1), ISTAK(1))
!     ..
!     .. External Functions ..
  integer ISTKGT
  external ISTKGT
!
!     count the nof boundary vertices
!
  NBPOIN = 0
  do 1 IPOIN = 0, NPOIN - 1
    if (ISTAK(LNODCOD + IPOIN) .gt. 0) NBPOIN = NBPOIN + 1
1   continue
    write (6, *) 'SetBndryNodePtr: Found ', NBPOIN, ' boundary points'
    LNODPTR = ISTKGT(3*NBPOIN, 2)
    LWKSP = ISTKGT(2*NBPOIN, 2)
    call IINIT(3*NBPOIN, 0, ISTAK(LNODPTR), 1)
    call IINIT(2*NBPOIN, 0, ISTAK(LWKSP), 1)
    call MYROUTINE(ISTAK(LNODCOD), ISTAK(LBNDFAC), NBFAC, ISTAK(LNODPTR),&
    &NPOIN, NBPOIN, ISTAK(LWKSP))
    call ISTKRL(1) ! release the workspace
    call FINDCOLOURS(ISTAK(LBNDFAC), ISTAK(LNODPTR),&
    &IC, CLOSED, MAXPATCHES, NBFAC, NBPOIN, NCLR)
!     call exit(0)
    write (6, *) 'SetBndryNodePtr: Found ', NCLR, ' colours or bndry patche&
    &s'
!
    LIAO = ISTKGT(NCLR + 1, 2) ! pointer in CSR format
    LICLR = ISTKGT(NCLR, 2) ! colour of the K-th boundary patch
!
!     We use two pointers IA(1:NCLR+1) and JA(1:NNZR) to address the bndry gridpoints
!     Gridpoints belonging to boundary (or patch) I, where 1<=I<=NCLR
!     are stored in JA(J), J=JBGN,JEND where
!     JBGN = IA(I), JEND = IA(I+1)-1
!     the total nof entries in JA is NNZR = IA(NCLR+1)-IA(1);
!     note that IA(1) = 1
!
!     therefore, in order to address all gridpoints of patch coloured 3, we do:
!
!     JBGN = IA(3)
!     JEND = IA(4) -1
!     DO J = JBGN,JEND
!        IPOIN = JA(J) ! global node number
!     ENDDO
!
!     Here we first setup the IA array
!
    ISTAK(LIAO) = 1
    K = 0
    do I = 0, MAXPATCHES
      if (IC(I) .gt. 0) then ! skip empty colours
        K = K + 1
        if (K .gt. NCLR) then
          write (6, *) 'SetBndryNodePtr: Found too many colours; check&
          & NCLR !'
          stop
        end if
        ISTAK(LIAO + K) = ISTAK(LIAO + K - 1) + IC(I)
        ISTAK(LICLR + K - 1) = I ! colour of the K-th boundary patch
      end if
    end do
    NNZR = ISTAK(LIAO + NCLR) - ISTAK(LIAO)
    write (6, *) 'SetBndryNodePtr: ', NNZR, ' entries expected in ja'
!
!     allocate JA
!
    LJAO = ISTKGT(NNZR, 2)
    call IINIT(NNZR, 0, ISTAK(LJAO), 1)
!
    write (6, *) 'SetBndryNodePtr: finished initializing'
    write (6, *) 'SetBndryNodePtr: now calling SetBndryNodeList'
    write (6, *)
!
!     Here we fill the JA array
!
    call SetBndryNodeList(ISTAK(LIAO), ISTAK(LJAO), ISTAK(LICLR),&
    &NCLR, CLOSED, ISTAK(LBNDFAC), NBFAC, ISTAK(LNODPTR), NBPOIN)
    return
    end subroutine SetBndryNodePtr
!
    subroutine MYROUTINE(NODCODE, IBNDPTR, NBFAC, INODPTR, NPOIN, NBPOIN,&
    &IWORK)
!
      implicit none(type, external)
      external BINSRC, ISORTRX
!
      integer NBFAC, NPOIN, NBPOIN
      integer NODCODE(*)
      integer IBNDPTR(3, NBFAC), INODPTR(NBPOIN, 3), IWORK(2*NBPOIN)
      integer IPOIN, IPOS, LAST, IFAIL, J, K, IFACE
      logical VERBOSE
      parameter(VERBOSE=.false.)
!     PARAMETER(VERBOSE=.TRUE.)
!
!     INODPTR is a nodal pointer for boundary nodes
!     INODPTR(1,*) addresses the global nodenumber
!     INODPTR(2,*) addresses one of the two edges it belongs to
!     INODPTR(3,*) addresses the other edge it belongs to
!
!     IBNDPTR is a nodal pointer for boundary edges
!     IBNDPTR(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
!     IBNDPTR(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
!     IBNDPTR(3,*) is the colour of the boundary face
!
!     fill the first entry of the pointer with the node number
!
      LAST = 0
      do 1 IPOIN = 1, NPOIN
        if (NODCODE(IPOIN) .gt. 0) then
          LAST = LAST + 1
          IWORK(LAST) = IPOIN
        end if
1       continue
        if (LAST .eq. NBPOIN) then
          write (8, *) 'Found ', NBPOIN, ' boundary points'
        else
          stop 'LAST .NE. NBPOIN'
        end if
!
!     IWORK(1:NBPOIN) stores the NBPOIN node numbers
!
        call ISORTRX(NBPOIN, IWORK, IWORK(NBPOIN + 1))
        do 2 IPOIN = 1, NBPOIN
          INODPTR(IPOIN, 1) = IWORK(IWORK(NBPOIN + IPOIN))
2         continue
!     write(6,*)(inodptr(ipoin,1),ipoin=1,nbpoin)
!     pause
!
!     CALL X04EAF('General',' ',NBPOIN,3,INODPTR,NBPOIN,
!    +            'Nodal Bndry pointer',IFAIL)
!     CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!    +            'Face Bndry pointer',IFAIL)
!
          IFAIL = 0
          do 4 IFACE = 1, NBFAC
            do 5 J = 1, 2
              IPOIN = IBNDPTR(J, IFACE)
              call BINSRC(IPOIN, INODPTR(1, 1), NBPOIN, IPOS, LAST)
              if (IPOS .eq. 0) then
                write (6, *) 'Subr. Myroutine: Entry NOT found for ',&
                &IPOIN
                stop
              end if
              do 6 K = 2, 3
                if (INODPTR(IPOS, K) .eq. 0) then
                  INODPTR(IPOS, K) = IFACE
                  goto 5
                end if
6               continue
!
!    the node seems to belong to more than 2 boundary faces
!
                IFAIL = IPOIN
                write (6, *) 'Subr. Myroutine: the node seems to belong to more&
                &    than 2 boundary faces'
                write (6, *) 'Subr. Myroutine: Face no. ', IFACE, (IBNDPTR(K, IFACE),&
                &K=1, 3)
                write (6, *) 'Subr. Myroutine: Node no. ', IPOIN, (INODPTR(IPOS, K), K&
                &=1, 3)
                goto 7
5               continue

4               continue
7               continue
!
                if (IFAIL .ne. 0 .or. VERBOSE) then
!         CALL X04EAF('General',' ',NBPOIN,3,INODPTR,NBPOIN,
!     +            'Nodal Bndry pointer',IFAIL)
!         CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!     +            'Face Bndry pointer',IFAIL)
                  if (IFAIL .ne. 0) stop 'Unrecoverable error in CheckBndryPntr'
                end if
!
                return
                end subroutine MYROUTINE

                subroutine FINDCOLOURS(IBNDPTR, INODPTR, IC, CLOSED, MAXPATCHES,&
                &NBFAC, NBPOIN, NC)
!
                  implicit none(type, external)
!
                  integer NBFAC, NBPOIN, NC, MAXPATCHES
                  integer IBNDPTR(3, NBFAC), INODPTR(NBPOIN, 3), IC(0:*)
                  logical CLOSED(0:*)
                  integer IPOIN, IPOS, I, J, IBC, IE
                  integer IFLG(2)

                  ! 13/06/2020 -- BugFix by Prof. Bonfiglioli
                  ! DO IBC = 1, MAXPATCHES
                  ! to avoid nasty initialization of IBC(0)
                  do IBC = 0, MAXPATCHES
                    CLOSED(IBC) = .true. ! TRUE if the boundary is closed (e.g. airfoil)
                    IC(IBC) = 0
                  end do
!
!     Here we count the nof boundary patches (NC) i.e. boundary surfaces
!     with different colours
!
                  NC = 0
                  do IPOS = 1, NBFAC
                    IBC = IBNDPTR(3, IPOS)
                    if ((IBC .lt. 0) .or. (IBC .gt. MAXPATCHES)) then
                      write (6, *) 'Subr. FindColours: Boundary colour ', IBC, ' is&
                      &          outside the range ', 0, MAXPATCHES
                      error stop 2
                    end if
                    IC(IBC) = IC(IBC) + 1
                  end do
                  do IBC = 1, MAXPATCHES
                    if (IC(IBC) .ne. 0) then
                      NC = NC + 1
                      write (6, *) 'Subr. FindColours: Boundary coloured ', IBC, ' ha&
                      &s ', IC(IBC), ' edges'
                    end if
                  end do
                  write (6, *)
                  write (6, *) 'Subr. FindColours: There are ', NC, ' different boundary&
                  & patches'
                  write (6, *)
!
!     check whether the patches are open or closed:
!     whenever a boundary gridpoint belongs to two edges coloured differently,
!     those two colours belong to an open boundary
!     otherwise the boundary is closed
!
                  do J = 1, NBPOIN ! loop over boundary points
                    IPOIN = INODPTR(J, 1) ! the global nodenumber of the J-th boundary point
                    do I = 2, 3 ! loop over the edges that share boundary point J
                      IE = INODPTR(J, I)
                      IFLG(I - 1) = IBNDPTR(3, IE) ! this is the colour
                    end do
                    if (IFLG(1) .ne. IFLG(2)) then ! identify the boundary gridpoints that belong to bndry edges of different colours
                      write (6, *) 'Subr. FindColours: Gridpoint ', IPOIN, ' belongs to&
                      & patches ', IFLG(1), ' and ', IFLG(2)
                      CLOSED(IFLG(1)) = .false.
                      CLOSED(IFLG(2)) = .false.
                    end if
                  end do ! end loop over bndry gridpoints
!
                  do IBC = 1, MAXPATCHES
                    if (IC(IBC) .ne. 0) then
                      write (6, *) 'Subr. FindColours: Bndry patch ', IBC, ' has ',&
                      &IC(IBC), ' bndry edges; closed is ', CLOSED(IBC)
!
!     if the boundary patch is closed, the nof boundary points equals the nof bndry edges
!     otherwise it equals the nof bndry edges+1
!     here we reset IC which is no longer the nof bndry edges, but the nof boundary points
!     with color IBC
!
                      if (CLOSED(IBC)) then
!                 IC(IBC) = IC(IBC)
                      else
                        IC(IBC) = IC(IBC) + 1 ! IC will return the nof bndry points lying on the bndry coloured IBC
                      end if
                    end if  !
                  end do
                  return
                end subroutine FINDCOLOURS
!
                subroutine SetBndryNodeList(IA, JA, ICLR, NCLR, CLOSED, IBNDPTR, NBFAC,&
                &INODPTR, NBPOIN)
!
!     this routine finds the list of bndry gridpoints belonging to bndry ICLR
!
                  implicit none(type, external)
                  external NEARBY
!
                  integer NBFAC, NBPOIN, NCLR
                  integer IA(*), JA(*), ICLR(NCLR)
                  integer IBNDPTR(3, NBFAC), INODPTR(NBPOIN, 3)
                  integer IPOIN, I, J, LAST, NOW, IEND, IBGN, ISTART, IPATCH, IBC, NC, IE
                  integer IFLG(2), NEIGHB(2)
                  logical CLOSED(0:*)
                  logical VERBOSE
!     PARAMETER(VERBOSE=.TRUE.)
                  parameter(VERBOSE=.false.)
!
!     Input
!     INODPTR is a nodal pointer for boundary nodes
!     INODPTR(1,*) addresses the global nodenumber
!     INODPTR(2,*) addresses one of the two edges it belongs to
!     INODPTR(3,*) addresses the other edge it belongs to
!
!     Input
!     IBNDPTR is a nodal pointer for boundary edges
!     IBNDPTR(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
!     IBNDPTR(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
!     IBNDPTR(3,*) is the colour of the boundary face
!
                  do 120 IPATCH = 1, NCLR  ! loo over all patches
                    IBC = ICLR(IPATCH)
                    write (6, *) 'SetBndryNodeList: Patch ', IPATCH, ' has colour ', IBC,&
                    &' expecting ', IA(IPATCH + 1) - IA(IPATCH), ' vertices'
                    if (CLOSED(IBC)) then
                      write (6, *) 'SetBndryNodeList: Patch ', IPATCH, ' is CLOSED ',&
                      &CLOSED(IBC)
!
!     Any bndry gridpoint that belongs to an adge coloured IBC is all right
!
                      do J = 1, NBPOIN ! loop over bndry gridpoints
                        IPOIN = INODPTR(J, 1) ! the global nodenumber of the J-th boundary point
                        do I = 2, 3 ! loop over the edges that share boundary point IPOIN
                          IE = INODPTR(J, I)
                          if (IBNDPTR(3, IE) .eq. IBC) then
                            ISTART = IPOIN
                            goto 100
                          end if
                        end do
                      end do
                      write (6, *) 'Cannot find node; un-recoverable error in SetBndr&
                      &yNodeList (1)'
                    else ! the boundary is open
                      write (6, *) 'SetBndryNodeList: Patch ', IPATCH, ' is OPEN ',&
                      &CLOSED(IBC)
!
!     identify one of the two endpoints lying on patch coloured IBC:
!     for an open boundary, this is the endpoint that is shared btw two
!     bndry faces of different colour
!
                      do J = 1, NBPOIN ! loop over bndry gridpoints
                        IPOIN = INODPTR(J, 1) ! the global nodenumber of the J-th boundary point
                        do I = 2, 3 ! loop over the edges that share boundary point IPOIN
                          IE = INODPTR(J, I)
                          IFLG(I - 1) = IBNDPTR(3, IE)
                        end do
                        if (IFLG(1) .ne. IFLG(2)) then ! we have found the boundary gridpoints that belong to bndry edges of different colours
                          write (6, *) 'SetBndryNodeList: Gridpoint ', IPOIN,&
                          &' belongs to patch ', IFLG(1), ' and ', IFLG(2)
                          if ((IFLG(1) .eq. IBC) .or. (IFLG(2) .eq. IBC)) then
                            ISTART = IPOIN ! set the starting point
                            goto 100
                          end if
                        end if
                      end do ! end loop over bndry gridpoints
!
                      write (6, *) 'Cannot find node; un-recoverable error in SetBndr&
                      &yNodeList (2)'
                      error stop 2
                    end if ! test on CLOSED(*)
100                 continue
                    IBGN = IA(IPATCH)
                    IEND = IA(IPATCH + 1) - 1
                    JA(IBGN) = ISTART
                    NOW = ISTART
                    LAST = -1
                    if (VERBOSE) write (6, 200) IPATCH, NOW, LAST, (NEIGHB(J), J=1, 2)
                    NC = IEND - IBGN + 1  ! number of vertices expected on current patch
2                   call NEARBY(NOW, IBNDPTR, INODPTR, NBPOIN, NEIGHB, IBC)
                    if (VERBOSE) write (6, 200) IPATCH, NOW, LAST, (NEIGHB(J), J=1, 2)
!
!     NEIGHB(1:2) gives the two vertices that surrount gridpoint NOW
!     AND have the same bndry colour
!
                    do 10 J = 1, 2
                      if ((NEIGHB(J) .ne. 0) .and. (NEIGHB(J) .ne. LAST)) then
                        LAST = NOW
                        NOW = NEIGHB(J)
                        if (NOW .eq. ISTART) then ! this should only occur when the patch is closed
                          goto 12
                        else
                          IBGN = IBGN + 1
                          JA(IBGN) = NOW
                          goto 2
                        end if
                      end if
10                    continue
!
12                    continue
                      if (IBGN .ne. IEND) then
                        write (6, *) 'Is ', IBGN, ' = ', IEND, ' ?????'
                        write (6, *) (ja(j), j=ibgn, iend)
                        error stop 4
                      end if
                      IBGN = IA(IPATCH)
                      write (6, *) 'Subr. SetBndryNodeList: Found ', IEND - IBGN + 1,&
                      &' vertices in patch ', IPATCH
                      write (6, *) 'Subr. SetBndryNodeList: vertices are: ',&
                      &(ja(j), j=ibgn, iend)
                      write (6, *)
120                   continue ! end the outermost loop over patches
200                   format(1x, 'Patch ', I2, ' curr, prev, vertices and neighb are ', 4(I3, 1x))
                      return
                      end subroutine SetBndryNodeList
!
                      subroutine NEARBY(INODE, IBNDPTR, INODPTR, NBPOIN, NEIGHB, IBC)
                        implicit none(type, external)
                        external BINSRC
                        integer INODE, NBPOIN, IBC ! Input
!     INODE is a GLOBAL nodenumber
                        integer NEIGHB(*) ! Output
                        integer IBNDPTR(3, *), INODPTR(NBPOIN, 3) ! Input
!
!     INODPTR is a nodal pointer for boundary nodes
!     INODPTR(1,*) addresses the global nodenumber
!     INODPTR(2,*) addresses one of the two edges it belongs to
!     INODPTR(3,*) addresses the other edge it belongs to
                        integer IPOS, N1, JBC, IE, I, J, LAST, K
                        logical VERBOSE
                        parameter(VERBOSE=.false.)
!     PARAMETER(VERBOSE=.TRUE.)
!
!     we need to find one of the two bndry vertices neighbouring INODE
!     that shares the same bndry code IBC; if the boundary is closed
!     there are two candidates, if it is open there must be only one
!
                        NEIGHB(1) = 0
                        NEIGHB(2) = 0
!     find the location IPOS in INODPTR where INODE is stored
                        call BINSRC(INODE, INODPTR(1, 1), NBPOIN, IPOS, LAST)
                        if (IPOS .eq. 0) then
                          write (6, *) 'Subr. Nearby(): Entry NOT found for ',&
                          &INODE
                          stop
                        else
                          if (VERBOSE)&
                          &write (6, *) 'Node ', INODE, ' found in entry ', IPOS,&
                          &INODPTR(IPOS, 1)
                        end if
!     write(6,*)
!     write(6,*)
!     write(6,*)
!
!                 INODE
!     +-------------+-------------+
!            ^             ^
!            |             |
!     INODPTR(IPOS,2) INODPTR(IPOS,3)
!
                        LAST = 0
                        do I = 2, 3 ! loop over the edges that share boundary point INODE
                          IE = INODPTR(IPOS, I)
                          JBC = IBNDPTR(3, IE) ! colour of the neighbouring boundary face
                          if (VERBOSE)&
                          &write (6, *) 'Node ', (INODPTR(IPOS, K), K=1, 3), ' edge', i - 1, ' is ', IE&
                          &, ' coloured ', JBC, ' with verts ', (IBNDPTR(J, IE), j=1, 2)
!
!        Pick up the edge that has colour IBC
!
                          if (JBC .eq. IBC) then ! the neighbouring face has the same colour
                            do J = 1, 2  ! loop over the two vertices of the boundary face
                              N1 = IBNDPTR(J, IE)
                              if ((N1 .ne. INODE)) then
                                LAST = LAST + 1
                                if (LAST .gt. 2) then ! A node on a boundary must not have more than two neighbours
                                  stop 'There is smthg very wrong'
                                end if
                                NEIGHB(LAST) = N1
                              end if
                            end do ! end loop over the two vertices of the neighbouring edges
                          else
                            if (VERBOSE)&
                            &write (6, *) 'Skipping face ', IE, ' has colour ', JBC,&
                            &'rather than ', IBC
!                    NEIGHB(J) = 0 ! J might be uninitialized
                          end if
                        end do ! end loop over the two edges that meet at a bndry gridpoint
                        if (VERBOSE)&
                        &write (6, *) 'Node ', INODE, ' has neighb ', (neighb(k), k=1, 2)
                        return
                      end subroutine NEARBY
!
                      subroutine CHECK(IA, JA, ICLR, NCLR, CORG, NDIM)
                        implicit none(type, external)
                        integer NDIM, NCLR
                        integer IA(*), JA(*), ICLR(NCLR)
                        double precision CORG(NDIM, *)
                        integer I, J, K, JBGN, JEND, L, NNZR
                        character*24 FNAME
                        FNAME = "bndry00.dat"
                        nnzr = ia(nclr + 1) - ia(1)
                        write (6, *) 'nnzr = ', nnzr
                        do I = 1, NCLR
                          write (FNAME(6:7), FMT="(I2.2)") I
                          open (10, FILE=FNAME)
                          write (10, *) '# patch ', I, ' has colour ', ICLR(I)
                          JBGN = IA(I)
                          JEND = IA(I + 1) - 1
                          write (6, *) 'jbgn, jend = ', jbgn, jend
                          do J = JBGN, JEND
                            K = JA(J)
                            write (10, *) (CORG(L, K), L=1, NDIM)
                          end do
                          close (10)
                        end do
                        return
                      end subroutine CHECK
