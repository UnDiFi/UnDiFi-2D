!> @LBNDFAC[in] bndry face pointer read from the poly file
!> @LNODCOD[in] nodal pointe read from the node or poly file
!> @NBFAC[in] nof bndry faces
!> @NPOIN[in] nof vertices
!> @NBPOIN[in] nof boundary vertices
!> @LNODPTR[out] pointer for boundary vertices
!> @LIAO[out] IA(1:NCLR+1)
!> @LJAO[out] JA(1:NNZR)
!> @LICLR[out] ICLR(1:NCLR) colour of the NCLR patches
      SUBROUTINE SetBndryNodePtr(LBNDFAC,LNODCOD,NBFAC,NPOIN,
     &NBPOIN,LNODPTR,LIAO,LJAO,LICLR,NCLR)
C
      IMPLICIT NONE
C
C     $Id: setbndrynodeptr.f,v 1.4 2018/08/06 09:18:23 abonfi Exp abonfi $
C
C     NDIM  is the space dimension =2
C     NVT = NDIM+1 is the number of vertices
C
c     NPSHMAX     : max number of shock points.
c     NESHMAX     : max number of shock elements
C
C     .. Parameters ..

C     .. Local Scalars ..
      INTEGER LBNDFAC,LNODCOD
      INTEGER LNODPTR,LWKSP
      INTEGER NBFAC,NPOIN,IPOIN,NBPOIN
      INTEGER IFAIL,I,K,NCLR,NNZR
      INTEGER LIAO,LJAO,LICLR

C     .. Local Arrays ..

      INTEGER MAXPATCHES
      PARAMETER(MAXPATCHES=50)
      INTEGER IC(0:MAXPATCHES) ! number of bndry gridpoints within patch coloured i, 0 <=i<= MAXPATCHES
      LOGICAL CLOSED(0:MAXPATCHES) ! .TRUE. if the boundary patch is closed (i.e. a profile)
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT,ISTKIN,ISTKRL
C     ..
C
c     CHARACTER*(*) FNAME
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
      COMMON/CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
C     .. External Functions ..
      INTEGER  ISTKGT
      EXTERNAL ISTKGT
C
C     count the nof boundary vertices
C
      NBPOIN = 0
      DO 1 IPOIN = 0, NPOIN-1
          IF( ISTAK(LNODCOD+IPOIN) .GT. 0 )NBPOIN = NBPOIN+1
    1 CONTINUE
      write(6,*)'SetBndryNodePtr: Found ',NBPOIN,' boundary points'
      LNODPTR = ISTKGT(3*NBPOIN,2)
      LWKSP   = ISTKGT(2*NBPOIN,2)
      CALL IINIT(3*NBPOIN,0,ISTAK(LNODPTR),1)
      CALL IINIT(2*NBPOIN,0,ISTAK(LWKSP),1)
      CALL MYROUTINE(ISTAK(LNODCOD),ISTAK(LBNDFAC),NBFAC,ISTAK(LNODPTR),
     &               NPOIN,NBPOIN,ISTAK(LWKSP))
      CALL ISTKRL(1) ! release the workspace
      CALL FINDCOLOURS(ISTAK(LBNDFAC),ISTAK(LNODPTR),
     &                 IC,CLOSED,MAXPATCHES,NBFAC,NBPOIN,NCLR)
!     call exit(0)
      write(6,*)'SetBndryNodePtr: Found ',NCLR,' colours or bndry patche
     &s'
C
      LIAO = ISTKGT(NCLR+1,2) ! pointer in CSR format
      LICLR = ISTKGT(NCLR,2) ! colour of the K-th boundary patch
C
C     We use two pointers IA(1:NCLR+1) and JA(1:NNZR) to address the bndry gridpoints
C     Gridpoints belonging to boundary (or patch) I, where 1<=I<=NCLR
C     are stored in JA(J), J=JBGN,JEND where
C     JBGN = IA(I), JEND = IA(I+1)-1
C     the total nof entries in JA is NNZR = IA(NCLR+1)-IA(1);
C     note that IA(1) = 1
C
C     therefore, in order to address all gridpoints of patch coloured 3, we do:
C
C     JBGN = IA(3)
C     JEND = IA(4) -1
C     DO J = JBGN,JEND
C        IPOIN = JA(J) ! global node number
C     ENDDO
C
C     Here we first setup the IA array
C
      ISTAK(LIAO) = 1
      K = 0
      DO I = 0, MAXPATCHES
         IF(IC(I).GT.0)THEN ! skip empty colours
            K = K + 1
            IF(K.GT.NCLR)THEN
               WRITE(6,*)'SetBndryNodePtr: Found too many colours; check
     & NCLR !'
               STOP
            ENDIF
            ISTAK(LIAO+K) = ISTAK(LIAO+K-1) + IC(I)
            ISTAK(LICLR+K-1) = I ! colour of the K-th boundary patch
         ENDIF
      ENDDO
      NNZR = ISTAK(LIAO+NCLR)-ISTAK(LIAO)
      write(6,*)'SetBndryNodePtr: ',NNZR,' entries expected in ja'
C
C     allocate JA
C
      LJAO = ISTKGT(NNZR,2)
      CALL IINIT(NNZR,0,ISTAK(LJAO),1)
!
      write(6,*)'SetBndryNodePtr: finished initializing'
      write(6,*)'SetBndryNodePtr: now calling SetBndryNodeList'
      write(6,*)
C
C     Here we fill the JA array
C
      CALL SetBndryNodeList(ISTAK(LIAO),ISTAK(LJAO),ISTAK(LICLR),
     &        NCLR,CLOSED,ISTAK(LBNDFAC),NBFAC,ISTAK(LNODPTR),NBPOIN)
      RETURN
      END
C
      SUBROUTINE MYROUTINE(NODCODE,IBNDPTR,NBFAC,INODPTR,NPOIN,NBPOIN,
     &IWORK)
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NPOIN,NBPOIN
      INTEGER NODCODE(*)
      INTEGER IBNDPTR(3,NBFAC),INODPTR(NBPOIN,3),IWORK(2*NBPOIN)
      INTEGER IPOIN,IPOS,LAST,IFAIL,J,K,IFACE
      LOGICAL VERBOSE
      PARAMETER(VERBOSE=.FALSE.)
!     PARAMETER(VERBOSE=.TRUE.)
C
C     INODPTR is a nodal pointer for boundary nodes
C     INODPTR(1,*) addresses the global nodenumber
C     INODPTR(2,*) addresses one of the two edges it belongs to
C     INODPTR(3,*) addresses the other edge it belongs to
C
C     IBNDPTR is a nodal pointer for boundary edges
C     IBNDPTR(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
C     IBNDPTR(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
C     IBNDPTR(3,*) is the colour of the boundary face
C
C     fill the first entry of the pointer with the node number
C
      LAST = 0
      DO 1 IPOIN = 1, NPOIN
          IF( NODCODE(IPOIN) .GT. 0 )THEN
              LAST = LAST+1
              IWORK(LAST) = IPOIN
          ENDIF
    1 CONTINUE
      IF(LAST.EQ.NBPOIN)THEN
         write(8,*)'Found ',NBPOIN,' boundary points'
      ELSE
         STOP 'LAST .NE. NBPOIN'
      ENDIF
C
C     IWORK(1:NBPOIN) stores the NBPOIN node numbers
C
      CALL ISORTRX(NBPOIN,IWORK,IWORK(NBPOIN+1))
      DO 2 IPOIN = 1, NBPOIN
         INODPTR(IPOIN,1) = IWORK(IWORK(NBPOIN+IPOIN))
    2 CONTINUE
!     write(6,*)(inodptr(ipoin,1),ipoin=1,nbpoin)
!     pause
C
!     CALL X04EAF('General',' ',NBPOIN,3,INODPTR,NBPOIN,
!    +            'Nodal Bndry pointer',IFAIL)
!     CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!    +            'Face Bndry pointer',IFAIL)
C
      IFAIL = 0
      DO 4 IFACE = 1, NBFAC
         DO 5 J = 1,2
         IPOIN = IBNDPTR(J,IFACE)
         CALL BINSRC(IPOIN,INODPTR(1,1),NBPOIN,IPOS,LAST)
         IF(IPOS.EQ.0)THEN
            write(6,*)'Subr. Myroutine: Entry NOT found for ',
     &IPOIN
            STOP
         ENDIF
         DO 6 K = 2,3
            IF(INODPTR(IPOS,K).EQ.0)THEN
               INODPTR(IPOS,K)=IFACE
               GOTO 5
            ENDIF
    6    CONTINUE
C
C    the node seems to belong to more than 2 boundary faces
C
         IFAIL = IPOIN
         write(6,*)'Subr. Myroutine: the node seems to belong to more
     & than 2 boundary faces'
         write(6,*)'Subr. Myroutine: Face no. ',IFACE,(IBNDPTR(K,IFACE),
     &K=1,3)
         write(6,*)'Subr. Myroutine: Node no. ',IPOIN,(INODPTR(IPOS,K),K
     &=1,3)
         GOTO 7
    5 CONTINUE

    4 CONTINUE
    7 CONTINUE
C
      IF(IFAIL.NE.0.OR.VERBOSE)THEN
!         CALL X04EAF('General',' ',NBPOIN,3,INODPTR,NBPOIN,
!     +            'Nodal Bndry pointer',IFAIL)
!         CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!     +            'Face Bndry pointer',IFAIL)
        IF(IFAIL.NE.0) STOP 'Unrecoverable error in CheckBndryPntr'
      ENDIF
C
      RETURN
      END

      SUBROUTINE FINDCOLOURS(IBNDPTR,INODPTR,IC,CLOSED,MAXPATCHES,
     &                       NBFAC,NBPOIN,NC)
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NBPOIN,NC,MAXPATCHES
      INTEGER IBNDPTR(3,NBFAC),INODPTR(NBPOIN,3),IC(0:*)
      LOGICAL CLOSED(0:*)
      INTEGER IPOIN,IPOS,I,J,IBC,IE
      INTEGER IFLG(2)

      ! 13/06/2020 -- BugFix by Prof. Bonfiglioli
      ! DO IBC = 1, MAXPATCHES
      ! to avoid nasty initialization of IBC(0)
      DO IBC = 0, MAXPATCHES
           CLOSED(IBC) = .TRUE. ! TRUE if the boundary is closed (e.g. airfoil)
           IC(IBC) = 0
      ENDDO
C
C     Here we count the nof boundary patches (NC) i.e. boundary surfaces
C     with different colours
C
      NC = 0
      DO IPOS = 1, NBFAC
         IBC = IBNDPTR(3,IPOS)
         IF((IBC.LT.0).OR.(IBC.GT.MAXPATCHES))THEN
               WRITE(6,*)'Subr. FindColours: Boundary colour ',IBC,' is
     &         outside the range ',0,MAXPATCHES
               CALL EXIT(2)
         ENDIF
         IC(IBC) = IC(IBC) + 1
      ENDDO
      DO IBC = 1, MAXPATCHES
           IF(IC(IBC).NE.0) THEN
              NC = NC + 1
              WRITE(6,*)'Subr. FindColours: Boundary coloured ',IBC,' ha
     &s ',IC(IBC),' edges'
           ENDIF
      ENDDO
      WRITE(6,*)
      WRITE(6,*)'Subr. FindColours: There are ',NC,' different boundary
     &patches'
      WRITE(6,*)
C
C     check whether the patches are open or closed:
C     whenever a boundary gridpoint belongs to two edges coloured differently,
C     those two colours belong to an open boundary
C     otherwise the boundary is closed
C
      DO J = 1, NBPOIN ! loop over boundary points
         IPOIN = INODPTR(J,1) ! the global nodenumber of the J-th boundary point
         DO I = 2,3 ! loop over the edges that share boundary point J
            IE = INODPTR(J,I)
            IFLG(I-1) = IBNDPTR(3,IE) ! this is the colour
         ENDDO
         IF(IFLG(1).NE.IFLG(2))THEN ! identify the boundary gridpoints that belong to bndry edges of different colours
            WRITE(6,*)'Subr. FindColours: Gridpoint ',IPOIN,' belongs to
     & patches ',IFLG(1),' and ',IFLG(2)
            CLOSED(IFLG(1)) = .FALSE.
            CLOSED(IFLG(2)) = .FALSE.
         ENDIF
      ENDDO ! end loop over bndry gridpoints
C
      DO IBC = 1, MAXPATCHES
           IF(IC(IBC).NE.0) THEN
              WRITE(6,*)'Subr. FindColours: Bndry patch ',IBC,' has ',
     &IC(IBC),' bndry edges; closed is ',CLOSED(IBC)
C
C     if the boundary patch is closed, the nof boundary points equals the nof bndry edges
C     otherwise it equals the nof bndry edges+1
C     here we reset IC which is no longer the nof bndry edges, but the nof boundary points
C     with color IBC
C
               IF(CLOSED(IBC))THEN
!                 IC(IBC) = IC(IBC)
               ELSE
                  IC(IBC) = IC(IBC)+1 ! IC will return the nof bndry points lying on the bndry coloured IBC
               ENDIF
           ENDIF  !
      ENDDO
      RETURN
      END
C
      SUBROUTINE SetBndryNodeList(IA,JA,ICLR,NCLR,CLOSED,IBNDPTR,NBFAC,
     &INODPTR,NBPOIN)
C
C     this routine finds the list of bndry gridpoints belonging to bndry ICLR
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NBPOIN,NCLR
      INTEGER IA(*),JA(*),ICLR(NCLR)
      INTEGER IBNDPTR(3,NBFAC),INODPTR(NBPOIN,3)
      INTEGER IPOIN,I,J,LAST,NOW,IEND,IBGN,ISTART,IPATCH,IBC,NC,IE
      INTEGER IFLG(2),NEIGHB(2)
      LOGICAL CLOSED(0:*)
      LOGICAL VERBOSE
!     PARAMETER(VERBOSE=.TRUE.)
      PARAMETER(VERBOSE=.FALSE.)
C
C     Input
C     INODPTR is a nodal pointer for boundary nodes
C     INODPTR(1,*) addresses the global nodenumber
C     INODPTR(2,*) addresses one of the two edges it belongs to
C     INODPTR(3,*) addresses the other edge it belongs to
C
C     Input
C     IBNDPTR is a nodal pointer for boundary edges
C     IBNDPTR(1,*) addresses the global nodenumber of one of the two vertices of the boundary edges
C     IBNDPTR(2,*) addresses the global nodenumber of the other      vertex   of the boundary edges
C     IBNDPTR(3,*) is the colour of the boundary face
C
      DO 120 IPATCH = 1,NCLR  ! loo over all patches
         IBC = ICLR(IPATCH)
         WRITE(6,*)'SetBndryNodeList: Patch ',IPATCH,' has colour ',IBC,
     &' expecting ',IA(IPATCH+1)-IA(IPATCH),' vertices'
         IF( CLOSED(IBC) )THEN
            WRITE(6,*)'SetBndryNodeList: Patch ',IPATCH,' is CLOSED ',
     &CLOSED(IBC)
C
C     Any bndry gridpoint that belongs to an adge coloured IBC is all right
C
            DO J = 1, NBPOIN ! loop over bndry gridpoints
               IPOIN = INODPTR(J,1) ! the global nodenumber of the J-th boundary point
               DO I = 2,3 ! loop over the edges that share boundary point IPOIN
                  IE = INODPTR(J,I)
                  IF( IBNDPTR(3,IE) .EQ. IBC )THEN
                      ISTART = IPOIN
                      GOTO 100
                  ENDIF
               ENDDO
            ENDDO
            WRITE(6,*)'Cannot find node; un-recoverable error in SetBndr
     &yNodeList (1)'
         ELSE ! the boundary is open
            WRITE(6,*)'SetBndryNodeList: Patch ',IPATCH,' is OPEN ',
     &CLOSED(IBC)
C
C     identify one of the two endpoints lying on patch coloured IBC:
C     for an open boundary, this is the endpoint that is shared btw two
C     bndry faces of different colour
C
            DO J = 1, NBPOIN ! loop over bndry gridpoints
               IPOIN = INODPTR(J,1) ! the global nodenumber of the J-th boundary point
               DO I = 2,3 ! loop over the edges that share boundary point IPOIN
                  IE = INODPTR(J,I)
                  IFLG(I-1) = IBNDPTR(3,IE)
               ENDDO
               IF(IFLG(1).NE.IFLG(2))THEN ! we have found the boundary gridpoints that belong to bndry edges of different colours
                  WRITE(6,*)'SetBndryNodeList: Gridpoint ',IPOIN,
     &' belongs to patch ',IFLG(1),' and ',IFLG(2)
                  IF( (IFLG(1).EQ.IBC) .OR. (IFLG(2).EQ.IBC) )THEN
                    ISTART = IPOIN ! set the starting point
                    GOTO 100
                  ENDIF
               ENDIF
            ENDDO ! end loop over bndry gridpoints
C
            WRITE(6,*)'Cannot find node; un-recoverable error in SetBndr
     &yNodeList (2)'
            CALL EXIT(2)
         ENDIF ! test on CLOSED(*)
  100 CONTINUE
         IBGN = IA(IPATCH)
         IEND = IA(IPATCH+1)-1
         JA(IBGN) = ISTART
         NOW = ISTART
         LAST = -1
         IF(VERBOSE)WRITE(6,200)IPATCH,NOW,LAST,(NEIGHB(J),J=1,2)
         NC = IEND-IBGN+1  ! number of vertices expected on current patch
    2    CALL NEARBY(NOW,IBNDPTR,INODPTR,NBPOIN,NEIGHB,IBC)
         IF(VERBOSE)WRITE(6,200)IPATCH,NOW,LAST,(NEIGHB(J),J=1,2)
C
C     NEIGHB(1:2) gives the two vertices that surrount gridpoint NOW
C     AND have the same bndry colour
C
         DO 10 J = 1,2
            IF( (NEIGHB(J) .NE. 0).AND.(NEIGHB(J) .NE. LAST) )THEN
               LAST = NOW
               NOW = NEIGHB(J)
               IF(NOW.EQ.ISTART)THEN ! this should only occur when the patch is closed
                  GOTO 12
               ELSE
                  IBGN = IBGN + 1
                  JA(IBGN) = NOW
                  GOTO 2
               ENDIF
            ENDIF
   10    CONTINUE
C
   12 CONTINUE
         IF(IBGN.NE.IEND)THEN
            WRITE(6,*)'Is ',IBGN,' = ',IEND,' ?????'
            write(6,*)(ja(j),j=ibgn,iend)
            call exit(4)
         ENDIF
         IBGN = IA(IPATCH)
         write(6,*)'Subr. SetBndryNodeList: Found ',IEND-IBGN+1,
     &' vertices in patch ',IPATCH
         write(6,*)'Subr. SetBndryNodeList: vertices are: ',
     &(ja(j),j=ibgn,iend)
         write(6,*)
  120 CONTINUE ! end the outermost loop over patches
  200 FORMAT(1X,'Patch ',I2,' curr, prev, vertices and neighb are ',
     &4(I3,1X))
      RETURN
      END
C
      SUBROUTINE NEARBY(INODE,IBNDPTR,INODPTR,NBPOIN,NEIGHB,IBC)
      IMPLICIT NONE
      INTEGER INODE,NBPOIN,IBC ! Input
!     INODE is a GLOBAL nodenumber
      INTEGER NEIGHB(*) ! Output
      INTEGER IBNDPTR(3,*),INODPTR(NBPOIN,3) ! Input
C
C     INODPTR is a nodal pointer for boundary nodes
C     INODPTR(1,*) addresses the global nodenumber
C     INODPTR(2,*) addresses one of the two edges it belongs to
C     INODPTR(3,*) addresses the other edge it belongs to
      INTEGER IPOS,N1,JBC,IE,I,J,LAST,K
      LOGICAL VERBOSE
      PARAMETER(VERBOSE=.FALSE.)
!     PARAMETER(VERBOSE=.TRUE.)
C
C     we need to find one of the two bndry vertices neighbouring INODE
C     that shares the same bndry code IBC; if the boundary is closed
C     there are two candidates, if it is open there must be only one
C
      NEIGHB(1) = 0
      NEIGHB(2) = 0
C     find the location IPOS in INODPTR where INODE is stored
      CALL BINSRC(INODE,INODPTR(1,1),NBPOIN,IPOS,LAST)
      IF(IPOS.EQ.0)THEN
            write(6,*)'Subr. Nearby(): Entry NOT found for ',
     &INODE
            STOP
      ELSE
            IF(VERBOSE)
     &      write(6,*)'Node ',INODE,' found in entry ',IPOS,
     &INODPTR(IPOS,1)
      ENDIF
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
      DO I = 2,3 ! loop over the edges that share boundary point INODE
         IE = INODPTR(IPOS,I)
         JBC = IBNDPTR(3,IE) ! colour of the neighbouring boundary face
         IF(VERBOSE)
     &   write(6,*)'Node ',(INODPTR(IPOS,K),K=1,3),' edge',i-1,' is ',IE
     &,' coloured ',JBC,' with verts ',(IBNDPTR(J,IE),j=1,2)
!
!        Pick up the edge that has colour IBC
!
         IF(JBC.EQ.IBC)THEN ! the neighbouring face has the same colour
            DO J = 1,2  ! loop over the two vertices of the boundary face
               N1 = IBNDPTR(J,IE)
               IF ( (N1.NE.INODE) )THEN
                     LAST = LAST + 1
                     IF(LAST.GT.2)THEN ! A node on a boundary must not have more than two neighbours
                         STOP 'There is smthg very wrong'
                     ENDIF
                     NEIGHB(LAST) = N1
               ENDIF
            ENDDO ! end loop over the two vertices of the neighbouring edges
         ELSE
         IF(VERBOSE)
     &      WRITE(6,*)'Skipping face ',IE,' has colour ',JBC,
     &'rather than ',IBC
!                    NEIGHB(J) = 0 ! J might be uninitialized
         ENDIF
      ENDDO ! end loop over the two edges that meet at a bndry gridpoint
      IF(VERBOSE)
     &write(6,*)'Node ',INODE,' has neighb ',(neighb(k),k=1,2)
      RETURN
      END
C
      SUBROUTINE CHECK(IA,JA,ICLR,NCLR,CORG,NDIM)
      IMPLICIT NONE
      INTEGER NDIM,NCLR
      INTEGER IA(*),JA(*),ICLR(NCLR)
      DOUBLE PRECISION CORG(NDIM,*)
      INTEGER I,J,K,JBGN,JEND,L,NNZR
      CHARACTER*24 FNAME
      FNAME = "bndry00.dat"
      nnzr = ia(nclr+1)-ia(1)
      write(6,*)'nnzr = ',nnzr
      DO I = 1, NCLR
         WRITE(FNAME(6:7),FMT="(I2.2)")I
         OPEN(10,FILE=FNAME)
         WRITE(10,*)'# patch ',I,' has colour ',ICLR(I)
         JBGN = IA(I)
         JEND = IA(I+1)-1
      write(6,*)'jbgn, jend = ',jbgn,jend
         DO J = JBGN, JEND
              K = JA(J)
              WRITE(10,*)(CORG(L,K),L=1,NDIM)
         ENDDO
         CLOSE(10)
      ENDDO
      RETURN
      END
