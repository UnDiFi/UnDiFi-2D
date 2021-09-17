      SUBROUTINE NODCOD(NODCODE,NPOIN,NBPOIN,ICELNOD,NOFVERT,NELEM,
     +                  IBNDFAC,NBFAC)
C
      IMPLICIT NONE
C
!     INCLUDE 'io.com'
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NBPOIN,NELEM,NOFVERT,NPOIN
C
C     On entry:
C     --------
C     NBFAC   no. of processor owned boundary faces/edges;
C             global number of boundary faces/edges 
C             in the uni-processor case
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C
C     Upon return:
C     -----------
C     ..
C     .. Array Arguments ..
      INTEGER IBNDFAC(3,NBFAC),ICELNOD(NOFVERT,NELEM),NODCODE(NPOIN)
C
C     On entry:
C     --------
C     IBNDFAC -- Integer IBNDFAC(1:3,1:NBFAC)
C            Boundary face pointer :
C
C          + IBNDFAC(1,ibfac) gives the global number of the
C            element the ibfac-th face belongs to
C
C          + IBNDFAC(2,ibfac) gives the local number of the
C            vertex opposite the ibfac-th boundary face, ie.
C            its global node number is given by:
C            ICELNOD(IBNDFAC(2,ibfac),IBNDFAC(1,ibfac))
C
C          + IBNDFAC(3,ibfac) gives the boundary color
C            of the ibfac-th boundary face
C
C     ICELNOD -- Integer ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            local node number (in the multi-proc. case) or the
C            global node number (in the uni-proc. case) of 
C            the i-th vertex of the ielem-th cell
C
C     Upon return:
C     -----------
C     NODCODE -- Integer NODCOD(1:NPOIN)
C         gives the number of boundary edges/faces meeting
C         in a given meshpoint, i.e. it is an integer flag 
C         used to distinguish between internal and boundary nodes
C
C         NODCODE(i) = 0 --> internal node
C         NODCODE(i) > 0 --> boundary node shared by nodcode(i)
C         boundary faces
C
C
C
C     ..
C     .. Local Scalars ..
      INTEGER IELEM,IFACE,IPOIN,IVERT,JVERT,KVERT
C     ..
C     .. External Functions ..
      INTEGER ICYCL
      EXTERNAL ICYCL
C     ..
      DO 100 IFACE = 1,NBFAC
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
          IF(IVERT.LT.1.OR.IVERT.GT.NOFVERT)THEN
             WRITE(6,*)'Wrong IVERT ',IVERT 
             CALL EXIT(1)
          ELSEIF(IELEM.LT.1.OR.IELEM.GT.NELEM)THEN
             WRITE(6,*)'Wrong IELEM ',IELEM 
             CALL EXIT(1)
          ENDIF
          DO 90 JVERT = 1,NOFVERT - 1
              KVERT = ICYCL(IVERT+JVERT,NOFVERT)
              IPOIN = ICELNOD(KVERT,IELEM)
              NODCODE(IPOIN) = NODCODE(IPOIN) + 1
   90     CONTINUE
  100 CONTINUE
C
C     ... Count the boundary nodes ...
C
      NBPOIN = 0
      DO 80 IPOIN = 1,NPOIN
          IF (NODCODE(IPOIN).NE.0) NBPOIN = NBPOIN + 1
   80 CONTINUE
      WRITE (6,FMT=300) NBPOIN
C
      RETURN

  300 FORMAT (I6,1X,'BOUNDARY POINTS (NBPOIN) found by NODCOD ()')

      END
