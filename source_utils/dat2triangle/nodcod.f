      SUBROUTINE NODCOD(NODCODE,NPOIN,NBPOIN,ICELNOD,NOFVERT,NELEM,
     +                  IBNDFAC,NBFAC)
C
      IMPLICIT NONE
C
C
C     ... NODCODE (ISTAK(LNODCOD)) is an integer flag used to
C         distinguish between internal and boundary nodes ...
C
C
C ICELNOD -- Integer ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
C IBNDFAC -- Integer IBNDFAC(1:3,1:NBFAC)
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
C     nodcode(i) = 0 --> internal node
C     nodcode(i) > 0 --> boundary node shared by nodcode(i)
C     boundary faces
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NBPOIN,NELEM,NOFVERT,NPOIN
C     ..
C     .. Array Arguments ..
      INTEGER IBNDFAC(3,NBFAC),ICELNOD(NOFVERT,NELEM),NODCODE(NPOIN)
C     ..
C     .. Local Scalars ..
      INTEGER IBC,IELEM,IFACE,IPOIN,IVERT,JVERT,KVERT
C     ..
C     .. External Functions ..
      INTEGER ICYCL
      EXTERNAL ICYCL
C     ..
      DO 100 IFACE = 1,NBFAC
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
          IBC = IBNDFAC(3,IFACE)
          DO 90 JVERT = 1,NOFVERT - 1
              KVERT = ICYCL(IVERT+JVERT,NOFVERT)
              IPOIN = ICELNOD(KVERT,IELEM)
              NODCODE(IPOIN) = NODCODE(IPOIN) + 1
   90     CONTINUE
  100 CONTINUE
C
C     ... Count the boundary nodes ...
C
      DO 80 IPOIN = 1,NPOIN
          IF (NODCODE(IPOIN).NE.0) NBPOIN = NBPOIN + 1
   80 CONTINUE
      WRITE (6,FMT=300) NBPOIN
C
      RETURN

  300 FORMAT (15X,'BOUNDARY POINTS (NBPOIN) = ',I6)

      END
