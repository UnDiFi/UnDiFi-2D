      SUBROUTINE WTRI(IBNDFAC,NBFAC,ICELNOD,NOFVERT,COOR,NDIM,ZROE,
     +                NOFVAR,NODCOD,NPOIN,FNAME,ICHOLE,NHOLE)

C
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVERT,NPOIN,NBFAC,NOFVAR,NHOLE
C     ..
C     .. Array Arguments ..
      CHARACTER*(*) FNAME
      CHARACTER     FWORK*255
      DOUBLE PRECISION COOR(NDIM,*),ZROE(NOFVAR,NPOIN)
      INTEGER IBNDFAC(3,*),ICELNOD(NOFVERT,*),NODCOD(*),ICHOLE(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DD,DUMMY,HH,XYSAVE
      DOUBLE PRECISION XYHOLE(2)
      INTEGER IA,K,IELEM,IFACE,IPOIN,IVERT,NELEM,IBODY,ITYPE
C     ..
C     .. External Functions ..
      DOUBLE PRECISION RAND
      INTEGER  ICYCL,LENSTR
      EXTERNAL ICYCL,LENSTR,RAND
C     ..
C     .. Intrinsic Functions ..
C     ..
      K = LENSTR(FNAME)
      write(6,*)' K = ',K
      FWORK(1:K+5) = FNAME(1:K)//".node"
      OPEN(UNIT=19,FILE=FWORK(1:K+5))
C     node files
C
C     * First line: <# of vertices> <dimension (must be 2)> <# of
C     * attributes> <# of boundary markers (0 or 1)>
C     * Remaining lines: <vertex #> <x> <y> [attributes]
C     * [boundary marker] 
C
C     Blank lines and comments prefixed by `#' may be placed
C     anywhere. Vertices must be numbered consecutively, starting from one or zero.
C
      WRITE(19,*)NPOIN,NDIM,NOFVAR,1
      HH = 1.d0/SQRT(REAL(NPOIN))
!     WRITE(6,*)'Typical mesh size is set to = ',HH
      DO 92 IPOIN = 1,NPOIN
!         WRITE(19,FMT=300)IPOIN,(COOR(IA,IPOIN),IA=1,NDIM),
!    &    (ZROE(IA,IPOIN),IA=1,NOFVAR),NODCOD(IPOIN)
          WRITE(19,FMT=*)IPOIN,(COOR(IA,IPOIN),IA=1,NDIM),
     &    (ZROE(IA,IPOIN),IA=1,NOFVAR),NODCOD(IPOIN)
C
   92 CONTINUE
      CLOSE(19)
Caldo
Caldo REM: non scrive il file poly!!!!!!!!
Caldo
Caldo
      RETURN
Caldo
      FWORK(1:K+5) = FNAME(1:K)//".poly"
      OPEN(UNIT=19,FILE=FWORK(1:K+5))
C     empty node list; this is just a poly file
      WRITE(19,*)0,NDIM,0,0
      WRITE(19,*)NBFAC,1
      DO 90 IFACE = 1,NBFAC
          IELEM = IBNDFAC(1,IFACE)
          IVERT = IBNDFAC(2,IFACE)
          ITYPE = IBNDFAC(3,IFACE)
          WRITE(19,FMT=*)IFACE,!-1
     &(ICELNOD(ICYCL(IVERT+IA,NOFVERT),IELEM),IA=1,NDIM),ITYPE
   90 CONTINUE
C     write nof holes
      WRITE(19,*)NHOLE
      DO K = 1,NHOLE
         CALL HOLEXY(ICHOLE(K),IBNDFAC,NBFAC,ICELNOD,NOFVERT,COOR,NDIM,
     &               XYHOLE)
         WRITE(19,*)K,(XYHOLE(IA),IA=1,NDIM)
      ENDDO
      CLOSE(19)
      RETURN
  300 FORMAT(I7,6(1X,F17.9),1X,I3)
      END
      SUBROUTINE HOLEXY(IBC,IBNDFAC,NBFAC,ICELNOD,NOFVERT,COOR,NDIM,
     &                  XYHOLE)
C
C     compute the barycenter of the hole; may not work for non-convex holes
C
      IMPLICIT NONE
      INTEGER IBC,NBFAC,NDIM,NOFVERT
      INTEGER IBNDFAC(3,NBFAC),ICELNOD(NOFVERT,*)
      DOUBLE PRECISION COOR(NDIM,*),XYHOLE(NDIM)
      INTEGER I,IELEM,IVERT,J,IPOIN,IC
      INTEGER ICYCL
      XYHOLE(1)  = 0.d0
      XYHOLE(2)  = 0.d0
      IC = 0
      DO I = 1,NBFAC
         IF(IBNDFAC(3,I).EQ.IBC)THEN
            IC = IC+1
            IELEM = IBNDFAC(1,I)
            IVERT = IBNDFAC(2,I)
            DO J = 1,2
               IPOIN = ICELNOD(ICYCL(IVERT+J,NOFVERT),IELEM)
               XYHOLE(1) = XYHOLE(1) + COOR(1,IPOIN)*0.5d0
               XYHOLE(2) = XYHOLE(2) + COOR(2,IPOIN)*0.5d0
            ENDDO
         ENDIF
      ENDDO
      XYHOLE(1) = XYHOLE(1) / REAL(IC)
      XYHOLE(2) = XYHOLE(2) / REAL(IC)
      RETURN
      END
