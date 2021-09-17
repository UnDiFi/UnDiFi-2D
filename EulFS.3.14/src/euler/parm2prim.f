      SUBROUTINE PARM_PRIM(ZROE,DVDZ,NOFVAR,NDIM)
C
C     $Id: parm2prim.f,v 1.4 2013/01/29 14:33:34 abonfi Exp $
C
      IMPLICIT NONE
C
C
C     .. Parameters ..
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DVDZ(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUM,DUMSQR
C     ..
      DUM = ONE/ZROE(1)
      DUMSQR = DUM*DUM
c
c       .. First row ..
c
      DVDZ(1,1) = TWO*ZROE(1)
      DVDZ(1,2) = ZERO
      DVDZ(1,3) = ZERO
      DVDZ(1,4) = ZERO
c
c       .. Second row ..
c
      DVDZ(2,1) = GM1OG*ZROE(2)
      DVDZ(2,2) = GM1OG*ZROE(1)
      DVDZ(2,3) = -GM1OG*ZROE(3)
      DVDZ(2,4) = -GM1OG*ZROE(4)
c
c       .. Third row ..
c
      DVDZ(3,1) = -ZROE(3)*DUMSQR
      DVDZ(3,2) = ZERO
      DVDZ(3,3) = DUM
      DVDZ(3,4) = ZERO
c
c       .. Fourth row ..
c
      DVDZ(4,1) = -ZROE(4)*DUMSQR
      DVDZ(4,2) = ZERO
      DVDZ(4,3) = ZERO
      DVDZ(4,4) = DUM
C
      IF (NDIM.EQ.2) RETURN
C
      DVDZ(1,5) = ZERO
      DVDZ(2,5) = -GM1OG*ZROE(5)
      DVDZ(3,5) = ZERO
      DVDZ(4,5) = ZERO
c
c       .. Fifth row ..
c
      DVDZ(5,1) = -ZROE(5)*DUMSQR
      DVDZ(5,2) = ZERO
      DVDZ(5,3) = ZERO
      DVDZ(5,4) = ZERO
      DVDZ(5,5) = DUM
C
      RETURN

      END
