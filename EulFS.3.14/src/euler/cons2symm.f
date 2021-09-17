      SUBROUTINE CONS2SYMM(ZROE,DVDU,NOFVAR,NDIM)
C
C     Watch out, there seem to be a bug somewhere here 
C
C     $Id: cons2symm.f,v 1.5 2013/01/29 14:33:34 abonfi Exp $
C
      IMPLICIT NONE
C
C
C     Transformation matrix from symmetrizing to conserved variables
C
C     .. Parameters ..
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DVDU(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AMACH2,ASQR,DENS,DENSINV,GM1ASQRINV,GM1RAINV,
     +                 TEMP,UX,UY,UZ,ZINV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      TEMP = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF (NDIM.EQ.3) TEMP = TEMP + ZROE(5)*ZROE(5)
      AMACH2 = TEMP
      TEMP = HALF*TEMP
C
C     TEMP = 0.5 * DENS * (U**2+V**2+W**2)
C
      ZINV = ONE/ZROE(1)
      DENS = ZROE(1)*ZROE(1)
      ASQR = GM1* (ZROE(2)-TEMP*ZINV)*ZINV
      AMACH2 = AMACH2/ASQR/DENS
C
      GM1ASQRINV = GM1/ASQR
      GM1RAINV = GM1/DENS/SQRT(ASQR)
      DENSINV = ZINV*ZINV
C
      UX = ZROE(3)*ZINV
      UY = ZROE(4)*ZINV
C
      DVDU(1,1) = ONE - HALF*GM1*AMACH2
      DVDU(2,1) = GM1RAINV*TEMP*DENSINV
      DVDU(3,1) = -UX*DENSINV
      DVDU(4,1) = -UY*DENSINV
C
      DVDU(1,2) = -GM1ASQRINV
      DVDU(2,2) = GM1RAINV
      DVDU(3,2) = ZERO
      DVDU(4,2) = ZERO
C
      DVDU(1,3) = GM1ASQRINV*UX
      DVDU(2,3) = -GM1RAINV*UX
      DVDU(3,3) = DENSINV
      DVDU(4,3) = ZERO
C
      DVDU(1,4) = GM1ASQRINV*UY
      DVDU(2,4) = -GM1RAINV*UY
      DVDU(3,4) = ZERO
      DVDU(4,4) = DENSINV
C
      IF (NDIM.EQ.2) RETURN
      UZ = ZROE(5)*ZINV
C
      DVDU(5,1) = -UZ*DENSINV
      DVDU(5,2) = ZERO
      DVDU(5,3) = ZERO
      DVDU(5,4) = ZERO
C
      DVDU(1,5) = GM1ASQRINV*UZ
      DVDU(2,5) = -GM1RAINV*UZ
      DVDU(3,5) = ZERO
      DVDU(4,5) = ZERO
      DVDU(5,5) = DENSINV
C
      RETURN

      END
