      SUBROUTINE SYMM2CONS(ZROE,DUDV,NOFVAR,NDIM)
C
C     $Id: symm2cons.f,v 1.4 2013/01/29 14:33:34 abonfi Exp $
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
      DOUBLE PRECISION DUDV(NOFVAR,NOFVAR),ZROE(NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AINV,ASQR,DENS,RUX,RUY,RUZ,TEMP,UX,UY,UZ,ZINV,
     +                 ZINVSQR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
      TEMP = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF (NDIM.EQ.3) TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF*TEMP
C
      ZINV = ONE/ZROE(1)
      DENS = ZROE(1)*ZROE(1)
      ZINVSQR = ZINV*ZINV
      ASQR = GM1* (ZROE(2)-TEMP*ZINV)*ZINV
      AINV = ONE/DSQRT(ASQR)
C
      UX = ZROE(3)*ZINV
      UY = ZROE(4)*ZINV
      RUX = ZROE(3)*ZROE(1)
      RUY = ZROE(4)*ZROE(1)
C
      DUDV(1,1) = ONE
      DUDV(2,1) = TEMP*ZINVSQR
      DUDV(3,1) = UX
      DUDV(4,1) = UY
C
      DUDV(1,2) = DENS*AINV
      DUDV(2,2) = ZROE(1)*ZROE(2)*AINV
      DUDV(3,2) = RUX*AINV
      DUDV(4,2) = RUY*AINV
C
      DUDV(1,3) = ZERO
      DUDV(2,3) = RUX
      DUDV(3,3) = DENS
      DUDV(4,3) = ZERO
C
      DUDV(1,4) = ZERO
      DUDV(2,4) = RUY
      DUDV(3,4) = ZERO
      DUDV(4,4) = DENS
C
      IF (NDIM.EQ.2) RETURN
      UZ = ZROE(5)*ZINV
      RUZ = ZROE(5)*ZROE(1)
C
      DUDV(5,1) = UZ
      DUDV(5,2) = RUZ*AINV
      DUDV(5,3) = ZERO
      DUDV(5,4) = ZERO
C
      DUDV(1,5) = ZERO
      DUDV(2,5) = RUZ
      DUDV(3,5) = ZERO
      DUDV(4,5) = ZERO
      DUDV(5,5) = DENS
C
      RETURN

      END
