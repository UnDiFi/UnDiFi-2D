      SUBROUTINE EVALTSA(IELEM,DUMMY1,DUMMY2,NDIM,NOFVERT,NOFVAR,NTURB,
     &                   DUMMY3,DUMMY4,DUMMY5,DUMMY6,LFLAG,SDUMMY1,
     &                   SDUMMY2,ZTURB,COMPRESSIBLE,VISCL,VISCT)
C
C     $Id: tsa1.f,v 1.3 2013/01/26 12:01:29 abonfi Exp $
C
C
C     RETURNS a cell-averaged value of the turbulent viscosity
C     computed using the Spalart-Allmaras model
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB
      DOUBLE PRECISION DENSIT,VISCL
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5,DUMMY6
      LOGICAL LFLAG,COMPRESSIBLE
      EXTERNAL SDUMMY1,SDUMMY2
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ZTURB(NTURB,*)
C     ..
C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION VISCT,TCHI
C
      DOUBLE PRECISION TFV1,SUTHERLAW
      EXTERNAL TFV1,SUTHERLAW
C     ..
      VISCT = ZERO
      DO 1 I = 1,NOFVERT
          VISCT = VISCT + ZTURB(1,I)
    1 CONTINUE
      VISCT = VISCT/NOFVERT
      IF(COMPRESSIBLE)THEN
         VISCL = SUTHERLAW(ZERO,ABAR,ASQR)
         DENSIT= UAVG(1)
      ELSE
         VISCL = ONE
         DENSIT= ONE
      ENDIF
      VISCT = MAX(ZERO,VISCT)
      TCHI = (DENSIT*VISCT)/VISCL
      VISCT = DENSIT*VISCT*TFV1(TCHI)

      RETURN

      END
