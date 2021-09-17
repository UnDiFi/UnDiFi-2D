      DOUBLE PRECISION FUNCTION EDDY(VISCT,NTURB,DENSIT,VISCL,NOFVERT)
C
C     RETURNS a cell-averaged value of the turbulent viscosity
C     computed using the Spalart-Allmaras model
C

      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER NOFVERT,NTURB
      DOUBLE PRECISION DENSIT,VISCL 
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VISCT(NTURB,*)
C     ..
C     .. Local Scalars ..
      INTEGER I
      DOUBLE PRECISION TVI,TCHI
C
      DOUBLE PRECISION TFV1
      EXTERNAL TFV1
C     ..
      TVI = 0.D0
      DO 1 I = 1,NOFVERT
          TVI = TVI + VISCT(1,I)
    1 CONTINUE
      TVI = TVI/NOFVERT
      TVI = MAX(0.d0,TVI)
      TCHI = (DENSIT*TVI)/VISCL
      EDDY = DENSIT*TVI*TFV1(TCHI)

      RETURN

      END
