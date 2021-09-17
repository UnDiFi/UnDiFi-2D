      DOUBLE PRECISION FUNCTION SUTHERLAW(DUMMY,ABAR,ASQR)
C
      IMPLICIT NONE
C
      INCLUDE 'suther.com'
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION ABAR,ASQR,DUMMY
C     ..
C     .. Local Scalars ..
C     ..
      SUTHERLAW = C1*ASQR*ABAR*C4/(C2*ASQR+C3)
C
      RETURN

      END
