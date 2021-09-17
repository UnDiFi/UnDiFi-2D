      DOUBLE PRECISION FUNCTION FUN_BETA_ALDO(X)
C
C     The argument (X) MUST be Mach^2 - 1.00
C
C


C
C     .. Parameters ..
      REAL*8 EPS_SONIC,EPSQR_SONIC
      PARAMETER (EPS_SONIC=0.10d0,EPSQR_SONIC=EPS_SONIC**2)
      REAL*8 EPS_STAG,EPSQR_STAG
      PARAMETER (EPS_STAG=0.10d0,EPSQR_STAG=EPS_STAG**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,DSQRT
C     ..
      A = 0.25D0* (EPS_SONIC** (-1.5D0))
      C = 0.75D0* (EPS_SONIC**0.5D0)
C
      IF (ABS(X).LT.EPS_SONIC) THEN
          FUN_BETA_ALDO = A*X*X + C

      ELSE
          FUN_BETA_ALDO = DSQRT(DABS(X))
      ENDIF
C
      RETURN

      END
