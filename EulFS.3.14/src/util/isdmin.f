      INTEGER FUNCTION ISDMIN( N , X , INCI )
      IMPLICIT NONE

      INTEGER I,N,INCI

      REAL*8 X(N)
      REAL*8 XMIN,ABS

      ISDMIN = 1
      XMIN = ABS(X(1))
      DO I = 2 , N , INCI
         IF (ABS(X(I)) .LT. XMIN)THEN
          ISDMIN = I
          XMIN = ABS(X(I))
         ENDIF
      ENDDO
      RETURN
      END
