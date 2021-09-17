C
C
      LOGICAL FUNCTION UNITMAT(A,N,M,LDA,TOLER)
C
      IMPLICIT NONE 
C
C     This function tests a matrix A(N,M) to see
C     whether it is the Identity matrix up to given
C     tolerance TOLER
C
C
C
C
C
C
C
C
C
C
C
C     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOLER
      INTEGER LDA,M,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S
      INTEGER I,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      DO 1 J = 1,M
          DO 1 I = 1,N
              IF (I.EQ.J) THEN
                  S = A(I,J) - ONE

              ELSE
                  S = A(I,J)
              ENDIF

              IF (ABS(S).GT.TOLER) THEN
                  UNITMAT = .FALSE.
                  RETURN

              ENDIF

    1 CONTINUE
      UNITMAT = .TRUE.
      RETURN

      END
