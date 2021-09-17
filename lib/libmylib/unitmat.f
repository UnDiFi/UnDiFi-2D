C
C ------------------------------ + ------------------------------
C
      LOGICAL FUNCTION UNITMAT( A, N, M, LDA, TOLER )
C
      IMPLICIT NONE
C
C     This function tests a matrix A(N,M) to see
C     whether it is the Identity matrix up to given 
C     tolerance TOLER 
C
      DOUBLE PRECISION ONE
      PARAMETER(ONE=1.D0)
C
C     .. Scalar Arguments ..
C
      INTEGER N,M,LDA
      DOUBLE PRECISION TOLER
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION A(LDA,*)
C
C     .. Local Scalars ..
C
      INTEGER I,J
      DOUBLE PRECISION S
C
C     .. Intrinsic Functions ..
C
      INTRINSIC ABS
C
C     .. Executable Statements ..
C
      DO 1 J = 1,M
	 DO 1 I = 1,N
            IF( I .EQ. J )THEN
               S = A(I,J) - ONE
            ELSE
               S = A(I,J)
            ENDIF 
            IF( ABS( S ) .GT. TOLER )THEN
                UNITMAT = .FALSE.
                RETURN
            ENDIF
    1 CONTINUE
      UNITMAT = .TRUE.
      RETURN
      END
