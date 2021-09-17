C
C ------------------------------ + ------------------------------
C
      LOGICAL FUNCTION NULLMAT( A, N, M, LDA, TOLER )
C
      IMPLICIT NONE
C
C     This function tests a matrix A(N,M) to see
C     whether it is the Null matrix within a given
C     tolerance TOLER 
C
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
C
C     .. Intrinsic Functions ..
C
      INTRINSIC ABS
C
C     .. Executable Statements ..
C
      DO 1 J = 1,M
	 DO 1 I = 1,N
            IF( ABS( A(I,J) ) .GT. TOLER )THEN
                NULLMAT = .FALSE.
                RETURN
            ENDIF
    1 CONTINUE
      NULLMAT = .TRUE.
      RETURN
      END
