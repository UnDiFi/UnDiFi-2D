C
      SUBROUTINE MATSUB(A,LDA,B,LDB,NR,NC)
C
      IMPLICIT NONE 
C
C     A := A - B
C
C     .. Scalar Arguments ..
      INTEGER LDA,LDB,NC,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
      DO 1 J = 1,NC
          DO 1 I = 1,NR
              A(I,J) = A(I,J) - B(I,J)
    1 CONTINUE
      RETURN

      END
