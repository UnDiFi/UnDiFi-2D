C
      SUBROUTINE MATSUM(A,LDA,B,LDB,NR,NC)
C
      IMPLICIT NONE
C
C     A := A + B
C
      INTEGER LDA,LDB,NR,NC
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
      INTEGER I,J
      DO 1 J = 1 , NC
         DO 1 I = 1 , NR
            A(I,J) = A(I,J) + B(I,J)
   1  CONTINUE
      RETURN
      END
C
