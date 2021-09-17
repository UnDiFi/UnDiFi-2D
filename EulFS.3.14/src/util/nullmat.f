      LOGICAL FUNCTION NULLMAT(A,N,M,LDA,TOLER)
C
C     $Id: nullmat.f,v 1.1 2013/01/24 08:37:50 abonfi Exp $
C
      IMPLICIT NONE
C
C     This function tests a matrix A(N,M) to see
C     whether it is the Null matrix within a given
C     tolerance TOLER
C
C
C
C
C
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOLER
      INTEGER LDA,M,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      DO 1 J = 1,M
          DO 1 I = 1,N
              IF (ABS(A(I,J)).GT.TOLER) THEN
                  NULLMAT = .FALSE.
                  RETURN

              ENDIF

    1 CONTINUE
      NULLMAT = .TRUE.
      RETURN

      END
