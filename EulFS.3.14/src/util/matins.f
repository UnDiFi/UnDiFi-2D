C
      SUBROUTINE MATINS(A,NOFVAR,B,NORDER,IMAX,JMAX,IOFF)
C
C
C     A := A + B
C

C     .. Scalar Arguments ..
      INTEGER IOFF,NOFVAR,IMAX,JMAX,NORDER
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NOFVAR,NOFVAR,IMAX,JMAX),
     +                 B(NORDER,NORDER,IMAX,JMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KOFF,L,LOFF
C     ..
      DO 1 L = 1,NORDER
          LOFF = L + IOFF
          DO 1 K = 1,NORDER
              KOFF = K + IOFF
              DO 1 J = 1,JMAX
                  DO 1 I = 1,IMAX
                      A(KOFF,LOFF,I,J) = B(K,L,I,J)
    1 CONTINUE
      RETURN

      END
