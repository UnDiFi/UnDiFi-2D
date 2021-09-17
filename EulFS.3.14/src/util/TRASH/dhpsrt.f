C**********************************************************************
C
      SUBROUTINE DHPSRT (N,X,Y,Z)
C
C     Purpose:
C     This subroutine sorts the  array X using HeapSort.  It rearranges
C     the elements of Y and Z at the same time.
C
C     Parameters:
C     N = the length of the arrays (input).
C     X = the primary array to be used in sorting (input/output).
C     Y = another array, sorted in the same order as X (input/output).
C     Z = another array, sorted in the same order as X (input/output).
C
C     Noel M. Nachtigal
C     October 4, 1990
C
C**********************************************************************
C
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), Z(N)
C
C     Local variables.
C
      INTEGER I, J, K, L
      DOUBLE PRECISION TMPX, TMPY, TMPZ
C
      IF (N.LE.1) RETURN
C
      L = N / 2 + 1
      K = N
 10   IF (L.GT.1) THEN
         L = L - 1
         TMPX = X(L)
         TMPY = Y(L)
         TMPZ = Z(L)
      ELSE
         TMPX = X(K)
         TMPY = Y(K)
         TMPZ = Z(K)
         X(K) = X(1)
         Y(K) = Y(1)
         Z(K) = Z(1)
         K = K - 1
         IF (K.LE.1) THEN
            X(1) = TMPX
            Y(1) = TMPY
            Z(1) = TMPZ
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN
         IF (J.LT.K) THEN
            IF (X(J).GT.X(J+1)) J = J + 1
         END IF
         IF (TMPX.GT.X(J)) THEN
            X(I) = X(J)
            Y(I) = Y(J)
            Z(I) = Z(J)
            I    = J
            J    = J + J
         ELSE
            J = K + 1
         END IF
         GO TO 20
      END IF
      X(I) = TMPX
      Y(I) = TMPY
      Z(I) = TMPZ
      GO TO 10
      END
C
C**********************************************************************
