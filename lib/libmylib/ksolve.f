C      ________________________________________________________
C     |                                                        |
C     |            SOLVE A GENERAL FACTORED SYSTEM             |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --KFACT'S OUTPUT                         |
C     |                                                        |
C     |         B     --RIGHT SIDE                             |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --SOLUTION (CAN BE IDENTIFIED WITH B     |
C     |                 ALTHOUGH THE RIGHT SIDE IS DESTROYED)  |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      SUBROUTINE KSOLVE(X,A,B)
C     .. Array Arguments ..
      DOUBLE PRECISION A(1),B(1),X(1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION T
      INTEGER H,I,J,K,L,M,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      T = A(1)
      IF (ABS(T).EQ.1236) GOTO 10
      WRITE (6,FMT=*) 'ERROR: MUST FACTOR WITH KFACT BEFORE SOLVING'
      STOP
C     -----------------------------
C     |*** FORWARD ELIMINATION ***|
C     -----------------------------
   10 N = A(2)
      H = N
      M = N + 1
      J = 4 - M
      IF (T.LT.0.D0) GOTO 100
      DO 20 I = 1,N
   20 X(I) = B(I)
      K = 1
   30 J = J + M
      IF (K.EQ.N) GOTO 50
      L = A(J)
      T = X(L)
      X(L) = X(K)
      X(K) = T
      K = K + 1
      IF (T.EQ.0.D0) GOTO 30
      DO 40 I = K,N
   40 X(I) = X(I) - T*A(I+J)
      GOTO 30
C     --------------------------------------
C     |*** BACK SUBSTITUTION BY COLUMNS ***|
C     --------------------------------------
   50 T = X(K)/A(J+K)
   60 X(K) = T
      IF (K.EQ.1) GOTO 80
      K = K - 1
      DO 70 I = 1,K
   70 X(I) = X(I) - T*A(I+J)
      J = J - M
      GOTO 50
C     -----------------------
C     |*** PIVOT SOLUTION***|
C     -----------------------
   80 L = 3 + M*N
   90 H = H - 1
      IF (H.EQ.0) RETURN
      I = A(H+L)
      T = X(H)
      X(H) = X(I)
      X(I) = T
      GOTO 90
C     -----------------------------
C     |*** COMPUTE NULL VECTOR ***|
C     -----------------------------
  100 K = 0
  110 K = K + 1
      J = J + M
      IF (A(J+K).NE.0.D0) GOTO 110
      DO 120 I = 1,N
  120 X(I) = 0.D0
      T = 1.D0
      H = K
      GOTO 60

      END
