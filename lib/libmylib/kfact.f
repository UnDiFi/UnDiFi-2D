C
C      ________________________________________________________
C     |                                                        |
C     |     FACTOR A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --ARRAY CONTAINING MATRIX                |
C     |                 (LENGTH AT LEAST 2 + N(N+2))           |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --MATRIX DIMENSION                       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --FACTORED MATRIX                        |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |    PACKAGE SUBROUTINES: PACK                           |
C     |________________________________________________________|
C
      SUBROUTINE KFACT(A,LA,N)
C     .. Scalar Arguments ..
      INTEGER LA,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R,S,T
      INTEGER B,C,D,E,F,G,H,I,J,K,L,M,O,P,Q
C     ..
C     .. External Subroutines ..
      EXTERNAL PACK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      IF (N.LT.LA) CALL PACK(A,LA,N)
      R = 0.D0
      O = N + 1
      P = O + 1
      L = 5 + N*P
      I = -N - 3
C     ---------------------------------------------
C     |*** INSERT PIVOT ROW AND COMPUTE 1-NORM ***|
C     ---------------------------------------------
   10 L = L - O
      IF (L.EQ.4) GOTO 30
      S = 0.D0
      DO 20 K = 1,N
          J = L - K
          T = A(I+J)
          A(J) = T
   20 S = S + ABS(T)
      IF (R.LT.S) R = S
      I = I + 1
      GOTO 10

   30 A(1) = 1236
      A(2) = N
      A(3) = R
      Q = 3 + N*O
      I = 5 - P
      B = 0
      K = 1
   40 I = I + P
      IF (K.EQ.N) GOTO 120
      E = N - K
      H = I
      DO 50 M = I,Q,O
          L = I + E
C     --------------------
C     |*** FIND PIVOT ***|
C     --------------------
          DO 50 J = M,L
   50 IF (ABS(A(J)).GT.ABS(A(H))) H = J
      C = (H-4)/O
      D = 4 + O*C + K
      G = H - D
      H = D - I
      L = I + E
      F = I - B
C     -----------------------------
C     |*** INTERCHANGE COLUMNS ***|
C     -----------------------------
      DO 60 J = F,L
          T = A(J)
          M = J + H
          A(J) = A(M)
   60 A(M) = T
      J = I - K
      A(J) = G + K
      H = G + I
      T = A(H)
      A(H) = A(I)
      A(I) = T
      B = K
      K = K + 1
      IF (T.EQ.0.D0) GOTO 120
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      M = I + 1
      DO 70 J = M,L
   70 A(J) = A(J)/T
      F = I + E*O
   80 J = K + L
      H = J + G
      T = A(H)
      A(H) = A(J)
      A(J) = T
      L = E + J
      IF (T.EQ.0.D0) GOTO 100
      H = I - J
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      M = J + 1
      DO 90 J = M,L
   90 A(J) = A(J) - T*A(J+H)
  100 IF (L.LT.F) GOTO 80
      A(L+B) = C + 1
      GOTO 40

  110 A(1) = -1236
      RETURN

  120 IF (A(I).EQ.0) GOTO 110
      RETURN

      END
C
C      ________________________________________________________
C     |                                                        |
C     |   REARRANGE THE ELEMENTS OF A REAL ARRAY SO THAT THE   |
C     |  ELEMENTS OF A SQUARE MATRIX ARE STORED SEQUENTIALLY   |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --REAL ARRAY CONTAINING SQUARE MATRIX    |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --MATRIX PACKED AT START OF ARRAY        |
C     |________________________________________________________|
C
      SUBROUTINE PACK(A,LA,N)
C     .. Scalar Arguments ..
      INTEGER LA,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(1)
C     ..
C     .. Local Scalars ..
      INTEGER H,I,J,K,L,O
C     ..
      H = LA - N
      IF (H.EQ.0) RETURN
      IF (H.GT.0) GOTO 10
      WRITE (6,FMT=*)
     +  'ERROR: LA ARGUMENT IN PACK MUST BE .GE. N ARGUMENT'
      STOP

   10 I = 0
      K = 1
      L = N
      O = N*N
   20 IF (L.EQ.O) RETURN
      I = I + H
      K = K + N
      L = L + N
      DO 30 J = K,L
   30 A(J) = A(I+J)
      GOTO 20

      END
