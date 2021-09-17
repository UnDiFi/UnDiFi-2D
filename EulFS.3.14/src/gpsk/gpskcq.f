      SUBROUTINE  GPSKCQ (N, INDEX, NVEC, DEGREE, ERROR)                GPSK2581
C
C     ==================================================================
C
C     I N S E R T I O N      S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         INDEX  -- AN ARRAY OF LENGTH  N  CONTAINING THE INDICES
C                 WHOSE DEGREES ARE TO BE SORTED
C         DEGREE -- AN  NVEC  VECTOR, GIVING THE DEGREES OF NODES
C                   WHICH ARE TO BE SORTED.
C
C     OUTPUT:
C         INDEX  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING
C                   ORDER
C         ERROR -- WILL BE ZERO UNLESS THE PROGRAM IS MALFUNCTIONING,
C                  IN WHICH CASE IT WILL BE EQUAL TO 1.
C
C     ==================================================================
C
      INTEGER     N, NVEC, ERROR
C
CIBM  INTEGER *2  INDEX(N), DEGREE(NVEC)
      INTEGER     INDEX(N), DEGREE(NVEC)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, V, INDEXI, INDXI1, INDEXJ, IP1, JM1
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ------------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE
C     ------------------------------------------------------------------
C
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     INDEXI = INDEX (I)
          INDXI1 = INDEX (IP1)
          IF ( DEGREE(INDEXI) .LE. DEGREE(INDXI1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              V = DEGREE (INDEXI)
              J = IP1
              JM1 = I
              INDEXJ = INDEX (J)
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR V FOUND'
C
 2600             INDEX (JM1) = INDEXJ
                  JM1 = J
                  J = J + 1
                  IF (J .GT. N)  GO TO 2700
                  INDEXJ = INDEX (J)
                  IF (DEGREE(INDEXJ) .LT. V)  GO TO 2600
C
 2700         INDEX (JM1) = INDEXI
C
 2800     IP1 = I
          I = I - 1
          IF ( I .GT. 0 )  GO TO 2500
C
 3000 RETURN
C
 6000 ERROR = 1
      GO TO 3000
C
      END
