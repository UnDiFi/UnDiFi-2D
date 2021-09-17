      SUBROUTINE   GPSKCO   (N, KEY, ERROR)                             GPSK2436
C
C     ==================================================================
C
C     I N S E R T I O N    S O R T
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES
C                 WHICH ARE TO BE SORTED
C
C     OUTPUT:
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN DECREASING
C                 ORDER
C
C     ==================================================================
C
      INTEGER     N, ERROR
C
CIBM  INTEGER *2  KEY(N)
      INTEGER     KEY(N)
C
C     ------------------------------------------------------------------
C
      INTEGER     I, J, K, IP1, JM1
C
C     ------------------------------------------------------------------
C
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     IF ( KEY (I) .GE. KEY (IP1) )  GO TO 2800
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              K = KEY (I)
              J = IP1
              JM1 = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR K FOUND'
C
 2600             KEY (JM1) = KEY (J)
                  JM1 = J
                  J = J + 1
                  IF  (J .GT. N)  GO TO 2700
                  IF (KEY (J) .GT. K)  GO TO 2600
C
 2700         KEY (JM1) = K
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
