      SUBROUTINE   QSORT2  (NN, KEY, ERROR)
C
C
C     ==================================================================
C
C     Q U I C K S O R T
C
C         IN THE STYLE OF THE CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C     INPUT:
C         N    -- NUMBER OF ELEMENTS TO BE SORTED (= NN)
C         KEY  -- AN ARRAY OF LENGTH  N  CONTAINING THE VALUES
C                 WHICH ARE TO BE SORTED
C
C     OUTPUT:
C         KEY  -- WILL BE ARRANGED SO THAT VALUES ARE IN INCREASING
C                 ORDER
C         ERROR -- WILL BE ZERO UNLESS YOUR INPUT FILE WAS OF TRULY
C                  ENORMOUS LENGTH, IN WHICH CASE IT WILL BE EQUAL TO 1.
C
C
C     ==================================================================
C
      INTEGER     NN, ERROR, KEY(NN)
C
C     ------------------------
C
      INTEGER     TOP, LEFT, RIGHT, I, J, TINY, V, K, IP1, JM1,
     1            LLEN, RLEN, N
      LOGICAL     DONE
      INTEGER     STKLEN, STACK(30)
C
      DATA        TINY, STKLEN / 9, 30 /
C
C     -----------------------------------
C
C     ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK S
C         PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C         LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C         PARTITIONING.
C
      N = NN
      IF (N .EQ. 1)  RETURN
      IF  (N .LE. 0)  GO TO 6000
C
      ERROR = 0
      TOP = 1
      LEFT = 1
      RIGHT = N
      DONE = (N .LE. TINY)
C
      IF(DONE) GO TO 2000
*     CALL IVFILL(STKLEN,STACK,0)
      CALL IINIT(STKLEN,0,STACK,1)
C
C     ===========================================================
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF
C     LENGTH GREATER THAN 'TINY'
C     ===========================================================
C
C     ... WHILE NOT DONE DO ...
C
  100 IF (DONE)  GO TO 2000
C
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT
C             SUBFILE, WHICH IS  KEY(LEFT), ..., KEY(RIGHT)
C
          LFRH2 = (LEFT + RIGHT)/2
          K = KEY (LFRH2)
          KEY (LFRH2) = KEY (LEFT)
          KEY (LEFT) = K
C
          IF ( KEY(LEFT+1) .LE. KEY(RIGHT) ) GO TO 200
              K = KEY (LEFT+1)
              KEY (LEFT+1) = KEY (RIGHT)
              KEY (RIGHT) = K
C
  200     IF ( KEY(LEFT) .LE. KEY(RIGHT) )  GO TO 300
              K = KEY (LEFT)
              KEY (LEFT) = KEY (RIGHT)
              KEY (RIGHT) = K
C
  300     IF ( KEY (LEFT+1) .LE. KEY (LEFT) )  GO TO 400
              K = KEY (LEFT+1)
              KEY (LEFT+1) = KEY (LEFT)
              KEY (LEFT) = K
C
  400     V = KEY (LEFT)
C
C         ... V IS NOW THE MEDIAN VALUE OF THE THREE KEYS.  NOW MOVE
C             FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING
C             KEYS UNTIL ALL KEYS LESS THAN  V  ARE PACKED TO
C             THE LEFT, ALL KEYS LARGER THAN  V  ARE PACKED TO THE
C             RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL KEY(I) >= V;
C             REPEAT J = J-1 UNTIL KEY(J) <= V;
C         EXIT IF J < I;
C             << EXCHANGE KEYS I AND J >>
C         END
C
  500     CONTINUE
  600         I  = I + 1
              IF ( KEY(I) .LT. V )  GO TO 600
C
  700         J = J - 1
              IF ( KEY(J) .GT. V )  GO TO 700
C
          IF (J .LT. I)  GO TO 800
              K = KEY (I)
              KEY (I) = KEY (J)
              KEY (J) = K
          GO TO 500
C
  800     K = KEY (LEFT)
          KEY (LEFT) = KEY (J)
          KEY (J) = K
C
C
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.
C
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF ( MAX0 (LLEN, RLEN) .GT. TINY )  GO TO 1100
C
C             ... BOTH SUBFILES ARE TINY, SO UNSTACK NEXT LARGER FILE
C
              IF (TOP .EQ. 1)  GO TO 900
                  TOP = TOP - 2
                  LEFT = STACK (TOP)
                  RIGHT = STACK (TOP+1)
                  GO TO 100
C
  900             DONE = .TRUE.
C
                  GO TO 100
C
C             ... ELSE ONE OR BOTH SUBFILES ARE LARGE
C
 1100     IF (MIN0 (LLEN, RLEN) .GT. TINY)  GO TO 1400
C
C             ... ONE SUBFILE IS SMALL, ONE LARGE.  IGNORE THE SMALL ONE
C
              IF ( LLEN .GT. RLEN )  GO TO 1200
                  LEFT = I
                  GO TO 100
C
 1200             RIGHT = J - 1
C
              GO TO 100
C
C         ... ELSE BOTH ARE LARGER THAN TINY.  ONE MUST BE STACKED.
C
 1400     IF  ( TOP .GE. STKLEN )  GO TO 6000
          IF ( LLEN .GT. RLEN )  GO TO 1500
              STACK (TOP) = I
              STACK (TOP+1) = RIGHT
              RIGHT = J-1
              GO TO 1600
C
 1500         STACK (TOP) = LEFT
              STACK (TOP+1) = J-1
              LEFT = I
C
 1600     TOP = TOP + 2
C
      GO TO 100
C
C     ------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.
C     ------------------------------------------------------------
C
C     ... FIRST, FIND LARGEST ELEMENT IN 'KEY'
C
 2000 I    = N - 1
      LEFT = MAX0 (0, N - TINY)
      K    = KEY (N)
      J    = N
C
 2100     IF  ( I .LE. LEFT )  GO TO 2300
              IF  ( KEY(I) .LE. K )  GO TO 2200
                  K = KEY(I)
                  J = I
C
 2200         I = I - 1
              GO TO 2100
C
 2300 IF  ( J .EQ. N )  GO TO 2400
C
C     ... LARGEST ELEMENT WILL BE IN  KEY(N)
C
          KEY(J)  = KEY(N)
          KEY(N)  = K
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2400 I = N - 1
      IP1 = N
C
 2500     IF ( KEY (I) .LE. KEY (IP1) )  GO TO 2800
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
                  IF (KEY (J) .LT. K)  GO TO 2600
C
              KEY (JM1) = K
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
