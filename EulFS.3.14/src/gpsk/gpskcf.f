      SUBROUTINE  GPSKCF  (N, ACTIVE, DEPTH, LVLLST, LVLPTR, LVLNUM,    GPSKC996
     1                     REVERS)
C
C     ==================================================================
C
C     CONVERT LEVEL STRUCTURE REPRESENTATION FROM A LIST OF NODES
C     GROUPED BY LEVEL TO A VECTOR GIVING LEVEL NUMBER FOR EACH NODE.
C
C     LVLLST, LVLPTR -- LIST OF LISTS
C
C     LVLNUM -- OUTPUT VECTOR OF LEVEL NUMBERS
C
C     REVERS -- IF .TRUE., NUMBER LEVEL STRUCTURE FROM BACK END
C               INSTEAD OF FROM FRONT
C
C     ==================================================================
C
      INTEGER     N, ACTIVE, DEPTH
C
CIBM  INTEGER *2  LVLLST(ACTIVE), LVLPTR(DEPTH), LVLNUM(N)
      INTEGER     LVLLST(ACTIVE), LVLPTR(DEPTH), LVLNUM(N)
      LOGICAL     REVERS
C
C     ------------------------------------------------------------------
C
      INTEGER     I, LEVEL, LSTART, LEND, XLEVEL, PLSTRT, LVLLSI
C
      IF  (ACTIVE .EQ. N)  GO TO 200
C
C         ... IF NOT ALL NODES OF GRAPH ARE ACTIVE, MASK OUT THE
C             NODES WHICH ARE NOT ACTIVE
C
          DO 100 I = 1, N
              LVLNUM(I) = 0
  100     CONTINUE
C
  200 DO 400 LEVEL = 1, DEPTH
          XLEVEL = LEVEL
          PLSTRT = DEPTH - LEVEL + 1
          IF (REVERS) XLEVEL = PLSTRT
          LSTART = LVLPTR (PLSTRT)
          LEND = LVLPTR (PLSTRT - 1) - 1
C
          DO 300 I = LSTART, LEND
              LVLLSI = LVLLST(I)
              LVLNUM (LVLLSI) = XLEVEL
  300     CONTINUE
  400 CONTINUE
C
      RETURN
      END
