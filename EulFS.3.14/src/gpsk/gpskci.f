      SUBROUTINE   GPSKCI   (N, ACTIVE, DEPTH, LSTRUC, LVLLST, LVLPTR,  GPSK1422
     1                       LTOTAL, ERROR, SPACE)
C
C     ==================================================================
C
C     TRANSITIONAL SUBROUTINE, ALGORITHM II TO IIIA OR IIIB.
C
C     CONVERT LEVEL STRUCTURE GIVEN AS VECTOR OF LEVEL NUMBERS FOR NODES
C     TO STRUCTURE AS LIST OF NODES BY LEVEL
C
C     N, ACTIVE, DEPTH -- PROBLEM SIZES
C     LSTRUC -- INPUT LEVEL STRUCTURE
C     LVLLST, LVLPTR -- OUTPUT LEVEL STRUCTURE
C     LTOTAL -- NUMBER OF NODES AT EACH LEVEL (PRECOMPUTED)
C
      INTEGER     N, ACTIVE, DEPTH, ERROR, SPACE
C
CIBM  INTEGER *2  LSTRUC(N), LVLLST(ACTIVE), LVLPTR(1), LTOTAL(DEPTH)
      INTEGER     LSTRUC(N), LVLLST(ACTIVE), LVLPTR(1), LTOTAL(DEPTH)
C
C     ===============================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C         INPUT (FROM COMBIN) ..
C
C     ------------------------------------------------------------------
C     :  NUMBERED  :  ..(N)..  :  TOTAL  :         ...        :  TREE  :
C     ------------------------------------------------------------------
C
C         OUTPUT (TO GPSKCJ OR GPSKCK) ..
C
C     ------------------------------------------------------------------
C     :  NUMBERED  :       ...             :  TLIST  :  TPTR  :  TREE  :
C     ------------------------------------------------------------------
C
C     HERE, NUMBERED IS THE SET OF NODES IN NUMBERED COMPONENTS
C         TOTAL IS A VECTOR OF LENGTH 'DEPTH' GIVING THE NUMBER
C         OF NODES IN EACH LEVEL OF THE 'TREE'.
C         TLIST, TPTR ARE LISTS OF NODES OF THE TREE, ARRANGED
C         BY LEVEL.  TLIST IS OF LENGTH 'ACTIVE', TPTR 'DEPTH+1'.
C
C     =================================================================
C
      INTEGER     I, ACOUNT, START, LEVEL, PLEVEL
C
C     ... ESTABLISH STARTING AND ENDING POINTERS FOR EACH LEVEL
C
      START = 1
      DO 100 I = 1, DEPTH
          LVLPTR(I) = START
          START = START + LTOTAL(I)
          LTOTAL(I) = START
  100 CONTINUE
      LVLPTR(DEPTH+1) = START
C
      ACOUNT = 0
      DO 300 I = 1, N
          IF (LSTRUC(I)) 200, 300, 6000
  200         LEVEL = -LSTRUC(I)
              LSTRUC(I) = LEVEL
              PLEVEL = LVLPTR (LEVEL)
              LVLLST (PLEVEL) = I
              LVLPTR (LEVEL) = LVLPTR (LEVEL) + 1
              ACOUNT = ACOUNT + 1
              IF (LVLPTR (LEVEL) .GT. LTOTAL (LEVEL))  GO TO 6100
  300 CONTINUE
C
C     ... RESET STARTING POINTERS
C
      LVLPTR(1) = 1
      DO 400 I = 1, DEPTH
          LVLPTR(I+1) = LTOTAL(I)
  400 CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C
 6000 ERROR = 40
      GO TO 6200
C
 6100 ERROR = 41
C
 6200 SPACE = -1
      RETURN
C
      END
