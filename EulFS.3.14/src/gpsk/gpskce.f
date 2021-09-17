      SUBROUTINE   GPSKCE   (N, AVAIL, ACTIVE, DEPTH, WRKLEN,           GPSKC855
     1                       LVLLST, LVLPTR, WORK, NXTNUM, TREE1, TREE2,
     2                       WIDTH1, WIDTH2, ONEIS1, ERROR, SPACE)
C
C     ==================================================================
C
C     TRANSITION BETWEEN ALGORITHM I AND ALGORITHM II OF
C     THE GIBBS-POOLE-STOCKMEYER PAPER.
C
C     IN THIS IMPLEMENTATION ALGORITHM I REPRESENTS LEVEL TREES AS
C     LISTS OF NODES ORDERED BY LEVEL.  ALGORITHM II APPEARS TO REQUIRE
C     LEVEL NUMBERS INDEXED BY NODE -- VECTORS FOR EFFICIENCY.
C     THIS SUBROUTINE CHANGES THE LEVEL TREE REPRESENTATION TO THAT
C     REQUIRED BY ALGORITHM II.  NOTE THAT THE FIRST ALGORITHM CAN BE
C     CARRIED OUT WITH THE LEVEL NUMBER VECTOR FORMAT, PROBABLY REQURING
C     MORE COMPUTATION TIME, BUT PERHAPS LESS STORAGE.
C
C     INPUT:  TWO LEVEL TREES, AS LEVEL LISTS AND LEVEL POINTERS,
C             FOUND IN TWO OF THE THREE COLUMNS OF THE ARRAYS 'LVLLST'
C             AND 'LVLPTR'
C
C     OUTPUT: TWO LEVEL TREES, AS VECTORS OF LEVEL NUMBERS,
C             ONE PACKED TO THE FRONT, ONE TO THE REAR OF THE WORKING
C             AREA 'WORK'.  NOTE THAT 'WORK', 'LVLLST' AND 'LVLPTR'
C             SHARE COMMON LOCATIONS.
C
C     ================================================================
C
C     ... STRUCTURE OF WORKSPACE
C
C         INPUT .. (OUTPUT FROM GPSKCB)
C
C     --------------------------------------------------------------
C     : NUMBERED : TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
C     --------------------------------------------------------------
C
C         OUTPUT .. (GOES TO COMBIN)
C
C     --------------------------------------------------------------
C     : NUMBERED :  TREE2  :           ...               :  TREE1  :
C     --------------------------------------------------------------
C
C     ==================================================================
C
      INTEGER     N, AVAIL, ACTIVE, DEPTH, WRKLEN, NXTNUM,
     1            WIDTH1, WIDTH2, TREE1, TREE2, ERROR, SPACE
C
CIBM  INTEGER *2  LVLLST(AVAIL,3), LVLPTR(AVAIL,3), WORK(WRKLEN)
      INTEGER     LVLLST(AVAIL,3), LVLPTR(AVAIL,3), WORK(WRKLEN)
C
      LOGICAL     ONEIS1
C
C     ------------------------------------------------------------------
C
      INTEGER     I, BTREE, FTREE, FWIDTH, BWIDTH
C
C
C     ... CHECK THAT WE HAVE ENOUGH ROOM TO DO THE NECESSARY UNPACKING
C
      IF (3*AVAIL .GT. WRKLEN)  GO TO 6000
      IF (AVAIL .LT. N)  GO TO 5100
C
C     ... INPUT HAS THREE POSSIBLE CASES:
C             LVLLST(*,1) IS EMPTY
C             LVLLST(*,2) IS EMPTY
C             LVLLST(*,3) IS EMPTY
C
      FTREE = TREE1
      BTREE = TREE2
      FWIDTH = WIDTH1
      BWIDTH = WIDTH2
C
      TREE1 = WRKLEN - N + 1
      TREE2 = NXTNUM
C
      IF ( (FTREE .EQ. 1) .OR. (BTREE .EQ. 1) )  GO TO 300
C
C         ... CASE 1:  1ST SLOT IS EMPTY.  UNPACK 3 INTO 1, 2 INTO 3
C
          IF (FTREE .NE. 2)  GO TO 100
              ONEIS1 = .TRUE.
              WIDTH2 = BWIDTH
              WIDTH1 = FWIDTH
              GO TO 200
C
  100         ONEIS1 = .FALSE.
              WIDTH1 = BWIDTH
              WIDTH2 = FWIDTH
C
  200     CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,3), LVLPTR(1,3),
     1                    WORK(TREE2), ONEIS1)
C
          CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,2), LVLPTR(1,2),
     1                    WORK(TREE1), .NOT. ONEIS1)
C
          GO TO 1000
C
C
  300 IF ( (FTREE .EQ. 2) .OR. (BTREE .EQ. 2) )  GO TO 600
C
C         ... CASE 2:  2ND SLOT IS EMPTY.  TO ENABLE COMPLETE
C              REPACKING, MOVE 3 INTO 2, THEN FALL INTO NEXT CASE
C
          DO 400 I = 1, ACTIVE
              LVLLST(I,2) = LVLLST(I,3)
  400     CONTINUE
C
          DO 500 I = 1, DEPTH
              LVLPTR(I,2) = LVLPTR(I,3)
  500     CONTINUE
C
C         ... CASE 3:  SLOT 3 IS EMPTY.  MOVE 1 INTO 3, THEN 2 INTO 1.
C
  600     IF (FTREE .EQ. 1) GO TO 700
              ONEIS1 = .FALSE.
              WIDTH1 = BWIDTH
              WIDTH2 = FWIDTH
              GO TO 800
C
  700         ONEIS1 = .TRUE.
              WIDTH1 = FWIDTH
              WIDTH2 = BWIDTH
C
  800     CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,1), LVLPTR(1,1),
     1                    WORK(TREE1), .NOT. ONEIS1)
C
          CALL GPSKCF (N, ACTIVE, DEPTH, LVLLST(1,2), LVLPTR(1,2),
     1                    WORK(TREE2), ONEIS1)
 1000 RETURN
C
C     ------------------------------------------------------------------
C
 5100 SPACE = 3 * (N - AVAIL)
      ERROR = 120
      RETURN
C
 6000 ERROR = 20
      SPACE = -1
      RETURN
C
      END
