      SUBROUTINE   GPSKCJ   (N, DEGREE, RSTART, CONNEC,                 GPSK1510
     1                       NCOMPN, INVNUM, SNODE1, SNODE2, REVRS1,
     2                       DEPTH, LVLLST, LVLPTR, LVLNUM, ERROR,
     3                       SPACE)
C
C     ==================================================================
C
C     NUMBER THE NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING
C     TO A GENERALIZATION OF THE CUTHILL MCKEE STRATEGY.
C
C     N      -- DIMENSION OF ORIGINAL PROBLEM
C     DEGREE, RSTART, CONNEC -- GIVE STRUCTURE OF SPARSE AND
C                               SYMMETRIC MATRIX
C
C     NCOMPN -- NUMBER OF NODES IN THIS COMPONENT OF MATRIX GRAPH
C
C     INVNUM -- WILL BECOME A LIST OF THE ORIGINAL NODES IN THE ORDER
C               WHICH REDUCES THE BANDWIDTH OF THE MATRIX.
C
C     NXTNUM -- THE NEXT INDEX TO BE ASSIGNED (1 FOR FIRST COMPONENT)
C
C     REVRS1 -- IF .TRUE., FIRST COMPONENT OF REDUCED GRAPH WAS NUMBERED
C               BACKWARDS.
C
C     LVLLST -- LIST OF NODES IN LEVEL TREE ORDERED BY LEVEL.
C
C     LVLPTR -- POSITION OF INITIAL NODE IN EACH LEVEL OF LVLLST.
C
C     LVLNUM -- LEVEL NUMBER OF EACH NODE IN COMPONENT
C
C
      INTEGER     N, RSTART(N), NCOMPN, SNODE1, SNODE2, DEPTH,
     1            ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), INVNUM(NCOMPN),
      INTEGER     DEGREE(N), CONNEC(1), INVNUM(NCOMPN),
     1            LVLLST(NCOMPN), LVLPTR(DEPTH), LVLNUM(N)
C
      LOGICAL     REVRS1
C
C
C     ==================================================================
C
C     NUMBERING REQUIRES TWO QUEUES, WHICH CAN BE BUILD IN PLACE
C     IN INVNUM.
C
C
C     ==================================================================
C     A L G O R I T H M    S T R U C T U R E
C     ==================================================================
C
C     << SET QUEUE1 TO BE THE SET CONTAINING ONLY THE START NODE. >>
C
C     FOR LEVEL = 1 TO DEPTH DO
C
C         BEGIN
C         LOOP
C
C             REPEAT
C                 BEGIN
C                 << CNODE <- FRONT OF QUEUE1                        >>
C                 << ADD UNNUMBERED NEIGHBORS OF CNODE TO THE BACK   >>
C                 << OF QUEUE1 OR QUEUE2 (USE QUEUE1 IF NEIGHBOR     >>
C                 << AT SAME LEVEL, QUEUE2 IF AT NEXT LEVEL).  SORT  >>
C                 << THE NEWLY QUEUED NODES INTO INCREASING ORDER OF >>
C                 << DEGREE.  NUMBER CNODE, DELETE IT FROM QUEUE1.   >>
C                 END
C             UNTIL
C                 << QUEUE1 IS EMPTY >>
C
C         EXIT IF << ALL NODES AT THIS LEVEL NUMBERED >>
C
C             BEGIN
C             << FIND THE UNNUMBERED NODE OF MINIMAL DEGREE AT THIS >>
C             << LEVEL, RESTART QUEUE1 WITH THIS NODE.              >>
C             END
C
C         END << LOOP LOOP >>
C
C         << PROMOTE QUEUE2 TO BE INITIAL QUEUE1 FOR NEXT ITERATION >>
C         << OF  FOR  LOOP.                                         >>
C
C         END <<FOR LOOP>>
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C     --------------------------------------------------------------
C     : NUMBERED :  QUEUE1  :  QUEUE2  : ... : TLIST : TPTR : TREE :
C     --------------------------------------------------------------
C
C     ON COMPLETION, WE HAVE ONLY A NEW, LONGER NUMBERED SET.
C
C     ==================================================================
      INTEGER     I, BQ1, BQ2, FQ1, INC, CPTR, CNODE,
     1            INODE, LEVEL, NLEFT, LSTART, LWIDTH, QUEUE1,
     2            QUEUE2, CDGREE, XLEVEL, STNODE, ILEVEL, SQ1, SQ2,
     3            NSORT, LOWDG, BPTR, LVLLSC, LVLLSB, INVNMI
C
      LOGICAL     FORWRD, RLEVEL
C
C     ------------------------------------------------------------------
C
C     ... GIBBS-POOLE-STOCKMEYER HEURISTIC CHOICE OF ORDER
C
      IF  (DEGREE(SNODE1) .GT. DEGREE(SNODE2))  GO TO 10
          FORWRD = REVRS1
          STNODE = SNODE1
          GO TO 20
C
   10     FORWRD = .NOT. REVRS1
          STNODE = SNODE2
C
C     ... SET UP INITIAL QUEUES AT FRONT OF 'INVNUM' FOR FORWRD ORDER,
C         AT BACK FOR REVERSED ORDER.
C
   20 IF (FORWRD) GO TO 100
          INC = -1
          QUEUE1 = NCOMPN
          GO TO 200
C
  100     INC = +1
          QUEUE1 = 1
C
  200 INVNUM (QUEUE1) = STNODE
      RLEVEL = (LVLNUM(STNODE) .EQ. DEPTH)
      LVLNUM (STNODE) = 0
      FQ1 = QUEUE1
      BQ1 = QUEUE1 + INC
C
C     -------------------------------
C     NUMBER NODES LEVEL BY LEVEL ...
C     -------------------------------
C
      DO 3000 XLEVEL = 1, DEPTH
          LEVEL = XLEVEL
          IF  (RLEVEL)  LEVEL = DEPTH - XLEVEL + 1
C
          LSTART = LVLPTR (LEVEL)
          LWIDTH = LVLPTR (LEVEL+1) - LSTART
          NLEFT = LWIDTH
          QUEUE2 = QUEUE1 + INC*LWIDTH
          BQ2 = QUEUE2
C
C         ==============================================================
C         ... 'LOOP' CONSTRUCT BEGINS AT STATEMENT 1000
C                 THE INNER 'REPEAT' WILL BE DONE AS MANY TIMES AS
C                 IS NECESSARY TO NUMBER ALL THE NODES AT THIS LEVEL.
C         ==============================================================
C
 1000     CONTINUE
C
C             ==========================================================
C             ... REPEAT ... UNTIL QUEUE1 BECOMES EMPTY
C                 TAKE NODE FROM FRONT OF QUEUE1, FIND EACH OF ITS
C                 NEIGHBORS WHICH HAVE NOT YET BEEN NUMBERED, AND
C                 ADD THE NEIGHBORS TO QUEUE1 OR QUEUE2 ACCORDING TO
C                 THEIR LEVELS.
C             ==========================================================
C
 1100             CNODE = INVNUM (FQ1)
                  FQ1 = FQ1 + INC
                  SQ1 = BQ1
                  SQ2 = BQ2
                  NLEFT = NLEFT - 1
C
                  CPTR = RSTART (CNODE)
                  CDGREE = DEGREE (CNODE)
                  DO 1300 I = 1, CDGREE
                      INODE = CONNEC (CPTR)
                      CPTR = CPTR + 1
                      ILEVEL = LVLNUM (INODE)
                      IF (ILEVEL .EQ. 0)  GO TO 1300
                          LVLNUM (INODE) = 0
                          IF ( ILEVEL .EQ. LEVEL ) GO TO 1200
C
                              IF  (IABS(LEVEL-ILEVEL) .NE. 1) GO TO 6400
                                  INVNUM (BQ2) = INODE
                                  BQ2 = BQ2 + INC
                                  GO TO 1300
C
 1200                             INVNUM (BQ1) = INODE
                                  BQ1 = BQ1 + INC
 1300             CONTINUE
C
C                 ==================================================
C                 ... SORT THE NODES JUST ADDED TO QUEUE1 AND QUEUE2
C                     SEPARATELY INTO INCREASING ORDER OF DEGREE.
C                 ==================================================
C
                  IF  (IABS (BQ1 - SQ1) .LE. 1)  GO TO 1500
                      NSORT = IABS (BQ1 - SQ1)
                      IF  (FORWRD)  GO TO 1400
                          CALL GPSKCP (NSORT, INVNUM(BQ1+1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
                          GO TO 1500
C
 1400                     CALL GPSKCQ (NSORT, INVNUM(SQ1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
C
 1500             IF  (IABS (BQ2 - SQ2) .LE. 1)  GO TO 1700
                      NSORT = IABS (BQ2 - SQ2)
                      IF  (FORWRD)  GO TO 1600
                          CALL GPSKCP (NSORT, INVNUM(BQ2+1), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
                          GO TO 1700
C
 1600                     CALL GPSKCQ (NSORT, INVNUM(SQ2), N, DEGREE,
     1                                 ERROR)
                          IF  (ERROR .NE. 0)  GO TO 6600
C
C                     ... END OF REPEAT LOOP
C
 1700             IF  (FQ1 .NE. BQ1)  GO TO 1100
C
C         ==============================================================
C         ... QUEUE1 IS NOW EMPTY ...
C             IF THERE ARE ANY UNNUMBERED NODES LEFT AT THIS LEVEL,
C             FIND THE ONE OF MINIMAL DEGREE AND RETURN TO THE
C             REPEAT LOOP ABOVE.
C         ==============================================================
C
 2000     IF  ((BQ1 .EQ. QUEUE2) .AND. (NLEFT .EQ. 0))  GO TO 2900
C
              IF ((NLEFT .LE. 0) .OR. (NLEFT .NE. INC * (QUEUE2 - BQ1)))
     1             GO TO 6200
C
              LOWDG = N + 1
              BPTR  = N + 1
              CPTR  = LSTART - 1
              DO 2800 I = 1, NLEFT
 2600             CPTR   = CPTR + 1
                  LVLLSC = LVLLST (CPTR)
                  IF (LVLNUM (LVLLSC) .EQ. LEVEL)  GO TO 2700
                      IF (LVLNUM (LVLLSC) .NE. 0)  GO TO 6300
                      GO TO 2600
C
 2700             IF  (DEGREE(LVLLSC) .GE. LOWDG)  GO TO 2800
                      LOWDG = DEGREE (LVLLSC)
                      BPTR  = CPTR
C
 2800         CONTINUE
C
C             ... MINIMAL DEGREE UNNUMBERED NODE FOUND ...
C
              IF  (BPTR .GT. N)  GO TO 6500
              LVLLSB = LVLLST (BPTR)
              INVNUM (BQ1) = LVLLSB
              LVLNUM (LVLLSB) = 0
              BQ1 = BQ1 + INC
              GO TO 1000
C
C             =============================================
C             ... ADVANCE QUEUE POINTERS TO MAKE QUEUE2 THE
C                 NEW QUEUE1 FOR THE NEXT ITERATION.
C             =============================================
C
 2900     QUEUE1 = QUEUE2
          FQ1 = QUEUE1
          BQ1 = BQ2
          IF  ((BQ1 .EQ. FQ1) .AND. (XLEVEL .LT. DEPTH))  GO TO 6100
C
 3000 CONTINUE
C
C     ... CHANGE SIGN OF DEGREE TO MARK THESE NODES AS 'NUMBERED'
C
      DO 3100 I = 1, NCOMPN
          INVNMI = INVNUM(I)
          DEGREE (INVNMI) = -DEGREE (INVNMI)
 3100 CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C
 6000 SPACE = -1
      RETURN
C
 6100 ERROR = 51
      GO TO 6000
C
 6200 ERROR = 52
      GO TO 6000
C
 6300 ERROR = 53
      GO TO 6000
C
 6400 ERROR = 54
      GO TO 6000
C
 6500 ERROR = 55
      GO TO 6000
C
 6600 ERROR = 56
      GO TO 6000
C
      END
