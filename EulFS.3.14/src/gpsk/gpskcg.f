      SUBROUTINE   GPSKCG   (N, DEGREE, RSTART, CONNEC, ACTIVE, WIDTH1, GPSK1047
     1                       WIDTH2, TREE1, TREE2, WORK, WRKLEN, DEPTH,
     2                       INC1, INC2, TOTAL, ONEIS1, REVRS1, ERROR,
     3                       SPACE)
C
C     ==================================================================
C
C     COMBINE THE TWO ROOTED LEVEL TREES INTO A SINGLE LEVEL STRUCTURE
C     WHICH MAY HAVE SMALLER WIDTH THAN EITHER OF THE TREES.  THE NEW
C     STRUCTURE IS NOT NECESSARILY A ROOTED STRUCTURE.
C
C     PARAMETERS:
C
C         N, DEGREE, RSTART, CONNEC -- GIVE THE DIMENSION AND STRUCTURE
C                                      OF THE SPARSE SYMMETRIC MATRIX
C
C         ACTIVE -- THE NUMBER OF NODES IN THIS CONNECTED COMPONENT OF
C                   THE MATRIX GRAPH
C
C         TREE1  -- ON INPUT, ONE OF THE INPUT LEVEL TREES.  ON
C                   OUTPUT, THE COMBINED LEVEL STRUCTURE
C
C         TREE2  -- THE SECOND INPUT LEVEL TREE
C
C         WIDTH1 -- THE MAXIMUM WIDTH OF A LEVEL IN TREE1
C
C         WIDTH2 -- THE MAXIMUM WIDTH OF A LEVEL IN TREE2
C
C         WORK   -- A WORKING AREA OF LENGTH 'WRKLEN'
C
C         INC1,  -- VECTORS OF LENGTH 'DEPTH'
C         INC2,
C         TOTAL
C
C         ONEIS1 -- INDICATES WHETHER TREE1 OR TREE2 REPRESENTS THE
C                   FORWARD TREE OR THE BACKWARDS TREE OF PHASE 1.
C                   USED TO MIMIC ARBITRARY TIE-BREAKING PROCEDURE OF
C                   ORIGINAL GIBBS-POOLE-STOCKMEYER CODE.
C
C         REVRS1 -- OUTPUT PARAMETER INDICATING WHETHER A BACKWARDS
C                   ORDERING WAS USED FOR THE LARGEST COMPONENT OF
C                   THE REDUCED GRAPH
C
C         ERROR  -- NON-ZERO ONLY IF FAILURE OF SPACE ALLOCATION OR
C                   DATA STRUCTURE ERROR FOUND
C
C         SPACE -- MINIMUM SPACE REQUIRED TO RERUN OR COMPLETE PHASE.
C
C     ------------------------------------------------------------------
C
      INTEGER     N, RSTART(N), ACTIVE, WIDTH1, WIDTH2, WRKLEN, DEPTH,
     2            ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), TREE1(N), TREE2(N),
      INTEGER     DEGREE(N), CONNEC(1), TREE1(N), TREE2(N),
     1            WORK(WRKLEN), INC1(DEPTH), INC2(DEPTH), TOTAL(DEPTH)
C
      LOGICAL     ONEIS1, REVRS1
C
C     ==================================================================
C
C     << REMOVE ALL NODES OF PSEUDO-DIAMETERS >>
C     << FIND CONNECTED COMPONENTS OF REDUCED GRAPH >>
C     << COMBINE LEVEL TREES, COMPONENT BY COMPONENT >>
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ...
C
C     ------------------------------------------------------------------
C     : NUMBERED : TREE2 : TOTAL : NODES : START : SIZE : INC1 : INC2 :
C     ------------------------------------------------------------------
C
C     --------
C      TREE1 :
C     --------
C
C         NUMBERED  IS THE SET OF  NUMBERED NODES (PROBABLY EMPTY)
C
C         TREE1 AND TREE1 ARE LEVEL TREES (LENGTH N)
C         TOTAL, INC1 AND INC2  ARE VECTORS OF NODE COUNTS PER LEVEL
C             (LENGTH 'DEPTH')
C         NODES IS THE SET OF NODES IN THE REDUCED GRAPH (THE NODES
C             NOT ON ANY SHORTEST PATH FROM ONE END OF THE
C             PSEUDODIAMETER TO THE OTHER)
C         START, SIZE ARE POINTERS INTO 'NODES', ONE OF EACH FOR
C         EACH CONNECTED COMPONENT OF THE REDUCED GRAPH.
C         THE SIZES OF NODES, START AND SIZE ARE NOT KNOWN APRIORI.
C
C     ==================================================================
      INTEGER     I, SIZE, AVAIL, CSTOP, START, COMPON, TREE1I, PCSTRT,
     1            CSTART, MXINC1, MXINC2, COMPNS, MXCOMP, OFFDIA,
     2            CSIZE, PCSIZE, WORKI, TWORKI
C
C     ------------------------------------------------------------------
C
C     ... FIND ALL SHORTEST PATHS FROM START TO FINISH.  REMOVE NODES ON
C         THESE PATHS AND IN OTHER CONNECTED COMPONENTS OF FULL GRAPH
C         FROM FURTHER CONSIDERATION.  SIGN OF ENTRIES IN TREE1 IS USED
C         AS A MASK.
C
      OFFDIA = ACTIVE
C
      DO 100 I = 1, DEPTH
          TOTAL(I) = 0
  100 CONTINUE
C
      DO 200 I = 1, N
          TREE1I = TREE1 (I)
          IF ((TREE1(I) .NE. TREE2(I)) .OR. (TREE1(I) .EQ. 0)) GO TO 200
              TOTAL (TREE1I) = TOTAL (TREE1I) + 1
              TREE1(I) = - TREE1(I)
              OFFDIA = OFFDIA - 1
  200 CONTINUE
C
      IF ( OFFDIA .EQ. 0 )  GO TO 1100
      IF ( OFFDIA .LT. 0 )  GO TO 6000
C
C     ... FIND CONNECTED COMPONENTS OF GRAPH INDUCED BY THE NODES NOT
C         REMOVED.  'MXCOMP' IS THE LARGEST NUMBER OF COMPONENTS
C         REPRESENTABLE IN THE WORKING SPACE AVAILABLE.
C
      AVAIL = WRKLEN - OFFDIA
      MXCOMP = AVAIL/2
      START = OFFDIA + 1
      SIZE = START + MXCOMP
C
      IF  (MXCOMP .LE. 0)  GO TO 5100
C
      CALL GPSKCH (N, DEGREE, RSTART, CONNEC, TREE1, OFFDIA, WORK,
     1             MXCOMP, WORK(START), WORK(SIZE), COMPNS, ERROR,
     2             SPACE)
      IF ( ERROR .NE. 0 )  GO TO 5000
C
C     ... RECORD SPACE ACTUALLY USED  (NOT INCLUDING  NUMBERED )
C
      SPACE = 2*N + 3*(DEPTH) + 2*COMPNS + OFFDIA
C
C     ... SORT THE COMPONENT START POINTERS INTO INCREASING ORDER
C         OF SIZE OF COMPONENT
C
      IF (COMPNS .GT. 1)
     1    CALL GPSKCN (COMPNS, WORK(SIZE), WORK(START), ERROR)
          IF  (ERROR .NE. 0)  GO TO 6200
C
C     ... FOR EACH COMPONENT IN TURN, CHOOSE TO USE THE ORDERING OF THE
C         'FORWARD' TREE1 OR OF THE 'BACKWARD' TREE2 TO NUMBER THE NODES
C         IN THIS COMPONENT.  THE NUMBERING IS CHOSEN TO MINIMIZE THE
C         MAXIMUM INCREMENT TO ANY LEVEL.
C
      DO 1000 COMPON = 1, COMPNS
          PCSTRT = START + COMPON - 1
          CSTART = WORK (PCSTRT)
          PCSIZE = SIZE + COMPON - 1
          CSIZE = WORK (PCSIZE)
          CSTOP  = CSTART + CSIZE - 1
          IF ( ( CSIZE .LT. 0 ) .OR. ( CSIZE .GT. OFFDIA ) )  GO TO 6100
C
          DO 300 I = 1, DEPTH
              INC1(I) = 0
              INC2(I) = 0
  300     CONTINUE
C
          MXINC1 = 0
          MXINC2 = 0
C
          DO 400 I = CSTART, CSTOP
              WORKI = WORK(I)
              TWORKI = -TREE1 (WORKI)
              INC1 (TWORKI) = INC1 (TWORKI) + 1
              TWORKI =  TREE2 (WORKI)
              INC2 (TWORKI) = INC2 (TWORKI) + 1
  400     CONTINUE
C
C         ... BAROQUE TESTS BELOW DUPLICATE THE GIBBS-POOLE-STOCKMEYER-
C             CRANE PROGRAM, *** NOT *** THE PUBLISHED ALGORITHM.
C
          DO 500 I = 1, DEPTH
              IF ((INC1(I) .EQ. 0) .AND. (INC2(I) .EQ. 0))  GO TO 500
                  IF  (MXINC1  .LT.  TOTAL(I) + INC1(I))
     1                 MXINC1 = TOTAL(I) + INC1(I)
                  IF  (MXINC2  .LT.  TOTAL(I) + INC2(I))
     1                 MXINC2 = TOTAL(I) + INC2(I)
  500     CONTINUE
C
C         ... USE ORDERING OF NARROWER TREE UNLESS IT INCREASES
C             WIDTH MORE THAN WIDER TREE.  IN CASE OF TIE, USE TREE 2!
C
          IF ( (MXINC1 .GT. MXINC2)  .OR.
     1         ( (MXINC1 .EQ. MXINC2) .AND. ( (WIDTH1 .GT. WIDTH2) .OR.
     2                                        ( (WIDTH1 .EQ. WIDTH2)
     3                                         .AND. ONEIS1) ) ) )
     4      GO TO 700
C
              IF ( COMPON .EQ. 1 )  REVRS1 = .NOT. ONEIS1
C
              DO 600 I = 1, DEPTH
                  TOTAL(I) = TOTAL(I) + INC1(I)
  600         CONTINUE
              GO TO 1000
C
  700         IF ( COMPON .EQ. 1 )  REVRS1 = ONEIS1
              DO 800 I = CSTART, CSTOP
                  WORKI = WORK(I)
                  TREE1 (WORKI) = - TREE2 (WORKI)
  800         CONTINUE
C
              DO 900 I = 1, DEPTH
                  TOTAL(I) = TOTAL(I) + INC2(I)
  900         CONTINUE
C
 1000 CONTINUE
      GO TO 2000
C
C     ... DEFAULT WHEN THE REDUCED GRAPH IS EMPTY
C
 1100 REVRS1 = .TRUE.
      SPACE = 2*N
C
 2000 RETURN
C
C     ------------------------------------------------------------------
C
C     ERROR FOUND ...
C
 5000 SPACE = -1
      GO TO 2000
C
 5100 SPACE = 2 - AVAIL
      ERROR = 131
      GO TO 2000
C
 6000 ERROR = 30
      GO TO 5000
C
 6100 ERROR = 31
      GO TO 5000
C
 6200 ERROR = 32
      GO TO 5000
C
      END
