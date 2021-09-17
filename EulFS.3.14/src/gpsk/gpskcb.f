      SUBROUTINE  GPSKCB  (N, DEGREE, RSTART, CONNEC, AVAIL, NLEFT,     GPSKC446
     1                     STNODE, RVNODE, WORK, FORWD, BESTBK, NNODES,
     2                     DEPTH, FWIDTH, BWIDTH, ERROR, SPACE)
C
C     ==================================================================
C
C     FIND A PSEUDO-DIAMETER OF THE MATRIX GRAPH ...
C
C         << BUILD A LEVEL TREE FROM STNODE >>
C         REPEAT
C             << BUILD A LEVEL TREE FROM EACH NODE 'BKNODE' IN THE
C                DEEPEST LEVEL OF  STNODE'S TREE >>
C             << REPLACE 'STNODE' WITH 'BKNODE' IF A DEEPER AND
C                NARROWER TREE WAS FOUND. >>
C         UNTIL
C             << NO FURTHER IMPROVEMENT MADE >>
C
C     ... HEURISTIC ABOVE DIFFERS FROM THE ALGORITHM PUBLISHED IN
C         SIAM J. NUMERICAL ANALYSIS, BUT MATCHES THE CODE
C         DISTRIBUTED BY TOMS.
C
C
C     PARAMETERS :
C
C         N, DEGREE, RSTART & CONNEC  DESCRIBE THE MATRIX STRUCTURE
C
C         WORK   -- WORKING SPACE, OF LENGTH  3*AVAIL, USED TO STORE
C         THREE LEVEL TREES.
C
C         STNODE IS INITIALLY THE NUMBER OF A NODE TO BE USED TO
C             START THE PROCESS, TO BE THE ROOT OF THE FIRST TREE.
C             ON OUTPUT, STNODE IS THE END OF THE PSEUDO-DIAMETER WHOSE
C             LEVEL TREE IS NARROWEST.
C
C         RVNODE WILL BE THE OTHER END OF THE PSEUDO-DIAMETER.
C
C         NNODES WILL BE THE NUMBER OF NODES IN THIS CONNECTED
C             COMPONNENT OF THE MATRIX GRAPH, I.E., THE LENGTH OF
C             THE LEVEL TREES.
C
C         DEPTH  -- THE DEPTH OF THE LEVEL TREES BEING RETURNED,
C                   I.E., THE LENGTH OF THE PSEUDO-DIAMETER.
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ...
C
C     ---------------------------------------------------------------
C     : NUMBERED :  TLIST1  PTR1  :  TLIST2  PTR2  :  TLIST3  PTR3  :
C     ---------------------------------------------------------------
C
C     TLISTI IS A LIST OF NODES OF LENGTH  'ACTIVE'
C     PTRI   IS A LIST OF POINTERS INTO TLISTI, OF LENGTH  'DEPTH+1'
C
C     ==================================================================
C
      INTEGER     N, RSTART(N), AVAIL, NLEFT,
     1            STNODE, RVNODE, FORWD, BESTBK, NNODES, DEPTH, FWIDTH,
     4            BWIDTH, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), WORK(AVAIL,3)
      INTEGER     DEGREE(N), CONNEC(1), WORK(AVAIL,3)
C
C     ----------------
C
      INTEGER     BACKWD, MXDPTH, WIDTH, FDEPTH, LSTLVL,
     1            NLAST, T, I, BKNODE, LSTLVI
C
      LOGICAL     IMPROV
C
C
C     ... BUILD INITIAL LEVEL TREE FROM 'STNODE'.  FIND OUT HOW MANY
C         NODES LIE IN THE CURRENT CONNECTED COMPONENT.
C
      FORWD = 1
      BACKWD = 2
      BESTBK = 3
C
      CALL GPSKCC (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL, NLEFT,
     1              WORK(1,FORWD), NNODES, DEPTH, WIDTH, ERROR,
     2              SPACE)
      IF ( ERROR .NE. 0 )  GO TO 5000
C
      MXDPTH = AVAIL - NNODES - 1
C
C     ==========================================
C     REPEAT UNTIL NO DEEPER TREES ARE FOUND ...
C     ==========================================
C
 1000     FWIDTH = WIDTH
          FDEPTH = DEPTH
          LSTLVL = AVAIL - DEPTH + 1
          NLAST = WORK (LSTLVL-1, FORWD) - WORK (LSTLVL, FORWD)
          LSTLVL = WORK (LSTLVL, FORWD)
          BWIDTH = N+1
C
C         ... SORT THE DEEPEST LEVEL OF 'FORWD' TREE INTO INCREASING
C             ORDER OF NODE DEGREE.
C
          CALL GPSKCQ (NLAST, WORK(LSTLVL,FORWD), N, DEGREE, ERROR)
          IF  (ERROR .NE. 0)  GO TO 6000
C
C         ... BUILD LEVEL TREE FROM NODES IN 'LSTLVL' UNTIL A DEEPER
C             AND NARROWER TREE IS FOUND OR THE LIST IS EXHAUSTED.
C
          IMPROV = .FALSE.
          DO 1200 I = 1, NLAST
              LSTLVI = LSTLVL + I - 1
              BKNODE = WORK (LSTLVI, FORWD)
              CALL GPSKCD (N, DEGREE, RSTART, CONNEC, BKNODE, AVAIL,
     1                     NNODES, MXDPTH, WORK(1,BACKWD), DEPTH, WIDTH,
     2                     BWIDTH, ERROR, SPACE)
              IF ( ERROR .NE. 0 )  GO TO 5000
C
              IF ( DEPTH .LE. FDEPTH )  GO TO 1100
C
C                 ... NEW DEEPER TREE ... MAKE IT NEW 'FORWD' TREE
C                     AND BREAK OUT OF 'DO' LOOP.
C
                  IMPROV = .TRUE.
                  T = FORWD
                  FORWD = BACKWD
                  BACKWD = T
                  STNODE = BKNODE
                  GO TO 1300
C
C                 ... ELSE CHECK FOR NARROWER TREE.
C
 1100             IF ( WIDTH .GE. BWIDTH )  GO TO 1200
                      T = BESTBK
                      BESTBK = BACKWD
                      BACKWD = T
                      BWIDTH = WIDTH
                      RVNODE = BKNODE
 1200     CONTINUE
C
C         ... END OF REPEAT LOOP
C         ----------------------
C
 1300     IF ( IMPROV )  GO TO 1000
C
      DEPTH = FDEPTH
      RETURN
C
C     ... IN CASE OF ERROR, SIMPLY RETURN ERROR FLAG TO USER.
C
 5000 RETURN
C
 6000 ERROR = 11
      SPACE = -1
      RETURN
C
      END
