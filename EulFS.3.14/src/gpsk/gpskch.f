      SUBROUTINE   GPSKCH   (N, DEGREE, RSTART, CONNEC, STATUS, NREDUC, GPSK1289
     1                       WORK, MXCOMP, START, SIZE, COMPNS, ERROR,
     2                       SPACE)
C
C     ==================================================================
C
C     FIND THE CONNECTED COMPONENTS OF THE GRAPH INDUCED BY THE SET
C     OF NODES WITH POSITIVE 'STATUS'.  WE SHALL BUILD THE LIST OF
C     CONNECTED COMPONENTS IN 'WORK', WITH A LIST OF POINTERS
C     TO THE BEGINNING NODES OF COMPONENTS LOCATED IN 'START'
C
C
      INTEGER     N, RSTART(N), NREDUC, MXCOMP, COMPNS, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), STATUS(N), WORK(NREDUC),
      INTEGER     DEGREE(N), CONNEC(1), STATUS(N), WORK(NREDUC),
     1            START(MXCOMP), SIZE(MXCOMP)
C
C
C     PARAMETERS ...
C
C         N      -- DIMENSION OF THE ORIGINAL MATRIX
C         DEGREE, RSTART, CONNEC -- THE STRUCTURE OF THE ORIGINAL MATRIX
C
C         STATUS -- DERIVED FROM A LEVEL TREE. POSITIVE ENTRIES INDICATE
C                   ACTIVE NODES.  NODES WITH STATUS <= 0 ARE IGNORED.
C
C         NREDUC -- THE NUMBER OF ACTIVE NODES
C
C         WORK   -- WORK SPACE, USED AS A QUEUE TO BUILD CONNECTED
C                   COMPONENTS IN PLACE.
C
C         MXCOMP -- MAXIMUM NUMBER OF COMPONENTS ALLOWED BY CURRENT
C                   SPACE ALLOCATION.  MUST NOT BE VIOLATED.
C
C         START  -- POINTER TO BEGINNING OF  I-TH  CONNECTED COMPONENT
C
C         SIZE   -- SIZE OF EACH COMPONENT
C
C         COMPNS -- NUMBER OF COMPONENTS ACTUALLY FOUND
C
C         ERROR  -- SHOULD BE ZERO ON RETURN UNLESS WE HAVE TOO LITTLE
C                   SPACE OR WE ENCOUNTER AN ERROR IN THE DATA STRUCTURE
C
C         SPACE  -- MAXIMUM AMOUNT OF WORKSPACE USED / NEEDED
C
C     ==================================================================
C
      INTEGER     I, J, FREE, JPTR, NODE, JNODE, FRONT, CDGREE, ROOT
C
C     ------------------------------------------------------------------
C
C
C     REPEAT
C         << FIND AN UNASSIGNED NODE AND START A NEW COMPONENT >>
C         REPEAT
C             << ADD ALL NEW NEIGHBORS OF FRONT NODE TO QUEUE, >>
C             << REMOVE FRONT NODE.                            >>
C         UNTIL <<QUEUE EMPTY>>
C     UNTIL << ALL NODES ASSIGNED >>
C
      FREE   = 1
      COMPNS = 0
      ROOT   = 1
C
C     ... START OF OUTER REPEAT LOOP
C
C         ... FIND AN UNASSIGNED NODE
C
  100     DO 200 I = ROOT, N
              IF (STATUS(I) .LE. 0) GO TO 200
                  NODE = I
                  GO TO 300
  200     CONTINUE
          GO TO 6100
C
C         ... START NEW COMPONENT
C
  300     COMPNS = COMPNS + 1
          ROOT   = NODE + 1
          IF (COMPNS .GT. MXCOMP)  GO TO 5000
          START (COMPNS) = FREE
          WORK (FREE) = NODE
          STATUS (NODE) = -STATUS (NODE)
          FRONT = FREE
          FREE = FREE + 1
C
C             ... INNER REPEAT UNTIL QUEUE BECOMES EMPTY
C
  400         NODE = WORK (FRONT)
              FRONT = FRONT + 1
C
              JPTR = RSTART (NODE)
              CDGREE = DEGREE (NODE)
              DO 500 J = 1, CDGREE
                  JNODE = CONNEC (JPTR)
                  JPTR = JPTR + 1
                  IF (STATUS(JNODE) .LT. 0) GO TO 500
                  IF (STATUS(JNODE) .EQ. 0) GO TO 6000
                      STATUS (JNODE) = -STATUS (JNODE)
                      WORK (FREE) = JNODE
                      FREE = FREE + 1
  500         CONTINUE
C
              IF (FRONT .LT. FREE) GO TO 400
C
C         ... END OF INNER REPEAT.  COMPUTE SIZE OF COMPONENT AND
C             SEE IF THERE ARE MORE NODES TO BE ASSIGNED
C
          SIZE (COMPNS) = FREE - START (COMPNS)
          IF (FREE .LE. NREDUC)  GO TO 100
C
      IF (FREE .NE. NREDUC+1)  GO TO 6200
      RETURN
C
C     ------------------------------------------------------------------
C
 5000 SPACE = NREDUC - FREE + 1
      ERROR = 130
      RETURN
C
 6000 ERROR = 33
      SPACE = -1
      RETURN
C
 6100 ERROR = 34
      SPACE = -1
      RETURN
C
 6200 ERROR = 35
      SPACE = -1
      RETURN
      END
