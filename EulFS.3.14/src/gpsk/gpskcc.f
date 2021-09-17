      SUBROUTINE   GPSKCC   (N, DEGREE, RSTART, CONNEC, STNODE, AVAIL,  GPSKC599
     1                       NLEFT, LIST, ACTIVE, DEPTH, WIDTH, ERROR,
     2                       SPACE)
C
C     ==================================================================
C     BUILD THE LEVEL TREE ROOTED AT 'STNODE' IN THE SPACE PROVIDED IN
C     LIST.  CHECK FOR OVERRUN OF SPACE ALLOCATION.
C     ==================================================================
C
      INTEGER     N, RSTART(N), STNODE, AVAIL, NLEFT,
     1            ACTIVE, DEPTH, WIDTH, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), LIST(AVAIL)
      INTEGER     DEGREE(N), CONNEC(1), LIST(AVAIL)
C
C     ... PARAMETERS:
C
C         INPUT ...
C
C             N, DEGREE, RSTART, CONNEC -- DESCRIBE THE MATRIX STRUCTURE
C
C             STNODE -- THE ROOT OF THE LEVEL TREE.
C
C             AVAIL  -- THE LENGTH OF THE WORKING SPACE AVAILABLE
C
C             NLEFT  -- THE NUMBER OF NODES YET TO BE NUMBERED
C
C             LIST   -- THE WORKING SPACE.
C
C         OUTPUT ...
C
C             ACTIVE -- THE NUMBER OF NODES IN THE COMPONENT
C
C             DEPTH  -- THE DEPTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             WIDTH  -- THE WIDTH OF THE LEVEL TREE ROOTED AT  STNODE.
C
C             ERROR  -- ZERO UNLESS STORAGE WAS INSUFFICIENT.
C
C     ------------------------------------------------------------------
C
      INTEGER         LSTART, NLEVEL, FRONT, J, NEWNOD, PTR, CDGREE,
     1                LFRONT, LISTJ
C
C     ... BUILD THE LEVEL TREE USING  LIST  AS A QUEUE AND LEAVING
C         THE NODES IN PLACE.  THIS GENERATES THE NODES ORDERED BY LEVEL
C         PUT POINTERS TO THE BEGINNING OF EACH LEVEL, BUILDING FROM
C         THE BACK OF THE WORK AREA.
C
      ACTIVE = 1
      DEPTH = 0
      WIDTH = 0
      ERROR = 0
      LSTART = 1
      FRONT = 1
      LIST (ACTIVE) = STNODE
      DEGREE (STNODE) = -DEGREE (STNODE)
      LIST (AVAIL)  = 1
      NLEVEL = AVAIL
C
C     ... REPEAT UNTIL QUEUE BECOMES EMPTY OR WE RUN OUT OF SPACE.
C     ------------------------------------------------------------
C
 1000     IF ( FRONT .LT. LSTART ) GO TO 1100
C
C         ... FIRST NODE OF LEVEL.  UPDATE POINTERS.
C
              LSTART = ACTIVE + 1
              WIDTH = MAX0 (WIDTH, LSTART - LIST(NLEVEL))
              NLEVEL = NLEVEL - 1
              DEPTH = DEPTH + 1
              IF ( NLEVEL .LE. ACTIVE )  GO TO 5000
                  LIST (NLEVEL) = LSTART
C
C         ... FIND ALL NEIGHBORS OF CURRENT NODE, ADD THEM TO QUEUE.
C
 1100     LFRONT = LIST (FRONT)
          PTR = RSTART (LFRONT)
          CDGREE = -DEGREE (LFRONT)
          IF (CDGREE .LE. 0)  GO TO 6000
          DO 1200 J = 1, CDGREE
              NEWNOD = CONNEC (PTR)
              PTR = PTR + 1
C
C             ... ADD TO QUEUE ONLY NODES NOT ALREADY IN QUEUE
C
              IF ( DEGREE(NEWNOD) .LE. 0 )  GO TO 1200
                  DEGREE (NEWNOD) = -DEGREE (NEWNOD)
                  ACTIVE = ACTIVE + 1
                  IF ( NLEVEL .LE. ACTIVE )  GO TO 5000
                  IF ( ACTIVE .GT. NLEFT  )  GO TO 6000
                      LIST (ACTIVE) = NEWNOD
 1200     CONTINUE
          FRONT = FRONT + 1
C
C         ... IS QUEUE EMPTY?
C         -------------------
C
          IF ( FRONT .LE. ACTIVE )  GO TO 1000
C
C     ... YES, THE TREE IS BUILT.  UNDO OUR MARKINGS.
C
      DO 1300 J = 1, ACTIVE
          LISTJ = LIST(J)
          DEGREE (LISTJ) = -DEGREE (LISTJ)
 1300 CONTINUE
C
      RETURN
C
C     ... INSUFFICIENT STORAGE ...
C
 5000 SPACE = 3 * ( (NLEFT+1-ACTIVE)*DEPTH / NLEFT + (NLEFT+1-ACTIVE) )
      ERROR = 110
      RETURN
C
 6000 ERROR = 12
      SPACE = -1
      RETURN
C
      END
