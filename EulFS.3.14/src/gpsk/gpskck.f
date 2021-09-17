      SUBROUTINE  GPSKCK  (N, DEGREE, RSTART, CONNEC, WRKLEN, NXTNUM,   GPSK1811
     1                     WORK, NCOMPN, DEPTH, LVLLST, LVLPTR, LVLNUM,
     2                     ERROR, SPACE)
C
      INTEGER     N, RSTART(N), WRKLEN, NXTNUM, NCOMPN, DEPTH, ERROR,
     1            SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), WORK(WRKLEN), LVLLST(N),
      INTEGER     DEGREE(N), CONNEC(1), WORK(WRKLEN), LVLLST(N),
     1            LVLPTR(DEPTH), LVLNUM(N)
C
C     ==================================================================
C
C     NUMBER NODES IN A GENERALIZED LEVEL STRUCTURE ACCORDING TO
C     A GENERALIZATION OF THE KING ALGORITHM, WHICH REDUCES
C     THE PROFILE OF THE SPARSE SYMMETRIC MATRIX.
C
C     ---------------------
C
C     CODE USES A PRIORITY QUEUE TO CHOOSE THE NEXT NODE TO BE NUMBERED
C     THE PRIORITY QUEUE IS REPRESENTED BY A SIMPLE LINEAR-LINKED LIST
C     TO SAVE SPACE.  THIS WILL REQUIRE MORE SEARCHING THAN A FULLY
C     LINKED REPRESENTATION, BUT THE DATA MANIPULATION IS SIMPLER.
C
C     -------------------
C
C     << ESTABLISH PRIORITY QUEUE 'ACTIVE' FOR LEVEL 1 NODES >>
C
C     FOR I = 1 TO DEPTH DO
C         << SET QUEUE 'QUEUED' TO BE EMPTY, LIST 'NEXT' TO BE >>
C         << SET OF NODES AT NEXT LEVEL.                       >>
C
C         FOR J = 1 TO 'NODES AT THIS LEVEL' DO
C             << FIND FIRST NODE IN ACTIVE WITH MINIMAL CONNECTIONS >>
C             << TO 'NEXT'.  NUMBER THIS NODE AND REMOVE HIM FROM   >>
C             << 'ACTIVE'.  FOR EACH NODE IN 'NEXT' WHICH CONNECTED >>
C             << TO THIS NODE, MOVE IT TO 'QUEUED' AND REMOVE IT    >>
C             << FROM 'NEXT'.                                       >>
C
C         << SET NEW QUEUE 'ACTIVE' TO BE 'QUEUED' FOLLOWED BY ANY >>
C         << NODES STILL IN 'NEXT'.                                >>
C
C     ==================================================================
C
C     DATA STRUCTURE ASSUMPTIONS:
C     THE FIRST 'NXTNUM-1' ELEMENTS OF  WORK  ARE ALREADY IN USE.
C     THE LEVEL STRUCTURE 'LVLLST' IS CONTIGUOUS WITH  WORK, THAT IS,
C     IT RESIDES IN ELEMENTS  WRKLEN+1, ...  OF  WORK.  'LVLPTR' AND
C     'LVLNUM' ARE ALSO EMBEDDED IN WORK, BEHIND 'LVLLST'.  THE
C     THREE VECTORS ARE PASSED SEPARATELY TO CLARIFY THE INDEXING,
C     BUT THE QUEUES DEVELOPED WILL BE ALLOWED TO OVERRUN 'LVLLST'
C     AS NEEDED.
C
C     ... BUILD THE FIRST 'ACTIVE' QUEUE STARTING W1 LOCATIONS FROM
C         THE FRONT OF THE CURRENT WORKING AREA  (W1 IS THE WIDTH OF THE
C         FIRST LEVEL).  BUILD THE FIRST 'QUEUED' QUEUE STARTING FROM
C         THE BACK OF WORK SPACE.  THE LIST 'NEXT' WILL BE REALIZED
C         IMPLICITLY IN 'LVLNUM' AS:
C                  LVLNUM(I) > 0   <== LEVEL NUMBER OF NODE.  'NEXT' IS
C                                      SET WITH LVLNUM(I) = LEVEL+1
C                  LVLNUM(I) = 0   <== I-TH NODE IS IN 'QUEUED' OR IS
C                                      NOT IN THIS COMPONENT OF GRAPH,
C                                      OR HAS JUST BEEN NUMBERED.
C                  LVLNUM(I) < 0   <== I-TH NODE IS IN 'ACTIVE' AND IS
C                                      CONNECTED TO -LVLNUM(I) NODES IN
C                                      'NEXT'.
C
C     ==================================================================
C
C     STRUCTURE OF WORKSPACE ..
C
C     --------------------------------------------------------------
C     : NUMBERED : DONE : ACTIVE : ALEVEL : ... : QUEUED : LVLLST :
C     --------------------------------------------------------------
C
C     -------------------
C       LVLPTR : LVLNUM :
C     -------------------
C
C     IN THE ABOVE,
C         NUMBERED IS THE SET OF NODES ALREADY NUMBERED FROM
C         PREVIOUS COMPONENTS AND EARLIER LEVELS OF THIS COMPONENT.
C         DONE, ACTIVE, ALEVEL  ARE VECTORS OF LENGTH THE WIDTH OF
C         THE CURRENT LEVEL.  ACTIVE IS A SET OF INDICES INTO
C         ALEVEL.  AS THE NODES IN ALEVEL ARE NUMBERED, THEY
C         ARE PLACED INTO 'DONE'.
C         QUEUED IS A QUEUE OF NODES IN THE 'NEXT' LEVEL, WHICH
C         GROWS FROM THE START OF THE 'NEXT' LEVEL IN LVLLST
C         FORWARDS TOWARD 'ALEVEL'.  QUEUED IS OF LENGTH NO MORE
C         THAN THE WIDTH OF THE NEXT LEVEL.
C         LVLLST IS THE LIST OF UNNUMBERED NODES IN THE TREE,
C         ARRANGED BY LEVEL.
C
C     ==================================================================
      INTEGER     I, J, K, PTR, JPTR, KPTR, LPTR, MPTR, PPTR, RPTR,
     1            MPPTR, JNODE, KNODE, CNODE, LEVEL, LOWDG, UNUSED,
     2            MXQUE, NNEXT, ASTART, MINDG, LSTART, LWIDTH, ACTIVE,
     2            QUEUEB, QUEUED, QCOUNT, NCONNC, NACTIV, CDGREE,
     3            LDGREE, NFINAL, JDGREE, STRTIC, ADDED, TWRKLN,
     4            LVLLSL, CONNEJ, CONNER, ASTPTR, ACTPTR, ACTIVI,
     5            ASTRTI, QUEUEI, ACPPTR
C
C     ------------------------------------------------------------------
C
      TWRKLN = WRKLEN + NCOMPN + N + DEPTH + 1
      UNUSED = TWRKLN
C
      ASTART = LVLPTR(1)
      LWIDTH = LVLPTR(2) - ASTART
      ASTART = WRKLEN  + 1
      ACTIVE = NXTNUM + LWIDTH + 1
      NACTIV = LWIDTH
      NFINAL = NXTNUM + NCOMPN
C
      NNEXT = LVLPTR(3) - LVLPTR(2)
      QUEUED = WRKLEN
      QUEUEB = QUEUED
      MXQUE = ACTIVE + LWIDTH
C
C     ... BUILD FIRST PRIORITY QUEUE 'ACTIVE'
C
      LOWDG = - (N + 1)
      LPTR = LVLPTR(1)
      DO 200 I = 1, LWIDTH
          NCONNC = 0
          LVLLSL= LVLLST (LPTR)
          JPTR = RSTART (LVLLSL)
          LDGREE = DEGREE(LVLLSL)
          DO 100 J = 1, LDGREE
              CONNEJ = CONNEC (JPTR)
              IF ( LVLNUM (CONNEJ) .EQ. 2 )  NCONNC = NCONNC - 1
              JPTR = JPTR + 1
  100     CONTINUE
C
          ACTIVI = ACTIVE + I - 1
          WORK (ACTIVI) = I
          LVLNUM (LVLLSL) = NCONNC
          LOWDG = MAX0 (LOWDG, NCONNC)
          LPTR = LPTR + 1
  200 CONTINUE
      WORK (ACTIVE-1) = 0
C
C     -----------------------------------
C     NOW NUMBER NODES LEVEL BY LEVEL ...
C     -----------------------------------
C
      DO 2000 LEVEL = 1, DEPTH
C
C         ... NUMBER ALL NODES IN THIS LEVEL
C
          DO 1100 I = 1, LWIDTH
              PPTR = -1
              PTR = WORK (ACTIVE-1)
              IF (NNEXT .EQ. 0)  GO TO 1000
C
C                 ... IF NODES REMAIN IN NEXT, FIND THE EARLIEST NODE
C                     IN ACTIVE OF MINIMAL DEGREE.
C
                  MINDG = -(N+1)
                  DO 400 J = 1, NACTIV
                      ASTPTR = ASTART + PTR
                      CNODE = WORK (ASTPTR)
                      IF ( LVLNUM (CNODE) .EQ. LOWDG )  GO TO 500
                      IF ( LVLNUM (CNODE) .LE. MINDG )  GO TO 300
                          MPPTR = PPTR
                          MPTR = PTR
                          MINDG = LVLNUM (CNODE)
  300                 PPTR = PTR
                      ACTPTR = ACTIVE + PTR
                      PTR = WORK (ACTPTR)
  400             CONTINUE
C
C                     ... ESTABLISH  PTR  AS FIRST MIN DEGREE NODE
C                         PPTR AS PREDECESSOR IN LIST.
C
                  PTR = MPTR
                  PPTR = MPPTR
C
  500             ASTPTR = ASTART + PTR
                  CNODE = WORK (ASTPTR)
                  LOWDG = LVLNUM (CNODE)
                  LVLNUM (CNODE) = 0
                  JPTR = RSTART (CNODE)
C
C                 ... UPDATE CONNECTION COUNTS FOR ALL NODES WHICH
C                     CONNECT TO  CNODE'S  NEIGHBORS IN  NEXT.
C
                  CDGREE = DEGREE(CNODE)
                  STRTIC = QUEUEB
C
                  DO 700 J = 1, CDGREE
                      JNODE = CONNEC (JPTR)
                      JPTR = JPTR + 1
                      IF (LVLNUM (JNODE) .NE. LEVEL+1 )  GO TO 700
                          IF (QUEUEB .LT. MXQUE)  GO TO 5000
                          WORK (QUEUEB) = JNODE
                          QUEUEB = QUEUEB - 1
                          NNEXT = NNEXT - 1
                          LVLNUM (JNODE) = 0
                          IF  (NACTIV .EQ. 1)  GO TO 700
                            KPTR = RSTART (JNODE)
                            JDGREE = DEGREE (JNODE)
                            DO 600 K = 1, JDGREE
                                KNODE = CONNEC (KPTR)
                                KPTR = KPTR + 1
                                IF (LVLNUM (KNODE) .GE. 0)  GO TO 600
                                    LVLNUM (KNODE) = LVLNUM (KNODE) + 1
                                    IF  (LOWDG .LT. LVLNUM(KNODE))
     1                                   LOWDG = LVLNUM(KNODE)
  600                       CONTINUE
  700             CONTINUE
C
C                 ... TO MIMIC THE ALGORITHM AS IMPLEMENTED BY GIBBS,
C                     SORT THE NODES JUST ADDED TO THE QUEUE INTO
C                     INCREASING ORDER OF ORIGINAL INDEX. (BUT, BECAUSE
C                     THE QUEUE IS STORED BACKWARDS IN MEMORY, THE SORT
C                     ROUTINE IS CALLED FOR DECREASING INDEX.)
C
C                     TREAT  0, 1 OR 2  NODES ADDED AS SPECIAL CASES
C
                  ADDED = STRTIC - QUEUEB
                  IF  (ADDED - 2)  1000, 800, 900
C
  800                 IF (WORK(STRTIC-1) .GT. WORK(STRTIC))  GO TO 1000
                          JNODE = WORK(STRTIC)
                          WORK(STRTIC) = WORK(STRTIC-1)
                          WORK(STRTIC-1) = JNODE
                          GO TO 1000
C
  900                 CALL GPSKCO (ADDED, WORK(QUEUEB+1), ERROR)
                      IF  (ERROR .NE. 0)  GO TO 5500
C
C
C                 ... NUMBER THIS NODE AND DELETE IT FROM 'ACTIVE'.
C                     MARK IT UNAVAILABLE BY CHANGING SIGN OF DEGREE
C
 1000         NACTIV = NACTIV - 1
              ASTPTR = ASTART + PTR
              CNODE = WORK (ASTPTR)
              WORK (NXTNUM) = CNODE
              DEGREE (CNODE) = -DEGREE (CNODE)
              NXTNUM = NXTNUM + 1
C
C             ... DELETE LINK TO THIS NODE FROM LIST
C
              ACPPTR = ACTIVE + PPTR
              ACTPTR = ACTIVE + PTR
              WORK (ACPPTR) = WORK (ACTPTR)
 1100     CONTINUE
C
C         ... NOW MOVE THE QUEUE 'QUEUED' FORWARD, AT THE SAME
C             TIME COMPUTING CONNECTION COUNTS FOR ITS ELEMENTS.
C             THEN DO THE SAME FOR THE REMAINING NODES IN 'NEXT'.
C
          UNUSED = MIN0 (UNUSED, QUEUEB - MXQUE)
          IF ( NXTNUM .NE. ACTIVE-1 )  GO TO 5100
          IF ( LEVEL .EQ. DEPTH ) GO TO 2000
              LSTART = LVLPTR (LEVEL+1)
              LWIDTH = LVLPTR (LEVEL+2) - LSTART
              ACTIVE = NXTNUM + LWIDTH + 1
              ASTART = ACTIVE + LWIDTH
              NACTIV = LWIDTH
              MXQUE = ASTART + LWIDTH
              IF ( MXQUE .GT. QUEUEB + 1 )  GO TO 5000
              UNUSED = MIN0 (UNUSED, QUEUEB - MXQUE + 1)
C
              QCOUNT = QUEUED - QUEUEB
              LOWDG = -N-1
              WORK (ACTIVE-1) = 0
C
              PTR = LSTART
              DO 1600 I = 1, LWIDTH
C
C                 ... CHOOSE NEXT NODE FROM EITHER 'QUEUED' OR 'NEXT'
C
                  IF (I .GT. QCOUNT )  GO TO 1200
                      QUEUEI = QUEUED + 1 - I
                      CNODE = WORK (QUEUEI)
                      GO TO 1300
C
 1200                 CNODE = LVLLST (PTR)
                      PTR = PTR + 1
                      IF ( PTR .GT. LVLPTR(LEVEL+2) )  GO TO 5200
                          IF (LVLNUM (CNODE) .GT. 0)  GO TO 1300
                              GO TO 1200
C
 1300             IF ( LEVEL+1 .EQ. DEPTH ) GO TO 1500
C
                      RPTR = RSTART (CNODE)
                      NCONNC = 0
                      JDGREE = DEGREE (CNODE)
                      DO 1400 J = 1, JDGREE
                          CONNER = CONNEC (RPTR)
                          IF ( LVLNUM (CONNER) .EQ. LEVEL+2 )
     1                        NCONNC = NCONNC - 1
                          RPTR = RPTR + 1
 1400                 CONTINUE
                      LVLNUM (CNODE) = NCONNC
                      LOWDG = MAX0 (LOWDG, NCONNC)
C
C             ... ADD CNODE TO NEW 'ACTIVE' QUEUE
C
 1500             ACTIVI = ACTIVE + (I - 1)
                  ASTRTI = ASTART + (I - 1)
                  WORK (ACTIVI) = I
                  WORK (ASTRTI) = CNODE
 1600         CONTINUE
C
              IF (DEPTH .EQ. LEVEL+1 ) GO TO 1700
                  NNEXT = LVLPTR (LEVEL+3) - LVLPTR (LEVEL+2)
                  QUEUED = LSTART - 1 + LWIDTH + WRKLEN
                  QUEUEB = QUEUED
                  GO TO 2000
C
 1700             NNEXT = 0
C
 2000 CONTINUE
C
      IF  (NXTNUM .NE. NFINAL)  GO TO 5300
      SPACE = MAX0 (SPACE, TWRKLN - UNUSED)
      RETURN
C
C
C     ------------------------------------------------------------------
C
 5000 SPACE = NACTIV + NNEXT
      ERROR = 160
      RETURN
C
 5100 ERROR = 61
      GO TO 5400
C
 5200 ERROR = 62
      GO TO 5400
C
 5300 ERROR = 63
C
 5400 RETURN
C
 5500 ERROR = 64
      GO TO 5400
C
      END
