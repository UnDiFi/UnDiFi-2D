      SUBROUTINE   GPSKCL   (N, DEGREE, RSTART, CONNEC, INVNUM, NEWNUM, GPSK2154
     1                       OLDNUM, BANDWD, PROFIL, ERROR, SPACE)
C
C
      INTEGER     N, RSTART(N), BANDWD, PROFIL, ERROR, SPACE
C
CIBM  INTEGER *2  DEGREE(N), CONNEC(1), INVNUM(N), NEWNUM(N), OLDNUM(N)
      INTEGER     DEGREE(N), CONNEC(1), INVNUM(N), NEWNUM(N), OLDNUM(N)
C
C     ==================================================================
C
C
C     COMPUTE THE BANDWIDTH AND PROFILE FOR THE RENUMBERING GIVEN
C     BY 'INVNUM' AND ALSO FOR THE RENUMBERING GIVEN BY 'OLDNUM'.
C     'NEWNUM' WILL BE A PERMUTATION VECTOR COPY OF THE NODE
C     LIST 'INVNUM'.
C
C     ==================================================================
C
      INTEGER     I, J, JPTR, IDGREE, OLDBND, OLDPRO, NEWBND, NEWPRO,
     1            OLDRWD, NEWRWD, OLDORG, NEWORG, JNODE, INVNMI
C
C     ------------------------------------------------------------------
C
C     ... CREATE NEWNUM AS A PERMUTATION VECTOR
C
      DO 100 I = 1, N
          INVNMI = INVNUM (I)
          NEWNUM (INVNMI) = I
  100 CONTINUE
C
C     ... COMPUTE PROFILE AND BANDWIDTH FOR BOTH THE OLD AND THE NEW
C         ORDERINGS.
C
      OLDBND = 0
      OLDPRO = 0
      NEWBND = 0
      NEWPRO = 0
C
      DO 300 I = 1, N
          IF (DEGREE(I) .EQ. 0)  GO TO 300
          IF (DEGREE(I) .GT. 0)  GO TO 6000
              IDGREE = -DEGREE(I)
              DEGREE(I) = IDGREE
              NEWORG = NEWNUM(I)
              OLDORG = OLDNUM(I)
              NEWRWD = 0
              OLDRWD = 0
              JPTR = RSTART (I)
C
C             ... FIND NEIGHBOR WHICH IS NUMBERED FARTHEST AHEAD OF THE
C                 CURRENT NODE.
C
              DO 200 J = 1, IDGREE
                  JNODE = CONNEC(JPTR)
                  JPTR = JPTR + 1
                  NEWRWD = MAX0 (NEWRWD, NEWORG - NEWNUM(JNODE))
                  OLDRWD = MAX0 (OLDRWD, OLDORG - OLDNUM(JNODE))
  200         CONTINUE
C
              NEWPRO = NEWPRO + NEWRWD
              NEWBND = MAX0 (NEWBND, NEWRWD)
              OLDPRO = OLDPRO + OLDRWD
              OLDBND = MAX0 (OLDBND, OLDRWD)
  300 CONTINUE
C
C     ... IF NEW ORDERING HAS BETTER BANDWIDTH THAN OLD ORDERING,
C         REPLACE OLD ORDERING BY NEW ORDERING
C
      IF  (NEWBND .GT. OLDBND)  GO TO 500
          BANDWD = NEWBND
          PROFIL = NEWPRO
          DO 400 I = 1, N
              OLDNUM(I) = NEWNUM(I)
  400     CONTINUE
          GO TO 600
C
C     ... RETAIN OLD ORDERING
C
  500     BANDWD = OLDBND
          PROFIL = OLDPRO
C
  600 RETURN
C
C     ------------------------------------------------------------------
C
 6000 SPACE = -1
      ERROR = 70
      RETURN
C
      END
