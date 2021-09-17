      SUBROUTINE BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB, MAP,INDX,
     +BAKFILE,VISCTFILE,OLDFILE)
C
C Subroutine for backing up DATA
C in the periodic case the solution in those periodic nodes
C that have been removed, i.e. those stored in
C NPOIN+NGHOST ......... NPOIN+NGHOST+NPNOD is reconstructed
C using the mapping MAP
C
C     $Id: backup.f,v 1.13 2009/06/11 09:55:05 abonfi Exp $
C
C VISCTFILE is used to store viscosity in the un-coupled approach 
C BAKFILE   is used to stores flow variables at the current time level
C OLDFILE   is used to stores flow variables at the previous time level
C
      INCLUDE 'paramt.h'
      INCLUDE 'time.h'
      INCLUDE 'nloc.com'
      INCLUDE 'periodic.com'
      INCLUDE 'time.com'
C
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NPOIN,NTURB
      CHARACTER BAKFILE* (*),VISCTFILE* (*), OLDFILE* (*)
C     .. Array Arguments ..
      INTEGER MAP(NPNOD),INDX(*)
C
C     INDX must be dimensioned MAX(NOFVAR,NTURB)*NPNOD
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
      INTEGER IPOIN,IVAR,LOC,I,NITEMS
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. External Subroutines ..
      EXTERNAL SOLZNE
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
C
      IF(NPNOD.GT.0)THEN
          LOC = 0
          DO 1 I = 1, NPNOD
              IPOIN = MAP(I)
              DO 1 IVAR = 1,NOFVAR
                  LOC = LOC + 1
                  INDX(LOC) = (IPOIN-1)*NOFVAR+IVAR
    1     CONTINUE
          IADD = LZROE+NOFVAR*(NPOIN+NGHOST)
          CALL DGTHR( NOFVAR*NPNOD, DSTAK(LZROE), DSTAK(IADD), INDX )
C
C    rotate velocities
C
          IF( PERIODIC_MESH .AND. ANNULAR )
     &    CALL ROTATE( DSTAK(IADD), NOFVAR, NPNOD )
          IF(NTURB.NE.0)THEN
              LOC = 0
              DO 3 I = 1, NPNOD
                  IPOIN = MAP(I)
                  DO 3 IVAR = 1,NTURB
                      LOC = LOC + 1
                      INDX(LOC) = (IPOIN-1)*NTURB+IVAR
    3         CONTINUE
              IADD = LTURB+NTURB*(NPOIN+NGHOST)
              CALL DGTHR( NTURB*NPNOD, DSTAK(LTURB), DSTAK(IADD), INDX )
          ENDIF
      ENDIF
c
c     Saves the solution ...
c
      NITEMS = NPOIN+NGHOST+NPNOD
      CALL SOLZNE(BAKFILE,DSTAK(LZROE),NOFVAR,NITEMS,'write')
      IF(LTIME)THEN ! save time level n
         LOC = LZROE + NOFVAR*NITEMS
         CALL SOLZNE(OLDFILE,DSTAK(LOC),NOFVAR,NITEMS,'write')
      ENDIF 
C
C     there is NO unsteady uncoupled approach, so we do NOT save
C     the previous time level
C
      IF (NTURB.NE.0)
     +CALL SOLZNE(VISCTFILE,DSTAK(LTURB),NTURB,NITEMS,'write')
C
      RETURN
 
      END
