      SUBROUTINE FDUMP                                                  FDMP0000
C  THIS IS A DUMMY ROUTINE TO BE SENT OUT ON
C  THE PORT SEDIT TAPE
C
      INCLUDE 'backup.com'
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
      COMMON /CSTAK/DSTAK
      EQUIVALENCE(ISTAK(1),DSTAK(1))
C
      LWORK = ISTKGT(MAX(NOFVAR,NTURB)*NPNOD,2)
      CALL BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB,ISTAK(LPMAP),
     +            ISTAK(LWORK),BAKFILE,VISCTFILE)
      CALL ISTKRL(1)
C
      RETURN
      END
