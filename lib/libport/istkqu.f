      INTEGER FUNCTION ISTKQU(ITYPE)                                    STKA0000
C
C  RETURNS THE NUMBER OF ITEMS OF TYPE ITYPE THAT REMAIN
C  TO BE ALLOCATED IN ONE REQUEST.
C
C  ERROR STATES -
C
C    1 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
C    2 - ITYPE .LE. 0 .OR. ITYPE .GE. 6
C
      COMMON /CSTAK/DSTAK
C
      DOUBLE PRECISION DSTAK(500)
      INTEGER ISTAK(1000)
      INTEGER ISIZE(5)
C
      LOGICAL INIT
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
      DATA INIT/.TRUE./
C
      IF (INIT) CALL I0TK00(INIT,500,4)
C
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
     1   (47HISTKQU - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN,
     2    47,1,2)
C
      IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
     1   (33HISTKQU - ITYPE.LE.0.OR.ITYPE.GE.6,33,2,2)
C
      ISTKQU = MAX0( ((LMAX-2)*ISIZE(2))/ISIZE(ITYPE)
     1             - (LNOW*ISIZE(2)-1)/ISIZE(ITYPE)
     2             - 1, 0 )
C
      RETURN
C
      END