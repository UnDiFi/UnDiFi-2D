      SUBROUTINE PERIODIC(NPOIN,NGHOST,NPNOD,NELEM,NOFVAR,
     +                    NOFVERT,ICELNOD,MAP)
C
      IMPLICIT NONE
C
C     $Id: periodic.F,v 1.12 2013/03/16 10:39:14 abonfi Exp $
C
C
      INCLUDE 'io.com' 
C
C     define GRAPHICS to see on screen the correspondance between 
C     periodic nodes
C     #define GRAPHICS
C
#ifdef GRAPHICS
      INCLUDE 'nloc' 
#endif
C
C
C     NPNOD is half the # of periodic nodes
C     Nodes must be ordered as follows in the datafile
C     the lists of periodic nodes that have been deleted must appear last
C     and there must be a correspondance between these
C     two lists; the correspondance is established by the array MAP
C
C
C     +-------------------------+-------+-------+
C     ^                        ^       ^       ^ 
C     |                        |       |       |
C     |                        |       |     NPOIN
C     1                        |       |
C                              |   NPOIN-NPNOD
C                              | 
C                        NPOIN-2*NPNOD
C
C     corresponding indexes of the two lists must be
C     the corresponding nodes
C
      INTEGER NPOIN,NGHOST,NPNOD,NELEM,NOFVAR,NOFVERT
      INTEGER ICELNOD(NOFVERT,NELEM),MAP(NPNOD)
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
      INTEGER IELEM,IVERT,IPOIN,IFAIL
#ifdef GRAPHICS
      integer INDGR,ipntr(2)
      external grper
#endif
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. External Functions ..
C     ..
C     .. External Subroutines ..
      EXTERNAL ISTKRL
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
C     ..
C
C     make a check on MAP
C     1 <= MAP(I) <= NPOIN+NGHOST
C
      IFAIL = 0
      do ivert = 1,npnod
         ipoin = MAP(IVERT)
         IF(.NOT.((1.LE.IPOIN).AND.(IPOIN.LE.(NPOIN+NGHOST))))THEN
             WRITE(6,*)'Smtgh. wrong with MAP(',ivert,') = ',
     &       ipoin 
             IFAIL = IFAIL + 1
         ENDIF
      enddo
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
C
C
C     loop over elements
C
      DO 200 IELEM = 1,NELEM
          DO 200 IVERT = 1, NOFVERT
             IPOIN = ICELNOD(IVERT,IELEM)
             IF( IPOIN.GT.(NPOIN+NGHOST) )THEN
                  ICELNOD(IVERT,IELEM) = MAP(IPOIN-(NPOIN+NGHOST))
             ENDIF
  200 CONTINUE
C
#ifdef GRAPHICS
      CALL GRINIT(5,6,'MESH')
      INDGR = 23
      CALL GRCTRL(GRPER,INDGR,'~x~y~TRPLOT',ISTAK(LCELNOD),
     +            NELEM,NPOIN,DSTAK(LCORG),1.d0,1.d0,1.d0,1.d0,
     +            istak(ipntr(1)),istak(ipntr(2)))
#endif
C

      RETURN
      END
