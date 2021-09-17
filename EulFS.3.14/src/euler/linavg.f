      SUBROUTINE LINAVG(VCZ,ZAVG,NDIM,NOFVERT,NOFVAR)
C
      IMPLICIT NONE 
C
C THIS SUBROUTINE COMPUTES:
C a) The averaged state ZAVG(1:NOFVAR) over the cell (in cart. coord.)
C
      INCLUDE 'constants.h'
C
C     .. Scalar Arguments ..
C
      INTEGER NDIM,NOFVAR,NOFVERT
      DOUBLE PRECISION VCZ(NOFVAR,NOFVERT),ZAVG(NOFVAR)
C
C     .. Local Scalars ..
C
      INTEGER I,IVAR,JVERT
      DOUBLE PRECISION DUM,HELP
C
C     .. Local Arrays ..
C
C
C     .. External Functions ..
C
C
C     .. Executable Statements ..
C
C *********************************************************************
C AVERAGED STATE OVER THE CELL (in cartesian coordinates)
C *********************************************************************
C
      HELP = ONE/REAL(NOFVERT)
      DO 10 IVAR = 1 , NOFVAR
       DUM = ZERO
            DO 12 JVERT = 1 , NOFVERT
               DUM = DUM + VCZ( IVAR , JVERT )
   12       CONTINUE
      ZAVG(IVAR) = DUM * HELP
   10 CONTINUE
C
      RETURN
      END
