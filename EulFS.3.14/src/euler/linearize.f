      SUBROUTINE LINEARIZE(IELEM,ALE,VCN,VCB,NDIM,NOFVERT,
     +VCZ,NOFVAR,VOLUME)
C
      IMPLICIT NONE 
C
C THIS SUBROUTINE COMPUTES:
C a) The averaged state ZAVG(1:NOFVAR) over the cell (in cart. coord.)
C b) The gradient of the parameter vector GRAD_PARM(1:NOFVAR,1:NOFVERT) 
C    (in cart. coord.)
C c) The gradient of the primitive variables GRAD_PRIM(1:NOFVAR,1:NOFVERT) 
C    (in cart. coord.)
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NDIM,NOFVAR,NOFVERT
      LOGICAL ALE
      DOUBLE PRECISION VCN(NDIM,NOFVERT),VCB(NDIM,NOFVERT),
     2VCZ(NOFVAR,NOFVERT),VOLUME
C
C     .. Local Scalars ..
C
      INTEGER I,IVAR,JVERT
      DOUBLE PRECISION HELP,TEMP
C
C     .. Local Arrays ..
C
C
C     .. External Functions ..
C
C
C
C     .. Executable Statements ..
C
C *********************************************************************
C AVERAGED STATE OVER THE CELL (in cartesian coordinates)
C *********************************************************************
C
      DO 10 IVAR = 1 , NOFVAR
         HELP = ZERO
            DO 12 JVERT = 1 , NOFVERT
               HELP = HELP + VCZ( IVAR , JVERT )
   12       CONTINUE
      ZAVG(IVAR) = HELP / NOFVERT
   10 CONTINUE
C
C *********************************************************************
C AVERAGED grid velocity OVER THE CELL (in cartesian coordinates)
C *********************************************************************
C
      IF(ALE)THEN
         DO 30 IVAR = 1 , NDIM
            HELP = ZERO
               DO 32 JVERT = 1 , NOFVERT
                  HELP = HELP + VCB( IVAR , JVERT )
   32          CONTINUE
         BAVG(IVAR) = HELP / NOFVERT
   30 CONTINUE
      ENDIF
C
C *********************************************************************
C COMPUTES THE GRADIENT OF THE PARAMETER VECTOR (in cartesian coordinates)
C *********************************************************************
C
      DO 20 IVAR = 1 , NOFVAR
         DO 22 I = 1 , NDIM
          HELP = ZERO
            DO 24 JVERT = 1 , NOFVERT

               TEMP = VCN( I , JVERT )
               HELP = HELP + VCZ( IVAR , JVERT ) * TEMP 

   24       CONTINUE
         GRAD_PARM( IVAR , I ) = HELP / NDIM / VOLUME
   22    CONTINUE
   20 CONTINUE
C
      RETURN
      END
