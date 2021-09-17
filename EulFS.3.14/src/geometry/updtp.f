      SUBROUTINE UPDTP( V1, V2, NDOF, NPOIN, NGHOST, MAP, INDX,
     & NPNOD, LROT )
C
C     $Id: updtp.f,v 1.1 2005/09/09 08:40:54 abonfi Exp $
C
C Subroutine for updating periodic nodes
C In the periodic case the solution in those periodic nodes
C that have been removed, i.e. those stored in
C NPOIN+NGHOST ......... NPOIN+NGHOST+NPNOD is reconstructed
C using the mapping MAP
C
      INCLUDE 'paramt.h'
C     INCLUDE 'periodic.com'
C
C     .. Scalar Arguments ..
      INTEGER NDOF,NPOIN,NGHOST,NPNOD
      LOGICAL LROT
C     .. Array Arguments ..
      DOUBLE PRECISION V1(*),V2(*)
C
C     INDX must be dimensioned MAX(NDOF,NTURB)*NPNOD
C
      INTEGER MAP(NPNOD),INDX(*)
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      INTEGER IPOIN,IDOF,LOC,I
C     ..
C     .. Local Arrays ..
C     ..
C     .. External Subroutines ..
C     ..
C     .. Common blocks ..
C     ..
C     .. Equivalences ..
C     ..
C
C     We check NPNOD since the mesh might be
C     periodic, BUT there might be no
C     periodic nodes on this processor
C
      IF(NPNOD.EQ.0)RETURN
      LOC = 0
      DO 1 I = 1, NPNOD
         IPOIN = MAP(I)
         DO 1 IDOF = 1,NDOF
            LOC = LOC + 1
            INDX(LOC) = (IPOIN-1)*NDOF+IDOF
    1 CONTINUE
      CALL DGTHR( NDOF*NPNOD, V1, V2, INDX )
C
C    rotate velocities, if required
C
      IF( LROT )CALL ROTATE( V2, NDOF, NPNOD )
C
      RETURN
      END
