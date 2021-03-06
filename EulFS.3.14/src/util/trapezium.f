      SUBROUTINE TRAPEZIUM(VNOR,Z,BCFLUX,FLUXN,NDIM,NOFVAR)
C
C    This routine applies the trapezium rule to compute the flux
C    through an edge ..
C
      IMPLICIT NONE 
C
C
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FLUXN(NOFVAR),VNOR(NDIM),Z(NOFVAR,*)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL BCFLUX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION H3
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F0(5),F1(5),F2(5),Z1(5)
C     ..
C
      CALL BCFLUX(NDIM,VNOR,Z(1,1),F0)
      CALL BCFLUX(NDIM,VNOR,Z(1,2),F2)
      IF(NDIM.EQ.3)CALL BCFLUX(NDIM,VNOR,Z(1,3),F1)
C
      H3 = 1.D0/REAL(NDIM)
      IF(NDIM.EQ.2)THEN
      DO 7 I = 1,NOFVAR
          FLUXN(I) = H3* (F0(I)+F2(I))
    7 CONTINUE
      ELSEIF(NDIM.EQ.3)THEN
      DO 9 I = 1,NOFVAR
          FLUXN(I) = H3* (F0(I)+F2(I)+F1(I))
    9 CONTINUE
      ENDIF
C
      RETURN

      END
