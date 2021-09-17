      SUBROUTINE GETNRM(VCORG,WORK,A,NDIM,NOFVAR,NPOIN,INNODE,NRMINF,
     +                  NRM2)
C
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NPOIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NOFVAR,NPOIN),NRM2(NOFVAR),NRMINF(NOFVAR),
     +                 VCORG(NDIM,NPOIN),WORK(3,*)
      INTEGER INNODE(NOFVAR)
C     ..
C     .. Local Scalars ..
      INTEGER IDIM,IPOIN,IVAR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER IDAMAX
      EXTERNAL DNRM2,IDAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      DO 5 IVAR = 1,NOFVAR
          IPOIN = IDAMAX(NPOIN,A(IVAR,1),NOFVAR)
          INNODE(IVAR) = IPOIN
          NRMINF(IVAR) = ABS(A(IVAR,IPOIN))
          DO 7 IDIM = 1,NDIM
              WORK(IDIM,IVAR) = VCORG(IDIM,IPOIN)
    7     CONTINUE
          NRM2(IVAR) = DNRM2(NPOIN,A(IVAR,1),NOFVAR)
    5 CONTINUE
C
      RETURN

      END
