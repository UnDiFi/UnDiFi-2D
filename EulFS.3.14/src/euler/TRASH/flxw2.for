C
      SUBROUTINE FLXW2(NDIM,NOFVAR,NOFVERT,STIFC,WORK,WORK2,VCZ,VCN,
     +                 NODRES,PICARD)



C

C
C     .. Parameters ..
      DOUBLE PRECISION ALPHA
      PARAMETER (ALPHA=0.75d0)
      DOUBLE PRECISION ZERO,HALF,ONE,TWO
      PARAMETER (ZERO=0.00d0,HALF=0.5d0,ONE=1.00d0,TWO=2.00d0)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT
      LOGICAL PICARD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION NODRES(NOFVAR,NDIM),
     +                 STIFC(NOFVAR,NOFVAR,NOFVERT,NOFVERT),VCN(NDIM),
     +                 VCZ(NOFVAR,NDIM),WORK(NOFVAR,NOFVAR,NOFVERT-1),
     +                 WORK2(NOFVAR,NOFVAR,NOFVERT-1,NOFVERT-1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BETA,CNST
      INTEGER I,IADD,IFAIL,IVERT,J,JVERT,K,L
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLUX(16)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DINIT,GETDF2CORRDU,INVWLLI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
      DO 2 IVERT = 1,NOFVERT
          IADD = (IVERT-1)*NOFVAR + 1
          CALL INVWLLI(NDIM,VCN,VCZ(1,IVERT),FLUX(IADD))
    2 CONTINUE
C
      CALL DINIT(NOFVAR* (NOFVERT-1),ZERO,NODRES,1)
      BETA = (1.-ALPHA)/REAL(NDIM-1)/REAL(NDIM)
      DO 5 IVERT = 1,NOFVERT - 1
          DO 7 JVERT = 1,NOFVERT - 1
              IADD = (JVERT-1)*NOFVAR + 1
              IF (JVERT.EQ.IVERT) THEN
                  CNST = ALPHA/NDIM

              ELSE
                  CNST = BETA
              ENDIF

              CALL DAXPY(NOFVAR,CNST,FLUX(IADD),1,NODRES(1,IVERT),1)
    7     CONTINUE
    5 CONTINUE
C
      IF (.NOT.PICARD) RETURN

      DO 3 IVERT = 1,NOFVERT - 1
          CALL GETDF2CORRDU(VCZ(2,IVERT),VCN,NDIM,NOFVAR,
     +                      WORK(1,1,IVERT))
    3 CONTINUE

      DO 8 I = 1,NOFVERT - 1
          DO 8 J = 1,NOFVERT - 1
              IF (J.EQ.I) THEN
                  CNST = ALPHA/REAL(NDIM)

              ELSE
                  CNST = BETA
              ENDIF

              DO 8 L = 1,NOFVAR
                  DO 8 K = 1,NOFVAR
                      STIFC(K,L,I,J) = 0.5d0*CNST*WORK(K,L,J)
    8 CONTINUE

      RETURN

  564 FORMAT ((E12.6,1X))

      END
