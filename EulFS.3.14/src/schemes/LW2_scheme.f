      SUBROUTINE LW2_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                      MATRIX_ASSEMBLY)
C
C     $Id: LW2_scheme.f,v 1.3 2012/04/17 10:39:29 abonfi Exp $
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.com'
C
      DOUBLE PRECISION CFLELEM
      PARAMETER (CFLELEM=HALF)
Caldo
Caldo This scheme should be time accurate
Caldo see De Palma et al JCP 208 (2005) 1-33
Caldo
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),DT(NOFVERT),NODRES(NOFVERT),
     +                 Q(NOFVERT),BETA(NOFVERT),STIFC(NOFVERT,NOFVERT),
     +                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S,HELP
      INTEGER I,IVERT,J,IADDR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
      CELRES = ZERO
C
      DO 10 IVERT = 1,NOFVERT
C
C Dotting advection speed with the face normal
C
          HELP = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = HELP/NDIM
C
          CELRES = CELRES + Q(IVERT)*K(IVERT)
C
          KSUM = KSUM + ABS(K(IVERT))
C
   10 CONTINUE
C
      CELRES = CELRES!+ SOURCE
C
C Loops over all nodes
C
      HELP = ONE/REAL(NOFVERT)
      KSUM = CFLELEM*DTVOL
      DO 20 I = 1,NOFVERT
          IF(K(I).GT.ZERO)DT(I) = DT(I) + K(I)
          BETA(I) = HELP + KSUM*K(I)
          NODRES(I) = -BETA(I)*CELRES-HELP*SOURCE
   20 CONTINUE
C
C
C
      IF (MATRIX_ASSEMBLY) THEN
          DO 40 J = 1,NOFVERT
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = -BETA(I)*K(J)
   40     CONTINUE
      ENDIF
C
C
      RETURN

      END
