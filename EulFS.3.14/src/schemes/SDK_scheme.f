!> \copydetails LDA_SCHEME()
      SUBROUTINE SDK_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     &                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     &                      MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C This routine computes Sidilkover's Q scheme on one triangle/tetrahedron
C
C     $Id: SDK_scheme.f,v 1.2 2013/01/24 07:46:33 abonfi Exp $ 
C
C
C
C
C
C
C
C
C     .. Parameters ..
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CELRES,SOURCE
      INTEGER IELEM,NDIM,NOFVERT
      LOGICAL MATRIX_ASSEMBLY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ADVECTION(NDIM),BETA(NOFVERT),
     &                 STIFC(NOFVERT,NOFVERT),
     &                 DT(NOFVERT),NODRES(NOFVERT),Q(NOFVERT),
     &                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DENOM,KNEGINV,KPOS,S,SUM,UIN,ABSPHI,HELP
      INTEGER I,I1,IROW,IVERT,J,JCOL,POSI,NEGI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT),PHI(MAXNOFVERT),
     &PHINL(MAXNOFVERT),P(MAXNOFVERT)
      INTEGER POS(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,MINMOD
      EXTERNAL DDOT,MINMOD
C     ..
      CELRES = ZERO
      ABSPHI = ZERO
      KPOS = ZERO
      UIN = ZERO
      POSI = 0
      NEGI = 0
C
C Loops over all vertices
C
      DO 10 IVERT = 1,NOFVERT
c
c Dotting advection speed with normal
c
          SUM = DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)
          K(IVERT) = SUM/NDIM
c
          HELP = Q(IVERT)*K(IVERT)
c
          CELRES = CELRES + HELP
          ABSPHI = ABSPHI + ABS(HELP)
c
          IF (K(IVERT).GT.ZERO) THEN
              KPOS = KPOS + K(IVERT)
              POSI = POSI+1
              POS(POSI) = IVERT
          ELSE
              NEGI = NEGI-1
              UIN = UIN - HELP
          ENDIF
c
   10 CONTINUE
C
      UIN = UIN/KPOS
C
      DO 20 I = 1,NOFVERT
         PHI(I) = MAX(K(I),ZERO)*(Q(I)-UIN) 
c
c Linear time-step for all the downstream nodes
c
         DT(I) = DT(I) + MAX(ZERO,K(I))
c Sidilkover's ratio
         P(I) = ONE-(TWO/(ONE+SIGN(ONE,PHI(I))*ABSPHI/CELRES)) 
   20 CONTINUE
C
C Target tracking
C
      DO 30 I = 1,NOFVERT
          PHINL(I) = PHI(I)
          DO 35 J = 1,NOFVERT
             IF(J.EQ.I)GOTO 30
             PHINL(I) = PHINL(I) + MINMOD(P(J),ONE)*PHI(J)    
   35     CONTINUE
          NODRES(I) = -PHINL(I)
          BETA(I) = PHINL(I)/CELRES
   30 CONTINUE
C
C
C  the element stiffness matrix is that of the N scheme
C
      IF (MATRIX_ASSEMBLY) THEN
          KNEGINV = -ONE/KPOS
C
C     The convection matrix has to be zeroth since in the
C     subsequent loops (28,30) not all vertices are touched
C
          DO 40 J = 1,NOFVERT
              DO 40 I = 1,NOFVERT
                  STIFC(I,J) = ZERO
   40     CONTINUE
C
          DO 32 I = 1,POSI
              IROW = POS(I)
              S = K(IROW)
              DO 28 J = NOFVERT,NEGI,-1
                  JCOL = POS(J)
                  STIFC(IROW,JCOL) = S*K(JCOL)*KNEGINV
   28         CONTINUE
              STIFC(IROW,IROW) = -S
   32     CONTINUE
      ENDIF
C

      RETURN

      END
C
