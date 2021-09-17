!> \copydetails LDA_SCHEME()
      SUBROUTINE NL_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     &                     NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                     MATRIX_ASSEMBLY)
C
      IMPLICIT NONE
C
C     $Id: NL_scheme.f,v 1.5 2013/01/24 07:46:33 abonfi Exp $
C
C
C This routine computes the NL scheme on one triangle/tetrahedron
C
C the limiting procedure is based on the L(x,y) limiter functions:
C
C \phi_i^{NL} = \phi_i^{N} + \sum_{j \neq i} L( -\phi_i, \phi_j )
C
C the limiter is only applied to the convective part
C while the source term (if any) is distributed using the LDA scheme
C
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
     2                 STIFC(NOFVERT,NOFVERT),
     +                 DT(NOFVERT),NODRES(NOFVERT),Q(NOFVERT),
     +                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DENOM,KNEGINV,KPOS,S,SUM,UIN
      INTEGER I,I1,IROW,IVERT,J,J1,J2,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DQ(MAXNOFVERT),K(MAXNOFVERT),PHI(MAXNOFVERT)
      INTEGER POS(MAXNOFVERT)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,MINMOD
      INTEGER ICYCL
      EXTERNAL DDOT,MINMOD,ICYCL
C     ..
      POSI = 0
      NEGI = NOFVERT + 1
      CELRES = ZERO
      KPOS = ZERO
      UIN = ZERO
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
          CELRES = CELRES + Q(IVERT)*K(IVERT)
c
          IF (K(IVERT).GT.ZERO) THEN
              POSI = POSI + 1
              POS(POSI) = IVERT
              KPOS = KPOS + K(IVERT)

          ELSE
              NEGI = NEGI - 1
              POS(NEGI) = IVERT
              UIN = UIN - K(IVERT)*Q(IVERT)
              NODRES(IVERT) = ZERO
          ENDIF
c
          BETA(IVERT) = ZERO
c
   10 CONTINUE
C
      IF (KPOS.EQ.ZERO) RETURN
      UIN = UIN/KPOS
C
C Target tracking
C
      S = ZERO
      DO 30 I = 1,POSI
          J = POS(I)
C N-scheme contribution for the convective part
          DQ(J) = K(J)* (Q(J)-UIN)
C LDA-scheme contribution for the source term
          PHI(J) = DQ(J) + K(J)*SOURCE/KPOS
          S = S + PHI(J)
c
c Linear time-step for all the downstream nodes
c
          DT(J) = DT(J) + K(J)
c
   30 CONTINUE
C
C Loop over downstream nodes (TARGET UPDATING)
C
      DO 20 I = 1,POSI
          J1 = POS(I)
          DO 22 I1 = 1,POSI - 1
              J2 = POS(ICYCL(I+I1,POSI))
              PHI(J1) = PHI(J1) - MINMOD(DQ(J1),-DQ(J2))
   22     CONTINUE
c
          NODRES(J1) = -PHI(J1)
          BETA(J1) = PHI(J1)/S
c
   20 CONTINUE
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
C
C     IF( MATRIX_ASSEMBLY )THEN
C
C        DENOM = CELRES+SOURCE
C
C        IF(DENOM.LT.1.D-15)RETURN
C
C        DO 40 J = 1 , NOFVERT
C           DO 40 I = 1 , NOFVERT
C              STIFC(I,J) = ZERO
C  40    CONTINUE
C        DO 34 I = 1 , POSI
C           IROW = POS(I)
C           S =-PHI(IROW)*DENOM
C           DO 28 JCOL = 1, NOFVERT
C              STIFC(IROW,JCOL) = S * K(JCOL)
C  28       CONTINUE
C  34    CONTINUE
C     ENDIF

      RETURN

      END
C
