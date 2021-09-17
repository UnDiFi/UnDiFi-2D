!> @copydetails LDA_SCHEME()
      SUBROUTINE NS_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     +                     NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                     MATRIX_ASSEMBLY)
C
C     $Id: NS_scheme.f,v 1.5 2013/01/24 07:46:33 abonfi Exp $
C
      IMPLICIT NONE
C
C This routine computes the N scheme on one triangle/tetrahedron
C
C the treatment of the source term with N scheme
C follows Mario Ricchiuto VKI PR ? 2001
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
      DOUBLE PRECISION FLUCT,SIGIN,ANVT
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
      SIGIN = ZERO
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
c reduced value of \Sigma_{in} assuming constant source term
              SIGIN = SIGIN - SOURCE/NOFVERT
          ELSE
              NEGI = NEGI - 1
              POS(NEGI) = IVERT
              UIN = UIN - K(IVERT)*Q(IVERT)
              NODRES(IVERT) = ZERO
          ENDIF
c
   10 CONTINUE
C
      IF (KPOS.EQ.ZERO) RETURN
      UIN = UIN/KPOS
      SIGIN = (SIGIN - SOURCE)/POSI
C
C Target tracking and updating
C
      FLUCT = ZERO
      DO 30 I = 1,POSI
          J = POS(I)
C N-scheme contribution for the convective part and source term
caldo     DQ(J) = K(J)* (Q(J)-UIN) - SOURCE/NOFVERT - SIGIN
C
C when the source term is constant it is simply
C equally split among the downstream nodes
C
          DQ(J) = K(J)* (Q(J)-UIN) + SOURCE/POSI
c
          FLUCT = FLUCT + DQ(J)
c
c Linear time-step for all the downstream nodes
c
          DT(J) = DT(J) + K(J)
          NODRES(J) = -DQ(J)
c
   30 CONTINUE

      IF (DABS(FLUCT-(CELRES+SOURCE)).GT.1.D-12) THEN
         write(6,*) "Check NwS.........",FLUCT,CELRES+SOURCE
      END IF

C
C  the element stiffness matrix is that of the N scheme
C  does not account for the source term
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
