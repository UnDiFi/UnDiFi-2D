!> \copydetails LDA_SCHEME()
      SUBROUTINE NL2_SCHEME(IELEM,VCN,ADVECTION,CELRES,SOURCE,Q,DT,
     2                      NODRES,BETA,STIFC,NDIM,NOFVERT,
     3                      MATRIX_ASSEMBLY)
C
C     $Id: NL2_scheme.f,v 1.8 2012/04/26 16:27:09 abonfi Exp $
C
      IMPLICIT NONE
C
C This routine computes the NL scheme on one triangle/tetrahedron
C
C the limiting procedure is based on the formula
C
C \beta_i = \max(0,\phi_i^N\phi)/[\sum_j \max(0,\phi_j^N\phi)]
C
C
C it SHOULD be the same as the PSI scheme as long as 
C there is not a source term
C
C the source term is treated as in
C Sidilkover and Roe ICASE Report No. 95-10
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
     2                 STIFC(NOFVERT,NOFVERT),DT(NOFVERT),
     3                 NODRES(NOFVERT),Q(NOFVERT),
     +                 VCN(NDIM,NOFVERT)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DENOM,KNEGINV,KPOS,S,HELP,UIN,PHIT
      INTEGER I,I1,IROW,IVERT,J,JCOL,NEGI,POSI
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION K(MAXNOFVERT),PHI(MAXNOFVERT),TMP(MAXNOFVERT)
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
          K(IVERT) =  DDOT(NDIM,VCN(1,IVERT),1,ADVECTION,1)/NDIM
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
              BETA(IVERT) = ZERO
          ENDIF
c
   10 CONTINUE
C
      UIN = UIN/KPOS
      PHIT=-CELRES-SOURCE
      DENOM=ZERO
!     IF(ABS(PHIT).LE.1.E-14)RETURN
C
C     we distinguish de
C
      IF(POSI.EQ.1)THEN ! 1-target case
         J = POS(1)
         BETA(J) = ONE
         NODRES(J)=PHIT
      ELSEIF(POSI.EQ.0)THEN ! 0-target case: should not occur
         WRITE(6,*) 'IE = ',IELEM
         WRITE(6,*) 'K = ',(K(IVERT),IVERT=1,NOFVERT)
         WRITE(6,*) 'U = ',(ADVECTION(J),J=1,NDIM)
         WRITE(6,*) 'Cannot handle this in NL2'
         DO I = 1, NOFVERT
            NODRES(I) = ZERO
         ENDDO 
      ELSE ! multiple target case
C
C Target tracking
C
         DO 30 I = 1,POSI
             J = POS(I)
C N-scheme contribution for the convective part and
C LDA-scheme contribution for the source term
             PHI(J) = -K(J)* ( (Q(J)-UIN) + SOURCE/KPOS )
c
c Linear time-step for all the downstream nodes
c
             DT(J) = DT(J) + K(J)
c
             TMP(J) = MAX(ZERO,PHIT*PHI(J))
             DENOM=DENOM+TMP(J)
c
   30    CONTINUE
C
C Loop over downstream nodes (TARGET UPDATING)
C
!        IF(ABS(DENOM).LE.1.D-16)then
!        write(6,*) ielem,denom,posi,(phi(pos(i)),i=1,posi)
!        DENOM=SIGN(ONE,DENOM)*1.D-16
!        ENDIF 
C
C        we avoid testing DENOM, since this tipically hamper convergence
C
         PHIT=PHIT/DENOM
         HELP = ZERO
         DO 20 I = 1,POSI
             J = POS(I)
             BETA(J) = TMP(J)/DENOM
             HELP = HELP+BETA(J) ! this is a check
             NODRES(J)=PHIT*TMP(J)
   20    CONTINUE
      ENDIF ! check on POSI
!     WRITE(6,*)posi,HELP,(BETA(POS(I)),i=1,posi)
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
