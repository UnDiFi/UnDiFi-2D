C
C
      SUBROUTINE TDIFF(IELEM,EPS,W,NODRES,DT,VCN,DMAT,
     &NDIM,NOFVERT,VOLUME)
C
      IMPLICIT NONE
C
C     NON conservative Diffusion term: TCB2 / TPR1 * (\nabla u)**2
C     from the Spalart-Allmaras model
C
      INCLUDE 'paramt.h'
      INCLUDE 'three.com'
C
      INTEGER IELEM,NDIM,NOFVERT
      DOUBLE PRECISION DT(NOFVERT),NODRES(NOFVERT),W(NOFVERT),
     +DMAT(NOFVERT,NOFVERT),VCN(NDIM,NOFVERT),VOLUME 
      DOUBLE PRECISION EPS
      DOUBLE PRECISION TEMP
      INTEGER I
      DOUBLE PRECISION GRADNRM
      DOUBLE PRECISION d1(4),d2
      common /diffune/d1,d2
C
C     Executable Statements ..
C
C     potrebbe essere piu` appropriato ricalcolare il
C     gradiente ed aggiungere la parte implicita
C
      GRADNRM = GRAD_PARM(1,1)**2 + GRAD_PARM(1,2)**2
     +        + GRAD_PARM(1,3)**2
C
      TEMP = EPS * GRADNRM * VOLUME / NOFVERT
      d2 = temp
C
      DO 2 i = 1 , NOFVERT
CXXXX    NODRES(I) = NODRES(I) - TEMP
         NODRES(I) = NODRES(I) + TEMP
    2 CONTINUE
C
      RETURN
      END
C
