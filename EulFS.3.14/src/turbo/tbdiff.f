      DOUBLE PRECISION FUNCTION TBDIFF(GRADVT,NDIM)
C
C     $Id: tbdiff.f,v 1.2 2013/01/26 12:01:29 abonfi Exp $
C
      IMPLICIT NONE
C
C     NON conservative Diffusion term: TCB2 / TPR1 * (\nabla u)**2
C     from the Spalart-Allmaras model
C
      INCLUDE 'paramt.h'
      INCLUDE 'three.com'
      INCLUDE 'turb.com'
      INCLUDE 'visco.com'
C
C     GRADVT is just a dummy, should be used to store the
C            gradient of turb. viscosity
C
      INTEGER NDIM
      DOUBLE PRECISION GRADVT(NDIM)
      DOUBLE PRECISION GRADNRM
C
C     Executable Statements ..
C
C     potrebbe essere piu` appropriato ricalcolare il
C     gradiente ed aggiungere la parte implicita
C
      GRADNRM = GRAD_PARM(1,1)**2 + GRAD_PARM(1,2)**2
      IF(NDIM.EQ.3)GRADNRM = GRADNRM + GRAD_PARM(1,3)**2
C
      TBDIFF = TCB2/TPR1*REINV * GRADNRM
C
      RETURN
      END
C
