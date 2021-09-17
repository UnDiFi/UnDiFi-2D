      SUBROUTINE GRADPRIM(IELEM,NDIM,NOFVAR)
C
      IMPLICIT NONE 
C
C THIS SUBROUTINE COMPUTES:
C c) The gradient of the primitive variables GRAD_PRIM(1:NOFVAR,1:NOFVERT) 
C    (in cart. coord.)
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
      INCLUDE 'pfcgas.com'
C
C     .. Scalar Arguments ..
C
      INTEGER IELEM,NDIM,NOFVAR
C
C     .. Local Scalars ..
C
      INTEGER I,IVAR,JVERT
      DOUBLE PRECISION SUM,DUM,DUMSQR,TEMP
C
C     .. Local Arrays ..
C
      DOUBLE PRECISION AMAT(5,5)
C
C     .. External Functions ..
C
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
      DATA (AMAT(1,I),I=1,5) / ONE,ZERO,ZERO,ZERO,ZERO /
      DATA (AMAT(2,I),I=1,5) / ZERO,ONE,ZERO,ZERO,ZERO /
      DATA (AMAT(3,I),I=1,5) / ZERO,ZERO,ONE,ZERO,ZERO /
      DATA (AMAT(4,I),I=1,5) / ZERO,ZERO,ZERO,ONE,ZERO /
      DATA (AMAT(5,I),I=1,5) / ZERO,ZERO,ZERO,ZERO,ONE /
C
C     .. Executable Statements ..
C
C     be careful because GRADPRIM must be
C     called with NOFVAR being equal to the
C     number of flow variables
C
C
C
C *********************************************************************
C ASSEMBLING THE LOWER TRIANGULAR MATRIX WHICH RELATES
C THE GRADIENT OF THE PARAMETER VECTOR TO THE GRADIENT
C OF THE PRIMITIVE VARIABLES( density,static pressure,velocities ).
C *********************************************************************
C
      DUM = ONE / ZAVG(1)
      DUMSQR = DUM * DUM
c
c	.. First row ..
c
      AMAT(1,1) = TWO * ZAVG(1)
      AMAT(1,2) = ZERO
      AMAT(1,3) = ZERO
      AMAT(1,4) = ZERO
      AMAT(1,5) = ZERO
c
c	.. Second row ..
c
      AMAT(2,1) = GM1OG * ZAVG(2)
      AMAT(2,2) = GM1OG * ZAVG(1)
      AMAT(2,3) =-GM1OG * ZAVG(3)
      AMAT(2,4) =-GM1OG * ZAVG(4)
      AMAT(2,5) =-GM1OG * ZAVG(5)
c
c	.. Third row ..
c
      AMAT(3,1) = - ZAVG(3) * DUMSQR
      AMAT(3,2) = ZERO
      AMAT(3,3) = DUM
      AMAT(3,4) = ZERO
      AMAT(3,5) = ZERO
c
c	.. Fourth row ..
c
      AMAT(4,1) = - ZAVG(4) * DUMSQR
      AMAT(4,2) = ZERO
      AMAT(4,3) = ZERO
      AMAT(4,4) = DUM
      AMAT(4,5) = ZERO
c
c	.. Fifth row ..
c
      AMAT(5,1) = - ZAVG(5) * DUMSQR
      AMAT(5,2) = ZERO
      AMAT(5,3) = ZERO
      AMAT(5,4) = ZERO
      AMAT(5,5) = DUM
C
C *********************************************************************
C PRODUCT OF THE MATRIX TIMES THE GRADIENT OF THE PARAMETER VECTOR
C *********************************************************************
C
   25 CONTINUE
C
      CALL DGEMM( 'N', 'N', NOFVAR, NDIM, NOFVAR, ONE, AMAT, 5,
     &GRAD_PARM, NMAX, ZERO, GRAD_PRIM, MAXNOFEQN )
C
      RETURN
      END
