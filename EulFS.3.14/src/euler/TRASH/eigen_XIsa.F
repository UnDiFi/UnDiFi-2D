      SUBROUTINE Eigen_XIsa(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_XI.F,v 1.1 2003/07/17 13:59:00 abonfi Exp $
C
C     Conserved variables + SA model
C
C     define FSPL_USE_DEBUG for debugging purposes, i.e. when using
C     -check
C
C
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants'
C
      INCLUDE 'three'
      INCLUDE 'flags.com'
C
C
      INTEGER LDA,NDIM,LDB
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB)
#ifdef FSPL_USE_DEBUG
      DOUBLE PRECISION WORK(NMAX*NMAX)
      LOGICAL NULLMAT
#endif
C
      DOUBLE PRECISION TMPA
      INTEGER IDIM,ifail
C
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C
C     Conserved variables ..
C
      IF(ICHECK.NE.0)THEN 
C
C     N.B. GRAD_CHAR  non serve, tranne quando ICHECK <> 0
C          cosi` la matrice, tranne che si utilizzi la
C          decomposizione numerica. 
C
C     here GRAD_CHAR are the conserved variables
C
      CALL PARM2CONSSA(ZAVG,DVDZ,LDB,NDIM)
      CALL DGEMM('NoTranspose','NoTranspose',LDB,NDIM,LDB,
     +           ONE,DVDZ,LDB,GRAD_PARM,NMAX,ZERO,GRAD_CHAR,LDW)
CREM  CALL X04CAF('G',' ',LDB,NDIM,GRAD_CHAR,LDW,
CREM +      'CONSV variables within eigenXI ',IFAIL)
CREM  pause
C
C     Jacobian matrix in conserved variables
C
C     the matrix is dF/dU ; dG/dU ; dH/dU
C
C
      TMPA = GM1*KINETIC-UAVG(2)
C
      IDIM = 1
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ONE
      Matrix(1,4,IDIM) = ZERO
      Matrix(1,5,IDIM) = ZERO
      Matrix(1,6,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = UAVG(3) * TMPA
      Matrix(2,2,IDIM) = UAVG(3) * GAM
      Matrix(2,3,IDIM) = UAVG(2) - GM1*UAVG(3)*UAVG(3)
      Matrix(2,4,IDIM) = -GM1 * UAVG(3)*UAVG(4)
      Matrix(2,5,IDIM) = -GM1 * UAVG(3)*UAVG(5)
      Matrix(2,6,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(3)+GM1*KINETIC
      Matrix(3,2,IDIM) = GM1
      Matrix(3,3,IDIM) = (3.d0-GAM)*UAVG(3)
      Matrix(3,4,IDIM) = -GM1*UAVG(4)
      Matrix(3,5,IDIM) = -GM1*UAVG(5)
      Matrix(3,6,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = -UAVG(3)*UAVG(4)
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = UAVG(4)
      Matrix(4,4,IDIM) = UAVG(3)
      Matrix(4,5,IDIM) = ZERO
      Matrix(4,6,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = -UAVG(3)*UAVG(5)
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = UAVG(5)
      Matrix(5,4,IDIM) = ZERO
      Matrix(5,5,IDIM) = UAVG(3)
      Matrix(5,6,IDIM) = ZERO
C
      Matrix(6,1,IDIM) = -UAVG(3)*UAVG(6)
      Matrix(6,2,IDIM) = ZERO
      Matrix(6,3,IDIM) = UAVG(6)
      Matrix(6,4,IDIM) = ZERO
      Matrix(6,5,IDIM) = ZERO
      Matrix(6,6,IDIM) = UAVG(3)
C
      IDIM = 2
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = ONE
      Matrix(1,5,IDIM) = ZERO
      Matrix(1,6,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = UAVG(4) * TMPA
      Matrix(2,2,IDIM) = UAVG(4) * GAM
      Matrix(2,3,IDIM) = -GM1 * UAVG(4)*UAVG(3)
      Matrix(2,4,IDIM) = UAVG(2) - GM1*UAVG(4)*UAVG(4)
      Matrix(2,5,IDIM) = -GM1 * UAVG(4)*UAVG(5)
      Matrix(2,6,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(4)
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(4)
      Matrix(3,4,IDIM) = UAVG(3)
      Matrix(3,5,IDIM) = ZERO
      Matrix(3,6,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = -UAVG(4)*UAVG(4)+GM1*KINETIC
      Matrix(4,2,IDIM) = GM1
      Matrix(4,3,IDIM) = -GM1*UAVG(3)
      Matrix(4,4,IDIM) = (3.d0-GAM)*UAVG(4)
      Matrix(4,5,IDIM) = -GM1*UAVG(5)
      Matrix(4,6,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = -UAVG(4)*UAVG(5)
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = ZERO
      Matrix(5,4,IDIM) = UAVG(5)
      Matrix(5,5,IDIM) = UAVG(4)
      Matrix(5,6,IDIM) = ZERO
C
      Matrix(6,1,IDIM) = -UAVG(4)*UAVG(6)
      Matrix(6,2,IDIM) = ZERO
      Matrix(6,3,IDIM) = ZERO
      Matrix(6,4,IDIM) = UAVG(6)
      Matrix(6,5,IDIM) = ZERO
      Matrix(6,6,IDIM) = UAVG(4)
C
C
C
      IDIM = 3
C
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = ZERO
      Matrix(1,5,IDIM) = ONE
      Matrix(1,6,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = UAVG(5) * TMPA
      Matrix(2,2,IDIM) = UAVG(5) * GAM
      Matrix(2,3,IDIM) = -GM1 * UAVG(5)*UAVG(3)
      Matrix(2,4,IDIM) = -GM1 * UAVG(5)*UAVG(4)
      Matrix(2,5,IDIM) = UAVG(2) - GM1*UAVG(5)*UAVG(5)
      Matrix(2,6,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(5)
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(5)
      Matrix(3,4,IDIM) = ZERO
      Matrix(3,5,IDIM) = UAVG(3)
      Matrix(3,6,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = -UAVG(4)*UAVG(5)
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = UAVG(5)
      Matrix(4,5,IDIM) = UAVG(4)
      Matrix(4,6,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = -UAVG(5)*UAVG(5)+GM1*KINETIC
      Matrix(5,2,IDIM) = GM1
      Matrix(5,3,IDIM) = -GM1*UAVG(3)
      Matrix(5,4,IDIM) = -GM1*UAVG(4)
      Matrix(5,5,IDIM) = (3.d0-GAM)*UAVG(5)
      Matrix(5,6,IDIM) = ZERO
C
      Matrix(6,1,IDIM) = -UAVG(5)*UAVG(6)
      Matrix(6,2,IDIM) = ZERO
      Matrix(6,3,IDIM) = ZERO
      Matrix(6,4,IDIM) = ZERO
      Matrix(6,5,IDIM) = UAVG(6)
      Matrix(6,6,IDIM) = UAVG(5)
C
C
C
C
      ENDIF
C
C
      RETURN
      END
