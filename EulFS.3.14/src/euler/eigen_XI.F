      SUBROUTINE Eigen_XI(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_XI.F,v 1.4 2020/03/28 09:51:14 abonfi Exp $
C
C     Conserved variables 
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
C
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
C
C
      INTEGER LDA,NDIM,LDB,I
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB)
C
      DOUBLE PRECISION DENS,DENSABARINV,DENSINVABARINV,TMPA,TRE
      PARAMETER(TRE=3.d0)
      INTEGER IDIM,ifail
C
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C     Conserved variables ..
C
      IF(ICHECK.EQ.0)THEN
         WRITE(6,*)' Calling eigen_XI with ICHECK = ',ICHECK 
         WRITE(6,*)' Doing nothing'
         RETURN
      ENDIF
C
      DENS = UAVG(1)
      DENSABARINV = DENS/ABAR
      DENSINVABARINV = ONE/(DENS*ABAR)
C
C     N.B. GRAD_CHAR  non serve, tranne quando ICHECK <> 0
C          cosi` la matrice, tranne che si utilizzi la
C          decomposizione numerica. 
C
C     here GRAD_CHAR is the gradient of the conserved variables
C
      CALL PARM2CONS(ZAVG,DVDZ,LDB,NDIM)
      CALL DGEMM('NoTranspose','NoTranspose',LDB,NDIM,LDB,
     +           ONE,DVDZ,LDB,GRAD_PARM,MAXNOFVAR,ZERO,GRAD_CHAR,LDW)
C
CREM  CALL R8Mat_Print('G',' ',LDB,NDIM,GRAD_CHAR,LDW,
CREM +      'CONSV variables within eigenXI ',IFAIL)
CREM  pause
C
C     Jacobian matrix in conserved variables
C
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
C
      Matrix(2,1,IDIM) = UAVG(3) * TMPA
      Matrix(2,2,IDIM) = UAVG(3) * GAM
      Matrix(2,3,IDIM) = UAVG(2) - GM1*UAVG(3)*UAVG(3)
      Matrix(2,4,IDIM) = -GM1 * UAVG(3)*UAVG(4)
      Matrix(2,5,IDIM) = -GM1 * UAVG(3)*UAVG(5)
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(3)+GM1*KINETIC
      Matrix(3,2,IDIM) = GM1
      Matrix(3,3,IDIM) = (TRE-GAM)*UAVG(3)
      Matrix(3,4,IDIM) = -GM1*UAVG(4)
      Matrix(3,5,IDIM) = -GM1*UAVG(5)
C
      Matrix(4,1,IDIM) = -UAVG(3)*UAVG(4)
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = UAVG(4)
      Matrix(4,4,IDIM) = UAVG(3)
      Matrix(4,5,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = -UAVG(3)*UAVG(5)
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = UAVG(5)
      Matrix(5,4,IDIM) = ZERO
      Matrix(5,5,IDIM) = UAVG(3)
C
      IDIM = 2
C
      Matrix(1,1,IDIM) = ZERO
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = ONE
      Matrix(1,5,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = UAVG(4) * TMPA
      Matrix(2,2,IDIM) = UAVG(4) * GAM
      Matrix(2,3,IDIM) = -GM1 * UAVG(4)*UAVG(3)
      Matrix(2,4,IDIM) = UAVG(2) - GM1*UAVG(4)*UAVG(4)
      Matrix(2,5,IDIM) = -GM1 * UAVG(4)*UAVG(5)
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(4)
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(4)
      Matrix(3,4,IDIM) = UAVG(3)
      Matrix(3,5,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = -UAVG(4)*UAVG(4)+GM1*KINETIC
      Matrix(4,2,IDIM) = GM1
      Matrix(4,3,IDIM) = -GM1*UAVG(3)
      Matrix(4,4,IDIM) = (TRE-GAM)*UAVG(4)
      Matrix(4,5,IDIM) = -GM1*UAVG(5)
C
      Matrix(5,1,IDIM) = -UAVG(4)*UAVG(5)
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = ZERO
      Matrix(5,4,IDIM) = UAVG(5)
      Matrix(5,5,IDIM) = UAVG(4)
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
C
      Matrix(2,1,IDIM) = UAVG(5) * TMPA
      Matrix(2,2,IDIM) = UAVG(5) * GAM
      Matrix(2,3,IDIM) = -GM1 * UAVG(5)*UAVG(3)
      Matrix(2,4,IDIM) = -GM1 * UAVG(5)*UAVG(4)
      Matrix(2,5,IDIM) = UAVG(2) - GM1*UAVG(5)*UAVG(5)
C
      Matrix(3,1,IDIM) = -UAVG(3)*UAVG(5)
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(5)
      Matrix(3,4,IDIM) = ZERO
      Matrix(3,5,IDIM) = UAVG(3)
C
      Matrix(4,1,IDIM) = -UAVG(4)*UAVG(5)
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = UAVG(5)
      Matrix(4,5,IDIM) = UAVG(4)
C
      Matrix(5,1,IDIM) = -UAVG(5)*UAVG(5)+GM1*KINETIC
      Matrix(5,2,IDIM) = GM1
      Matrix(5,3,IDIM) = -GM1*UAVG(3)
      Matrix(5,4,IDIM) = -GM1*UAVG(4)
      Matrix(5,5,IDIM) = (TRE-GAM)*UAVG(5)
C
C     When doing ALE, the components of the average grid velocity
C     have to be substracted from the diagonal entries of the
C     NDIM Jacobian matrices
C
      IF(LALE)THEN
         DO IDIM = 1,3
            TMPA = BAVG(IDIM)
            DO I = 1, 5
                  Matrix(I,I,IDIM) = Matrix(I,I,IDIM)-TMPA
            ENDDO ! I
         ENDDO ! IDIM
      ENDIF
C
      RETURN
      END
