      SUBROUTINE Eigen_VII(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_VII.F,v 1.10 2020/03/28 09:51:14 abonfi Exp $
C
C     Symmetrizing variables (dS,dp/(r*a),velocity) variables
C        as observed by E. van der Weide any choice of variables,
C        at least when using linear schemes, gives the same result
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
      INCLUDE 'constants.h'
C
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
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
      DOUBLE PRECISION DENS,DENSABARINV,DENSINVABARINV,Z1SQRINV
      INTEGER IDIM,ifail
C
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
      DENS = UAVG(1)
      DENSABARINV = DENS/ABAR
      DENSINVABARINV = ONE/(DENS*ABAR)
C
C     Advection speed for entropy ..
C
      R_SPEED(1,1) = UAVG(3)
      R_SPEED(2,1) = UAVG(4)
      IF(NDIM.EQ.3)R_SPEED(3,1) = UAVG(5)
C
C     Symmetrizing variables ..
C
      IF(ICHECK.NE.0)THEN 
C
C     N.B. GRAD_CHAR  non serve, tranne quando ICHECK <> 0
C          cosi` la matrice, tranne che si utilizzi la
C          decomposizione numerica. 
C
          CALL GRADPRIM(-10,NDIM,NDIM+2)
C
          DO 11 IDIM = 1 , NDIM
C
             GRAD_CHAR(1,IDIM) = GRAD_PRIM(1,IDIM) -
     &   GRAD_PRIM(2,IDIM)/ASQR
             GRAD_CHAR(2,IDIM) = GRAD_PRIM(2,IDIM)*DENSINVABARINV
             GRAD_CHAR(3,IDIM) = GRAD_PRIM(3,IDIM)
             GRAD_CHAR(4,IDIM) = GRAD_PRIM(4,IDIM)
             GRAD_CHAR(5,IDIM) = GRAD_PRIM(5,IDIM)
caldo    R_SPEED(IDIM,1) = UAVG(2+IDIM)
C
   11 CONTINUE
C
          IDIM = 1
C
          Matrix(1,1,IDIM) = UAVG(3)
          Matrix(1,2,IDIM) = ABAR
          Matrix(1,3,IDIM) = ZERO
          Matrix(1,4,IDIM) = ZERO
C
          Matrix(2,1,IDIM) = ABAR
          Matrix(2,2,IDIM) = UAVG(3)
          Matrix(2,3,IDIM) = ZERO
          Matrix(2,4,IDIM) = ZERO
C
          Matrix(3,1,IDIM) = ZERO
          Matrix(3,2,IDIM) = ZERO
          Matrix(3,3,IDIM) = UAVG(3)
          Matrix(3,4,IDIM) = ZERO
C
          Matrix(4,1,IDIM) = ZERO
          Matrix(4,2,IDIM) = ZERO
          Matrix(4,3,IDIM) = ZERO
          Matrix(4,4,IDIM) = UAVG(3)
C
          IDIM = 2
C
          Matrix(1,1,IDIM) = UAVG(4)
          Matrix(1,2,IDIM) = ZERO
          Matrix(1,3,IDIM) = ABAR
          Matrix(1,4,IDIM) = ZERO
C
          Matrix(2,1,IDIM) = ZERO
          Matrix(2,2,IDIM) = UAVG(4)
          Matrix(2,3,IDIM) = ZERO
          Matrix(2,4,IDIM) = ZERO
C
          Matrix(3,1,IDIM) = ABAR
          Matrix(3,2,IDIM) = ZERO
          Matrix(3,3,IDIM) = UAVG(4)
          Matrix(3,4,IDIM) = ZERO
C
          Matrix(4,1,IDIM) = ZERO
          Matrix(4,2,IDIM) = ZERO
          Matrix(4,3,IDIM) = ZERO
          Matrix(4,4,IDIM) = UAVG(4)
C
          IDIM = 3
C
          Matrix(1,1,IDIM) = UAVG(5)
          Matrix(1,2,IDIM) = ZERO
          Matrix(1,3,IDIM) = ZERO
          Matrix(1,4,IDIM) = ABAR
C
          Matrix(2,1,IDIM) = ZERO
          Matrix(2,2,IDIM) = UAVG(5)
          Matrix(2,3,IDIM) = ZERO
          Matrix(2,4,IDIM) = ZERO
C
          Matrix(3,1,IDIM) = ZERO
          Matrix(3,2,IDIM) = ZERO
          Matrix(3,3,IDIM) = UAVG(5)
          Matrix(3,4,IDIM) = ZERO
C
          Matrix(4,1,IDIM) = ABAR
          Matrix(4,2,IDIM) = ZERO
          Matrix(4,3,IDIM) = ZERO
          Matrix(4,4,IDIM) = UAVG(5)
C
      ENDIF
C
C     Transformation matrix from symmetrizing to conserved variables ..
C
      DUDV(1,1) = ONE
      DUDV(2,1) = KINETIC
      DUDV(3,1) = UAVG(3)
      DUDV(4,1) = UAVG(4)
C
      DUDV(1,2) = DENSABARINV
      DUDV(2,2) = DENSABARINV*UAVG(2)
      DUDV(3,2) = DENSABARINV*UAVG(3)
      DUDV(4,2) = DENSABARINV*UAVG(4)
C
      DUDV(1,3) = ZERO
      DUDV(2,3) = DENS*UAVG(3)
      DUDV(3,3) = DENS
      DUDV(4,3) = ZERO
C
      DUDV(1,4) = ZERO
      DUDV(2,4) = DENS*UAVG(4)
      DUDV(3,4) = ZERO
      DUDV(4,4) = DENS
C
      IF(NDIM.EQ.3)THEN
         DUDV(5,1) = UAVG(5)
         DUDV(5,2) = DENSABARINV*UAVG(5)
         DUDV(5,3) = ZERO
         DUDV(5,4) = ZERO
C
         DUDV(1,5) = ZERO
         DUDV(2,5) = DENS*UAVG(5)
         DUDV(3,5) = ZERO
         DUDV(4,5) = ZERO
         DUDV(5,5) = DENS
      ENDIF
C
      Z1SQRINV = ONE/(ZAVG(1)*ZAVG(1))
C
C     Transformation matrix from parameter vector to symmetrizing variables
C
      DVDZ(2,1) = GM1OG * DENSINVABARINV * ZAVG(2)
      DVDZ(1,1) = TWO*ZAVG(1) - GM1OG * ZAVG(2) / ASQR
      DVDZ(3,1) = - ZAVG(3) * Z1SQRINV
      DVDZ(4,1) = - ZAVG(4) * Z1SQRINV
C
      DVDZ(2,2) = GM1OG * ZAVG(1) * DENSINVABARINV
      DVDZ(1,2) = - GM1OG * ZAVG(1) / ASQR
      DVDZ(3,2) = ZERO
      DVDZ(4,2) = ZERO
C
      DVDZ(2,3) = - GM1OG * ZAVG(3) * DENSINVABARINV
      DVDZ(1,3) = GM1OG * ZAVG(3) / ASQR
      DVDZ(3,3) = ONE/ZAVG(1)
      DVDZ(4,3) = ZERO
C
      DVDZ(2,4) = - GM1OG * ZAVG(4) * DENSINVABARINV
      DVDZ(1,4) = GM1OG * ZAVG(4) / ASQR
      DVDZ(3,4) = ZERO
      DVDZ(4,4) = ONE/ZAVG(1)
C
      IF(NDIM.EQ.3)THEN
      DVDZ(5,1) = - ZAVG(5) * Z1SQRINV
      DVDZ(5,2) = ZERO
      DVDZ(5,3) = ZERO
      DVDZ(5,4) = ZERO
C
      DVDZ(2,5) = - GM1OG * ZAVG(5) * DENSINVABARINV
      DVDZ(1,5) = GM1OG * ZAVG(5) / ASQR
      DVDZ(3,5) = ZERO
      DVDZ(4,5) = ZERO
      DVDZ(5,5) = ONE/ZAVG(1)
      ENDIF
C
#ifdef FSPL_USE_DEBUG
C
C    Tests two ways of computing the gradient of simmetrysing 
C    variables 
C
            CALL DGEMM('N','N',NDIM+2,NDIM,NDIM+2,ONE,
     +      DVDZ,NDIM+2,GRAD_PARM,NMAX,0.D0,WORK,NDIM+2)
            CALL MATSUB(WORK,NDIM+2,GRAD_CHAR,LDW,
     +      NDIM+2,NDIM)
            IF( .NOT. NULLMAT(WORK,NDIM+2,NDIM+2,NDIM+2,1.D-14) )THEN
            CALL R8Mat_Print('G',' ',NDIM+2,NDIM,WORK,NDIM+2,
     +      'dZdU*dZ matrix ????',IFAIL)
            CALL R8Mat_Print('G',' ',NDIM+2,NDIM,GRAD_CHAR,LDW,
     +      'dW matrix ????',IFAIL)
            PAUSE
            ENDIF
#endif
C
      RETURN
      END
