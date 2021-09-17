      SUBROUTINE Eigen_VIIsa(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_VII.F,v 1.6 2002/09/13 08:42:01 abonfi Exp $
C
C     Symmetrizing variables (dS,dp/(r*a),velocity,dn_t) variables
C     Spalart-Allmaras model
C
C
C     define FSPL_USE_DEBUG for debugging purposes
C
C#define FSPL_USE_DEBUG
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
      INTEGER NVAR,LDC,I
      PARAMETER(LDC=6)
      DOUBLE PRECISION WORK(NMAX*NMAX),AMAT(LDC,LDC,3),
     +Jacobian(LDC,LDC,3)
      LOGICAL NULLMAT
#endif
C
      DOUBLE PRECISION DENS,DENSABARINV,DENSINVABARINV,Z1SQRINV
      INTEGER IDIM,ifail
C
      INTEGER ITURB
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C     location where the turbulence variable is stored
C
      ITURB = NDIM + 3 
      UAVG(ITURB) = ZAVG(ITURB)/ZAVG(1)
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
         IF(NDIM.EQ.3)GRAD_CHAR(5,IDIM) = GRAD_PRIM(5,IDIM)
c
c    dn_t = -z_1^{-2}z_6 \nabla z_1 + z_1^{-1} \nabla z_6 
c
         GRAD_CHAR(ITURB,IDIM) = -ZAVG(ITURB)*GRAD_PARM(1,IDIM)/DENS +
     >   GRAD_PARM(ITURB,IDIM)/ZAVG(1)
C
   11 CONTINUE
C
         IDIM = 1
C
C     matrix for (dp/ra,du,dv,dw)
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
#ifdef FSPL_USE_DEBUG
C
C     matrix for (dS,dp/ra,du,dv,dw,dn_t)
C
         NVAR = NDIM + 3
         CALL DINIT(3*LDC*LDC,ZERO,AMAT,1)
         DO IDIM = 1,NDIM
            DO I = 1,NVAR
               Amat(I,I,IDIM) = UAVG(IDIM+2)
            ENDDO
         ENDDO
         DO IDIM = 1,NDIM
            I = 2 + IDIM
            Amat(2,I,IDIM) = ABAR
            Amat(I,2,IDIM) = ABAR
         ENDDO
C
c
c     JACOBIAN IN Parameter Vector (Compressible)
c
c	Matrix dF/dZ
c
      Jacobian(1,1,1) = ZAVG(3)
      Jacobian(1,2,1) = ZERO
      Jacobian(1,3,1) = ZAVG(1)
      Jacobian(1,4,1) = ZERO
c
      Jacobian(2,1,1) = ZERO
      Jacobian(2,2,1) = ZAVG(3)
      Jacobian(2,3,1) = ZAVG(2)
      Jacobian(2,4,1) = ZERO
c
      Jacobian(3,1,1) = GM1OG * ZAVG(2)
      Jacobian(3,2,1) = GM1OG * ZAVG(1)
      Jacobian(3,3,1) = GP1OG * ZAVG(3)
      Jacobian(3,4,1) =-GM1OG * ZAVG(4)
c
      Jacobian(4,1,1) = ZERO
      Jacobian(4,2,1) = ZERO
      Jacobian(4,3,1) = ZAVG(4)
      Jacobian(4,4,1) = ZAVG(3)
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,5,1) = ZERO
         Jacobian(2,5,1) = ZERO
         Jacobian(3,5,1) =-GM1OG * ZAVG(5)
         Jacobian(4,5,1) = ZERO
c
         Jacobian(5,1,1) = ZERO
         Jacobian(5,2,1) = ZERO
         Jacobian(5,3,1) = ZAVG(5)
         Jacobian(5,4,1) = ZERO
         Jacobian(5,5,1) = ZAVG(3)
      ENDIF
c
      Jacobian(1,NVAR,1) = ZERO
      Jacobian(2,NVAR,1) = ZERO
      Jacobian(3,NVAR,1) = ZERO
      Jacobian(4,NVAR,1) = ZERO
      IF(NDIM.EQ.3)Jacobian(5,NVAR,1) = ZERO
c
      Jacobian(NVAR,1,1) = ZERO
      Jacobian(NVAR,2,1) = ZERO
      Jacobian(NVAR,3,1) = ZAVG(NVAR)
      Jacobian(NVAR,4,1) = ZERO
      IF(NDIM.EQ.3)Jacobian(NVAR,5,1) = ZERO
      Jacobian(NVAR,NVAR,1) = ZAVG(3)
c
c	MATRIX dG/dZ
c
      Jacobian(1,1,2) = ZAVG(4)
      Jacobian(1,2,2) = ZERO
      Jacobian(1,3,2) = ZERO
      Jacobian(1,4,2) = ZAVG(1)
c
      Jacobian(2,1,2) = ZERO
      Jacobian(2,2,2) = ZAVG(4)
      Jacobian(2,3,2) = ZERO
      Jacobian(2,4,2) = ZAVG(2)
c
      Jacobian(3,1,2) = ZERO
      Jacobian(3,2,2) = ZERO
      Jacobian(3,3,2) = ZAVG(4)
      Jacobian(3,4,2) = ZAVG(3)
c
      Jacobian(4,1,2) = GM1OG * ZAVG(2)
      Jacobian(4,2,2) = GM1OG * ZAVG(1)
      Jacobian(4,3,2) =-GM1OG * ZAVG(3)
      Jacobian(4,4,2) = GP1OG * ZAVG(4)
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,5,2) = ZERO
         Jacobian(2,5,2) = ZERO
         Jacobian(3,5,2) = ZERO
         Jacobian(4,5,2) =-GM1OG * ZAVG(5)
         Jacobian(5,1,2) = ZERO
         Jacobian(5,2,2) = ZERO
         Jacobian(5,3,2) = ZERO
         Jacobian(5,4,2) = ZAVG(5)
         Jacobian(5,5,2) = ZAVG(4)
      ENDIF
c
      Jacobian(1,NVAR,2) = ZERO
      Jacobian(2,NVAR,2) = ZERO
      Jacobian(3,NVAR,2) = ZERO
      Jacobian(4,NVAR,2) = ZERO
      IF(NDIM.EQ.3)Jacobian(5,NVAR,2) = ZERO
c
      Jacobian(NVAR,1,2) = ZERO
      Jacobian(NVAR,2,2) = ZERO
      Jacobian(NVAR,3,2) = ZERO
      Jacobian(NVAR,4,2) = ZAVG(NVAR)
      IF(NDIM.EQ.3)Jacobian(NVAR,5,2) = ZERO
      Jacobian(NVAR,NVAR,2) = ZAVG(4)
c
c	Matrix dH/dZ
c
      IF(NDIM.EQ.3)THEN
         Jacobian(1,1,3) = ZAVG(5)
         Jacobian(1,2,3) = ZERO
         Jacobian(1,3,3) = ZERO
         Jacobian(1,4,3) = ZERO
         Jacobian(1,5,3) = ZAVG(1)
c
         Jacobian(2,1,3) = ZERO
         Jacobian(2,2,3) = ZAVG(5)
         Jacobian(2,3,3) = ZERO
         Jacobian(2,4,3) = ZERO
         Jacobian(2,5,3) = ZAVG(2)
c
         Jacobian(3,1,3) = ZERO
         Jacobian(3,2,3) = ZERO
         Jacobian(3,3,3) = ZAVG(5)
         Jacobian(3,4,3) = ZERO
         Jacobian(3,5,3) = ZAVG(3)
c
         Jacobian(4,1,3) = ZERO
         Jacobian(4,2,3) = ZERO
         Jacobian(4,3,3) = ZERO
         Jacobian(4,4,3) = ZAVG(5)
         Jacobian(4,5,3) = ZAVG(4)
c
         Jacobian(5,1,3) = GM1OG * ZAVG(2)
         Jacobian(5,2,3) = GM1OG * ZAVG(1)
         Jacobian(5,3,3) =-GM1OG * ZAVG(3)
         Jacobian(5,4,3) =-GM1OG * ZAVG(4)
         Jacobian(5,5,3) = GP1OG * ZAVG(5)
c
         Jacobian(1,NVAR,3) = ZERO
         Jacobian(2,NVAR,3) = ZERO
         Jacobian(3,NVAR,3) = ZERO
         Jacobian(4,NVAR,3) = ZERO
         Jacobian(5,NVAR,3) = ZERO
c
         Jacobian(NVAR,1,3) = ZERO
         Jacobian(NVAR,2,3) = ZERO
         Jacobian(NVAR,3,3) = ZERO
         Jacobian(NVAR,4,3) = ZERO
         Jacobian(NVAR,5,3) = ZAVG(NVAR)
         Jacobian(NVAR,NVAR,3) = ZAVG(5)
C
      ENDIF
#endif
C
C     Transformation matrix from symmetrizing to conserved variables ..
C
      DUDV(1,1) = ONE
      DUDV(2,1) = KINETIC
      DUDV(3,1) = UAVG(3)
      DUDV(4,1) = UAVG(4)
      DUDV(ITURB,1) = UAVG(ITURB)
C
      DUDV(1,2) = DENSABARINV
      DUDV(2,2) = DENSABARINV*UAVG(2)
      DUDV(3,2) = DENSABARINV*UAVG(3)
      DUDV(4,2) = DENSABARINV*UAVG(4)
      DUDV(ITURB,2) = DENSABARINV*UAVG(ITURB)
C
      DUDV(1,3) = ZERO
      DUDV(2,3) = DENS*UAVG(3)
      DUDV(3,3) = DENS
      DUDV(4,3) = ZERO
      DUDV(ITURB,3) = ZERO
C
      DUDV(1,4) = ZERO
      DUDV(2,4) = DENS*UAVG(4)
      DUDV(3,4) = ZERO
      DUDV(4,4) = DENS
      DUDV(ITURB,4) = ZERO
C
      DUDV(1,ITURB) = ZERO
      DUDV(2,ITURB) = ZERO
      DUDV(3,ITURB) = ZERO
      DUDV(4,ITURB) = ZERO
      DUDV(ITURB,ITURB) = DENS
C
C
      IF(NDIM.EQ.3)THEN
         DUDV(5,1) = UAVG(5)
         DUDV(5,2) = DENSABARINV*UAVG(5)
         DUDV(5,3) = ZERO
         DUDV(5,4) = ZERO
         DUDV(5,5) = ZERO
C
         DUDV(1,5) = ZERO
         DUDV(2,5) = DENS*UAVG(5)
         DUDV(3,5) = ZERO
         DUDV(4,5) = ZERO
         DUDV(5,5) = DENS
C
         DUDV(6,5) = ZERO
         DUDV(5,6) = ZERO
      ENDIF
C
C
      Z1SQRINV = ONE/(ZAVG(1)*ZAVG(1))
C
C     Transformation matrix from parameter vector to 
C     symmetrizing variables
C
      DVDZ(1,1) = TWO*ZAVG(1) - GM1OG * ZAVG(2) / ASQR
      DVDZ(2,1) = GM1OG * DENSINVABARINV * ZAVG(2)
      DVDZ(3,1) = - ZAVG(3) * Z1SQRINV
      DVDZ(4,1) = - ZAVG(4) * Z1SQRINV
      DVDZ(ITURB,1) = - ZAVG(ITURB) * Z1SQRINV
C
      DVDZ(1,2) = - GM1OG * ZAVG(1) / ASQR
      DVDZ(2,2) = GM1OG * ZAVG(1) * DENSINVABARINV
      DVDZ(3,2) = ZERO
      DVDZ(4,2) = ZERO
      DVDZ(ITURB,2) = ZERO
C
      DVDZ(1,3) = GM1OG * ZAVG(3) / ASQR
      DVDZ(2,3) = - GM1OG * ZAVG(3) * DENSINVABARINV
      DVDZ(3,3) = ONE/ZAVG(1)
      DVDZ(4,3) = ZERO
      DVDZ(ITURB,3) = ZERO
C
      DVDZ(1,4) = GM1OG * ZAVG(4) / ASQR
      DVDZ(2,4) = - GM1OG * ZAVG(4) * DENSINVABARINV
      DVDZ(3,4) = ZERO
      DVDZ(4,4) = ONE/ZAVG(1)
      DVDZ(ITURB,4) = ZERO
C
      DVDZ(1,ITURB) = ZERO
      DVDZ(2,ITURB) = ZERO
      DVDZ(3,ITURB) = ZERO
      DVDZ(4,ITURB) = ZERO
      DVDZ(ITURB,ITURB) = ONE/ZAVG(1)
C
      IF(NDIM.EQ.3)THEN
         DVDZ(5,1) = - ZAVG(5) * Z1SQRINV
         DVDZ(5,2) = ZERO
         DVDZ(5,3) = ZERO
         DVDZ(5,4) = ZERO
C
         DVDZ(1,5) = GM1OG * ZAVG(5) / ASQR
         DVDZ(2,5) = - GM1OG * ZAVG(5) * DENSINVABARINV
         DVDZ(3,5) = ZERO
         DVDZ(4,5) = ZERO
         DVDZ(5,5) = ONE/ZAVG(1)
C
         DVDZ(6,5) = ZERO
         DVDZ(5,6) = ZERO
      ENDIF
C
#ifdef FSPL_USE_DEBUG
C
C    Tests two ways of computing the gradient of simmetrysing 
C    variables 
C
      NVAR = NDIM + 3
caldo CALL X04CAF('G',' ',NVAR,NDIM,GRAD_PARM,NMAX,
caldo+      'grad(Z) matrix ',IFAIL)
C
            CALL DGEMM('N','N',NVAR,NDIM,NVAR,ONE,
     +      DVDZ,LDB,GRAD_PARM,NMAX,0.D0,WORK,NVAR)
caldo       CALL X04CAF('G',' ',NVAR,NDIM,WORK,NVAR,
caldo+      'dVdZ*grad(Z) matrix ????',IFAIL)
            CALL MATSUB(WORK,NVAR,GRAD_CHAR,LDW,
     +      NVAR,NDIM)
            IF( .NOT. NULLMAT(WORK,NVAR,NDIM,NVAR,1.D-14) )THEN
            CALL X04CAF('G',' ',NVAR,NDIM,WORK,NVAR,
     +      'This should be the null matrix !!!!',IFAIL)
            CALL X04CAF('G',' ',NVAR,NDIM,GRAD_CHAR,LDW,
     +      'dW matrix ????',IFAIL)
            PAUSE
            ENDIF
C
C     dUdV \tilde{A} dVdZ = F_x
C
            DO IDIM = 1,NDIM
            CALL DGEMM('N','N',NVAR,NVAR,NVAR,ONE,
     +      DUDV,LDB,AMAT(1,1,IDIM),LDC,0.D0,WORK,NVAR)
            CALL DGEMM('N','N',NVAR,NVAR,NVAR,ONE,
     +      WORK,NVAR,DVDZ,LDB,0.D0,AMAT(1,1,IDIM),LDC)
      CALL X04CAF('G',' ',NVAR,NVAR,AMAT(1,1,IDIM),LDC,
     +      'F,G,H matrix ',IFAIL)
      CALL X04CAF('G',' ',NVAR,NVAR,Jacobian(1,1,IDIM),LDC,
     +      'F_Z,G_Z,H_Z matrix ',IFAIL)
      ENDDO
      pause
#endif
C
      RETURN
      END
