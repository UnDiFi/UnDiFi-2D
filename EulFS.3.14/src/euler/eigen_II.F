      SUBROUTINE Eigen_II(Matrix,LDA,DVDZ,DUDV,NDIM,NOFVAR)
C
C#define DEBUG
C
C     van Leer, Lee, Roe preconditioner in conservative form ..
C
      IMPLICIT NONE
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'three.com'
      INCLUDE 'blkEulerII.com'
      INCLUDE 'frame.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
C
      INTEGER LDA,NDIM,NOFVAR
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(NOFVAR,NOFVAR),
     +DUDV(NOFVAR,NOFVAR)
C
      DOUBLE PRECISION MSQRM1,CNST,QBAR,ZDOTN,ZDOTT,ZDOTS
      DOUBLE PRECISION DensAbar,TEMP1,TEMP2,TEMP3,TEMP4,
     +TEMP0,TEMPA,TEMPX,TEMPY,TEMPZ,DENSABARINV,DENSINVABARINV
      INTEGER IDIM,JDIM,IFAIL
      logical nullmat
C
      DOUBLE PRECISION U_dot_DU(3),S_dot_DU(3),WORK(4,4)
C
C
      DOUBLE PRECISION DDOT,fun_BETA_aldo
      EXTERNAL         DDOT,fun_BETA_aldo
C
C
C              .. Mach related variables ..
C
      MSQRM1       = MACHSQR - ONE
caldo BETA        = SQRT( MAX( EPSQR_SONIC , ABS( MSQRM1 ) ) )
      BETA       = fun_BETA_aldo( MSQRM1 )
      BETASQR       = BETA*BETA
      X              = BETA / MAX( MACH , ONE )
      nu_p       = HALF * ( MSQRM1/BETASQR + ONE )
      nu_m       = HALF * ( MSQRM1/BETASQR - ONE )
C
      QBAR = ONE / QINV
      DensAbar = UAVG(1)*ABAR
      TEMP1 = BETA/DensAbar
      DENSABARINV = UAVG(1)/ABAR
      DENSINVABARINV = ONE/DENSABAR
C
C     unroll !!!
C
      CALL GRADPRIM(-10,NDIM,NDIM+2)
C
         U_dot_dU(1) = DDOT(NDIM,GRAD_PRIM(3,1),1,
     1                         RotationMatrix(1,1),1)
         U_dot_dU(2) = DDOT(NDIM,GRAD_PRIM(3,2),1,
     1                         RotationMatrix(1,1),1)
         S_dot_dU(1) = DDOT(NDIM,GRAD_PRIM(3,1),1,
     1                         RotationMatrix(1,2),1)
         S_dot_dU(2) = DDOT(NDIM,GRAD_PRIM(3,2),1,
     1                         RotationMatrix(1,2),1)
C
C     Transformation matrix from parameter vector 
C     to characteristic variables
C
      ZDOTN =(RotationMatrix(1,1) * ZAVG(3) +
     1        RotationMatrix(2,1) * ZAVG(4) )/UAVG(1)
C
      ZDOTS =(RotationMatrix(1,2) * ZAVG(3) +
     1        RotationMatrix(2,2) * ZAVG(4) )/UAVG(1)
C
C
      TEMP0 = GM1OG/ASQR
      TEMP2 = GM1OG*DENSINVABARINV
      TEMPX = TEMP2 * ZAVG(3)
      TEMPY = TEMP2 * ZAVG(4)
      TEMPZ = TEMP2 * ZAVG(5)
      TEMPA = TWO*ZAVG(1) - TEMP0 * ZAVG(2)
C
      dVdZ(1,1) = TWO*ZAVG(1)-TEMP0*ZAVG(2)
      dVdZ(1,2) = -TEMP0*ZAVG(1)
      dVdZ(1,3) =  TEMP0*ZAVG(3)
      dVdZ(1,4) =  TEMP0*ZAVG(4)
C
      dVdZ(2,1) = TEMP2 * ZAVG(2) - MACH * ZDOTN
      dVdZ(2,2) = TEMP2 * ZAVG(1)
      dVdZ(2,3) =-TEMPX + MACH*RotationMatrix(1,1)/ZAVG(1)
      dVdZ(2,4) =-TEMPY + MACH*RotationMatrix(2,1)/ZAVG(1)
C
      dVdZ(3,1) = BETA * TEMP2 * ZAVG(2) - MACH * ZDOTS
      dVdZ(3,2) = BETA * dVdZ(2,2)
      dVdZ(3,3) =-BETA * TEMPX + MACH * RotationMatrix(1,2) / ZAVG(1)
      dVdZ(3,4) =-BETA * TEMPY + MACH * RotationMatrix(2,2) / ZAVG(1)
C
      dVdZ(4,1) = BETA * TEMP2 * ZAVG(2) + MACH * ZDOTS
      dVdZ(4,2) = dVdZ(3,2)
      dVdZ(4,3) =-BETA * TEMPX - MACH * RotationMatrix(1,2) / ZAVG(1)
      dVdZ(4,4) =-BETA * TEMPY - MACH * RotationMatrix(2,2) / ZAVG(1)
C
      dVdZ(1,5) =  TEMP0*ZAVG(5)
      dVdZ(2,5) =-TEMPZ + MACH*RotationMatrix(3,1)/ZAVG(1)
      dVdZ(3,5) =-BETA * TEMPZ + MACH * RotationMatrix(3,2) / ZAVG(1)
      dVdZ(4,5) =-BETA * TEMPZ - MACH * RotationMatrix(3,2) / ZAVG(1)
C
C     Advection speeds ..
C
         R_SPEED(1,1) = RotationMatrix(1,1)
         R_SPEED(1,2) = RotationMatrix(1,1)
         R_SPEED(1,3) = X*(nu_p*RotationMatrix(1,1)+
     1                          RotationMatrix(1,2)/BETA)
         R_SPEED(1,4) = X*(nu_p*RotationMatrix(1,1)-
     1                          RotationMatrix(1,2)/BETA)
C
         R_SPEED(2,1) = RotationMatrix(2,1)
         R_SPEED(2,2) = RotationMatrix(2,1)
         R_SPEED(2,3) = X*(nu_p*RotationMatrix(2,1)+
     1                          RotationMatrix(2,2)/BETA)
         R_SPEED(2,4) = X*(nu_p*RotationMatrix(2,1)-
     1                          RotationMatrix(2,2)/BETA)
C
      IF(ICHECK.NE.0)THEN
C
C       The Characteristic variables ..
C
C       1. Entropy ..
C
         GRAD_CHAR(1,1)=-GRAD_PRIM(2,1)/ASQR + GRAD_PRIM(1,1)
         GRAD_CHAR(1,2)=-GRAD_PRIM(2,2)/ASQR + GRAD_PRIM(1,2)
C
C       2. Total Enthalpy ..
C
         GRAD_CHAR(2,1)= GRAD_PRIM(2,1)/DensAbar + Mach*U_dot_DU(1)
         GRAD_CHAR(2,2)= GRAD_PRIM(2,2)/DensAbar + Mach*U_dot_DU(2)
C
C       3. Acoustic variables ..
C
         GRAD_CHAR(3,1)= TEMP1*GRAD_PRIM(2,1)+Mach*S_dot_DU(1)
         GRAD_CHAR(3,2)= TEMP1*GRAD_PRIM(2,2)+Mach*S_dot_DU(2)
         GRAD_CHAR(4,1)= TEMP1*GRAD_PRIM(2,1)-Mach*S_dot_DU(1)
         GRAD_CHAR(4,2)= TEMP1*GRAD_PRIM(2,2)-Mach*S_dot_DU(2)
C
C       Jacobian matrix of the acoustic+shear subsystem ..
C       
C
         Matrix(1,1,1) = X*(nu_p*RotationMatrix(1,1)+
     &                           RotationMatrix(1,2)/BETA)
         Matrix(1,2,1) = X*nu_m* RotationMatrix(1,1)
C
         Matrix(2,1,1) = Matrix(1,2,1)
         Matrix(2,2,1) = X*(nu_p*RotationMatrix(1,1)-
     &                           RotationMatrix(1,2)/BETA)
C
         Matrix(1,1,1) = X*(nu_p*RotationMatrix(1,1)+
     &                           RotationMatrix(1,2)/BETA)
         Matrix(1,2,1) = X*nu_m* RotationMatrix(1,1)
C
         Matrix(2,1,1) = Matrix(1,2,1)
         Matrix(2,2,1) = X*(nu_p*RotationMatrix(1,1)-
     &                           RotationMatrix(1,2)/BETA)
C
C
C
         Matrix(1,1,2) = X*(nu_p*RotationMatrix(2,1)+
     &                           RotationMatrix(2,2)/BETA)
         Matrix(1,2,2) = X*nu_m* RotationMatrix(2,1)
C
         Matrix(2,1,2) = Matrix(1,2,2)
         Matrix(2,2,2) = X*(nu_p*RotationMatrix(2,1)-
     &                           RotationMatrix(2,2)/BETA)
C
         Matrix(1,1,2) = X*(nu_p*RotationMatrix(2,1)+
     &                           RotationMatrix(2,2)/BETA)
         Matrix(1,2,2) = X*nu_m* RotationMatrix(2,1)
C
         Matrix(2,1,2) = Matrix(1,2,2)
         Matrix(2,2,2) = X*(nu_p*RotationMatrix(2,1)-
     &                           RotationMatrix(2,2)/BETA)
C
      ENDIF
C
C     Right eigenvector matrix in conserved variables ..
C
      CNST = QBAR
C
      TEMP1 = HALF + ONE / (GM1*MACHSQR)
      TEMP2 = CNST * HALF * DensAbar / KINETIC
      TEMP3 = HALF * TEMP2 * BETA / X
      TEMP4 = CNST * DensAbar / X * QINV
C
      dUdV(1,1) = CNST
      dUdV(2,1) = CNST * KINETIC
C
      dUdV(1,2) = TEMP2
      dUdV(2,2) = CNST * DensAbar * ( TEMP1 + ONE )
C
      dUdV(1,3) = TEMP3
      dUdV(2,3) = TEMP3 * UAVG(2)
C
      dUdV(1,4) = dUdV(1,3)
      dUdV(2,4) = dUdV(2,3)
C
      dUdV(1,5) = ZERO
      dUdV(2,5) = ZERO
C
      DO 31 IDIM = 1 , NDIM
         JDIM = IDIM + 2
C
         dUdV(JDIM,1) = CNST         * UAVG(JDIM)
         dUdV(JDIM,2) = TWO  * TEMP2 * UAVG(JDIM) 
         dUdV(JDIM,3) = HALF * TEMP4 * (BETA*RotationMatrix(IDIM,1)
     &                                      +RotationMatrix(IDIM,2))
         dUdV(JDIM,4) = HALF * TEMP4 * (BETA*RotationMatrix(IDIM,1)
     &                                      -RotationMatrix(IDIM,2))
         dUdV(JDIM,5) =        TEMP4 *       RotationMatrix(IDIM,3)
C
   31 CONTINUE
C
#ifdef DEBUG
C
      CALL DGEMM('N','N',4,2,4,ONE,dVdZ,NOFVAR,GRAD_PARM,NMAX,ZERO,
     +           WORK,4) 
C
      do jdim = 1,4
      do idim = 1,4
         work(idim,jdim) = work(idim,jdim) -GRAD_CHAR(idim,jdim)
      enddo
      enddo
      if(.not. nullmat(work,4,2,4,1.D-12) )then
      CALL R8Mat_Print('General',' ',4,2,GRAD_CHAR,LDW,
     +            'DW ',IFAIL)
      CALL R8Mat_Print('General',' ',4,2,WORK,4,
     +            'DW - dW/dZ * grad(Z) ',IFAIL)
      write(6,*)'M^2, beta, abar ',machsqr,beta,abar
      CALL R8Mat_Print('General',' ',NDIM+2,NDIM,GRAD_PARM,NMAX,
     +            'grad(parametr) ',IFAIL)
      CALL R8Mat_Print('General',' ',NDIM+2,NDIM,GRAD_PRIM,NMAX,
     +            'grad(primitive) ',IFAIL)
      CALL R8Mat_Print('General',' ',NDIM+2,NDIM+2,dVdZ,NOFVAR,
     +            'grad(primitive) ',IFAIL)

      pause
      endif
#endif
C
      RETURN
      END
