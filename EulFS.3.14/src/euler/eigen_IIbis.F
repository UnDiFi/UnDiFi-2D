      SUBROUTINE Eigen_IIbis(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_IIbis.F,v 1.9 2020/03/28 09:51:14 abonfi Exp $
C
C     van Leer, Lee, Roe preconditioner in conservative form ..
C
C     compute the transformation matrices from parameter vector
C     to conserved and from characteristic to conserved
C
C     if OLD_CODE is defined, the preconditioned residual PF_x
C     is recovered, rather than F_x
C
C#define OLD_CODE
C#define FSPL_USE_DEBUG
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'blkEulerII.com'
      INCLUDE 'frame.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
      INCLUDE 'three.com'
C
      INTEGER LDA,NDIM,LDB
C
      DOUBLE PRECISION Matrix(LDA,LDA,1),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB)
C
      DOUBLE PRECISION MSQRM1,CNST,QBAR,ZDOTN,ZDOTS,ZDOTT
      DOUBLE PRECISION DENSABAR,TEMP1,TEMP2,TEMP3,TEMP4,TEMPA,
     +DENSABARINV,DENSINVABARINV,TEMP0,TEMPX,TEMPY,TEMPZ
      INTEGER IDIM,JDIM
C
      DOUBLE PRECISION U_dot_DU(3),S_dot_DU(3),T_dot_DU(3)
#ifdef FSPL_USE_DEBUG
      DOUBLE PRECISION WORK(15),TOLER
      PARAMETER(TOLER=1.d-15)
      LOGICAL NULLMAT
      INTEGER IFAIL
#endif
C
C
      DOUBLE PRECISION DDOT,fun_BETA_aldo
      EXTERNAL         DDOT,fun_BETA_aldo
C
C
C     Mach related variables ..
C
      MSQRM1= MACHSQR - ONE
caldo BETA 	= SQRT( MAX( EPSQR_SONIC , ABS( MSQRM1 ) ) )
      BETA= fun_BETA_aldo( MSQRM1 )
      BETASQR= BETA*BETA
      X= BETA / MAX( MACH , ONE )
      nu_p= HALF * ( MSQRM1/BETASQR + ONE )
      nu_m= HALF * ( MSQRM1/BETASQR - ONE )
C
      QBAR = ONE / QINV
      DENSABAR = UAVG(1)*ABAR
      TEMP1 = BETA/DENSABAR
C
C
      DENSABARINV = UAVG(1)/ABAR
      DENSINVABARINV = ONE/DENSABAR
C
C     Transformation matrix from parameter vector 
C     to characteristic variables
C
      ZDOTN = RotationMatrix(1,1) * ZAVG(3) +
     1        RotationMatrix(2,1) * ZAVG(4)  
caldo1        RotationMatrix(3,1) * ZAVG(5) )/UAVG(1)
C
      ZDOTS = RotationMatrix(1,2) * ZAVG(3) +
     1        RotationMatrix(2,2) * ZAVG(4)  
caldo1        RotationMatrix(3,2) * ZAVG(5) )/UAVG(1)
C
C
      IF(NDIM .EQ. 3)THEN 
        ZDOTN = ZDOTN + RotationMatrix(3,1) * ZAVG(5)
        ZDOTS = ZDOTS + RotationMatrix(3,2) * ZAVG(5)
        ZDOTT =(RotationMatrix(1,3) * ZAVG(3) +
     1          RotationMatrix(2,3) * ZAVG(4) +
     1          RotationMatrix(3,3) * ZAVG(5) )/UAVG(1)
      ENDIF
      ZDOTN = ZDOTN/UAVG(1)
      ZDOTS = ZDOTS/UAVG(1)
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
C     advection velocity for entropy and total enthalpy
C
      R_SPEED(1,1) = RotationMatrix(1,1)
      R_SPEED(2,1) = RotationMatrix(2,1)
      R_SPEED(1,2) = RotationMatrix(1,1)
      R_SPEED(2,2) = RotationMatrix(2,1)
C
      IF(NDIM.EQ.3)THEN
C
      dVdZ(1,5) =  TEMP0*ZAVG(5)
      dVdZ(2,5) =-TEMPZ + MACH*RotationMatrix(3,1)/ZAVG(1)
      dVdZ(3,5) =-BETA * TEMPZ + MACH * RotationMatrix(3,2) / ZAVG(1)
      dVdZ(4,5) =-BETA * TEMPZ - MACH * RotationMatrix(3,2) / ZAVG(1)
C
      dVdZ(5,1) = -MACH * ZDOTT
      dVdZ(5,2) = ZERO
      dVdZ(5,3) = MACH * RotationMatrix(1,3) / ZAVG(1)
      dVdZ(5,4) = MACH * RotationMatrix(2,3) / ZAVG(1)
      dVdZ(5,5) = MACH * RotationMatrix(3,3) / ZAVG(1)
C
      R_SPEED(3,1) = RotationMatrix(3,1)
      R_SPEED(3,2) = RotationMatrix(3,1)
      ENDIF
C
      IF(ICHECK.NE.0)THEN
C
      CALL GRADPRIM(-10,NDIM,NDIM+2)
C
      DO 5 IDIM = 1 , NDIM
C
         U_dot_dU(IDIM) = DDOT(NDIM,GRAD_PRIM(3,IDIM),1,
     1                         RotationMatrix(1,1),1)
         S_dot_dU(IDIM) = DDOT(NDIM,GRAD_PRIM(3,IDIM),1,
     1                         RotationMatrix(1,2),1)
         T_dot_dU(IDIM) = DDOT(NDIM,GRAD_PRIM(3,IDIM),1,
     1                         RotationMatrix(1,3),1)
C
C	Advection speeds ..
C
         R_SPEED(IDIM,1) = RotationMatrix(IDIM,1)
         R_SPEED(IDIM,2) = RotationMatrix(IDIM,1)
C
         R_SPEED(IDIM,5) = X * RotationMatrix(IDIM,1)
         R_SPEED(IDIM,3) = X*(nu_p*RotationMatrix(IDIM,1)+
     1                             RotationMatrix(IDIM,2)/BETA)
         R_SPEED(IDIM,4) = X*(nu_p*RotationMatrix(IDIM,1)-
     1                             RotationMatrix(IDIM,2)/BETA)
C
C	The Characteristic variables ..
C
C	Entropy ..
C
         GRAD_CHAR(1,IDIM)=-GRAD_PRIM(2,IDIM)/ASQR + GRAD_PRIM(1,IDIM)
C
C	Total Enthalpy ..
C
         GRAD_CHAR(2,IDIM)= GRAD_PRIM(2,IDIM)/DENSABAR 
     1   + Mach*U_dot_DU(IDIM)
C
C	Acoustic variables and shear ..
C
         GRAD_CHAR(5,IDIM)= MACH * T_dot_DU(IDIM)
         GRAD_CHAR(3,IDIM)= TEMP1*GRAD_PRIM(2,IDIM)+Mach*S_dot_DU(IDIM)
         GRAD_CHAR(4,IDIM)= TEMP1*GRAD_PRIM(2,IDIM)-Mach*S_dot_DU(IDIM)
C
C	Jacobian matrix of the acoustic+shear subsystem ..
C
         Matrix(1,1,IDIM) = X*(nu_p*RotationMatrix(IDIM,1)+
     &                              RotationMatrix(IDIM,2)/BETA)
         Matrix(1,2,IDIM) = X*nu_m* RotationMatrix(IDIM,1)
         Matrix(1,3,IDIM) = X*RotationMatrix(IDIM,3)/BETA
C
         Matrix(2,1,IDIM) = Matrix(1,2,IDIM)
         Matrix(2,2,IDIM) = X*(nu_p*RotationMatrix(IDIM,1)-
     &                              RotationMatrix(IDIM,2)/BETA)
         Matrix(2,3,IDIM) = Matrix(1,3,IDIM)
C
         Matrix(3,1,IDIM) = HALF*Matrix(1,3,IDIM)
         Matrix(3,2,IDIM) = Matrix(3,1,IDIM)
         Matrix(3,3,IDIM) = X*RotationMatrix(IDIM,1)
C
caldo    CALL R8Mat_Print('General',' ',NDIM,NDIM,Matrix(1,1,IDIM),LDA,
caldo+            'Jacobian matrix ',IFAIL)
    5 CONTINUE
C
      ENDIF
C
C
#ifdef OLD_CODE
C
C     Non conservative version, the preconditioned residual
C        is recovered (more stable than the conservative version) ..
C
      TEMP1 = HALF * UAVG(1)
      TEMP2 = TEMP1 / (BETA*MACH)
      TEMP3 = UAVG(1) / MACH
      TEMP4 = HALF * TEMP3
      CNST  = HALF * MACHSQR - (GAM - TWO) / GM1
C
C	Right eigenvectors ..
C
C	Entropy ..
C
      dUdV(1,1) = ONE
      dUdV(2,1) = KINETIC
      dUdV(3,1) = UAVG(3)
      dUdV(4,1) = UAVG(4)
C
C	Entropy+Total Enthalpy ..
C
      dUdV(1,2) = ZERO
      dUdV(2,2) = DENSABAR
      dUdV(3,2) = TEMP3*RotationMatrix(1,1)
      dUdV(4,2) = TEMP3*RotationMatrix(2,1)
C
C	First acoustic ..
C
      dUdV(1,3) = TEMP1/(ABAR*BETA)
      dUdV(2,3) = TEMP1*ABAR/BETA*CNST
      dUdV(3,3) = TEMP4*(MSQRM1/BETA*RotationMatrix(1,1)+
     &                               RotationMatrix(1,2))
      dUdV(4,3) = TEMP4*(MSQRM1/BETA*RotationMatrix(2,1)+
     &                               RotationMatrix(2,2))
C
C	Second acoustic ..
C
      dUdV(1,4) = dUdV(1,3)
      dUdV(2,4) = dUdV(2,3)
      dUdV(3,4) = TEMP4*(MSQRM1/BETA*RotationMatrix(1,1)-
     &                            RotationMatrix(1,2))
      dUdV(4,4) = TEMP4*(MSQRM1/BETA*RotationMatrix(2,1)-
     &                            RotationMatrix(2,2))
      IF(NDIM.EQ.3)THEN
C
C	Second shear ..
C
      dUdV(5,1) = UAVG(5)
      dUdV(5,2) = TEMP3*RotationMatrix(3,1)
      dUdV(5,3) = TEMP4*(MSQRM1/BETA*RotationMatrix(3,1)+
     &                               RotationMatrix(3,2))
      dUdV(5,4) = TEMP4*(MSQRM1/BETA*RotationMatrix(3,1)-
     &                            RotationMatrix(3,2))
C
      dUdV(1,5) = ZERO
      dUdV(2,5) = ZERO
      dUdV(3,5) = TEMP3*RotationMatrix(1,3)
      dUdV(4,5) = TEMP3*RotationMatrix(2,3)
      dUdV(5,5) = TEMP3*RotationMatrix(3,3)
      ENDIF
C
#endif
#ifndef OLD_CODE
C
C     Right eigenvector matrix in conserved variables ..
C
      CNST = QBAR
C
      TEMP1 = HALF + ONE / (GM1*MACHSQR) ! H / q^2
      TEMP2 = CNST * HALF * DENSABAR / KINETIC
      TEMP3 = HALF * TEMP2 * BETA / X
      TEMP4 = CNST * DENSABAR / X * QINV
C
      dUdV(1,1) = CNST
      dUdV(2,1) = CNST * KINETIC
C
      dUdV(1,2) = TEMP2
      dUdV(2,2) = CNST * DENSABAR * ( TEMP1 + ONE )
C
      dUdV(1,3) = TEMP3
      dUdV(2,3) = TEMP3 * UAVG(2)
C
      dUdV(1,4) = dUdV(1,3)
      dUdV(2,4) = dUdV(2,3)
C
      dUdV(3,1) = CNST         * UAVG(3)
      dUdV(3,2) = TWO  * TEMP2 * UAVG(3) 
      dUdV(3,3) = HALF * TEMP4 * (BETA*RotationMatrix(1,1)
     &                                +RotationMatrix(1,2))
      dUdV(3,4) = HALF * TEMP4 * (BETA*RotationMatrix(1,1)
     &                                -RotationMatrix(1,2))
C
C
      dUdV(4,1) = CNST         * UAVG(4)
      dUdV(4,2) = TWO  * TEMP2 * UAVG(4) 
      dUdV(4,3) = HALF * TEMP4 * (BETA*RotationMatrix(2,1)
     &                                +RotationMatrix(2,2))
      dUdV(4,4) = HALF * TEMP4 * (BETA*RotationMatrix(2,1)
     &                                -RotationMatrix(2,2))
C
      IF(NDIM.EQ.3)THEN
C
        dUdV(1,5) = ZERO
        dUdV(2,5) = ZERO
        dUdV(3,5) = TEMP4 *       RotationMatrix(1,3)
        dUdV(4,5) = TEMP4 *       RotationMatrix(2,3)
C
         dUdV(5,1) = CNST         * UAVG(5)
         dUdV(5,2) = TWO  * TEMP2 * UAVG(5) 
         dUdV(5,3) = HALF * TEMP4 * (BETA*RotationMatrix(3,1)
     &                                   +RotationMatrix(3,2))
         dUdV(5,4) = HALF * TEMP4 * (BETA*RotationMatrix(3,1)
     &                                   -RotationMatrix(3,2))
         dUdV(5,5) =        TEMP4 *       RotationMatrix(3,3)
C
      ENDIF
C
#endif
C
C
#ifdef FSPL_USE_DEBUG
C
C     check that grad(W) == dVdZ grad(Z)
C
      CALL DGEMM('N','N',NDIM+2,NDIM,NDIM+2,ONE,dVdZ,NDIM+2,
     +           GRAD_PARM,NMAX,ZERO,WORK,NDIM+2)
C
      CALL MATSUB(WORK,NDIM+2,GRAD_CHAR,LDW,NDIM+2,NDIM)
      IF(.NOT. NULLMAT(WORK,NDIM+2,NDIM,NDIM+2,TOLER) )then
      CALL R8Mat_Print('General',' ',NDIM+2,NDIM,GRAD_CHAR,LDW,
     +            'DW ',IFAIL)
      CALL R8Mat_Print('General',' ',NDIM+2,NDIM,WORK,NDIM+2,
     +            'grad(W) - dW/dZ * grad(Z) ',IFAIL)
      PAUSE
      ENDIF
C
#endif
C
      RETURN
      END
