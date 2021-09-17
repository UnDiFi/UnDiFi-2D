      SUBROUTINE parm2charHE(ZAVG,dVdZ,dUdV,NDOF,NDIM)
C
C     $Id: parm2charHE.f,v 1.2 2013/01/29 14:33:34 abonfi Exp $
C
C     van Leer, Lee, Roe preconditioner in conservative form ..
C
C     compute the transformation matrices from parameter vector
C     to conserved and vice-versa
C
C
C#define FSPL_USE_DEBUG
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C
      INTEGER NDOF,NDIM
C
      DOUBLE PRECISION ZAVG(NDOF),DVDZ(NDOF,*),DUDV(NDOF,*)
C
      DOUBLE PRECISION MSQRM1,CNST,ZDOTN,ZDOTS,ZDOTT
      DOUBLE PRECISION DENSABAR,TEMP1,TEMP2,TEMP3,TEMP4,TEMPA,
     +DENSABARINV,DENSINVABARINV,TEMP0,TEMPX,TEMPY,TEMPZ,
     2COSG,SING,QBAR,QINV,BETA,BETASQR,NU_M,NU_P,X
      INTEGER k
C
      DOUBLE PRECISION UAVG(MAXNOFEQN),ROTMAT(3,3)
      DOUBLE PRECISION KINETIC,ASQR,ABAR,MACH,MACHSQR
C
      DOUBLE PRECISION DDOT,fun_BETA_aldo
      EXTERNAL         DDOT,fun_BETA_aldo
C
C     Mach related variables ..
C
      IF(NDIM.EQ.3)
     &UAVG(5) = ZAVG(5)/ZAVG(1) ! z component of the velocity vector
      UAVG(4) = ZAVG(4)/ZAVG(1) ! y component of the velocity vector
      UAVG(3) = ZAVG(3)/ZAVG(1) ! x component of the velocity vector
      UAVG(2) = ZAVG(2)/ZAVG(1) ! Total Enthalpy
      UAVG(1) = ZAVG(1)*ZAVG(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5) 
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(6,FMT=*)'Negative averaged sound speed ',ASQR,-100
         WRITE(6,FMT=*)'d = ',ndim
         WRITE(6,FMT=*)'k = ',kinetic
         WRITE(6,FMT=*)(zavg(k),k=1,NDOF)
         WRITE(6,FMT=*)(uavg(k),k=1,NDOF)
         STOP
      ENDIF
      QBAR = SQRT(TWO*KINETIC)
      QINV = ONE/QBAR
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )
C
      MSQRM1= MACHSQR - ONE
caldo BETA 	= SQRT( MAX( EPSQR_SONIC , ABS( MSQRM1 ) ) )
      BETA= fun_BETA_aldo( MSQRM1 )
      BETASQR= BETA*BETA
      X= BETA / MAX( MACH , ONE )
      nu_p= HALF * ( MSQRM1/BETASQR + ONE )
      nu_m= HALF * ( MSQRM1/BETASQR - ONE )
C
      DENSABAR = UAVG(1)*ABAR
      TEMP1 = BETA/DENSABAR
C
      DENSABARINV = UAVG(1)/ABAR
      DENSINVABARINV = ONE/DENSABAR
C
      IF(NDIM.EQ.2)THEN
C
         ROTMAT(1,1) = UAVG(3)*QINV
         ROTMAT(2,1) = UAVG(4)*QINV
         ROTMAT(3,1) = ZERO
C
         ROTMAT(1,2) =-UAVG(4)*QINV
         ROTMAT(2,2) = UAVG(3)*QINV
         ROTMAT(3,2) = ZERO
C
         ROTMAT(1,3) = ZERO
         ROTMAT(2,3) = ZERO
         ROTMAT(3,3) = ONE
C
      ELSE
C
         ROTMAT(1,1) = UAVG(3) * QINV
         ROTMAT(2,1) = UAVG(4) * QINV
         ROTMAT(3,1) = UAVG(5) * QINV
C
         TEMP2 = ONE / SQRT( UAVG(3)*UAVG(3)+UAVG(4)*UAVG(4) )
C
         COSG = ONE
         SING = ZERO
C
         ROTMAT(1,2) =(-UAVG(4)*COSG-UAVG(3)*UAVG(5)*SING*
     1                  QINV)*TEMP2
      ROTMAT(2,2) =( UAVG(3)*COSG-UAVG(4)*UAVG(5)*SING*
     1                  QINV)*TEMP2
         ROTMAT(3,2) =  SING*QINV/TEMP2
C
         ROTMAT(1,3) =( UAVG(4)*SING-UAVG(3)*UAVG(5)*COSG*
     1                  QINV)*TEMP2
         ROTMAT(2,3) =(-UAVG(3)*SING-UAVG(4)*UAVG(5)*COSG*
     1                  QINV)*TEMP2
         ROTMAT(3,3) =  COSG*QINV/TEMP2
      ENDIF
C
C     Transformation matrix from parameter vector 
C     to characteristic variables
C
      ZDOTN = ROTMAT(1,1) * ZAVG(3) +
     1        ROTMAT(2,1) * ZAVG(4)  
C
      ZDOTS = ROTMAT(1,2) * ZAVG(3) +
     1        ROTMAT(2,2) * ZAVG(4)  
C
C
      IF(NDIM .EQ. 3)THEN 
        ZDOTN = ZDOTN + ROTMAT(3,1) * ZAVG(5)
        ZDOTS = ZDOTS + ROTMAT(3,2) * ZAVG(5)
        ZDOTT =(ROTMAT(1,3) * ZAVG(3) +
     1          ROTMAT(2,3) * ZAVG(4) +
     1          ROTMAT(3,3) * ZAVG(5) )/UAVG(1)
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
      dVdZ(2,3) =-TEMPX + MACH*ROTMAT(1,1)/ZAVG(1)
      dVdZ(2,4) =-TEMPY + MACH*ROTMAT(2,1)/ZAVG(1)
C
      dVdZ(3,1) = BETA * TEMP2 * ZAVG(2) - MACH * ZDOTS
      dVdZ(3,2) = BETA * dVdZ(2,2)
      dVdZ(3,3) =-BETA * TEMPX + MACH * ROTMAT(1,2) / ZAVG(1)
      dVdZ(3,4) =-BETA * TEMPY + MACH * ROTMAT(2,2) / ZAVG(1)
C
      dVdZ(4,1) = BETA * TEMP2 * ZAVG(2) + MACH * ZDOTS
      dVdZ(4,2) = dVdZ(3,2)
      dVdZ(4,3) =-BETA * TEMPX - MACH * ROTMAT(1,2) / ZAVG(1)
      dVdZ(4,4) =-BETA * TEMPY - MACH * ROTMAT(2,2) / ZAVG(1)
C
C
      IF(NDIM.EQ.3)THEN
C
      dVdZ(1,5) =  TEMP0*ZAVG(5)
      dVdZ(2,5) =-TEMPZ + MACH*ROTMAT(3,1)/ZAVG(1)
      dVdZ(3,5) =-BETA * TEMPZ + MACH * ROTMAT(3,2) / ZAVG(1)
      dVdZ(4,5) =-BETA * TEMPZ - MACH * ROTMAT(3,2) / ZAVG(1)
C
      dVdZ(5,1) = -MACH * ZDOTT
      dVdZ(5,2) = ZERO
      dVdZ(5,3) = MACH * ROTMAT(1,3) / ZAVG(1)
      dVdZ(5,4) = MACH * ROTMAT(2,3) / ZAVG(1)
      dVdZ(5,5) = MACH * ROTMAT(3,3) / ZAVG(1)
C
      ENDIF
C
C     Right eigenvector matrix in conserved variables ..
C
      CNST = ONE/QINV
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
      dUdV(3,3) = HALF * TEMP4 * (BETA*ROTMAT(1,1)
     &                                +ROTMAT(1,2))
      dUdV(3,4) = HALF * TEMP4 * (BETA*ROTMAT(1,1)
     &                                -ROTMAT(1,2))
C
C
      dUdV(4,1) = CNST         * UAVG(4)
      dUdV(4,2) = TWO  * TEMP2 * UAVG(4) 
      dUdV(4,3) = HALF * TEMP4 * (BETA*ROTMAT(2,1)
     &                                +ROTMAT(2,2))
      dUdV(4,4) = HALF * TEMP4 * (BETA*ROTMAT(2,1)
     &                                -ROTMAT(2,2))
C
      IF(NDIM.EQ.3)THEN
C
        dUdV(1,5) = ZERO
        dUdV(2,5) = ZERO
        dUdV(3,5) = TEMP4 *       ROTMAT(1,3)
        dUdV(4,5) = TEMP4 *       ROTMAT(2,3)
C
         dUdV(5,1) = CNST         * UAVG(5)
         dUdV(5,2) = TWO  * TEMP2 * UAVG(5) 
         dUdV(5,3) = HALF * TEMP4 * (BETA*ROTMAT(3,1)
     &                                   +ROTMAT(3,2))
         dUdV(5,4) = HALF * TEMP4 * (BETA*ROTMAT(3,1)
     &                                   -ROTMAT(3,2))
         dUdV(5,5) =        TEMP4 *       ROTMAT(3,3)
C
      ENDIF
C
C
      RETURN
      END
