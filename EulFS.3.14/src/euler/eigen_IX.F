      SUBROUTINE Eigen_IX(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     Turkel preconditioner in conservative form ..
C
C#define DEBUG
C#define INVERSE
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'chorin.com'
      INCLUDE 'three.com'
      INCLUDE 'frame.com'
C
      INTEGER LDA,NDIM,LDB
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB)
C
      DOUBLE PRECISION TEMP1,TMP,TMPSQR,VELMAG
      INTEGER IDIM,JDIM
C
#ifdef DEBUG
      DOUBLE PRECISION n_dot_du(3),s_dot_du(3),t_dot_du(3),three
      INTEGER LWORK
      PARAMETER(LWORK=20,THREE=3.d0)
      DOUBLE PRECISION WKSP(4,3),DIFF(4,3),WORK(LWORK),PRE(4,4),
     +RIGHT(4,4),C(4,4),PINV(4,4)
      LOGICAL NULLMAT,UNITMAT
      INTEGER IFAIL,I,J,IPIV(10)
#endif
C
C
      DOUBLE PRECISION DDOT,DNRM2
      EXTERNAL         DDOT,DNRM2
C
C     Mach related variables ..
C
      VELMAG  = DNRM2(NDIM,ZAVG(2),1)
      MACHSQR = VELMAG*VELMAG/BETA
C
C
C     Transformation matrix from primitive
C     to characteristic variables (i.e. the left eigenvectors matrix)
C
C
      dVdZ(1,1) = ONE
      dVdZ(1,2) = ZAVG(2)
      dVdZ(1,3) = ZAVG(3)
C
      dVdZ(2,1) = ZERO
      dVdZ(2,2) = VELMAG * RotationMatrix(1,2)
      dVdZ(2,3) = VELMAG * RotationMatrix(2,2)
C
      dVdZ(3,1) = ONE 
      dVdZ(3,2) = ZERO
      dVdZ(3,3) = ZERO
C
      R_SPEED(1,1) = RotationMatrix(1,1)
      R_SPEED(2,1) = RotationMatrix(2,1)
C
      IF (NDIM.EQ.3)THEN
      dVdZ(1,4) = ZAVG(4)
      dVdZ(2,4) = VELMAG * RotationMatrix(3,2)
      dVdZ(3,4) = ZERO
      dVdZ(4,1) = ZERO
      dVdZ(4,2) = VELMAG * RotationMatrix(1,3)
      dVdZ(4,3) = VELMAG * RotationMatrix(2,3)
      dVdZ(4,4) = VELMAG * RotationMatrix(3,3)
C
      R_SPEED(3,1) = RotationMatrix(3,1)
C
      ENDIF
C
#ifdef DEBUG
C
      DO 5 IDIM = 1 , NDIM
C
C	Advection speed ..
C
         R_SPEED(IDIM,1) = RotationMatrix(IDIM,1)
C
C	Jacobian matrix of the acoustic+shear subsystem ..
C
         Matrix(1,1,IDIM) = RotationMatrix(IDIM,1)
         Matrix(1,2,IDIM) = RotationMatrix(IDIM,2)
         Matrix(1,3,IDIM) = ZERO
C
         Matrix(2,1,IDIM) = Matrix(1,2,IDIM)
         Matrix(2,2,IDIM) =-Matrix(1,1,IDIM)
         Matrix(2,3,IDIM) = RotationMatrix(IDIM,3)
C
         Matrix(3,1,IDIM) = ZERO
         Matrix(3,2,IDIM) = Matrix(2,3,IDIM)
         Matrix(3,3,IDIM) = Matrix(1,1,IDIM)
C
         n_dot_du(idim) = ddot(ndim,grad_parm(2,idim),1,
     >                    rotationmatrix(1,1),1) 
         s_dot_du(idim) = ddot(ndim,grad_parm(2,idim),1,
     >                    rotationmatrix(1,2),1) 
         t_dot_du(idim) = ddot(ndim,grad_parm(2,idim),1,
     >                    rotationmatrix(1,3),1) 
C
         grad_char(1,idim) = grad_parm(1,idim) + velmag*
     >                       n_dot_du(idim)
         grad_char(2,idim) = velmag*s_dot_du(idim)
         grad_char(3,idim) = grad_parm(1,idim)
         grad_char(4,idim) = velmag*t_dot_du(idim)
C
    5 CONTINUE
#endif
C
C     Right eigenvector matrix in primitive variables ..
C
      TMP = ONE/VELMAG
      TMPSQR = TMP*TMP
C
#ifdef INVERSE
C
C     Total pressure
C
      pinv(1,1) = TWO/MACHSQR*velmag
      pinv(2,1) = Three*ZAVG(2)*tmpsqr*velmag
      pinv(3,1) = Three*ZAVG(3)*tmpsqr*velmag
      pinv(4,1) = Three*ZAVG(4)*tmpsqr*velmag
C
      pinv(1,2) = zavg(2)/MACHSQR*velmag
      pinv(2,2) = (one+zavg(2)*zavg(2)*TMPSQR)  * velmag
      pinv(3,2) = (zavg(3)*zavg(2)*TMPSQR    )  * velmag
      pinv(4,2) = (zavg(4)*zavg(2)*TMPSQR    )  * velmag
C
      pinv(1,3) = zavg(3)/MACHSQR*velmag
      pinv(2,3) = (zavg(3)*zavg(2)*TMPSQR    )  * velmag
      pinv(3,3) = (one+zavg(3)*zavg(3)*TMPSQR    )  * velmag
      pinv(4,3) = (zavg(4)*zavg(3)*TMPSQR    )  * velmag
C
      pinv(1,4) = zavg(4)/MACHSQR*velmag
      pinv(2,4) =  (zavg(4)*zavg(2)*TMPSQR    )  * velmag
      pinv(3,4) = (zavg(4)*zavg(3)*TMPSQR    )  * velmag
      pinv(4,4) = (one+zavg(4)*zavg(4)*TMPSQR    )  * velmag
C
#endif
C
C     Right eigenvector matrix in primitive variables ..
C
C
C     Total pressure
C
      dUdV(1,1) = VELMAG/MACHSQR
      dUdV(2,1) = TWO*RotationMatrix(1,1)
      dUdV(3,1) = TWO*RotationMatrix(2,1)
C
      dUdV(1,2) = ZERO
      dUdV(2,2) = RotationMatrix(1,2)
      dUdV(3,2) = RotationMatrix(2,2)
C
      dUdV(1,3) = VELMAG/MACHSQR
      dUdV(2,3) = RotationMatrix(1,1)
      dUdV(3,3) = RotationMatrix(2,1)
C
      IF(NDIM.EQ.3)THEN
      dUdV(4,1) = TWO*RotationMatrix(3,1)
      dUdV(4,2) = RotationMatrix(3,2)
      dUdV(4,3) = RotationMatrix(3,1)
      dUdV(1,4) = ZERO
      dUdV(2,4) = RotationMatrix(1,3)
      dUdV(3,4) = RotationMatrix(2,3)
      dUdV(4,4) = RotationMatrix(3,3)
      ENDIF
C
#ifdef WRONG_STUFF
C
C     Right eigenvector matrix in primitive variables ..
C
      TMP = ONE/VELMAG
      TMPSQR = TMP*TMP
      TEMP1 = -TMPSQR*(3.d0*MACHSQR+TWO)
C
C     Total pressure
C
      dUdV(1,1) =-TMP
      dUdV(2,1) = TWO*TMPSQR*RotationMatrix(1,1)
      dUdV(3,1) = TWO*TMPSQR*RotationMatrix(2,1)
      dUdV(4,1) = TWO*TMPSQR*RotationMatrix(3,1)
C
      dUdV(1,2) = ZERO
      dUdV(2,2) = TMPSQR*RotationMatrix(1,2)
      dUdV(3,2) = TMPSQR*RotationMatrix(2,2)
      dUdV(4,2) = TMPSQR*RotationMatrix(3,2)
C
      dUdV(1,3) = TMP*(2.d0*MACHSQR+1.d0)
      dUdV(2,3) = TEMP1*RotationMatrix(1,1)
      dUdV(3,3) = TEMP1*RotationMatrix(2,1)
      dUdV(4,3) = TEMP1*RotationMatrix(3,1)
C
      dUdV(1,4) = ZERO
      dUdV(2,4) = TMPSQR*RotationMatrix(1,3)
      dUdV(3,4) = TMPSQR*RotationMatrix(2,3)
      dUdV(4,4) = TMPSQR*RotationMatrix(3,3)
C
#endif
C
#ifdef DEBUG
C
      CALL DGEMM('No Transpose','No Transpose',ndim+1,ndim,
     +ndim+1,ONE,dVdZ,LDB,GRAD_PARM,NMAX,ZERO,wksp,4)
      do 3 j = 1,ndim
         do 3 i = 1,ndim+1
            diff(i,j)=grad_char(i,j)-wksp(i,j)
    3 continue
      if( .not. nullmat(diff,ndim+1,ndim,4,1.D-14) )THEN
            CALL R8Mat_Print('General',' ',ndim+1,ndim,grad_char,ldw,
     +      'grad_char ',IFAIL)
            CALL R8Mat_Print('General',' ',ndim+1,ndim,work,4,
     +      'dVdZ * grad(U)',IFAIL)
      pause
      endif
C
      PRE(1,1) = TMP * TWO * MACHSQR
      PRE(1,2) = - TMP * ZAVG(2)
      PRE(1,3) = - TMP * ZAVG(3)
      PRE(1,4) = - TMP * ZAVG(4)
C
      PRE(2,1) = - 3.d0* TMP * ZAVG(2) / BETA
      PRE(3,1) = - 3.d0* TMP * ZAVG(3) / BETA
      PRE(4,1) = - 3.d0* TMP * ZAVG(4) / BETA
C
      PRE(2,2) = TMP * (ONE + ZAVG(2)*ZAVG(2)*TMPSQR)
      PRE(2,3) = TMP * (ZAVG(2)*ZAVG(3)*TMPSQR)
      PRE(2,4) = TMP * (ZAVG(2)*ZAVG(4)*TMPSQR)
C
      PRE(3,2) = TMP * (ZAVG(3)*ZAVG(2)*TMPSQR)
      PRE(3,3) = TMP * (ONE + ZAVG(3)*ZAVG(3)*TMPSQR)
      PRE(3,4) = TMP * (ZAVG(3)*ZAVG(4)*TMPSQR)
C
      PRE(4,2) = TMP * (ZAVG(2)*ZAVG(4)*TMPSQR)
      PRE(4,3) = TMP * (ZAVG(3)*ZAVG(4)*TMPSQR)
      PRE(4,4) = TMP * (ONE + ZAVG(4)*ZAVG(4)*TMPSQR)
C
      CALL DGETRF(NDIM+1,NDIM+1,PRE,4,IPIV,IFAIL)
      CALL DGETRI(NDIM+1,PRE,4,IPIV,WORK,LWORK,IFAIL)
C
      RIGHT(1,1) = ZERO
      RIGHT(2,1) = TMP * RotationMatrix(1,1)
      RIGHT(3,1) = TMP * RotationMatrix(2,1)
      RIGHT(4,1) = TMP * RotationMatrix(3,1)
C
C
      RIGHT(1,2) = ZERO
      RIGHT(2,2) = TMP * RotationMatrix(1,2)
      RIGHT(3,2) = TMP * RotationMatrix(2,2)
      RIGHT(4,2) = TMP * RotationMatrix(3,2)
C
C
      RIGHT(1,3) = ONE
      RIGHT(2,3) =-TMP * RotationMatrix(1,1)
      RIGHT(3,3) =-TMP * RotationMatrix(2,1)
      RIGHT(4,3) =-TMP * RotationMatrix(3,1)
C
C
      RIGHT(1,4) = ZERO
      RIGHT(2,4) = TMP * RotationMatrix(1,3)
      RIGHT(3,4) = TMP * RotationMatrix(2,3)
      RIGHT(4,4) = TMP * RotationMatrix(3,3)
C
      CALL DGEMM('N','N',NDIM+1,NDIM+1,NDIM+1,ONE,PRE,4,
     +           RIGHT,4,ZERO,C,4)
C
C     test if left*right = Identity
C
C     CALL DGEMM('N','N',NDIM+1,NDIM+1,NDIM+1,ONE,
C    +           RIGHT,4,dVdZ,LDB,ZERO,C,4)
C
C     test inv(P) with its analitical expr.
C
      CALL R8Mat_Print('General',' ',ndim+1,ndim+1,pinv,4,
     +      'inverse of P analytical ',IFAIL)
      CALL R8Mat_Print('General',' ',ndim+1,ndim+1,pre,4,
     +      'inverse of P by dgetri ',IFAIL)
C
      CALL R8Mat_Print('General',' ',ndim+1,ndim+1,dUdV,nmax,
     +      'Analitical ',IFAIL)
      CALL R8Mat_Print('General',' ',ndim+1,ndim+1,c,4,
     +      'through DGEMM ',IFAIL)
      pause
C
c     if(.not.unitmat(c,ndim+1,ndim+1,4,1.d-14))then
c     CALL R8Mat_Print('General',' ',ndim+1,ndim+1,c,4,
c    +      'Identity ',IFAIL)
c     endif
c     do j = 1,ndim+1
c     do i = 1,ndim+1
c        dUdv(i,j)=c(i,j)
c     enddo
c     enddo
c     pause
C
#endif
      RETURN
      END
