!> \copydoc MatSplitIX()
      SUBROUTINE MatSplitII(IELEM,NDIM,NORDER,VECN,JacobianMatrix,LDJ,
     +KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplII.F,v 1.11 2020/03/28 09:51:15 abonfi Exp $
C
C#define FSPL_USE_DEBUG
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'blkEulerII.com'
      INCLUDE 'frame.com'
C
C
      INTEGER   IELEM,NDIM,NORDER,LDA,LDJ
      LOGICAL   EIGENDECO
C
      DOUBLE PRECISION  VECN(NDIM),JacobianMatrix(LDJ,*)
      DOUBLE PRECISION  KMAT(LDA,*),  KPOS(LDA,*),KNEG(LDA,*),
     1                  VLEFT(LDA,*),VRIGHT(LDA,*),
     2                  WR(*),LPOS(*),LNEG(*)
C
      INTEGER I,J,K,IFBRANCH
      DOUBLE PRECISION N,S,T,NSIGN,SUM1,SUM2
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3,TEMP4,TEMP5
C
#ifdef FSPL_USE_DEBUG
      INTEGER IFAIL
      DOUBLE PRECISION WORK1(3,3),WORK2(3,3)
      DOUBLE PRECISION WKSPA(3,3),WKSPB(3,3)
      DOUBLE PRECISION WKSPC(3,3),WKSPD(3,3)
      LOGICAL LFLAG
      LOGICAL NULLMAT
#endif
C
C
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
      INTRINSIC DSQRT
C
C
C       N,S and T are the dot products of the face normal
C       with the streamline and two crossflow directions
C       divided by the space dimension  ..
C
      N = RotationMatrix(1,1)*VECN(1) +
     >    RotationMatrix(2,1)*VECN(2)
      S = RotationMatrix(1,2)*VECN(1) +
     >    RotationMatrix(2,2)*VECN(2)
      T = ZERO
      IF(NDIM.EQ.3)THEN 
          N = N + RotationMatrix(3,1)*VECN(3)
          S = S + RotationMatrix(3,2)*VECN(3)
          T = RotationMatrix(1,3)*VECN(1) +
     >        RotationMatrix(2,3)*VECN(2) +
     >        RotationMatrix(3,3)*VECN(3)
      ENDIF 
      N = N/NDIM
      S = S/NDIM
      T = T/NDIM
C
C     Jacobian Matrix ..
C
      TEMP1 = X*nu_p*N
      TEMP2 = X*nu_m*N
      TEMP3 = X*S/BETA
      TEMP4 = X*T/BETA
      KMAT(1,1) = TEMP1+TEMP3
      KMAT(1,2) = TEMP2
      KMAT(2,1) = TEMP2
      KMAT(2,2) = TEMP1-TEMP3
      IF(NDIM.EQ.3)THEN
         KMAT(1,3) = TEMP4
         KMAT(2,3) = TEMP4
         KMAT(3,1) = HALF*TEMP4
         KMAT(3,2) = HALF*TEMP4
         KMAT(3,3) = X*N
      ENDIF
C
      IF( EIGENDECO .EQV. .FALSE. )RETURN
C
C     Eigenvalues ..
C
      TEMP5 = T*T+S*S
      TEMP3 = DSQRT( (nu_m*nu_m*BETASQR*N*N+TEMP5 ) )
      TEMP2 = TEMP3 / BETA
C
C
      WR(1) = X*(nu_p*N - TEMP2)
      WR(2) = X*(nu_p*N + TEMP2)
      WR(3) = X*N
C
      IF( TEMP5 .GT. 1.D-15 )THEN
         IFBRANCH = 1
         TEMP1 = BETA*nu_m*N
         TEMP2 = HALF * TEMP3 / BETA
         VLEFT(1,1) = - S - TEMP1 + TEMP3
         VLEFT(1,2) =   S - TEMP1 + TEMP3
         VLEFT(2,1) = - S - TEMP1 - TEMP3
         VLEFT(2,2) =   S - TEMP1 - TEMP3
         IF(NDIM.EQ.3)THEN
            VLEFT(1,3) = - TWO * T
            VLEFT(2,3) = - TWO * T
            VLEFT(3,1) = - T
            VLEFT(3,2) =   T
            VLEFT(3,3) = TWO * S
         ENDIF
*
         TEMP2 = ONE / (4.D0 * TEMP3 * TEMP5)
         TEMP4 = T*T
         VRIGHT(1,1) = TEMP2 * ( S * VLEFT(2,2) + TEMP4 )
         VRIGHT(1,2) = TEMP2 * (-S * VLEFT(1,2) - TEMP4 )
         VRIGHT(2,1) = TEMP2 * (-S * VLEFT(2,1) + TEMP4 )
         VRIGHT(2,2) = TEMP2 * ( S * VLEFT(1,1) - TEMP4 )
         IF(NDIM.EQ.3)THEN
            VRIGHT(1,3) = TEMP2 * (-TWO * T * TEMP3 )
            VRIGHT(2,3) = TEMP2 * ( TWO * T * TEMP3 )
            VRIGHT(3,1) = TEMP2 * ( T * (VLEFT(2,2)-S) )
            VRIGHT(3,2) = TEMP2 * (-T * (VLEFT(1,1)+S) )
            VRIGHT(3,3) = TEMP2 * ( TWO * S * TEMP3 )
         ENDIF
      ELSE
         IFBRANCH = 2
         NSIGN = SIGN( ONE , N )
         VLEFT(1,1) = ONE
         VLEFT(1,2) = NSIGN
         VLEFT(2,1) = ONE
         VLEFT(2,2) =-NSIGN
         IF(NDIM.EQ.3)THEN
            VLEFT(1,3) = ZERO
            VLEFT(2,3) = ZERO
            VLEFT(3,1) = ZERO
            VLEFT(3,2) = ZERO
            VLEFT(3,3) = ONE
         ENDIF
*
         VRIGHT(1,1) = HALF
         VRIGHT(1,2) = HALF
         VRIGHT(2,1) = HALF*NSIGN
         VRIGHT(2,2) =-HALF*NSIGN
         IF(NDIM.EQ.3)THEN
            VRIGHT(1,3) = ZERO
            VRIGHT(2,3) = ZERO
            VRIGHT(3,1) = ZERO
            VRIGHT(3,2) = ZERO
            VRIGHT(3,3) = ONE
         ENDIF
      ENDIF
C
C     Positive and negative eigenvalues ..
C
      LPOS(1) = HALF * (WR(1) + DABS(WR(1)) )
      LNEG(1) = HALF * (WR(1) - DABS(WR(1)) )
      LPOS(2) = HALF * (WR(2) + DABS(WR(2)) )
      LNEG(2) = HALF * (WR(2) - DABS(WR(2)) )
      IF(NDIM.EQ.3)THEN
         LPOS(3) = HALF * (WR(3) + DABS(WR(3)) )
         LNEG(3) = HALF * (WR(3) - DABS(WR(3)) )
      ENDIF
C
C     compute VRIGHT * LAMBDA(+/-) * VLEFT
C
      DO 15 J = 1 , NORDER
         DO 15 I = 1 , NORDER
         SUM1 = ZERO
         SUM2 = ZERO
            DO 13 K = 1 , NORDER
               TEMP1= VRIGHT(I,K) * VLEFT(K,J)
               SUM1 = SUM1 + TEMP1* LPOS(K)
               SUM2 = SUM2 + TEMP1* LNEG(K)
   13       CONTINUE
         KPOS(I,J) = SUM1
         KNEG(I,J) = SUM2
   15 CONTINUE
C
#ifdef FSPL_USE_DEBUG
         SUM1 = ZERO
         DO 30 I =  1, NORDER
            DO 30 J =  1, NORDER
               WORK1(I,J) = (KPOS(I,J)+KNEG(I,J))
               WORK2(I,J) = WORK1(I,J)
               SUM1 = MAX(ABS(WORK1(I,J)),SUM1)
   30    CONTINUE
C
C        normalize by the maximum value
C
         CALL MATSUB(WORK2,3,KMAT,3,NORDER,NORDER)
         LFLAG = NULLMAT(WORK2,NORDER,NORDER,3,1.D-14*SUM1)
         LFLAG = .FALSE.
C
         IF( LFLAG .EQV. .FALSE. )THEN
 
            WRITE(6,*)'If branch, (N,S,T), TEMP5 '
            WRITE(6,*)IFBRANCH,N,S,T,TEMP5
            WRITE(6,*)'nu+ ',nu_p,' nu- ',nu_m,' ',nu_p+nu_m
C
            CALL R8Mat_Print('G',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Jacobian matrix',IFAIL)
            WRITE(6,*)(WR(I),I=1,NORDER)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA, 
     +                  'LEFT eigenvectors',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA, 
     +                  'RIGHT eigenvectors',IFAIL)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VLEFT,LDA,
     +                 KMAT,LDA,ZERO,WKSPA,3)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0, KMAT,LDA,
     +                 VRIGHT,LDA,ZERO,WKSPB,3)
            DO 6 I =  1, NORDER
               DO 6 J =  1, NORDER
                  WKSPC(I,J) = WR(I)*VLEFT(I,J)
                  WKSPD(I,J) = WR(J)*VRIGHT(I,J)
    6       CONTINUE
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,3, 
     >                  'LEFT*A ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,3, 
     >                  'lambda*LEFT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,3, 
     >                  'A*RIGHT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,3, 
     >                  'RIGHT*lambda',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     +                 'Original matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,3,
     +                 'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,3,
     +                 'Error matrix',IFAIL)
C
C      Test if VRIGHT*VLEFT = Identity ..
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LDA,
     >                 VLEFT,LDA,0.D0,WORK1,3)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,3,
     >                  'Identity matrix ????',IFAIL)
C           PAUSE
         ENDIF
#endif
C
      RETURN
      END
