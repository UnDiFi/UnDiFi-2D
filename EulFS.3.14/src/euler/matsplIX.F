!> \brief \b MatSplitIX
!
!> \par Purpose
!>
!> Compute the split matrices for the
!> INcompressible preconditioned equations
!>
!>
!>
!> @param[in] IELEM the current simplicial element
!> @param[in] NDIM the dimension of the space
!> @param[in] NORDER the order of the system to be solved
!> @param[in] VECN the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] JacobianMatrix is the jacobian matrix
!> @param[in] LDJ leading dimension of the jacobian matrix
!> @param[out] KMAT  the jacobian matrix times the face normal 
!> @param[out] KPOS is the "positive" part of the jacobian matrix
!> @param[out] KNEG is the "negative" part of the jacobian matrix
!> @param[out] VLEFT is the left eigenvector matrix of the jacobian matrix
!> @param[out] VRIGHT is the right eigenvector matrix of the jacobian matrix
!> @param[in] LDA is the leading dimension of the matrices KMAT, KPOS, KNEG, VLEFT, VRIGHT
!> @param[out] WR are the eigenvalues of the jacobian matrix
!> @param[out] LPOS are the positive eigenvalues of the jacobian matrix
!> @param[out] LNEG are the negative eigenvalues of the jacobian matrix
!> @param[in] EIGENDECO .TRUE. if the eigenvector decomposition of the matrix is needed, if .FALSE. only the KMAT matrix is computed
C
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2020/03/28 09:51:15 $
!> \warning Define the preprocessor macro DEBUG when modifying this piece of code
C
      SUBROUTINE MatSplitIX(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C************************************************************
C
C     $Id: matsplIX.F,v 1.5 2020/03/28 09:51:15 abonfi Exp $     
C     
C************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
C
C#define DEBUG
C
      INCLUDE 'constants.h'
      INCLUDE 'frame.com'
C
      INTEGER IELEM,NDIM,NORDER,LDA,LDJ
      LOGICAL   EIGENDECO
C
      DOUBLE PRECISION VECN(NDIM),JacobianMatrix(LDJ,*),
     1                 KMAT(LDA,*),KPOS(LDA,*),KNEG(LDA,*),
     2                 VLEFT(LDA,*),VRIGHT(LDA,*),
     3                 WR(*),LPOS(*),LNEG(*)
C
      INTEGER I,J,K,IFAIL,IFBRANCH
      LOGICAL LFLAG
      DOUBLE PRECISION N,S,T,SUM1,SUM2,TEMP1,AREA,TMP,NM1,NP1
C
#ifdef DEBUG
      DOUBLE PRECISION WORK1(3,3),WORK2(3,3)
      DOUBLE PRECISION WKSPA(3,3),WKSPB(3,3)
      DOUBLE PRECISION WKSPC(3,3),WKSPD(3,3)
#endif
C
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
      INTRINSIC DSQRT
      logical nullmat
      EXTERNAL nullmat
C
      DATA LFLAG / .FALSE. /
C
C       N,S and T are the dot products of the face normal
C       with the streamline and two crossflow directions
C       divided by the space dimension  ..
C
      AREA = VECN(1)*VECN(1) + VECN(2)*VECN(2)
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
         AREA = AREA + VECN(3)*VECN(3)
      ENDIF
      AREA = SQRT(AREA)
      N = N/AREA
      S = S/AREA
      T = T/AREA
C
      TMP = AREA/NDIM 
C
C      Jacobian Matrix ..
C
      KMAT(1,1) = N      * TMP
      KMAT(1,2) = S      * TMP
      KMAT(2,1) = S      * TMP
      KMAT(2,2) =-N      * TMP
      IF(NDIM.EQ.3)THEN
         KMAT(1,3) = ZERO
         KMAT(2,3) = T      * TMP
         KMAT(3,1) = ZERO
         KMAT(3,2) = T      * TMP
         KMAT(3,3) = N      * TMP
      ENDIF
C
      IF( EIGENDECO .EQV. .FALSE. )RETURN
C
C      Eigenvalues ..
C
      WR(1) = TMP
      WR(2) =-TMP
      IF(NDIM.EQ.3)WR(3) = N*TMP
C
      IF( ABS(N-ONE) .LE. 1.D-15 )THEN
         IFBRANCH = 1
         VLEFT(1,1) = ONE
         VLEFT(1,2) = ZERO
         VLEFT(2,1) = ZERO
         VLEFT(2,2) = ONE
         IF(NDIM.EQ.3)THEN 
             VLEFT(1,3) = ZERO
             VLEFT(2,3) = ZERO
             VLEFT(3,1) = ZERO
             VLEFT(3,2) = ZERO
             VLEFT(3,3) = ONE
         ENDIF
C
         VRIGHT(1,1) = ONE
         VRIGHT(1,2) = ZERO
         VRIGHT(2,1) = ZERO
         VRIGHT(2,2) = ONE
         IF(NDIM.EQ.3)THEN 
             VRIGHT(1,3) = ZERO
             VRIGHT(2,3) = ZERO
             VRIGHT(3,1) = ZERO
             VRIGHT(3,2) = ZERO
             VRIGHT(3,3) = ONE
         ENDIF
      ELSEIF( ABS(N+ONE) .LE. 1.D-15 )THEN
         IFBRANCH = 2
         VLEFT(1,1) = ZERO
         VLEFT(1,2) = ONE
         VLEFT(2,1) = ONE
         VLEFT(2,2) = ZERO
         IF(NDIM.EQ.3)THEN 
             VLEFT(1,3) = ZERO
             VLEFT(2,3) = ZERO
             VLEFT(3,1) = ZERO
             VLEFT(3,2) = ZERO
             VLEFT(3,3) = ONE
         ENDIF
C
         VRIGHT(1,1) = ZERO
         VRIGHT(1,2) = ONE
         VRIGHT(2,1) = ONE
         VRIGHT(2,2) = ZERO
         IF(NDIM.EQ.3)THEN 
             VRIGHT(1,3) = ZERO
             VRIGHT(2,3) = ZERO
             VRIGHT(3,1) = ZERO
             VRIGHT(3,2) = ZERO
             VRIGHT(3,3) = ONE
         ENDIF
      ELSE
C
         NM1 = N-ONE
         NP1 = N+ONE
C
         IFBRANCH = 3
         VLEFT(1,1) =-S
         VLEFT(1,2) = N-ONE
         VLEFT(2,1) =-S
         VLEFT(2,2) = N+ONE
         IF(NDIM.EQ.3)THEN 
             VLEFT(1,3) =-T
             VLEFT(2,3) =-T
             VLEFT(3,1) =-T
             VLEFT(3,2) = ZERO
             VLEFT(3,3) = S
         ENDIF
C
         VRIGHT(1,1) = HALF*S / NM1
         VRIGHT(1,2) =-HALF*S / NP1
         VRIGHT(2,1) =-HALF
         VRIGHT(2,2) = HALF
         IF(NDIM.EQ.3)THEN 
             VRIGHT(1,3) = T / (NM1*NP1)
             VRIGHT(2,3) = ZERO
             VRIGHT(3,1) = HALF*T/NM1
             VRIGHT(3,2) =-HALF*T/NP1
             VRIGHT(3,3) =-S / (NM1*NP1)
         ENDIF
      ENDIF
C
C        Positive and negative eigenvalues ..
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
C     WRITE(6,*)'If branch ',IFBRANCH,'in elmt ',ielem
#ifdef DEBUG
         DO 30 I =  1, NORDER
            DO 30 J =  1, NORDER
               WORK1(I,J) = (KPOS(I,J)+KNEG(I,J))
               WORK2(I,J) = WORK1(I,J)-KMAT(I,J)
c              IF( DABS(WORK2(I,J)) .GT. 1.D- 14)LFLAG = .TRUE.
   30    CONTINUE
         LFLAG = NULLMAT(WORK2,NORDER,NORDER,3,1.D-15)
C
         IF( LFLAG .NEQV. .TRUE. )THEN
 
            WRITE(6,*)'If branch ',IFBRANCH,N,S,T
caldo
caldo       CALL R8Mat_Print('G',' ',3,3,RotationMatrix,3, 'Q-matrix',
caldo+      IFAIL)
caldo
C
            LFLAG = .FALSE.
C
            CALL R8Mat_Print('G',' ',NORDER,NORDER,KMAT,LDA, 
     >'Jacobian matrix',IFAIL)
            WRITE(6,*)(WR(I),I=1,NORDER)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA, 
     >'LEFT eigenvectors',
     +      IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA, 
     >'RIGHT eigenvectors',
     +      IFAIL)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VLEFT,LDA, 
     >KMAT,LDA,ZERO,
     +      WKSPA,3)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0, KMAT,LDA,
     >VRIGHT,LDA,ZERO,
     +      WKSPB,3)
            DO 6 I =  1, NORDER
               DO 6 J =  1, NORDER
                  WKSPC(I,J) = WR(I)*VLEFT(I,J)
                  WKSPD(I,J) = WR(J)*VRIGHT(I,J)
    6       CONTINUE
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,3, 
     >'LEFT*A ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,3, 
     >'lambda*LEFT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,3, 
     >'A*RIGHT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,3, 
     >'RIGHT*lambda',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     +      'Original matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,3,
     +      'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,3,
     +      'Error matrix',IFAIL)
C
C        Test if VRIGHT*VLEFT = Identity ..
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LDA, 
     >                 VLEFT,LDA,0.D0,WORK1,3)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,3, 
     >                  'Identity matrix ????',IFAIL)
            PAUSE
         ENDIF
#endif
C
      RETURN
      END
