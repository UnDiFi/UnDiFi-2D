!> \par Purpose
!>
!> Compute positive and negative inflow parameters when the equations are
!> written in Euler symmetrising variables; the split matrices. i.e. \c KPOS and \c KNEG are:
!> \f[
!> \tilde{K}_{\mathbf{n}}^{\pm} = \frac{1}{d}
!> \tilde{R}_{\mathbf{n}} \Lambda_{\mathbf{n}}^{\pm}  \tilde{L}_{\mathbf{n}} = \frac{1}{d}
!>  \left( \begin{array}{ccc}
!>  {\lambda_{0}}^{\pm} & 0 & \mathbf{0} \\\
!>  0 &
!>  \frac{1}{2} \left( \lambda_+^{\pm} + \lambda_-^{\pm} \right) &
!>  \frac{1}{2} \left( \lambda_+^{\pm} - \lambda_-^{\pm} \right) \mathbf{n}^t \\\
!>  0 &
!>  \frac{1}{2} \left( \lambda_+^{\pm} - \lambda_-^{\pm} \right) \mathbf{n} &
!>  \lambda_0^{\pm} I +
!>  \left[ \frac{\lambda_+^{\pm} + \lambda_-^{\pm}}{2} -\lambda_0^{\pm} \right] \mathbf{n} \mathbf{n}
!>  \end{array} \right)
!> \f]
!>
!> Matrix \c KMAT \f$ = \tilde{K}_{\mathbf{n}} = \tilde{K}_{\mathbf{n}}^+ + \tilde{K}_{\mathbf{n}}^-  \f$ is:
!> \f[
!> \tilde{K}_{\mathbf{n}} = \frac{1}{d}
!> \left(   \begin{array}{ccc}
!>    u_n & 0 & \mathbf{0}^t \\\
!>    0 & u_n & a \mathbf{n}^t \\\
!>    \mathbf{0} & a \mathbf{n} & u_n I_{d \times d}
!>    \end{array} \right)
!> \f]
!> The left \c RIGHT and \c LEFT eigenvector matrices are:
!> \f[
!> \tilde{R}_{\mathbf{n}} = \left( \begin{array}{ccccc}
!> 1 &0 & 0 & 0 & 0 \\\
!> 0& 0 & 1 & 1 & 0 \\\
!> 0& \mathbf{s} & \mathbf{n} & -\mathbf{n} & \mathbf{t} \\\
!> \end{array} \right)
!> \quad \mbox{and} \quad
!> \tilde{L}_{\mathbf{n}} = \left( \begin{array}{ccc}
!> 1& 0 & 0 \\\
!> 0& 0 & \mathbf{s}^t \\\
!> 0& \frac{1}{2} & \frac{1}{2} \mathbf{n}^t  \\\
!> 0& \frac{1}{2} &-\frac{1}{2} \mathbf{n}^t  \\\
!> 0& 0 & \mathbf{t}^t
!> \end{array} \right).
!> \f]
!
!> @param[in] IELEM the current simplicial element
!> @param[in] NDIM the dimension of the space
!> @param[in] NORDER the order of the system to be solved
!> @param[in] VECN the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] JacobianMatrix is the jacobian matrix
!> @param[in] LDJ leading dimension of the jacobian matrix
!> @param[out] KMAT  the jacobian matrix times the face normal 
!> @param[out] KPOS is the "positive" part of the jacobian matrix
!> @param[out] KNEG is the "negative" part of the jacobian matrix
!> @param[out] VLEFT is the left eigenvector matrix of the jacobian matrix, i.e. \f$ \tilde{L}_{\mathbf{n}} \f$
!> @param[out] VRIGHT is the right eigenvector matrix of the jacobian matrix, i.e. \f$ \tilde{R}_{\mathbf{n}} \f$
!> @param[in] LDA is the leading dimension of the matrices KMAT, KPOS, KNEG, VLEFT, VRIGHT
!> @param[out] WR are the eigenvalues of the jacobian matrix: \f$ \lambda_0, \lambda_0, \lambda_+, \lambda_-, \lambda_0 \f$
!> @param[out] LPOS are the positive eigenvalues of the jacobian matrix: \f$ \lambda_k^+ \f$
!> @param[out] LNEG are the negative eigenvalues of the jacobian matrix: \f$ \lambda_k^- \f$
!> @param[in] EIGENDECO .TRUE. if the eigenvector decomposition of the matrix is needed, if .FALSE. only the KMAT matrix is computed
!
!> \author $Author: abonfi $
!> \version $Revision: 1.9 $
!> \date $Date: 2020/03/28 09:51:15 $
!> \warning When the pre-processor flag FSPL_USE_DEBUG is defined,
!> the analyical expression for \f$ K^{\pm} \f$ is compared with
!> the one obtained as \f$ R \Lambda^{\pm} L \f$. Should only be used for debugging purposes
!> \warning The eigenvectors VLEFT and VRIGHT are only computed when the pre-processor flag FSPL_USE_DEBUG is defined 
!> \bug Does NOT support moving grids
!
      SUBROUTINE MatSplitVII(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplVII.F,v 1.9 2020/03/28 09:51:15 abonfi Exp $
C
C
C#define FSPL_USE_DEBUG
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'pfcgas.com'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
C
C
      INTEGER IELEM,NDIM,NORDER,LDA,LDJ
      LOGICAL EIGENDECO
C
      DOUBLE PRECISION VECN(NDIM)
      DOUBLE PRECISION JacobianMatrix(LDJ,*)
      DOUBLE PRECISION KMAT(LDA,*),KPOS(LDA,*),KNEG(LDA,*),
     +VLEFT(LDA,*),VRIGHT(LDA,*),WR(*),LPOS(*),LNEG(*)
C
      INTEGER I,J,K,IDIM,JDIM,IFAIL
      DOUBLE PRECISION l,m,n,U_dot_M,AREAINV,AREA
      DOUBLE PRECISION FSPSLP,FSPSLN,FSMSLP,FSMSLN,APLUS,AMINS,NX,NY,NZ
C
#ifdef FSPL_USE_DEBUG
C
      DOUBLE PRECISION VECM(3),VECS(3),VECT(3)
      DOUBLE PRECISION WKSPA(4,4),WKSPB(4,4)
      DOUBLE PRECISION WKSPC(4,4),WKSPD(4,4)
      DOUBLE PRECISION WORK1(4,4),WORK2(4,4)
      DOUBLE PRECISION KPLUS(4,4),KMINS(4,4)
      INTEGER TABLE(3,3),INDX
      DOUBLE PRECISION SURF,SUM1,SUM2,TEMP1,Nxyz,Sxyz,Txyz
      LOGICAL LFLAG
#endif
C
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER  ISDMIN
      LOGICAL  NULLMAT,UNITMAT
      EXTERNAL ISDMIN,DDOT,DNRM2,NULLMAT,UNITMAT
C
      INTRINSIC DSQRT
C
      DATA l,m,n / 3*ZERO /
#ifdef FSPL_USE_DEBUG
      DATA (table(1,j),j=1,3) / 1,0,0 /
      DATA (table(2,j),j=1,3) / 0,1,0 /
      DATA (table(3,j),j=1,3) / 0,0,1 /
      DATA LFLAG / .FALSE. /
#endif
C
C     Executable Statements ..
C
C     U_dot_M is the dot product between the velocity vector
C          and the face normal ..
C
      l = VECN(1)/NDIM
      m = VECN(2)/NDIM
      U_dot_M = UAVG(3) * l + UAVG(4) * m
      IF( NDIM .EQ. 3 )THEN
           n = VECN(3)/NDIM
           U_dot_M = U_dot_M + UAVG(5) * n
      ENDIF
C
C     Jacobian Matrix for symmetrizing variables ..
C
      KMAT(1,1) = U_dot_M
      KMAT(1,2) = ABAR*l
      KMAT(1,3) = ABAR*m
C
      KMAT(2,1) = KMAT(1,2)
      KMAT(2,2) = U_dot_M
      KMAT(2,3) = ZERO
C
      KMAT(3,1) = KMAT(1,3)
      KMAT(3,2) = ZERO
      KMAT(3,3) = U_dot_M
C
      IF( NDIM .NE. 3 )GOTO 7
      KMAT(1,4) = ABAR*n
      KMAT(2,4) = ZERO
      KMAT(3,4) = ZERO
      KMAT(4,1) = KMAT(1,4)
      KMAT(4,2) = ZERO
      KMAT(4,3) = ZERO
      KMAT(4,4) = U_dot_M
C
    7 IF( EIGENDECO .EQV. .FALSE. )RETURN
C
C     Some useful constants ..
C
      AREA = VECN(1)*VECN(1)+VECN(2)*VECN(2)
      IF( NDIM .EQ. 3 ) AREA=AREA + VECN(3)*VECN(3)
      AREA = SQRT(AREA)
      AREAINV = ONE/AREA
      AREA = AREA/NDIM
C
C     Eigenvalues ..
C
      U_dot_M = U_dot_M / AREA
C
      WR(1) = AREA * U_dot_M
      WR(2) = AREA * ( U_dot_M + ABAR )
      WR(3) = AREA * ( U_dot_M - ABAR )
      WR(4) = AREA * U_dot_M
C
C    Positive and negative eigenvalues ..
C
      DO 12 IDIM = 1 , NORDER
         LPOS(IDIM) = HALF * (WR(IDIM) + ABS(WR(IDIM)) )
         LNEG(IDIM) = HALF * (WR(IDIM) - ABS(WR(IDIM)) )
   12 CONTINUE
C
#ifdef FSPL_USE_DEBUG
C
      IF( NDIM .EQ. 2 )THEN
         VECM(1) = VECN(1)*AREAINV
         VECM(2) = VECN(2)*AREAINV
         VECM(3) = ZERO
         VECS(1) = - VECM(2)
         VECS(2) = VECM(1)
         VECS(3) = ZERO
         VECT(1) = ZERO
         VECT(2) = ZERO
         VECT(3) =  ONE
      ELSE
         VECM(1) = VECN(1)*AREAINV
         VECM(2) = VECN(2)*AREAINV
         VECM(3) = VECN(3)*AREAINV
         INDX = ISDMIN(NDIM,VECM,1)
         VECS(1) = table(2,INDX) * VECM(3) - table(3,INDX) * VECM(2)
         VECS(2) = table(3,INDX) * VECM(1) - table(1,INDX) * VECM(3)
         VECS(3) = table(1,INDX) * VECM(2) - table(2,INDX) * VECM(1)
         SURF = ONE/DNRM2(NDIM,VECS,1)
         VECS(1) = VECS(1)*SURF
         VECS(2) = VECS(2)*SURF
         VECS(3) = VECS(3)*SURF
         CALL CROSS_PROD( VECM , VECS , VECT )
      ENDIF
C
C     Left Eigenvectors ..
C
      VLEFT(1,1) = ZERO
      VLEFT(2,1) = HALF
      VLEFT(3,1) = HALF
      VLEFT(4,1) = ZERO
C
C     Right Eigenvectors ..
C
      VRIGHT(1,1) = ZERO
      VRIGHT(1,2) = ONE
      VRIGHT(1,3) = ONE
      VRIGHT(1,4) = ZERO
C
      DO 10 JDIM = 1 , NDIM
C
         IDIM = JDIM + 1
         Nxyz = VECM(JDIM)
         Sxyz = VECS(JDIM)
         Txyz = VECT(JDIM)
C
C     Left Eigenvectors ..
C
         VLEFT(1,IDIM) = Sxyz
         VLEFT(2,IDIM) = HALF * Nxyz
         VLEFT(3,IDIM) =-HALF * Nxyz
         VLEFT(4,IDIM) = Txyz
C
C     Right Eigenvectors ..
C
         VRIGHT(IDIM,1) = Sxyz
         VRIGHT(IDIM,2) = Nxyz
         VRIGHT(IDIM,3) =-Nxyz
         VRIGHT(IDIM,4) = Txyz
C
   10 CONTINUE
C
C    Test if the product of the two eigenvectors equals the
C    identity matrix: VRIGHT*VLEFT = Identity
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,VRIGHT,LDA,
     +                 VLEFT,LDA,0.D0,WORK1,4)
            IF( UNITMAT(WORK1,NORDER,NORDER,4,1.D-15) .EQV. .FALSE. )
     +      THEN
                CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,4,
     +                      'L*R = Identity matrix ????',IFAIL)
            PAUSE
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
            KPLUS(I,J) = SUM1
            KMINS(I,J) = SUM2
   15 CONTINUE
C
#endif
C
      NX = VECN(1) * AREAINV
      NY = VECN(2) * AREAINV
      IF( NDIM .EQ. 3 )NZ = VECN(3) * AREAINV
C
      FSPSLP = LPOS(2)+LPOS(3)
CXXX  FSMSLP = LPOS(2)+LPOS(3) !!! previous bug
      FSMSLP = LPOS(2)-LPOS(3)
      FSPSLN = LNEG(2)+LNEG(3)
      FSMSLN = LNEG(2)-LNEG(3)
      APLUS = HALF*FSPSLP - LPOS(1)
      AMINS = HALF*FSPSLN - LNEG(1)
C
      KPOS(1,1) = HALF * FSPSLP
      KPOS(1,2) = HALF * FSMSLP * NX
      KPOS(1,3) = HALF * FSMSLP * NY
C
      KPOS(2,1) = KPOS(1,2)
      KPOS(2,2) = LPOS(1) + APLUS * NX * NX
      KPOS(2,3) = APLUS * NX * NY
C
      KPOS(3,1) = KPOS(1,3)
      KPOS(3,2) = KPOS(2,3)
      KPOS(3,3) = LPOS(1) + APLUS * NY * NY
C
C
C
      KNEG(1,1) = HALF * FSPSLN
      KNEG(1,2) = HALF * FSMSLN * NX
      KNEG(1,3) = HALF * FSMSLN * NY
C
      KNEG(2,1) = KNEG(1,2)
      KNEG(2,2) = LNEG(1) + AMINS * NX * NX
      KNEG(2,3) = AMINS * NX * NY
C
      KNEG(3,1) = KNEG(1,3)
      KNEG(3,2) = KNEG(2,3)
      KNEG(3,3) = LNEG(1) + AMINS * NY * NY
C
C
      IF(NDIM.EQ.3)THEN
C
          KPOS(1,4) = HALF * FSMSLP * NZ
          KPOS(2,4) = APLUS * NX * NZ
          KPOS(3,4) = APLUS * NY * NZ
C
          KPOS(4,1) = KPOS(1,4)
          KPOS(4,2) = KPOS(2,4)
          KPOS(4,3) = KPOS(3,4)
          KPOS(4,4) = LPOS(1) + APLUS * NZ * NZ
C
          KNEG(1,4) = HALF * FSMSLN * NZ
          KNEG(2,4) = AMINS * NX * NZ
          KNEG(3,4) = AMINS * NY * NZ
C
          KNEG(4,1) = KNEG(1,4)
          KNEG(4,2) = KNEG(2,4)
          KNEG(4,3) = KNEG(3,4)
          KNEG(4,4) = LNEG(1) + AMINS * NZ * NZ
C
      ENDIF
C
#ifdef FSPL_USE_DEBUG
C
      DO 32 I =  1, NORDER
         DO 32 J =  1, NORDER
            WORK1(I,J) = KPOS(I,J)-KPLUS(I,J)
            WORK2(I,J) = KNEG(I,J)-KMINS(I,J)
   32 CONTINUE
      IF( NULLMAT(WORK1,NORDER,NORDER,4,1.D-15) .EQV. .FALSE. .OR.
     +    NULLMAT(WORK2,NORDER,NORDER,4,1.D-15) .EQV. .FALSE. )THEN
         WRITE(6,*)ielem
         WRITE(6,*)NULLMAT(WORK1,NORDER,NORDER,4,1.D-15)
         WRITE(6,*)NULLMAT(WORK2,NORDER,NORDER,4,1.D-15)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS,LDA,
     +               'true positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPLUS,4,
     +               'analytical positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,4,
     +               'error positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG,LDA,
     +               'true negative matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KMINS,4,
     +               'analytical negative matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK2,4,
     +               'error negative matrix',IFAIL)
         PAUSE
      ENDIF
C
C     test if the positive and negative matrices
C     sum up to the jacobian matrix 
C
      DO 30 I =  1, NORDER
            DO 30 J =  1, NORDER
               WORK1(I,J) = (KPOS(I,J)+KNEG(I,J))
               WORK2(I,J) = WORK1(I,J)-KMAT(I,J)
   30 CONTINUE
      LFLAG = ( .NOT. NULLMAT(WORK2,NORDER,NORDER,4,1.D-15) )
C
      IF( LFLAG .EQV. .TRUE. )THEN
C
            LFLAG = .FALSE.
C
            WRITE(6,*)'Element # ',IELEM
            CALL R8Mat_Print('G',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Jacobian matrix',IFAIL)
            WRITE(6,*)(WR(I),I=1,NORDER)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA,
     +                  'LEFT eigenvectors',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA,
     +                  'RIGHT eigenvectors',IFAIL)
C
C     compute L*K and K*R
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,VLEFT,LDA,
     +                 KMAT,LDA,ZERO,WKSPA,4)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,KMAT,LDA,
     +                 VRIGHT,LDA,ZERO,WKSPB,4)
C
C     compute R*\Lambda and L*\Lambda 
C
            DO 8 I =  1, NORDER
               DO 8 J =  1, NORDER
                  WKSPC(I,J) = WR(I)*VLEFT(I,J)
                  WKSPD(I,J) = WR(J)*VRIGHT(I,J)
    8       CONTINUE
C
C     test if L*K = \Lambda*L and K*R = R*\Lambda 
C
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,4,'LEFT*A ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,4,'lambda*LEFT',
     +                  IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,4, 'A*RIGHT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,4, 'RIGHT*lambda',
     +                  IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Original matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,4,
     +                  'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,4,
     +                  'Error matrix',IFAIL)
            PAUSE
      ENDIF
#endif
C
      RETURN
      END
