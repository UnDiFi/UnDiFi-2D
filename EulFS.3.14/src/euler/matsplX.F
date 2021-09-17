!> \copydoc MatSplitIX()
      SUBROUTINE MatSplitX(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplX.F,v 1.10 2020/03/28 09:51:15 abonfi Exp $
C
C     generalized inflow parameters for the Choi-Merkle preconditioner in
c     primitive variables 
C
C     when the pre-processor flag PETSC_USE_DEBUG is defined
C     compare the analyical expresion for \K^{\pm} with
C     the one obtained as R \Lambda^{\pm} L.
C
CCCCCCCC#define FSPL_USE_DEBUG
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'pfcgas.com'
      INCLUDE 'flags.com'
      INCLUDE 'merkle.com'
      INCLUDE 'stream.com'
      INCLUDE 'three.com'
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
      DOUBLE PRECISION l,m,n,U_dot_M,NX,NY,NZ,DELTA,
     &AREAINV,AREA,DENSINV,HELP,T3,T4,D43,D34,TEMP
C
C
#ifdef FSPL_USE_DEBUG
      DOUBLE PRECISION VECM(3),VECS(3),VECT(3)
      DOUBLE PRECISION WKSPA(5,5),WKSPB(5,5)
      DOUBLE PRECISION WKSPC(5,5),WKSPD(5,5)
      DOUBLE PRECISION WORK1(5,5),WORK2(5,5)
      DOUBLE PRECISION KPLUS(5,5),KMINS(5,5)
      INTEGER TABLE(3,3),INDX,IV1(5)
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
C
C     REM: 
C     U_dot_M is the dot product between the velocity vector
C          and the face normal divided by the spatial dimension (d) ..
C
      l = VECN(1)/NDIM
      m = VECN(2)/NDIM
      U_dot_M = UAVG(3) * l + UAVG(4) * m
      IF( NDIM .EQ. 3 )THEN
           n = VECN(3)/NDIM
           U_dot_M = U_dot_M + UAVG(5) * n
      ENDIF
!     HELP = GM1*M_INFTY*M_INFTY
      HELP = GM1OG/RSTAR
C
C     Jacobian Matrix for primitive (p,T,u,v,w) variables ..
C
      KMAT(1,1) = U_dot_M * AMPSQR
      KMAT(1,2) = ZERO
      KMAT(1,3) = UAVG(1)*ASQR*AMPSQR*l
      KMAT(1,4) = UAVG(1)*ASQR*AMPSQR*m
C
      KMAT(2,1) = HELP*(AMPSQR-ONE)*U_dot_M/UAVG(1)
      KMAT(2,2) = U_dot_M
      KMAT(2,3) = HELP*AMPSQR*ASQR*l
      KMAT(2,4) = HELP*AMPSQR*ASQR*m
C
      KMAT(3,1) = l/UAVG(1)
      KMAT(3,2) = ZERO
      KMAT(3,3) = U_dot_M
      KMAT(3,4) = ZERO
C
      KMAT(4,1) = m/UAVG(1)
      KMAT(4,2) = ZERO
      KMAT(4,3) = ZERO
      KMAT(4,4) = U_dot_M
C
      IF( NDIM .NE. 3 )GOTO 7
C
      KMAT(1,5) = UAVG(1)*ASQR*AMPSQR*n
      KMAT(2,5) = HELP*AMPSQR*ASQR*n
      KMAT(3,5) = ZERO
      KMAT(4,5) = ZERO
C
      KMAT(5,1) = n/UAVG(1)
      KMAT(5,2) = ZERO
      KMAT(5,3) = ZERO
      KMAT(5,4) = ZERO
      KMAT(5,5) = U_dot_M
C
C     return if there is no need to split the
C     KMAT matrix into its positive and negative parts
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
      U_dot_M = U_dot_M / AREA
C
C     Now U_dot_M is velocity \cdot unit normal vector
C
      DELTA = 4.d0*ASQR*AMPSQR+((AMPSQR-ONE)*U_dot_M)**2
      DELTA = SQRT(DELTA)
C
C     REMARK that at this stage the Eigenvalues do NOT
C     contain the area
C
      WR(1) = U_dot_M
      WR(2) = U_dot_M
      WR(3) = ( U_dot_M*(ONE+AMPSQR) + DELTA )* HALF
      WR(4) = ( U_dot_M*(ONE+AMPSQR) - DELTA )* HALF
      WR(5) = U_dot_M
C
#if 0
      call rg(LDA,4,KMAT,WKSPA,WKSPB,0,VRIGHT,IV1,WKSPD,IFAIL)
      if (IFAIL.EQ.0)then
      write(6,*)(WKSPA(j,1),j=1,4)
      else
         stop 'rg has failed'
      endif
#endif
C
C    Positive and negative eigenvalues ..
C
      DO 12 IDIM = 1 , NORDER
         LPOS(IDIM) = HALF * (WR(IDIM) + ABS(WR(IDIM)) ) * AREA
         LNEG(IDIM) = HALF * (WR(IDIM) - ABS(WR(IDIM)) ) * AREA
   12 CONTINUE
C
      DENSINV = ONE/UAVG(1)
C
      D43 = WR(4)-WR(3)
      T3 = WR(3)-AMPSQR*WR(1)!U_dot_M
      T4 = WR(4)-AMPSQR*WR(1)
!     HELP = GM1*M_INFTY*M_INFTY/UAVG(1)
      HELP = GM1OG/(RSTAR*UAVG(1))
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
C     Left Eigenvectors (1st col) ..
C
      VLEFT(1,1) =-HELP
      VLEFT(2,1) = ZERO
      VLEFT(3,1) = ONE
      VLEFT(4,1) = ONE
      IF(NDIM.EQ.3) VLEFT(5,1) = ZERO
C
C     Left Eigenvectors (4th col) ..
C
      VLEFT(1,2) = ONE
      VLEFT(2,2) = ZERO
      VLEFT(3,2) = ZERO
      VLEFT(4,2) = ZERO
      IF(NDIM.EQ.3) VLEFT(5,2) = ZERO
C
C     Right Eigenvectors (1st row) ..
C
      VRIGHT(1,1) = ZERO
      VRIGHT(1,2) = ZERO
      VRIGHT(1,3) = T4/D43
      VRIGHT(1,4) =-T3/D43
      IF(NDIM.EQ.3) VRIGHT(1,5) = ZERO
C
C     Right Eigenvectors (2nd row) ..
C
      VRIGHT(2,1) = ONE
      VRIGHT(2,2) = ZERO
      VRIGHT(2,3) = HELP*T4/D43
      VRIGHT(2,4) =-HELP*T3/D43
      IF(NDIM.EQ.3) VRIGHT(2,5) = ZERO
C
      DO 10 JDIM = 1 , NDIM
C
         IDIM = JDIM + 2
         Nxyz = VECM(JDIM)
         Sxyz = VECS(JDIM)
         Txyz = VECT(JDIM)
C
C     Left Eigenvectors ..
C
         VLEFT(1,IDIM) = ZERO
         VLEFT(2,IDIM) = Sxyz
         VLEFT(3,IDIM) = Nxyz*UAVG(1)*T3
         VLEFT(4,IDIM) = Nxyz*UAVG(1)*T4
         IF(NDIM.EQ.3) VLEFT(5,IDIM) = Txyz
C
C     Right Eigenvectors ..
C
         VRIGHT(IDIM,1) = ZERO
         VRIGHT(IDIM,2) = Sxyz
         VRIGHT(IDIM,3) =-DENSINV*Nxyz/D43
         VRIGHT(IDIM,4) = DENSINV*Nxyz/D43 
         VRIGHT(IDIM,5) = Txyz ! not sure
C
   10 CONTINUE
C
C
C    Test if the product of the two eigenvectors equals the
C    identity matrix: VRIGHT*VLEFT = Identity
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LDA,
     +                 VLEFT,LDA,0.D0,WORK1,5)
            IF( UNITMAT(WORK1,NORDER,NORDER,5,1.D-15) .EQV. .FALSE. )
     +      THEN
                CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,5,
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
CREM  VECN is not divided by d ?!?! 
C
C
      NX = VECN(1) * AREAINV
      NY = VECN(2) * AREAINV
      IF( NDIM .EQ. 3 )NZ = VECN(3) * AREAINV
C
      D34 = -ONE/D43
C
      KPOS(1,1) = (LPOS(4)*T3-LPOS(3)*T4)*D34
      KPOS(1,2) = ZERO
      KPOS(1,3) = UAVG(1)*T3*T4*(LPOS(4)-LPOS(3))*NX*D34
      KPOS(1,4) = UAVG(1)*T3*T4*(LPOS(4)-LPOS(3))*NY*D34
C
      KPOS(2,1) = HELP*(KPOS(1,1)-LPOS(1))
      KPOS(2,2) = LPOS(1)
      KPOS(2,3) = HELP*KPOS(1,3)
      KPOS(2,4) = HELP*KPOS(1,4)
C
      KPOS(3,1) = (LPOS(3)-LPOS(4))*D34*DENSINV*NX
      KPOS(3,2) = ZERO
      TEMP = (LPOS(3)*T3-LPOS(4)*T4)*D34 - LPOS(1)
      KPOS(3,3) = LPOS(1) + TEMP*NX*NX
      KPOS(3,4) = TEMP*NX*NY
C
      KPOS(4,1) = (LPOS(3)-LPOS(4))*D34*DENSINV*NY
      KPOS(4,2) = ZERO
      KPOS(4,3) = KPOS(3,4)
      KPOS(4,4) = LPOS(1) + TEMP*NY*NY
C
      IF(NDIM.EQ.3)THEN
C
         KPOS(1,5) = UAVG(1)*T3*T4*(LPOS(4)-LPOS(3))*NZ*D34
         KPOS(2,5) = HELP*KPOS(1,5)
         KPOS(5,1) = (LPOS(3)-LPOS(4))*D34*DENSINV*NZ
         KPOS(5,2) = ZERO
         KPOS(3,5) = TEMP*NX*NZ
         KPOS(4,5) = TEMP*NY*NZ
         KPOS(5,5) = LPOS(1) + TEMP*NZ*NZ
         KPOS(5,3) = KPOS(3,5)
         KPOS(5,4) = KPOS(4,5)
C
      ENDIF
C
C
      KNEG(1,1) = (LNEG(4)*T3-LNEG(3)*T4)*D34
      KNEG(1,2) = ZERO
      KNEG(1,3) = UAVG(1)*T3*T4*(LNEG(4)-LNEG(3))*NX*D34
      KNEG(1,4) = UAVG(1)*T3*T4*(LNEG(4)-LNEG(3))*NY*D34
C
      KNEG(2,1) = HELP*(KNEG(1,1)-LNEG(1))
      KNEG(2,2) = LNEG(1)
      KNEG(2,3) = HELP*KNEG(1,3)
      KNEG(2,4) = HELP*KNEG(1,4)
C
      KNEG(3,1) = (LNEG(3)-LNEG(4))*D34*DENSINV*NX
      KNEG(3,2) = ZERO
      TEMP = (LNEG(3)*T3-LNEG(4)*T4)*D34 - LNEG(1)
      KNEG(3,3) = LNEG(1) + TEMP*NX*NX
      KNEG(3,4) = TEMP*NX*NY
C
      KNEG(4,1) = (LNEG(3)-LNEG(4))*D34*DENSINV*NY
      KNEG(4,2) = ZERO
      KNEG(4,3) = KNEG(3,4)
      KNEG(4,4) = LNEG(1) + TEMP*NY*NY
C
      IF(NDIM.EQ.3)THEN
C
         KNEG(1,5) = UAVG(1)*T3*T4*(LNEG(4)-LNEG(3))*NZ*D34
         KNEG(2,5) = HELP*KNEG(1,5)
         KNEG(5,1) = (LNEG(3)-LNEG(4))*D34*DENSINV*NZ
         KNEG(5,2) = ZERO
         KNEG(3,5) = TEMP*NX*NZ
         KNEG(4,5) = TEMP*NY*NZ
         KNEG(5,5) = LNEG(1) + TEMP*NZ*NZ
         KNEG(5,3) = KNEG(3,5)
         KNEG(5,4) = KNEG(4,5)
C
      ENDIF
C
C
C
#ifdef FSPL_USE_DEBUG
C
      DO 32 I =  1, NORDER
         DO 32 J =  1, NORDER
            WORK1(I,J) = KPOS(I,J)-KPLUS(I,J)
            WORK2(I,J) = KNEG(I,J)-KMINS(I,J)
   32 CONTINUE
      IF( NULLMAT(WORK1,NORDER,NORDER,5,1.D-15) .EQV. .FALSE. .OR.
     +    NULLMAT(WORK2,NORDER,NORDER,5,1.D-15) .EQV. .FALSE. )THEN
         WRITE(6,*)'Errors on Inflow parameters for IELEM =',ielem
         WRITE(6,*)'Positive = ',NULLMAT(WORK1,NORDER,NORDER,5,1.D-15)
         WRITE(6,*)'Negative = ',NULLMAT(WORK2,NORDER,NORDER,5,1.D-15)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS,LDA,
     +               'true positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPLUS,5,
     +               'analytical positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,5,
     +               'error positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG,LDA,
     +               'true negative matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KMINS,5,
     +               'analytical negative matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK2,5,
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
CREM           WORK1(I,J) = (KPLUS(I,J)+KMINS(I,J))
               WORK2(I,J) = WORK1(I,J)-KMAT(I,J)
   30 CONTINUE
      LFLAG = ( .NOT. NULLMAT(WORK2,NORDER,NORDER,5,1.D-15) )
C
      IF( LFLAG .EQV. .TRUE. )THEN
C
            LFLAG = .FALSE.
C
            WRITE(6,*)'Element # ',IELEM
            WRITE(6,*)'Eigenvalues ',(WR(I),I=1,NORDER)
C
C     skip the following ...
C
#if 0
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA,
     +                  'LEFT eigenvectors',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA,
     +                  'RIGHT eigenvectors',IFAIL)
C
C     compute L*K and K*R
C     REM: we divide by Area/d
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE/AREA,VLEFT,LDA,
     +                 KMAT,LDA,ZERO,WKSPA,5)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE/AREA,KMAT,LDA,
     +                 VRIGHT,LDA,ZERO,WKSPB,5)
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
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,5,'LEFT*A ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,5,'lambda*LEFT',
     +                  IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,5, 'A*RIGHT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,5, 'RIGHT*lambda',
     +                  IFAIL)
#endif
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Analytical inflow parameter matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,5,
     +                  'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,5,
     +                  'Error matrix',IFAIL)
C           PAUSE
      ENDIF
#endif
C
      RETURN
      END
