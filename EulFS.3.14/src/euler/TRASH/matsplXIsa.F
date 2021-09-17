      SUBROUTINE MatSplitXISA(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplXI.F,v 1.1 2003/07/17 13:58:36 abonfi Exp $
C
C     generalized inflow parameters in conserved variables 
C     + SA model
C
C     when the pre-processor flag PETSC_USE_DEBUG is defined
C     compare the analyical expresion for \K^{\pm} with
C     the one obtained as R \Lambda^{\pm} L.
C
#define FSPL_USE_DEBUG
C
C************************************************************
C
C     NORDER: is the order of the system to be solved (INPUT)
C     VECN:   is the face normal (INPUT)
C     KMAT:   the jacobian matrix times the face normal (OUTPUT)
C     KPOS:   is the "positive" part of the jacobian matrix (OUTPUT)
C     KNEG:   is the "negative" part of the jacobian matrix (OUTPUT)
C     VLEFT:  is left eigenvector matrix of the jacobian matrix (OUTPUT)
C     VRIGHT: is right eigenvector matrix of the jacobian matrix (OUTPUT)
C     LD:     is the leading dimension of the previous matrices
C     WR:     eigenvalues of the jacobian matrix (OUTPUT)
C     LPOS:   positive eigenvalues of the jacobian matrix (OUTPUT)
C     LNEG:   negative eigenvalues of the jacobian matrix (OUTPUT)
C
C************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
C
      INCLUDE 'three'
      INCLUDE 'flags.com'
C
      INCLUDE 'constants'
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
      DOUBLE PRECISION l,m,n,U_dot_M,U_DOT_S,U_DOT_T,
     &AREAINV,AREA,TRATIO,DELTA,DENSINV,ABARINV
      DOUBLE PRECISION UPLUS,UMINS,AMACHN,UPOS,APOS,UNEG,ANEG
      DOUBLE PRECISION FSPSLP,FSPSLN,FSMSLP,FSMSLN,APLUS,AMINS,NX,NY,NZ
C
#ifdef FSPL_USE_DEBUG
C
      INTEGER LDC
      PARAMETER(LDC=6)
      DOUBLE PRECISION VECM(3),VECS(3),VECT(3)
      DOUBLE PRECISION WKSPA(LDC,LDC),WKSPB(LDC,LDC)
      DOUBLE PRECISION WKSPC(LDC,LDC),WKSPD(LDC,LDC)
      DOUBLE PRECISION WORK1(LDC,LDC),WORK2(LDC,LDC)
      DOUBLE PRECISION KPLUS(LDC,LDC),KMINS(LDC,LDC)
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
C     Jacobian Matrix for conserved variables ..
C
      KMAT(1,1) = ZERO
      KMAT(1,2) = ZERO
      KMAT(1,3) = l
      KMAT(1,4) = m
C
      KMAT(2,1) = (-UAVG(2)+GM1*KINETIC)*U_dot_M
      KMAT(2,2) = GAM*U_dot_M
      KMAT(2,3) = UAVG(2)*l-GM1*U_dot_M*UAVG(3)
      KMAT(2,4) = UAVG(2)*m-GM1*U_dot_M*UAVG(4)
C
      KMAT(3,1) = -U_DOT_M*UAVG(3)+GM1*KINETIC*l
      KMAT(3,2) = GM1 * l
      KMAT(3,3) = UAVG(3)*l - GM1*l*UAVG(3) + U_dot_M
      KMAT(3,4) = UAVG(3)*m - GM1*l*UAVG(4)
C
      KMAT(4,1) = -U_DOT_M*UAVG(4)+GM1*KINETIC*m
      KMAT(4,2) = GM1 * m
      KMAT(4,3) = UAVG(4)*l - GM1*m*UAVG(3)
      KMAT(4,4) = UAVG(4)*m - GM1*m*UAVG(4) + U_dot_M
C
C     turbulent stuff
C
      KMAT(1,NORDER) = ZERO
      KMAT(2,NORDER) = ZERO
      KMAT(3,NORDER) = ZERO
      KMAT(4,NORDER) = ZERO
      KMAT(NORDER,1) = ZERO
      KMAT(NORDER,2) = ZERO
      KMAT(NORDER,3) = ZERO
      KMAT(NORDER,4) = ZERO
      KMAT(NORDER,NORDER) = U_dot_M
C
      IF( NDIM .NE. 3 )GOTO 7
      KMAT(1,5) = n
C
      KMAT(2,5) = UAVG(2)*n-GM1*U_dot_M*UAVG(5)
      KMAT(3,5) = UAVG(3)*n - GM1*l*UAVG(5)
      KMAT(4,5) = UAVG(4)*n - GM1*m*UAVG(5)
C
      KMAT(5,1) = -U_DOT_M*UAVG(5)+GM1*KINETIC*n
      KMAT(5,2) = GM1 * n
      KMAT(5,3) = UAVG(5)*l - GM1*n*UAVG(3) 
      KMAT(5,4) = UAVG(5)*m - GM1*n*UAVG(4)
      KMAT(5,5) = UAVG(5)*n - GM1*n*UAVG(5) + U_dot_M
C
      KMAT(5,6) = ZERO
      KMAT(6,5) = ZERO
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
C     Eigenvalues are divided by d ..
C
      U_dot_M = U_dot_M / AREA
C
      WR(1) = AREA * U_dot_M
      WR(2) = AREA * U_dot_M
      WR(3) = AREA * ( U_dot_M + ABAR )
      WR(4) = AREA * ( U_dot_M - ABAR )
      WR(5) = AREA * U_dot_M
      WR(6) = AREA * U_dot_M
C
C    Positive and negative eigenvalues ..
C
      DO 12 IDIM = 1 , NORDER
         LPOS(IDIM) = HALF * (WR(IDIM) + ABS(WR(IDIM)) )
         LNEG(IDIM) = HALF * (WR(IDIM) - ABS(WR(IDIM)) )
   12 CONTINUE
C
#if 0
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
      U_DOT_S = DDOT(NDIM,UAVG(3),1,VECS,1)
      U_DOT_T = DDOT(NDIM,UAVG(3),1,VECT,1)
      DENSINV = ONE/UAVG(1)
      ABARINV = ONE/ABAR
C
C     Left Eigenvectors (1st row) ..
C
      VLEFT(1,1) = ONE-HALF*GM1*MACHSQR
      VLEFT(2,1) = -U_DOT_S*DENSINV
      VLEFT(3,1) = HALF*DENSINV*(GM1*ABARINV*KINETIC-U_DOT_M)
      VLEFT(4,1) = HALF*DENSINV*(GM1*ABARINV*KINETIC+U_DOT_M)
      IF(NDIM.EQ.3) VLEFT(5,1) = -U_DOT_T*DENSINV
C
C     Left Eigenvectors (2nd row) ..
C
      VLEFT(1,2) = -GM1/ASQR
      VLEFT(2,2) = ZERO
      VLEFT(3,2) = HALF*DENSINV*GM1*ABARINV
      VLEFT(4,2) = HALF*DENSINV*GM1*ABARINV
      IF(NDIM.EQ.3) VLEFT(5,2) = ZERO
C
C     Right Eigenvectors (1st row) ..
C
      VRIGHT(1,1) = ONE
      VRIGHT(1,2) = ZERO
      VRIGHT(1,3) = UAVG(1)*ABARINV
      VRIGHT(1,4) = UAVG(1)*ABARINV
      IF(NDIM.EQ.3) VRIGHT(1,5) = ZERO
C
C     Right Eigenvectors (1st row) ..
C
      VRIGHT(2,1) = KINETIC
      VRIGHT(2,2) = UAVG(1) * U_DOT_S
      VRIGHT(2,3) = UAVG(1) * (UAVG(2)/ABAR + U_DOT_M )
      VRIGHT(2,4) = UAVG(1) * (UAVG(2)/ABAR - U_DOT_M )
      IF(NDIM.EQ.3) VRIGHT(2,5) = UAVG(1) * U_DOT_T
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
         VLEFT(1,IDIM) = GM1/ASQR*UAVG(IDIM)
         VLEFT(2,IDIM) = DENSINV*Sxyz
         VLEFT(3,IDIM) = HALF*DENSINV*(-GM1*ABARINV*UAVG(IDIM)+Nxyz)
         VLEFT(4,IDIM) = HALF*DENSINV*(-GM1*ABARINV*UAVG(IDIM)-Nxyz)
         IF(NDIM.EQ.3) VLEFT(5,IDIM) = DENSINV*Txyz
C
C     Right Eigenvectors ..
C
         VRIGHT(IDIM,1) = UAVG(IDIM)
         VRIGHT(IDIM,2) = UAVG(1)*Sxyz
         VRIGHT(IDIM,3) = UAVG(1)*(UAVG(IDIM)/ABAR+Nxyz)
         VRIGHT(IDIM,4) = UAVG(1)*(UAVG(IDIM)/ABAR-Nxyz)
         VRIGHT(IDIM,5) = UAVG(1)*Txyz
C
   10 CONTINUE
C
C    Test if the product of the two eigenvectors equals the
C    identity matrix: VRIGHT*VLEFT = Identity
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LDA,
     +                 VLEFT,LDA,0.D0,WORK1,5)
            IF( UNITMAT(WORK1,NORDER,NORDER,5,1.D-15) .EQV. .FALSE. )
     +      THEN
                CALL X04CAF('G',' ',NORDER,NORDER,WORK1,5,
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
#endif
C
CREM  VECN is not divided by d ?!?! 
C
C
      NX = VECN(1) * AREAINV
      NY = VECN(2) * AREAINV
      IF( NDIM .EQ. 3 )NZ = VECN(3) * AREAINV
C
      FSPSLP = LPOS(3)+LPOS(4)
      FSMSLP = LPOS(3)-LPOS(4)
      FSPSLN = LNEG(3)+LNEG(4)
      FSMSLN = LNEG(3)-LNEG(4)
      UPLUS = HALF*FSPSLP - LPOS(1)
      UMINS = HALF*FSPSLN - LNEG(1)
      APLUS = HALF*FSMSLP
      AMINS = HALF*FSMSLN
      AMACHN = U_DOT_M/ABAR
      DELTA = GM1 * HALF
      TRATIO = ONE+DELTA*MACHSQR
C
      KPOS(1,1) = LPOS(1) + DELTA*MACHSQR*UPLUS-AMACHN*APLUS
      KPOS(1,2) = GM1/ASQR*UPLUS
      KPOS(1,3) = -GM1*UAVG(3)/ASQR*UPLUS+NX/ABAR*APLUS
      KPOS(1,4) = -GM1*UAVG(4)/ASQR*UPLUS+NY/ABAR*APLUS
C
      KPOS(2,1) = UPLUS*(DELTA*MACHSQR*UAVG(2)-U_DOT_M*U_DOT_M)
     &          + APLUS*AMACHN*(GM1*KINETIC-UAVG(2))
      KPOS(2,2) = HALF*FSPSLP + DELTA*MACHSQR*UPLUS+GM1*AMACHN*APLUS
      KPOS(2,3) = (U_DOT_M*NX-TRATIO*UAVG(3))*UPLUS
     &          + APLUS * (UAVG(2)/ABAR*NX-GM1*AMACHN*UAVG(3))
      KPOS(2,4) = (U_DOT_M*NY-TRATIO*UAVG(4))*UPLUS
     &          + APLUS * (UAVG(2)/ABAR*NY-GM1*AMACHN*UAVG(4))
C
      KPOS(3,1) = (DELTA*MACHSQR*UAVG(3)-U_DOT_M*NX)*UPLUS
     &          + (DELTA*MACHSQR*ABAR*NX-AMACHN*UAVG(3))*APLUS
      KPOS(4,1) = (DELTA*MACHSQR*UAVG(4)-U_DOT_M*NY)*UPLUS
     &          + (DELTA*MACHSQR*ABAR*NY-AMACHN*UAVG(4))*APLUS
C
      KPOS(3,2) = GM1/ABAR*(UAVG(3)/ABAR*UPLUS+NX*APLUS)
      KPOS(4,2) = GM1/ABAR*(UAVG(4)/ABAR*UPLUS+NY*APLUS)
C
      KPOS(3,3) = (NX*NX-GM1*UAVG(3)*UAVG(3)/ASQR)*UPLUS
     &          + (UAVG(3)*NX-GM1*NX*UAVG(3))*APLUS/ABAR + LPOS(1)
      KPOS(3,4) = (NX*NY-GM1*UAVG(3)*UAVG(4)/ASQR)*UPLUS
     &          + (UAVG(3)*NY-GM1*NX*UAVG(4))*APLUS/ABAR
C
      KPOS(4,3) = (NY*NX-GM1*UAVG(4)*UAVG(3)/ASQR)*UPLUS
     &          + (UAVG(4)*NX-GM1*NY*UAVG(3))*APLUS/ABAR
      KPOS(4,4) = (NY*NY-GM1*UAVG(4)*UAVG(4)/ASQR)*UPLUS
     &          + (UAVG(4)*NY-GM1*NY*UAVG(4))*APLUS/ABAR + LPOS(1)
C
C     turbulent stuff
C
      KPOS(1,NORDER) = ZERO
      KPOS(2,NORDER) = ZERO
      KPOS(3,NORDER) = ZERO
      KPOS(4,NORDER) = ZERO
C
      KPOS(NORDER,1) = (DELTA*MACHSQR*UPLUS-AMACHN*APLUS)*UAVG(NORDER)
      KPOS(NORDER,2) = GM1*ABARINV*UAVG(NORDER)*UPLUS
      KPOS(NORDER,3) = UAVG(NORDER)*(-GM1*UAVG(3)*ABARINV*UPLUS
     &                +NX*APLUS/ABAR)
      KPOS(NORDER,4) = UAVG(NORDER)*(-GM1*UAVG(4)*ABARINV*UPLUS
     &                +NY*APLUS/ABAR)
      KPOS(NORDER,NORDER) = LPOS(1)
C
      KNEG(1,1) = LNEG(1) + DELTA*MACHSQR*UMINS-AMACHN*AMINS
      KNEG(1,2) = GM1/ASQR*UMINS
      KNEG(1,3) = -GM1*UAVG(3)/ASQR*UMINS+NX/ABAR*AMINS
      KNEG(1,4) = -GM1*UAVG(4)/ASQR*UMINS+NY/ABAR*AMINS
C
      KNEG(2,1) = UMINS*(DELTA*MACHSQR*UAVG(2)-U_DOT_M*U_DOT_M)
     &          + AMINS*AMACHN*(GM1*KINETIC-UAVG(2))
      KNEG(2,2) = HALF*FSPSLN + DELTA*MACHSQR*UMINS+GM1*AMACHN*AMINS
      KNEG(2,3) = (U_DOT_M*NX-TRATIO*UAVG(3))*UMINS
     &          + AMINS * (UAVG(2)/ABAR*NX-GM1*AMACHN*UAVG(3))
      KNEG(2,4) = (U_DOT_M*NY-TRATIO*UAVG(4))*UMINS
     &          + AMINS * (UAVG(2)/ABAR*NY-GM1*AMACHN*UAVG(4))
C
      KNEG(3,1) = (DELTA*MACHSQR*UAVG(3)-U_DOT_M*NX)*UMINS
     &          + (DELTA*MACHSQR*ABAR*NX-AMACHN*UAVG(3))*AMINS
      KNEG(4,1) = (DELTA*MACHSQR*UAVG(4)-U_DOT_M*NY)*UMINS
     &          + (DELTA*MACHSQR*ABAR*NY-AMACHN*UAVG(4))*AMINS
C
      KNEG(3,2) = GM1/ABAR*(UAVG(3)/ABAR*UMINS+NX*AMINS)
      KNEG(4,2) = GM1/ABAR*(UAVG(4)/ABAR*UMINS+NY*AMINS)
C
      KNEG(3,3) = (NX*NX-GM1*UAVG(3)*UAVG(3)/ASQR)*UMINS
     &          + (UAVG(3)*NX-GM1*NX*UAVG(3))*AMINS/ABAR + LNEG(1)
      KNEG(3,4) = (NX*NY-GM1*UAVG(3)*UAVG(4)/ASQR)*UMINS
     &          + (UAVG(3)*NY-GM1*NX*UAVG(4))*AMINS/ABAR
C
      KNEG(4,3) = (NY*NX-GM1*UAVG(4)*UAVG(3)/ASQR)*UMINS
     &          + (UAVG(4)*NX-GM1*NY*UAVG(3))*AMINS/ABAR
      KNEG(4,4) = (NY*NY-GM1*UAVG(4)*UAVG(4)/ASQR)*UMINS
     &          + (UAVG(4)*NY-GM1*NY*UAVG(4))*AMINS/ABAR + LNEG(1)
C
C     turbulent stuff
C
      KNEG(1,NORDER) = ZERO
      KNEG(2,NORDER) = ZERO
      KNEG(3,NORDER) = ZERO
      KNEG(4,NORDER) = ZERO
C
      KNEG(NORDER,1) = (DELTA*MACHSQR*UMINS-AMACHN*AMINS)*UAVG(NORDER)
      KNEG(NORDER,2) = GM1*ABARINV*UAVG(NORDER)*UMINS
      KNEG(NORDER,3) = UAVG(NORDER)*(-GM1*UAVG(3)*ABARINV*UMINS
     &                +NX*AMINS/ABAR)
      KNEG(NORDER,4) = UAVG(NORDER)*(-GM1*UAVG(4)*ABARINV*UMINS
     &                +NY*AMINS/ABAR)
      KNEG(NORDER,NORDER) = LNEG(1)
C
C
      IF(NDIM.EQ.3)THEN
C
      KPOS(1,5) = -GM1*UAVG(5)/ASQR*UPLUS+NZ/ABAR*APLUS
      KPOS(2,5) = (U_DOT_M*NZ-TRATIO*UAVG(5))*UPLUS
     &          + APLUS * (UAVG(2)/ABAR*NZ-GM1*AMACHN*UAVG(5))
      KPOS(3,5) = (NX*NZ-GM1*UAVG(3)*UAVG(5)/ASQR)*UPLUS
     &          + (UAVG(3)*NZ-GM1*NX*UAVG(5))*APLUS/ABAR
      KPOS(4,5) = (NY*NZ-GM1*UAVG(4)*UAVG(5)/ASQR)*UPLUS
     &          + (UAVG(4)*NZ-GM1*NY*UAVG(5))*APLUS/ABAR
C
      KPOS(5,1) = (DELTA*MACHSQR*UAVG(5)-U_DOT_M*NZ)*UPLUS
     &          + (DELTA*MACHSQR*ABAR*NZ-AMACHN*UAVG(5))*APLUS
      KPOS(5,2) = GM1/ABAR*(UAVG(5)/ABAR*UPLUS+NZ*APLUS)
C
      KPOS(5,3) = (NZ*NX-GM1*UAVG(5)*UAVG(3)/ASQR)*UPLUS
     &          + (UAVG(5)*NX-GM1*NZ*UAVG(3))*APLUS/ABAR
      KPOS(5,4) = (NZ*NY-GM1*UAVG(5)*UAVG(4)/ASQR)*UPLUS
     &          + (UAVG(5)*NY-GM1*NZ*UAVG(4))*APLUS/ABAR
      KPOS(5,5) = (NZ*NZ-GM1*UAVG(5)*UAVG(5)/ASQR)*UPLUS
     &          + (UAVG(5)*NZ-GM1*NZ*UAVG(5))*APLUS/ABAR + LPOS(1)
C
C     turbulent stuff
C
      KPOS(5,NORDER) = ZERO
      KPOS(NORDER,5) = UAVG(NORDER)*(-GM1*UAVG(5)*ABARINV*UPLUS
     &                +NZ*APLUS/ABAR)
C
      KNEG(1,5) = -GM1*UAVG(5)/ASQR*UMINS+NZ/ABAR*AMINS
      KNEG(2,5) = (U_DOT_M*NZ-TRATIO*UAVG(5))*UMINS
     &          + AMINS * (UAVG(2)/ABAR*NZ-GM1*AMACHN*UAVG(5))
      KNEG(3,5) = (NX*NZ-GM1*UAVG(3)*UAVG(5)/ASQR)*UMINS
     &          + (UAVG(3)*NZ-GM1*NX*UAVG(5))*AMINS/ABAR
      KNEG(4,5) = (NY*NZ-GM1*UAVG(4)*UAVG(5)/ASQR)*UMINS
     &          + (UAVG(4)*NZ-GM1*NY*UAVG(5))*AMINS/ABAR
C
      KNEG(5,1) = (DELTA*MACHSQR*UAVG(5)-U_DOT_M*NZ)*UMINS
     &          + (DELTA*MACHSQR*ABAR*NZ-AMACHN*UAVG(5))*AMINS
      KNEG(5,2) = GM1/ABAR*(UAVG(5)/ABAR*UMINS+NZ*AMINS)
C
      KNEG(5,3) = (NZ*NX-GM1*UAVG(5)*UAVG(3)/ASQR)*UMINS
     &          + (UAVG(5)*NX-GM1*NZ*UAVG(3))*AMINS/ABAR
      KNEG(5,4) = (NZ*NY-GM1*UAVG(5)*UAVG(4)/ASQR)*UMINS
     &          + (UAVG(5)*NY-GM1*NZ*UAVG(4))*AMINS/ABAR
      KNEG(5,5) = (NZ*NZ-GM1*UAVG(5)*UAVG(5)/ASQR)*UMINS
     &          + (UAVG(5)*NZ-GM1*NZ*UAVG(5))*AMINS/ABAR + LNEG(1)
C
C     turbulent stuff
C
      KNEG(5,NORDER) = ZERO
      KNEG(NORDER,5) = UAVG(NORDER)*(-GM1*UAVG(5)*ABARINV*UMINS
     &                +NZ*AMINS/ABAR)
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
      IF( NULLMAT(WORK1,NORDER,NORDER,LDC,1.D-15) .EQV. .FALSE. .OR.
     +    NULLMAT(WORK2,NORDER,NORDER,LDC,1.D-15) .EQV. .FALSE. )THEN
         WRITE(6,*)'Errors on Inflow parameters for IELEM =',ielem
         WRITE(6,*)NULLMAT(WORK1,NORDER,NORDER,LDC,1.D-15)
         WRITE(6,*)NULLMAT(WORK2,NORDER,NORDER,LDC,1.D-15)
         CALL X04CAF('G',' ',NORDER,NORDER,KPOS,LDA,
     +               'true positive matrix',IFAIL)
         CALL X04CAF('G',' ',NORDER,NORDER,KPLUS,LDC,
     +               'analytical positive matrix',IFAIL)
         CALL X04CAF('G',' ',NORDER,NORDER,WORK1,LDC,
     +               'error positive matrix',IFAIL)
         CALL X04CAF('G',' ',NORDER,NORDER,KNEG,LDA,
     +               'true negative matrix',IFAIL)
         CALL X04CAF('G',' ',NORDER,NORDER,KMINS,LDC,
     +               'analytical negative matrix',IFAIL)
         CALL X04CAF('G',' ',NORDER,NORDER,WORK2,LDC,
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
      LFLAG = ( .NOT. NULLMAT(WORK2,NORDER,NORDER,LDC,1.D-15) )
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
            CALL X04CAF('G',' ',NORDER,NORDER,VLEFT,LDA,
     +                  'LEFT eigenvectors',IFAIL)
            CALL X04CAF('G',' ',NORDER,NORDER,VRIGHT,LDA,
     +                  'RIGHT eigenvectors',IFAIL)
C
C     compute L*K and K*R
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VLEFT,LDA,
     +                 KMAT,LDA,ZERO,WKSPA,LDC)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,KMAT,LDA,
     +                 VRIGHT,LDA,ZERO,WKSPB,LDC)
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
            CALL X04CAF('G',' ',NORDER,NORDER,WKSPA,LDC,'LEFT*A ',IFAIL)
            CALL X04CAF('G',' ',NORDER,NORDER,WKSPC,LDC,'lambda*LEFT',
     +                  IFAIL)
            CALL X04CAF('G',' ',NORDER,NORDER,WKSPB,LDC, 'A*RIGHT',IFAIL)
            CALL X04CAF('G',' ',NORDER,NORDER,WKSPD,LDC, 'RIGHT*lambda',
     +                  IFAIL)
#endif
            CALL X04CAF('General',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Original matrix',IFAIL)
            CALL X04CAF('General',' ',NORDER,NORDER,WORK1,LDC,
     +                  'Reassembled matrix',IFAIL)
            CALL X04CAF('General',' ',NORDER,NORDER,WORK2,LDC,
     +                  'Error matrix',IFAIL)
C           PAUSE
      ENDIF
#endif
C
      RETURN
      END
