      SUBROUTINE MatSplitXII(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplXII.F,v 1.10 2020/03/28 09:51:15 abonfi Exp $
C
C     generalized inflow parameters in conserved variables 
C
C     when the pre-processor flag PETSC_USE_DEBUG is defined
C     compare the analyical expresion for \K^{\pm} with
C     the one obtained as R \Lambda^{\pm} L.
C
CC#define FSPL_USE_DEBUG
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
C     NSP:    number of chemical species            
C
C************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'plasma.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'three.com'
      INCLUDE 'four.com'
      INCLUDE 'dofs.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
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
      DOUBLE PRECISION l,m,n,U_dot_M,U_DOT_S,U_DOT_T,
     &AREAINV,AREA,TRATIO,DELTA,ABARINV,CE,CD,CDQ
      DOUBLE PRECISION UPLUS,UMINS,AMACHN,UPOS,APOS,UNEG,ANEG
      DOUBLE PRECISION FSPSLP,FSPSLN,FSMSLP,FSMSLN,APLUS,AMINS,NX,NY,NZ
C
#ifdef FSPL_USE_DEBUG
C
      DOUBLE PRECISION VECM(3),VECS(3),VECT(3)
      DOUBLE PRECISION WKSPA(MAXNOFVAR,MAXNOFVAR),
     2WKSPB(MAXNOFVAR,MAXNOFVAR)
      DOUBLE PRECISION WKSPC(MAXNOFVAR,MAXNOFVAR),
     2WKSPD(MAXNOFVAR,MAXNOFVAR)
      DOUBLE PRECISION WORK1(MAXNOFVAR,MAXNOFVAR),
     2WORK2(MAXNOFVAR,MAXNOFVAR)
      DOUBLE PRECISION WORKL(MAXNOFVAR,MAXNOFVAR),
     2WORKR(MAXNOFVAR,MAXNOFVAR)
      DOUBLE PRECISION KPLUS(MAXNOFVAR,MAXNOFVAR),
     2KMINS(MAXNOFVAR,MAXNOFVAR)
      INTEGER TABLE(3,3),INDX
      DOUBLE PRECISION SURF,SUM1,SUM2,TEMP1,Nxyz,Sxyz,Txyz
      LOGICAL LFLAG
#endif
C
      DOUBLE PRECISION DDOT,DNRM2
C
C     PRESSURE DERIVATIVE (Conservative Variables)          
      DOUBLE PRECISION UU(3) ! crpepe aggIXnto UU
      DOUBLE PRECISION PIR,PIE,PIM
C
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
C     U_dot_M is the dot product between the velocity vector
C          and the face normal ..
C
      l = VECN(1)/NDIM
      m = VECN(2)/NDIM
C
C     UAVG(RHO1,...,RHON,H,U,V,W)
      U_dot_M = UAVG(IX) * l + UAVG(IY) * m ! rpepe 
      IF( NDIM .EQ. 3 )THEN
           n = VECN(3)/NDIM
           U_dot_M = U_dot_M + UAVG(IZ) * n
      ENDIF

C
C     Jacobian Matrix for conserved variables ..
C
C     (Matrix (1...Ns) x (1...Ns)) ..
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            KMAT(I,J) = - ALPHA(I) * U_dot_M
            IF (I.EQ.J)THEN
                KMAT(I,J) = KMAT(I,J) + U_dot_M
            ENDIF
        ENDDO
      ENDDO

C     (Matrix (1...Ns) x (Ns+1...Ns+3)) ..

      DO I = 1 , NSP
        KMAT(I,IE) = ZERO
        KMAT(I,IX) = ALPHA(I) * l
        KMAT(I,IY) = ALPHA(I) * m
      ENDDO

C     (Matrix (Ns+1...Ns+3) x (1...Ns)) ..

      DO J = 1 , NSP
        KMAT(IE,J) =  U_dot_M * (-UAVG(IE) + DR(J))
        KMAT(IX,J) = -U_DOT_M * UAVG(IX) + DR(J) * l
        KMAT(IY,J) = -U_DOT_M * UAVG(IY) + DR(J) * m
      ENDDO

C     (Matrix (Ns+1...Ns+3) x (Ns+1...Ns+3)) ..

      KMAT(IE,IE) = (DE+1) * U_dot_M
      KMAT(IE,IX) = UAVG(IE) * l + U_dot_M * DM(1)
      KMAT(IE,IY) = UAVG(IE) * m + U_dot_M * DM(2)
C
      KMAT(IX,IE) = DE * l
      KMAT(IX,IX) = (UAVG(IX) + DM(1)) * l + U_dot_M
      KMAT(IX,IY) = UAVG(IX) * m + DM(2) * l
C
      KMAT(IY,IE) = DE * m
      KMAT(IY,IX) = UAVG(IY) * l + DM(1) * m
      KMAT(IY,IY) = (UAVG(IY) + DM(2)) * m + U_dot_M
C
      IF( NDIM .NE. 3 )GOTO 7
C      
      DO I = 1 , NSP
        KMAT(I,IZ) = ALPHA(I) * n
      ENDDO
C
      DO J = 1 , NSP
        KMAT(IZ,J) = -U_DOT_M * UAVG(IZ) + DR(J) * n
      ENDDO
C
      KMAT(IE,IZ) = UAVG(IE) * n + U_dot_M * DM(3)
      KMAT(IX,IZ) = UAVG(IX) * n + DM(3) * l
      KMAT(IY,IZ) = UAVG(IY) * n + DM(3) * m
C
      KMAT(IZ,IE) = DE * n
      KMAT(IZ,IX) = UAVG(IZ) * l + DM(1) * n
      KMAT(IZ,IY) = UAVG(IZ)*m + DM(2) * n
      KMAT(IZ,IZ) = (UAVG(IZ) + DM(3)) * n + U_dot_M
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
C     Speed of Sound
C
C
C     Eigenvalues are divided by d ..
C
      U_dot_M = U_dot_M / AREA
C
      DO I = 1 , NSP
        WR(I) = AREA * U_dot_M
      ENDDO

      WR(IE) = AREA * U_dot_M
      WR(IX) = AREA * ( U_dot_M + ABAR )
      WR(IY) = AREA * ( U_dot_M - ABAR )
      WR(IZ) = AREA * U_dot_M
C
C    PositIYe and negative eigenvalues ..
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
      U_DOT_S = DDOT(NDIM,UAVG(IX),1,VECS,1)
      U_DOT_T = DDOT(NDIM,UAVG(IX),1,VECT,1)    

      ABARINV = ONE/ABAR
C
C     Left Eigenvectors (Matrix (1...Ns) x (1...Ns)) ..
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            VLEFT(I,J) = - ALPHA(I) * DR(J) / ASQR
            IF (I.EQ.J)THEN
                VLEFT(I,J) = 1 + VLEFT(I,J)
            ENDIF
        ENDDO
      ENDDO

C     Left Eigenvectors (Matrix (1...Ns) x (Ns+1...Ns+4)) ..
      
      DO I = 1 , NSP
        VLEFT(I,IE) = - ALPHA(I) * DE/ASQR
        VLEFT(I,IX) = ALPHA(I) * DE/ASQR * UAVG(IX)
        VLEFT(I,IY) = ALPHA(I) * DE/ASQR * UAVG(IY)
        IF(NDIM.EQ.3) VLEFT(I,IZ) = ALPHA(I) * DE/ASQR * UAVG(IZ)
      ENDDO

C     Left Eigenvectors (Matrix (Ns+1...Ns+4) x (1...Ns)) ..

      DO J = 1 , NSP
        VLEFT(IE,J) = -U_DOT_S * DENSINV
        VLEFT(IX,J) = HALF * DENSINV * (DR(J)*ABARINV - U_DOT_M)
        VLEFT(IY,J) = HALF * DENSINV * (DR(J)*ABARINV + U_DOT_M)
        IF(NDIM.EQ.3) VLEFT(IZ,J) = -U_DOT_T * DENSINV
      ENDDO

C     Left Eigenvectors (Matrix (Ns+1...Ns+4) x (Ns+1)) ..
      
      VLEFT(IE,IE) = ZERO
      VLEFT(IX,IE) = HALF*DENSINV*DE*ABARINV
      VLEFT(IY,IE) = HALF*DENSINV*DE*ABARINV
      IF(NDIM.EQ.3) VLEFT(IZ,IE) = ZERO
C
C     Right Eigenvectors (Matrix (1...Ns) x (1...Ns)) ..
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            VRIGHT(I,J) = ZERO
            IF (I.EQ.J)THEN
                VRIGHT(I,J) = ONE
            ENDIF
        ENDDO
      ENDDO

C     Right Eigenvectors (Matrix (1...Ns) x (Ns+1 ... Ns+4)) ..
      
      DO I = 1 , NSP
        VRIGHT(I,IE) = ZERO
        VRIGHT(I,IX) = DENS * ALPHA(I) * ABARINV
        VRIGHT(I,IY) = DENS * ALPHA(I) * ABARINV
        IF(NDIM.EQ.3) VRIGHT(I,IZ) = ZERO
      ENDDO

C     Right Eigenvectors (Matrix (Ns+1...Ns+4) x (1...Ns)) ..

      DO J = 1 , NSP
        VRIGHT(IE,J) = 2 * KINETIC - DR(J) / DE
        VRIGHT(IX,J) = UAVG(IX)
        VRIGHT(IY,J) = UAVG(IY)
        IF(NDIM.EQ.3) VLEFT(IZ,J) = UAVG(IZ)
      ENDDO

C     Right Eigenvectors (Matrix (Ns+1) x (Ns+1...Ns+4)) ..

      VRIGHT(IE,IE) = DENS * U_DOT_S
      VRIGHT(IE,IX) = DENS * (UAVG(IE)/ABAR + U_DOT_M )
      VRIGHT(IE,IY) = DENS * (UAVG(IE)/ABAR - U_DOT_M )
      IF(NDIM.EQ.3) VRIGHT(IE,IZ) = DENS * U_DOT_T
C
C
      DO 10 JDIM = 1 , NDIM
C
         IDIM = JDIM + 2
         Nxyz = VECM(JDIM)
         Sxyz = VECS(JDIM)
         Txyz = VECT(JDIM)
C
         IDIM = JDIM + NSP + 1
C
C     Left Eigenvectors (Matrix (Ns+1...Ns+4) x (Ns+2...Ns+4)) ..
C
         VLEFT(IE,IDIM) = DENSINV * Sxyz
         VLEFT(IX,IDIM) = HALF*DENSINV*(-DE*ABARINV*UAVG(IDIM)+Nxyz)
         VLEFT(IY,IDIM) = HALF*DENSINV*(-DE*ABARINV*UAVG(IDIM)-Nxyz)
         IF(NDIM.EQ.3) VLEFT(IZ,IDIM) = DENSINV * Txyz
C     Right Eigenvectors (Matrix (Ns+2...Ns+4) x (Ns+1...Ns+4)) ..
             VRIGHT(IDIM,IE) = DENS * Sxyz
             VRIGHT(IDIM,IX) = DENS * (UAVG(IDIM)/ABAR+Nxyz)
             VRIGHT(IDIM,IY) = DENS * (UAVG(IDIM)/ABAR-Nxyz)
         IF(NDIM.EQ.3) VRIGHT(IDIM,IZ) = DENS * Txyz
C
C
   10 CONTINUE
#if 0
C
C    Test if the product of the two eigenvectors equals the
C    identity matrix: VRIGHT*VLEFT = Identity
C
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA,
     +                      'Left ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA,
     +                      'Right ',IFAIL)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,VRIGHT,LDA,
     +                 VLEFT,LDA,ZERO,WORK1,MAXNOFVAR)
            IF( UNITMAT(WORK1,NORDER,NORDER,MAXNOFVAR,1.D-15) .EQV. 
     &          .FALSE. ) THEN
                CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,MAXNOFVAR,
     +                      'L*R = Identity matrix ????',IFAIL)
                PAUSE
            ENDIF
#endif
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
!           KPOS(I,J) = SUM1 ! comment out once the analytical ones will be all right
!           KNEG(I,J) = SUM2 ! comment out once the analytical ones will be all right
   15 CONTINUE
C
#endif
C
CREM  VECN is not dIYided by d ?!?! 
C
C
      NX = VECN(1) * AREAINV
      NY = VECN(2) * AREAINV
      IF( NDIM .EQ. 3 )NZ = VECN(3) * AREAINV
C
      FSPSLP = LPOS(IX)+LPOS(IY)
      FSMSLP = LPOS(IX)-LPOS(IY)
      FSPSLN = LNEG(IX)+LNEG(IY)
      FSMSLN = LNEG(IX)-LNEG(IY)
      UPLUS = HALF*FSPSLP - LPOS(1)
      UMINS = HALF*FSPSLN - LNEG(1)
      APLUS = HALF*FSMSLP
      AMINS = HALF*FSMSLN
      AMACHN = U_DOT_M/ABAR
      DELTA = GM1 * HALF
      TRATIO = ONE+DELTA*MACHSQR
      CE = DE*UAVG(IE)/ASQR
!      CD = ZERO
!      CDQ = ZERO
!      DO I = 1 , NSP
!        CD = ALPHA(I) * DR(I) + CD
!        CDQ = ALPHA(I) * DR(I)**2 + CDQ
!      ENDDO
!      CD = KINETIC * TWO * CD / ASQR
!      CDQ = CDQ / (DE * ASQR) 
C
C
C     (Matrix (1...Ns) x (1...Ns)) ..
      DO I = 1 , NSP
        DO J = 1 , NSP
            KPOS(I,J) = ALPHA(I)* (UPLUS*DR(J)/ASQR-AMACHN*APLUS)
            IF (I.EQ.J) THEN
                KPOS(I,J) = KPOS(I,J) + LPOS(1)
            ENDIF
        ENDDO
      ENDDO
C
C     (Matrix (1...Ns) x (Ns+1...Ns+3)) ..
      DO I = 1 , NSP
       KPOS(I,IE) = ALPHA(I)*DE/ASQR*UPLUS
       KPOS(I,IX) = ALPHA(I)*(-DE*UAVG(IX)/ASQR*UPLUS+NX/ABAR*APLUS)
       KPOS(I,IY) = ALPHA(I)*(-DE*UAVG(IY)/ASQR*UPLUS+NY/ABAR*APLUS)
      ENDDO
C
C     (Matrix (Ns+1...Ns+3) x (1...Ns)) .. 
      DO J = 1 , NSP
        KPOS(IE,J) = (DR(J)*UAVG(IE)/ASQR - U_DOT_M*U_DOT_M)*UPLUS 
     &          + (DR(J) - UAVG(IE))*APLUS*AMACHN    
!rpepe        KPOS(IE,J) = HALF * FSPSLP*(DR(J)*UAVG(IE)/ASQR-U_DOT_M*U_DOT_M)
!rpepe     &          + APLUS * AMACHN * (DR(J) - UAVG(IE)) 
!rpepe     &          + LPOS(1) * (U_DOT_M*U_DOT_M - CD + CDQ - DR(J)/DE)
        KPOS(IX,J) = (DR(J)/ASQR*UAVG(IX)-U_DOT_M*NX)*UPLUS
     &          + (DR(J)/ABAR*NX-AMACHN*UAVG(IX))*APLUS
        KPOS(IY,J) = (DR(J)/ASQR*UAVG(IY)-U_DOT_M*NY)*UPLUS
     &          + (DR(J)/ABAR*NY-AMACHN*UAVG(IY))*APLUS
      ENDDO        
C
C     (Matrix (Ns+1...Ns+3) x (Ns+1...Ns+3)) ..     
      KPOS(IE,IE) = LPOS(1)+CE*UPLUS+APLUS*DE*AMACHN
C
      KPOS(IE,IX) = (U_DOT_M * NX-CE * UAVG(IX)) * UPLUS
     &          + APLUS * (UAVG(IE)/ABAR*NX-DE*AMACHN*UAVG(IX))
      KPOS(IE,IY) = (U_DOT_M*NY-CE*UAVG(IY))*UPLUS
     &          + APLUS * (UAVG(IE)/ABAR*NY-DE*AMACHN*UAVG(IY))
C
      KPOS(IX,IE) = DE/ABAR*(UAVG(IX)/ABAR*UPLUS+NX*APLUS)
      KPOS(IY,IE) = DE/ABAR*(UAVG(IY)/ABAR*UPLUS+NY*APLUS)
C
      KPOS(IX,IX) = (NX*NX-DE*UAVG(IX)*UAVG(IX)/ASQR)*UPLUS
     &          + (UAVG(IX)*NX-DE*NX*UAVG(IX))*APLUS/ABAR + LPOS(1)
      KPOS(IX,IY) = (NX*NY-DE*UAVG(IX)*UAVG(IY)/ASQR)*UPLUS
     &          + (UAVG(IX)*NY-DE*NX*UAVG(IY))*APLUS/ABAR
C
      KPOS(IY,IX) = (NY*NX-DE*UAVG(IY)*UAVG(IX)/ASQR)*UPLUS
     &          + (UAVG(IY)*NX-DE*NY*UAVG(IX))*APLUS/ABAR
      KPOS(IY,IY) = (NY*NY-DE*UAVG(IY)*UAVG(IY)/ASQR)*UPLUS
     &          + (UAVG(IY)*NY-DE*NY*UAVG(IY))*APLUS/ABAR + LPOS(1)

C     KNEG

C     (Matrix (1...Ns) x (1...Ns)) ..
      DO I = 1 , NSP
        DO J = 1 , NSP
            KNEG(I,J) = ALPHA(I)*(UMINS*DR(J)/ASQR-AMACHN*AMINS)
            IF (I.EQ.J) THEN
                KNEG(I,J) = KNEG(I,J) + LNEG(1)
            ENDIF
        ENDDO
      ENDDO
C
C     (Matrix (1...Ns) x (Ns+1...Ns+3)) ..
      DO I = 1 , NSP
       KNEG(I,IE) = ALPHA(I)*DE/ASQR*UMINS
       KNEG(I,IX) = ALPHA(I)*(-DE*UAVG(IX)/ASQR*UMINS+NX/ABAR*AMINS)
       KNEG(I,IY) = ALPHA(I)*(-DE*UAVG(IY)/ASQR*UMINS+NY/ABAR*AMINS)
      ENDDO
C
C     (Matrix (Ns+1...Ns+3) x (1...Ns)) .. 
      DO J = 1 , NSP
        KNEG(IE,J) = (DR(J)*UAVG(IE)/ASQR - U_DOT_M*U_DOT_M)*UMINS
     &          + (DR(J) - UAVG(IE))*AMINS*AMACHN
!rpepe        KNEG(IE,J) = HALF*FSPSLN*(DR(J)*UAVG(IE)/ASQR - U_DOT_M*U_DOT_M)
!rpepe     &          + AMINS * AMACHN * (DR(J) - UAVG(IE)) 
!rpepe     &          + LNEG(1) * (U_DOT_M * U_DOT_M - CD + CDQ - DR(J)/DE)
        KNEG(IX,J) = (DR(J)/ASQR*UAVG(IX)-U_DOT_M*NX)*UMINS
     &          + (DR(J)/ABAR*NX-AMACHN*UAVG(IX))*AMINS
        KNEG(IY,J) = (DR(J)/ASQR*UAVG(IY)-U_DOT_M*NY)*UMINS
     &          + (DR(J)/ABAR*NY-AMACHN*UAVG(IY))*AMINS
      ENDDO
C
C     (Matrix (Ns+1...Ns+3) x (Ns+1...Ns+3)) ..
      KNEG(IE,IE) = LNEG(1)+CE*UMINS+AMINS*DE*AMACHN
      KNEG(IE,IX) = (U_DOT_M*NX-CE*UAVG(IX))*UMINS
     &          + AMINS * (UAVG(IE)/ABAR*NX-DE*AMACHN*UAVG(IX))
      KNEG(IE,IY) = (U_DOT_M*NY-CE*UAVG(IY))*UMINS
     &          + AMINS * (UAVG(IE)/ABAR*NY-DE*AMACHN*UAVG(IY))
C
      KNEG(IX,IE) = DE/ABAR*(UAVG(IX)/ABAR*UMINS+NX*AMINS)
      KNEG(IY,IE) = DE/ABAR*(UAVG(IY)/ABAR*UMINS+NY*AMINS)
C
      KNEG(IX,IX) = (NX*NX-DE*UAVG(IX)*UAVG(IX)/ASQR)*UMINS
     &          + (UAVG(IX)*NX-DE*NX*UAVG(IX))*AMINS/ABAR + LNEG(1)
      KNEG(IX,IY) = (NX*NY-DE*UAVG(IX)*UAVG(IY)/ASQR)*UMINS
     &          + (UAVG(IX)*NY-DE*NX*UAVG(IY))*AMINS/ABAR
C
      KNEG(IY,IX) = (NY*NX-DE*UAVG(IY)*UAVG(IX)/ASQR)*UMINS
     &          + (UAVG(IY)*NX-DE*NY*UAVG(IX))*AMINS/ABAR
      KNEG(IY,IY) = (NY*NY-DE*UAVG(IY)*UAVG(IY)/ASQR)*UMINS
     &          + (UAVG(IY)*NY-DE*NY*UAVG(IY))*AMINS/ABAR + LNEG(1)
C
      IF(NDIM.EQ.3)THEN
C
      DO I = 1 , NSP
       KPOS(I,IZ) = ALPHA(I)*(-DE*UAVG(IZ)/ASQR*UPLUS+NZ/ABAR*APLUS)
      ENDDO
C
      KPOS(IE,IZ) = (U_DOT_M*NZ-CE*UAVG(IZ))*UPLUS
     &          + APLUS * (UAVG(IE)/ABAR*NZ-DE*AMACHN*UAVG(IZ))
      KPOS(IX,IZ) = (NX*NZ-DE*UAVG(IX)*UAVG(IZ)/ASQR)*UPLUS
     &          + (UAVG(IX)*NZ-DE*NX*UAVG(IZ))*APLUS/ABAR
      KPOS(IY,IZ) = (NY*NZ-DE*UAVG(IY)*UAVG(IZ)/ASQR)*UPLUS
     &          + (UAVG(IY)*NZ-DE*NY*UAVG(IZ))*APLUS/ABAR
C
      DO J = 1 , NSP
        KPOS(IZ,J) = (DR(J)/ASQR*UAVG(IZ)-U_DOT_M*NZ)*UPLUS
     &          + (DR(J)/ABAR*NZ-AMACHN*UAVG(IZ))*APLUS
      ENDDO  
C
      KPOS(IZ,IE) = DE/ABAR*(UAVG(IZ)/ABAR*UPLUS+NZ*APLUS)
C
      KPOS(IZ,IX) = (NZ*NX-DE*UAVG(IZ)*UAVG(IX)/ASQR)*UPLUS
     &          + (UAVG(IZ)*NX-DE*NZ*UAVG(IX))*APLUS/ABAR
      KPOS(IZ,IY) = (NZ*NY-DE*UAVG(IZ)*UAVG(IY)/ASQR)*UPLUS
     &          + (UAVG(IZ)*NY-DE*NZ*UAVG(IY))*APLUS/ABAR
      KPOS(IZ,IZ) = (NZ*NZ-DE*UAVG(IZ)*UAVG(IZ)/ASQR)*UPLUS
     &          + (UAVG(IZ)*NZ-DE*NZ*UAVG(IZ))*APLUS/ABAR + LPOS(1)
C
      DO I = 1 , NSP
       KNEG(I,IZ) = ALPHA(I)*(-DE*UAVG(IZ)/ASQR*UMINS+NZ/ABAR*AMINS)
      ENDDO
C
      KNEG(IE,IZ) = (U_DOT_M*NZ-CE*UAVG(IZ))*UMINS
     &          + AMINS * (UAVG(IE)/ABAR*NZ-DE*AMACHN*UAVG(IZ))
      KNEG(IX,IZ) = (NX*NZ-DE*UAVG(IX)*UAVG(IZ)/ASQR)*UMINS
     &          + (UAVG(IX)*NZ-DE*NX*UAVG(IZ))*AMINS/ABAR
      KNEG(IY,IZ) = (NY*NZ-DE*UAVG(IY)*UAVG(IZ)/ASQR)*UMINS
     &          + (UAVG(IY)*NZ-DE*NY*UAVG(IZ))*AMINS/ABAR
C
      DO J = 1 , NSP
        KNEG(IZ,J) = (DR(J)/ASQR*UAVG(IZ)-U_DOT_M*NZ)*UMINS
     &          + (DR(J)/ABAR*NZ-AMACHN*UAVG(IZ))*AMINS
      ENDDO
C      
      KNEG(IZ,IE) = DE/ABAR*(UAVG(IZ)/ABAR*UMINS+NZ*AMINS)
C
      KNEG(IZ,IX) = (NZ*NX-DE*UAVG(IZ)*UAVG(IX)/ASQR)*UMINS
     &          + (UAVG(IZ)*NX-DE*NZ*UAVG(IX))*AMINS/ABAR
      KNEG(IZ,IY) = (NZ*NY-DE*UAVG(IZ)*UAVG(IY)/ASQR)*UMINS
     &          + (UAVG(IZ)*NY-DE*NZ*UAVG(IY))*AMINS/ABAR
      KNEG(IZ,IZ) = (NZ*NZ-DE*UAVG(IZ)*UAVG(IZ)/ASQR)*UMINS
     &          + (UAVG(IZ)*NZ-DE*NZ*UAVG(IZ))*AMINS/ABAR + LNEG(1)
C
C
      ENDIF ! NDIM
C
#ifdef FSPL_USE_DEBUG
C
      DO 32 I =  1, NORDER
         DO 32 J =  1, NORDER
            WORK1(I,J) = KPOS(I,J)-KPLUS(I,J)
            WORK2(I,J) = KNEG(I,J)-KMINS(I,J)
   32 CONTINUE
      IF( NULLMAT(WORK1,NORDER,NORDER,MAXNOFVAR,1.D-15) .EQV. .FALSE. 
     2    .OR.
     +    NULLMAT(WORK2,NORDER,NORDER,MAXNOFVAR,1.D-15) .EQV. .FALSE. 
     4   )THEN
         WRITE(6,*)'Errors on Inflow parameters for IELEM =',ielem
         WRITE(6,*)NULLMAT(WORK1,NORDER,NORDER,MAXNOFVAR,1.D-15)
         WRITE(6,*)NULLMAT(WORK2,NORDER,NORDER,MAXNOFVAR,1.D-15)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS,LDA,
     +               'A+ analytical positive matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPLUS,MAXNOFVAR,
     +               'R diag(+) L matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,MAXNOFVAR,
     +               'error positive matrix (A(+) - Rdiag(+)L)',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG,LDA,
     +               'R diag(-) L matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KMINS,MAXNOFVAR,
     +               'A(-) analytical negative matrix',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK2,MAXNOFVAR,
     +               'error negative matrix (A(-) - Rdiag(-)L)',IFAIL)
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
      LFLAG = ( .NOT. NULLMAT(WORK2,NORDER,NORDER,MAXNOFVAR,1.D-15) )
!     lflag = .true.
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
#if 1
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VLEFT,LDA,
     +                  'LEFT eigenvectors',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,VRIGHT,LDA,
     +                  'RIGHT eigenvectors',IFAIL)
C
C     compute L*K and K*R
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,VLEFT,LDA,
     +                 KMAT,LDA,ZERO,WKSPA,MAXNOFVAR)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,ONE,KMAT,LDA,
     +                 VRIGHT,LDA,ZERO,WKSPB,MAXNOFVAR)
C
C     compute R*\Lambda and L*\Lambda 
C
            DO 8 I =  1, NORDER
               DO 8 J =  1, NORDER
                  WKSPC(I,J) = WR(I)*VLEFT(I,J)
                  WKSPD(I,J) = WR(J)*VRIGHT(I,J)
                  WORKL(I,J) = WKSPA(I,J) - WKSPC(I,J)
                  WORKR(I,J) = WKSPB(I,J) - WKSPD(I,J)  
    8       CONTINUE
C
C     test if L*K = \Lambda*L and K*R = R*\Lambda 
C
!            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,MAXNOFVAR,'LEFT*A ',
!     2                  IFAIL)
!            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,MAXNOFVAR,
!     +                  'lambda*LEFT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WORKL,MAXNOFVAR,
     +                  'error LEFT (LEFT*A - labda*LEFT)',IFAIL)
           PAUSE
!            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,MAXNOFVAR,'A*RIGHT',
!     6                  IFAIL)
!            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,MAXNOFVAR,
!     +                  'RIGHT*lambda',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WORKR,MAXNOFVAR,
     +                  'error RIGHT (A*RIGHT - RIGHT*labda)',IFAIL)
           PAUSE
#endif
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     +                  'Analytical Jacobian matrix (A)',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,MAXNOFVAR,
     +                  'Reassembled Jacobian matrix (A(+)+A(-))',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,MAXNOFVAR,
     +                  'Error matrix ((A(+)+A(-))- A)',IFAIL)
           PAUSE
      ENDIF
#endif
C
      RETURN
      END


C     PRESSURE DERIVATIVE (CHEMICAL SPECIES DENSITY)
      DOUBLE PRECISION FUNCTION PIR(u,NDIM,chii,kappa)
      IMPLICIT NONE
      INCLUDE 'constants.h'
      INTEGER NDIM
      DOUBLE PRECISION u(NDIM),chii,kappa

      PIR = chii + HALF * kappa * (u(1)*u(1) + u(2)*u(2))
      
      IF(NDIM.EQ.3)THEN
         PIR = PIR + HALF * kappa * u(3)*u(3)
      ENDIF
      RETURN
      END

C     PRESSURE DERIVATIVE (TOTAL ENERGY)
      DOUBLE PRECISION FUNCTION PIE(kappa)
      IMPLICIT NONE
      DOUBLE PRECISION kappa

      PIE = kappa
      RETURN
      END

C     PRESSURE DERIVATIVE (MOMENTUM)
      DOUBLE PRECISION FUNCTION PIM(u,kappa)
      IMPLICIT NONE
      DOUBLE PRECISION u,kappa

      PIM = - u * kappa
            RETURN
      END
