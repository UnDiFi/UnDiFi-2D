!> \copydoc MatSplitIX()
      SUBROUTINE MatSplitVIII(IELEM,NDIM,NORDER,VECN,
     +JacobianMatrix,LDJ,KMAT,KPOS,KNEG,VLEFT,VRIGHT,LDA,
     +WR,LPOS,LNEG,EIGENDECO)
C
C     $Id: matsplVIII.F,v 1.6 2020/03/28 09:51:15 abonfi Exp $
C
C
C     Subroutine to compute the generalized inflow parameters
C     for the UNpreconditioned pseudo-compressible (incompressible)
C     Euler equations
C     Positive and negative inflow parameters are normally
C     computed using their analytical expressions. 
C     if the flag DEBUG is defined they are also computed through
C     their eigenvector decomposition K = R \Lambda L
C     and the two are compared for debugging purposes
C
C
      IMPLICIT NONE
C
C#define DEBUG
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'three.com'
      INCLUDE 'chorin.com'
      INCLUDE 'flags.com'
C
C
      INTEGER   IELEM,NDIM,NORDER,LDA,LDJ
      LOGICAL   EIGENDECO
C
      DOUBLE PRECISION VECN(*)
      DOUBLE PRECISION JacobianMatrix(LDJ,*)
      DOUBLE PRECISION KMAT(LDA,*),KPOS(LDA,*),KNEG(LDA,*),
     >       VLEFT(LDA,*),VRIGHT(LDA,*),WR(*),LPOS(*),LNEG(*)
C
      DOUBLE PRECISION UDOTM,AREAINV,AREA,CSQRINV,CNST,
     >                 CSQR,C,LPLUS,LMINS,ANX,ANY,ANZ
C
      DOUBLE PRECISION VPOSX,VPOSY,VPOSZ,WX,WY,WZ
      DOUBLE PRECISION VNEGX,VNEGY,VNEGZ,L0POS,L0NEG,
     +                 CINV,TEMPX,TEMPY,TEMPZ,TEMPA,TEMPB,TEMP
C
#ifdef DEBUG
      INTEGER I,J,K,IDIM,IFAIL
      LOGICAL LFLAG,WARNA,WARNB
      INTEGER TABLE(3,3),INDX
      DOUBLE PRECISION VECM(3),VECS(3),VECT(3)
      DOUBLE PRECISION WKSPA(NMAX,NMAX),WKSPB(NMAX,NMAX)
      DOUBLE PRECISION WKSPC(NMAX,NMAX),WKSPD(NMAX,NMAX)
      DOUBLE PRECISION WORK1(NMAX,NMAX),WORK2(NMAX,NMAX)
      DOUBLE PRECISION KPLUS(4,4),KMINS(4,4)
      DOUBLE PRECISION SURF,SUM1,SUM2,UDOTS,UDOTT,TEMP1
C
      DATA (table(1,j),j=1,3) / 1,0,0 /
      DATA (table(2,j),j=1,3) / 0,1,0 /
      DATA (table(3,j),j=1,3) / 0,0,1 /
      DATA LFLAG / .FALSE. /
C
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER   ISDMIN
      LOGICAL   UNITMAT,NULLMAT
      EXTERNAL  ISDMIN,DDOT,DNRM2,UNITMAT,NULLMAT
#endif
C
      INTRINSIC DSQRT
C
      DATA ANX,ANY,ANZ / 3*ZERO /
C
C     UDOTM is the dot product between the velocity vector
C           and the face normal
C
      ANX = VECN(1)/NDIM
      ANY = VECN(2)/NDIM
      UDOTM = ZAVG(2) * ANX + ZAVG(3) * ANY
      IF( NDIM .EQ. 3 )THEN
         ANZ = VECN(3)/NDIM
         UDOTM = UDOTM + ZAVG(4) * ANZ
      ENDIF
C
C     Jacobian Matrix for primitive variables ..
C
      KMAT(1,1) = ZERO
      KMAT(1,2) = BETA*ANX
      KMAT(1,3) = BETA*ANY
C
      KMAT(2,1) = ANX
      KMAT(2,2) = ZAVG(2) * ANX + UDOTM
      KMAT(2,3) = ZAVG(2) * ANY
C
      KMAT(3,1) = ANY
      KMAT(3,2) = ZAVG(3) * ANX
      KMAT(3,3) = ZAVG(3) * ANY + UDOTM
C
      IF( NDIM .EQ. 3 )THEN
          KMAT(1,4) = BETA*ANZ
          KMAT(2,4) = ZAVG(2) * ANZ
          KMAT(3,4) = ZAVG(3) * ANZ
          KMAT(4,1) = ANZ
          KMAT(4,2) = ZAVG(4) * ANX
          KMAT(4,3) = ZAVG(4) * ANY
          KMAT(4,4) = ZAVG(4) * ANZ + UDOTM
      ENDIF
C
      IF( EIGENDECO .EQV. .FALSE. )RETURN
C
      AREA = VECN(1)*VECN(1)+VECN(2)*VECN(2)
      IF( NDIM .EQ. 3 )AREA = AREA + VECN(3)*VECN(3)
      AREA = SQRT(AREA)
      AREAINV = ONE/AREA
C
      AREA = AREA/NDIM
C
C     Some useful constants ..
C
      UDOTM = UDOTM / AREA
      CSQR = UDOTM * UDOTM + BETA
      C = DSQRT( CSQR )
      CSQRINV = ONE/CSQR
      CINV = ONE/C
C
      LPLUS = UDOTM + C
      LMINS = UDOTM - C
C
C     Eigenvalues ..
C
      WR(1) = AREA * UDOTM
      WR(2) = AREA * LMINS
      WR(3) = AREA * LPLUS
      WR(4) = AREA * UDOTM
C
#ifdef DEBUG
      IF( NDIM .EQ. 2 )THEN
         VECM(1) = VECN(1) * AREAINV
         VECM(2) = VECN(2) * AREAINV
         VECM(3) = ZERO
         VECS(1) = - VECM(2)
         VECS(2) = VECM(1)
         VECS(3) = ZERO
         VECT(1) = ZERO
         VECT(2) = ZERO
         VECT(3) =  ONE
      ELSE
         VECM(1) = VECN(1) * AREAINV
         VECM(2) = VECN(2) * AREAINV
         VECM(3) = VECN(3) * AREAINV
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
      UDOTM = DDOT(NDIM,ZAVG(2),1,VECM,1)
      UDOTS = DDOT(NDIM,ZAVG(2),1,VECS,1)
      UDOTT = DDOT(NDIM,ZAVG(2),1,VECT,1)
C
C     Left Eigenvectors ..
C
      VLEFT(1,1) =-UDOTS*CSQRINV
      VLEFT(2,1) = HALF * CSQRINV
      VLEFT(3,1) = HALF * CSQRINV
      IF(NDIM.EQ.3)VLEFT(4,1) =-UDOTT*CSQRINV
C
C     Right Eigenvectors ..
C
      VRIGHT(1,1) = ZERO
      VRIGHT(1,2) = LPLUS * C
      VRIGHT(1,3) =-LMINS * C
      IF(NDIM.EQ.3)VRIGHT(1,4) = ZERO
C
C     Left Eigenvectors ..
C
         VLEFT(1,2) = VECS(1) - UDOTM*CSQRINV*UDOTS*VECM(1)
         VLEFT(2,2) = HALF * CSQRINV * LMINS * VECM(1)
         VLEFT(3,2) = HALF * CSQRINV * LPLUS * VECM(1)
C
         VLEFT(1,3) = VECS(2) - UDOTM*CSQRINV*UDOTS*VECM(2)
         VLEFT(2,3) = HALF * CSQRINV * LMINS * VECM(2)
         VLEFT(3,3) = HALF * CSQRINV * LPLUS * VECM(2)
C
         IF( NDIM .EQ. 3 )THEN
            VLEFT(4,2) = VECT(1) - UDOTM*CSQRINV*UDOTT*VECM(1)
            VLEFT(4,3) = VECT(2) - UDOTM*CSQRINV*UDOTT*VECM(2)
            VLEFT(1,4) = VECS(3) - UDOTM*CSQRINV*UDOTS*VECM(3)
            VLEFT(2,4) = HALF * CSQRINV * LMINS * VECM(3)
            VLEFT(3,4) = HALF * CSQRINV * LPLUS * VECM(3)
            VLEFT(4,4) = VECT(3) - UDOTM*CSQRINV*UDOTT*VECM(3)
         ENDIF
C
C     Right Eigenvectors ..
C
         VRIGHT(2,1) = VECS(1)
         VRIGHT(2,2) = ZAVG(2) - LPLUS * VECM(1)
         VRIGHT(2,3) = ZAVG(2) - LMINS * VECM(1)
C
         VRIGHT(3,1) = VECS(2)
         VRIGHT(3,2) = ZAVG(3) - LPLUS * VECM(2)
         VRIGHT(3,3) = ZAVG(3) - LMINS * VECM(2)
C
         IF( NDIM .EQ. 3 )THEN
             VRIGHT(2,4) = VECT(1)
             VRIGHT(3,4) = VECT(2)
C
             VRIGHT(4,1) = VECS(3)
             VRIGHT(4,2) = ZAVG(4) - LPLUS * VECM(3)
             VRIGHT(4,3) = ZAVG(4) - LMINS * VECM(3)
             VRIGHT(4,4) = VECT(3)
         ENDIF 
C
C     Test if VRIGHT*VLEFT = Identity ..
C
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VRIGHT,LDA,
     +                 VLEFT,LDA,0.D0,WORK1,NMAX)
C
         WARNA = UNITMAT( WORK1, NORDER, NORDER, NMAX, 1.D-15 )
C
         IF( .NOT. WARNA )THEN
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WORK1,NMAX,
     +                  'Left * Right eigenvectors',IFAIL)
            PAUSE
         ENDIF
C
C     Positive and negative eigenvalues ..
C
      DO 12 IDIM = 1 , NORDER
         LPOS(IDIM) = HALF * (WR(IDIM) + DABS(WR(IDIM)) )
         LNEG(IDIM) = HALF * (WR(IDIM) - DABS(WR(IDIM)) )
   12 CONTINUE
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
C     here we need a unit normal
C
      ANX = VECN(1) * AREAINV
      ANY = VECN(2) * AREAINV
      IF( NDIM .EQ. 3 )ANZ = VECN(3) * AREAINV
      UDOTM = ANX*ZAVG(2)+ANY*ZAVG(3)
      IF( NDIM .EQ. 3 )UDOTM = UDOTM + ANZ*ZAVG(4)
C
      WR(1) = UDOTM
      WR(2) = LMINS
      WR(3) = LPLUS
      WR(4) = UDOTM
C
      LPOS(1) =HALF* (WR(1)+ABS(WR(1)))
      LPOS(2) =HALF* (WR(2)+ABS(WR(2)))
      LPOS(3) =HALF* (WR(3)+ABS(WR(3)))
      LPOS(4) =HALF* (WR(4)+ABS(WR(4)))
      LNEG(1) =HALF* (WR(1)-ABS(WR(1)))
      LNEG(2) =HALF* (WR(2)-ABS(WR(2)))
      LNEG(3) =HALF* (WR(3)-ABS(WR(3)))
      LNEG(4) =HALF* (WR(4)-ABS(WR(4)))
C
      CNST = AREA * CSQRINV 
C
C     Auxiliary arrays VPOS,VNEG,W
C
      VPOSX = ZAVG(2) * LPLUS + BETA * ANX
      VPOSY = ZAVG(3) * LPLUS + BETA * ANY
      VPOSZ = ZAVG(4) * LPLUS + BETA * ANZ
      VNEGX = ZAVG(2) * LMINS + BETA * ANX
      VNEGY = ZAVG(3) * LMINS + BETA * ANY
      VNEGZ = ZAVG(4) * LMINS + BETA * ANZ
      WX = ZAVG(2) - UDOTM * ANX
      WY = ZAVG(3) - UDOTM * ANY
      WZ = ZAVG(4) - UDOTM * ANZ
C
      TEMP = HALF * BETA * CINV * AREA
C
      KPOS(1,1) = TEMP
      KNEG(1,1) =-TEMP
C
      TEMPX = TEMP * LPLUS
      KPOS(1,2) = TEMPX * ANX
      KPOS(1,3) = TEMPX * ANY
      IF( NDIM.EQ.3 ) KPOS(1,4) = TEMPX * ANZ
C
      TEMPX =-TEMP * LMINS
      KNEG(1,2) = TEMPX * ANX
      KNEG(1,3) = TEMPX * ANY
      IF( NDIM.EQ.3 ) KNEG(1,4) = TEMPX * ANZ
C
      KPOS(2,1) = CNST * (HALF * VPOSX - LPOS(1) * WX )
      KPOS(3,1) = CNST * (HALF * VPOSY - LPOS(1) * WY )
      IF( NDIM .EQ. 3 )
     >KPOS(4,1) = CNST * (HALF * VPOSZ - LPOS(1) * WZ )
C
      KNEG(2,1) = CNST * (HALF * VNEGX - LNEG(1) * WX )
      KNEG(3,1) = CNST * (HALF * VNEGY - LNEG(1) * WY )
      IF( NDIM .EQ. 3 )
     >KNEG(4,1) = CNST * (HALF * VNEGZ - LNEG(1) * WZ )
C
      TEMPX = ( CSQR * (ONE - ANX * ANX) - UDOTM * WX * ANX )
      TEMPY = ( CSQR * (ONE - ANY * ANY) - UDOTM * WY * ANY )
      TEMPZ = ( CSQR * (ONE - ANZ * ANZ) - UDOTM * WZ * ANZ )
C
      L0POS = LPOS(1) * CNST 
      L0NEG = LNEG(1) * CNST 
      TEMPA = HALF * LPLUS * CNST
      TEMPB = HALF * LMINS * CNST
C
      KPOS(2,2) = TEMPA * VPOSX * ANX + L0POS * TEMPX
      KPOS(3,3) = TEMPA * VPOSY * ANY + L0POS * TEMPY
C
      KNEG(2,2) = TEMPB * VNEGX * ANX + L0NEG * TEMPX
      KNEG(3,3) = TEMPB * VNEGY * ANY + L0NEG * TEMPY
C
      IF( NDIM .EQ. 3 )THEN
          KPOS(4,4) = TEMPA * VPOSZ * ANZ + L0POS * TEMPZ
          KNEG(4,4) = TEMPB * VNEGZ * ANZ + L0NEG * TEMPZ
      ENDIF
C
      TEMPY = ( -CSQR * ANX * ANY - UDOTM * WX * ANY )
C
      KPOS(2,3) = TEMPA * VPOSX * ANY + L0POS * TEMPY
      KNEG(2,3) = TEMPB * VNEGX * ANY + L0NEG * TEMPY
C
      IF( NDIM .EQ. 3 )THEN
          TEMPZ = ( -CSQR * ANX * ANZ - UDOTM * WX * ANZ )
          KPOS(2,4) = TEMPA * VPOSX * ANZ + L0POS * TEMPZ
          KNEG(2,4) = TEMPB * VNEGX * ANZ + L0NEG * TEMPZ
      ENDIF
C
      TEMPX = ( -CSQR * ANY * ANX - UDOTM * WY * ANX ) 
      TEMPZ = ( -CSQR * ANY * ANZ - UDOTM * WY * ANZ ) 
C
      KPOS(3,2) = TEMPA * VPOSY * ANX + L0POS * TEMPX
      KNEG(3,2) = TEMPB * VNEGY * ANX + L0NEG * TEMPX
C
      IF( NDIM .EQ. 3 )THEN
          KPOS(3,4) = TEMPA * VPOSY * ANZ + L0POS * TEMPZ
          KNEG(3,4) = TEMPB * VNEGY * ANZ + L0NEG * TEMPZ
      ENDIF
C
      IF( NDIM .EQ. 3 )THEN
          TEMPX = ( -CSQR * ANZ * ANX - UDOTM * WZ * ANX ) 
          TEMPY = ( -CSQR * ANZ * ANY - UDOTM * WZ * ANY )
          KPOS(4,2) = TEMPA * VPOSZ * ANX + L0POS * TEMPX
          KPOS(4,3) = TEMPA * VPOSZ * ANY + L0POS * TEMPY
          KNEG(4,2) = TEMPB * VNEGZ * ANX + L0NEG * TEMPX
          KNEG(4,3) = TEMPB * VNEGZ * ANY + L0NEG * TEMPY
      ENDIF
C
C     the following ensures a correct computation of the timestep
C
      LPOS(1) = LPOS(1) * AREA
      LPOS(2) = LPOS(2) * AREA
      LPOS(3) = LPOS(3) * AREA
      LPOS(4) = LPOS(4) * AREA
C      
#ifdef DEBUG
C
      DO 32 I =  1, NORDER
         DO 32 J =  1, NORDER
            WORK1(I,J) = (KPOS(I,J)-KPLUS(I,J))
            WORK2(I,J) = (KNEG(I,J)-KMINS(I,J))
   32 CONTINUE
      WARNA = NULLMAT( WORK1, NORDER, NORDER, NMAX, 1.D-14 )
      WARNB = NULLMAT( WORK2, NORDER, NORDER, NMAX, 1.D-14 )
      IF( .NOT. WARNA )THEN
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPOS,LDA,'A+ ',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KPLUS,4,'R*K*L ',IFAIL)
         PAUSE
      ENDIF
      IF( .NOT. WARNB )THEN
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KNEG,LDA,'A-',IFAIL)
         CALL R8Mat_Print('G',' ',NORDER,NORDER,KMINS,4,'R*K*L  ',IFAIL)
         PAUSE
      ENDIF
C
C        WARNA = UNITMAT(....1.D-14)
C
         DO 30 I =  1, NORDER
            DO 30 J =  1, NORDER
               WORK1(I,J) = (KPOS(I,J)+KNEG(I,J))
               WORK2(I,J) = WORK1(I,J)-KMAT(I,J)
   30    CONTINUE
         LFLAG = .NOT. NULLMAT(WORK2,NORDER,NORDER,NMAX,1.D-14)
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
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,VLEFT,LDA,
     +                  KMAT,LDA,ZERO,WKSPA,NMAX)
            CALL DGEMM('N','N',NORDER,NORDER,NORDER,1.D0,KMAT,LDA,
     +                  VRIGHT,LDA,ZERO,WKSPB,NMAX)
            DO 16 I =  1, NORDER
               DO 16 J =  1, NORDER
                  WKSPC(I,J) = WR(I)*VLEFT(I,J)
                  WKSPD(I,J) = WR(J)*VRIGHT(I,J)
   16       CONTINUE
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPA,NMAX,
     &                  'LEFT*A ',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPC,NMAX,
     &                  'lambda*LEFT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPB,NMAX,
     &                  'A*RIGHT',IFAIL)
            CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSPD,NMAX,
     &                  'RIGHT*lambda',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,KMAT,LDA,
     &                  'Original matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK1,NMAX,
     +                  'Reassembled matrix',IFAIL)
            CALL R8Mat_Print('General',' ',NORDER,NORDER,WORK2,NMAX,
     +                  'Error matrix',IFAIL)
      ENDIF
C
#endif
C
      RETURN
      END
