      SUBROUTINE Eigen_X(Matrix,LDA,DVDZ,DUDV,DVDU,NDIM,LDB)
C
C     $Id: eigen_X.F,v 1.7 2020/03/28 09:51:14 abonfi Exp $
C
C      This is Merkle's preconditioner
C
C     define FSPL_USE_DEBUG for debugging purposes, i.e. when using
C     -check
C
CCCCCCCC#define FSPL_USE_DEBUG
C
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
C
      INCLUDE 'merkle.com'
      INCLUDE 'stream.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'flags.com'
      INCLUDE 'pfcgas.com'
C
      INTEGER LDA,NDIM,LDB
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB),DVDU(LDB,LDB)
#ifdef FSPL_USE_DEBUG
      DOUBLE PRECISION WORK(MAX_NOFVAR_SQR),WKSP(MAX_NOFVAR_SQR)
      INTEGER IPIV(MAXNOFVAR)
      LOGICAL UNITMAT
#endif
C
      DOUBLE PRECISION DENS,DENSINV,HELP,XI,TEMP,DUMMY
      DOUBLE PRECISION GP11,GP12,GP13,GP14,GP21,GP22,GP31,GP33,GP41,
     2GP42,GP43,GP44,BOX,CIRC
      INTEGER IDIM,ifail
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
      DENS = UAVG(1)
      DENSINV = ONE/UAVG(1) 
!     MINFSQR = M_INFTY*M_INFTY
      DUMMY = ONE/(GAM*RSTAR) ! DUMMY replaces MINFSQR
      HELP = GM1*DUMMY
      BOX = HELP*(AMPSQR-ONE)*DENSINV
      CIRC = HELP*AMPSQR*ASQR
C
C
      IDIM = 1
C
      Matrix(1,1,IDIM) = UAVG(3)*AMPSQR
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = DENS*ASQR*AMPSQR
      Matrix(1,4,IDIM) = ZERO
      Matrix(1,5,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = BOX*UAVG(3)
      Matrix(2,2,IDIM) = UAVG(3)
      Matrix(2,3,IDIM) = CIRC
      Matrix(2,4,IDIM) = ZERO
      Matrix(2,5,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = DENSINV
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(3)
      Matrix(3,4,IDIM) = ZERO
      Matrix(3,5,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = ZERO
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = UAVG(3)
      Matrix(4,5,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = ZERO
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = ZERO
      Matrix(5,4,IDIM) = ZERO
      Matrix(5,5,IDIM) = UAVG(3)
C
      IDIM = 2
C
      Matrix(1,1,IDIM) = UAVG(4)*AMPSQR
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = DENS*AMPSQR*ASQR
      Matrix(1,5,IDIM) = ZERO
C
      Matrix(2,1,IDIM) = BOX*UAVG(4) 
      Matrix(2,2,IDIM) = UAVG(4)
      Matrix(2,3,IDIM) = ZERO
      Matrix(2,4,IDIM) = CIRC
      Matrix(2,5,IDIM) = ZERO
C
      Matrix(3,1,IDIM) = ZERO
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(4)
      Matrix(3,4,IDIM) = ZERO
      Matrix(3,5,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = DENSINV
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = UAVG(4)
      Matrix(4,5,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = ZERO
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = ZERO
      Matrix(5,4,IDIM) = ZERO
      Matrix(5,5,IDIM) = UAVG(4)
C
      IDIM = 3
C
      Matrix(1,1,IDIM) = UAVG(5)*AMPSQR
      Matrix(1,2,IDIM) = ZERO
      Matrix(1,3,IDIM) = ZERO
      Matrix(1,4,IDIM) = ZERO
      Matrix(1,5,IDIM) = DENS*AMPSQR*ASQR
C
      Matrix(2,1,IDIM) = BOX*UAVG(5)
      Matrix(2,2,IDIM) = UAVG(5)
      Matrix(2,3,IDIM) = ZERO
      Matrix(2,4,IDIM) = ZERO
      Matrix(2,5,IDIM) = CIRC
C
      Matrix(3,1,IDIM) = ZERO
      Matrix(3,2,IDIM) = ZERO
      Matrix(3,3,IDIM) = UAVG(5)
      Matrix(3,4,IDIM) = ZERO
      Matrix(3,5,IDIM) = ZERO
C
      Matrix(4,1,IDIM) = ZERO
      Matrix(4,2,IDIM) = ZERO
      Matrix(4,3,IDIM) = ZERO
      Matrix(4,4,IDIM) = UAVG(5)
      Matrix(4,5,IDIM) = ZERO
C
      Matrix(5,1,IDIM) = DENSINV
      Matrix(5,2,IDIM) = ZERO
      Matrix(5,3,IDIM) = ZERO
      Matrix(5,4,IDIM) = ZERO
      Matrix(5,5,IDIM) = UAVG(5)
C
C     Transformation matrix from primitive to conserved variables ..
C     this is matrix \Gamma_p^{-1} of MSC's notes
C
      TEMP = ONE/AMPSQR+GM1
      HELP =-DENS/(ASQR*DUMMY)
C
      DUDV(1,1) = TEMP/ASQR
      DUDV(1,2) = HELP
      DUDV(1,3) = ZERO
      DUDV(1,4) = ZERO
C
      DUDV(2,1) = TEMP*(ONE/GM1+HALF*MACHSQR)-ONE
      DUDV(2,2) = HELP*KINETIC
      DUDV(2,3) = DENS*UAVG(3)
      DUDV(2,4) = DENS*UAVG(4)
C
      DUDV(3,1) = TEMP*UAVG(3)/ASQR
      DUDV(3,2) = HELP*UAVG(3)
      DUDV(3,3) = UAVG(1)
      DUDV(3,4) = ZERO
C
      DUDV(4,1) = TEMP*UAVG(4)/ASQR
      DUDV(4,2) = HELP*UAVG(4)
      DUDV(4,3) = ZERO
      DUDV(4,4) = UAVG(1)
C
      IF(NDIM.EQ.3)THEN
         DUDV(5,1) = TEMP*UAVG(5)/ASQR
         DUDV(5,2) = HELP*UAVG(5)
         DUDV(5,3) = ZERO
         DUDV(5,4) = ZERO
C
         DUDV(1,5) = ZERO
         DUDV(2,5) = DENS*UAVG(5)
         DUDV(3,5) = ZERO
         DUDV(4,5) = ZERO
         DUDV(5,5) = DENS
      ENDIF
C
      IF(.NOT.LTIME)GOTO 120
C
C     Transformation matrix from primitive to conserved variables ..
C
      GP14 = GM1*AMPSQR
      GP11 = KINETIC*GP14
      GP12 =-GP14*UAVG(3)
      GP13 =-GP14*UAVG(4)
C
      TEMP = (ONE+GM1*AMPSQR)*DUMMY 
C
      GP44 = GM1*TEMP*DENSINV
      GP42 =-GP44*UAVG(3)
      GP43 =-GP44*UAVG(4)
      GP41 =-ASQR*DENSINV*DUMMY+GP44*KINETIC
C
      GP21 =-UAVG(3)*DENSINV
      GP22 = DENSINV
C
      GP31 =-UAVG(4)*DENSINV
      GP33 = DENSINV
C
      DVDU(1,1) = GP11
      DVDU(1,2) = GP14
      DVDU(1,3) = GP12
      DVDU(1,4) = GP13
C
      DVDU(3,1) = GP21
      DVDU(3,2) = ZERO
      DVDU(3,3) = GP22
      DVDU(3,4) = ZERO
C
      DVDU(4,1) = GP31
      DVDU(4,2) = ZERO
      DVDU(4,3) = ZERO
      DVDU(4,4) = GP33
C
      DVDU(2,1) = GP41
      DVDU(2,2) = GP44
      DVDU(2,3) = GP42
      DVDU(2,4) = GP43
C
      IF(NDIM.EQ.3)THEN
         DVDU(5,1) =-UAVG(5)*DENSINV
         DVDU(5,2) = ZERO
         DVDU(5,3) = ZERO
         DVDU(5,4) = ZERO
C
         DVDU(1,5) = -UAVG(5)*GP14
         DVDU(2,5) = -GP44*UAVG(5)
         DVDU(3,5) = ZERO
         DVDU(4,5) = ZERO
         DVDU(5,5) = DENSINV
      ENDIF
C
#ifdef FSPL_USE_DEBUG
C
C    Test if the product of the two preconditioners equals the
C    identity matrix: DUDV*DVDU = Identity
C
      CALL DCOPY(LDB*LDB,dUdV,1,WKSP,1)
      CALL DGETRF(LDB,LDB,WKSP,LDB,IPIV,IFAIL)
      CALL DGETRI(LDB,WKSP,LDB,IPIV,WORK,MAXNOFVAR,IFAIL)
                CALL R8Mat_Print('G',' ',LDB,LDB,WKSP,LDB,
     +                      'Inverse via DGETRI ????',IFAIL)
                CALL R8Mat_Print('G',' ',LDB,LDB,DVDU,LDB,
     +                      'Inverse via analitical ????',IFAIL)

            CALL DGEMM('N','N',LDB,LDB,LDB,ONE,DUDV,LDB,
     +                 DVDU,LDB,ZERO,WORK,LDB)
            IF( UNITMAT(WORK,LDB,LDB,LDB,1.D-15) .EQV. .FALSE. )
     +      THEN
                CALL R8Mat_Print('G',' ',LDB,LDB,WORK,LDB,
     +                      'P*P^{-1} = Identity matrix ????',IFAIL)
                PAUSE
            ENDIF
#endif
C
  120 CONTINUE
C
C     Transformation matrix from parameter vector to symmetrizing variables
C
      DVDZ(1,1) = GM1OG*ZAVG(2) ! dp/dz1
      DVDZ(1,2) = GM1OG*ZAVG(1) ! dp/dz2
      DVDZ(1,3) =-GM1OG*ZAVG(3) ! dp/dz3
      DVDZ(1,4) =-GM1OG*ZAVG(4) ! dp/dz4
C
      DVDZ(3,1) =-ZAVG(3)*DENSINV  ! du/dz1
      DVDZ(3,2) = ZERO            ! du/dz2
      DVDZ(3,3) = ONE/ZAVG(1)     ! du/dz3
      DVDZ(3,4) = ZERO            ! du/dz4
C
      DVDZ(4,1) =-ZAVG(4)*DENSINV  ! dv/dz1
      DVDZ(4,2) = ZERO
      DVDZ(4,3) = ZERO
      DVDZ(4,4) = ONE/ZAVG(1)
C
      HELP = ZAVG(3)*ZAVG(3)+ZAVG(4)*ZAVG(4) 
      IF(NDIM.EQ.3)HELP = HELP + ZAVG(5)*ZAVG(5)
      HELP = HELP/ZAVG(1)
      XI = GM1*DUMMY
      TEMP = XI*DENSINV 
C
      DVDZ(2,1) = TEMP*(HELP-ZAVG(2))        ! dT/dz1
      DVDZ(2,2) = XI/ZAVG(1)
      DVDZ(2,3) =-TEMP*ZAVG(3)
      DVDZ(2,4) =-TEMP*ZAVG(4)
C
      IF (NDIM.EQ.2) RETURN
C
      DVDZ(1,5) =-GM1OG*ZAVG(5)
      DVDZ(2,5) =-TEMP*ZAVG(5) 
      DVDZ(3,5) = ZERO
      DVDZ(4,5) = ZERO
C
      DVDZ(5,1) =-ZAVG(5)*DENSINV
      DVDZ(5,2) = ZERO
      DVDZ(5,3) = ZERO
      DVDZ(5,4) = ZERO
      DVDZ(5,5) = ONE/ZAVG(1)
C
      RETURN 
C
      END
