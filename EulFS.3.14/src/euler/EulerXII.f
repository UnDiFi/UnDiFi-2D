!> \copydetails EulerIX()
      SUBROUTINE EulerXII(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Unsteady equations in conserved variables ..
C
C     $Id: EulerXII.f,v 1.7 2020/03/28 09:51:15 abonfi Exp $
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'bnd.h'
      INCLUDE 'plasma.h'
C
C     NEQMAX is the max. no. of equations (5 in 3D)
C            for the matrix scheme (solves for conserved vars.)
C     MAXNOFEQN is the max. no. of mean flow equations (5 in 3D)
C
C
      DOUBLE PRECISION SOURCECHEM(NSP),PRESS,MAXDAM,SOURCEEN
      INTEGER NEQMAX,LNNVV
      DOUBLE PRECISION TOLER
      PARAMETER (NEQMAX=MAXNOFEQN,TOLER=1.D-15)
      PARAMETER (LNNVV=NMAX*NMAX*MAXNOFVERT*MAXNOFVERT)
      INTEGER IWRK,FRSTEQ
      PARAMETER(IWRK=10,FRSTEQ=1)
      INTEGER MOVE(IWRK)
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
      INCLUDE 'bodyf.com'
      INCLUDE 'commonchem.inc'
      INCLUDE 'tauchem.com'
C
C
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,ISP
C
C
      EXTERNAL ScalarScheme,MatrixScheme
C
C
      INTEGER IVAR,IVERT,JVERT,NTURB,IDIM,IADD,JADD,KADD
      INTEGER NORDER,ifail,M,N,N4,JCOL
      DOUBLE PRECISION FLUCT,SCALRES,ESQR
      LOGICAL LFLAG,PICARD
C
C
      DOUBLE PRECISION VCZ(*),VCN(*),VOLUME(*),
     +                 STIFEL(*),NODRES(*),TSTEP(*),VCB(*)
C
C
C
      DOUBLE PRECISION CONSV(MAXNOFEQN*MAXNOFVERT),
     2                DCONSV(MAXNOFEQN*MAXNOFVERT),
     3       TAUX(MAXNOFEQN*MAXNOFVERT),SOURCE(MAXNOFVAR),
     4       Jacobian(NEQMAX,NEQMAX*3),
     5       TEMPA((MAXNOFEQN*MAXNOFVERT)**2),
     6       TEMPB((MAXNOFEQN*MAXNOFVERT)**2),
     7       BETA(MAX_NOFVERT_SQR*MAXNOFEQN*MAXNOFEQN),
     8       PHI(NEQMAX),RESIDUAL(2*NEQMAX),DAMK(NSP),
     9       STIFC(MAXNOFVERT*MAXNOFVERT*NEQMAX*NEQMAX)
         DOUBLE PRECISION HELP,SUMK,ULOCAL,HEIGHT 
C
C     NOFEQN (= DIM+2) is actual no. of mean flow equations
C     NORDER (= DIM+2) is actual no. of equations being solved
C                      with the system scheme
C
      INTEGER NOFEQN,ORDSQR
C
C     RESIDUAL[1:NORDER] stores the residual computed by
C                        the Matrix scheme as \sum K_j U_j
C     RESIDUAL[NORDER+1:2*NORDER]
C                        stores the residual computed by
C                        the Matrix scheme as \sum C_{ij} U_j
C     it is used just for debugging purposes, to be compared with
C     the residual computed as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C
C
      DOUBLE PRECISION DDOT,DNRM2,PRESS4Ar,OHMHEAT
      EXTERNAL DNRM2,DDOT,OHMHEAT
      EXTERNAL PRESS4Ar
!      INTEGER CHEM
C
      EXTERNAL MatSplitNum,MatSplitXII
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
C
C     Statement function
C
      INTEGER JADDR
      JADDR(IVERT,JVERT,N) = (((JVERT-1)*NOFVERT+IVERT-1)*N*N) + 1
C
C
      NORDER = NDIM + NSP + 1
      NOFEQN = NDIM + NSP + 1
      ORDSQR = NOFEQN*NOFEQN
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
!     CALL R8Mat_Print('G',' ',NOFVAR,Nofvert,VCZ,NOFVAR,
!    +      'Nodal values of the Z variables ',IFAIL)
      CALL PARM2PRIM4Ar(NDIM,IELEM)
!     write(6,*)'<Z>',ielem,(ZAVG(jcol),jcol=1,NORDER)
!     write(6,*)'<U>',ielem,(UAVG(jcol),jcol=1,NORDER)
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
C
C     set local residual and timestep to zero (should maybe bring it
C     in the calling routine)
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DCONSV,1)
      CALL DINIT(NOFVERT*NORDER,ZERO,TAUX,1)
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.EQ.0) GOTO 7
C
C    Some initializations ....
C
      CALL DINIT(NORDER,ZERO,PHI,1)
      CALL CHECK(IELEM,NDIM,NOFEQN)
C
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_XII(Jacobian,NEQMAX,dVdZ,dUdV,NDIM,NOFEQN)
C
C
C     COMPUTES THE RESIDUAL/VOLUME as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C     for debugging purposes
C
      DO 11 idim = 1,NDIM
          jcol = (idim-1)*NEQMAX + 1
          CALL DGEMV('N',NORDER,NORDER,ONE,Jacobian(1,jcol),NEQMAX,
     +               GRAD_CHAR(FrstEq,idim),1,ONE,PHI,1)
    
   11 CONTINUE
C
C     passing this test ensures that the gradient of the
C     conserved variables and the jacobian matrix
C     in conserved variables are computed correctly
C     withing subr. eigen_XII
C
C     CALL TEST(DivFlux,PHI,TOLER,IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
    7 CONTINUE
C
C
C     The quasi linear form we discretize is in
C     conservation variables CONSV = (r,rE,ru,rvrw)
C
C
C     conserved variables are computed through the relation
C     \U_j = M Z_j
C
C     where M = \frac{\partial U}{\partial Z}
C
C
      CALL PARM2CONS4Ar(ZAVG,DVDZ,NOFEQN,NDIM)
      CALL DGEMM('NoTranspose','NoTranspose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dVdZ,NOFEQN,VCZ,NOFVAR,ZERO,CONSV,NOFEQN)
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,dVdZ,NOFEQN,
!    +      'Matrice di trasformazione',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,Nofvert,CONSV,NOFEQN,
!    +      'Nodal values of the CONS variables ',IFAIL)
!     pause
C
C     ---------- Matrix scheme ----------
C
C     If CHEM=1 the Source term is ON
C     If CHEM=~1 the Source term is OFF
C
!     *************************************
!     N.B. Inserire if(Plasma) esterno 
!          a if(Chem.eq.1) 
!     ************************************

!      IF(PLASMA)THEN
C
      IF (CHEM.EQ.1)THEN
C
          PRESS = PRESS4Ar(NDIM,ZAVG)      
C
          CALL CHEMSOURCE(ZAVG,PRESS,SOURCECHEM)       
C
          DO ISP = 1 , NSP
              SOURCE(ISP) = -SOURCECHEM(ISP)*VOLUME(1)
!              write(6,*)'SOURCE=',ISP,SOURCE(ISP)
          ENDDO
          IF (OHM.EQ.1)THEN              
              ESQR = DDOT(NDIM,GRAD_PARM(NOFVAR,1),NOFVAR,
     &                         GRAD_PARM(NOFVAR,1),NOFVAR)
              SOURCEEN = OHMHEAT(ZAVG,PRESS,ESQR,NDIM,NOFVAR)
              SOURCE(IE) = -SOURCEEN*VOLUME(1)
          ELSE         
              SOURCE(IE) = ZERO
          ENDIF         
          SOURCE(IX) = ZERO
          SOURCE(IY) = ZERO
          IF(NDIM.EQ.3)SOURCE(IZ) = ZERO
C
      ELSE ! CHEM = 0 
          DO ISP = 1 , NOFVAR
              SOURCE(ISP) = ZERO
          ENDDO
      ENDIF  
C
C
      CALL MatrixScheme(MatSplitXII,CONSV(1),DCONSV(1),TAUX(1),
     +                  BETA,STIFEL,NORDER,NORDER,NOFVERT,VCN,NDIM,
     +                  Jacobian,NEQMAX,RESIDUAL,SOURCE,IELEM,PICARD)
C
C
C     compute the timestep
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         IADD = (IVERT-1)*NORDER
         HELP = ZERO 
C
C     sum over the dofs
C
         DO IVAR = 1,NORDER
            JADD = IADD + IVAR
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD)) ! take the maximum
            ELSE
               HELP = HELP + TAUX(JADD) ! sum up
            ENDIF
         ENDDO
         DO IVAR = 1,NORDER
            JADD = IADD + IVAR
            TAUX(JADD) = HELP
         ENDDO
      ENDDO ! loop over vertices
      IF( CHAR_TIMESTEPPING )THEN  
  
!      OPEN(12,FILE='fort.12')
C
!      cell lenght scale
C      h = (Vol * |u|)/(sum k_j^+)
         SUMK = ZERO
         DO IVERT = 1, NOFVERT
            HELP = DDOT(NDIM,VCN((IVERT-1)*NDIM+1),1,UAVG(IX),1)
            SUMK = SUMK + MAX(ZERO,HELP)
         ENDDO
         HELP = REAL(NDIM)
         SUMK = SUMK/HELP ! SUMK = (u \cdot n_j )/d
         ULOCAL = DNRM2(NDIM,UAVG(IX),1)
         HEIGHT = VOLUME(1) * ULOCAL / SUMK
!      write(13,*)IELEM,'cell size=',HEIGHT
!      pause
!      HEIGHT = ONE

!     Damkohler number  (Da_i = S_i*L/(u*\rho_i))
!        MAXDAM=10.0d0
         DO ISP = 1, NSP
!           DAMK(ISP)=SOURCECHEM(ISP)/(ULOCAL*UAVG(ISP))        
!           DAMK(ISP)=abs(DAMK(ISP))*HEIGHT ! /(1.0d0 - 1.0d0/exp(1.0d0)
!           WRITE(6,*)'DAMK(',ISP,') before=',DAMK(ISP)           
           DAMK(ISP)=HEIGHT/ULOCAL/TAU(ISP)    
!           WRITE(6,*)'DAMK(',ISP,') after=',DAMK(ISP)                     
         ENDDO

!         pause
!       WRITE(12,*)IELEM,HEIGHT,MAXDAM
!      WRITE(12,*)IELEM,(DAMK(ISP),ISP=1,NSP)
      
C
C     Time step scaling      
C
         DO IVERT = 1, NOFVERT
            IADD = (IVERT-1)*NOFEQN
            DO ISP = 1,NSP 
               IF (DAMK(ISP).GT.ONE)THEN
                  TAUX(IADD+ISP) = TAUX(IADD+ISP)*DAMK(ISP)
               ENDIF
            ENDDO
         ENDDO
      ENDIF ! test on char_timestepping

!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,TAUX,
!    +            NOFEQN,'Nodal timestep ',IFAIL)
C
C     add DCONSV to NODRES and
C     copy the timestep from TAUX into TSTEP (update)
C     could we call MatrixScheme passing NODRES ?????
C     should we really add or just insert ?!?!
C
      DO 16 IVERT = 1, NOFVERT
         IADD = (IVERT-1)*NOFEQN
         JADD = (IVERT-1)*NOFVAR
         DO 16 IVAR = 1, NOFEQN
            TSTEP(JADD+IVAR) = TSTEP(JADD+IVAR) + TAUX(IADD+IVAR)
            NODRES(JADD+IVAR) = NODRES(JADD+IVAR) + DCONSV(IADD+IVAR)
   16 CONTINUE
C
      IF(LTIME)THEN
         DTVOL = DELT/VOLUME(1)
         DO IVERT = 1, NOFVERT
             IADD = (IVERT-1)*NOFVAR+1 ! must be NOFVAR for RANS when NOFVAR != NOFEQN
             JADD = (IVERT-1)*ORDSQR+1 
             CALL PARM2CONS4Ar(VCZ(IADD),TEMPA(JADD),NOFEQN,NDIM)
         ENDDO
C
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz ',IFAIL)
!          pause
C
         IF(DUALTS)
     &   CALL UNSTEADY4(TEMPA,BETA,VCZ,NOFVAR,NODRES,STIFEL,VOLUME,
     &                  NORDER,NDIM,NOFVERT,PICARD)
C
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz ',IFAIL)
C
      ENDIF
C
!     write(6,*)'ielem = ',ielem
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,beta,Norder,
!    +      'beta(U) ',IFAIL)
C
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +            NOFEQN,'Nodal update in U',IFAIL)
      IF(PICARD)CALL DSCAL(N4,-ONE,STIFEL,1)
C
C     CALL R8Mat_Print('G',' ',Norder,2,residual,Norder,
C    +      'residual ',IFAIL)
C     CALL R8Mat_Print('G',' ',Nofvar,Nofvert,consv,Nofvar,
C    +      'cons. vars. ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvert,dconsv,Nofvar,
!    +      'update ',IFAIL)
!     pause
C
C --------------- Debugging code starts here ---------------
C
C    Checks the decomposition ..
C
      IF (ICHECK.NE.0) THEN
          CALL DSCAL(2*NORDER,ONE/VOLUME(1),RESIDUAL,1)
          LFLAG = .TRUE.
          DO 18 IVAR = 1,NORDER
              IF (DABS(PHI(IVAR)-RESIDUAL(IVAR)).GT.
     +            TOLER) LFLAG = .FALSE.
   18     CONTINUE
          IF (LFLAG .EQV. .FALSE.) THEN
              WRITE (6,99999) IELEM
              DO 22 IVAR = 1,NORDER
                  WRITE (6,*) PHI(IVAR),RESIDUAL(IVAR),
     +              DABS(PHI(IVAR)-RESIDUAL(IVAR))
   22         CONTINUE
C           PAUSE 'ehi ! you!'
C
          ENDIF
C
C     test the residual as computed by the "explicit" scheme
C
          CALL TEST(DivFlux,RESIDUAL,TOLER,IELEM,NOFEQN)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     Add the element stiffness matrix to the global stiffness matrix
C
C
C     transform the element stiffness matrix into conserved variables
C
      IF(ICHECK.NE.0)THEN
caldo IF( NOFVAR .EQ. NOFEQN )THEN
caldo     CALL DGEMM('No Transpose','No Transpose',NOFVAR,
caldo+               NOFVAR*NOFVERT*NOFVERT,NOFVAR,-ONE,dUdV,
caldo+               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
caldo ELSE
          CALL DINIT(NOFEQN*NOFVAR,ZERO,TEMPB,1)
          DO 14 IVERT = 1, NOFVERT
                  IADD = (IVERT-1)*NOFEQN
          DO 13 JVERT = 1, NOFVERT
                  IF(JVERT.EQ.IVERT)GOTO 13
                  JADD = (JVERT-1)*NOFEQN
                  KADD = JADDR(IVERT,JVERT,NOFVAR)
                  DO IVAR = 1,NOFEQN
                     TEMPA(IADD+IVAR) = CONSV(JADD+IVAR)-
     &                                  CONSV(IADD+IVAR)
                  ENDDO
          CALL DGEMV('NoTranspose',NOFEQN,NOFEQN,ONE,STIFEL(KADD),
     +               NOFVAR,TEMPA(IADD+1),1,ONE,TEMPB(IADD+1),1)
   13     CONTINUE
caldo     do ivar = 1,nofeqn
caldo        write(6,*)ivar,TEMPB(IADD+IVAR),dconsv(iadd+ivar)
caldo     enddo
   14     CONTINUE
caldo ENDIF
      CALL DINIT(NOFEQN*NOFVAR,ZERO,TEMPA,1)
      DO 15 IVERT = 1, NOFVERT
         IADD = (IVERT-1)*NOFEQN
         CALL DAXPY(NOFEQN,MONE/VOLUME(1),TEMPB(IADD+1),1,TEMPA,1)
   15 CONTINUE
caldo do ivar = 1,nofeqn
caldo   write(6,*)ivar,TEMPA(IVAR),residual(ivar),divflux(ivar)
caldo enddo
      CALL TEST(DivFlux,TEMPA,TOLER,-100*IELEM,NOFEQN)
      PAUSE 'check passed'
      GOTO 747
C
C
caldo CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
caldo+      'TEMPB ',IFAIL)
caldo pause
C
          DO 35 IVERT = 1,NOFVERT
                  IADD = JADDR(IVERT,JVERT,NOFVAR)
          DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFEQN*NOFEQN + 1
                  KADD = JADDR(IVERT,JVERT,NOFEQN)
C
C
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
!    +      'C(i,j) ',IFAIL)
C
   35     CONTINUE
!     pause
C
  747 continue
      ENDIF
C
      IF (ICHECK.EQ.0) RETURN
C
C --------------- Debugging code starts here ---------------
C
C     test the residual as computed by the "implicit" scheme
C
      CALL TEST(DivFlux,RESIDUAL(NORDER+1),TOLER,-IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
      RETURN

99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerVII')

      END
