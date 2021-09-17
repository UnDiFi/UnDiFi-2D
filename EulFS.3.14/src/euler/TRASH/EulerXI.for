      SUBROUTINE EulerXI(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     +                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     +                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Unsteady equations in conserved variables ..
C
C     $Id: EulerXI.f,v 1.6 2011/12/30 11:26:48 abonfi Exp abonfi $
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'bnd.h'
C
C     NEQMAX is the max. no. of equations (5 in 3D)
C            for the matrix scheme (solves for conserved vars.)
C     MAXNOFEQN is the max. no. of mean flow equations (5 in 3D)
C
      INTEGER NEQMAX,LNNVV
      DOUBLE PRECISION TOLER
      PARAMETER (NEQMAX=5,TOLER=1.D-15)
      PARAMETER (LNNVV=NMAX*NMAX*MAXNOFVERT*MAXNOFVERT)
      INTEGER FRSTEQ
      PARAMETER(FRSTEQ=1)
C
C
      INCLUDE 'three.com'
      INCLUDE 'time.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
      INCLUDE 'bodyf.com'
C
C
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR
C
C
      EXTERNAL ScalarScheme,MatrixScheme
C
C
      INTEGER IVAR,IVERT,JVERT,NTURB,IDIM,IADD,JADD,KADD
      INTEGER NORDER,ifail,M,N,MN,N4,JCOL,I,J
      DOUBLE PRECISION FLUCT,SCALRES,ALPHa,HELP
      LOGICAL LFLAG,PICARD
C
C
      DOUBLE PRECISION VCZ(*),VCN(*),VOLUME,
     +                 STIFEL(*),NODRES(*),TSTEP(*),VCB(*)
C
C
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the RHS PETSc vector
C
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the DT PETSc vector 
C
C     CONSV(1:NOFEQN,1:NOFVERT) is used to store the vector
C         of symmetrizing variables
C
C     DCONSV(1:NOFEQN,1:NOFVERT) is used to store the change
C         in the vector of symmetrizing variables
C
      DOUBLE PRECISION CONSV(MAXNOFEQN*MAXNOFVERT),
     2                DCONSV(MAXNOFEQN*MAXNOFVERT),
     3       TAUX(MAXNOFEQN*MAXNOFVERT),SOURCE(MAXNOFVAR),
     4       Jacobian(NEQMAX,NEQMAX*3),
     5       TEMPA((MAXNOFEQN*MAXNOFVERT)**2),
     6       TEMPB((MAXNOFEQN*MAXNOFVERT)**2),
     7       BETA(MAX_NOFVERT_SQR*MAXNOFEQN*MAXNOFEQN),
     8       PHI(NEQMAX),RESIDUAL(2*NEQMAX),
     9       STIFC(MAXNOFVERT*MAXNOFVERT*NEQMAX*NEQMAX)
C
C     NOFEQN (= DIM+2) is actual no. of mean flow equations
C     NORDER (= DIM+2) is actual no. of equations being solved
C                      with the system scheme
C
      INTEGER NOFEQN,ORDSQR
      INTEGER INFO
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
      DOUBLE PRECISION DDOT,dnrm2
      EXTERNAL DDOT,dnrm2
C
      EXTERNAL MatSplitNum,MatSplitXI
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
C
C     Statement function
C
      INTEGER JADDR
      JADDR(IVERT,JVERT,N) = (((JVERT-1)*NOFVERT+IVERT-1)*N*N) + 1
C
      NORDER = NDIM + 2
      NOFEQN = NDIM + 2
      ORDSQR = NOFEQN*NOFEQN
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
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
      CALL Eigen_XI(Jacobian,NEQMAX,dVdZ,dUdV,NDIM,NOFEQN)
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
C     withing subr. eigen_XI
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
      CALL PARM2CONS(ZAVG,DVDZ,NOFEQN,NDIM)
      CALL DGEMM('NoTranspose','NoTranspose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dVdZ,NOFEQN,VCZ,NOFVAR,ZERO,CONSV,NOFEQN)
C
C     CALL X04CAF('G',' ',NOFEQN,Nofvert,CONSV,NOFEQN,
C    +      'Nodal values of the CONS variables ',IFAIL)
C     pause
C
      IF(LALE)THEN
         HELP = ONE/REAL(NOFVERT)
         CALL DINIT(NOFEQN*NDIM,ZERO,TEMPB(1),1)
         DO IVERT = 1, NOFVERT
C
C     copy GRADZ into TEMPA (this is repeated NOFVERT times, could it be avoided?)
C
            DO J = 1, NDIM
               JADD = (J-1)*NOFEQN
               DO I = 1, NOFEQN
                  IADD = JADD+I
                  TEMPA(IADD) = HELP * GRAD_PARM(I,J) ! TEMPA := \nabla_Z/(d+1) is a (NOFEQN x NDIM) matrix
               ENDDO
            ENDDO
!
!        IF(IVERT.EQ.1)
!    &   CALL X04CAF('G',' ',NOFEQN,NDIM,GRAD_PARM,MAXNOFVAR,'nabla Z ',
!    &IFAIL)
C
            IADD = NOFVAR*(IVERT-1)+1
            JADD = NDIM*(IVERT-1)+1
            ALPHA = -HALF/(NDIM*VOLUME) 
            CALL DGER(NOFEQN,NDIM,ALPHA,VCZ(IADD),1, ! TEMPA := \nabla Z + alpha * Z_ivert n_ivert^t
     &                      VCN(JADD),1,TEMPA,NOFEQN)
!
!        CALL X04CAF('G',' ',NOFEQN,NDIM,TEMPA,NOFEQN,
!    +      'local matrix (after DGER)',IFAIL)
!
            CALL PARM2CONS(VCZ(IADD),DVDZ,NOFEQN,NDIM) ! compute (dUdZ)_j
!
!        CALL X04CAF('G',' ',NOFEQN,NDIM,DVDZ,NOFEQN,
!    +      'dUdZ matrix (locale)',IFAIL)
!
            CALL DGEMM('N','N',NOFEQN,NDIM,NOFEQN,ONE,DVDZ,NOFEQN,
     &                       TEMPA,NOFEQN,ONE,TEMPB,NOFEQN) ! sum up within TEMPB
!
!        CALL X04CAF('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,
!    +      'ddu matrix (parziale)',IFAIL)
         ENDDO ! end loop over vertices
!
!        CALL X04CAF('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,
!    +      'ddu matrix ',IFAIL)
!     SOURCE(1) = ZERO
!     SOURCE(2) = ZERO
!     SOURCE(3) = ZERO
!     SOURCE(4) = ZERO
!     IF(NDIM.EQ.3)SOURCE(5) = ZERO
C
C     store ddu * b in SOURCE and multiply by the volume
C
         CALL DGEMV('No',NOFEQN,NDIM,VOLUME,TEMPB(1),NOFEQN,BAVG,1,
     &              ZERO,SOURCE(1),1)
            INFO = 0
caldo if(dnrm2(NOFEQN,SOURCE,1).GT.1.d-6)then
caldo    write(6,*)'ie = ',ielem,(bavg(i),i=1,ndim)
!
caldo    CALL X04CAF('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,
caldo+      'ddu matrix ',IFAIL)
caldo    CALL X04CAF('G',' ',NOFEQN,1,SOURCE,NOFEQN,
caldo+      'b*ddu values ',IFAIL)
caldo       INFO = 1
caldo pause
caldo endif
      ENDIF ! Arbitrary Lagrangean Eulerian
C
C     ---------- Matrix scheme ----------
C
!     SOURCE(2) = -GRAV(1)*VOLUME
!     SOURCE(3) = -GRAV(2)*VOLUME
!     IF(NDIM.EQ.3)SOURCE(4) = -GRAV(3)*VOLUME
C
C REM: calling with TAUX(1) will add new contributions
C      to those already computed when solving entropy
C
      CALL MatrixScheme(MatSplitXI,CONSV(1),DCONSV(1),TAUX(1),
     +                  BETA,STIFEL,NORDER,NORDER,NOFVERT,VCN,
     +                  NDIM,Jacobian,NEQMAX,RESIDUAL,SOURCE,IELEM,
     +                  PICARD)
caldo if(info.NE.0)then
caldo    CALL X04CAF('G',' ',2*NOFEQN,1,RESIDUAL,2*NOFEQN,
caldo+      'residual values ',IFAIL)
caldo    CALL X04CAF('G',' ',NOFEQN,1,source,NOFEQN,
caldo+      'source values ',IFAIL)
caldo  endif
C     
C     copy the timestep from TAUX into TSTEP (update)
C
      DO 12 IVERT = 1, NOFVERT
         IADD = (IVERT-1)*NOFVAR+1
         JADD = (IVERT-1)*NORDER+1
         TSTEP(IADD) = TSTEP(IADD) + TAUX(JADD)
   12 CONTINUE
C
!     CALL X04CAF('G',' ',NOFEQN,nofvert,DCONSV,
!    +            NOFEQN,'Nodal update in U (before)',IFAIL)
C
C     add DCONSV to NODRES
C     could we call MatrixScheme passing NODRES ?????
C     should we really add or just insert ?!?!
C
      DO 16 IVERT = 1, NOFVERT
         IADD = (IVERT-1)*NOFEQN
         JADD = (IVERT-1)*NOFVAR
         DO 16 IVAR = 1, NOFEQN
            NODRES(JADD+IVAR) = NODRES(JADD+IVAR) + DCONSV(IADD+IVAR)
   16 CONTINUE
C
      IF(LTIME)THEN
         DTVOL = DELT/VOLUME
         DO IVERT = 1, NOFVERT
             IADD = (IVERT-1)*NOFVAR+1 ! must be NOFVAR for RANS when NOFVAR != NOFEQN
             JADD = (IVERT-1)*ORDSQR+1 
             CALL PARM2CONS(VCZ(IADD),TEMPA(JADD),NOFEQN,NDIM) ! might re-use when LALE = .TRUE.
         ENDDO
C
!     CALL X04CAF('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz ',IFAIL)
!          pause
C
         IF(DUALTS)
     &   CALL UNSTEADY4(TEMPA,BETA,VCZ,NOFVAR,NODRES,STIFEL,NORDER,NDIM,
     &                  NOFVERT,PICARD)
C
!     CALL X04CAF('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz ',IFAIL)
C
      ENDIF
C
!     write(6,*)'ielem = ',ielem
!     CALL X04CAF('G',' ',Norder,norder*nofvert,beta,Norder,
!    +      'beta(U) ',IFAIL)
C
!     CALL X04CAF('G',' ',NOFEQN,nofvert,NODRES,
!    +            NOFEQN,'Nodal update in U',IFAIL)
      IF(PICARD)CALL DSCAL(N4,-ONE,STIFEL,1)
C
C     CALL X04CAF('G',' ',Norder,2,residual,Norder,
C    +      'residual ',IFAIL)
C     CALL X04CAF('G',' ',Nofvar,Nofvert,consv,Nofvar,
C    +      'cons. vars. ',IFAIL)
!     CALL X04CAF('G',' ',Nofvar,Nofvert,dconsv,Nofvar,
!    +      'update ',IFAIL)
!     pause
C
C --------------- Debugging code starts here ---------------
C
C    Checks the decomposition ..
C
      IF (ICHECK.NE.0) THEN
          CALL DSCAL(2*NORDER,ONE/VOLUME,RESIDUAL,1)
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
         CALL DAXPY(NOFEQN,-ONE/VOLUME,TEMPB(IADD+1),1,TEMPA,1)
   15 CONTINUE
caldo do ivar = 1,nofeqn
caldo   write(6,*)ivar,TEMPA(IVAR),residual(ivar),divflux(ivar)
caldo enddo
      CALL TEST(DivFlux,TEMPA,TOLER,-100*IELEM,NOFEQN)
      PAUSE 'check passed'
      GOTO 747
C
C
caldo CALL X04CAF('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
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
!     CALL X04CAF('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
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

89999 FORMAT (5X,'Scalar residual in Element ',I6,' Wave # ',I1)
99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerVII')

      END 
