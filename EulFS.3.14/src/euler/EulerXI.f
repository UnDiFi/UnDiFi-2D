!> \par Purpose
!>
!> Unsteady equations in conserved variables \f$\rho\left(1,E,\mathbf{u}\right)\f$
!> 
!> @param[in] IELEM the current simplicial element
!> @param[in] VCN the NDIM cartesian components of the inward face normal, scaled by its measure in the NOFVERT vertices
!> @param[in] VCB the NDIM cartesian components of the grid velocity
!> @param[in] VCZ the NOFVAR dofs in the NOFVERT vertices
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVERT the nof vertices of the current simplex (=NDIM+1)
!> @param[in] NOFVAR the nof dependent variables
!> @param[in] NTURB the nof turbulent variables
!> @param[out] NODRES the NOFVAR components of the residual vector in the NOFVERT vertices
!> @param[out] TSTEP is the contribution of the current cell to the global array 
!> @param[out] STIFEL is the contribution of the current cell to the global jacobian matrix
!> @param[in] VOLUME is the area/volume of the current cell
!> @param[in] PICARD is .TRUE. if the (approximate) jacobian matrix has to be assembled
!> @param[in] ScalarScheme the subroutine implementing the FS scheme for scalar equations
!> @param[in] MatrixScheme the subroutine implementing the FS scheme for hyperbolic system
C
      SUBROUTINE EulerXI(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C
C     $Id: EulerXI.f,v 1.16 2020/03/28 09:51:15 abonfi Exp $
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
      INTEGER IWRK,FRSTEQ
      PARAMETER(IWRK=10,FRSTEQ=1)
      INTEGER MOVE(IWRK)
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
      DOUBLE PRECISION FLUCT,SCALRES,ALPHA,HELP,DIVB
      LOGICAL LFLAG,PICARD
C
C
      DOUBLE PRECISION VCZ(*),VCN(*),VOLUME(*),
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
     4       Jacobian(NEQMAX,NEQMAX*3),UMEAN(MAXNOFVAR),
     5       TEMPA((MAXNOFEQN*MAXNOFVERT)**2),
     6       TEMPB((MAXNOFEQN*MAXNOFVERT)**2),
     6       TEMPC((MAXNOFEQN*MAXNOFVERT)**2),
     6       TEMPD((MAXNOFEQN*MAXNOFVERT)**2),
     7       BETA(MAX_NOFVERT_SQR*MAXNOFEQN*MAXNOFEQN),
     8       PHI(NEQMAX),RESIDUAL(2*NEQMAX),
     9       STIFC(MAXNOFVERT*MAXNOFVERT*NEQMAX*NEQMAX)
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
      DOUBLE PRECISION DDOT,dnrm2,DIV
      EXTERNAL DDOT,dnrm2,DIV
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
      IF (ICHECK.NE.0)THEN
C
         WRITE(6,*)'*******************************************' 
         WRITE(6,*)' IE = ',IELEM
         WRITE(6,*)'*******************************************' 
C
C    Some initializations ....
C
         CALL DINIT(NORDER,ZERO,PHI,1)
         CALL CHECK(IELEM,NDIM,NOFEQN) ! here we compute DivFlux
C
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C     when doing ALE the averaged grid velocity must be subtracted
C     from the diagonal entries of the Jacobian (this is done inside Eigen_XI) 
C
C     REM: dUdV is NOT set within Eigen_XI (it should equal the Identity matrix)
C
         CALL Eigen_XI(Jacobian,NEQMAX,dVdZ,dUdV,NDIM,NOFEQN)
C
C
C     COMPUTES THE RESIDUAL/VOLUME as:
C     PHI := dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C          - <b_x> * dU/dX - <b_y> * dU/dy - [ <b_z> * dU/dz ] 
C     for debugging purposes
C
         DO 11 idim = 1,NDIM
             jcol = (idim-1)*NEQMAX + 1
             CALL DGEMV('N',NORDER,NORDER,ONE,Jacobian(1,jcol),NEQMAX,
     +               GRAD_CHAR(FrstEq,idim),1,ONE,PHI,1)
    
   11    CONTINUE
C
C     Here we subtract from the Eulerian flux (DivFlux) the term
C     dUdx (NOFEQN x NDIM) b (NDIM x 1 ) where dUdx is the gradient of the conserved
C     variables and b the averaged grid velocity
C
C                                      |
C                                      |
C                                      V 
         IF(LALE)THEN
            CALL DGEMV('No',NOFEQN,NDIM,MONE,GRAD_CHAR(FrstEq,1),
     &              LDW,BAVG,1,ONE,DIvFlux,1)
         ENDIF
C
C     passing the following test ensures that the gradient of the
C     conserved variables and the jacobian matrix
C     in conserved variables are computed correctly
C     withing subr. eigen_XI
C
         write(6,*)'Here we test that (A-bI)gradU = divflux - gradU*b'
         CALL TEST(DivFlux,PHI,TOLER,IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
      ENDIF
C
C
C     The quasi linear form we discretize is in
C     conservation variables CONSV = (r,rE,ru,rvrw)
C
C     nodal values of the
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
      IF(LTIME)DTVOL = DELT/VOLUME(1)
      IF(LALE)THEN
         GOTO 76
!
!        there is smthg still wrong with this approach
!
         CALL ALEFLUX(NDIM,NOFVERT,NOFEQN,VCN,VCB,VCZ,TEMPA,.FALSE.) ! this is the ALE flux, which is stored in TEMPA
!
!        once the "exact" ALE flux is available, we add b_i \frac{\partial U}{\partial x_i} and
!        distribute
!
         CALL DGEMM('NoTranspose','NoTranspose',NORDER,NDIM,NORDER,
     +              ONE,DVDZ,NORDER,GRAD_PARM,MAXNOFVAR,ZERO,TEMPB,
     3              NORDER) ! store \nabla U in TEMPB
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,GRAD_CHAR,LDW,'nabla U ',
!    &IFAIL)
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,'nabla U(1) ',
!    &IFAIL)
         CALL DGEMV('No',NOFEQN,NDIM,VOLUME(1),TEMPB,NOFEQN,BAVG,1,
     &              ZERO,SOURCE,1)
         CALL DAXPY(NORDER,MONE,TEMPA,1,SOURCE,1)
         IF(ICHECK.NE.0)THEN
            CALL DGEMV('No',NOFEQN,NDIM,ONE,TEMPB,NOFEQN,BAVG,1,
     &              ONE,PHI,1)
            CALL DAXPY(NORDER,MONE/VOLUME(1),TEMPA,1,PHI,1)
         ENDIF
C
         GOTO 67
   76    CONTINUE
C
C        Deconinck's stuff
C
         HELP = ONE/REAL(NOFVERT)
         ALPHA = -HALF/(NDIM*VOLUME(1)) 
         CALL DINIT(NOFEQN*NDIM,ZERO,TEMPB,1)
         CALL DINIT(NOFEQN,ZERO,UMEAN,1)
C
         DO IVERT = 1, NOFVERT ! loop over the vertices
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
!    &   CALL R8Mat_Print('G',' ',NOFEQN,NDIM,GRAD_PARM,MAXNOFVAR,'nabla Z ',
!    &IFAIL)
C
            IADD = NOFVAR*(IVERT-1)+1
            JADD = NDIM*(IVERT-1)+1
            CALL DGER(NOFEQN,NDIM,ALPHA,VCZ(IADD),1, ! TEMPA := \nabla Z + alpha * Z_ivert n_ivert^t
     &                      VCN(JADD),1,TEMPA,NOFEQN)
!
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,TEMPA,NOFEQN,
!    +      'local matrix (after DGER)',IFAIL)
!
            CALL PARM2CONS(VCZ(IADD),DVDZ,NOFEQN,NDIM) ! compute (dUdZ)_j
C
            CALL DGEMV('No',NOFEQN,NOFEQN,ONE,DVDZ,NOFEQN,VCZ(IADD),1, ! compute <<U>> = [\sum_j (1/2) (dUdZ)_j x Z_j]/(d+1)
     &                 ONE,UMEAN(1),1)
!
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,DVDZ,NOFEQN,
!    +      'dUdZ matrix (locale)',IFAIL)
!
            CALL DGEMM('N','N',NOFEQN,NDIM,NOFEQN,ONE,DVDZ,NOFEQN,
     &                       TEMPA,NOFEQN,ONE,TEMPB,NOFEQN) ! sum up within TEMPB
!
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,
!    +      'ddu matrix (parziale)',IFAIL)
         ENDDO ! end loop over vertices
C
!
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,TEMPD,NOFEQN,
!    +      'ddu(1) matrix ',IFAIL)
!        CALL R8Mat_Print('G',' ',NOFEQN,NDIM,TEMPB,NOFEQN,
!    +      'ddu matrix ',IFAIL)
!      pause
C
C     store ddu * b in SOURCE and multiply by the volume
C
         CALL DGEMV('No',NOFEQN,NDIM,VOLUME(1),TEMPB(1),NOFEQN,BAVG,1,
     &              ZERO,SOURCE(1),1)
!        write(6,*)'ie (b) = ',ielem,(source(i),i=1,nofeqn)
C
C     add the term which depends upon div(b)
C
         IF(.NOT.DUALTS)THEN
            DIVB = DIV(NDIM,NOFVERT,VCN,VCB) ! computes the divergence of the grid velocity (divb is multiplied by the volume)
            ALPHA = -HALF*DIVB/REAL(NOFVERT) ! 1/(2*(d+1)) is part of <<U>>
!
!        write(6,*)'ie (U) = ',ielem,(umean(i),i=1,nofeqn)
!
            CALL DAXPY(NOFEQN,ALPHA,UMEAN,1,SOURCE,1)
         ENDIF
!
!        write(6,*)'ie (a) = ',ielem,(source(i),i=1,nofeqn)
C
         IF(ICHECK.NE.0)THEN
            CALL DGEMV('No',NOFEQN,NDIM,ONE,TEMPB(1),NOFEQN,BAVG,1,
     &              ONE,PHI,1)
            CALL DAXPY(NOFEQN,ALPHA/VOLUME(1),UMEAN,1,PHI,1)
         ENDIF
c
   67 CONTINUE
      ELSE
         SOURCE(1) = ZERO
         SOURCE(2) = ZERO
         SOURCE(3) = ZERO
         SOURCE(4) = ZERO
         IF(NDIM.EQ.3)SOURCE(5) = ZERO
      ENDIF ! on LALE
C
C     CALL R8Mat_Print('G',' ',NOFEQN,Nofvert,CONSV,NOFEQN,
C    +      'Nodal values of the CONS variables ',IFAIL)
C     pause
C
C     ---------- Matrix scheme ----------
C
C
      CALL MatrixScheme(MatSplitXI,CONSV(1),DCONSV(1),TAUX(1),
     +                  BETA,STIFEL,NORDER,NORDER,NOFVERT,VCN,NDIM,
     +                  Jacobian,NEQMAX,RESIDUAL,SOURCE,IELEM,PICARD)
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
               HELP = MAX(HELP,TAUX(JADD))
            ELSE
               HELP = HELP + TAUX(JADD)
            ENDIF
         ENDDO
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NORDER
            JADD = IADD + IVAR
            TSTEP(JADD) = HELP
         ENDDO
      ENDDO
C
C
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,DCONSV,
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
         DO IVERT = 1, NOFVERT
             IADD = (IVERT-1)*NOFVAR+1 ! must be NOFVAR for RANS when NOFVAR != NOFEQN
             JADD = (IVERT-1)*ORDSQR+1 
             CALL PARM2CONS(VCZ(IADD),TEMPA(JADD),NOFEQN,NDIM)
         ENDDO
      ENDIF
C
      IF(PICARD)THEN
C
C     compute the transformation matrices from
C     conserved to parameter variables in the vertices
C
          DO 33 IVERT = 1,NOFVERT
              IADD = (IVERT-1)*NOFVAR+1
              JADD = (IVERT-1)*NOFEQN*NOFEQN+1
c
c a bit of care here: Z and dZdU are dimensioned NOFVAR in
c subr MatdZdU(), but this should not be a problem if called
c with NOFEQN
c
          CALL MatdZdU(VCZ(IADD),dZdU(JADD),NDIM,NOFEQN)
   33     CONTINUE
      ENDIF
C
C
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz before',IFAIL)
!          pause
C
C
      IF(LTIME.AND.DUALTS)THEN
          DTVOL = DELT/VOLUME(1)
          CALL UNSTEADY4(TEMPA,BETA,VCZ,NOFVAR,NODRES,STIFEL,VOLUME,
     &                   NORDER,NDIM,NOFVERT,PICARD)
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,tempa,Norder,
!    +      'dudz after ',IFAIL)
      ENDIF 
C
C
C
!     write(6,*)'ielem = ',ielem
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,beta,Norder,
!    +      'beta(U) ',IFAIL)
C
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +            NOFEQN,'Nodal update in U',IFAIL)
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
          IF(LALE)THEN
              CALL CHECK(IELEM,NDIM,NOFEQN) ! here we re-compute DivFlux (it is NOT volume-multiplied)
              CALL ALEFLUX(NDIM,NOFVERT,NOFEQN,VCN,VCB,VCZ,TEMPA,.TRUE.)! this is the ALE flux
              ALPHA = MONE/VOLUME(1)
              CALL DAXPY(NOFEQN,ALPHA,TEMPA,1,DivFlux,1) 
          ENDIF
!
!     The following test is NOT going to work when using Deconinck's ALE formulation
!
          write(6,*)'Here we test the ALE fluxes in cell ',IELEM
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
C     Compute: TEMPB := 2. * avg(dU/dZ) * (dZ/dU)_j
C
C     attenzione perche' dVdZ e' stata sovrascritta
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,TWO,dVdZ,NOFEQN,dZdU,
     +               NOFEQN,ZERO,TEMPB,NOFEQN)
c
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
!    +      'TEMPB ',IFAIL)
!     pause
C
C     Here we re-use beta
C
      CALL DCOPY(N4,STIFEL,1,BETA,1)
C
          DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFEQN*NOFEQN + 1
          DO 35 IVERT = 1,NOFVERT
                  IADD = JADDR(IVERT,JVERT,NOFVAR)
!                 KADD = JADDR(IVERT,JVERT,NOFEQN)
C
               CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +              NOFEQN,NOFEQN,ONE,BETA(IADD),NOFVAR,
     +              TEMPB(JADD),NOFEQN,ZERO,STIFEL(IADD),NOFVAR)
C
C     now STIFEL contains the convection stiffness matrix in
C     conserved variables
C
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,beta(iadd),Nofvar,
!    +      'C(i,j) ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
!    +      'Cp(i,j) ',IFAIL)
C
   35     CONTINUE
      CALL DSCAL(N4,MONE,STIFEL,1)
!     pause
C
C
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

99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerXI')

      END
