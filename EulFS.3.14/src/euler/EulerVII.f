!> \copydetails EulerIX()
      SUBROUTINE EulerVII(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Unsteady equations in symmetrizing variables ..
C
C     $Id: EulerVII.f,v 1.43 2020/03/28 09:51:15 abonfi Exp $
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'time.h'
C
C     NEQMAX is the max. no. of equations (4 in 3D)
C            for the matrix scheme (solves for dp/ra,du,dv,dw)
C     MAXNOFEQN is the max. no. of mean flow equations (5 in 3D)
C
      INTEGER NEQMAX,LNNVV
      DOUBLE PRECISION TOLER
      PARAMETER (NEQMAX=4,TOLER=1.D-15)
      PARAMETER (LNNVV=MAX_NOFVAR_SQR*MAX_NOFVERT_SQR)
      INTEGER IWRK,FRSTEQ
      PARAMETER(IWRK=10,FRSTEQ=2)
      INTEGER MOVE(IWRK)
C
C
      INCLUDE 'time.com'
      INCLUDE 'three.com'
      INCLUDE 'transf.com'
      INCLUDE 'flags.com'
      INCLUDE 'bodyf.com'
C
C     .. Scalar Arguments ..
      INTEGER IELEM,NDIM,NOFVERT,NOFVAR,NTURB
      DOUBLE PRECISION VOLUME(*)
      LOGICAL PICARD
C
C     .. Array Arguments ..
      DOUBLE PRECISION VCZ(*),VCN(*),VCB(*),STIFEL(*),NODRES(*),TSTEP(*)
C
C     .. External Arguments ..
      EXTERNAL ScalarScheme,MatrixScheme
C
C     .. Local Scalar ..
      INTEGER IVAR,IVERT,JVERT,IDIM,IADD,JADD,KADD
      INTEGER NORDER,ifail,M,N,MN,N4,JCOL
      DOUBLE PRECISION FLUCT,SCALRES,HELP
      LOGICAL LFLAG
      LOGICAL unitmat
C
C
C
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the RHS PETSc vector
C
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the DT PETSc vector 
C
C     SYMMV(1:NOFEQN,1:NOFVERT) is used to store the vector
C         of symmetrizing variables
C
C     DSYMMV(1:NOFEQN,1:NOFVERT) is used to store the change
C         in the vector of symmetrizing variables
C
      DOUBLE PRECISION SYMMV(MAXNOFEQN*MAXNOFVERT),PHI(NEQMAX)
      DOUBLE PRECISION DSYMMV(MAXNOFEQN*MAXNOFVERT),WKSP(MAXNOFEQN),
     1                 TAUX(MAXNOFEQN*MAXNOFVERT),SOURCE(MAXNOFVAR),
     2                 Jacobian(NEQMAX,NEQMAX*3),RESIDUAL(2*NEQMAX),
     3                 TEMPA(3*MAX_NOFVERT_SQR*MAXNOFEQN**2),
     4                 TEMPB(3*MAX_NOFVERT_SQR*MAXNOFEQN**2)
      DOUBLE PRECISION BETA1(NEQMAX*NEQMAX*MAXNOFVERT),
     1                 BETA(2*MAXNOFEQN*MAXNOFEQN*MAXNOFVERT),
     2                 STIFC(MAX_NOFVERT_SQR*NEQMAX*NEQMAX)
C
C     Need to dimension BETA with 2*MAXNOF etc. otherwise
C     run-time error in 3D; maybe a bug somewhere
C
C     NOFEQN (= DIM+2) is actual no. of mean flow equations
C     NORDER (= DIM+1) is actual no. of equations being solved
C                      with the system scheme
C     NOFVAR (=NOFEQN) when the are no turbulence models
C
      INTEGER NOFEQN
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
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C
      EXTERNAL MatSplitNum,MatSplitVII
C
      DATA SOURCE/MAXNOFVAR*ZERO/
C
C     Statement function
C
      INTEGER JADDR
      JADDR(IVERT,JVERT,N) = (((JVERT-1)*NOFVERT+IVERT-1)*N*N) + 1
C
      NORDER = NDIM + 1 ! we solve for dp/ra,du,dv,dw
      NOFEQN = NDIM + 2 ! this is for r,rE,ru,rv,rw
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
      IF (LTIME) THEN
         DTVOL = DELT/VOLUME(1)
         IF(DUALTS)THEN 
            CALL DINIT(NOFEQN*NOFEQN*NOFVERT,ZERO,BETA,1)
            CALL DINIT(NORDER*NORDER*NOFVERT,ZERO,BETA1,1)
         ENDIF
      ENDIF
C
C     set local residual and timestep to zero (should maybe bring it
C     in the calling routine)
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DSYMMV,1)
      CALL DINIT(NOFVERT*NOFEQN,ZERO,TAUX,1)
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C    Some initializations ....
C
          CALL DINIT(NORDER,ZERO,PHI,1)
          CALL DINIT(NOFEQN,ZERO,WKSP,1)
          CALL CHECK(IELEM,NDIM,NOFEQN)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_VII(Jacobian,NEQMAX,dVdZ,dUdV,NDIM,NOFEQN)
C
      IF (ICHECK.EQ.0) GOTO 7
C
C --------------- Debugging code starts here ---------------
C
C
C     COMPUTES THE RESIDUAL/VOLUME as:
C     dF/dU * dU/dX + dG/dU * dU/dy + [ dH/dU * dU/dz ]
C     for debugging purposes
C
      DO 12 idim = 1,NDIM
          jcol = (idim-1)*NEQMAX + 1
          CALL DGEMV('N',NORDER,NORDER,ONE,Jacobian(1,jcol),NEQMAX,
     +               GRAD_CHAR(FrstEq,idim),1,ONE,PHI,1)

   12 CONTINUE
C
C --------------- Debugging code ends here ---------------
C
    7 CONTINUE
C
C     The quasi linear form we discretize is in
C     symmetrizing variables SYMMV = (dS,dp/(r*a),du,dv,dw)
C
C
C     symmetrizing variables are computed through the relation
C     \tilde{U}_i = M Z_j
C     where M = \frac{\partial \tilde{U}}{\partial Z}
C     in what follows we actually compute
C      \tilde{U}^t = Z^t M^t where ' denotes transposition
C     so that symmetrizing variables are stored as
C     SYMMV(1:NOFVERT,1:NOFEQN) and entropy
C     appears in the first NOFVERT locations, i.e.
C
C     dS(1) dp/ra(1) du(1) dv(1) dw(1)
C     ....    ....    ..    ..    ..
C     ....    ....    ..    ..    ..
C     dS(4) dp/ra(4) du(4) dv(4) dw(4)    
C
      CALL DGEMM('Transpose','Transpose',NOFVERT,NOFEQN,NOFEQN,
     +           ONE,VCZ,NOFVAR,dVdZ,NOFEQN,ZERO,SYMMV,NOFVERT)
C
!     CALL R8Mat_Print('G',' ',Nofvar,NOFvert,vcz,Nofvar,
!    +      'Z variables ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +      NOFEQN,'Nodal update in U on ENTRY ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvert,NOFEQN,symmv,Nofvert,
!    +      'SYMM variables ',IFAIL)
C     pause
C
C     ---------- Scalar scheme on the entropy wave ----------
C
C     Entropy wave ..
C
      IVAR = 1
C
C The entropy eqn. is solved using an upwind scheme
C
      CALL ScalarScheme(IELEM,VCN,R_SPEED(1,IVAR),SCALRES,ZERO,
     +                  SYMMV(1),TAUX(1),DSYMMV(1),BETA1,STIFC,NDIM,
     +                  NOFVERT,PICARD)
!     write(6,*)ielem
!     CALL R8Mat_Print('G',' ',Nofvert,Nofeqn,symmv,Nofvert,
!    +      'symm. vars. after entropy ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvert,Nofeqn,dsymmv,Nofvert,
!    +      'residual in symm. vars. after entropy ',IFAIL)
C
C     the first NOFVERT entries of TAUX now keep
C     the timestep computed by the Scalar Scheme
C
C
C
C    Copy the convective jacobian (STIFC) into STIFEL and
C    reset it to 0. since it will be reused by the matrix scheme
C
      IF (LTIME.AND.DUALTS) THEN
!     CALL R8Mat_Print('G',' ',Nofvert,1,beta1,Nofvert,
!    +      'vettore beta per entropia',IFAIL)
! bug somewhere: unless the following dcopy is commented out, one of
!                the entries in SYMMV is set to 0.d0
          CALL DCOPY(NOFVERT*NOFVERT,BETA1,1,BETA,NOFVAR*NOFVAR)
!         CALL MATINS(BETA,NOFEQN,BETA1,1,NOFVERT,NOFVERT,0)
          CALL DINIT(NOFVERT*NOFVERT,ZERO,BETA1,1) ! maybe useless
      ENDIF
      IF (PICARD) THEN
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL,NOFVAR*NOFVAR)
          CALL DINIT(NOFVERT*NOFVERT,ZERO,STIFC,1)
      ENDIF
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C    Checks for the scalar residual ..
C
          FLUCT = DDOT(NDIM,R_SPEED(1,IVAR),1,GRAD_CHAR(IVAR,1),LDW)
C
          SCALRES = SCALRES/VOLUME(1)
          CALL DAXPY(NOFEQN,SCALRES,dUdV((IVAR-1)*NOFEQN+1),1,WKSP,1)
C
          IF (DABS(FLUCT-SCALRES).GT.TOLER) THEN
              WRITE (6,89999) IELEM,IVAR
              WRITE (6,FMT=*)FLUCT,SCALRES,ABS(FLUCT-SCALRES)
          ENDIF

      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C     the matrix of symmetrizing variables (except entropy) is now
C     transposed so that these are stored (starting at SYMMV(M+1))
C     as SYMMV(1:NORDER,1:NOFVERT)
C
      M = NOFVERT
      N = NORDER
      MN = NOFVERT*NORDER
C
!     CALL R8Mat_Print('G',' ',Nofvert,Norder,symmv(M+1),Nofvert,
!    +      'last symm. vars before calling TRANS ',IFAIL)
C
      CALL TRANS(SYMMV(M+1),NOFVERT,NORDER,MN,MOVE,IWRK,IFAIL)
!
!     CALL R8Mat_Print('G',' ',NOFVERT,NOFEQN,TAUX,Nofvert,
!    +      'vettore DT before TRANS',IFAIL)
!
C
C     at this stage the matrix SYMM is as follows:
C
C     dS(1) dp/ra(1) ...  dp/ra(4)
C     ....   du(1)   ...   du(4)
C            dv(1)   ...   dv(4) 
C     dS(4)  dw(1)   ...   dw(4)
C
C
C     ---------- Matrix scheme ----------
C
      SOURCE(2) = -GRAV(1)*VOLUME(1)
      SOURCE(3) = -GRAV(2)*VOLUME(1)
      IF(NDIM.EQ.3)SOURCE(4) = -GRAV(3)*VOLUME(1)
C
C REM: calling with TAUX(NOFVERT+1) will add new contributions
C      in the correct locations
C
!     CALL R8Mat_Print('G',' ',Norder,Nofvert,symmv(m+1),Norder,
!    +      'symm. vars before calling the matrixscheme ',IFAIL)
!     write(6,*)'source ',(source(ivar),ivar=1,norder)
!
      CALL MatrixScheme(MatSplitVII,SYMMV(M+1),DSYMMV(M+1),TAUX(M+1),
     +                  BETA1,STIFC,NORDER,NORDER,NOFVERT,VCN,NDIM,
     +                  Jacobian,NEQMAX,RESIDUAL,SOURCE,IELEM,PICARD)
!
!     write(6,*)'residual ',(residual(ivar),ivar=1,2*norder)
!     write(6,*)ielem,' bp (2)'
!     CALL R8Mat_Print('G',' ',Nofvert,Nofeqn,dsymmv,Nofvert,
!    +      'residual in symm. vars. ',IFAIL)
!     CALL R8Mat_Print('G',' ',Norder,norder*nofvert,beta1,Norder,
!    +      'beta1 in EulerVI ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOrder,NOrder*nofvert,BETA1,
!    +      NOrder,'Beta in symm only last d+1 eqns',IFAIL)
C
C     at this stage the matrix DSYMM is as follows (same for TAUX)
C
C     dS(1) dp/ra(1) ...  dp/ra(4)
C     ....   du(1)   ...   du(4)
C            dv(1)   ...   dv(4) 
C     dS(4)  dw(1)   ...   dw(4)
C
C     the block involving symm. variables other than entropy
C     is transposed, same thing is done for the timestep
C
      CALL TRANS(DSYMMV(M+1),NORDER,NOFVERT,MN,MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)THEN
         WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
         CALL EXIT(IFAIL)
      ENDIF
      CALL TRANS(TAUX(M+1),NORDER,NOFVERT,MN,MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)THEN
         WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
         CALL EXIT(IFAIL)
      ENDIF
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFVERT,TAUX,NoFEQN,
!    +      'vettore DT dopo TRANS(2)',IFAIL)
C
C     so that now it looks like (same for TAUX):
C
C     dS(1) dp/ra(1)  du(1)  dv(1)  dw(1)
C     ....    ...      ...    ...    ...
C     ....    ...      ...    ...    ...
C           
C     dS(4) dp/ra(4)  du(4)  dv(4)  dw(4)
C
!     CALL TRANS(TAUX(1),NOFVERT,NOFEQN,NOFEQN*NOFVERT,MOVE,IWRK,IFAIL)
!     IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     compute the timestep: TAUX is transposed, see the table above for DSYMM
C
      DO IVERT = 1, NOFVERT ! loop over the vertices
         JADD = IVERT
         HELP = ZERO 
C
C     sum or maximum over the dofs
C
         DO IVAR = 1,NOFEQN ! loop over the NDIM+2 dofs
            IF( CHAR_TIMESTEPPING )THEN
               HELP = MAX(HELP,TAUX(JADD))
            ELSE
               HELP = HELP + TAUX(JADD)
            ENDIF
            JADD = JADD + NOFVERT
         ENDDO ! end loop over the dofs
         IADD = (IVERT-1)*NOFVAR
         DO IVAR = 1,NOFEQN
            JADD = IADD + IVAR
            TSTEP(JADD) = HELP
         ENDDO
      ENDDO
C
C
C Transform the nodal residual into conserved variables
C           note that DSYMMV is transposed during the MM product
C
      CALL DGEMM('No Transpose','Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dUdV,NOFEQN,DSYMMV(1),NOFVERT,ZERO,NODRES,
     +           NOFVAR)
C
C     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,dUdV,nofeqn,
C    +      'dUdV variables ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFvar,NOFvert,NODRES,nofvar,
!    +      'residual (before) ',IFAIL)
!     pause
C
C Insert the element stiffness matrix corresponding to the last
C d+1 eqns. into the element stiffness matrix of dimension (d+2)
C the offset equals the number of decoupled eqns.
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,BETA,
!    +      NOFEQN,'Betas ',IFAIL)
C
!     write(6,*)ielem,' bp (3)'
      IF (PICARD) THEN
         CALL MATINS(STIFEL,NOFVAR,STIFC,NORDER,NOFVERT,NOFVERT,1)
      ENDIF
C
      IF (LTIME.AND.DUALTS) 
     &CALL MATINS(BETA,NOFEQN,BETA1,NORDER,NOFVERT,1,1)
!     write(6,*)ielem,' bp (4)'
C
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,BETA,
!    +      NOFEQN,'Beta after insertion',IFAIL)
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
C
          ENDIF
C
C     transforms the residual into conserved variables
C
          CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV((FrstEq-1)*NOFEQN+1),
     +               NOFEQN,RESIDUAL,1,ONE,WKSP,1)
C
C     test the residual as computed by the "explicit" scheme
C     WKSP := dUdV(1) * \phi_{entropy} + dUdV(1,FrstEq) * \Phi
C
          CALL TEST(DivFlux,WKSP,TOLER,IELEM,NOFEQN)
      ENDIF ! on ICHECK
C
C --------------- Debugging code ends here ---------------
C
      IF (PICARD) THEN
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
      ENDIF ! on picard
C
      IF (LTIME) THEN
          DO IVERT = 1,NOFVERT
              IADD = (IVERT-1)*NOFVAR+1
              JADD = (IVERT-1)*NOFEQN*NOFEQN+1
              CALL PARM2CONS(VCZ(IADD),TEMPB(JADD),NOFEQN,NDIM)
          ENDDO
      ENDIF
C
!     CALL R8Mat_Print('G',' ',NOFVAR,nofvert,NODRES,
!    +      NOFVAR,'Nodal update in U before td',IFAIL)
!
      IF(LTIME)THEN 
c
c     compute transformation matrix from 
c
         CALL CONS2SYMM(ZAVG,BETA1,NOFEQN,NDIM)
!
!           CALL DGEMM('No Transpose','No Transpose',NOFEQN,
!    +               NOFEQN,NOFEQN,ONE,BETA1,NOFEQN,dUdV,NOFEQN,
!    2               ZERO,TEMPA(1),NOFEQN) 
!           IF(.NOT.UNITMAT(TEMPA,NOFEQN,NOFEQN,NOFEQN,1.D-12))THEN
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,BETA1,
!    +      NOFEQN,'dVdS ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,dUdV,
!    +      NOFEQN,'dUdV ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,TEMPB,
!    +      NOFEQN,'TEMPA ',IFAIL)
!           ENDIF 
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,BETA,
!    +      NOFEQN,'Beta in symm',IFAIL)
C
cold     DO IVERT = 1,NOFVERT ! si puo forse ridurre il nof di DGEMM !?!?!?!?!
cold        JADD = (IVERT-1)*NOFEQN*NOFEQN+1
cold        CALL DGEMM('No Transpose','No Transpose',NOFEQN,
cold +               NOFEQN,NOFEQN,ONE,dUdV,NOFEQN,BETA(JADD),NOFEQN,
cold 2               ZERO,TEMPA(JADD),NOFEQN) 
cold        CALL DGEMM('No Transpose','No Transpose',NOFEQN,
cold +               NOFEQN,NOFEQN,ONE,TEMPA(JADD),NOFEQN,BETA1,NOFEQN,
cold 2               ZERO,BETA(JADD),NOFEQN) 
cold     ENDDO
C
         CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,ONE,dUdV,NOFEQN,
     2               BETA,NOFEQN,ZERO,TEMPA,NOFEQN) 
         DO IVERT = 1,NOFVERT
            JADD = (IVERT-1)*NOFEQN*NOFEQN+1
            CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,ONE,TEMPA(JADD),NOFEQN,BETA1,NOFEQN,
     2               ZERO,BETA(JADD),NOFEQN) 
         ENDDO
!     write(6,*)ielem
!     pause
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +      NOFEQN,'Nodal update in U before the unstdy term',IFAIL)
C
         IF(DUALTS)
     &   CALL UNSTEADY4(TEMPB,BETA,VCZ,NOFVAR,NODRES,STIFEL,VOLUME,
     &                  NOFEQN,NDIM,NOFVERT,PICARD)
      ENDIF
C
!     write(6,*)'ielem = ',ielem
!     CALL R8Mat_Print('G',' ',Nofeqn,nofeqn*nofvert,beta,Nofeqn,
!    +      'beta(U) in EulerVI ',IFAIL)
!     CALL R8Mat_Print('G',' ',Nofvar,NOFvert*NTIMLEVS,vcz,Nofvar,
!    +      'Z variables ',IFAIL)
!     CALL R8Mat_Print('G',' ',NOFEQN,nofvert,NODRES,
!    +      NOFEQN,'Nodal update in U',IFAIL)
!     pause
C
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     Add the element stiffness matrix to the global stiffness matrix
C
C
C     transform the element stiffness matrix into conserved variables
C
      IF( NOFVAR .EQ. NOFEQN )THEN
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,MONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
      ELSE
          DO 13 JVERT = 1, NOFVERT
             DO 13 IVERT = 1, NOFVERT
                JADD = JADDR(IVERT,JVERT,NOFEQN)
                IADD = JADDR(IVERT,JVERT,NOFVAR)
                CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,MONE,dUdV,
     +               NOFEQN,STIFEL(IADD),NOFVAR,ZERO,TEMPA(JADD),NOFEQN)
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN,tempa(jadd),NOFEQN,
!    +      'TEMPA ',IFAIL)
!     pause
   13     CONTINUE
      ENDIF
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,TWO,dVdZ,NOFEQN,dZdU,
     +               NOFEQN,ZERO,TEMPB,NOFEQN)
!     CALL R8Mat_Print('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
!    +      'TEMPB ',IFAIL)
!     pause
C
          DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFEQN*NOFEQN + 1
          DO 35 IVERT = 1,NOFVERT
                  IADD = JADDR(IVERT,JVERT,NOFVAR)
                  KADD = JADDR(IVERT,JVERT,NOFEQN)
C
               CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +              NOFEQN,NOFEQN,ONE,TEMPA(KADD),NOFEQN,
     +              TEMPB(JADD),NOFEQN,ZERO,STIFEL(IADD),NOFVAR)
C
C     now STIFEL contains the convection stiffness matrix in
C     conserved variables
C
!     write(6,*)ivert,jvert
!     CALL R8Mat_Print('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
!    +      'C(i,j) ',IFAIL)
C
   35     CONTINUE
!     pause
C
C
      IF (ICHECK.EQ.0) RETURN
C
C --------------- Debugging code starts here ---------------
C
C     test the residual as computed by the "implicit" scheme
C     WKSP := dUdV(1) * \phi_{entropy} + dUdV(1,FrstEq) * \Phi
C
      CALL DINIT(NOFEQN,ZERO,WKSP,1)
      CALL DAXPY(NOFEQN,SCALRES,dUdV,1,WKSP,1)
      CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV((FrstEq-1)*NOFEQN+1),
     +           NOFEQN,RESIDUAL(NORDER+1),1,ONE,WKSP,1)
      CALL TEST(DivFlux,WKSP,TOLER,-IELEM,NOFEQN)
C
C --------------- Debugging code ends here ---------------
C
      RETURN

89999 FORMAT (5X,'Scalar residual in Element ',I6,' Wave # ',I1)
99999 FORMAT (5X,'Vector residual in Element ',I6,' EulerVII')

      END
