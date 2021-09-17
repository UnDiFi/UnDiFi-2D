!> \copydetails EulerIX()
      SUBROUTINE EulerVIIsa(IELEM,VCN,VCB,VCZ,NDIM,NOFVERT,NOFVAR,
     &                    NTURB,NODRES,TSTEP,STIFEL,VOLUME,PICARD,
     &                    ScalarScheme,MatrixScheme)
C
      IMPLICIT NONE 
C
C     Unsteady equations in symmetrizing variables ..
C
C     $Id: EulerVII.f,v 1.30 2002/09/14 09:10:40 abonfi Exp abonfi $
C
C
C
      INCLUDE 'paramt.h'
C
C     NEQMAX is the max. no. of equations (4 in 3D)
C            for the matrix scheme (solves for dp/ra,du,dv,dw)
C     MAXNOFEQN is the max. no. of mean flow equations (5 in 3D)
C
      INTEGER NEQMAX,LNNVV
      DOUBLE PRECISION TOLER
      PARAMETER (NEQMAX=4,TOLER=1.D-15)
      PARAMETER (LNNVV=NMAX*NMAX*VMAX*VMAX)
      INTEGER IWRK,FRSTEQ
      PARAMETER(IWRK=10,FRSTEQ=2)
      INTEGER MOVE(IWRK)
C
C
      INCLUDE 'constants'
      INCLUDE 'bnd.h'
      INCLUDE 'three'
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
      INTEGER NORDER,ifail,M,N,MN,N4,JCOL
      DOUBLE PRECISION FLUCT,SCALRES(2)
      LOGICAL LFLAG,PICARD
C
C
      DOUBLE PRECISION VCZ(*),VCN(*),VCB(*),VOLUME,
     +                 STIFEL(*),NODRES(*),TSTEP(*)
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
      DOUBLE PRECISION SYMMV(MAXNOFVAR*VMAX),DSYMMV(MAXNOFVAR*VMAX),
     +                 TAUX(MAXNOFEQN*VMAX),SOURCE(MAXNOFVAR)
      DOUBLE PRECISION Jacobian(NEQMAX,NEQMAX*3),
     +                 TEMPA((MAXNOFVAR*VMAX)**2),TMPV(MAXNOFVAR),
     2                 TEMPB((MAXNOFVAR*VMAX)**2)
      DOUBLE PRECISION PHI(NEQMAX),WKSP(MAXNOFVAR),RESIDUAL(2*NEQMAX),
     +                 STIFC(VMAX*VMAX*NEQMAX*NEQMAX)
C
C     NOFEQN (= DIM+2) is actual no. of mean flow equations
C     NORDER (= DIM+1) is actual no. of equations being solved
C                      with the system scheme
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
      INTEGER IADDR
      IADDR(IVERT,JVERT,N) = (((JVERT-1)*NOFVERT+IVERT-1)*N*N) + 1
C
      NORDER = NDIM + 1
C
C     solving also for turbulent viscosity
C
      NOFEQN = NDIM + 3
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
C
C     The element stiffness matrix is initialized to 0.d0
C
      IF (PICARD) THEN
          CALL DINIT(N4,ZERO,STIFEL,1)
          CALL DINIT((NORDER*NOFVERT)**2,ZERO,STIFC,1)
      ENDIF
C
C     set local residual and timestep to zero
C
      CALL DINIT(NOFVERT*NOFEQN,ZERO,DSYMMV,1)
      CALL DINIT(NOFVERT*NORDER,ZERO,TAUX,1)
C
C --------------- Debugging code starts here ---------------
C
cxxxx CALL X04CAF('G',' ',NOFVAR,NOFVERT,VCZ,NOFVAR,
cxxxx+      'VCZ ',IFAIL)
      IF (ICHECK.NE.0) THEN
C
C    Some initializations ....
C
          CALL DINIT(NORDER,ZERO,PHI,1)
          CALL DINIT(NOFEQN,ZERO,WKSP,1)
          CALL DINIT(NOFEQN,ZERO,TMPV,1)
          CALL CHECK2(IELEM,NDIM,NOFEQN)
cxxxx CALL X04CAF('G',' ',NOFVAR,1,DivFlux,NOFVAR,
cxxxx+      'Flux divergence ',IFAIL)
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C     The Jacobian Matrix of the subsystem is assembled and
C         the eigenvectors computed ..
C
      CALL Eigen_VIIsa(Jacobian,NEQMAX,dVdZ,dUdV,NDIM,NOFEQN)
C
      IF (ICHECK.NE.0) THEN
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
     +                  GRAD_CHAR(FrstEq,idim),1,ONE,PHI,1)

   12    CONTINUE
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
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
C     dS(1) dp/ra(1) du(1) dv(1) dw(1) dn_t(1)
C     ....    ....    ..    ..    ..
C     ....    ....    ..    ..    ..
C     dS(4) dp/ra(4) du(4) dv(4) dw(4) dn_t(4)
C
      CALL DGEMM('Transpose','Transpose',NOFVERT,NOFEQN,NOFEQN,
     +           ONE,VCZ,NOFVAR,dVdZ,NOFEQN,ZERO,SYMMV,NOFVERT)
C     CALL X04CAF('G',' ',Nofvert,NOFEQN,symmv,Nofvert,
C    +      'SYMM variables (transposed) dS, dp/ra, du, dn_t ',IFAIL)
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
      CALL ScalarScheme(IELEM,VCN,R_SPEED(1,IVAR),SCALRES(1),ZERO,
     +                  SYMMV(1),TAUX(1),DSYMMV(1),STIFC,NDIM,NOFVERT,
     +                  PICARD)
C
C     the first NOFVERT entries of TAUX now keep
C     the timestep computed by the Scalar Scheme
C
C     a little trick so that we can directly feed
C     TAUX to the matrix scheme:
C
      CALL TRANS(TAUX,NOFVERT,NORDER,NORDER*NOFVERT,MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C    Copy the convective jacobian (STIFC) into STIFEL and
C    reset it to 0. since it will be reused by the matrix scheme
C
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
          SCALRES(1) = SCALRES(1)/VOLUME
          CALL DAXPY(NOFEQN,SCALRES(1),dUdV((IVAR-1)*NOFEQN+1),1,WKSP,1)
C
          IF (DABS(FLUCT-SCALRES(1)).GT.TOLER) THEN
              WRITE (6,89999) IELEM,IVAR
              WRITE (6,FMT=*)FLUCT,SCALRES(1),ABS(FLUCT-SCALRES(1))
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
      MN = M*N
      CALL TRANS(SYMMV(M+1),M,N,MN,MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     ---------- Matrix scheme ----------
C
      SOURCE(2) = -GRAV(1)*VOLUME
      SOURCE(3) = -GRAV(2)*VOLUME
      IF(NDIM.EQ.3)SOURCE(4) = -GRAV(3)*VOLUME
C
C REM: calling with TAUX(1) will add new contributions
C      to those already computed when solving entropy
C
      CALL MatrixScheme(MatSplitVII,SYMMV(M+1),DSYMMV(M+1),TAUX(1),
     +                  STIFC,NORDER,NORDER,NOFVERT,VCN,NDIM,Jacobian,
     +                  NEQMAX,RESIDUAL,SOURCE,IELEM,PICARD)
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
C           PAUSE
C
          ENDIF
C
C     transforms the residual into conserved variables
C
          CALL DGEMV('N',NOFEQN,NORDER,ONE,dUdV((FrstEq-1)*NOFEQN+1),
     +               NOFEQN,RESIDUAL,1,ONE,WKSP,1)
C
      ENDIF
C
C --------------- Debugging code ends here ---------------
C
C
C     copy the timestep from TAUX into TSTEP (update)
C
      DO 14 IVERT = 1, NOFVERT
         IADD = (IVERT-1)*NOFVAR+1
         JADD = (IVERT-1)*NORDER+1
         TSTEP(IADD) = TSTEP(IADD) + TAUX(JADD)
   14 CONTINUE
C
C Insert the element stiffness matrix corresponding to the last
C d+1 eqns. into the element stiffness matrix of dimension (d+2)
C the offset equals the number of decoupled eqns.
C
      IF (PICARD) THEN
         CALL MATINS(STIFEL,NOFVAR,STIFC,NORDER,NOFVERT,NOFVERT,1)
C        it is enough to clear the first NOFVERT**2 entries
         CALL DINIT(NOFVERT*NOFVERT,ZERO,STIFC,1)
      ENDIF
C
C     at this stage the matrix DSYMM is as follows
C
C     dS(1) dp/ra(1) ...  dp/ra(4) dn_t(1)
C     ....   du(1)   ...   du(4)   ...
C            dv(1)   ...   dv(4)   ...
C     dS(4)  dw(1)   ...   dw(4)   dn_t(4)
C
C     the block involving symm. variables other than entropy
C     is transposed
C
      CALL TRANS(DSYMMV(M+1),N,M,MN,MOVE,IWRK,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,*)'TRANS HAS RETURNED IFAIL = ',IFAIL
C
C     so that now it looks like:
C
C     dS(1) dp/ra(1)  du(1)  dv(1)  dw(1) dn_t(1)
C     ....    ...      ...    ...    ...   ...
C     ....    ...      ...    ...    ...   ...
C           
C     dS(4) dp/ra(4)  du(4)  dv(4)  dw(4) dn_t(4)
C
C     ---------- Scalar scheme on the turbulent viscosity ----------
C
      IADD = (NOFVAR-1)*NOFVERT+1
C
C The turbulent viscosity eqn. is solved using an upwind scheme
C
      CALL DINIT(NOFVERT,ZERO,TAUX,1)
      CALL NS_scheme(IELEM,VCN,R_SPEED(1,1),SCALRES(2),ZERO,
     +                  SYMMV(IADD),TAUX(1),DSYMMV(IADD),STIFC,
     +                  NDIM,NOFVERT,PICARD)
      CALL DCOPY(NOFVERT,TAUX,1,TSTEP(NOFVAR),NOFVAR)
!     CALL X04CAF('G',' ',Nofvar,Nofvert,tstep,Nofvar,
!    +      'time-step ',IFAIL)
!     pause
C
      IF (PICARD) THEN
          IADD = NOFVAR*NOFVAR
          CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL(IADD),NOFVAR*NOFVAR)
      ENDIF
C
C Transform the nodal residual into conserved variables
C           note that DSYMMV is transposed during the MM product
C
      CALL DGEMM('No Transpose','Transpose',NOFEQN,NOFVERT,NOFEQN,
     +           ONE,dUdV,NOFEQN,DSYMMV(1),NOFVERT,ZERO,NODRES,
     +           NOFVAR)
C
C     write(6,*)ielem
C     CALL X04CAF('G',' ',Nofvar,Nofvert,nodres,Nofvar,
C    +      'residual ',IFAIL)
caldo do ivert = 1,nofvert
caldo CALL DAXPY(NOFEQN,-ONE/volume,NODRES((IVERT-1)*NOFVAR+1),1,TMPV,1)
caldo enddo
C
C --------------- Debugging code starts here ---------------
C
      IF (ICHECK.NE.0) THEN
C
C    Checks for the scalar residual ..
C
          IVAR = NOFVAR 
          FLUCT = DDOT(NDIM,R_SPEED(1,1),1,GRAD_CHAR(IVAR,1),LDW)
C
          SCALRES(2) = SCALRES(2)/VOLUME
          CALL DAXPY(NOFEQN,SCALRES(2),dUdV((IVAR-1)*NOFEQN+1),1,WKSP,1)
C
          IF (DABS(FLUCT-SCALRES(2)).GT.TOLER) THEN
              WRITE (6,89999) IELEM,IVAR
              WRITE (6,FMT=*)FLUCT,SCALRES(2),ABS(FLUCT-SCALRES(2))
          ENDIF
C
C     test the residual as computed by the "explicit" scheme
C
caldo     CALL TEST(DivFlux,TMPV,TOLER,IELEM,NOFVAR)
          CALL TEST(DivFlux,WKSP,TOLER,IELEM,NOFVAR)

      ENDIF
C
C --------------- Debugging code ends here ---------------
C
!     CALL X04CAF('G',' ',Nofvert,Nofvar,dsymmv,Nofvert,
!    +      'residual ',IFAIL)
!     do ivert=1,nofvert
!     do jvert=1,nofvert
!                 IADD = IADDR(IVERT,JVERT,NOFVAR)
!     CALL X04CAF('G',' ',NOFEQN,NOFEQN,STIFEL(iadd),nofvar,
!    +      'C(i,j) matrix ',IFAIL)
!     pause
!     enddo
!     enddo
C
C     --------------- If explicit, return now  ---------------
C
      IF (.NOT.PICARD) RETURN
C
C     compute the transformation matrices from
C     conserved to parameter variables in the vertices
C
          DO 33 IVERT = 1,NOFVERT
             IADD = (IVERT-1)*NOFVAR+1
             JADD = (IVERT-1)*NOFEQN*NOFEQN+1
             CALL Cons2ParmSA(VCZ(IADD),dZdU(JADD),NDIM,NOFEQN)
   33     CONTINUE
C
C     Add the element stiffness matrix to the global stiffness matrix
C
C
C     transform the element stiffness matrix into conserved variables
C
      IF( NOFVAR .EQ. NOFEQN )THEN
          CALL DGEMM('No Transpose','No Transpose',NOFVAR,
     +               NOFVAR*NOFVERT*NOFVERT,NOFVAR,-ONE,dUdV,
     +               NOFVAR,STIFEL,NOFVAR,ZERO,TEMPA,NOFVAR)
      ELSE
          DO 13 JVERT = 1, NOFVERT
          DO 13 IVERT = 1, NOFVERT
                  JADD = IADDR(IVERT,JVERT,NOFEQN)
                  IADD = IADDR(IVERT,JVERT,NOFVAR)
          CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN,NOFEQN,-ONE,dUdV,
     +               NOFEQN,STIFEL(IADD),NOFVAR,ZERO,TEMPA(JADD),NOFEQN)
!     write(6,*)ivert,jvert
!     CALL X04CAF('G',' ',NOFEQN,NOFEQN,tempa(jadd),NOFEQN,
!    +      'TEMPA ',IFAIL)
!     pause
   13     CONTINUE
      ENDIF
C
!         DO 37 JVERT = 1,NOFVERT
!         DO 37 IVERT = 1,NOFVERT
!                 IADD = (((JVERT-1)*NOFVERT+IVERT-1)*NOFEQN*NOFEQN) + 1
C
C
!  37     CONTINUE
C
      CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +               NOFEQN*NOFVERT,NOFEQN,TWO,dVdZ,NOFEQN,dZdU,
     +               NOFEQN,ZERO,TEMPB,NOFEQN)
!     CALL X04CAF('G',' ',NOFEQN,NOFEQN*nofvert,tempb,NOFEQN,
!    +      'TEMPB ',IFAIL)
!     pause
C
          DO 35 JVERT = 1,NOFVERT
                  JADD = (JVERT-1)*NOFEQN*NOFEQN + 1
          DO 35 IVERT = 1,NOFVERT
                  IADD = IADDR(IVERT,JVERT,NOFVAR)
                  KADD = IADDR(IVERT,JVERT,NOFEQN)
C
               CALL DGEMM('No Transpose','No Transpose',NOFEQN,
     +              NOFEQN,NOFEQN,ONE,TEMPA(KADD),NOFEQN,
     +              TEMPB(JADD),NOFEQN,ZERO,STIFEL(IADD),NOFVAR)
C
C     now STIFEL contains the convection stiffness matrix in
C     conserved variables
C
!     write(6,*)ivert,jvert
!     CALL X04CAF('G',' ',Nofvar,Nofvar,stifel(iadd),Nofvar,
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
C     WKSP := dUdV(1) * \phi_{entropy} + 
C             dUdV(n) * \phi_{viscosity} + 
C             dUdV(1,FrstEq) * \Phi
C
      CALL DINIT(NOFEQN,ZERO,WKSP,1)
      CALL DAXPY(NOFEQN,SCALRES(1),dUdV,1,WKSP,1)
      IADD = (NOFVAR-1)*NOFVAR+1
      CALL DAXPY(NOFEQN,SCALRES(2),dUdV(IADD),1,WKSP,1)
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
