!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] VFACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] XYZDOT the cartesian components of the grid velocity
!> @param[in] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] ZROE Nodal values of the dependent variable;  compressible equations: \f$ Z = \sqrt{\rho} \left( 1, H, \mathbf{u} \right) \f$; incompressible equations: \f$ Z = \left( p, \mathbf{u} \right) \f$
!> @param[in] ZTURB Nodal values of the turbulent variable(s); for both the compressible and incompressible equations: \f$ \tilde{\nu} \f$
!> @param[in,out] STIFEL elemental stiffness matrix
!> @param[in,out] NODRES 
!> @param[in,out] NEWRES 
!> @param[out] VCZ returns the \c NOFVAR dofs of the NOFVERT vertices of cell \c IELEM; VCZ(1:NOFVAR,*,1) stores the values at time level \c n+1,k ; \c VCZ(1:NOFVAR,*,2) those at tme level n and \c VCZ(1:NOFVAR,*,3) those at time level \c n-1
!> @param[out] VCDZ returns the \c NOFVAR dofs of the NOFVERT vertices of cell \c IELEM; VCZ(1:NOFVAR,*,1) stores the values at time level \c n+1,k ; \c VCZ(1:NOFVAR,*,2) those at tme level n and \c VCZ(1:NOFVAR,*,3) those at time level \c n-1
!> @param[out] VCB returns the \c NDIM Cartesian components of the nodal grid velocities at time \c n+1/2 of the NOFVERT vertices of cell \c IELEM; only if \c LALE is \c .TRUE.
!> @param[in] NELEM is the no. of processor owned elements (triangles/tetrahedra); global number of elements in the uni-processor case
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in] NGHOST nof ghost nodes in the mesh
!> @param[in] NPNOD nof periodic nodes in the mesh
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT number of vertices per element (=NDIM+1, since only triangles or tetrahedra are allowed)
!> @param[in] NOFVAR number of variables (degrees of freedom) in each meshpoint
!> @param[in] NTURB nof turbulent dofs =1 for SA
!> @param[in,out] PetscVecDT Petsc Vec where time-step is stored
!> @param[in,out] PetscVecRHS Petsc Vec where the nodal residual is stored
!> @param[in,out] PetscMatA Petsc Mat where the implicit matrix is stored
!> @param[in] ScalarScheme is the SUBROUTINE used to discretize scalar equation
!> @param[in] MatrixScheme is the SUBROUTINE used to discretize hyperbolic systems of equations
!> @param[in] TModelScheme is the SUBROUTINE used to discretize the turbulence transport equations
!> @param[in] EulerModel is the SUBROUTINE used to discretize the inviscid mean-flow equations
!> @param[in] NSModel is the SUBROUTINE used to discretize the viscous terms of the  mean-flow equations
!> @param[in] TurbulenceModel is the SUBROUTINE used to discretize the turbulence model
      SUBROUTINE SetupRHS(ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE,ZTURB,
     +                    STIFEL,NODRES,NEWRES,VCZ,VCDZ,VCB,
     3                    NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     4                    NTURB,PetscVecDT,PetscVecRHS,PetscMatA,
     +                    ScalarScheme,MatrixScheme,TModelScheme,
     6                    EulerModel,NSModel,TurbulenceModel)
C
C     $Id: setupRHS.F,v 1.35 2021/02/02 10:23:18 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.h'
      INCLUDE 'constants.h'
C
      Vec PetscVecDT,PetscVecRHS
      Mat PetscMatA
C
C
C#define TEST_JACOBIAN
C#define PRINT_MATRICES
C
C     On Aug. 21, 2002 I have tested the diffusive terms of 
C     the incompressible NS eqns.
C
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'periodic.com'
      INCLUDE 'nloc.com'
      INCLUDE 'sa.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'time.com'
C
C
      INTEGER NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NELEM,NTURB
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
C
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM), ZROE(NOFVAR,*),
     1                 ZTURB(*),XYZDOT(NDIM,*),
     +STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT), NODRES(NOFVAR,NOFVERT),
     +NEWRES(NOFVAR,NOFVERT), VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT),
     +VCB(NDIM,NOFVERT)
 
C
C     FACNOR(1:NDIM,1:NFACE)  cartesian components
C                             of the NFACE faces/edges
C     NODRES(1:NOFVAR,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the PetscVecRHS PETSc vector
C     NEWRES(1:NOFVAR,1:NOFVERT) same as NODRES, but for the
C         perturbed dependent variable (only Newton iteration)
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C         and then scatter it to the PETSc vector PetscVecDT
C
C     ICN stores the vertices of the current element (0-based indexing)
C
C
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      EXTERNAL ScalarScheme,MatrixScheme,
     +         EulerModel,NSModel,TurbulenceModel,TModelScheme
C
C
      INTEGER JVAR,IVAR,IVERT,JVERT,JV,IELEM,NP,IADDR,IOFF,IFAIL,N4
      INTEGER I,J,M
      DOUBLE PRECISION EPS
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C     ..
C     .. Local Arrays ..
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION TSTEP(MAXNOFVAR*MAXNOFVERT)
      DOUBLE PRECISION VCN(3*MAXNOFVERT),VOLUME(MAXTIMLEVS+1),
     &DUMMY(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR),VISCT(MAXNOFVERT)
      LOGICAL COMPRESSIBLE,VISCOUS,PLASMA
C
#ifdef TEST_JACOBIAN
      DOUBLE PRECISION ERRMAT((MAXNOFVAR*MAXNOFVERT)**2)
      DOUBLE PRECISION STIFD((MAXNOFVAR*MAXNOFVERT)**2)
      DOUBLE PRECISION S
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
      IADDR(I,J) = ((J-1)*NOFVERT+(I-1))*(NOFVAR**2)+1
#endif
C
C     Some initializations ....
C
      N4 = NOFVAR*NOFVAR*NOFVERT*NOFVERT
      CALL DINIT(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR,ZERO,DUMMY,1)
C
C
      COMPRESSIBLE = ABS(KAN) .EQ. 4
      PLASMA = ABS(KAN) .EQ. 3
      VISCOUS = (KAN .GT. 0)
      NP = NPOIN + NGHOST + NPNOD
C
C
      DO 2000 IELEM = 1,NELEM
C
C     The element stiffness matrix is initialized to 0.d0
C
          IF (TIMEIMPL) CALL DINIT(N4,ZERO,STIFEL,1)
C
C     set local residual and timestep to zero
C
          CALL DINIT(NOFVERT*NOFVAR,ZERO,NODRES,1)
          CALL DINIT(NOFVERT*NOFVAR,ZERO,TSTEP,1)
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
!         write(6,*)ielem,(icn(j)+1,j=1,nofvert)
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert*NTIMLEVS,VCZ,
!    +    Nofvar,'Z array just after celptr ',IFAIL)
C
C     in the segregated approach, VISCT is not altered,
C     so we can compute it only once.
C
          IF(TURBULENT.AND.SEGREGATED)THEN
              DO 10 IVERT = 1, NOFVERT
                  VISCT(IVERT) = ZTURB(ICN(IVERT)+1)
   10         CONTINUE
          ELSEIF(TURBULENT.AND.COUPLED)THEN
C
C     here we compute wall and trip point distance
C     as required by the Spalart-Allmaras model
C     this is NOT very elegant (ugly!!)
C     also requires DSTAK
C
              TD  = ZERO
              TTD = ZERO
              DO 12 IVERT = 1, NOFVERT
                  TD  =  TD + DSTAK(LTD+ICN(IVERT))
                  TTD = TTD + DSTAK(LTTD+ICN(IVERT))
                  visct(ivert) = ZERO
   12         CONTINUE
              TD  = TD / NOFVERT
              TTD = TTD / NOFVERT
          ENDIF
C
C     test for periodic i-nodes
C
      INCLUDE 'ptest.f'
C
C     Compute nodal residual and timestep
C
#ifdef TEST_JACOBIAN
          CALL NONAME(IELEM,VCN,VCB,VCZ,VISCT,NDIM,NOFVERT,NOFVAR,
     1    NTURB,VOLUME(1),.TRUE.,VISCOUS,COMPRESSIBLE,
     2    EulerModel,TurbulenceModel,NSModel,
     3    ScalarScheme,MatrixScheme,TModelScheme,
     4    NODRES,TSTEP,STIFD)
#else
!         write(6,*)ielem
!         write(6,*)(vcn(i),i=1,nofvert*ndim)
!         do i = 1,nofvert
!         write(6,*) (vcz(j,i),j=1,nofvar)
!         enddo
!         write(6,*)'visct = ',(visct(i),i=1,nofvert)
!         write(6,*) ndim,nofvert,nofvar,
!    &    NTURB,VOLUME(1),PICARD,VISCOUS,COMPRESSIBLE
caldo
          CALL NONAME(IELEM,VCN,VCB,VCZ,VISCT,NDIM,NOFVERT,NOFVAR,
     1                NTURB,VOLUME(1),PICARD,VISCOUS,COMPRESSIBLE,
     2                EulerModel,TurbulenceModel,NSModel,
     3                ScalarScheme,MatrixScheme,TModelScheme,
     4                NODRES,TSTEP,STIFEL)
C
!         write(6,*)'Beyond noname outside Newton loop ',ielem
!         write(6,*)' ielem = ',ielem,' ivert = ',ivert,' ivar = ',ivar
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert,Nodres,
!    +    Nofvar,'REsidual Matrix before Newton',IFAIL)
#endif
C
          IF(.NOT.TIMEIMPL)GOTO 1500
C
          IF(NEWTON)THEN
C
C     Compute Jacobian by finite differences
C
             DO IVERT = 1, NOFVERT
                DO IVAR = 1, NOFVAR
                   CALL DCOPY(NOFVAR*NOFVERT*NTIMLEVS,VCZ,1,VCDZ,1)
C
C     For compressible flows, it is necessary to transform to
C     conserved variables, perturb and then transform back to
C     parameter vector
C
                  IF(COMPRESSIBLE)THEN
                    CALL PARM_TO_CONS(VCZ (1,IVERT),NDIM,NOFVAR,1,LTSQR,
     &                                IFAIL)
                    CALL PARM_TO_CONS(VCDZ(1,IVERT),NDIM,NOFVAR,1,LTSQR,
     &                                IFAIL)
                  ELSEIF(PLASMA)THEN
                    CALL PARM_TO_CONS4Ar(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &                                .FALSE.,IFAIL)
                    CALL PARM_TO_CONS4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                .FALSE.,IFAIL)
                  ENDIF
C
                  EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),ONE)*
     +            SIGN(ONE,VCZ(IVAR,IVERT))
                  VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
C
                  IF(COMPRESSIBLE)THEN
                    CALL CONS_TO_PARM(VCZ (1,IVERT),NDIM,NOFVAR,1,LTSQR,
     &                                IFAIL)
                  ELSEIF(PLASMA)THEN
                    CALL CONS_TO_PARM4Ar(VCZ(1,IVERT),NDIM,NOFVAR,1,
     &                                .FALSE.,IFAIL)
                  ELSE ! incompressible
                      IFAIL = 0
                  ENDIF ! compressible
                  IF(IFAIL.NE.0)THEN
                       WRITE(ERRMSG(1:64),FMT=105)'setupRHS__',MY_PE,
     &                                            "(1)"
                       CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                  ENDIF   ! ifail
                  IF(COMPRESSIBLE)THEN
                    CALL CONS_TO_PARM(VCDZ(1,IVERT),NDIM,NOFVAR,1,LTSQR,
     &                                IFAIL)
                  ELSEIF(PLASMA)THEN
                    CALL CONS_TO_PARM4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &                                .FALSE.,IFAIL)
                  ELSE ! incompressible
                      IFAIL = 0
                  ENDIF ! compressible
                  IF(IFAIL.NE.0)THEN
                     WRITE(ERRMSG(1:64),FMT=105)'setupRHS__',MY_PE,
     &                                            "(2)"
                     CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                  ENDIF   ! ifail
C
C         we set the "perturbed" residual to 0.d0
C         since this is not any more done within
C         the EulerModel routines. Note, however,
C         that most scheme would do this (to be checked)
C
                  CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,NEWRES,1)
C
!         write(6,*)' before noname in Newton loop ',ielem
!         write(6,*)ielem
!         write(6,*)(vcn(i),i=1,nofvert*ndim)
!         do i = 1,nofvert
!         write(6,*) (vcdz(j,i),j=1,nofvar)
!         enddo
!         write(6,*)'visct = ',(visct(i),i=1,nofvert)
!         write(6,*) ndim,nofvert,nofvar,
!    &    NTURB,VOLUME(1),PICARD,VISCOUS,COMPRESSIBLE
!         call dcopy(nofvar*nofvert,VCDZ,1,JUMMY,1)
!         iflg = 13
                  CALL NONAME(IELEM,VCN,VCB,VCDZ,VISCT,NDIM,NOFVERT,
     1                        NOFVAR,NTURB,VOLUME(1),.FALSE.,VISCOUS,
     2                        COMPRESSIBLE,EulerModel,TurbulenceModel,
     3                        NSModel,ScalarScheme,MatrixScheme,
     3                        TModelScheme,NEWRES,DUMMY,DUMMY)
C
!         write(6,*)' ielem = ',ielem,' ivert = ',ivert,' ivar = ',ivar
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert*NTIMLEVS,VCZ,
!    +    Nofvar,'Z array ',IFAIL)
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert*NTIMLEVS,VCDZ,
!    +    Nofvar,'Perturbed Z array ',IFAIL)
!         CALL R8Mat_Print('General',' ',Nofvar,Nofvert,NEWRES,
!    +    Nofvar,'REsidual Matrix ',IFAIL)
!         pause
C
C     construct the jacobian entries using one-sided FD
C
                  DO JVAR = 1, NOFVAR
                      DO JVERT = 1, NOFVERT
                          STIFEL(JVAR,IVAR,JVERT,IVERT) = -(NEWRES
     +                    (JVAR,JVERT)-NODRES(JVAR,JVERT))/EPS
C
                      ENDDO ! JVERT
                   ENDDO ! JVAR
                ENDDO ! loop over IVAR
             ENDDO ! loop over IVERT
          ENDIF ! Newton
C
C   Insert element matrices into the global one
!
!   STIFEL è una matrice quadridimensionale
!   double precision STIFEL (bs,bs,nvt,nvt) 
!   dove (nell'esempio qui sotto)
!   bs = NOFVAR = 2, nvt = NOFVERT = 3
!   -------------------------------------
!   |  1.   3. |   13.  15.|   25.   27.|
!   |  2.   4. |   14.  16.|   26.   28.| 
!   -------------------------------------
!   |  5.   7. |   17.  19.|   29.   31.| 
!   |  6.   8. |   18.  20.|   30.   32.| 
!   -------------------------------------
!   |  9.   11.|   21.  23.|   33.   35.| 
!   |  10.  12.|   22.  24.|   34.   36.| 
!   -------------------------------------
!   Vista come un vettore 1d, gli elementi
!   di questo sono nell'ordine degli interi,
!   perché riempio la matrice con:
!   do jvt = 1, nvt
!      do ivt = 1, nvt
!         do j = 1, bs
!            do i = 1, bs
!               C(i,j,ivt,jvt) = ...
!            enddo
!         enddo
!      enddo
!   enddo
!   Sull' ordine degli elementi di C non posso intervenire
!   senza stravolgere troppe cose, pero':
!   Fino alla versione 1.31 inserivo uno alla volta i blocchi nella matrice
!   PETsc, cioe' chiamavo nvt**2 volte MatSetValuesBlocked()
!
!   Nella 1.32 faccio un'unica chiamata a MatSetValuesBlocked(),
!   ma devo passare a Petsc una matrice 2d
!   che "vede" C come un'unica matrice, cioe' 
!   A(bs*nvt,bs*nvt) =
!   -----------------------------------
!   |  1.   3.    13.  15.   25.   27.|
!   |  2.   4.    14.  16.   26.   28.| 
!   |  5.   7.    17.  19.   29.   31.| 
!   |  6.   8.    18.  20.   30.   32.| 
!   |  9.   11.   21.  23.   33.   35.| 
!   |  10.  12.   22.  24.   34.   36.| 
!   -----------------------------------
!   oppure, passandogli il vettore 1d
!   y = 1,2,5,6,9,10,3, etc.
!   
!   C'è un modo più clever di fare questa cosa,
!   rispetto alla copia che faccio adesso?
C
          IF(NOFVAR.EQ.1)THEN
C
C     Scalar problems ******************
C
#ifdef MPI
             CALL MatSetValuesLocal(PetscMatA,NOFVERT,ICN(1),NOFVERT,
     +                ICN(1),STIFEL(1,1,1,1),ADD_VALUES,
     +                IFAIL)
#else
             CALL MatSetValues(PetscMatA,NOFVERT,ICN(1),NOFVERT,
     +                ICN(1),STIFEL(1,1,1,1),ADD_VALUES,
     +                IFAIL)
#endif
          ELSE ! NOFVAR
C
              IF( PERIODIC_MESH .AND. ANNULAR )
     &        CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NOFVERT,NOFVERT)
C
C     Euler/NS problems ******************
C
              DO 37 IVERT = 1,NOFVERT
                  DO 37 JVERT = 1,NOFVERT
                      CALL DCOPY(NOFVAR*NOFVAR,STIFEL(1,1,IVERT,JVERT),
     &                           1,DUMMY,1)
#ifdef MPI
                      CALL MatSetValuesBlockedLocal(PetscMatA,1,
     +                ICN(IVERT),1,ICN(JVERT),DUMMY,ADD_VALUES, IFAIL)
#else
                      CALL MatSetValuesBlocked(PetscMatA,1,ICN(IVERT),1,
     +                ICN(JVERT),DUMMY,ADD_VALUES,IFAIL)
#endif
#ifdef PRINT_MATRICES
                      ERRMSG(1:32) = 'Jacobian Matrix in setupRHS (1) '
                      write(ERRMSG(33:59),FMT=120)IELEM,IVERT,JVERT 
                      CALL R8Mat_Print('General',' ',Nofvar,Nofvar,
     +                stifel(1,1,ivert,jvert),Nofvar,ERRMSG(1:59),
     &                IFAIL)
#endif
   37         CONTINUE
C
          ENDIF ! Check on NOFVAR
C
#ifdef TEST_JACOBIAN
C
              CALL DCOPY(N4,STIFD,1,ERRMAT,1)
              CALL DAXPY(N4,MONE,STIFEL,1,ERRMAT,1)
              DO 39 IVERT = 1,NOFVERT
                  DO 39 JVERT = 1,NOFVERT
                     I= IADDR(IVERT,JVERT)
                     S = DNRM2(NOFVAR*NOFVAR,ERRMAT(I),1)
                     IF(ABS(S).GT.1.E-5)THEN
                        write(6,*)ielem,icn(ivert),icn(jvert),s
                 CALL R8Mat_Print('General',' ',Nofvar,Nofvar,stifel(1,
     +           1,ivert,jvert), Nofvar,
     &           'Jacobian Matrix in setupRHS (1) ',IFAIL)
                 CALL R8Mat_Print('General',' ',Nofvar,Nofvar,STIFD(i),
     +           Nofvar,'Picard Matrix ',IFAIL)
                 pause
                     ENDIF
   39         CONTINUE
#endif
C
 1500 CONTINUE
C
C     Insert nodal residual and timestep into PETSc vectors
C
       INCLUDE 'rotaterhs.f'
C
caldo     CALL R8Mat_Print('General',' ',Nofvar,Nofvert,NODRES,
caldo+    Nofvar,'REsidual Matrix ',IFAIL)
caldo     pause
C
         
         CALL DCOPY(NOFVAR*NOFVERT,NODRES,1,DUMMY,1)
         IF(NOFVAR.EQ.1)THEN
              CALL VecSetValues(PetscVecDT,NOFVERT,ICN,TSTEP,
     +        ADD_VALUES,IFAIL)
              CALL VecSetValues(PetscVecRHS,NOFVERT,ICN,DUMMY,
     +        ADD_VALUES,IFAIL)
         ELSE
              CALL VecSetValuesBlocked(PetscVecDT,NOFVERT,ICN,TSTEP,
     +        ADD_VALUES,IFAIL)
              CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT,ICN,DUMMY,
     +        ADD_VALUES,IFAIL)
         ENDIF
 2000 CONTINUE ! end loop over elements
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
  120 FORMAT(1X,'ie = ',I6,' iv = ',I1,' jv = ',I1)
      END
