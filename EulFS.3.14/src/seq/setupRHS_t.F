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
!> @param[in] Model is the SUBROUTINE used to discretize the turbulence transport equations
      SUBROUTINE SetupRHS_t(ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE,
     +                    ZTURB,STIFEL,NODRES,NEWRES,VCZ,VCDZ,VCB,NELEM,
     3                    NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NTURB,
     &                    PetscVecDT,PetscVecRHS,PetscMatA,
     +                    ScalarScheme,MatrixScheme,Model)
C
C     $Id: setupRHS_t.F,v 1.22 2021/01/29 08:53:23 abonfi Exp $
C
C
C#define BUDGET
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
C
      Vec PetscVecDT,PetscVecRHS
      Mat PetscMatA
C
      INCLUDE 'constants.h'
      INCLUDE 'three.com'
      INCLUDE 'flags.com'
      INCLUDE 'turb.com'
      INCLUDE 'newtoneps.com'
C
      INTEGER NDIM,NOFVERT,NTURB,NELEM,NPOIN,NGHOST,NPNOD,NOFVAR
C
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
C
C
      DOUBLE PRECISION VFACNOR(NDIM,*),XYZDOT(NDIM,*),VOL(NELEM),
     1 ZROE(NOFVAR,*), ZTURB(NTURB,*),
     +STIFEL(NTURB,NTURB,NOFVERT,NOFVERT), NODRES(NTURB,NOFVERT),
     +NEWRES(NTURB,NOFVERT), VCZ(NTURB,NOFVERT),VCDZ(NTURB,NOFVERT),
     4VCB(NDIM,NOFVERT)
C
C     NODRES(1:NTURB,1:NOFVERT) is used to accumulate the
C         nodal residual in conserved variables and scatter
C         it to the RHS PETSc vector
C     NEWRES(1:NTURB,1:NOFVERT) same as NODRES, but for the
C         perturbed dependent variable (only Newton iteration)
C     TSTEP(1:NOFVERT) is used to accumulate the timestep
C
C     ICN stores the vertices of the current element (0-based indexing)
C         and then scatter it to the PETSc vector PetscVecDT
C
      EXTERNAL ScalarScheme,MatrixScheme,Model
C
      INTEGER JVAR,IVAR,IVERT,JVERT,IELEM,NP,M,IOFF,IADDR
      INTEGER IFAIL,N4
      DOUBLE PRECISION EPS
C
      DOUBLE PRECISION TSTEP(MAXNOFVERT)
C
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION VCN(12),VOLUME(MAXTIMLEVS+1),DUMMY(400),RWORK(20)
      LOGICAL COMPRESSIBLE,MATRIX_ASSEMBLY
C
#ifdef BUDGET
      INCLUDE 'tdebug.com'
      INCLUDE 'nloc'
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
      INTEGER istkgt,loc,npoin
#endif
C
C     Some initializations ....
C
C     the "implicit" matrix is only computed by FDs
C
      MATRIX_ASSEMBLY = .FALSE.
C
      N4 = NTURB*NTURB*NOFVERT*NOFVERT
      NP = NPOIN + NGHOST + NPNOD 
C
      COMPRESSIBLE = ABS(KAN) .EQ. 4
C
C     set vorticity at trip point
C
      IF(TTFLAG.EQ.1)THEN
          CALL TSETUP(TTELEM,NELEM,ZROE,NDIM,NOFVERT,NOFVAR,
     +                (NPOIN+NGHOST+NPNOD),COMPRESSIBLE,RWORK)
          TST = RWORK(3)
          write(6,*)'vorticity at trip point is ',tst
      ENDIF
C
#ifdef BUDGET
      MATRIX_ASSEMBLY = .TRUE.
#endif
C
      DO 2000 IELEM = 1,NELEM
C
#ifdef BUDGET
      write(6,*)ielem
      call icopy(nofvert,icn,1,kcn,1)
#endif
C
C     The element stiffness matrix is initialized to 0.d0
C
          IF (TIMEIMPL) CALL DINIT(N4,ZERO,STIFEL,1)
C
C     set local residual and timestep to zero
C
          CALL DINIT(NOFVERT*NTURB,ZERO,NODRES,1)
          CALL DINIT(NOFVERT,ZERO,TSTEP,1)
C
C**********************************************
C
C     eddy viscosity
C
C**********************************************
C
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZTURB,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NTURB, NP, ICN, 
     3                VCZ, VCN, VCB, VOLUME)
C
!         CALL R8Mat_Print('General',' ',NDIM,NofvERT,VCB,
!    +                NDIM,
!    &                'grid velocity ',IFAIL)
C
          CALL TSETUP(IELEM,NELEM, ZROE, NDIM,NOFVERT, NOFVAR,
     +                (NPOIN+NGHOST+NPNOD),COMPRESSIBLE,RWORK)
C
C
C**********************************************
C
C     Compute nodal residual and timestep
C
C**********************************************
C
          CALL MODEL(IELEM,VCN,VCZ,VCB,NDIM,NOFVERT,
     1    NTURB,NODRES,TSTEP,STIFEL,VOLUME,MATRIX_ASSEMBLY,
     2    ScalarScheme,MatrixScheme,RWORK)
C
C
C     Insert nodal residual and timestep into PETSc vectors
C
          CALL VecSetValues(PetscVecDT,NOFVERT,ICN,TSTEP,ADD_VALUES,
     &                      IFAIL)
C
C     This ugly copy is needed because VecSetValuesBlocked only
C     allows 1D arrays
C
          CALL DCOPY(NTURB*NOFVERT,NODRES,1,DUMMY,1)
          IF(NTURB.EQ.1)THEN
              CALL VecSetValues(PetscVecRHS,NOFVERT,ICN,DUMMY,
     +                          ADD_VALUES,IFAIL)
          ELSE
              CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT,ICN,DUMMY,
     +                                 ADD_VALUES,IFAIL)
          ENDIF
C
          IF(.NOT.TIMEIMPL)GOTO 2000
C
          IF(MATRIX_ASSEMBLY)GOTO 1000
C
C     Compute Jacobian by finite differences
C
          DO IVERT = 1, NOFVERT
              DO IVAR = 1, NTURB
                  CALL DCOPY(NTURB*NOFVERT,VCZ,1,VCDZ,1)
C
                  EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),ONE)*
     +            SIGN(ONE,VCZ(IVAR,IVERT))
                  VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
C
                  CALL MODEL(IELEM,VCN,VCDZ,VCB,NDIM,NOFVERT,
     +            NTURB,NEWRES,TSTEP,DUMMY,VOLUME(1),.FALSE.,
     +            ScalarScheme,MatrixScheme,RWORK )
 
C
C     construct the jacobian entries
C
                  DO JVAR = 1, NTURB
                      DO JVERT = 1, NOFVERT
                          STIFEL(JVAR,IVAR,JVERT,IVERT) = -(NEWRES
     +                    (JVAR,JVERT)-NODRES(JVAR,JVERT))/EPS
C
                      ENDDO ! JVERT
                  ENDDO ! JVAR
             ENDDO ! IVAR
          ENDDO ! IVERT
 1000     CONTINUE
C
C     Insert element matrices into the global one
C
          IF(NTURB.EQ.1)THEN
C
C     Scalar problems, i.e. 1-eqn turbulence model ******************
C
#ifdef MPI
             CALL MatSetValuesLocal(PetscMatA,NOFVERT,ICN(1),NOFVERT,
     +                ICN(1),STIFEL(1,1,1,1),ADD_VALUES,IFAIL)
 
#else
             CALL MatSetValues(PetscMatA,NOFVERT,ICN,NOFVERT,ICN,
     +                STIFEL(1,1,1,1),ADD_VALUES,IFAIL)
#endif
          ELSE ! NTURB
C
C     2-eqn (or more) turbulence model
C
C     this part has NEVER been tested, because only SA is available
C
              DO JVERT = 1,NOFVERT
                 m = (jvert-1)*NOFVERT*NTURB*NTURB
                 DO IVERT = 1,NOFVERT
                    ioff = m + (IVERT-1)*NTURB
                    do jvar = 1, NTURB
                       iaddr = ioff + (jvar-1)*NTURB*NOFVERT
                       do ivar = 1, NTURB
                          DUMMY(iaddr+ivar) = 
     &                         STIFEL(ivar,jvar,IVERT,JVERT)
                       enddo ! i
                    enddo ! j
#ifdef PRINT_MATRICES
                     write(6,*)ivert,jvert
                     CALL R8Mat_Print('General',' ',NTURB,NTURB,
     +               stifel(1,1,ivert,jvert), NTURB,
     &               'Jacobian Matrix in setupRHS_t (1) ',IFAIL)
#endif
                 enddo ! IVERT
              enddo ! JVERT
C
C Insert all at once
C
#ifdef MPI
              CALL MatSetValuesBlockedLocal(PetscMatA,NOFVERT,
     +                        ICN,NOFVERT,ICN,DUMMY,
     +                         ADD_VALUES, IFAIL)
#else
              CALL MatSetValuesBlocked(PetscMatA,NOFVERT,ICN,
     +                               NOFVERT,ICN,DUMMY,ADD_VALUES,
     +                            IFAIL)
#endif
         ENDIF ! NTURB
C
 2000 CONTINUE
caldo
#ifdef BUDGET
      call vecgetsize(dt,npoin,ifail)
      write(6,*)ifail,npoin,lvol,lcelnod,nofvert,nelem
      loc = istkgt(npoin,4)
      call median(dstak(loc),npoin,dstak(lvol),istak(lcelnod),
     +            nofvert,nelem)
      write(6,*)'beyond median'
C
      call dcopy(npoin,tmpdif1,1,tmpdiff,1)
      call daxpy(npoin,ONE,tmpdif2,1,tmpdiff,1)
C
      call dcopy(npoin,tmpdiff,1,tmpsum,1)
      call daxpy(npoin,ONE,tmpsou1,1,tmpsum,1)
      call daxpy(npoin,ONE,tmpsou2,1,tmpsum,1)
C
C     SUM = PROD - DESTR + DIFF
C
      call scala(npoin,tmpconv,dstak(loc))
      call scala(npoin,tmpdif1,dstak(loc))
      call scala(npoin,tmpdif2,dstak(loc))
      call scala(npoin,tmpsou1,dstak(loc))
      call scala(npoin,tmpsou2,dstak(loc))
      call scala(npoin,tmpdiff,dstak(loc))
      call scala(npoin,tmpsum,dstak(loc))
C
      call scala(npoin,tmpsou3,dstak(loc))
C
      call istkrl(1)
      call SOLZNE('file030.dat',tmpconv,1,NPOIN,'w')
      call SOLZNE('file031.dat',tmpdif1,1,NPOIN,'w')
      call SOLZNE('file032.dat',tmpdif2,1,NPOIN,'w')
      call SOLZNE('file033.dat',tmpsou1,1,NPOIN,'w')
      call SOLZNE('file034.dat',tmpsou2,1,NPOIN,'w')
      call SOLZNE('file036.dat',tmpdiff,1,NPOIN,'w')
      call SOLZNE('file037.dat',tmpsum,1,NPOIN,'w')
      call SOLZNE('file035.dat',tmpsou3,1,NPOIN,'w')
      stop
#endif
      RETURN
      END

C
      SUBROUTINE SCALA(NPOIN,V,DT)


C     .. Scalar Arguments ..
      INTEGER NPOIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 1 I = 1,NPOIN
          V(I) = V(I)/DT(I)
    1 CONTINUE
      RETURN

      END

