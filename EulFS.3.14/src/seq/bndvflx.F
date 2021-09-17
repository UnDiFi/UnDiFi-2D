!> @param[in] NSModel is the SUBROUTINE used to discretize the viscous terms of the  mean-flow equations
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
!> @param[in] IBNDFAC
!> @param[in] NBFAC
!> @param[in] TurbulenceModel is the SUBROUTINE used to discretize the turbulence model
!> @param[in,out] PetscVecDT Petsc Vec where time-step is stored
!> @param[in,out] PetscVecRHS Petsc Vec where the nodal residual is stored
!> @param[in,out] PetscMatA Petsc Mat where the implicit matrix is stored
      SUBROUTINE BNDVFLX(NSModel,ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,
     1                   ZROE,ZTURB,STIFEL,NODRES,NEWRES,VCZ,VCDZ,VCB,
     2                   NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     3                   NTURB,IBNDFAC,NBFAC,TurbulenceModel,PetscVecDT,
     4                   PetscVecRHS,PetscMatA)
C
C     Purpose: driver for computing viscous fluxes on the boundaries
C
C     $Id: bndvflx.F,v 1.38 2020/04/23 09:51:45 abonfi Exp $
C
C
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      Mat PetscMatA
      Vec PetscVecRHS,PetscVecDT 
C
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'implicit.h'
C
C
      DOUBLE PRECISION DSTAK
      COMMON /CSTAK/DSTAK(1)
      INTEGER ISTAK(1)
      EQUIVALENCE(ISTAK(1),DSTAK(1))
C
C
      INCLUDE 'dofs.com'
      INCLUDE 'bnd.com'
      INCLUDE 'flags.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'nloc.com'
      INCLUDE 'periodic.com'
      INCLUDE 'stream.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
      EXTERNAL NSModel,TurbulenceModel,N_scheme,NSYS_scheme
C
C
      INTEGER NBFAC,NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NTURB
      INTEGER IBC,IELEM,IVERT,IFACE,LOC,ICLR,JV,
     +LOCB,IVAR,JVAR,JVERT,NP,N,IADDR,I,K,L,II,JJ,IFAIL
      LOGICAL COMPRESSIBLE,LFLG
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
      INTEGER IBNDFAC(3,NBFAC)
C
      INTEGER ICN(MAXNOFVERT),KCN(MAXNOFVERT)
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM),XYZDOT(NDIM,*),
     1                 ZROE(NOFVAR,*),ZTURB(*),
     2                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     3                 NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT),
     4                 VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT),
     5                 VCB(NDIM,NOFVERT)

      DOUBLE PRECISION VCN(3*MAXNOFVERT),VOLUME(MAXTIMLEVS+1),
     2                 VCT(MAXNOFVERT),
     +                 TSTEP(MAXNOFVERT*MAXNOFVAR),
     3                 DUMMY(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR)
C
C
C     LOC is a pointer to the position of the current face
C         in the array SKINF; it assumes that boundary faces
C         are ordered for increasing color
C
      DOUBLE PRECISION VISCL,VISCT
      INTEGER ICYCL
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION PRESSI,PRESSC
      EXTERNAL PRESSI,PRESSC
C
      DOUBLE PRECISION EPS
C
      DATA VISCL,VISCT/ONE,ZERO/


      COMPRESSIBLE = KAN .EQ. 4
      NP = NPOIN + NGHOST + NPNOD
      LFLG = .FALSE.
C
      LOC = LSKINF - 1
      LOCB = LHEAT - 1
C
      IF (TIMEIMPL)
     1    CALL DINIT(NOFVAR*NOFVAR*NOFVERT*NOFVERT,ZERO,STIFEL,1)
C
C     Loop over all boundary faces
C
      DO 100 IFACE = 1,NBFAC
C
         ICLR = IBNDFAC(3,IFACE)
         IBC = ICOLOR(ICLR,1) ! colours for the mean flow equations
C
C     Skip inviscid walls AND periodic boundaries
C
         IF (IBC.EQ.BC_TYPE_SLIP_FREE.OR.
     &       IBC.EQ.BC_TYPE_PERIODIC) GOTO 100
C
         IELEM = IBNDFAC(1,IFACE)
         IVERT = IBNDFAC(2,IFACE)
C
         CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
C
      INCLUDE 'ptest.f'
C
         IF(TURBULENT.AND.SEGREGATED)THEN
            DO 30 II = 1, NOFVERT
               VCT(II) = ZTURB(ICN(II)+1)
   30       CONTINUE
         ENDIF
C
C     here we do not use the trick of swapping the
C     vertex opposite the bnd. face into the nofvert position....
C
         IF ( .NOT. COMPRESSIBLE ) THEN
C
C     Compute pressure force
C
            IF (IBC.EQ.BC_TYPE_NO_SLIP) CALL PFORCE(ICLR,IVERT,VCN,
     &        NDIM,VCZ,NOFVAR,NOFVERT,PRESSI)
C
         ELSEIF ( COMPRESSIBLE ) THEN
C
C     Compute pressure force
C
            IF (IBC.EQ.BC_TYPE_NO_SLIP) THEN
               IF(NDIM.EQ.2)THEN
                  CALL PFORCESimpson(ICLR,IVERT,VCN,
     &                    NDIM,VCZ,NOFVAR,NOFVERT,PRESSC)
               ELSE
                  CALL PFORCE(ICLR,IVERT,VCN,
     &                    NDIM,VCZ,NOFVAR,NOFVERT,PRESSC)
               ENDIF
            ENDIF
C
         ELSE ! COMPRESSIBLE
C
              STOP 'Invalid KAN in Subroutine BNDVFLX'
C
         ENDIF ! COMPRESSIBLE
C
         IF (IBC.EQ.BC_TYPE_NO_SLIP) THEN
            LOC = LOC + 1
            LOCB = LOCB + 1
         ENDIF
C
C     reset LOCAL timestep and nodal residual to 0.d0
C
         CALL DINIT(NOFVAR*NOFVERT,ZERO,NODRES,1)
         CALL DINIT(NOFVAR*NOFVERT,ZERO,TSTEP,1)
C---should change below
         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +               VCZ,NOFVAR,VOLUME(1))
         IF(COMPRESSIBLE)CALL PARM2PRIM(NDIM,-IELEM)
C---should change above
caldo     CALL TurbulenceModel(VCT,NTURB,NOFVERT,arg,arg,COMPRESSIBLE,
caldo1                     viscl,visct)
C
cxxxx call TurbulenceModel(IELEM,VCN,VCT,NDIM,NOFVERT,NOFVAR,NTURB,
         CALL LAminar(IELEM,VCN,VCT,NDIM,NOFVERT,NOFVAR,NTURB,
     &                  NodRes,TSTEP,STIFEL,VOLUME(1),PICARD,
     &                  N_scheme,NSYS_scheme,DUMMY,COMPRESSIBLE,
     &                  VISCL,VISCT)
C
C
         CALL NSModel(ICLR,IELEM,IVERT,VCZ,NODRES,TSTEP,NOFVAR,NOFVERT,
     +                VCN,NDIM,VOLUME(1),STIFEL,VISCL,VISCT,DSTAK(LOC),
     +                DSTAK(LOCB),(IBC.EQ.BC_TYPE_NO_SLIP),PICARD)
C
C     Scatter and add nodal update and timestep to global arrays
C
      INCLUDE 'rotaterhs.f'
C
         CALL VecSetValuesBlocked(PetscVecDT,NOFVERT,ICN,TSTEP,
     +                      ADD_VALUES,IFAIL)
C
C     This ugly copy is needed because VecSetValuesBlocked only
C     allows 1D arrays
C
         CALL DCOPY(NOFVERT*NOFVAR,NODRES,1,DUMMY,1)
         CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT,ICN,DUMMY,
     +                      ADD_VALUES,IFAIL)
#if 0
         if(ibc.eq.BC_TYPE_NO_SLIP)then
            do i = 1,nofvert-1 
            II = ICYCL(IVERT+I,NOFVERT)
            do jj = 1,nofvert
               write(6,*)'ielem,ivert,ii,jj ',ielem,ivert,ii,jj
               CALL R8Mat_Print('General',' ',NOFVAR,NOFVert,
     +        stifd(1,1,ii,jj),nofvar,'Picard matrix ',ifail)
               enddo
            enddo
         endif
#endif
C
         IF(.NOT.TIMEIMPL)GOTO 100
         IF(PICARD)GOTO 45
C
C     loop over ALL vertices including the one opposite
C     the boundary face
C
         DO I = 1,NOFVERT 
            DO IVAR = 1, NOFVAR
               CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,VCDZ,1)
               IF(COMPRESSIBLE)THEN
                  CALL PARM_TO_CONS(VCZ(1,I),NDIM,NOFVAR,1,LFLG,IFAIL)
                  CALL PARM_TO_CONS(VCDZ(1,I),NDIM,NOFVAR,1,LFLG,IFAIL)
               ENDIF
               EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,I)),ONE)*
     +               SIGN(ONE,VCZ(IVAR,I))
               VCDZ(IVAR,I) = VCZ(IVAR,I) + EPS
               IF(COMPRESSIBLE)THEN
                  CALL CONS_TO_PARM(VCZ(1,I),NDIM,NOFVAR,1,LFLG,IFAIL)
                  IF(IFAIL.NE.0)THEN
                     WRITE(ERRMSG(1:64),FMT=105)'BNDVFLX___',MY_PE,"(1)"
                     CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                  ENDIF   ! ifail
                  CALL CONS_TO_PARM(VCDZ(1,I),NDIM,NOFVAR,1,LFLG,IFAIL)
                  IF(IFAIL.NE.0)THEN
                     WRITE(ERRMSG(1:64),FMT=105)'BNDVFLX___',MY_PE,"(2)"
                     CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                  ENDIF   ! ifail
               ENDIF
C
caldo
            CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     &                  VCDZ,NOFVAR,VOLUME(1))
            IF(COMPRESSIBLE)CALL PARM2PRIM(NDIM,-IELEM)
caldo    CALL TurbulenceModel(VCT,NTURB,NOFVERT,arg,arg,COMPRESSIBLE,
caldo1                     viscl,visct)
C
cxxx  call TurbulenceModel(IELEM,VCN,VCT,NDIM,NOFVERT,NOFVAR,NTURB,
            CALL Laminar(IELEM,VCN,VCT,NDIM,NOFVERT,NOFVAR,NTURB,
     &                  NodRes,TSTEP,STIFEL,VOLUME(1),PICARD,
     &                  N_scheme,NSYS_scheme,DUMMY,COMPRESSIBLE,
     &                  VISCL,VISCT)
C
caldo
C
C     Set NEWRES to 0.d0 since the "viscous" routines
C     ADD the viscous residual
C
            CALL DINIT(NOFVAR*NOFVERT,ZERO,NEWRES,1)
C
            CALL NSModel(ICLR,IELEM,IVERT,VCDZ,NEWRES,TSTEP,NOFVAR,
     &                NOFVERT,VCN,NDIM,VOLUME(1),DUMMY,VISCL,VISCT,
     &                DSTAK(LOC),DSTAK(LOCB),.FALSE.,.FALSE.)
C
#if 0
            CALL R8Mat_Print('General',' ',NOFVAR,NOFVERT,
     +            NEWRES,nofvar,'nodal residual ',ifail)
            pause
#endif
C
            DO JVAR = 1, NOFVAR
                DO JVERT = 1, NOFVERT-1
                    II = ICYCL(IVERT+JVERT,NOFVERT)
                    STIFEL(JVAR,IVAR,II,I) =
     +           -(NEWRES(JVAR,II)-NODRES(JVAR,II))/EPS
C
                ENDDO ! JVERT
            ENDDO ! JVAR
            ENDDO ! IVAR
         ENDDO ! I
C
C
   45 CONTINUE ! Jumps here when PICARD .EQV. .TRUE.
C     INCLUDE 'rotatecij.f'
         IF( PERIODIC_MESH .AND. ANNULAR )
     &CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NOFVERT-1,NOFVERT)
C
         N = (NOFVERT-1)*NOFVAR
         DO I = 1, NOFVERT-1 
            II = ICYCL(IVERT+I,NOFVERT)
            KCN(I) = ICN(II)
            DO JJ = 1, NOFVERT
               DO JVAR = 1, NOFVAR
                  DO IVAR = 1, NOFVAR
                     K = ( I-1)*NOFVAR+IVAR
                     L = (JJ-1)*NOFVAR+JVAR
                     IADDR = (L-1)*N+K
                     DUMMY(IADDR) = STIFEL(IVAR,JVAR,II,JJ)
                  ENDDO ! IVAR
               ENDDO ! JVAR
#if 0
               IF(IBC.EQ.IBC)THEN
                  WRITE(6,*)II,JJ
                  CALL R8Mat_Print('General',' ',NOFVAR,NOFVAR,
     +            STIFEL(1,1,II,JJ),NOFVAR,'STIFEL ',ifail)
               ENDIF
#endif
            ENDDO ! JJ
         ENDDO ! I
C
C     Add the element diffusion stiffness matrix to the global one
C
C     call MatSetValuesBlocked only once per boundary face
C
#ifdef MPI
         CALL MatSetValuesBlockedLocal(PetscMatA,NOFVERT-1,KCN,NOFVERT,
     +                     ICN,DUMMY,ADD_VALUES,IFAIL)
#else
         CALL MatSetValuesBlocked(PetscMatA,NOFVERT-1,KCN,NOFVERT,ICN,
     +                     DUMMY,ADD_VALUES,IFAIL)
#endif

  100 CONTINUE ! IFACE
C
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
