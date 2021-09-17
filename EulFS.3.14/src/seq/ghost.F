      SUBROUTINE GHOST(BCSUB,BC_TYPE,ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,
     +                 ZROE,STIFEL,NODRES,NEWRES,VCZ,VCDZ,VCB,
     +                 NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +                 IBNDFAC,NBFAC,PetscVecDT,PetscVecRHS,PetscMatA)
C
C     Purpose: far-field and inviscid wall boundary conditions;
C              this routine creates a dummy element for each
C              vertex of a boundary face
C
C
C     $Id: ghost.F,v 1.41 2020/04/23 09:51:45 abonfi Exp $
C
C
C
CCCCCC#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
C
      INCLUDE 'dofs.com'
      INCLUDE 'nloc.com'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
      INCLUDE 'periodic.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE 
C
C
      Vec PetscVecDT,PetscVecRHS
      Mat PetscMatA
C
C#define PRINT
C
C
C     .. Scalar Arguments ..
      INTEGER BC_TYPE,NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN,NGHOST,NPNOD
C     ..
C     .. Subroutine Arguments ..
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),XYZDOT(NDIM,*),VOL(NELEM),
     +                 ZROE(NOFVAR,*),
     3                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     +                 NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT),
     +                 VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT),
     5                 VCB(NDIM,NOFVERT)
      INTEGER IBNDFAC(3,NBFAC) 
C     ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLUME(MAXTIMLEVS+1),EPS
      INTEGER IVERT,IBC,ICLR,IELEM,IFACE,N,IV,JV,
     +        IFAIL,I,J,LOCA,LOCB,IVAR,JVAR
      LOGICAL COMPRESSIBLE,PLASMA,LWK
      INTEGER NERR,IOPT,NP
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
      EXTERNAL BCSUB
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION VCN(3*MAXNOFVERT),WORK(MAXNOFVERT*MAXNOFVAR),
     +TSTEP(MAXNOFVERT*MAXNOFVAR),DUMMY(MAX_NOFVAR_SQR*MAX_NOFVERT_SQR)
C
C     VCN(1:NDIM,1:NOFVERT)
C     VCZ(1:NDIM,1:NOFVERT)
C     TSTEP(1:NOFVERT-1)
C     NODRES(1:NOFVAR,1:NOFVERT-1)
C     STIFEL(1:NOFVAR,1:NOFVAR,1:NOFVERT-1)
C
C
      INTEGER ICN(VMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,PRESSC
      INTEGER ICYCL
      EXTERNAL DDOT,DNRM2,ICYCL,PRESSC
C     ..
C     .. External Subroutines ..
      EXTERNAL CELPTR,DAXPY,DCOPY,DGEMV,DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,REAL
C     ..
C     ..
C     .. Statement Functions ..
      INTEGER IDX
C     ..
C     ..
C     .. Statement Function definitions ..
C
C
      IDX(I,J,N) = (I-1)*N + J
C     ..
C
C
      COMPRESSIBLE = ABS(KAN) .EQ. 4
      PLASMA = ABS(KAN) .EQ. 3
      NP = NPOIN + NGHOST + NPNOD
C
      DO 100 IFACE = 1,NBFAC
C
         IELEM = IBNDFAC(1,IFACE)
         IVERT = IBNDFAC(2,IFACE)
         ICLR = IBNDFAC(3,IFACE)
C
         IBC = ICOLOR(ICLR,1) ! colours for the mean flow equations
C
C     Do NOT apply (inviscid) boundary conditions if:
C
C     IBC = 1 supersonic inlet
C     IBC = 3 supersonic outlet
C     IBC = 6 viscous wall
C     IBC = 7 prescribed velocity profile (maybe I should)
C
         IF (IBC.NE.BC_TYPE)GOTO 100
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
C
      INCLUDE 'ptest.f'
C
         IF(IVERT.NE.NOFVERT)THEN
            CALL DSWAP(NOFVAR,VCZ(1,IVERT),1,VCZ(1,NOFVERT),1)
            LOCA = IDX(IVERT,1,NDIM)
            LOCB = IDX(NOFVERT,1,NDIM)
            CALL DSWAP(NDIM,VCN(LOCA),1,VCN(LOCB),1)
            IF(LALE)CALL DSWAP(NDIM,VCB(1,IVERT),1,VCB(1,NOFVERT),1)
            LOCA = ICN(NOFVERT)
            ICN(NOFVERT) = ICN(IVERT)
            ICN(IVERT) = LOCA
            LWK = PFLAG(NOFVERT)
            PFLAG(NOFVERT) = PFLAG(IVERT)
            PFLAG(IVERT) = LWK
         ENDIF
C
         IF (BC_TYPE.EQ.BC_TYPE_FAR_FIELD)
     &   CALL DCOPY(NOFVAR,U_INFTY,1,VCZ(1,NOFVERT),1)
C
         CALL DINIT(NOFVAR*(NOFVERT-1),ZERO,TSTEP,1)
C
         CALL BCSUB(IELEM,VCN,VCZ,VCB,STIFEL,NODRES,TSTEP,
     +             NDIM,NOFVERT,NOFVAR)
C
C     Here we compute the forces on the solid wall
C
         IF(NDIM.EQ.2)THEN
             CALL PFORCESimpson(ICLR,NOFVERT,VCN,NDIM,VCZ,
     +                          NOFVAR,NOFVERT,PRESSC)
         ELSE
             CALL PFORCE(ICLR,NOFVERT,VCN,NDIM,VCZ,NOFVAR,
     +                   NOFVERT,PRESSC)
      ENDIF
C
         IF(.NOT.TIMEIMPL)GOTO 50
#if 0
      do ivert=1,nofvert-1
      write(6,*)ielem,ivert,ivert
      CALL R8Mat_Print('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
     +            Nofvar,'Picard matrix ',IFAIL)
      enddo
#endif
         IF(NEWTON)THEN
C
C     Compute Jacobian by finite differences
C
           DO 15 IVERT = 1, NOFVERT-1
            DO 15 IVAR = 1, NOFVAR
             CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,VCDZ,1)
             IF(COMPRESSIBLE)THEN
                CALL PARM_TO_CONS(VCZ (1,IVERT),NDIM,NOFVAR,1,.FALSE.,
     &          IFAIL)
                CALL PARM_TO_CONS(VCDZ(1,IVERT),NDIM,NOFVAR,1,.FALSE.,
     &          IFAIL)     
             ELSEIF(PLASMA)THEN
                CALL PARM_TO_CONS4Ar(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &          .FALSE.,IFAIL)
                CALL PARM_TO_CONS4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &          .FALSE.,IFAIL)
             ENDIF
             EPS = ROOT_MACHINE_EPS*MAX(ABS(VCZ(IVAR,IVERT)),1.d0)*
     +             SIGN(1.d0,VCZ(IVAR,IVERT))
             VCDZ(IVAR,IVERT) = VCZ(IVAR,IVERT) + EPS
             IF(COMPRESSIBLE)THEN
                CALL CONS_TO_PARM(VCZ (1,IVERT),NDIM,NOFVAR,1,.FALSE.,
     &          IFAIL)
                IF(IFAIL.NE.0)THEN
                    WRITE(ERRMSG(1:64),FMT=105)'GHOST_____',MY_PE,"(1)"
                    CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                ENDIF  
                CALL CONS_TO_PARM(VCDZ(1,IVERT),NDIM,NOFVAR,1,.FALSE.,
     &          IFAIL)
                IF(IFAIL.NE.0)THEN
                    WRITE(ERRMSG(1:64),FMT=105)'GHOST_____',MY_PE,"(2)"
                    CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                ENDIF   ! ifail
             ELSEIF(PLASMA)THEN
                CALL CONS_TO_PARM4Ar(VCZ (1,IVERT),NDIM,NOFVAR,1,
     &           .FALSE.,IFAIL)
                IF(IFAIL.NE.0)THEN
                    WRITE(ERRMSG(1:64),FMT=105)'GHOST_____',MY_PE,"(1)"
                    CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                ENDIF
                CALL CONS_TO_PARM4Ar(VCDZ(1,IVERT),NDIM,NOFVAR,1,
     &          .FALSE.,IFAIL)
                IF(IFAIL.NE.0)THEN
                    WRITE(ERRMSG(1:64),FMT=105)'GHOST_____',MY_PE,"(2)"
                    CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                ENDIF   ! ifail
              ENDIF ! compressible or plasma
C
#if 0
      CALL R8Mat_Print('General',' ',Nofvar,Nofvert,vcz(1,IVERT),
     +            Nofvar,'Z ',IFAIL)
      CALL R8Mat_Print('General',' ',Nofvar,Nofvert,vcdz(1,IVERT),
     +            Nofvar,'DZ ',IFAIL)
#endif
C
             CALL BCSUB(IELEM,VCN,VCDZ,VCB,DUMMY,NEWRES,WORK,
     +                   NDIM,NOFVERT,NOFVAR)
C
C     The residual sent to vertex IVERT only depends
C     on the node itself
C
             DO 15 JVAR = 1, NOFVAR
                 STIFEL(JVAR,IVAR,IVERT,1) =
     +           -(NEWRES(JVAR,IVERT)-NODRES(JVAR,IVERT))/EPS
C
   15      CONTINUE
      ENDIF ! Newton
C
      IF( PERIODIC_MESH .AND. ANNULAR )
     &CALL ROTATECIJ(STIFEL,NOFVAR,NOFVERT,NOFVERT-1,1)
C
      DO 25 IVERT = 1,NOFVERT-1 
C
C     This ugly copy is needed because MatSetValuesBlocked only
C     allows 1D arrays or 2D matrices, not 4D matrices such as STIFEL
C
         CALL DCOPY(NOFVAR*NOFVAR,STIFEL(1,1,IVERT,1),1,DUMMY,1)
#ifdef MPI
         CALL MatSetValuesBlockedLocal(PetscMatA,1,ICN(IVERT),1,
     +                                 ICN(IVERT),DUMMY,
     +                                 ADD_VALUES,IFAIL)
#else
         CALL MatSetValuesBlocked(PetscMatA,1,ICN(IVERT),1,ICN(IVERT),
     +                            DUMMY,ADD_VALUES,IFAIL)
#endif
#ifdef PRINT
      write(6,*)ielem,ivert,ivert
      CALL R8Mat_Print('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
     +            Nofvar,'Jacobian ',IFAIL)
#endif
C
   25 CONTINUE
C
C     insert rhs into PETSc vector
C
   50 CONTINUE
C
      INCLUDE 'rotaterhs2.f'
C
         CALL VecSetValuesBlocked(PetscVecDT,NOFVERT-1,ICN,TSTEP,
     +                      ADD_VALUES,IFAIL)
C
C     This ugly copy is needed because VecSetValuesBlocked only
C     allows 1D arrays
C
         CALL DCOPY(NOFVAR*(NOFVERT-1),NODRES,1,DUMMY,1)
         CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT-1,ICN,DUMMY,
     +                      ADD_VALUES,IFAIL)
C
  100 CONTINUE ! End loop over boundary faces
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
