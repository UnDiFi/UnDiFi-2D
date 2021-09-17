      SUBROUTINE GHOST4(BC5,ICELNOD,ICELFAC,VFACNOR,VOL,ZROE,
     +                 STIFEL,NODRES,NEWRES,VCZ,VCDZ,
     +                 NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +                 IBNDFAC,NBFAC,PetscVecDT,PetscVecRHS,PetscMatA)
C
C     Purpose: inviscid wall boundary conditions;
C              this routine creates a dummy element for each
C              vertex of a boundary face
C
C
C     $Id: ghost.F,v 1.28 2009/06/11 13:10:58 abonfi Exp $
C
C
      IMPLICIT NONE
C
C
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
C
      Vec PetscVecDT,PetscVecRHS
      Mat PetscMatA
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants'
      INCLUDE 'implicit.h'
C
      INCLUDE 'nloc'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
      INCLUDE 'periodic.com'
      INCLUDE 'newtoneps.com'
C
C#define PRINT
C
C
      INCLUDE 'three'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE 
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN,NGHOST,NPNOD
C     ..
C     .. Subroutine Arguments ..
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM),
     +                 ZROE(NOFVAR,*),
     +                 STIFEL(NOFVAR,NOFVAR,NOFVERT,NOFVERT),
     +                 NODRES(NOFVAR,NOFVERT),NEWRES(NOFVAR,NOFVERT),
     +                 VCZ(NOFVAR,NOFVERT),VCDZ(NOFVAR,NOFVERT)
      INTEGER IBNDFAC(3,NBFAC) 
C     ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLUME,EPS
      INTEGER IVERT,IBC,ICLOR,IELEM,IFACE,N,IV,JV,
     +        IFAIL,I,J,LOCA,LOCB,IVAR,JVAR
      LOGICAL COMPRESSIBLE,LWK
      INTEGER NERR,IOPT,NP
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
      EXTERNAL BC5
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
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL
      EXTERNAL DDOT,DNRM2,ICYCL
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
      NP = NPOIN + NGHOST + NPNOD
C
      DO 100 IFACE = 1,NBFAC
C
         IELEM = IBNDFAC(1,IFACE)
         IVERT = IBNDFAC(2,IFACE)
         ICLOR = IBNDFAC(3,IFACE)
C
         IBC = ICOLOR(ICLOR)
C
C     Do NOT apply (inviscid) boundary conditions if:
C
C     IBC = 1 supersonic inlet
C     IBC = 3 supersonic outlet
C     IBC = 6 viscous wall
C     IBC = 7 prescribed velocity profile (maybe I should)
C
         IF (IBC.NE.BC_TYPE_SLIP_FREE)GOTO 100
C
         CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL,
     +   ZROE, VFACNOR, NDIM, NOFVERT, NOFVAR, NP, ICN,
     +   VCZ, VCN, VOLUME)
C
      INCLUDE 'ptest.f'
C
         IF(IVERT.NE.NOFVERT)THEN
            CALL DSWAP(NOFVAR,VCZ(1,IVERT),1,VCZ(1,NOFVERT),1)
            LOCA = IDX(IVERT,1,NDIM)
            LOCB = IDX(NOFVERT,1,NDIM)
            CALL DSWAP(NDIM,VCN(LOCA),1,VCN(LOCB),1)
            LOCA = ICN(NOFVERT)
            ICN(NOFVERT) = ICN(IVERT)
            ICN(IVERT) = LOCA
            LWK = PFLAG(NOFVERT)
            PFLAG(NOFVERT) = PFLAG(IVERT)
            PFLAG(IVERT) = LWK
         ENDIF
C
         CALL DCOPY(NOFVAR,U_INFTY,1,VCZ(1,NOFVERT),1)
C
         CALL DINIT(NOFVAR*(NOFVERT-1),ZERO,TSTEP,1)
C
         CALL BC4VII(IELEM,VCN,VCZ,STIFEL,NODRES,TSTEP,
     +             NDIM,NOFVERT,NOFVAR)
C
         IF(.NOT.TIMEIMPL)GOTO 50
#if 0
      do ivert=1,nofvert-1
      write(6,*)ielem,ivert,ivert
      CALL X04CAF('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
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
             ENDIF ! compressible
C
#if 0
      CALL X04CAF('General',' ',Nofvar,Nofvert,vcz(1,IVERT),
     +            Nofvar,'Z ',IFAIL)
      CALL X04CAF('General',' ',Nofvar,Nofvert,vcdz(1,IVERT),
     +            Nofvar,'DZ ',IFAIL)
#endif
C
             CALL BC4VII(IELEM,VCN,VCDZ,DUMMY,NEWRES,WORK,
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
#ifdef MPI
          CALL MatSetValuesBlockedLocal(PetscMatA,1,ICN(IVERT),1,
     +                                  ICN(IVERT),STIFEL(1,1,IVERT,1),
     +                                  ADD_VALUES,IFAIL)
#else
          CALL MatSetValuesBlocked(PetscMatA,1,ICN(IVERT),1,ICN(IVERT),
     +                             STIFEL(1,1,IVERT,1),
     +                             ADD_VALUES,IFAIL)
#endif
#ifdef PRINT
      write(6,*)ielem,ivert,ivert
      CALL X04CAF('General',' ',Nofvar,Nofvar,STIFEL(1,1,IVERT,1),
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
         CALL VecSetValuesBlocked(PetscVecRHS,NOFVERT-1,ICN,NODRES,
     +                      ADD_VALUES,IFAIL)
C
  100 CONTINUE ! End loop over boundary faces
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
      END
