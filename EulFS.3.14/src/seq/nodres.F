      SUBROUTINE NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     1                  NELEM,NFACE,NBFAC,ZROE,PetscDT,PetscRHS,PetscA,
     2                  SCALARSCHEME,MATRIXSCHEME,TModelScheme,
     3                  EulerModel,NSModel,TurbulenceModel)
C
C     This routine computes the NODal RESidual, i.e. it
C     assembles the r.h.s. (PetscRHS), the nodal timestep (PetscDT)
C     and the matrix (A)
C
C
C     $Id: nodres.F,v 1.61 2020/04/23 09:37:34 abonfi Exp $
C
C#define DEBUG
C#define SKIP_BNDRY_VFLUX
C#define SKIP_BNDRY_WFLUX
C#define SKIP_BNDRY_FFLUX
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscmat
      use petscvec
C
      IMPLICIT NONE
C
C     Scalar Arguments ..
C
      INTEGER NDIM,NOFVERT,NOFVAR,NPOIN,NGHOST,NPNOD,NELEM,NFACE,NBFAC,
     &NTURB
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NOFVAR  number of variables (degrees of freedom)
C             in each meshpoint
C     NPOIN   no. of interior (processor owned) meshpoints; 
C             global number of meshpoints in the uni-processor case
C     NGHOST  no. of ghost meshpoints on the current processor; 
C             0 the uni-processor case
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C     NFACE   number of edges/faces in the mesh for the current 
C             processor (multi-processor case) or global number 
C             of edges/faces in the mesh (uni-processor case).
C     NBFAC   no. of processor owned boundary faces/edges;
C             global number of boundary faces/edges
C             in the uni-processor case.
C
C     PETSc stuff
C
      Mat PetscA
      Vec PetscDT,PetscRHS
      DOUBLE PRECISION ZROE(NOFVAR*NPOIN)
      PetscOffset VEC_I
      PetscLogDouble t1beg,t1elapsed,t2beg,t2elapsed
      PetscLogDouble t3beg,t3elapsed 
      PetscScalar VEC_V(1),         s1
      integer ivar
#ifdef DEBUG
      PetscTruth flg
      PetscScalar s1,s2
#endif
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bctype.h'
      INCLUDE 'datatype.h'
      INCLUDE 'time.h'
      INCLUDE 'implicit.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'visco.com'
      INCLUDE 'time.com'
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
      INCLUDE 'bnd.com'
      INCLUDE 'nloc.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'flags.com'
      INCLUDE 'stream.com'
      INCLUDE 'conv.com'
      INCLUDE 'io.com'
      INTEGER MY_PE 
      COMMON/MPICOM/MY_PE
C
C     External Arguments ..
C
      EXTERNAL EulerModel,NSModel,TurbulenceModel,
     >         ScalarScheme,MatrixScheme,TModelScheme
C
C     External Subroutines ..
C
      EXTERNAL BVSFLX2,BVSFLX4,BC5VII,BC5VIII,BC28VII,BC28VIII,BC4VII
      EXTERNAL BC5VII4Ar,BC28VII4Ar
C
C     Local Scalars ..
C
      INTEGER ifail,ICLR,IUNIT,IDIM
      INTEGER ISTKGT
      INTEGER LWORK(7)
      LOGICAL LFLAG(2)
      double precision s
C
C     Sets residual and local timestep to ZERO
C
#ifdef DEBUG
      call VecValid(PetscDT,flg,iFAIL)
       if (flg .ne. PETSC_TRUE) then
          if (MY_PE .eq. 0) then
             write(6,*) 'Corrupted vector DT!'
          endif
          STOP
          SETERRQ(1,' ',IFAIL)
       endif
      CALL VecSet(PetscDT,ZERO,IFAIL)
      call VecValid(PetscRHS,flg,iFAIL)
       if (flg .ne. PETSC_TRUE) then
          if (MY_PE .eq. 0) then
             write(6,*) 'Corrupted vector RHS!'
          endif
          STOP
          SETERRQ(1,' ',IFAIL)
       endif
#endif
C
      CALL VecSet(PetscRHS,ZERO,IFAIL)
C
C     set ghost values to zero
C     note that VecSet only sets values belonging to the
C          owning processor 
C     it is necessary to set the ghost point values using
C     the LOCAL representation
C
      CALL VecGetArray(PetscRHS,vec_v,vec_i,IFAIL) 
      CALL DINIT(NGHOST*NOFVAR,0.d0,vec_v(vec_i+NPOIN*NOFVAR+1),1)
      CALL VecRestoreArray(PetscRHS,vec_v,vec_i,IFAIL) 
C
C     get local representation of the time-step vector
C     and set ghost values to 0.d0
C
      CALL VecGetArray(PetscDT,vec_v,vec_i,IFAIL)
      CALL DINIT(NGHOST*NOFVAR,ZERO,vec_v(vec_i+NPOIN*NOFVAR+1),1)
      CALL VecRestoreArray(PetscDT,vec_v,vec_i,IFAIL)
C
C     Set entries of the matrix to 0.d0
C     observe that at the first iteration the memory has not been allocated yet !
C
      IF( TIMEIMPL ) THEN
          CALL MatZeroEntries(PetscA,IFAIL)
          IF(IFAIL.NE.0)THEN
             WRITE(NOUT,*)'MatZeroEntries has ret = ',IFAIL
          ENDIF
      ENDIF
c
#ifdef DEBUG
         CALL VecNorm(PetscRHS,NORM_2,s1,IFAIL)
      IF(MY_PE.EQ.0)THEN
         write(32,*)'RHS norm before all calculations ',s1
!        write(32,*)'A   norm before all calculations ',s2
      ENDIF
#endif
c
      CALL DINIT(3*(MBODIES+1),ZERO,PRESF,1)
      CALL DINIT(3*(MBODIES+1),ZERO,VISCF,1)
C     
      CALL PetscTime(t1beg,IFAIL)
C
C     the hybrid model has been removed
C
      LWORK(1) = ISTKGT((NOFVAR*NOFVERT)**2,KIND_REAL8) ! jacobian matrix
      LWORK(2) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8) ! nodal residual
      LWORK(3) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8) ! perturbed nodal residual
      LWORK(4) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8) ! nodal values
      LWORK(5) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8) ! perturbed nodal values
      LWORK(6) = ISTKGT(NDIM*NOFVERT*NTIMLEVS,KIND_REAL8) ! grid velocities
      CALL DINIT((NOFVAR*NOFVERT)**2,ZERO,DSTAK(LWORK(1)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(2)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(3)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(4)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(5)),1) 
      CALL DINIT(NDIM*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(6)),1) 
C
      CALL SetupRHS(ISTAK(LCELNOD),ISTAK(LCELFAC),DSTAK(LFACNOR),
     +              DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,DSTAK(LTURB),
     +              DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +              DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +              NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NTURB,
     +              PetscDT,PetscRHS,PetscA,
     +              ScalarScheme,MatrixScheme,TModelScheme,
     +              EulerModel,NSModel,TurbulenceModel)
C
      CALL ISTKRL(6)
C
      CALL PetscTime(t2beg,IFAIL)
      t1elapsed = t2beg-t1beg
C
C     If a scalar problem goto 50
C
      IF(ABS(KAN).EQ.1)GOTO 50
C
      IF(RAD_EQUI)THEN
         LWORK(1) = ISTKGT((NOFVAR+NCL)/2,KIND_INTEGER)
         LWORK(2) = ISTKGT(NOFVAR*NCL,KIND_REAL8)
         CALL CIRCUM_AV(DSTAK(LCLA),ISTAK(LCLJA),ISTAK(LCLIA),
     &                  ISTAK(LCLDEG),DSTAK(LCLRC),DSTAK(LZROE),
     &                  DSTAK(LCLZB),DSTAK(LWORK(2)),NOFVAR,NCL,LCLDA)
         CALL ISRE(DSTAK(LCLRC),DSTAK(LCLZB),NOFVAR,NCL)
         CALL TRANS(DSTAK(LCLZB),NOFVAR,NCL,NOFVAR*NCL,
     &              ISTAK(LWORK(1)),(NOFVAR+NCL)/2,IFAIL)
         CALL ISTKRL(2)
      ENDIF
C
C     Apply Boundary Conditions for the inviscid fluxes
C
      LWORK(1) = ISTKGT((NOFVAR*NOFVERT)**2,KIND_REAL8)
      LWORK(2) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8)
      LWORK(3) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8)
      LWORK(4) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8)
      LWORK(5) = ISTKGT(NOFVAR*NOFVERT*NTIMLEVS,KIND_REAL8)
      LWORK(6) = ISTKGT(NDIM*NOFVERT*NTIMLEVS,KIND_REAL8) ! grid velocities
C
C     it is important to set the nodal residual to ZERO
C     otherwise we run into troubles with coupled RANS: this
C     is because non-zero values previously computed would be re-inserted
C     into the RHS vector
C
      CALL DINIT((NOFVAR*NOFVERT)**2,ZERO,DSTAK(LWORK(1)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(2)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(3)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(4)),1) 
      CALL DINIT(NOFVAR*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(5)),1) 
      CALL DINIT(NDIM*NOFVERT*NTIMLEVS,ZERO,DSTAK(LWORK(6)),1) 
C
#ifdef SKIP_BNDRY_FFLUX
      write(6,*)"skipping inviscid far-field fluxes  !!!!!"
#else
C
C     Far-field boundary conditions
C
      IF    ( ABS(KAN) .EQ. 2 )THEN ! INcompressible flow equations
         CALL GHOST(BC5VIII,BC_TYPE_FAR_FIELD,ISTAK(LCELNOD),
     +           ISTAK(LCELFAC),DSTAK(LFACNOR),DSTAK(LXYZDOT),
     3           DSTAK(LVOL),ZROE,
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA)
      ELSEIF( ABS(KAN) .EQ. 4 )THEN ! Compressible flow equations
         CALL GHOST(BC5VII,BC_TYPE_FAR_FIELD,ISTAK(LCELNOD),
     +           ISTAK(LCELFAC),DSTAK(LFACNOR),DSTAK(LXYZDOT),
     3           DSTAK(LVOL),ZROE,
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA)
      ELSEIF( KAN .EQ. -3)THEN !Plasma flow equations
         CALL GHOST(BC5VII4AR,BC_TYPE_FAR_FIELD,ISTAK(LCELNOD),
     +           ISTAK(LCELFAC),DSTAK(LFACNOR),DSTAK(LXYZDOT),
     3           DSTAK(LVOL),ZROE,
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA)          
      ELSE 
          STOP 'INVALID KAN IN SUBROUTINE NODRES'
      ENDIF
caldo
C
C     Inflow-outflow boundary conditions
C
      IF(IBCTYPE.EQ.BC_TYPE_WEAK)THEN
caldo
         IF    ( ABS(KAN) .EQ. 2 )THEN ! INCompressible flow equations
           CALL GHOST2(BC28VIII,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +           DSTAK(LFACNOR),DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     &           DSTAK(LCORG),DSTAK(LCLRC),DSTAK(LCLZB),
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,ISTAK(LKLIST),DSTAK(LVLIST),
     +           PetscDT,PetscRHS,PetscA)
         ELSEIF( ABS(KAN) .EQ. 4 )THEN ! Compressible flow equations
           CALL GHOST2(BC28VII,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +           DSTAK(LFACNOR),DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     &           DSTAK(LCORG),DSTAK(LCLRC),DSTAK(LCLZB),
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,ISTAK(LKLIST),DSTAK(LVLIST),
     +           PetscDT,PetscRHS,PetscA)
         ELSEIF( KAN .EQ. -3)THEN !Plasma flow equations
           CALL GHOST2(BC28VII4Ar,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +           DSTAK(LFACNOR),DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     &           DSTAK(LCORG),DSTAK(LCLRC),DSTAK(LCLZB),
     +           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +           ISTAK(LBNDFAC),NBFAC,ISTAK(LKLIST),DSTAK(LVLIST),
     +           PetscDT,PetscRHS,PetscA)
         ELSE 
           STOP 'INVALID KAN IN SUBROUTINE NODRES'
         ENDIF
caldo
      ENDIF ! IBCTYPE
caldo
#endif
#ifdef SKIP_BNDRY_WFLUX
      write(6,*)"skipping inviscid wall fluxes  !!!!!"
#else
C
C     weak type inviscid wall boundary conditions
C
C     mirror type inviscid wall boundary conditions
C     currently implemented for compressible flows only
C
      IF(SLIP_FREE_BC_TYPE.EQ.BC_TYPE_MIRROR)THEN
         IF(ABS(KAN).EQ.4)THEN 
         CALL GHOST(BC4VII,BC_TYPE_SLIP_FREE,ISTAK(LCELNOD),
     2           ISTAK(LCELFAC),DSTAK(LFACNOR),DSTAK(LXYZDOT),
     3           DSTAK(LVOL),ZROE,
     3           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     5           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     6           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA)
         ELSEIF(ABS(KAN).EQ.2)THEN 
            STOP 'Mirror bcs not yet available for Incompressible flows'
!        CALL GHOST(BC4VII,BC_TYPE_SLIP_FREE,ISTAK(LCELNOD),
!    2           ISTAK(LCELFAC),DSTAK(LFACNOR),DSTAK(LXYZDOT),
!    3           DSTAK(LVOL),ZROE,
!    3           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
!    +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
!    5           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
!    6           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA)
         ELSEIF( KAN .EQ. -3)THEN 
            STOP 'Mirror bcs not yet available for Plasma'
         ELSE
            STOP 'DO NOT KNOW WHAT TO DO (1)'
         ENDIF
C
C     weak wall boundary conditions for inviscid walls
C
      ELSEIF(SLIP_FREE_BC_TYPE.EQ.BC_TYPE_FLUX)THEN
         LFLAG(1) = (ABS(KAN).EQ.4) ! Compressible perfect gas
         LFLAG(2) = (ABS(KAN).EQ.3) ! Compressible plasma (Argon)
         CALL WEAKBC(ISTAK(LCELNOD),ISTAK(LCELFAC),DSTAK(LFACNOR),
     2           DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     3           DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +           DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     5           NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     6           ISTAK(LBNDFAC),NBFAC,PetscDT,PetscRHS,PetscA,
     7           LFLAG(1))
      ELSE
            WRITE(6,*)SLIP_FREE_BC_TYPE
            STOP 'DO NOT KNOW WHAT TO DO (2)'
      ENDIF
C
      CALL PetscTime(t3beg,IFAIL)
      t2elapsed = t3beg-t2beg
#endif
C
C     Apply boundary conditions for the viscous terms
C
#ifdef SKIP_BNDRY_VFLUX
      write(6,*)"skipping boundary viscous fluxes  !!!!!"
#else
      IF    ( KAN .EQ. 2 )THEN
         CALL BNDVFLX(BVSFLX2,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +                DSTAK(LFACNOR),DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     3                DSTAK(LTURB),
     +                DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +                DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +                NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +                NTURB,ISTAK(LBNDFAC),NBFAC,TurbulenceModel,
     +                PetscDT,PetscRHS,PetscA)
      ELSEIF( KAN .EQ. 4 )THEN
         CALL BNDVFLX(BVSFLX4,ISTAK(LCELNOD),ISTAK(LCELFAC),
     +                DSTAK(LFACNOR),DSTAK(LXYZDOT),DSTAK(LVOL),ZROE,
     3                DSTAK(LTURB),
     +                DSTAK(LWORK(1)),DSTAK(LWORK(2)),DSTAK(LWORK(3)),
     +                DSTAK(LWORK(4)),DSTAK(LWORK(5)),DSTAK(LWORK(6)),
     +                NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     +                NTURB,ISTAK(LBNDFAC),NBFAC,TurbulenceModel,
     +                PetscDT,PetscRHS,PetscA)
      ENDIF 
#endif
C
C     release workarrays
C
      CALL ISTKRL(6)
C
      CALL PetscTime(t3elapsed,IFAIL)
      t3elapsed = t3elapsed-t3beg
C
C
C     write aerodynamic coefficients
C
      DO 33 ICLR = 0, NCOLOR
C
C     Skip the current color unless slip/no-slip body
C
         IF( MCOLOR(ICLR) .NE. 0 .AND.
     +     ( ICOLOR(ICLR,1) .EQ. BC_TYPE_SLIP_FREE .OR. 
     2       ICOLOR(ICLR,1) .EQ. BC_TYPE_NO_SLIP ) )THEN
C
      IUNIT = IMUNIT(ICLR)
C
C     Write Fx(pressure) Fx(skin) Fy(pressure) Fy(skin)
C
      WRITE(IUNIT,2100)NITER,(PRESF(IDIM,ICLR),VISCF(IDIM,ICLR),
     +              IDIM=1,3),ITIM
C
      ENDIF
   33 CONTINUE
   50 CONTINUE
C
C     write timings:
C
      WRITE(ITIM1,*)NITER,T1ELAPSED,T2ELAPSED,T3ELAPSED
C
      RETURN
 2100 FORMAT(I6,1X,3(1X,E20.14,1X,E12.6),1X,I6)
C
C
      END
