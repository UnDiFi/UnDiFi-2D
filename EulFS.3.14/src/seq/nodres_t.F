      SUBROUTINE NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,
     +           NGHOST,NPNOD,NELEM,NFACE,NBFAC,DT,ZTURB,RHS,A,
     +           Model,ScalarScheme,MatrixScheme)
C
C     This routine computes the NODal RESidual, i.e. it
C     assembles the r.h.s. (RHS), the nodal timestep (DT)
C     and the matrix (A) for the turbulence model
C
C
C     $Id: nodres_t.F,v 1.20 2020/04/23 09:37:34 abonfi Exp $
C
CC#define DEBUG
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec 
      use petscmat 
C
      IMPLICIT NONE
C
C     PETSc stuff
C
      Mat A
      Vec DT,RHS
      DOUBLE PRECISION ZTURB(*)
      PetscLogDouble t1beg,t1elapsed,t2beg
      PetscScalar VEC_V(1),S
      PetscOffset VEC_I
C
      INCLUDE 'paramt.h'
      INCLUDE 'datatype.h'
      INCLUDE 'bnd.h'
      INCLUDE 'implicit.h'
      INCLUDE 'constants.h'
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
      INCLUDE 'nloc.com'
      INCLUDE 'flags.com'
      INCLUDE 'stream.com'
      INCLUDE 'conv.com'
      INCLUDE 'nboun.com'
      INCLUDE 'io.com'
      integer my_pe 
      common/MPICOM/my_pe
C
C     Scalar Arguments ..
C
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     2NBFAC
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
C     External Arguments ..
C
      EXTERNAL Model,ScalarScheme,MatrixScheme
C
C     External Subroutines ..
C
C
C     Local Scalars ..
C
      INTEGER IFAIL
      INTEGER ISTKGT
      INTEGER LWORKA,LWORKB,LWORKC,LWORKD,LWORKE ,LWORKF
C
C     Sets residual and local timestep to ZERO
C
      CALL VecSet(DT,ZERO,IFAIL)
      CALL VecSet(RHS,ZERO,IFAIL)
C
C     set ghost values to 0.d0
C
      CALL VecGetArray(DT,vec_v,vec_i,IFAIL)
      CALL DINIT(NGHOST,ZERO,vec_v(vec_i+NPOIN+1),1)
      CALL VecRestoreArray(DT,vec_v,vec_i,IFAIL)
C
C     set ghost values to zero
C     note that VecSet only sets values belonging to the
C          owning processor 
C     it is necessary to set the ghost point values using
C     the LOCAL representation
C
      CALL VecGetArray(RHS,vec_v,vec_i,IFAIL) 
      CALL DINIT(NGHOST*NTURB,ZERO,vec_v(vec_i+NPOIN*NTURB+1),1)
      CALL VecRestoreArray(RHS,vec_v,vec_i,IFAIL) 
C
C     Set entries of the matrix to 0.d0
C
      IF( TIMEIMPL )THEN
         CALL MatZeroEntries(A,IFAIL)
      ENDIF
C     
      CALL PetscTime(t1beg,IFAIL)
C
C     can change nofvar into nturb
C
      LWORKA = ISTKGT((NOFVAR*NOFVERT)**2,KIND_REAL8)
      LWORKB = ISTKGT(NOFVAR*NOFVERT,KIND_REAL8)
      LWORKC = ISTKGT(NOFVAR*NOFVERT,KIND_REAL8)
      LWORKD = ISTKGT(NOFVAR*NOFVERT,KIND_REAL8)
      LWORKE = ISTKGT(NOFVAR*NOFVERT,KIND_REAL8)
      LWORKF = ISTKGT(NDIM*NOFVERT,KIND_REAL8)
C
C
      CALL SetupRHS_t(ISTAK(LCELNOD),ISTAK(LCELFAC),DSTAK(LFACNOR),
     +              DSTAK(LXYZDOT),
     +              DSTAK(LVOL),DSTAK(LZROE),ZTURB,
     +              DSTAK(LWORKA),DSTAK(LWORKB),DSTAK(LWORKC),
     +              DSTAK(LWORKD),DSTAK(LWORKE),DSTAK(LWORKF),
     +              NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NTURB,
     +              DT,RHS,A,
     +              ScalarScheme,MatrixScheme,Model)
C
      CALL ISTKRL(6)
C
#ifdef DEBUG
      CALL VecNorm(RHS,NORM_2,S,IFAIL)
      WRITE(6,*)'L-2 norm of the rhs is ',S
      CALL VecNorm(DT,NORM_2,S,IFAIL)
      WRITE(6,*)'L-2 norm of the DT is ',S
#endif
C
      CALL PetscTime(t2beg,IFAIL)
      t1elapsed = t2beg-t1beg
C
C
      RETURN
      END
