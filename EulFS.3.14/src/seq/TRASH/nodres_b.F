      SUBROUTINE NODRES_B(NDIM,NOFVERT,NPOIN,
     +           NGHOST,NPNOD,NELEM,NBFAC,NFACE,
     3           PetscMatC)
C
C     This routine computes the NODal RESidual, i.e. it
C     assembles the r.h.s. (PetscVecRHS),
C     and the matrix (A) for the Motion Solver
C
C
      IMPLICIT NONE
C
C     $Id: nodres_b.F,v 1.1 2013/06/28 09:09:15 abonfi Exp $
C
CC#define DEBUG
C
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
C
C     PETSc stuff
C
      Mat PetscMatC
      PetscLogDouble t1beg,t1elapsed,t2beg
!     PetscScalar VEC_V(1),S
!     PetscOffset VEC_I
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
      INTEGER NDIM,NOFVERT,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     2NBFAC
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
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
      INTEGER LWORKA,NBPOIN
C
C
C
C     Set entries of the matrix to 0.d0
C
      CALL MatZeroEntries(PetscMatC,IFAIL)
C     
ctime CALL PetscTime(t1beg,IFAIL)
C

!     CALL NODCOD(ISTAK(LNODCOD),NPOIN+NGHOST+NPNOD,NBPOIN,
!    +            ISTAK(LCELNOD),NOFVERT,NELEM,ISTAK(LBNDFAC),NBFAC)
C
C
      LWORKA = ISTKGT(NELEM,KIND_REAL8) ! <--- distance from the nearest boundary
C
C     We compute the distance of the cell centres from the closest boundary
C
      CALL C_DISTIJ(DSTAK(LCORG),NDIM,ISTAK(LNODCOD),
     &              NPOIN+NGHOST+NPNOD,DSTAK(LWORKA),
     &              ISTAK(LCELNOD),NOFVERT,NELEM)
C
      CALL SetupLHS_b(ISTAK(LCELNOD),ISTAK(LCELFAC),DSTAK(LFACNOR),
     +              DSTAK(LCORG),DSTAK(LXYZDOT),DSTAK(LVOL),
!                   ^
!                   |
!                   |
!          we pass the grid velocity because we have there 3 time levels in the un-steady case
!          otherwise there will be an access out-of-bounds in celptr
!
     +              DSTAK(LWORKA), ! <---- distance from the nearest boundary
     +              NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,
     4              PetscMatC)
C
      CALL ISTKRL(1)
C
C
ctime CALL PetscTime(t2beg,IFAIL)
ctime t1elapsed = t2beg-t1beg
C
C
      RETURN
      END
