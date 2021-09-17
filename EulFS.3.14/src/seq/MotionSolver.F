!>  \brief
!>  the Motion solver: given the grid velocities at all boundaries,
!>  solves Laplace's equation to get the grid velocites at all interior gridpoints
!>
!> @param [in] MotionSolver
!> @param [in] PetscMatC
!> @param [in] NodalBcs
!> @param [in] NOFVERT number of vertices per element (=NDIM+1), since only triangles or tetrahedra are allowed)
!> @param [in] NPOIN   nof of interior (processor owned) meshpoints; global number of meshpoints in the uni-processor case
!> @param [in] NGHOST  nof of ghost meshpoints on the current processor; 0 the uni-processor case
!> @param [in] NPNOD  nof of periodic meshpoints on the current processor
!> @param [in] NELEM   nof of processor owned elements (triangles/tetrahedra); global number of elements in the uni-processor case
!> @param [in] NBFAC   nof of processor owned boundary faces/edges; global number of boundary faces/edges in the uni-processor case.
!> @param [in] NFACE   number of edges/faces in the mesh for the current processor (multi-processor case) or global number of edges/faces in the mesh (uni-processor case).
!> \author $Author: abonfi $
!> \version $Revision: 1.8 $
!> \date $Date: 2020/04/24 06:04:41 $
!> \warning The subroutine C_DISTIJ is un-optimized and will be VERY slow on large meshes
!> \warning This subroutine is un-tested in 3D
!> \bug This subroutine is NOT yet working in parallel
      SUBROUTINE GridDeformation(MotionSolver,PetscMatC,NodalBCS,
     &                           NDIM,NOFVERT,NPOIN,NGHOST,NPNOD,NELEM,
     3                           NBFAC,NFACE)
!
!  $Id: MotionSolver.F,v 1.8 2020/04/24 06:04:41 abonfi Exp $
!
C
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscpc.h"
      use petscpc
      use petscksp
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
      INCLUDE 'datatype.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
      INCLUDE 'nloc.com'
      INCLUDE 'stream.com'
      INCLUDE 'conv.com'
      INCLUDE 'nboun.com'
      INCLUDE 'io.com'
      integer my_pe 
      common/MPICOM/my_pe
C
      Mat PetscMatC
      Vec              X,RHS,dummy
      IS, dimension(0:*) :: NodalBcs
      KSP Motionsolver
      PetscScalar x_array(1)
      PetscOffset i_x

cxxx  COMMON/MYTSCOM/A,B,PetscMatC,dummy
C
      INTEGER NDIM,NOFVERT,NOFVAR,NPOIN,NGHOST,NPNOD,
     &NELEM,NFACE,NBFAC
      integer IFAIL,I,LIWORK(3),NR,NZR,BS
      double precision s 
C
C
C     Executable Statements ..
C
C
C     PETSc stuff
C
      PetscLogDouble t1beg,t1elapsed,t2beg
C
C     FORTRAN stack
C
      DOUBLE PRECISION DSTAK(1)
      COMMON /CSTAK/ DSTAK
      INTEGER ISTAK(1)
      EQUIVALENCE(DSTAK(1),ISTAK(1))
C
C
C     External Arguments ..
C
C
C     External Subroutines ..
C
      INTEGER ISTKGT
C
C     Local Scalars ..
C
      INTEGER LWORKA
C
C
C
C     Set entries of the matrix to 0.d0
C
      CALL MatZeroEntries(PetscMatC,IFAIL)
C     
ctime CALL PetscTime(t1beg,IFAIL)
C
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
!     CALL DINIT(NELEM,ONE,DSTAK(LWORKA),1) ! put 1 into the diffusion coefficient: just 4 debugging
C
      CALL SetupLHS_b(ISTAK(LCELNOD),ISTAK(LCELFAC),DSTAK(LFACNOR),
     +              DSTAK(LCORG),DSTAK(LXYZDOT),DSTAK(LVOL),
C                   ^
C                   |
C                   |
C          we pass the grid velocity because we have there 3 time levels in the un-steady case
C          otherwise there will be an access out-of-bounds in celptr
     +              DSTAK(LWORKA),
     +              NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,
     4              PetscMatC)
C
!     CALL R8Mat_Print('General',' ',NDIM,NPOIN,DSTAK(LXYZDOT),
!    +            NDIM,'Grid velocity ',IFAIL)
C
      CALL ISTKRL(1)
C
C     The matrix has now been assembled
C
      CALL KSPSetOperators(MotionSolver,PetscMatC,PetscMatC,IFAIL)
      CALL KSPSetFromOptions(MotionSolver,IFAIL) 
C
      CALL VecCreate(PETSC_COMM_WORLD,X,IFAIL)
#ifdef MPI
      CALL VecSetType(X,VECMPI,IFAIL)
#else
      CALL VecSetType(X,VECSEQ,IFAIL)
#endif
      CALL VecSetSizes(X,NPOIN,PETSC_DECIDE,IFAIL)
      CALL VecDuplicate(X,RHS,IFAIL)
C
      DO I = 0,NDIM-1
         CALL VecSet(RHS,ZERO,IFAIL)
         call VecGetArray(X,x_array,i_x,ifail)
         CALL DCOPY(NPOIN,DSTAK(LXYZDOT+I),NDIM,x_array(i_x+1),1)
         call VecRestoreArray(X,x_array,i_x,ifail)
!        CALL VecNorm(X,NORM_2,s,IFAIL)
!        write(6,*)'X norm before solving ',s
C
         CALL MatZeroRowsIS(PetscMatC,NodalBcs(MotionSolverBCS),ONE,
     &                   X,RHS,IFAIL)
C
      IF( .FALSE. )THEN
!     IF( .TRUE. )THEN
          LIWORK(1) = 1
          LIWORK(2) = 1
          LIWORK(3) = 1
          NR = NPOIN
          BS = 1
          CALL PrintMatMM(PetscMatC,RHS,X,
     3                    ISTAK(LIWORK(3)),NR,NZR,BS,ITER,0)
          LIWORK(3) = ISTKGT(NZR,KIND_INTEGER)! storage for ir
          CALL PrintMatMM(PetscMatC,RHS,X,
     3                    ISTAK(LIWORK(3)),NR,NZR,BS,ITER,100)
          write(6,*)'Beyond PrintMatMM'
          CALL ISTKRL(1)
          CALL EXIT(0)
      ENDIF
C
         CALL KSPSolve(MotionSolver,RHS,X,IFAIL)
!        CALL VecNorm(X,NORM_2,s,IFAIL)
!        write(6,*)'X norm after solving ',s
         call VecGetArray(X,x_array,i_x,ifail)
         CALL DCOPY(NPOIN,x_array(i_x+1),1,DSTAK(LXYZDOT+I),NDIM)
         call VecRestoreArray(X,x_array,i_x,ifail)
      ENDDO
C
      CALL VecDestroy(X,IFAIL)
      CALL VecDestroy(RHS,IFAIL)
#ifdef MPI
      STOP 'Motion solver has not been checked in parallel'
#endif
C
      RETURN
      END
