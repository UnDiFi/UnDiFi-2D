      subroutine UserInitializeLinearSolver(m,n,userctx,ierr)

      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

      PetscInt m,n
      PetscErrorCode ierr
      PetscFortranAddr userctx(*)

      common /param/ hx2,hy2
      double precision hx2,hy2

!  Local variable declararions
      Mat     A
      Vec     b,x
      KSP    ksp
      PetscInt Ntot,five,one


!  Here we assume use of a grid of size m x n, with all points on the
!  interior of the domain, i.e., we do not include the points corresponding
!  to homogeneous Dirichlet boundary conditions.  We assume that the domain
!  is [0,1]x[0,1].

      hx2 = (m+1)*(m+1)
      hy2 = (n+1)*(n+1)
      Ntot = m*n

      five = 5
      one = 1

!  Create the sparse matrix. Preallocate 5 nonzeros per row.

      call MatCreateSeqAIJ(PETSC_COMM_SELF,Ntot,Ntot,five,              &
     &     PETSC_NULL_INTEGER,A,ierr)
!
!  Create vectors. Here we create vectors with no memory allocated.
!  This way, we can use the data structures already in the program
!  by using VecPlaceArray() subroutine at a later stage.
!
      call VecCreateSeqWithArray(PETSC_COMM_SELF,one,Ntot,              &
     &     PETSC_NULL_SCALAR,b,ierr)
      call VecDuplicate(b,x,ierr)

!  Create linear solver context. This will be used repeatedly for all
!  the linear solves needed.

      call KSPCreate(PETSC_COMM_SELF,ksp,ierr)

      userctx(1) = x
      userctx(2) = b
      userctx(3) = A
      userctx(4) = ksp
      userctx(5) = m
      userctx(6) = n

      return
      end
