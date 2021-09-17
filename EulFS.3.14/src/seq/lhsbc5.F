C
      SUBROUTINE LHSBC5( B, NodalBcs )
C
C     This routine applies the Dirichlet boundary conditions to the
C     stiffness matrix for the turbulence one-eqn. model.
C     For a no-slip wall node (NODCODE=1000) A_ij j \neq i is set to 0.d0
C     and the diagonal entry is left unchanged
C
C     $Id: lhsbc5.F,v 1.20 2020/04/23 09:46:19 abonfi Exp $
C
C
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscis.h"
      use petscmat
      use petscis
C
      IMPLICIT NONE
C
      INCLUDE "constants.h"
      INCLUDE "bnd.h"
      INCLUDE "iset.h"
C
C
      IS, dimension(0:*) :: NodalBcs
      Mat B
C
      INTEGER IFAIL
C
C     zero-out rows of nodes where no-slip and free-stream nodal bcs
C     apply
C
      CALL MatZeroRowsIS(B,NodalBcs(NoSlipNodes),ONE,PETSC_NULL_VEC,
     &                   PETSC_NULL_VEC,IFAIL)
      CALL MatZeroRowsIS(B,NodalBcs(FreestreamTurbulence),ONE,
     &                   PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
C
      CALL MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,IFAIL)
      CALL MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,IFAIL)
C
C     This will remove the entries from the matrix; but when 
C     you assemble the matrix for the next time-step it will not 
C     insert those unneeded values (the values that would be zeroed
C     in the next call to MatZeroRows()). This is a very efficient
C     way to do multiple linear systems BUT assumes
C
C     1) the rows you want to zero remain the same for each successive
C        linear system in the series
C     2) you are not adding new nonzeros in other parts of the matrix
C        at later timesteps (because the above option will cause those 
C        new values to be ignored).
C
      CALL MatSetOption(B,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE,IFAIL)
C
      RETURN
      END
C
