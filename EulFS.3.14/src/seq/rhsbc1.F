      SUBROUTINE RHSBC1( RHS, NodalBcs )
C
C     This routine sets the boundary conditions in the
C     "EXPLICIT" residual for scalar problems ..
C
C     $Id: rhsbc1.F,v 1.18 2020/04/23 09:46:19 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
      use petscvec
      use petscis
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
C
C
      INTEGER MY_PE 
      COMMON/MPICOM/MY_PE
C
      Vec RHS
      IS, dimension(0:*) :: NodalBcs
C
      INTEGER IPOIN,I,IFAIL,NI
      PetscScalar RHS_V(1)
      PetscOffset RHS_I,IDX_I
      INTEGER IDX_V(1),ROWBGN,ROWEND
C
C
C    EXPLICIT Treatment of dirichlet nodes
C
      CALL VecGetOwnerShipRange(RHS,ROWBGN,ROWEND,IFAIL)
      CALL VecGetArray(RHS,RHS_V,RHS_I,IFAIL)
C
      CALL ISGetSize(NodalBcs(SupersonicNodes),NI,IFAIL)
      CALL ISGetIndices(NodalBcs(SupersonicNodes),IDX_V,IDX_I,IFAIL)
C
      DO 1000 I= 1 , NI
C
          IPOIN = IDX_V(IDX_I+I)
          IPOIN = IPOIN-ROWBGN
          RHS_V(RHS_I+IPOIN+1) = ZERO
C
 1000 CONTINUE
C
      CALL VecRestoreArray(RHS,RHS_V,RHS_I,IFAIL)
      CALL ISRestoreIndices(NodalBcs(SupersonicNodes),IDX_V,IDX_I,IFAIL)
C
C
      RETURN
      END 
