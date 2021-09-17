C
      SUBROUTINE RHSBC5( RHS, NodalBcs )
C
C     This routine sets the boundary conditions in the
C     "EXPLICIT" residual for the turbulence equation ..
C
C     $Id: rhsbc5.F,v 1.19 2020/04/23 09:46:19 abonfi Exp $
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
C
C
      Vec RHS
      IS, dimension(0:*) :: NodalBcs
      PetscOffset RHS_I,IDX_I
C
      INTEGER I
C
C
      INTEGER IPOIN,NI,IFAIL,ROWBGN,ROWEND
      PetscScalar RHS_V(1)
      INTEGER IDX_V(1)
C
C
C    EXPLICIT Treatment of no-slip wall nodes
C
C
      CALL VecGetOwnerShipRange(RHS,ROWBGN,ROWEND,IFAIL)
      CALL VecGetArray(RHS,RHS_V,RHS_I,IFAIL)
      CALL ISGetIndices(NodalBcs(NoSlipNodes),IDX_V,IDX_I,IFAIL)
      CALL ISGetSize(NodalBcs(NoSlipNodes),NI,IFAIL)
C
      DO 1000 I= 1 , NI
C
          IPOIN = IDX_V( IDX_I+I )-ROWBGN
          RHS_V(RHS_I+IPOIN+1) = ZERO
C
 1000 CONTINUE
C
      CALL ISRestoreIndices(NodalBcs(NoSlipNodes),IDX_V,IDX_I,IFAIL)
C
C
C
      CALL ISGetIndices(NodalBcs(FreestreamTurbulence),IDX_V,IDX_I,
     &                  IFAIL)
      CALL ISGetSize(NodalBcs(FreestreamTurbulence),NI,IFAIL)
C
      DO 2000 I= 1 , NI
C
          IPOIN = IDX_V( IDX_I+I )-ROWBGN
          RHS_V(RHS_I+IPOIN+1) = ZERO
C
 2000 CONTINUE
C
      CALL ISRestoreIndices(NodalBcs(FreestreamTurbulence),IDX_V,IDX_I,
     &                      IFAIL)
      CALL VecRestoreArray(RHS,RHS_V,RHS_I,IFAIL)
C
      RETURN
      END 
