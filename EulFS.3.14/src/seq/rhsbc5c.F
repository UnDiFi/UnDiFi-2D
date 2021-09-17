      SUBROUTINE RHSBC5C( RHS, ROWBGN, NodalBcs )
C
C     $Id: rhsbc5c.F,v 1.7 2020/04/23 09:46:19 abonfi Exp $
C
C
C     This routine sets the boundary conditions in the
C     "EXPLICIT" residual for the turbulence equation ..
C
#include "petsc/finclude/petscis.h"
      use petscis
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
C
      IS, dimension(0:*) :: NodalBcs
C
      PetscOffset RHS_I,IDX_I
C
      INTEGER I
C
      DOUBLE PRECISION RHS(*)
C
      INTEGER IPOIN,NI,IFAIL,ROWBGN
      PetscScalar RHS_V(1)
      INTEGER IDX_V(1)
C
C
C    EXPLICIT Treatment of no-slip wall nodes
C
C
      CALL ISGetIndices(NodalBcs(NoSlipNodes),IDX_V,IDX_I,IFAIL)
      CALL ISGetSize(NodalBcs(NoSlipNodes),NI,IFAIL)
C
      DO 1000 I= 1 , NI
C
          IPOIN = IDX_V( IDX_I+I )-ROWBGN  +1
c         write(6,*)ipoin,RHS(IPOIN)
          RHS(IPOIN) = ZERO
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
          IPOIN = IDX_V( IDX_I+I )-ROWBGN  +1
          RHS(IPOIN) = ZERO
C
 2000 CONTINUE
C
      CALL ISRestoreIndices(NodalBcs(FreestreamTurbulence),IDX_V,IDX_I,
     &                      IFAIL)
C
      RETURN
      END
