!> \brief
!> This routine destroys PETSc Objects
!>
!> @param[in] ZROE PETSc Vec where Roe's parameter vector is stored
!> @param[in] VISCT PETSc Vec where turbulent viscosity is stored
!> @param[in] A PETSc Mat where the Jacobian matrix for the mean flow equations is stored
!> @param[in] B PETSc Mat where the Jacobian matrix for the turbulent transport equations is stored
!> @param[in] C PETSc Mat where the implicit matrix for Laplace's equations is stored
!> @param[in] NodalBcs PETSc Index Sets where the nodal boundary conditions are stored
!> @param[in] NOFVAR number of dofs and leading dimension of ZROE
!> \author $Author: abonfi $
!> \version $Revision: 1.14 $
!> \date $Date: 2020/04/24 06:04:16 $
C
      SUBROUTINE CLEARMEM( ZROE, VISCT, PetscMatrix, NodalBcs, NOFVAR )
C
C     $Id: clearmem.F,v 1.14 2020/04/24 06:04:16 abonfi Exp $
C
C
C
#include "petsc/finclude/petscpc.h"
      use petscpc
C
      IMPLICIT NONE
C
      INCLUDE 'bnd.h' 
      INCLUDE 'implicit.h'
      INCLUDE 'iset.h'
      INCLUDE 'flags.com'
      INCLUDE 'bnd.com'
      INCLUDE 'time.com'
C
#include "mapping.com"
C
      INTEGER NOFVAR 
      INTEGER IFAIL , I
      Mat, dimension(*) :: PetscMatrix
      Vec ZROE,VISCT
      IS, dimension(0:*) :: NodalBcs
C
C
      CALL VecDestroy(ZROE,IFAIL)
      IF(TURBULENT.AND.SEGREGATED) CALL VecDestroy(VISCT,IFAIL)
C
      IF( TIMEIMPL )THEN
          CALL MatDestroy(PetscMatrix(1),IFAIL)
          IF(TURBULENT.AND.SEGREGATED)THEN
              CALL MatDestroy(PetscMatrix(2),IFAIL)
          ENDIF
      ENDIF
C
      IF( ALE_LAPLACE_GRID )THEN
         CALL MatDestroy(PetscMatrix(3),IFAIL)
         CALL ISDestroy(NodalBcs(MotionSolverBCS),IFAIL)
      ENDIF 
C
C     clear Index sets
C
      IF(NOFVAR.EQ.1)THEN
         CALL ISDestroy(NodalBcs(SupersonicNodes),IFAIL)
      ELSE
         CALL ISDestroy(NodalBcs(SupersonicVariables),IFAIL)
      ENDIF
      CALL ISDestroy(NodalBcs(NoSlipVelocities),IFAIL)
      CALL ISDestroy(NodalBcs(Isothermal),IFAIL)
      IF(TURBULENT.AND.SEGREGATED)THEN
          CALL ISDestroy(NodalBcs(NoSlipNodes),IFAIL)
          CALL ISDestroy(NodalBcs(FreestreamTurbulence),IFAIL)
      ENDIF

      DO I = 0,NCOLOR
         IF( MCOLOR(I) .NE. 0 )CALL ISDestroy(NodalBcs(I),IFAIL)
      ENDDO
C
C
#ifdef MPI
      IF(TIMEIMPL)THEN
          CALL ISLocalToGlobalMappingDestroy(mapping(1),IFAIL)
          IF(TURBULENT.AND.SEGREGATED)
     +    CALL ISLocalToGlobalMappingDestroy(mapping(2),IFAIL)
      ENDIF 
#endif
C
      RETURN
 
      END
