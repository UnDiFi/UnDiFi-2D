!> \details
!> @param[in] ts the PETSc timestepping context, currently a dummy argument
!> @param[in] U is the PETSc vector of CONSERVED variables
!> @param[in] A is the PETSc jacobian matrix
!> @param[in] X is the PETSc unknown vector
!> @param[in] RHS is the PETSc rhs vector
!> @param[in] NodalBcs is the PETSc IS that addresses nodal bcs
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR is the nof dofs in each gridpoint
      SUBROUTINE JacobianBoundaryConditions(ts,U,A,X,RHS,NodalBcs,
     &                                      NDIM,NOFVAR)
C
C
C     N.B.: upon entry U contains the vector of CONSERVED
C           variables
C
C
C     $Id: JacobianBoundaryConditions.F,v 1.27 2020/04/23 09:46:19 abonfi Exp $
C
#include "petsc/finclude/petscsys.h"
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscmat.h"
CCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscts.h"
      use petscts
C
      IMPLICIT NONE
C
      Mat A
      Vec rhs
      Vec X ! vector of the unkwnouns
      Vec U
      TS ts
      PetscLogDouble telapsed,tbeg,tend
      PetscScalar x_array(1)
      PetscOffset i_x
#ifdef DEBUG
      integer      is_array(1)
      PetscOffset i_is
      integer maxcols
      parameter(maxcols=130)
      double precision vals(maxcols)
      INTEGER i,j,k,icol(maxcols),irow(5),lda,nrows,ncols
#endif
      INTEGER NDIM,NOFVAR
      INTEGER IFAIL,IBGN,IEND
!     INTEGER ISTAK(1)
!     DOUBLE PRECISION DSTAK(1)
!     COMMON/CSTAK/DSTAK
!     EQUIVALENCE (ISTAK(1),DSTAK(1))
C     ..
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bctype.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'nloc.com'
      INCLUDE 'ibc8.com'
C
      IS, dimension(0:*) :: NodalBcs 
C
C
      INTEGER ISTKGT
      EXTERNAL ISTKGT
C
      CALL PetscTime(tbeg,IFAIL)
C
C     the stiffness matrix has to be modified
C     to account for boundary conditions .. 
C
C     Note that the Index Sets address rows by GLOBAL number
C
      CALL MatZeroRowsIS(A,NodalBcs(SupersonicVariables),ONE,
     &                   PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
      CALL MatZeroRowsIS(A,NodalBcs(NoSlipVelocities),ONE,
     &                   PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
      IF(TURBULENT.AND.COUPLED)THEN
         CALL MatZeroRowsIS(A,NodalBcs(NoSlipNodes),ONE,
     &                      PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
         CALL MatZeroRowsIS(A,NodalBcs(FreestreamTurbulence),ONE,
     &                      PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
      ENDIF
C
C     ``strong'' boundary conditions for isothermal wall
C
      IF(IADIA.NE.0)THEN
          CALL MatZeroRowsIS(A,NodalBcs(Isothermal),ONE,
     &                       PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
          call VecGetOwnerShipRange(U,ibgn,iend,ifail)
          call VecGetArray(U,x_array,i_x,ifail)
C
C     Flush assembly (because of MatSetValues inside Isothwall)
C
          CALL MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,IFAIL)
          IF(IFAIL.NE.0)THEN
             WRITE(6,*)' MatAssemblyBegin has ret ',IFAIL
             CALL EXIT(IFAIL)
          ENDIF
          CALL MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY,IFAIL)
          IF(IFAIL.NE.0)THEN
             WRITE(6,*)' MatAssemblyEnd has ret ',IFAIL
             CALL EXIT(IFAIL)
          ENDIF
C
          CALL Isothwall(A,x_array(i_x+1),NodalBcs,IBGN,NDIM,NOFVAR,
     $         (ABS(KAN).EQ.4))
          call VecRestoreArray(U,x_array,i_x,ifail)
C
C     Final assembly (because of MatSetValues inside lhsbc)
C
          CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL)
          CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL)
      ENDIF
C
C
      IF( LAPLACE ) THEN 
          CALL SETBC4LAPLACE(A,X,RHS,U,NodalBcs)
      ENDIF
caldo
      CALL PetscTime(tend,IFAIL)
      telapsed=tend-tbeg
C
C
      RETURN
      END
