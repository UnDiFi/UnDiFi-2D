      subroutine lhsbcX(A,X,IBGN,NDIM,NOFVAR,COMPRESSIBLE)
C
C     $Id: lhsbc6.F,v 1.5 2000/10/19 16:10:34 aldo Exp $
C
C     set b.c. for adiabatic wall nodes in the Jacobian
C
      implicit none
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscsles.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscts.h"
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants'
      INCLUDE 'visco.com'
      INCLUDE 'stream.com'
C
#include "iset.com"
C
C  Input/output parameters:
      double precision x(*)
      integer     NDIM,NOFVAR,IBGN
      Mat              A
      PetscOffset i_isA,i_isB
      LOGICAL COMPRESSIBLE

C  Local variables:
      integer maxcols
      parameter(maxcols=5)
      double precision alpha,y(maxcols)
      integer          i,IFAIL,nrows,irow(maxcols),icol(maxcols)
      integer     is_arrayA(1),is_arrayB(1)
      double precision vals(maxcols)
C
C
      call ISGetIndices(PeriodicNodesA,is_arrayA,i_isA,IFAIL)
      call ISGetIndices(PeriodicNodesB,is_arrayB,i_isB,IFAIL)
      call ISGetSize(PeriodicNodesA,Nrows,IFAIL)
C
      DO 1 I = 1, Nrows
         irow(1) = is_arrayA(i_isA + I)-IBGN
         icol(1) = is_arrayB(i_isB + I)-IBGN
         y(1) = -1.d0
#ifdef MPI
         CALL MatSetValuesLocal(A,1,irow,1,icol(1),y(1),INSERT_VALUES,
     &                          IFAIL)
#else
         CALL MatSetValues(A,1,irow,1,icol(1),y(1),INSERT_VALUES,IFAIL)
#endif
    1 CONTINUE 
C
C
      CALL ISRestoreIndices(PeriodicNodesA,is_arrayA,i_isA,IFAIL)
      CALL ISRestoreIndices(PeriodicNodesB,is_arrayB,i_isB,IFAIL)
C
C
      RETURN
      END
