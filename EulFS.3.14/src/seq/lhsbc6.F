!> \details
!> @param[in] A is the PETSc jacobian matrix
!> @param[in] X is the FORTRAN vector of CONSERVED variables
!> @param[in] NodalBcs is the PETSc IS where nodal bndry conditions are stored
!> @param[in] IBGN the order of the system to be solved
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR nof degrees of freedom within each gridpoint
!> @param[in] COMPRESSIBLE is currently unused
C
      subroutine isothwall(A,X,NodalBcs,IBGN,NDIM,NOFVAR,COMPRESSIBLE)
C
C     $Id: lhsbc6.F,v 1.15 2020/04/23 09:46:19 abonfi Exp $
C
C     set b.c. for isothermal wall nodes in the Jacobian
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscis.h"
      use petscvec
      use petscmat
      use petscis
C
      implicit none
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
      INCLUDE 'pfcgas.com'
      INCLUDE 'stream.com'
      INCLUDE 'visco.com'
C
C
C  Input/output parameters:
      double precision x(*)
      integer     NDIM,NOFVAR,IBGN
      Mat A
      PetscOffset i_is
      IS, dimension(0:*) :: NodalBcs
      LOGICAL COMPRESSIBLE

C  Local variables:
      integer maxcols
      parameter(maxcols=5)
      double precision alpha,y(maxcols),help
      integer     iloc,i,IFAIL,nrows,irow(maxcols),icol(maxcols)
      integer     is_array(1)
C
C
      ALPHA = RSTAR*TWALL/GM1
C
      call ISGetIndices(NodalBcs(Isothermal),is_array,i_is,IFAIL)
      call ISGetSize(NodalBcs(Isothermal),Nrows,IFAIL)
C
      DO 1 I = 1, Nrows
         irow(1) = is_array(i_is + I)
         icol(1) = irow(1)-1
         iloc = irow(1) - IBGN + 1! local node numbering: adds one because FORTRAN based
         help = -X(iloc)/X(iloc-1)
C
C        help equals alpha as long as the initial temperature on the walls is twall
C        this is not any longer true when restarting from a solution with a different twall
C
!        y(1) = -alpha
         y(1) =  help
!        write(6,*)I,'i = ',irow(1),'j = ',icol(1),y(1),help,ibgn
         CALL MatSetValues(A,1,irow,1,icol,y,INSERT_VALUES,IFAIL)
    1 CONTINUE 
C
C
      CALL ISRestoreIndices(NodalBcs(Isothermal),is_array,i_is,IFAIL)
C
C
      RETURN
      END
