C
C   $Id: PCPierluigi.F,v 1.1 2001/07/19 09:59:33 abonfi Exp $
C

!/***********************************************************************/
!/*          Routines for a user-defined shell preconditioner           */
!/***********************************************************************/

!
!   SampleShellPCSetUp - This routine sets up a user-defined
!   preconditioner context.  
!
!   Input Parameters:
!   pmat  - preconditioner matrix
!   x     - vector
!
!   Output Parameter:
!   ierr  - error code (nonzero if error has been detected)
!
!   Notes:
!   In this example, we define the shell preconditioner to be Jacobi
!   method.  Thus, here we create a work vector for storing the reciprocal
!   of the diagonal of the preconditioner matrix; this vector is then
!   used within the routine SampleShellPCApply().
!
      subroutine PCMyASM(pmat,pmat1,x,ierr)

      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"

#include "iset.com"

      Vec     x
      Mat     pmat,pmat1
      IS      isolap
      integer ierr
      double  precision info(MAT_INFO_SIZE)
      PetscTruth flg

!  Common block to store data for user-provided preconditioner 

      PC  asm
      integer MY_PE
      common /MPICOM/ MY_PE
C
C
caldo CALL MPI_Barrier(PETSC_COMM_WORLD,ierr) 
caldo CALL MPI_Abort(PETSC_COMM_WORLD,-12,ierr) 
C
      call PCCreate(PETSC_COMM_WORLD,asm,ierr)
      call PCSetType(asm,PCASM,ierr)
      call PCASMSetLocalSubdomains(asm,1,isolap,ierr)
      call PCSetOperators(asm,pmat,pmat,SAME_NONZERO_PATTERN,
     &                     ierr)
      call PCSetUp(asm,ierr)

      end

