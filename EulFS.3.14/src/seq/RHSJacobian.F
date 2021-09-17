      subroutine RHSJacobian(ts,t,global_in,A,B,str,ctx)
C
C     $Id: RHSJacobian.F,v 1.7 2020/04/23 09:37:34 abonfi Exp $
C
CCCCC#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
CCCCC#include "petsc/finclude/petscksp.h"
CCCCC#include "petsc/finclude/petscpc.h"
CCCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscts.h"
              use petscvec
              use petscmat
              use petscts
C
      implicit none
C
      INCLUDE 'implicit.h'
      INCLUDE 'flags.com'
      INCLUDE 'visco.com'
C
C  Input/output parameters:
      TS               ts
      double precision t
      Vec              global_in
      Mat              A,B
      MatStructure     str
      integer          ctx

C  Local variables:
      integer          IFAIL
C
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C  Complete the matrix assembly process and set some options
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C
      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL) 
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL) 
C
      RETURN
      END
