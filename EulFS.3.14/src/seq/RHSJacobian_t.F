      subroutine RHSJacobian_t(ts,t,global_in,A,B,str,ctx)
C
C     $Id: RHSJacobian_t.F,v 1.7 2020/04/23 09:37:34 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscts.h"
              use petscts
              use petscvec
              use petscmat
C
      implicit none
C
C  Input/output parameters:
      TS               ts
      double precision t
      Vec              global_in
      Mat              A,B
      MatStructure     str
      integer          ctx

      integer          IFAIL
C
      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL) 
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL) 
C

      return
      end
