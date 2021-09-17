
      subroutine petsc_kernel_dgefa(n,a)
      implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

      integer n
      integer ierr
      double precision a(*)
      PetscReal shift
      PetscInt ipvt(15)
      MatScalar work(15)
      PetscErrorCode Kernel_A_gets_inverse_A_2
      PetscErrorCode Kernel_A_gets_inverse_A_3
      PetscErrorCode Kernel_A_gets_inverse_A_4
      PetscErrorCode Kernel_A_gets_inverse_A_5
      PetscErrorCode Kernel_A_gets_inverse_A_6
      PetscErrorCode Kernel_A_gets_inverse_A_7
      PetscErrorCode Kernel_A_gets_inverse_A_9
      PetscErrorCode Kernel_A_gets_inverse_A_15
      PetscErrorCode LINPACKdgedi

      goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)n
   1  goto 777
   2  ierr = Kernel_A_gets_inverse_A_2(a,%val(shift))
      return
   3  ierr = Kernel_A_gets_inverse_A_3(a,shift)
      return
   4  ierr = Kernel_A_gets_inverse_A_4(a,shift)
      return
   5  ierr = Kernel_A_gets_inverse_A_5(a,ipvt,work,shift)
      return
   6  ierr = Kernel_A_gets_inverse_A_6(a,shift)
      return
   7  ierr = Kernel_A_gets_inverse_A_7(a,shift)
      return
   8  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
   9  ierr = Kernel_A_gets_inverse_A_9(a,shift)
      return
  10  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
  11  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
  12  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
  13  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
  14  ierr = LINPACKdgedi(a,ipvt,work,shift)
      return
  15  ierr = Kernel_A_gets_inverse_A_15(a,ipvt,work,shift)
      return
CKernel_A_gets_inverse_A_9(MatScalar *a,PetscReal shift)
CKernel_A_gets_inverse_A_15(MatScalar *a,PetscInt *ipvt,MatScalar *work,PetscReal shift)
CKernel_A_gets_inverse_A_5(MatScalar *a,PetscInt *ipvt,MatScalar *work,PetscReal shift)
CKernel_A_gets_inverse_A_6(MatScalar *a,PetscReal shift)
CKernel_A_gets_inverse_A_7(MatScalar *a,PetscReal shift)
  777 write(6,*)'Invalid ',n,' in petsc_kernel '
      call exit(n) 
      end


