/*
  $Id: petsc_kernel.c,v 1.2 2002/02/22 15:35:22 abonfi Exp $
*/

/*
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define solven_ SOLVEN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define solven_ solven
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define setsolven_ SETSOLVEN
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define setsolven_ setsolven
#endif
*/
#include <stdlib.h>
#include "petsc.h"

extern PetscErrorCode Kernel_A_gets_inverse_A_2(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_3(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_4(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_5(MatScalar *, PetscInt *, MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_6(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_7(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_9(MatScalar *, PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_15(MatScalar *, PetscInt *, MatScalar *, PetscReal);

static PetscErrorCode (*a[8])(MatScalar *mat, PetscReal);

int solven_(int *n, MatScalar *mat)
{
    switch(*n) {
		case 5:	
		    {
			PetscInt ipvt[5]:
			MatScalar work[5*5];
			return (int)Kernel_A_gets_inverse_A_5(mat, ipvt, work, (PetscReal)0.0);
			}
		case 15:
		    {
			PetscInt ipvt[15]:
			MatScalar work[15*15];
			return (int)Kernel_A_gets_inverse_A_15(mat, ipvt, work, (PetscReal)0.0);
			}
		default:
			return (int)a[(*n-2)](mat, (PetscReal) 0.0);
	}
}


void setsolven_(int *ierr)
{
    a[0] = Kernel_A_gets_inverse_A_2;
    a[1] = Kernel_A_gets_inverse_A_3;
    a[2] = Kernel_A_gets_inverse_A_4;
    a[3] = NULL;
    a[4] = Kernel_A_gets_inverse_A_6;
    a[5] = Kernel_A_gets_inverse_A_7;
	a[6] = NULL;
	a[7] = Kernel_A_gets_inverse_A_9;
    *ierr = 0;
}
