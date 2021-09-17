
/* */
#include <stdlib.h>
#include <stdio.h>
#include "petsc.h"
/* */

extern PetscErrorCode PetscKernel_A_gets_inverse_A_2(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_3(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_4(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_5(MatScalar *, PetscInt *, MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_6(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_7(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_9(MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscKernel_A_gets_inverse_A_15(MatScalar *, PetscInt *, MatScalar *, PetscReal, PetscBool, PetscBool *);
extern PetscErrorCode PetscLINPACKgefa(MatScalar*,PetscInt,PetscInt*, PetscBool, PetscBool *);
extern PetscErrorCode PetscLINPACKgedi(MatScalar*,PetscInt,PetscInt*,MatScalar*);

static PetscErrorCode (*a[8])(MatScalar *mat, PetscReal, PetscBool, PetscBool *);

int solven_(int *n, MatScalar *mat)
{
PetscBool      allowzeropivot=PETSC_FALSE,zeropivotdetected;
    switch(*n) {   
		case 5:	
		    {
			PetscInt ipvt[5];
			MatScalar work[5*5];
			return (int)PetscKernel_A_gets_inverse_A_5(mat, ipvt, work, (PetscReal)0.0, allowzeropivot,&zeropivotdetected);
			}
		case 15:
		    {
			PetscInt ipvt[15];
			MatScalar work[15*15];
			return (int)PetscKernel_A_gets_inverse_A_15(mat, ipvt, work, (PetscReal)0.0, allowzeropivot,&zeropivotdetected);
			}
		case 2:
		case 3:
		case 4:
		case 6:
		case 7:
		case 9:
			return (int)a[(*n-2)](mat, (PetscReal) 0.0, allowzeropivot,&zeropivotdetected);
		default:
			printf("Cannot deal with this bs");
			return *n;
/*
			{
			PetscInt ipvt[*n];
			MatScalar work[*n * *n];
                        int q = (int)PetscLINPACKgefa(mat, (PetscInt)*n, ipvt, allowzeropivot,&zeropivotdetected);
                        if (q) return q;                                                                                                                                                                           
                        return (int)PetscLINPACKgedi(mat, (PetscInt)*n, ipvt, work);
			}
*/
	}
}


void setsolven_(int *ierr)
{
    a[0] = PetscKernel_A_gets_inverse_A_2;
    a[1] = PetscKernel_A_gets_inverse_A_3;
    a[2] = PetscKernel_A_gets_inverse_A_4;
    a[3] = NULL;
    a[4] = PetscKernel_A_gets_inverse_A_6;
    a[5] = PetscKernel_A_gets_inverse_A_7;
	a[6] = NULL;
	a[7] = PetscKernel_A_gets_inverse_A_9;
    *ierr = 0;
}
