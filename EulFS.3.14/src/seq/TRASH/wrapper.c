#include "petsc.h"
#include "petscfix.h"
#include "HYPRE.h"
#include "IJ_mv.h"
/* Fortran interface file */


#ifdef PETSC_USE_POINTER_CONVERSION
#if defined(__cplusplus)
extern "C" { 
#endif 
extern void *PetscToPointer(void*);
extern int PetscFromPointer(void *);
extern void PetscRmPointer(void*);
#if defined(__cplusplus)
} 
#endif 

#else

#define PetscToPointer(a) (*(long *)(a))
#define PetscFromPointer(a) (long)(a)
#define PetscRmPointer(a)
#endif

#include "petscmat.h"
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mathypre_ijmatrixcreate_ MATHYPRE_IJMATRIXCREATE
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mathypre_ijmatrixcreate_ mathypre_ijmatrixcreate
#endif
#ifdef PETSC_HAVE_FORTRAN_CAPS
#define mathypre_ijmatrixcopy_ MATHYPRE_IJMATRIXCOPY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define mathypre_ijmatrixcopy_ mathypre_ijmatrixcopy
#endif


/* Definitions of Fortran Wrapper routines */
#if defined(__cplusplus)
extern "C" {
#endif
void PETSC_STDCALL   mathypre_ijmatrixcreate_(Mat v,HYPRE_IJMatrix *ij, int *__ierr ){
*__ierr = MatHYPRE_IJMatrixCreate( v, &ij);
}
void PETSC_STDCALL   mathypre_ijmatrixcopy_(Mat v,HYPRE_IJMatrix ij, int *__ierr ){
*__ierr = MatHYPRE_IJMatrixCopy( v, ij);
}
#if defined(__cplusplus)
}
#endif
