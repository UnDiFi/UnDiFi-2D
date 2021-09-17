#define PETSCMAT_DLL

#include "petscmat.h"
#include "src/mat/order/order.h"

EXTERN_C_BEGIN
/*
 *     MatOrdering_GPSK - Find the Reverse Cuthill-McKee ordering of a given matrix.
 *     */
#undef __FUNCT__
#define __FUNCT__ "MatOrdering_RCM"
PetscErrorCode PETSCMAT_DLLEXPORT MatOrdering_GPSK(Mat mat,const MatOrderingType type,IS *row,IS *col)
{
  int        wrklen,nnz,ierr,i,nrow,*work,*degree,*rstart,*connec,*ia,*ja,*perm;
/*
  *work,*degree, etc. sono puntatori ad interi 
*/    
  PetscTruth done;

  PetscFunctionBegin;
  ierr = MatGetRowIJ(mat,1,PETSC_TRUE,&nrow,&ia,&ja,&done);CHKERRQ(ierr);
  if (!done) SETERRQ(PETSC_ERR_SUP,"Cannot get rows for matrix");

  nnz = ia[nrow]-ia[0];
  wrklen = 6*nrow + 3;
  ierr = PetscMalloc((3*nrow+1+wrklen+nnz) * sizeof(int),&work);CHKERRQ(ierr);
  perm = work + wrklen;
  degree  = perm + nrow;
  rstart  = degree + nrow;
  connec  = rstart + nrow+1;

  printf("%d %d\n",nnz,nrow);

  exrcm_(&nrow,ia,ja,perm,degree,rstart,connec,work,&wrklen);
  ierr = MatRestoreRowIJ(mat,1,PETSC_TRUE,&nrow,&ia,&ja,&done);CHKERRQ(ierr);

  /* shift because Sparsepack indices start at one */
  for (i=0; i<nrow; i++) perm[i]--;

  ierr = ISCreateGeneral(PETSC_COMM_SELF,nrow,perm,row);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,nrow,perm,col);CHKERRQ(ierr);
  ierr = PetscFree(work);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
EXTERN_C_END
