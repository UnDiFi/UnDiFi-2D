
static char help[] = "Reads a PETSc matrix and vector from a file and reorders it.\n\
  -f0 <input_file> : first file to load (small system)\n\
  -f1 <input_file> : second file to load (larger system)\n\n";

/*T
   Concepts: Mat^ordering a matrix - loading a binary matrix and vector;
   Concepts: Mat^loading a binary matrix and vector;
   Concepts: Vectors^loading a binary vector;
   Concepts: PetscLog^preloading executable
   Processors: 1
T*/

/* 
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h    - vectors
     petscsys.h    - system routines       petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers               
*/
#include "petscmat.h"
/* 
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h    - vectors
     petscsys.h    - system routines       petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers               
*/
#include "src/mat/matimpl.h"
/*#include "petscmat.h"*/

EXTERN_C_BEGIN
EXTERN int MatOrdering_GPSK(Mat,MatOrderingType,IS*,IS*);
EXTERN_C_END


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Mat                   A;                /* matrix */
  PetscViewer           fd;               /* viewer */
  char                  file[2][PETSC_MAX_PATH_LEN];     /* input file name */
  IS                    isrow,iscol;      /* row and column permutations */
  PetscErrorCode        ierr;
  const MatOrderingType rtype = MATORDERING_RCM;
  const MatOrderingType mytype = "my_order";
  PetscTruth            flg,PreLoad = PETSC_FALSE;

  PetscInitialize(&argc,&args,(char *)0,help);


  /* 
     Determine files from which we read the two linear systems
     (matrix and right-hand-side vector).
  */
  ierr = PetscOptionsGetString(PETSC_NULL,"-f0",file[0],PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate binary file with the -f0 option");
  ierr = PetscOptionsGetString(PETSC_NULL,"-f1",file[1],PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) PreLoad = PETSC_TRUE;
  MatOrderingRegisterDynamic(mytype,"/home/aldo/src/petsc-2.3.0/lib/linux_gnu/libgpsk.a","MatOrdering_GPSK",MatOrdering_GPSK);CHKERRQ(ierr);

  /* -----------------------------------------------------------
                  Beginning of loop
     ----------------------------------------------------------- */
  /* 
     Loop through the reordering 2 times.  
      - The intention here is to preload and solve a small system;
        then load another (larger) system and solve it as well.
        This process preloads the instructions with the smaller
        system so that more accurate performance monitoring (via
        -log_summary) can be done with the larger one (that actually
        is the system of interest). 
  */
  PreLoadBegin(PreLoad,"Load");

    /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
                           Load system i
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /* 
       Open binary file.  Note that we use PETSC_FILE_RDONLY to indicate
       reading from this file.
    */
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[PreLoadIt],PETSC_FILE_RDONLY,&fd);CHKERRQ(ierr);

    /*
       Load the matrix; then destroy the viewer.
    */
    ierr = MatLoad(fd,MATSEQAIJ,&A);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);


    /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    PreLoadStage("Reordering");
    /* - - - - - - - - - - - New Stage - - - - - - - - - - - - -
    ierr = MatGetOrdering(A,rtype,&isrow,&iscol);CHKERRQ(ierr);
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = MatGetOrdering(A,mytype,&isrow,&iscol);CHKERRQ(ierr);

    /* 
       Free work space.  All PETSc objects should be destroyed when they
       are no longer needed.
    */
    ierr = MatDestroy(A);CHKERRQ(ierr);
    ierr = ISDestroy(isrow);CHKERRQ(ierr);
    ierr = ISDestroy(iscol);CHKERRQ(ierr);
  PreLoadEnd();

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

