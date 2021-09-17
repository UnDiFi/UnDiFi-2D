      SUBROUTINE PrintMat(A,RHS,X,cols,values,maxcols,ITER,JOB)
C
CC#define HYPER_PARCSR 5555
C
CCCC#define VIEWER_ASCII
CC#define HAVE_HYPRE
#define VIEWER_BINARY
C
C
C
C     $Id: printmat.F,v 1.9 2020/04/23 09:37:34 abonfi Exp $
C
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
CCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscviewer.h"
              use petscksp
              use petscpc
C
      IMPLICIT NONE
C
#ifdef HAVE_HYPRE
CCC#include "HYPRE.h"
CCC#include "IJ_mv.h"
#endif

C
      Mat A
      Vec RHS,X
      INTEGER ITER,JOB
      INTEGER rstart,rend,ncols,i,k
      INTEGER HYPRE_PARCSR
      parameter(HYPRE_PARCSR=5555)
      INTEGER maxcols
      integer cols(maxcols)
      double precision values(maxcols)
      double  precision info(MAT_INFO_SIZE)
      INTEGER nz_a
C
C
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,MY_PE
C     ..
C     .. Local Arrays ..
      PetscViewer MyOpenMindedViewer
      INTEGER*8 ij
      CHARACTER* 6 matfile,rhsfile,solfile
      PetscScalar x_array(1)
      PetscOffset i_x
      PetscInt n
!     PetscErrorCode IFAIL

C     ..
C     .. External Functions ..
C     ..
C     .. External Subroutines ..
C     ..
C
C     .. Intrinsic Functions ..
C     ..
C     .. Common blocks ..
      COMMON /MPICOM/MY_PE
C     ..
C     .. Equivalences ..
C     ..
C     .. Data statements ..
      DATA matfile,rhsfile,solfile/"matXXX","rhsXXX","solXXX"/
C     ..
C     Executable statements 
C
      IF(JOB.EQ.0)THEN
         call MatGetInfo(A,MAT_LOCAL,info,IFAIL)
         nz_a = info(MAT_INFO_NZ_ALLOCATED)
!        maxcols = info(MAT_INFO_COLUMNS_LOCAL)
         maxcols = 0
         CALL VecGetLocalSize(rhs,n,IFAIL)
         maxcols = max(maxcols,n)
         write(6,*)'PE # ',MY_PE,nz_a,maxcols,n
         RETURN
      ENDIF
C
      IF( ITER. NE. 11 )RETURN
C
      WRITE(matfile(4:6),FMT="(I3.3)")0
      WRITE(rhsfile(4:6),FMT="(I3.3)")0
      WRITE(solfile(4:6),FMT="(I3.3)")0
!
!    #     #  #     #  ######   ######   #######  
!    #     #   #   #   #     #  #     #  #        
!    #     #    # #    #     #  #     #  #        
!    #######     #     ######   ######   #####    
!    #     #     #     #        #   #    #        
!    #     #     #     #        #    #   #        
!    #     #     #     #        #     #  #######  
!    
#ifdef HAVE_HYPRE
      CALL MatGetOwnershipRange(A,rstart,rend,IFAIL)
      CALL HYPRE_IJMatrixCreate(PETSC_COMM_WORLD,rstart,rend-1,rstart,
     &rend-1,ij,IFAIL)
      write(6,*)'HYPRE_IJMatrixCreate has returned ifail = ',ifail
      CALL HYPRE_IJMatrixSetObjectType(ij,HYPRE_PARCSR,IFAIL) 
      write(6,*)'HYPRE_IJMatrixSetObjectType has returned ifail = ',
     &ifail
      CALL  HYPRE_IJMatrixInitialize(ij,IFAIL)
      write(6,*)'HYPRE_IJMatrixInitialize has returned ifail = ',
     &ifail
      do i= rstart, rend-1
        CALL MatGetRow(A,i,ncols,cols,values,IFAIL)
        if(ncols.GT.maxcols)STOP 'Must increase maxcols'
        CALL HYPRE_IJMatrixSetValues(ij,1,ncols,i,cols,values,IFAIL)
        CALL MatRestoreRow(A,i,ncols,cols,values,IFAIL)
      enddo
      CALL HYPRE_IJMatrixAssemble(ij,IFAIL)
      CALL HYPRE_IJMatrixPrint(ij,matfile,IFAIL)
      write(6,*)'HYPRE_IJMatrixPrint has returned ifail = ',
     &ifail
      CALL HYPRE_IJMatrixDestroy(ij,IFAIL)
      write(6,*)'HYPRE_IJMatrixDestroy has returned ifail = ',
     &ifail
C
C     now dump rhs
C
      CALL VecGetOwnerShipRange(rhs,rstart,rend,IFAIL)
      CALL HYPRE_IJVectorCreate(PETSC_COMM_WORLD,rstart,rend-1,
     &ij,IFAIL)
      write(6,*)'HYPRE_IJVectorCreate has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorSetObjectType(ij,HYPRE_PARCSR,IFAIL)
      write(6,*)'HYPRE_IJSetObjectType has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorInitialize(ij,IFAIL)
      write(6,*)'HYPRE_IJVectorInitialize has returned ifail = ',
     &ifail
      CALL VecGetLocalSize(rhs,n,IFAIL)
      if(n.GT.maxcols)STOP 'Must increase maxcols'
      CALL VecGetArray(rhs,x_array,i_x,IFAIL)
      k = 0 
      do i = rstart,rend-1
         k = k + 1
         cols(k) = i
      enddo
!     write(6,*)rend-rstart,n
      CALL HYPRE_IJVectorSetValues(ij,n,cols,x_array(i_x + 1),IFAIL)
      write(6,*)'HYPRE_IJVectorSetValues has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorAssemble(ij,IFAIL)
      CALL VecRestoreArray(rhs,x_array,i_x,IFAIL)
      CALL HYPRE_IJVectorPrint(ij,rhsfile,IFAIL)
      write(6,*)'HYPRE_IJVectorPrint has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorDestroy(ij,IFAIL)
c
c     dump solution
c
      CALL VecGetOwnerShipRange(x,rstart,rend,IFAIL)
      CALL HYPRE_IJVectorCreate(PETSC_COMM_WORLD,rstart,rend-1,
     &ij,IFAIL)
      write(6,*)'HYPRE_IJVectorCreate has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorSetObjectType(ij,HYPRE_PARCSR,IFAIL)
      write(6,*)'HYPRE_IJSetObjectType has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorInitialize(ij,IFAIL)
      write(6,*)'HYPRE_IJVectorInitialize has returned ifail = ',
     &ifail
      CALL VecGetLocalSize(x,n,IFAIL)
      if(n.GT.maxcols)STOP 'Must increase maxcols'
      CALL VecGetArray(x,x_array,i_x,IFAIL)
      k = 0 
      do i = rstart,rend-1
         k = k + 1
         cols(k) = i
      enddo
      CALL HYPRE_IJVectorSetValues(ij,n,cols,x_array(i_x + 1),IFAIL)
      write(6,*)'HYPRE_IJVectorSetValues has returned ifail = ',
     &ifail
      CALL VecRestoreArray(x,x_array,i_x,IFAIL)
      CALL HYPRE_IJVectorAssemble(ij,IFAIL)
      CALL HYPRE_IJVectorPrint(ij,solfile,IFAIL)
      write(6,*)'HYPRE_IJVectorPrint has returned ifail = ',
     &ifail
      CALL HYPRE_IJVectorDestroy(ij,IFAIL)

      RETURN
#endif
C
C     This is the non-HYPRE portion
C
#ifdef VIEWER_ASCII
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,matfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                  PETSC_VIEWER_ASCII_DEFAULT,IFAIL)
#endif
#ifdef VIEWER_MATLAB
C
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,matfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_MATLAB,IFAIL)
C
#endif
#ifdef VIEWER_BINARY
C
      CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,matfile//'.bin',
     &                          FILE_MODE_WRITE,MyOpenMindedViewer,
     &                          IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_DEFAULT,IFAIL)
C
#endif
!
C
      write(6,*)"Dumping matrix at iteration on proc. no. ",my_pe
      CALL MatView(A,MyOpenMindedViewer,IFAIL)
      write(6,*)"I just finished Dumping matrix on proc. no. ",
     &my_pe,' IFAIL = ',IFAIL
C
      CALL PetscViewerDestroy(MyOpenMindedViewer,IFAIL)
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
#ifdef VIEWER_ASCII
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,rhsfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_DEFAULT,IFAIL)
#endif
#ifdef VIEWER_MATLAB
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,rhsfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_MATLAB,IFAIL)
#endif
#ifdef VIEWER_BINARY
C
      CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,rhsfile//'.bin',
     &                          FILE_MODE_WRITE,MyOpenMindedViewer,
     &                          IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_DEFAULT,IFAIL)
C
#endif
C
      write(6,*)"Dumping rhs on proc. no. ",my_pe
      CALL VecView(RHS,MyOpenMindedViewer,IFAIL)
      write(6,*)"I just finished Dumping the rhs on proc. no. ",
     &my_pe,' IFAIL = ',IFAIL
C
      CALL PetscViewerDestroy(MyOpenMindedViewer,IFAIL)
C
C     dump solution
C
#ifdef VIEWER_BINARY
C
      CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,solfile//'.bin',
     &                          FILE_MODE_WRITE,MyOpenMindedViewer,
     &                          IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_DEFAULT,IFAIL)
C
#endif
#ifdef VIEWER_ASCII
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,solfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_DEFAULT,IFAIL)
#endif
#ifdef VIEWER_MATLAB
      CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,solfile,
     &                          MyOpenMindedViewer,IFAIL)
      CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_MATLAB,IFAIL)
#endif
C
      write(6,*)"Dumping solution on proc. no. ",my_pe
      CALL VecView(X,MyOpenMindedViewer,IFAIL)
      write(6,*)"I just finished with the solution vector on proc. no. "
     &,my_pe
C
      CALL PetscViewerDestroy(MyOpenMindedViewer,IFAIL)
C
      RETURN
      END
