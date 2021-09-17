      SUBROUTINE PrintMatMM(A,RHS,X,ir,nr,nz_a,ordsqr,
     &ITER,JOB)
C
CCCC#define VIEWER_MATLAB
#define VIEWER_ASCII
#define ASCII
C
C     $Id: printmatmm.F,v 1.9 2020/04/23 09:48:01 abonfi Exp $
C
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscviewer.h"
      use petscpc
      use petscksp
C
      IMPLICIT NONE
C
#include "io.com"
C
      Mat A
      Vec RHS,X
      INTEGER ITER,JOB
      INTEGER nr,nz_a,nnz,ix
      integer ia(1),ja(1),ir(*)
C     local vars
      double  precision info(MAT_INFO_SIZE)
      double  precision dnrm2
!     double  precision smax,sl2
      integer ival(1)
      complex cval(1)
      INTEGER jbgn,jend,ncols,i,j,k,iout,shift,iwrite,ordsqr,jobid
      INTEGER rows,cols,entries,bs
      PetscScalar mat_array(1),x_array(1),s
      PetscOffset i_mat,i_ia,i_ja,i_x
      PetscErrorCode ierr
      PetscBool done,blockcompressed
      MatType type
      PARAMETER(blockcompressed=PETSC_TRUE) ! it looks like we can only use PETSC_TRUE when the matrix is BAIJ 
!     PARAMETER(blockcompressed=PETSC_FALSE,ifmt=4)
      PetscViewer MyOpenMindedViewer
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
      CHARACTER* 10 matfile,rhsfile,solfile
      character*10 rep
      character*7 field 
      character*19 symm
      PetscInt n
!     PetscErrorCode IFAIL
C     ..
C     .. External Functions ..
C     ..
      integer idamax
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
      DATA matfile,rhsfile,solfile/"matXXX.bin","rhsXXX.bin",
     &"solXXX.bin"/
C     ..
C     Executable statements 
C
      IF(JOB.EQ.0)THEN
         write(6,*)'Welcome inside printmat MM; JOB = ',JOB
         call MatGetInfo(A,MAT_LOCAL,info,IFAIL)
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         nz_a = info(MAT_INFO_NZ_USED)
         ordsqr = info(MAT_INFO_BLOCK_SIZE)
         call MatGetType(A,type,IFAIL)
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         write(6,*)'PE # ',MY_PE,' nnz = ',nz_a,' nr = ',nr,
     &' bs = ',ordsqr,' Mat type is ',type
         RETURN
      ELSE
         write(6,*)'PE # ',MY_PE,' nnz = ',nz_a,' nr = ',nr,
     &' bs = ',ordsqr
      ENDIF
      iwrite = 99
C
      IF( ITER. EQ. 1 .OR. ITER .EQ. 8 )THEN
          WRITE(6,*)'ITER = ',ITER,' JOB = ',JOB
      ELSE
          WRITE(6,*)'ITER = ',ITER,' JOB = ',JOB
          RETURN
      ENDIF
C
      WRITE(matfile(4:6),FMT="(I3.3)")ITER
      WRITE(rhsfile(4:6),FMT="(I3.3)")ITER
      WRITE(solfile(4:6),FMT="(I3.3)")ITER
C
#ifndef MPI
      done = PETSC_FALSE
c     use Fortran based indexing
      shift = 1
      write(6,*)'blockcompressed = ',blockcompressed,
     &' PETSC_TRUE = ',PETSC_TRUE
!      call MatGetRowIJ(mat,shift,symmetric,inodecompressed,n,ia,iia,ja,jja,done,ierr)
         CALL MatGetRowIJ(A,shift,PETSC_FALSE,blockcompressed,nr,
     &   ia,i_ia,ja,i_ja,done,IFAIL)
      IF((done.NEQV.PETSC_TRUE) .OR. (IFAIL.NE.0) )then
         write(6,*)' MatGetRowIJ has ret done = ',done,PETSC_FALSE
         write(6,*)' MatGetRowIJ has ret IFAIL = ',IFAIL
         STOP
      endif
!
      nnz = ia(i_ia+nr+1)-ia(i_ia+1)
      write(NOUT,*)' Nof rows = ',nr
      write(NOUT,*)' Nof non-zero block entries is = ',nnz
      IF(blockcompressed.EQV.PETSC_TRUE)nz_a = nz_a/ordsqr
      if( nnz .NE. nz_a )then
             write(6,*)'There is a mismatch btw nnz,nz_a ',
     &       nnz,nz_a,ordsqr
             STOP
      else
             write(6,*)' nnz and nz_a match ',nnz,nz_a,ordsqr
      endif
#endif
#if 0
caldo
      call MatSeqAIJGetArray(A,mat_array,i_mat,IFAIL)
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
      iout = 100
!     call dump (1,nr,.FALSE.,mat_array(i_mat+1),ja(i_ja+1),
!    &ia(i_ia+1),iout)
!     call exit(0)
#endif
#ifndef MPI
         jobid = 1
c
c job   = integer serving as a job indicator. 
c         if job = 1 fill in only the array ir, ignore jc, and ao.
c         if job = 2 fill in ir, and jc but not ao 
c         if job = 3 fill in everything.
c         The reason why these options are provided is that on return 
c         ao and jc are the same as a, ja. So when job = 3, a and ja are
c         simply copied into ao, jc.  When job=2, only jc and ir are
c         returned. With job=1 only the array ir is returned. Moreover,
c         the algorithm is in place:
c            call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
c         will write the output matrix in coordinate format on a, ja,ia.
c
         write(6,*)'Calling csrcoo with job = ',jobid
         call csrcoo (nr,jobid,nz_a,mat_array(i_mat+1),ja(i_ja+1),
     2                ia(i_ia+1),nnz,mat_array(i_mat+1),ir,ja(i_ja+1),
     3                ifail)
         write(6,*)'Done! '
c
c        writing pattern
c
         iwrite = 99
         open(iwrite,FILE="pattern.mm")
         rep = "coordinate"
         field = "real"
         field = "pattern"
         symm = "general"
         write(6,*)'Calling mmwrite'
         call mmwrite(iwrite,rep,field,symm,nr,nr,nnz,
     *                    ir,ja(i_ja+1),ival,mat_array(i_mat+1),cval)
!        template
!        call mmwrite(ounit,rep,field,symm,rows,cols,nnz,
!    *                    indx,jndx,ival,rval,cval)
         close(iwrite)
         write(6,*)'Done! '
         WRITE(6,*)'nr before MatRestoreRowIJ = ',nr
         CALL MatRestoreRowIJ(A,shift,PETSC_FALSE,blockcompressed,nr,
     &   ia,i_ia,ja,i_ja,done,IFAIL)
         WRITE(6,*)'nr after MatRestoreRowIJ = ',nr
         write(6,*)' MatRestoreRowIJ has ret IFAIL = ',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
#endif
#if 0
c
c        write matrix entries
c
#ifdef ASCII
         write(6,*)'Writing matrix to ASCII file'
         open(iwrite,FILE=matfile,FORM="FORMATTED")
         write(iwrite,*)(mat_array(i_mat+i),i=1,nz_a*ordsqr)
         close(iwrite)
#else
         write(6,*)'Writing matrix to binary file'
         open(iwrite,FILE=matfile,FORM="UNFORMATTED")
         write(iwrite)(mat_array(i_mat+i),i=1,nz_a*ordsqr)
         close(iwrite)
         write(6,*)'Done! '
#endif
#endif
#ifdef ASCII
C
C     This is where the MAtlab section beginENDs +++++++++++++++++++++++++++++++++++++++
C
         write(matfile(8:8),FMT="('m')")
         CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matrix.ascii",
     &                          MyOpenMindedViewer,IFAIL)
!        CALL PetscViewerASCIIOpen(PETSC_COMM_SELF,matfile(1:8),
!    &                          MyOpenMindedViewer,IFAIL)
         write(6,*)'PetscViewerASCIIOpen IFAIL = ',IFAIL,' on ',my_pe
#ifdef VIEWER_ASCII
         CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                  PETSC_VIEWER_ASCII_IMPL,IFAIL)
#else
C
!     CALL PetscViewerBinaryOpen(PETSC_COMM_SELF,matfile,  !"mat.bin",
!    &FILE_MODE_WRITE,MyOpenMindedViewer,IFAIL)
         CALL PetscViewerSetFormat(MyOpenMindedViewer,
     &                          PETSC_VIEWER_ASCII_MATLAB,IFAIL)
C
#endif
         write(6,*)'PetscViewerSetFormat IFAIL = ',IFAIL,' on ',my_pe
!
C
         write(6,*)"Start dumping matrix at iteration on proc. no. ",
     &my_pe
         CALL MatView(A,MyOpenMindedViewer,IFAIL)
!     CALL MatView(A,PETSC_VIEWER_STDOUT_SELF,IFAIL)
         write(6,*)"I just finished Dumping matrix on proc. no. ",
     &my_pe,' MatView says IFAIL = ',IFAIL
C
C
         CALL PetscViewerDestroy(MyOpenMindedViewer,IFAIL)
C
C     This is where the MAtlab section ENDs +++++++++++++++++++++++++++++++++++++++
C
#endif
#if 0
         call MatSeqAIJRestoreArray(A,mat_array,i_mat,ifail)
         write(6,*)' MatRestoreArray has ret ierr = ',ifail
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
#endif
#ifndef MPI
         open(iwrite,FILE="pattern.mm")
         write(6,*)'Reading info from ',"pattern.mm"
         write(6,*)'Calling mminfo'
         call mminfo(iwrite,rep,field,symm,rows,cols,entries) 
         close(iwrite)
         write(6,*)rep
         write(6,*)field
         write(6,*)symm
         write(6,*)'rows are ',rows,' cols are ',cols,' nnzr = ',
     &entries 
         close(iwrite)
#endif
C
C        writing the solution
C
         write(6,*)'Writing solution to binary file'
         call VecGetArray(x,x_array,i_x,ifail)
         call VecGetSize(x,rows,ifail)
         call VecGetBlockSize(x,bs,ifail)
         write(6,*)'Vec block size is = ',bs
         nr = rows/bs
         write(6,*)'Vec block size of X is = ',nr
#ifdef ASCII
         write(6,*)'Writing solution to ASCII file'
         open(iwrite,FILE=solfile,FORM="FORMATTED")
         write(iwrite,*)nr*bs
         do i = 1,nr*bs
         write(iwrite,*)i,x_array(i_x+i)
         enddo
         close(iwrite)
#else
         write(6,*)'Writing solution to binary file'
         open(iwrite,FILE=solfile,FORM="UNFORMATTED")
         write(iwrite)nr*bs
         write(iwrite)(x_array(i_x+i),i=1,nr*bs)
         close(iwrite)
#endif
         CALL VecNorm(X,NORM_2,S,IFAIL)
         write(6,*)'L-2 norm of X is ',s
         write(6,*)'Done! '
         call VecRestoreArray(x,x_array,i_x,ifail)
C
C        writing the rhs
C
         call VecGetArray(rhs,x_array,i_x,ifail)
#ifdef ASCII
         write(6,*)'Writing rhs to ASCII file'
         open(iwrite,FILE=rhsfile,FORM="FORMATTED")
         write(iwrite,*)nr*bs
         do i = 1,nr*bs
         write(iwrite,*)i,x_array(i_x+i)
         enddo
         close(iwrite)
#else
         write(6,*)'Writing rhs to binary file'
         open(iwrite,FILE=rhsfile,FORM="UNFORMATTED")
         write(iwrite)nr*bs
         write(iwrite)(x_array(i_x+i),i=1,nr*bs)
         close(iwrite)
#endif
         CALL VecNorm(rhs,NORM_2,S,IFAIL)
         write(6,*)'L-2 norm of RHS is ',s
         call VecRestoreArray(x,x_array,i_x,ifail)
         write(6,*)'Done! '
         write(6,*)'Now leaving printmatmm'
         call exit(0)
C
      RETURN
      END
