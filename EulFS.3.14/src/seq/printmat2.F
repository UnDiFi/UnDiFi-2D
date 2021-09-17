      SUBROUTINE PrintMatCSR(A,RHS,X,ia,ja,mat_array,nr,nz_a,ordsqr,
     &ITER,JOB)
C
C     $Id: printmat2.F,v 1.9 2020/04/23 09:37:34 abonfi Exp $
C
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
CCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscviewer.h"
              use petscpc
              use petscksp
C
      IMPLICIT NONE
#include "io.com"
C
      Mat A
      Vec RHS,X
      INTEGER ITER,JOB
      INTEGER nr,nz_a,nnz,ix
      integer ia(*),ja(*)
C     local vars
      double  precision info(MAT_INFO_SIZE)
!     double  precision dnrm2
!     double  precision smax,sl2
      integer idamax
      INTEGER jbgn,jend,ncols,i,j,k,iout,shift,iwrite,ordsqr
      PetscScalar mat_array(1)
      PetscOffset i_mat,i_ia,i_ja
      PetscErrorCode ierr
      PetscBool done,blockcompressed
      MatType type
      PARAMETER(blockcompressed=PETSC_TRUE)
!     PARAMETER(blockcompressed=PETSC_FALSE,ifmt=4)
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
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         nz_a = info(MAT_INFO_NZ_USED)
         ordsqr = info(MAT_INFO_BLOCK_SIZE)
!        CALL VecGetLocalSize(rhs,n,IFAIL)
!        maxcols = max(maxcols,n)
         call MatGetType(A,type,IFAIL)
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         write(6,*)'PE # ',MY_PE,' nnz = ',nz_a,' nr = ',nr,
     &' bs = ',ordsqr,' type = ',type
         RETURN
      ENDIF
C
C     IF( ITER. NE. 3 )RETURN
C
      WRITE(matfile(4:6),FMT="(I3.3)")MY_PE
!     WRITE(rhsfile(4:6),FMT="(I3.3)")0
!     WRITE(solfile(4:6),FMT="(I3.3)")0
      done = PETSC_FALSE
c     use Fortran based indexing
      shift = 1
      write(6,*)'blockcompressed = ',blockcompressed,
     &' TRUE = ',PETSC_TRUE
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
!     write(NOUT,*)(ia(i+i_ia),i=1,nr+1)
!     do i = 1,nr
!        jbgn = ia(i_ia+i)
!        jend = ia(i_ia+i+1)-1
!     write(NOUT,*)'row = ',i,(ja(j+i_ja),j=jbgn,jend)
!     enddo
c
!     call dump (1,nr,.FALSE.,mat_array(1),ja(i_ja+1),
!    &ia(i_ia+1),iout)
         call MatSeqAIJGetArray(A,mat_array,i_mat,IFAIL)
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         nz_a = nz_a/ordsqr
         if( nnz .NE. nz_a )then
             write(6,*)'There is a mismatch btw nnz,nz_a ',
     &       nnz,nz_a,ordsqr
             STOP
         else
             write(6,*)' nnz,nz_a match ',nnz,nz_a,ordsqr
         endif
         iwrite = 99
cnag     call x04abf(1,iwrite) ! removed in release 3.7.1 when the NAG library has been removed
         do i = 1,nr
         call printrow(mat_array(i_mat+1),ia(i_ia+1),ja(i_ja+1),i,
     &ordsqr,iwrite)
         enddo
         CALL EXIT(0)
         return
!        ix = idamax(nnz,mat_array(i_mat+1),1)
!        smax = mat_array(i_mat+ix)
!        sl2 = dnrm2(nnz,mat_array(i_mat+1),1)
!        write(NOUT,*)' Inf/L2 norms of ILU on PE # ',MY_PE,
!    &   'at Newton step ',ITER,' are : ',smax,sl2
!        call dump (1,nr,.TRUE.,mat_array(i_mat+1),ja(i_ja+1),
!    &ia(i_ia+1),iout)
         iout = 10+MY_PE
         open(iout,FILE=matfile)
!        call skit (nr,mat_array(i_mat+1),ja(i_ja+1),
!    &ia(i_ia+1),ifmt,iout,IFAIL)
         close(iout)
         call MatSeqAIJRestoreArray(A,mat_array,i_mat,ierr)
!     ENDIF
      CALL MatRestoreRowIJ(A,shift,PETSC_FALSE,blockcompressed,nr,
     &   ia,i_ia,ja,i_ja,done,IFAIL)
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
C
c
      RETURN
      END
      subroutine printrow(a,ia,ja,ir,ordsqr,iwrite)
      implicit none
      include 'paramt.h'
      double precision a(*) 
      integer ia(*),ja(*)
      integer ir,iwrite
      integer ordsqr
      integer jbgn,jend,iaddr,IFAIL,j,jcol,n,i
      character*25 label
      double precision deta
      double precision work(MAX_NOFVAR_SQR)
      integer ipiv(MAXNOFVAR)
      jbgn = ia(ir) 
      jend = ia(ir+1) -1
      n = int(sqrt(real(ordsqr)))
      write(iwrite,*)'*********************************************'
      write(iwrite,*)' block row i = ',ir
      write(iwrite,*)'*********************************************'
      do j = jbgn,jend
         jcol = ja(j)
         iaddr = (j-1)*ordsqr +1 
         write(label,FMT=100)ir,jcol
         CALL R8Mat_Print('General',' ',n,n,a(iaddr),
     +   n,label,IFAIL)
         if(ir.EQ.jcol)THEN
            deta = 1.d0
            CALL DCOPY(ordsqr,a(iaddr),1,work,1)
            CALL DGETRF( n, N, work, n, IPIV, IFAIL )
            IF(IFAIL.NE.0)THEN
                WRITE(iwrite,*)'DGETRF has ret IFAIL = ',IFAIL
                CALL EXIT(IFAIL)
            ELSE
                do i = 1,ordsqr,n+1
                   deta = deta * work(i) 
                enddo
            ENDIF
         CALL R8Mat_Print('General',' ',n,n,work,
     +   n,'LU factor',IFAIL)
            WRITE(iwrite,*)'det(A) = ',deta
         endif
      enddo
  100 FORMAT('row = ',I6,' col = ',I6)
      return
      end
