      SUBROUTINE PREFATT(A,N,NNZ)
C
      IMPLICIT NONE
C
C     This routine computes the blocks of the matrix A with an overlapping
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscsles.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/pcpetsc.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscviewer.h"
C
      Mat A
      Mat B
      PetscTruth done 
C
      IS subind

C
C
C
C     .. Parameters ..
C      INTEGER NDNM
C      PARAMETER (NDNM=3*NMAX)
C     ..
C     .. Scalar Arguments ..
C      INTEGER NDIM,NOFVAR,NPOIN
      INTEGER NNZ
C     ..
C     .. Arrays in Common ..
C      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
C      DOUBLE PRECISION CNST,s,norm
      INTEGER N,IFAIL,MY_PE,ov,overlap
      INTEGER nrows,NNZ,lrow,lcol,gfirst,glastpo,i,ii,ll,rmin,rmax,ir,ic
      INTEGER nc,k,kk
C     ..
C     .. Local Arrays ..
C      DOUBLE PRECISION WKSP1(3,NMAX),WKSP2(3,NMAX)
C      INTEGER ISTAK(1)
      INTEGER IA(N+1),JA(NNZ+1),list(NNZ+1)
C      INTEGER IA(0:N),JA(0:NNZ),list(0:NNZ)
C     ..
C     .. External Functions ..
C      DOUBLE PRECISION DNRM2
C      INTEGER ISTKGT,ISTKST
C      EXTERNAL DNRM2,ISTKGT,ISTKST
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
C      INTRINSIC DLOG10,MIN
C     ..
C     .. Common blocks ..
C      COMMON /CSTAK/DSTAK
      COMMON /MPICOM/MY_PE
C     ..
C     .. Equivalences ..
C      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
C     .. Data statements ..
C      DATA WKSP1,WKSP2,ITS/NDNM*ZERO,NDNM*ZERO,0/
C     ..

C

      overlap = 1
 
      CALL MatGetLocalSize(A,lrow,lcol,IFAIL)
      write(6,*)'IFAIL=',IFAIL
      write(6,*)'lrow=',lrow
      write(6,*)'lcol=',lcol
      CALL MatGetOwnershipRange(A,gfirst,glastpo,IFAIL)
      write(6,*)'IFAIL=',IFAIL
      write(6,*)'gfirst=',gfirst
      write(6,*)'glastpo=',glastpo
      CALL MatGetRowIJ(A,0,PETSC_FALSE,nrows,IA,JA,done,IFAIL)
      write(6,*)'IFAIL=',IFAIL
      DO i = 1,N+1
          write(6,*)'IA(',i,')=',IA(i)
      END DO
      DO i = 1,NNZ
          write(6,*)'JA(',i,')=',JA(i)
      END DO

      DO i = 1,lrow
          list(i) = gfirst+i
      END DO
      ll = lrow
      write(6,*)'ll=',ll
      DO ov = 1,overlap
          rmin = list(1)
          write(6,*)'rmin=',rmin
          rmax = list(ll)
          write(6,*)'rmax=',rmax
          DO i = 1,lrow
 	      ir = gfirst+i
              ic = JA(IA(ir))
              nc = IA(ir+1)-IA(ir)
              k = 1
              DO WHILE (k.LE.nc)
                  IF (ic.LT.rmin) THEN
                      DO ii = ll,1,-1
                          list(ii+1) = list(ii)
                      END DO
                      list(1) = ic
                      ll = ll+1
                  ELSE IF (ic.GT.rmax) THEN
                      ll = ll+1
                      list(ll) = ic
                  ELSE
                      ii = 1
                      DO WHILE ((ic.GT.list(ii)).AND.(ii.LT.ll))
                          IF (ic.LT.list(ii+1)) THEN
                              DO kk = ll,ii+1,-1
                                  list(kk+1) = list(kk)
                              END DO
                              list(ii+1) = ic
                              ll = ll+1
                          ELSE
                              ii = ii+1
                          ENDIF
                      END DO
                  ENDIF
                  ic = JA(IA(ir)+k)
                  k = k+1
              END DO 
          END DO
      END DO

      CALL MatRestoreRowIJ(A,1,PETSC_FALSE,nrows,IA,JA,done,IFAIL)
      write(6,*)'IFAIL=',IFAIL

      CALL ISCreateGeneral(MPI_COMM_SELF,ll,list,subind,IFAIL)
      write(6,*)'IFAIL=',IFAIL
      CALL ISSorted(subind,IFAIL)
      write(6,*) 'list is sorted? IFAIL=',IFAIL

      CALL MatGetSubmatrices(A,1,subind,subind,MAT_INITIAL_MATRIX,B)
      write(6,*)'IFAIL=',IFAIL




      END




