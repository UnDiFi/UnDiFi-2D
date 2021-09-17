      SUBROUTINE JacobianBoundaryConditions(ts,x,A,NDIM,NOFVAR)
C
      IMPLICIT NONE
C
C     N.B.: upon entry X contains the vector of CONSERVED
C           variables
C
C
C     $Id: JacobianBoundaryConditions.F,v 1.6 2001/05/24 10:24:09 abonfi Exp abonfi $
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscts.h"
C
      Mat A
      Vec X
      TS ts
      PLogDouble telapsed,tbeg,tend
      Scalar      x_array(1)
      integer      is_array(1)
      integer      iarray(30,30)
      PetscOffset i_x
      PetscOffset i_is
      integer maxcols
      parameter(maxcols=130)
      double precision vals(maxcols)
      double precision amt(5,maxcols)
      integer nmt(5,maxcols),ipiv(maxcols)
      INTEGER i,j,k,icol(maxcols),irow(5),lda,nrows,ncols(maxcols),i1
      INTEGER ic,jc
      INTEGER NDIM,NOFVAR
      INTEGER IFAIL,IBGN,IEND
      DOUBLE PRECISION DSTAK(1)
      COMMON/CSTAK/DSTAK
C     ..
C
#include "iset.com"
C
C
      INCLUDE 'constants'
      INCLUDE 'paramt.h'
      INCLUDE 'visco.com'
      INCLUDE 'flags.com'
      INCLUDE 'nloc'
C
      CALL PetscGetTime(tbeg,IFAIL)
C
C     the stiffness matrix has to be modified
C     to account for boundary conditions .. 
C
C     Note that the Index Sets address rows by GLOBAL number
C
      CALL MatZeroRows(A,SupersonicVariables,ONE,IFAIL)
      CALL MatZeroRows(A,NoSlipVelocities,ONE,IFAIL)
      CALL MatZeroRows(A,Internal_inflow,ONE,IFAIL)
      CALL MatZeroRows(A,CnstPressure,ONE,IFAIL)
      CALL MatZeroRows(A,Isothermal,ONE,IFAIL)
caldo
caldo CALL MatZeroRows(A,IsDebug,ONE,IFAIL)
caldo
C
      call VecGetOwnerShipRange(x,ibgn,iend,ifail)
      call VecGetArray(x,x_array,i_x,ifail)
      IF(IADIA.NE.0)
     &CALL Isothwall(A,x_array(i_x+1),IBGN,NDIM,NOFVAR,(ABS(KAN).EQ.4))
      CALL    lhsbc2(A,x_array(i_x+1),IBGN,NDIM,NOFVAR,(ABS(KAN).EQ.4))
      CALL    lhsbc0(A,x_array(i_x+1),DSTAK(LPTOT),
     +               IBGN,NDIM,NOFVAR,(ABS(KAN).EQ.4))
      call VecRestoreArray(X,x_array,i_x,ifail)
      CALL PetscGetTime(tend,IFAIL)
      telapsed=tend-tbeg
C
C     Final assembly
C
      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL)
C
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL)
C
C
C     some debugging stuff follows
C
#if 1
C
      call ISGetIndices(Internal_inflow,is_array,i_is,IFAIL)
      call ISGetSize(Internal_inflow,Nrows,IFAIL)
C
C     outermost loop on inflow nodes
C
C     lda equals the no. of incoming characteristics, hence
C     the number of rows that need to be modified
C     once again, this is currently implemented ONLY for compressible flows
C
      lda = ndim + 1
      ic=0
      DO 1 I = 1, Nrows, lda
      ic=ic+1
C
C     row (and column) pointers (0-based indexing) for the NOFVAR eqns.
C
         irow(1) = is_array(i_is + I) -IBGN
caldo    irow(2) = is_array(i_is + I+1) -IBGN
caldo    irow(3) = is_array(i_is + I+2) -IBGN
caldo    IF(NDIM.EQ.3)irow(4) = is_array(i_is + I+3) -IBGN
         irow(2) = irow(1)+1
         irow(3) = irow(2)+1
         irow(4) = irow(3)+1
         irow(5) = irow(4)+1
C
C        1-based (FORTRAN) indexing
C
C
C        write(6,*)irow(1),i1,i2,i3,i4
C        write(6,*)x(i1)
C        write(6,*)x(i2)
C        write(6,*)x(i3)
C        write(6,*)x(i4)
C        pause
C
         jc=2
         do j = 1,nofvar
         call MatGetRow(A,irow(j),ncols(j),icol,vals,IFAIL)
         call dcopy(ncols(j),vals,1,amt(j,1),5)
         call icopy(ncols(j),icol,1,nmt(j,1),5)
         if(j.eq.3)then
              iarray(1,ic) = irow(j)
              do k = 1, ncols(j)
                if( icol(k) .LT. irow(1) .OR. icol(k) .GT. irow(4) )then
                       jc = jc+1
                       iarray(jc,ic) = icol(k) 
                endif
              enddo
              iarray(2,ic) = jc-2
         endif
         call MatRestoreRow(A,irow(j),ncols(j),icol,vals,IFAIL)
         enddo 
         write(6,*)'Row numbers ',(irow(j),j=1,nofvar)
         write(6,*)'Col. number ',(ncols(j),j=1,nofvar)
         CALL X04EAF('General',' ',NOFVAR,ncols(1),nmt,
     +            5,'matrix for inflow b.c. ',IFAIL)
         CALL X04CAF('General',' ',NOFVAR,ncols(1),amt,
     +            5,'matrix for inflow b.c. ',IFAIL)
         do j = 1, ncols(1)
              if( nmt(1,j) .EQ. irow(1) )then
                    i1 = j
                    goto 33
              endif
         enddo
         stop ' should not get here'
   33    CALL DGETRF(nofvar,nofvar,amt(1,i1),5,ipiv,ifail)
         if( ifail .GT. 0 )then
             write(6,*)i,IFAIL
             pause
         endif
c        pause
C
    1 CONTINUE 
C
      CALL ISRestoreIndices(Internal_inflow,is_array,i_is,IFAIL)
#if 0
caldo
caldo
caldo
      call dinit(maxcols,0.d0,vals,1)
#if 0
      CALL X04EAF('General',' ',22,ic,iarray,
     +            30,'matrix for inflow b.c. ',IFAIL)
#endif
      DO I = 1, ic
caldo write(6,*)i,iarray(1,i),iarray(2,I)
              CALL MatSetValues(A,1,iarray(1,i),iarray(2,i),iarray(3,i),
     +        vals,INSERT_VALUES, IFAIL)
      END DO
C
C     Final assembly
C
      CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL)
C
      CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL)
C
#if 0
      call ISGetIndices(Internal_inflow,is_array,i_is,IFAIL)
      call ISGetSize(Internal_inflow,Nrows,IFAIL)
      DO 2 I = 1, Nrows, lda
C
C     row (and column) pointers (0-based indexing) for the NOFVAR eqns.
C
         irow(1) = is_array(i_is + I) -IBGN
         irow(2) = irow(1)+1
         irow(3) = irow(2)+1
         irow(4) = irow(3)+1
         irow(5) = irow(4)+1
C
         do j = 1,nofvar
         call MatGetRow(A,irow(j),ncols(j),icol,vals,IFAIL)
         call dcopy(ncols(j),vals,1,amt(j,1),5)
         call icopy(ncols(j),icol,1,nmt(j,1),5)
         call MatRestoreRow(A,irow(j),ncols(j),icol,vals,IFAIL)
         enddo 
         write(6,*)'Row numbers ',(irow(j),j=1,nofvar)
         write(6,*)'Col. number ',(ncols(j),j=1,nofvar)
         CALL X04EAF('General',' ',NOFVAR,ncols(1),nmt,
     +            5,'matrix for inflow b.c. ',IFAIL)
         CALL X04CAF('General',' ',NOFVAR,ncols(1),amt,
     +            5,'matrix for inflow b.c. ',IFAIL)
         pause
C
    2 CONTINUE 
      CALL ISRestoreIndices(Internal_inflow,is_array,i_is,IFAIL)
#endif
C
#endif
#endif
C
      RETURN
      END
