      subroutine sillysub(A,X,VRATIO,IBGN,NDIM,NOFVAR,COMPRESSIBLE)
C
C     $Id: lhsbc0.F,v 1.9 2001/05/24 10:25:07 abonfi Exp abonfi $
C
C     set inflow b.c. for internal flows:
C     total enthalpy, entropy and flow angle in the compressible case
C     total pressure and flow angle in the INcompressible case
C
      implicit none
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscts.h"
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants'
      INCLUDE 'visco.com'
      INCLUDE 'stream.com'
C
#include "iset.com"
C
C  Input/output parameters:
      integer ndim,nofvar
      double precision x(*),vratio(3,*)
      Mat              A
      logical compressible
 
C  Local variables:
      PetscOffset i_is
      integer     is_array(1)
      integer maxcols
      parameter(maxcols=50)
      double precision amt(25)
      double precision cosb1,sinb1,kinetic,asqr,hh,denom,vals(maxcols)
      integer          i,IFAIL,nrows,irow(maxcols),jcol(maxcols),j
      integer          i1,i2,i3,i4,i5,ibgn
      integer          ipiv(maxcols)
      integer lda
      double precision pressc
      external pressc
C
      integer iadd
      iadd(i,j,lda) = (j-1)*lda + i
C
C
C     set diagonal element to 0.d0
C
      call ISGetIndices(Internal_inflow,is_array,i_is,IFAIL)
      call ISGetSize(Internal_inflow,Nrows,IFAIL)
C
C     outermost loop on inflow nodes
C
C     lda equals the no. of incoming characteristics, hence
C     the number of rows that need to be modified;
C
      j = 0
      IF(COMPRESSIBLE)THEN
C
          lda = ndim + 1
          DO 1 I = 1, Nrows, lda
              j = j+1
C
C     row (and column) pointers (0-based indexing) for the NOFVAR eqns.
C     irow(1) and irow(2) always address (r,rE) while
C     irow(3) and irow(4) can address any two of (ru,rv,rw)
C     jcol(*) addresses the whole vector of conserved variables
C     irow,jcol address rows and colums by global number,
C     while i1,..,i5 address local entries
C
              irow(1) = is_array(i_is + I+1)+1
              jcol(1) = irow(1)
              jcol(2) = jcol(1)+1
              jcol(3) = jcol(2)+1
              jcol(4) = jcol(3)+1
              jcol(5) = jcol(4)+1
C
C        1-based (FORTRAN) indexing and local numbering (-IBGN)
C
              i1 = is_array(i_is + I)
              i2 = i1 + nofvar -1
c
              call MatGetRow(A,irow(1),ncols,jcol,vals,IFAIL)
              do j = 1, ncols
                   if( jcol(j) .LT. i1 .OR. jcol(j) .GT. i2 )then
                       nc = nc + 1 
                       icol(j,nc) = jcol(j)
                   endif
              enddo
              call MatRestoreRow(A,irow(1),ncols,jcol,vals,IFAIL)
C
C     amt is:
C
C        r |rE|ru|rv|rw|
C       +---------------
C     r |  |  |  |  |  |
C    rE |  |  |  |  |  |
C    rv |  |  |  |  |  |
C    rw |  |  |  |  |  |
C
C
C
C     compute kinetic energy from the conserved variables
C
              kinetic = x(i3)*x(i3)+x(i4)*x(i4)
              if(ndim.eq.3)kinetic=kinetic+x(i5)*x(i5)
              kinetic = half*kinetic /(x(i1)*x(i1))
C
C        first row, entropy equation
C
              denom = 1.d0/(x(i1)**GAM)
              asqr = gam*gm1*(x(i2)/x(i1)-kinetic)
              amt(iadd(1,1,lda)) =(gm1*kinetic-asqr)*denom
              amt(iadd(1,2,lda)) =(gm1)*denom
              amt(iadd(1,3,lda)) =(-gm1*x(i3)/x(i1))*denom
              amt(iadd(1,4,lda)) =(-gm1*x(i4)/x(i1))*denom
              if(ndim.eq.3)amt(iadd(1,5,lda)) =(-gm1*x(i5)/x(i1))*denom
C
C        second row, energy equation
C
              hh = asqr / gm1 + kinetic
              amt(iadd(2,1,lda)) = (gm1*kinetic-hh)/x(i1)
              amt(iadd(2,2,lda)) = gam/x(i1)
              amt(iadd(2,3,lda)) =-gm1*x(i3)/x(i1)/x(i1)
              amt(iadd(2,4,lda)) =-gm1*x(i4)/x(i1)/x(i1)
              if(ndim.eq.3)amt(iadd(2,5,lda)) =-gm1*x(i5)/x(i1)/x(i1)
C
C        third or fourth row, flow angle
C
              amt(iadd(3,1,lda)) = (x(i4)*cosb1-x(i3)*sinb1)/(x(i1)*x
     +        (i1))
              amt(iadd(3,2,lda)) = 0.d0
              amt(iadd(3,3,lda)) = sinb1/x(i1)
              amt(iadd(3,4,lda)) = -cosb1/x(i1)
              if(ndim.eq.3)amt(iadd(3,5,lda)) = 0.d0
C
C        in 3D the b.c. is \Delta (w) = 0.
C
              if(ndim.eq.3)then
                  amt(iadd(4,1,lda)) = 0.d0
                  amt(iadd(4,2,lda)) = 0.d0
                  amt(iadd(4,3,lda)) = 0.d0
                  amt(iadd(4,4,lda)) = 0.d0
                  amt(iadd(4,5,lda)) = 1.d0/x(i1)
              endif
C
C        CALL X04CAF('General',' ',lda,NOFVAR,amt,
C    +   lda,'matrix for inflow b.c. ',IFAIL)
C        pause
C
              CALL MatSetValues(A,lda,irow,nofvar,jcol,amt,
     +        INSERT_VALUES, IFAIL)
    1     CONTINUE
C
      ELSE
C
C     Incompressible case
C
C
          lda = ndim
          DO 2 I = 1, Nrows, lda
              j = j+1
C
C     row (and column) pointers (0-based indexing) for the NOFVAR eqns.
C     irow(1) always address (p) while
C     irow(3) and irow(4) can address any two of (u,v,w)
C     jcol(*) addresses the whole vector of primitive variables
C             for the current meshpoint
C     irow,jcol address rows and colums by global number,
C     while i1,..,i4 address local entries
C
              irow(1) = is_array(i_is + I)
              irow(2) = is_array(i_is + I+1)
              IF(NDIM.EQ.3)irow(3) = is_array(i_is + I+2)
              jcol(1) = irow(1)
              jcol(2) = jcol(1)+1
              jcol(3) = jcol(2)+1
              jcol(4) = jcol(3)+1
C
C        1-based (FORTRAN) indexing and local numbering (-IBGN)
C
              i1 = jcol(1)+ 1-IBGN
              i2 = jcol(2)+ 1-IBGN
              i3 = jcol(3)+ 1-IBGN
              i4 = jcol(4)+ 1-IBGN
C
C        write(6,*)irow(1),i1,i2,i3,i4
C        write(6,*)x(i1)
C        write(6,*)x(i2)
C        write(6,*)x(i3)
C        write(6,*)x(i4)
C        pause
C
C     amt is:
C
C        dp|du|dv|dw|
C       +------------
C     p | 1| u| v| w|
C     u |  |  |  |  |
C     v |  |  |  |  |
C     w |  |  |  |  |
C
C
C
C     compute kinetic energy from the conserved variables
C
              kinetic = x(i2)*x(i2)+x(i3)*x(i3)
              if(ndim.eq.3)kinetic=kinetic+x(i4)*x(i4)
              kinetic = half*kinetic
C
C        first row, total pressure equation
C
              amt(iadd(1,1,lda)) = vratio(3,j)**2
c      write(6,*)j,vratio(3,j)**2
              amt(iadd(1,2,lda)) = x(i2)
              amt(iadd(1,3,lda)) = x(i3)
              if(ndim.eq.3)amt(iadd(1,4,lda)) = x(i4)
C
C        third or fourth row, flow angle
C
              amt(iadd(2,1,lda)) = zero
              amt(iadd(2,2,lda)) = sinb1
              amt(iadd(2,3,lda)) = -cosb1
              if(ndim.eq.3)amt(iadd(2,4,lda)) = 0.d0
C
C        in 3D the b.c. is \Delta (w) = 0.
C
              if(ndim.eq.3)then
                  amt(iadd(3,1,lda)) = 0.d0
                  amt(iadd(3,2,lda)) = 0.d0
                  amt(iadd(3,3,lda)) = 0.d0
                  amt(iadd(3,4,lda)) = 1.d0
              endif
C    +   lda,'matrix for inflow b.c. ',IFAIL)
C        pause
C
              CALL MatSetValues(A,lda,irow,nofvar,jcol,amt,
     +        INSERT_VALUES, IFAIL)
C
C        CALL X04CAF('General',' ',lda,NOFVAR,amt,
C    +   lda,'matrix for inflow b.c. ',IFAIL)
C        CALL DGETRF(3,3,amt,3,ipiv,ifail)
C if( ifail .GT. 0 )then
C            write(6,*)i,IFAIL
C     pause
C        endif
    2     CONTINUE
C
      ENDIF
C
      CALL ISRestoreIndices(Internal_inflow,is_array,i_is,IFAIL)
C
      RETURN
      END
