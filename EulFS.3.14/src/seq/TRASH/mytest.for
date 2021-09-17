      subroutine Mytest(X,B)
C
C     $Id: update3.F,v 1.24 2000/07/21 15:49:03 aldo Exp aldo $
C
C     set b.c. for adiabatic wall nodes in the Jacobian
C
      implicit none
C
#include "include/finclude/petsc.h"
#include "include/finclude/vec.h"
#include "include/finclude/mat.h"
#include "include/finclude/sles.h"
#include "include/finclude/ksp.h"
C#include "include/finclude/pc.h"
#include "include/finclude/is.h"
#include "include/finclude/ts.h"
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
      double precision t,dtmin
      integer     is_array(1)
      PetscOffset i_is

C  Local variables:
      integer maxcols
      parameter(maxcols=50)
      double precision alpha,y(maxcols),x(*),b(*)
      integer          i,IFAIL,nrows,ncols,irow(maxcols),icol(maxcols),j
      PetscOffset      i_x
      Scalar           x_array(1)
      PLogDouble       tbeg,tend
      integer cols(maxcols)
      double precision vals(maxcols),pressc
      external pressc
C
C
C
      ALPHA = TWALL/GAM/GM1/(M_infty*M_infty)
C
C     diagonal element corresponding to the energy eqn.
C
C
      call ISGetIndices(CnstPressure,is_array,i_is,IFAIL)
      call ISGetSize(CnstPressure,Nrows,IFAIL)
C
      DO 1 I = 1, Nrows
         irow(2) = is_array(i_is + I)+1
         irow(1) = irow(2)-1
         write(6,*).4*(x(irow(2))-0.5*(x(irow(2)+1)**2+
     &                                 x(irow(2)+2)**2))
c        write(6,*)X(irow(1)),X(irow(2))
    1 CONTINUE 
      CALL ISRestoreIndices(CnstPressure,is_array,i_is,IFAIL)
C
      RETURN
      END
