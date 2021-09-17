      SUBROUTINE AddDttoLhs(DT,VMEDIAN,DELT,NOFVAR,NPOIN)
C
C     $Id: JacobianBoundaryConditions.F,v 1.12 2004/02/27 15:46:46 abonfi Exp abonfi $
C
C     Add real time step in the lhs
C     
C
      IMPLICIT NONE
C
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
C
      Vec DT
      PetscScalar  x_array(1)
      PetscOffset i_x
      INTEGER NOFVAR,NPOIN
      INTEGER IFAIL
      DOUBLE PRECISION VMEDIAN(NPOIN),DELT
      INTEGER ioff,ivar,ipoin
      DOUBLE PRECISION temp
C
C
C
      INCLUDE 'constants'
      INCLUDE 'paramt.h'
C
      CALL VecGetArray(DT,x_array,i_x,IFAIL)
C
C     add the term dx/dt to the lhs
C
      DO 100 IPOIN = 1,NPOIN
         temp = vmedian(IPOIN)/DELT
!        write(6,*)ipoin,temp
C
C     Add the real time-step to the momentum eqns only 
C
         DO 200 IVAR = 1,NOFVAR
            ioff = (IPOIN-1)*NOFVAR+IVAR 
            x_array(i_x+ioff) = x_array(i_x+ioff)+temp
  200    CONTINUE 
  100 CONTINUE 
      CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
C
C
      RETURN
      END
