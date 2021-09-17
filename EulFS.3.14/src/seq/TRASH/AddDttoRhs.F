      SUBROUTINE AddDttoRhs(Rhs,ZRoe,ZOld,VMEDIAN,DELT,NOFVAR,NPOIN)
C
C     $Id: JacobianBoundaryConditions.F,v 1.12 2004/02/27 15:46:46 abonfi Exp abonfi $
C
C
      IMPLICIT NONE
C
C     Add the time derivative term in the rhs
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
C
      Vec Rhs
      PetscScalar  x_array(1),y_array(1),z_array(1)
      PetscOffset i_x,i_y,i_z
      INTEGER NOFVAR,NPOIN
      INTEGER IFAIL
      DOUBLE PRECISION DELT
      DOUBLE PRECISION VMEDIAN(NPOIN),ZROE(NOFVAR,*),ZOLD(NOFVAR,*)
      INTEGER ioff,ivar,ipoin
      DOUBLE PRECISION temp
C
C
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'flags.com'
C
      CALL VecGetArray(Rhs ,z_array,i_z,IFAIL)
C
C     subtract the term (u^{n+1,m}-u^n)*dx/dt from the rhs
C
      DO 100 IPOIN = 1,NPOIN
         temp = -vmedian(IPOIN)/DELT
         DO 200 IVAR = 1,NOFVAR
            ioff = (IPOIN-1)*NOFVAR+IVAR 
            z_array(i_z+ioff) = z_array(i_z+ioff)+
     &      (ZROE(IVAR,IPOIN)-ZOLD(IVAR,IPOIN)) * temp
!    &      (x_array(i_x+ioff)-y_array(i_y+ioff))
!           write(6,*)'rhs = ',ioff,z_array(i_z+ioff)
  200    CONTINUE 
  100 CONTINUE 
      CALL VecRestoreArray(Rhs ,z_array,i_z,IFAIL)
C
C
      RETURN
      END
