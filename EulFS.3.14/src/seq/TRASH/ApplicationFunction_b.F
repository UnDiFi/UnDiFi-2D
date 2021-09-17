      SUBROUTINE ApplicationFunction_b(PetscMatC,X,PetscVecRHS,
     &                                 ctx,IFAIL)
C ---------------------------------------------------------------------
C
C  ApplicationFunction_b - Computes nonlinear function
C  for the Motion solver
C  This function is called by
C  the higher level routines RHSFunction_b() and RHSJacobian_b()
C
C  $Id: ApplicationFunction_t.F,v 1.9 2011/09/16 07:19:48 abonfi Exp $
C
C  Input Parameter:
C  x - local vector data as a FORTRAN array
C  A - the jacobian matrix
C
C  Output Parameters:
C  PetscVecRHS - local vector data, PetscVecRHS(x)
C  PetscMatC - the jacobian matrix
C
C  Notes:
C
      IMPLICIT NONE
C
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
C
      Vec PetscVecRHS
      Mat PetscMatC
C
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
C
      double precision x(*)
      integer ctx(*),IFAIL
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     &NELEM,NFACE,NBFAC,IOPT
C
      CHARACTER*38 ERRMSG
      DATA ERRMSG(1:33)/'SUBR. TURBSCH, NO SUCH KSCHEME = '/
C
      NDIM   =ctx(1)
      NOFVERT=ctx(2)
      NOFVAR =ctx(3)
      NTURB  =ctx(4)
      NPOIN  =ctx(5)
      NGHOST =ctx(6)
      NELEM  =ctx(7)
      NFACE  =ctx(8)
      NBFAC  =ctx(9)
      NPNOD  =ctx(10)
C
C
C     Executable Statements ..
C
C
      CALL NODRES_B(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +              NELEM,NFACE,NBFAC,X,
     3              PetscVecRHS,PetscMatC)
C
      RETURN
      END
