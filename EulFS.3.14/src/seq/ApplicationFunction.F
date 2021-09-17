      SUBROUTINE ApplicationFunction(PetscMatA,X,PetscVecRHS,
     &PetscVecDT,ctx,IERR)
C
C  $Id: ApplicationFunction.F,v 1.15 2020/04/23 09:37:34 abonfi Exp $
C
C ---------------------------------------------------------------------
C
C  ApplicationFunction - Computes nonlinear function
C  This function is called by
C  the higher level routines RHSFunction() and RHSJacobian()
C  Input Parameter:
C  x - local vector data
C  A - the jacobian matrix
C
C  Output Parameters:
C  PetscVecRHS - local vector data, PetscVecRHS(x)
C  A - the jacobian matrix
C  ierr - error code
C
C  Notes:
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
      IMPLICIT NONE
C
      Vec PetscVecRHS,PetscVecDT
      Mat PetscMatA
C
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
C
      double precision x(*)
      integer ierr
      integer ctx(*)
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     &NBFAC
C
      EXTERNAL SCALAR,EULERII,EULERIIBIS,EULERVII,EULERVIII,EULERIX,
     +EULERX,EULERXI,EULERXII,VSFLX2,VSFLX4,VISCO
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
C
C **************************************************************
C         Scalar case (ABS(KAN)=1;NOFVAR=1)
C **************************************************************
C
      IF (IABS(KAN).EQ.1) THEN
C
          CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                NELEM,NFACE,NBFAC,
     2                X,PetscVecDT,PetscVecRHS,PetscMatA,
     3                SCALAR,VISCO)
C
C **************************************************************
C         Euler and Navier-Stokes
C **************************************************************
C
      ELSEIF (IABS(KAN).EQ. 4.OR. IABS(KAN).EQ.2) THEN
C
C         Compute NODAL residuals ..
C
C         HE splitting for supersonic 2D flows
C
          IF (DECOMP.EQ.2) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERII, VSFLX4)
 
C
C         HE splitting
C
 
          ELSEIF (DECOMP.EQ.3) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERIIBIS, VSFLX4)
 
C
C         Hybrid
C
          ELSEIF (DECOMP.EQ.4) THEN
cold         CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
cold +                   NELEM,NFACE,NBFAC,
cold 2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
cold +                   EULERHYB, VSFLX4)
             STOP ' DECOMP == 4 is not implemented anymore'
C
C         Unsteady compressible in conserved variables
C
 
          ELSEIF (DECOMP.EQ.5) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERXI, VSFLX4)
C
C         Preconditioned Unsteady compressible using the Merkle preconditioner
C
 
          ELSEIF (DECOMP.EQ.6) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERX, VSFLX4)
C
C         Unsteady compressible
C
          ELSEIF (DECOMP.EQ.7) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERVII, VSFLX4)
 
C
C         Unsteady incompressible
C
          ELSEIF (DECOMP.EQ.8) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERVIII, VSFLX2)
 
C
C         Pseudo-Unsteady incompressible
C
          ELSEIF (DECOMP.EQ.9) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERIX, VSFLX2)
C
          ELSE
             STOP 'INVALID DECOMP in ApplicationFunction'
C
          ENDIF
C
      ELSEIF ( KAN .EQ. -3 )THEN
C
C         Plasma: Ar, Ar+, Ar++, e
C
          IF (DECOMP.EQ.10) THEN
             CALL SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                   NELEM,NFACE,NBFAC,
     2                   X,PetscVecDT,PetscVecRHS,PetscMatA,
     +                   EULERXII, VSFLX4)
          ELSE
             STOP 'INVALID DECOMP in ApplicationFunction'
C
          ENDIF
 
 
      ELSE
          STOP 'INVALID KAN in ApplicationFunction'
C
      ENDIF
C
C
      RETURN
      END
