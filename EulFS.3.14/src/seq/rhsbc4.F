      SUBROUTINE RHSBC4(X,RHS,NodalBcs,ibgn,NDIM,PRATIO,COMPRESSIBLE)
C
C     This routine sets the boundary conditions in the
C     nodal residual ..
C
C     $Id: rhsbc4.F,v 1.33 2020/04/23 09:46:19 abonfi Exp $
C
C
#include "petsc/finclude/petscvec.h"
      use petscvec
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bctype.h'
      INCLUDE 'iset.h'
C
      DOUBLE PRECISION X(*),RHS(*),PRATIO(3,*)
      LOGICAL COMPRESSIBLE
C
C
      INCLUDE 'stream.com'
      INCLUDE 'flags.com'
      INCLUDE 'visco.com'
      INCLUDE 'pfcgas.com'
C
C
      Integer IDX_V(1),IDX_X(1)
      PetscScalar X_V(1)
C
      PetscOffset RHS_I,IDX_I
      IS, dimension(0:*) :: NodalBcs
C
      INTEGER IPOIN,I,NI,IFAIL,I1,I2,I3,I4,I5,IBGN,NDIM,LDA,J
      DOUBLE PRECISION ALPHA,DUM,U1,U2
!     DOUBLE PRECISION KINETIC,ASQR,MACHSQR,P0,TEMP,DENS
      DOUBLE PRECISION PRESSC
      EXTERNAL PRESSC 
C
C    EXPLICIT Treatment of supersonic and/or 
C             viscous nodes
C
C    note that the indices for are
C         given in global numbering (and 0-based indexing)
C         the reason for using global numbering is that
C         these index sets are among the arguments
C         of the call to MatZeroRows
C
C ********************************************************************
C Supersonic variables: a Dirichlet condition (IBC=1)
C this is in the form of NO increment
C d(r,rE,ru,rv,rw) = 0
C ********************************************************************
      CALL ISGetSize(NodalBcs(SupersonicVariables),NI,IFAIL)
      IF(NI.NE.0)THEN
         CALL ISGetIndices(NodalBcs(SupersonicVariables),IDX_V,IDX_I,
     &                     IFAIL)
         DO 100 I = 1,NI
             IPOIN = IDX_V(IDX_I+I)+1-IBGN
             RHS(IPOIN) = ZERO
  100    CONTINUE
         CALL ISRestoreIndices(NodalBcs(SupersonicVariables),IDX_V,
     &                         IDX_I,IFAIL)
      ENDIF
C ********************************************************************
C No-slip & no-penetration boundary condition (IBC=6)
C this is in the form of NO increment
C d(ru,rv,rw) = 0
C so a NON zero velocity would be preserved
C ********************************************************************
      CALL ISGetSize(NodalBcs(NoSlipVelocities),NI,IFAIL)
      IF(NI.NE.0)THEN
         CALL ISGetIndices(NodalBcs(NoSlipVelocities),IDX_V,IDX_I,IFAIL)
         DO 110 I = 1,NI
             IPOIN = IDX_V(IDX_I+I)+1-IBGN
             RHS(IPOIN) = ZERO
  110    CONTINUE
         CALL ISRestoreIndices(NodalBcs(NoSlipVelocities),IDX_V,IDX_I,
     &                         IFAIL)
      ENDIF
C
C ********************************************************************
C Fixed temperature wall boundary condition (IBC=6,IADIA=1)
C this is in the form
C d(T) = T(required) - T^n
C so can also be applied with the ``wrong'' initial condition
C
C ************* W A R N I N G ******************
C this is currently implemented only for 
C -nondimensionalisation external
C
C
C ********************************************************************
C
      CALL ISGetSize(NodalBcs(Isothermal),NI,IFAIL)
      IF(NI.NE.0)THEN
         CALL ISGetIndices(NodalBcs(Isothermal),IDX_V,IDX_I,IFAIL)
c
c
C        ALPHA = TWALL/GAM/GM1/(M_infty*M_infty)
         ALPHA = RSTAR*TWALL/GM1
C
         DO 120 I = 1,NI
C
C N.B. X is the parameter vector 
C
C
C I2  1-based address of the location where \sqrt(\rho) H is stored
C     (+1) because of the 0-based indexing in the index set
             I2 = IDX_V(IDX_I+I)+1-IBGN
             I1 = I2-1
             I3 = I2+1
             I4 = I3+1
             I5 = I4+1
             DUM = X(I3)*X(I3)+X(I4)*X(I4)
             IF(NDIM.EQ.3)DUM = DUM + X(I5)*X(I5)
             U2 = GINV * X(I1)*X(I2)+HALF*GM1OG*DUM
             U1 = X(I1)*X(I1)
             RHS(I2) = ALPHA*U1-U2
!            write(6,*)i2,DUM,alpha*U1,U2
  120    CONTINUE
         CALL ISRestoreIndices(NodalBcs(Isothermal),IDX_V,IDX_I,IFAIL)
      ENDIF
C
C ********************************************************************
C Dirichlet condition (IBC=1)
C this is in the form of NO increment
C d(phi) = 0
C ********************************************************************
      IF( LAPLACE )THEN
         CALL ISGetSize(NodalBcs(Dirichlet4Poisson),NI,IFAIL)
         IF(NI.NE.0)THEN
            CALL ISGetIndices(NodalBcs(Dirichlet4Poisson),IDX_V,IDX_I,
     &                        IFAIL)
            DO 104 I = 1,NI
               IPOIN = IDX_V(IDX_I+I)+1-IBGN
               RHS(IPOIN) = ZERO
  104       CONTINUE
            CALL ISRestoreIndices(NodalBcs(Dirichlet4Poisson),IDX_V,
     &                            IDX_I,IFAIL)
         ENDIF
      ENDIF ! test on Poisson
C
      RETURN
      END
