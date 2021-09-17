      SUBROUTINE TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                  NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,
     3                  PetscMatA,EulerModel,NSModel,ScalarScheme,
     4                  MatrixScheme,TModelScheme)
C
C     $Id: tmodel.F,v 1.14 2020/04/23 09:37:34 abonfi Exp $
C
C     Purpose:
C     This routine selects the turbulence model to be used.
C     Only Spalart-Allmaras is currently implemented
C     For laminar flow a dummy routine is used, which returns
C     the (cell averaged) laminar viscosity
C     The idea is to pass the turbulence model as a subroutine
C     at the uppermost level (i.e. at the non-linear iteration
C     loop) to avoid multiple IFs inside the innermost loop
C     over the elements.
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
C
      IMPLICIT NONE
C
C
      Vec  PetscVecDT,PetscVecRHS
      Mat  PetscMatA
C
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     &NBFAC
      INTEGER IFAIL,IOPT
      DOUBLE PRECISION X(*)
C
      INCLUDE 'flags.com'
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NOFVAR  number of variables (degrees of freedom)
C             in each meshpoint
C     NTURB   number of variables (degrees of freedom)
C             in each meshpoint for the turbulence model
C     NPOIN   no. of interior (processor owned) meshpoints;
C             global number of meshpoints in the uni-processor case
C     NGHOST  no. of ghost meshpoints on the current processor;
C             0 the uni-processor case
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C     NFACE   number of edges/faces in the mesh for the current
C             processor (multi-processor case) or global number
C             of edges/faces in the mesh (uni-processor case).
C     NBFAC   no. of processor owned boundary faces/edges;
C             global number of boundary faces/edges
C             in the uni-processor case.
C
      EXTERNAL EulerModel,NSModel,ScalarScheme,MatrixScheme,
     2         TModelScheme,SA8,SA7,SA7b,SA
C
      EXTERNAL Laminar,EVALTSA,Tdummy,EulerVIIsa
C
      CHARACTER*38 ERRMSG
      LOGICAL COMPRESSIBLE
C
C     Executable Statements ..
C
      COMPRESSIBLE = (ABS(KAN).EQ.4)
C
C     Scalar problems
C
      IF(IABS(KAN).EQ.1)THEN
          CALL NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     2         NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     4         Scalarscheme,MatrixScheme,TModelScheme,
     5         EulerModel,NSModel,Tdummy)
      ELSE
C
C     Euler Navier-Stokes
C
          IF(TURBULENT.AND.SEGREGATED)THEN
C
              CALL NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     2             NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     4             Scalarscheme,MatrixScheme,TModelScheme,
     5             EulerModel,NSModel,EVALTSA)
C
          ELSEIF(TURBULENT.AND.COUPLED)THEN
C
            IF(COMPRESSIBLE)THEN
              CALL NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     2             NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     3             Scalarscheme,MatrixScheme,TModelScheme,EulerModel,
!    4             NSModel,SA)
     4             NSModel,SA7)
            ELSE
              CALL NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     2             NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     3             Scalarscheme,MatrixScheme,TModelScheme,EulerModel,
     4             NSModel,SA8)
            ENDIF
C             ********laminar********
          ELSE
              CALL NODRES(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     2             NELEM,NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     3             Scalarscheme,MatrixScheme,TModelScheme,
     4             EulerModel,NSModel,LAMINAR)
C
          ENDIF
      ENDIF
      RETURN
      END
