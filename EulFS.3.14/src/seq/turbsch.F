      SUBROUTINE TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                  NELEM,NFACE,NBFAC,
     2                  X,PetscVecDT,PetscVecRHS,PetscMatA,
     2                  EulerModel,NSModel,ScalarScheme,MatrixScheme)
C
C     $Id: turbsch.F,v 1.3 2020/04/23 09:37:34 abonfi Exp $
C
C     Purpose:
C     This routine selects the scalar scheme for turbulence modeling to be used,
C     depending on the flag KSCHEME.
C     The idea is to pass the schemes as subroutines
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
      INCLUDE 'flags.com'
C
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     &NBFAC
      INTEGER IFAIL,IOPT
C
      Vec  PetscVecDT,PetscVecRHS
      Mat  PetscMatA
      DOUBLE PRECISION X(*)
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
      EXTERNAL EulerModel,NSModel,ScalarScheme,MatrixScheme
C
C     On entry:
C     --------
C     EulerModel    
C              subroutine of the model equation being solved, e.g.
C              EulerII,EulerVII,SCALAR,etc.
C     NSModel    
C              subroutine of the model equation being solved, e.g.
C              VSFLX2,VSFLX4
C     ScalarScheme 
C              subroutine of the scalar scheme being used, e.g.
C              N_scheme, NL_scheme, etc. 
C     MatrixScheme 
C              subroutine of the matrix scheme being used, e.g.
C              N_scheme, NL_scheme, etc. 
C
      CHARACTER*38 ERRMSG
C
      EXTERNAL N_scheme,PSI_scheme,LDA_scheme,NL_scheme,FV_scheme,
     +FCT_scheme,NL2_scheme,NS_scheme,LW2_scheme,LDA2_scheme,LW_scheme,
     +LDASqr_scheme,SUPG_scheme,SDK_scheme
C
      DATA ERRMSG(1:33)/'SUBR. TURBSCH, NO SUCH KSCHEME = '/
C
C     Executable Statements ..
C
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14) KSCHEME

      GOTO 9999
C
    1 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            N_scheme)
      RETURN
    2 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            PSI_scheme)
      RETURN
    3 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            LDA_scheme)
      RETURN
    4 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            NL_scheme)
      RETURN
    5 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            FV_scheme)
      RETURN
    6 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            FCT_scheme)
      RETURN
    7 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            NL2_scheme)
      RETURN
    8 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            LW_scheme)
      RETURN
    9 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            NS_scheme)
      RETURN
   10 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            SUPG_scheme)
      RETURN
   11 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            LDA2_scheme)
      RETURN
   12 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            LW2_scheme)
      RETURN
   13 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            LDASqr_scheme)
      RETURN
   14 CALL TMODEL(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,MatrixScheme,
     3            SDK_scheme)
      RETURN
C
 9999 WRITE(ERRMSG(33:38),FMT="(I6.6)")KSCHEME
      IFAIL = 2
      IOPT = 1
      CALL SETERR(ERRMSG,38,IFAIL,IOPT)
C
      END
