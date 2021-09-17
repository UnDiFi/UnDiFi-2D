      SUBROUTINE MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                  NELEM,NFACE,NBFAC,
     2                  X,PetscVecDT,PetscVecRHS,PetscMatA,
     2                  EulerModel,NSModel,ScalarScheme)
C
C     $Id: matsch.F,v 1.23 2020/04/23 09:37:34 abonfi Exp $
C
C     Purpose:
C     This routine selects the MATRIX scheme to be used,
C     depending on the flag JSCHEME.
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
      EXTERNAL EulerModel,NSModel,ScalarScheme
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
C
      CHARACTER*38 ERRMSG
C
      EXTERNAL NSys_scheme,LDASys_scheme,LDA2Sys_scheme,LWSys_scheme,
     +LDAS_fast,CDSSys_scheme,LW2Sys_scheme,VKIB_scheme,
     +LSSys_scheme,FVSys_scheme,LWTSys_scheme,UCVSys_scheme
C
      DATA ERRMSG(1:32)/'SUBR. MATSCH, NO SUCH JSCHEME = '/
C
C     Executable Statements ..
C
C     Select the PETSc kernel used to invert small matrices
C
      CALL setsolven( IFAIL )
C
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12) JSCHEME

      GOTO 9999
C
    1 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,NSys_scheme)
      RETURN

    2 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LDAS_fast)

      RETURN

    3 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LDASys_scheme)
      RETURN

    4 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LW2Sys_scheme)

      RETURN

    5 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,FVSys_scheme)
      RETURN

    6 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LDA2Sys_scheme)
      RETURN

    7 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,CDSSys_scheme)

      RETURN

    8 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LWSys_scheme)
      RETURN
    9 STOP 'JSCHEME = 9 Not available any more'
      RETURN
   10 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,VKIB_scheme)
      RETURN
   11 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,LWTSys_scheme)
      RETURN
   12 CALL TURBSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,PetscVecDT,PetscVecRHS,PetscMatA,
     2            EulerModel,NSModel,ScalarScheme,UCVSys_scheme)
      RETURN
C
 9999 WRITE(ERRMSG(33:38),FMT="(I6.6)")JSCHEME
      IFAIL = 2
      IOPT = 1
      CALL SETERR(ERRMSG,38,IFAIL,IOPT)
C
      END
