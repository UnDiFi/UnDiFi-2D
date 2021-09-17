      SUBROUTINE SCLSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     +                  NELEM,NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel)
C
C     $Id: sclsch.F,v 1.23 2020/04/23 09:37:34 abonfi Exp $
C
C     Purpose:
C     This routine selects the SCALAR scheme to be used,
C     depending on the flag ISCHEME.
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
      Vec  DT,RHS
      Mat  A
C
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,NFACE,
     &NBFAC
      INTEGER NERR,IOPT
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
      EXTERNAL EulerModel,NSModel
C
      EXTERNAL N_scheme,LDA_scheme,PSI_scheme,NL_Scheme,FCT_scheme,
     + LW_scheme,LW2_scheme,NL2_scheme,NS_scheme,SUPG_scheme,FV_scheme,
     & LDA2_scheme,LDASqr_scheme,SDK_scheme
C
      CHARACTER*38 ERRMSG
C
      DATA ERRMSG(1:32)/'SUBR. SCLSCH, NO SUCH ISCHEME = '/
C
C     Executable Statements ..
C
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14),ISCHEME
      GOTO 9999

C
    1 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,N_scheme)
      RETURN

    2 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,PSI_scheme)
      RETURN
    3 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,LDA_scheme)
      RETURN

    4 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,NL_scheme)
      RETURN

    5 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,FV_scheme)
      RETURN

    6 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,FCT_scheme)
      RETURN

    7 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,NL2_scheme)
      RETURN
    8 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,LW_scheme)
      RETURN
    9 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,NS_scheme)
      RETURN
   10 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,SUPG_scheme)
      RETURN
   11 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,LDA2_scheme)
      RETURN
   12 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,LW2_scheme)
      RETURN
   13 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,
     +            LDASqr_scheme)
      RETURN 
   14 CALL MATSCH(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,X,DT,RHS,A,EulerModel,NSModel,
     +            SDK_scheme)
      RETURN
C
 9999 WRITE(ERRMSG(33:38),FMT="(I6.6)")ISCHEME
      NERR = 1
      IOPT = 1
      CALL SETERR(ERRMSG,38,NERR,IOPT)
C
      END
