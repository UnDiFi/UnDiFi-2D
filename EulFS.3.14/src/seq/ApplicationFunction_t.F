      SUBROUTINE ApplicationFunction_t(PetscMatB,X,PetscVecRHS,
     &                                 PetscVecDT,ctx,IFAIL)
C ---------------------------------------------------------------------
C
C  ApplicationFunction_t - Computes nonlinear function
C  for the turbulence model (segregated approach)
C  This function is called by
C  the higher level routines RHSFunction_t() and RHSJacobian_t()
C
C  $Id: ApplicationFunction_t.F,v 1.11 2020/04/23 09:37:34 abonfi Exp $
C
C  Input Parameter:
C  x - local vector data as a FORTRAN array
C  A - the jacobian matrix
C
C  Output Parameters:
C  PetscVecRHS - local vector data, PetscVecRHS(x)
C  PetscMatB - the jacobian matrix
C
C  Notes:
C
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec 
      use petscmat 
C
      IMPLICIT NONE
      Vec PetscVecRHS,PetscVecDT
      Mat PetscMatB
C
      INCLUDE 'flags.com'
      INCLUDE 'implicit.h'
C
      double precision x(*)
      integer ctx(*),IFAIL
      INTEGER NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
     &NELEM,NFACE,NBFAC,IOPT
C
      EXTERNAL TURBO,NSYS_SCHEME
C
      EXTERNAL N_scheme,PSI_scheme,LDA_scheme,NL_scheme,FV_scheme,
     +FCT_scheme,NL2_scheme,NS_scheme,LW2_scheme,LDA2_scheme,LW_scheme,
     +LDASqr_scheme,SUPG_scheme
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
!     CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,
!    +              NELEM,NFACE,NBFAC,PetscVecDT,X,
!    3              PetscVecRHS,PetscMatB,TURBO,
!    +              NL2_SCHEME,NSYS_SCHEME)
!    
C
C
C     Executable Statements ..
 
C
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13) KSCHEME

      GOTO 9999
C
    1 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            N_scheme,NSYS_SCHEME)
      RETURN
    2 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            PSI_scheme,NSYS_SCHEME)
      RETURN
    3 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            LDA_scheme,NSYS_SCHEME)
      RETURN
    4 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            NL_scheme,NSYS_SCHEME)
      RETURN
    5 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            FV_scheme,NSYS_SCHEME)
      RETURN
    6 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            FCT_scheme,NSYS_SCHEME)
      RETURN
    7 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            NL2_scheme,NSYS_SCHEME)
      RETURN
    8 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            LW_scheme,NSYS_SCHEME)
      RETURN
    9 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            NS_scheme,NSYS_SCHEME)
      RETURN
   10 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            SUPG_scheme,NSYS_SCHEME)
      RETURN
   11 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            LDA2_scheme,NSYS_SCHEME)
      RETURN
   12 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            LW2_scheme,NSYS_SCHEME)
      RETURN
   13 CALL NODRES_T(NDIM,NOFVERT,NOFVAR,NTURB,NPOIN,NGHOST,NPNOD,NELEM,
     +            NFACE,NBFAC,PetscVecDT,X,PetscVecRHS,PetscMatB,
     2            TURBO,
     3            LDASqr_scheme,NSYS_SCHEME)
      RETURN
C
 9999 WRITE(ERRMSG(33:38),FMT="(I6.6)")KSCHEME
      IFAIL = 2
      IOPT = 1
      CALL SETERR(ERRMSG,38,IFAIL,IOPT)
C
C
      RETURN
      END
