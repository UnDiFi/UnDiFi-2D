      SUBROUTINE POISSON(IELEM,VCN,VCZ,NDIM,NOFVERT,NOFVAR,NDUMMY,
     +                  NODRES,TSTEP,STIFEL,VOLUME,PICARD)
C
      IMPLICIT NONE
C
C     $Id: laplace.f,v 1.3 2013/06/25 14:17:33 abonfi Exp $
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
C
C
C     FORTRAN stack
C
!     DOUBLE PRECISION DSTAK(1)
!     COMMON /CSTAK/ DSTAK
!     INTEGER ISTAK(1)
!     EQUIVALENCE(DSTAK(1),ISTAK(1))
C
C     ..
C     .. Common blocks ..
C     ..
C
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'nloc.com'
      INCLUDE 'three.com'
      INCLUDE 'time.com'
C
C
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NOFVERT,NDUMMY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TSTEP(NOFVAR,*),NODRES(NOFVAR,*),
     +STIFEL(NOFVERT,NOFVERT),VCN(*),VCZ(NOFVAR,*)
C     ..
C     .. Subroutine Arguments ..
!     EXTERNAL SCALARSCHEME,MATRIXSCHEME
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ERR,RESIDUAL,S,SOURCE,VOLUME,DIVB
      INTEGER IELEM,IFAIL,IOFF,I,J,IPOIN,IVERT
      LOGICAL PICARD
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION PHI(MAXNOFVERT),DUMMY(MAX_NOFVAR_SQR),
     2STIFC(MAX_NOFVERT_SQR),DT(MAXNOFVERT),RHS(MAXNOFVERT)
C     ..
C     .. External Functions ..
!     DOUBLE PRECISION DDOT,DNRM2,DIV
!     EXTERNAL DDOT,DNRM2,DIV
C     ..
C     .. External Subroutines ..
      EXTERNAL ADVECT,DINIT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN
C
C
C     Sets residual and local timestep to zero
C
      CALL DINIT(NOFVERT,ZERO,DT,1)
      CALL DINIT(NOFVERT,ZERO,RHS,1)
C
      DO 12 IVERT = 1, NOFVERT
         PHI(IVERT) = VCZ(NOFVAR,IVERT)
   12 CONTINUE
C
c     computes
c
      CALL VISCO(IELEM,PHI,RHS,DT,1,VCN,NDIM,NOFVERT,VOLUME,
     &STIFC,ONE,DUMMY,.FALSE.)
c
      DO 10 IVERT = 1, NOFVERT
         NODRES(NOFVAR,IVERT) = NODRES(NOFVAR,IVERT) + RHS(IVERT)
!        TSTEP(NOFVAR,IVERT) = TSTEP(NOFVAR,IVERT) + DT(IVERT)
         TSTEP(NOFVAR,IVERT) = ZERO
   10 CONTINUE
c
C     Assembling the Stiffness matrix ...
C
      IF (.NOT.PICARD) RETURN
      STOP 'Should NOT use Picard with Laplace'
C
      CALL DSCAL(NOFVERT*NOFVERT,MONE,STIFC,1) 
      CALL DCOPY(NOFVERT*NOFVERT,STIFC,1,STIFEL,1) 
C
      RETURN

! 200 FORMAT (5X,'Error on scalar residual in ELEM # ',I6,/,12X,'true',
!    +       17X,'computed',14X,'error',/,3 (10X,D12.5))
C
      END
