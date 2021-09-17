      SUBROUTINE TURBINI(NodalBcs,NOFVERT,NOFVAR,NDIM,NP,
     &                   TURB,ZROE,TD,TTD,RESTART,ROWBGN,
     &                   VISCTFILE,WDISTFILE)
C  ==============================================================
C
C
C
C      $Id: turbini.F,v 1.27 2020/04/23 09:46:19 abonfi Exp $
C
C      DESCRIPTION:
C      -----------
C
C      initialize turbulence computation
C      (Spalart & Allmaras model)
C
C      REM: NP is the number of meshpoints in the uni-processor case
C           in the multi-procs case it is the NPOIN+NGHOST+NPNOD
C           
C
C  ==============================================================
C
C
#include "petsc/finclude/petscis.h"
      use petscis
C
      IMPLICIT NONE
c2345678
C
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
      INCLUDE 'iset.h'
      INCLUDE 'turb.com'
      INCLUDE 'nloc.com'
      INCLUDE 'nboun.com' 
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
C
C
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVERT,NOFVAR,NP,ROWBGN
      LOGICAL RESTART
      CHARACTER*(*) VISCTFILE,WDISTFILE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TURB(NP),ZROE(NP),TD(NP),TTD(NP)
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      INTEGER IROW,NTURB,NI,IFAIL,I,IOPT,NERR
      CHARACTER TRIPFILE*11,ERRMSG*72
      DOUBLE PRECISION TVAL
C
C     ..
C     .. Petsc stuff
C
      PetscBool IFLAG
      PetscOffset IDX_I
      IS, dimension(0:*) :: NodalBcs 
C
      INTEGER IDX_V(1)
C     ..
C     .. Local Arrays ..
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,SOLZNE
C     ..
C     .. Common blocks ..
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C     ..
C     .. Equivalences ..
C     ..
C     .. Data statements ..

      DATA NTURB/1/
      DATA TRIPFILE/'file014.dat'/
C
      IF(MY_PE.EQ.0)WRITE(IWUNIT,1000)
C
C************************************************************
C
C Set trip point variables
C
C************************************************************
C23456
C
C**************************************************************
C
C     when restarting, read turbulent viscosity
C
      IF (RESTART.AND.SEGREGATED) 
     &CALL SOLZNE(VISCTFILE,TURB,NTURB,NP,'r')
C
C**************************************************************
C
C Turbulent viscosity initialization
C
C*************************************************************
C
      CALL PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,
     +        '-inflow_turbulent_viscosity',TVAL,IFLAG,IFAIL)
      IF(IFLAG.EQV.PETSC_FALSE)THEN
          WRITE(ERRMSG,300)'-inflow_turbulent_viscosity [value]'
          NERR = 14
          IOPT = 1
          CALL SETERR(ERRMSG,69,NERR,IOPT)
      ELSE
          IF(SEGREGATED.AND.(.NOT.RESTART))THEN
              CALL DINIT(NP,TVAL,TURB,1)
          ELSEIF(COUPLED.AND.(.NOT.RESTART))THEN 
              CALL DINIT(NP,TVAL,ZROE(NOFVAR),NOFVAR)
          ENDIF
      ENDIF
C
C Retrieve from the Index Set the gridpoints where no-slip
C boundary conditions are applied so that turbulent viscosity
C can be set to 0. there
C
      CALL ISGetSize(NodalBcs(NoSlipNodes),NI,IFAIL)
      IF(NI.NE.0)THEN
         CALL ISGetIndices(NodalBcs(NoSlipNodes),IDX_V,IDX_I,IFAIL)
C
C     set turbulent viscosity = 0.0 on no-slip walls
C
         DO 110 I = 1,NI
             IROW = IDX_V(IDX_I+I)+1-ROWBGN
             IF(SEGREGATED)THEN
                 TURB(IROW) = ZERO
             ELSE 
                 ZROE(IROW) = ZERO
             ENDIF
  110    CONTINUE
         CALL ISRestoreIndices(NodalBcs(NoSlipNodes),IDX_V,IDX_I,IFAIL)
      ENDIF
C
C**************************************************************
C
C Read Wall distance from a file
C
C*************************************************************

          I = 1
          CALL SOLZNE(WDISTFILE,TD,I,NP,'r')
C
C read trip point distance from file011.dat
C
          IF(TTFLAG.EQ.1)CALL SOLZNE('file011.dat',TTD,I,NP,'r')
C
      RETURN
  300 FORMAT(5X,'Need to specify option: ',A35)
 1000 FORMAT(//5X,'TURBULENCE MODELING: Spalart-Allmaras',/,5X,38("="))
      END
