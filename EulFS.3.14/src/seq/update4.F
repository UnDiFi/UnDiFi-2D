      SUBROUTINE UPDATE4( NDIM, NTURB, NPOIN, solver, B, RHS, ViscT, DT,
     +                    NIT, IFLAG, NodalBcs )
C
C     $Id: update4.F,v 1.39 2020/04/23 09:46:19 abonfi Exp $
C
CCC#define DO_NOT_UPDATE
C
C     Implicit timestepping for the one-eqn. turbulence model ..
C
CCC#include "petsc/finclude/petscvec.h"
CCC#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
CCC#include "petsc/finclude/petscis.h"
!     use petscvec
!     use petscmat
!     use petscis
      use petscksp
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
      INTEGER NDNM
      PARAMETER(NDNM=3*NMAX)
      INCLUDE 'conv.com'
      INCLUDE 'nboun.com'
      INCLUDE 'verbose.com'
      INCLUDE 'io.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C
C
      Mat B
      Vec RHS,DT,ViscT
      Vec X
      KSP solver
      IS, dimension(0:*) :: NodalBcs
C
      INTEGER NDIM,NTURB,NPOIN,IFAIL,ITS,NIT
C
      INTEGER IVAR,I,IFLAG,N
      DOUBLE PRECISION CNST,RES0,RATIO2,RATIO1,EPS
      PARAMETER(EPS=1.d-20)
      SAVE RES0,RATIO2,RATIO1,CNST
C
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP1(3,MAXNOFVAR),WKSP2(3,MAXNOFVAR),
     &                 WKSP3(3,MAXNOFVAR),WKSP4(3,MAXNOFVAR),
     &                 WKSP5(3,MAXNOFVAR), 
     &                 AL2(MAXNOFVAR,5),AMAX(MAXNOFVAR,5)
      INTEGER WHEREMAX(MAXNOFVAR,5)
C
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
C
      DATA WKSP1,WKSP2,IFAIL,ITS/NDNM*ZERO,NDNM*ZERO,0,-1/
C
      IFLAG = 100
C
C     store the residual of the turbulence eqn. at the previous
C     time level
C
      IF( NIT .EQ. 1 ) RATIO1 = RESL2(1,2)
C
!     CALL VecNorm(RHS,NORM_2,RESL2(1,2),IFAIL)
!     CALL VecNorm(RHS,NORM_INFINITY,RESMAX(1,2),IFAIL)
C
      CALL FindVecStrideAbsMinMax(RHS,RESL2(1,2),RESMAX(1,2),INMAX(1,2),
     &WKSP1,NDIM)
C
      CALL FindVecStrideMinMax(DT,AL2(1,4),AMAX(1,4),WHEREMAX(1,4),
     &WKSP4,NDIM)
C
C
      IF( NIT .EQ. 1 ) THEN
C
C     RATIO1 = RES(ITER-1,NIT)/RES(ITER,1) where ITER refers to the outer
C             iteration counter 1 <= ITER <= NITER while NIT refers
C             to the inner iteration counter 1 <= NIT <= NSUBIT
C
C     store the residual at the first sub-iteration
C
          RES0 = RESL2(1,2)
          RATIO1 = RATIO1/RES0
      ENDIF
      IF (ITER.EQ.1) RATIO1 = 1.d0
C
C     RATIO2 = RES(ITER,1)/RES(ITER,NIT)
C
      RATIO2 = RES0/RESL2(1,2) 
C
      CNST = MAX(1.d0,RATIO1)*MAX(1.d0,RATIO2)
      CNST = MIN(CFL(2)*CNST,CFLMAX(2))
C
C     divide V_i/Dt by the CFL number ...
C
      CALL VecScale(DT,ONE/CNST,IFAIL)
C
C     Adds V_i/Dt to the diagonal elements of B ...
C
      CALL MatDiagonalSet(B,DT,ADD_VALUES,IFAIL)
C
C     We do not need DT any more
C
      CALL VecDestroy(DT,IFAIL) 
C
C     modify the stiffness matrix structure to account for
C     no-slip wall b.c. this is done only for the first iteration
C
caldo IF( ITER .EQ. 1 )THEN
          CALL LHSBC5( B, NodalBcs )
caldo ENDIF
C
      N=NPOIN*NTURB
C
C     create a vector to store the solution
C
#ifdef MPI
      CALL VecCreateMPI(PETSC_COMM_WORLD,N,PETSC_DECIDE,X,IFAIL)
#else
      CALL VecCreateSeq(PETSC_COMM_SELF,N,X,IFAIL)
#endif
#ifdef DO_NOT_UPDATE
      WRITE(6,*)'Not updating turbulence'
      goto 656
#endif
C
C     Sets the matrix associated with the linear system
C
      CALL KSPSetOperators(solver,B,B,IFAIL)
      CALL KSPSetup(solver,IFAIL)
C
C     Solve THE linear system
C
      CALL KSPSolve(solver,RHS,X,IFAIL)
      CALL KSPGetIterationNumber(solver,ITS,IFAIL)
C
C     Update the nodal unknown vector by forming ViscT := X + ViscT
C
      CALL VecAXPY(ViscT,ONE,X,IFAIL)
  656 continue
C
C     Monitor the norms of the update .. 

      CALL VecNorm(X,NORM_2,DELL2(1,2),IFAIL)
      CALL VecNorm(X,NORM_INFINITY,DELMAX(1,2),IFAIL)
C
C     Monitor the L2 and L infinity norms of the update ..
C
      CALL FindVecStrideAbsMinMax(X,DELL2(1,2),DELMAX(1,2),INDEL(1,2),
     &WKSP2,NDIM)
      CALL FindVecStrideMinMax(ViscT,AL2(1,5),AMAX(1,5),WHEREMAX(1,5),
     &WKSP5,NDIM)
C
C
      IF(MY_PE.EQ.0)THEN
          IF ((ITER/ISTMP)*ISTMP.EQ.ITER) THEN
             WRITE(IWUNIT,200)ITER,ITS,CNST
C
C     Print out the convergence history ..
C
             WRITE(IWUNIT,215)
             DO 20 IVAR = 1 , NTURB
                WRITE(IWUNIT,210)IVAR,DLOG10(MAX(EPS,RESL2(IVAR,2))),
     &          DLOG10(MAX(EPS,RESMAX(IVAR,2))),INMAX(IVAR,2),
     2          (WKSP1(I,IVAR),I=1,3)
   20        CONTINUE
             WRITE(IWUNIT,225)
             DO 30 IVAR = 1 , NTURB
                WRITE(IWUNIT,210)IVAR,DLOG10(MAX(EPS,DELL2(IVAR,2))),
     &          DLOG10(MAX(EPS,DELMAX(IVAR,2))),INDEL(IVAR,2),
     2          (WKSP2(I,IVAR),I=1,3)
   30        CONTINUE
              WRITE (IWUNIT,FMT=240)
              DO 34 IVAR = 1,NTURB
                  WRITE (IWUNIT,FMT=212) IVAR,
     &            AL2(IVAR,4), AMAX(IVAR,4),
     &            WHEREMAX(IVAR,4), (WKSP4(I,IVAR),I=1,3)
   34         CONTINUE
              WRITE (IWUNIT,FMT=230)
              DO 36 IVAR = 1,NTURB
                  WRITE (IWUNIT,FMT=212) IVAR,
     &            AL2(IVAR,5), AMAX(IVAR,5),
     &            WHEREMAX(IVAR,5), (WKSP5(I,IVAR),I=1,3)
   36         CONTINUE
         ENDIF
C
C     Writing convergence history to file ...
C
      WRITE(IHST3,235)NITER,ITS,RESL2(1,2),RESMAX(1,2),CNST
      ENDIF
       
C     clear memory allocated for the solution vector
C
      CALL VecDestroy(X,IFAIL)
C
C     this is the Zingg stuff
C
!     IF(RESL2(1,2).LE.RESL2(IVCNVG,1))IFLAG = 0
C
      RETURN
  200 FORMAT(5X,70('-'),/,14X,'TURBULENCE ITERATION # ',I4,
     +' (',I4,') CFL = ',E10.4/,5X,70('-'),/,
     +5X,70('-'),/,5X,'Var.',4X,'L2-norm',3X,'L_infty',3X,
     +'node #',3X,'(',4X,'x',7X,'y',7X,'z',4X,')',/,5X,70('-'))
  210 FORMAT(5X,I1,5X,F10.5,1X,F10.5,2X,I5,3X,
     +'(',2(F8.5,','),F8.5,')')
  212 FORMAT (5X,I1,5X,E10.3,1X,E10.3,1X,I6,3X,'(', 2(F8.4,','),F8.4,
     +')')
  215 FORMAT(5X,'Nodal Residual'/)
  225 FORMAT(5X,'Nodal Update'/)
  230 FORMAT (5X,'Nodal Values',/)
  235 FORMAT (I5,1X,I4,3 (1X,E10.4))
  240 FORMAT (5X,'Timestep ',/)
C
      END
