      SUBROUTINE UPDATE3(NDIM,NOFVAR,NPOIN,FlowSolver,A,RHS,DT,ZROE,
     &                   VMEDIAN)
C
C     $Id: update3.F,v 1.57 2009/06/11 13:10:58 abonfi Exp abonfi $
C
      IMPLICIT NONE
C
C#define INITIAL_GUESS_NONZERO
C#define PRINT_KSP
CCCC#define DEBUG
CC#define DO_NOT_UPDATE
C
C     This routine updates the nodal solution by implicit
C     time integration when KAN = 2,4
C
C     N.B.: upon entry ZROE contains the vector of CONSERVED
C           variables
C
C
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscis.h"
#include "finclude/petscviewer.h"
#include "finclude/petscpc.h"
C
      Mat A
      Vec RHS,DT,ZROE,X,LEFT
      SAVE X
      KSP FlowSolver
      PetscLogDouble TBEGALL,telapsed,tbeg,tend
      DOUBLE PRECISION VMEDIAN(*)
C
      Mat AijMat
      MatType mt
C
caldo
      PC mypc
      PC my_sub_pc
      KSP my_sub_ksp
      Mat factored_mat
      PetscTruth matflag
      DOUBLE PRECISION INFO(MAT_INFO_SIZE)
caldo
#include "iset.com"
C
      PetscTruth flag,flg
C
      COMMON /TIMING/TBEGALL
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.h'
      INCLUDE 'datatype.h'
      INCLUDE 'time.com'
      INCLUDE 'conv.com'
      INCLUDE 'visco.com'
      INCLUDE 'nloc.com'
      INCLUDE 'verbose.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
C
C     .. Parameters ..
      INTEGER NDNM
      PARAMETER (NDNM=3*MAXNOFVAR)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NPOIN
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CNST,DTMAX,eps,XMIN,XMAX,S,THRESH
      parameter(eps=1.d-20) 
      INTEGER I,IVAR,N,ITS,IFAIL,MY_PE,ROWBGN,ROWEND,IPOIN,NI,INDX,
     &IADDR,INDXMIN,INDXMAX,LIWORK,LJWORK,LDWORK,NITEMS,NR,NZR
      integer nofeqn
      CHARACTER*255 string
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP1(3,MAXNOFVAR),WKSP2(3,MAXNOFVAR),
     &  WKSP3(3,MAXNOFVAR),WKSP4(3,MAXNOFVAR),
     &  WKSP5(3,MAXNOFVAR), AL2(MAXNOFVAR,5),AMAX(MAXNOFVAR,5)
      INTEGER WHEREMAX(MAXNOFVAR,5)
      INTEGER ISTAK(1)
      PetscScalar x_array(1)
      INTEGER idx_v(1),i_x
      PetscOffset idx_i
#ifdef DEBUG
      PetscScalar b_array(1)
      PetscOffset i_b
      integer j,k
#endif
C
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISTKGT,ISTKST
      EXTERNAL DNRM2,ISTKGT,ISTKST
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DLOG10,MIN
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON /MPICOM/MY_PE
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
C     .. Data statements ..
      DATA WKSP1,WKSP2,WKSP3,WKSP4,WKSP5,ITS/NDNM*ZERO,NDNM*ZERO,
     &NDNM*ZERO,NDNM*ZERO,NDNM*ZERO,0/
C
!     SAVE OMEGA
C     ..
#ifdef INITIAL_GUESS_NONZERO
      PARAMETER(flag = PETSC_TRUE)
#else
      PARAMETER(flag = PETSC_FALSE)
#endif
C 
C     Executable statements 
C
#ifdef DEBUG
      do ivar = 1,nofvar
      CALL VecStrideNorm(ZROE,ivar-1,NORM_2,s,IFAIL)
      IF(MY_PE.EQ.0)write(6,*)'U(',ivar,')',s
      enddo
#endif
      CALL VecGetOwnerShipRange(RHS,ROWBGN,ROWEND,IFAIL)
C
#ifdef DEBUG
      call VecNorm(RHS,NORM_2,s,ifail)
         write(6,*)'||RHS|| before  is ',s
         do ivar = 1,nofvar
         CALL VecStrideNorm(RHS,ivar-1,NORM_2,s,IFAIL)
      if(MY_PE.EQ.0)then
         write(6,*)'rhs(',ivar,')',s
      endif
         enddo
#endif
      CALL FindVecStrideAbsMinMax(RHS,RESL2(1,1),RESMAX(1,1),INMAX(1,1),
     &WKSP1,NDIM)
C
      IF (ITER.EQ.1) THEN
          RESMAX0(1) = RESMAX(IVCNVG,1)
          RESL20(1) = RESL2(IVCNVG,1)
      ENDIF
C
      IF( TIMEIMPL )THEN
          IF( CFL_RAMP .EQ. CFL_SER )THEN
!            write(6,*)'Omega is = ',omega
             CNST = RESL20(1)/RESL2(IVCNVG,1) * OMEGA
             CNST = MIN(CFLMAX(1),CFL(1)*CNST)
          ELSEIF( CFL_RAMP .EQ. CFL_EXP )THEN
             CNST = REAL(ITER-1)/REAL(ITMAX-1)
             CNST = CFL(1)+((EXP(PE*CNST)-ONE)/(EXP(PE)-ONE))*
     &              (CFLMAX(1)-CFL(1))
          ENDIF
      ELSE
!         CNST=ONE/CFL(1)
          CNST=CFL(1)
      ENDIF
!     CALL VecGetArray(rhs,x_array,i_x,IFAIL)
!     CALL solzne("res010.dat",x_array(i_x+1),nofvar,npoin,"w")
!     CALL VecRestoreArray(rhs,x_array,i_x,IFAIL)
C
C     For pure advection problems, the inverse of
C         the local timestep in the Dirichlet Nodes
C         is 0., so we set it to the maximum DT
C
C to be changed
C
      IF(KAN.EQ.-4.OR.KAN.EQ.-2)THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,x_array,i_x,IFAIL)
          CALL ISGetIndices(SupersonicVariables,IDX_V,IDX_I,IFAIL)
          CALL ISGetSize(SupersonicVariables,NI,IFAIL)
          DO 10 I = 1, NI, NOFVAR
C     supersonic nodes are 0-based indexed
              IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
C             write(6,*)x_array(i_x+IPOIN) , DTMAX
              x_array(i_x+IPOIN) = DTMAX
   10     CONTINUE
          CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          CALL ISRestoreIndices(SupersonicVariables,IDX_V,IDX_I,IFAIL)
      ENDIF
C
#ifdef DEBUG
      IF(TIMEIMPL)THEN
         call MatNorm(A,NORM_FROBENIUS,s,ifail)
         if(MY_PE.EQ.0)write(6,*)'||A|| before is ',s
      ENDIF
#endif
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||dt|| before myroutine is ',s
#endif
C
C     Handles hanging nodes
C
      CALL ISGetSize(HangingNodes,NI,IFAIL)
      IF( NI .GT. 0 .AND. IGLOB .EQ. 1 )THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,x_array,i_x,IFAIL)
          CALL ISGetIndices(HangingNodes,IDX_V,IDX_I,IFAIL)
          DO 12 I = 1, NI, NOFVAR
C     hanging nodes are 0-based indexed
              IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
              x_array(i_x+IPOIN) = DTMAX
   12     CONTINUE
          CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          CALL ISRestoreIndices(HangingNodes,IDX_V,IDX_I,IFAIL)
      ENDIF
C
      IF( ABS(KAN) .EQ. 4 )THEN
          NOFEQN = NDIM + 2
      ELSEIF( ABS(KAN) .EQ. 2 )THEN
          NOFEQN = NDIM + 1
      ELSE
          NOFEQN = NOFVAR
      ENDIF 
C
      CALL VecGetArray(DT,x_array,i_x,IFAIL)
      call mysillyroutine(x_array(i_x+1),NOFVAR,NOFEQN,NPOIN)
C     CALL solzne("vdt010.dat",x_array(i_x+1),nofvar,npoin,"w")
      CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||dt|| after myroutine is ',s
#endif
C
      CALL FindVecStrideMinMax(DT,AL2(1,4),AMAX(1,4),WHEREMAX(1,4),
     &WKSP4,NDIM)
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
       if(MY_PE.EQ.0)write(6,*)'||dt|| after FindVecStrideMinMax is ',s
#endif
C
C     Global time-stepping (if required)
C
      IF( IGLOB .EQ. 0 )THEN
          write(6,*)'Hello! IGLOB = ',IGLOB,PETSC_TRUE
          IF( LTIME )THEN ! Time-accurate
             CALL VecGetArray(DT,x_array,i_x,IFAIL)
             DO 14 IPOIN = 1, NPOIN ! May not be ok in mpi !
                ROWBGN = (IPOIN-1)*NOFVAR
                DO 14 IVAR = 1,NOFVAR
                   x_array(i_x+IVAR+ROWBGN) = VMEDIAN(IPOIN)/DELT
   14     CONTINUE
             CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          ELSE
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecSet(DT,DTMAX,IFAIL)
          ENDIF
      ENDIF
C
C     add V/(CNST*Dt) to the diagonal block of the stiffness matrix
C
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||dt|| before VecScale is ',s
#endif
      CALL VecScale(DT,ONE/CNST,IFAIL)
C
      IF( TIMEIMPL )THEN
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||dt|| before MatDiagonalSet is ',s
#endif
          CALL MatDiagonalSet(A,DT,ADD_VALUES,IFAIL)
#ifdef DEBUG
          call MatNorm(A,NORM_FROBENIUS,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||A|| after MatDiagonalSet is ',s
#endif
      ELSE
          CALL VecPointwiseDivide(RHS,RHS,DT,IFAIL)
      ENDIF
C
      IF( TIMEIMPL )THEN
C
          CALL JacobianBoundaryConditions(-1,ZROE,A,NDIM,NOFVAR)
#ifdef DEBUG
          call MatNorm(A,NORM_FROBENIUS,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||A|| after JacBCs is ',s
#endif
C
      ENDIF
C
#ifndef INITIAL_GUESS_NONZERO
C
      IF( TIMEIMPL )THEN
C
C     create a vector X to store the solution
C
          N = NPOIN*NOFVAR
#ifdef MPI
          CALL VecCreateMPI(PETSC_COMM_WORLD,N,PETSC_DECIDE,X,IFAIL)
#else
          CALL VecCreateSeq(PETSC_COMM_SELF,N,X,IFAIL)
#endif
          CALL VecSetBlockSize(X,NOFVAR,IFAIL)
C
          call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-rowScale',
     &         flg,IFAIL)
          if (flg .EQV. PETSC_TRUE) then
             N = NPOIN*NOFVAR
#ifdef MPI
             CALL VecCreateMPI(PETSC_COMM_WORLD,N,PETSC_DECIDE,LEFT,
     &                         IFAIL)
#else
             CALL VecCreateSeq(PETSC_COMM_SELF,N,LEFT,IFAIL)
#endif
             CALL VecSetBlockSize(LEFT,NOFVAR,IFAIL)
!         call MatGetRowMax(A,LEFT,IFAIL)
!         CALL VecReciprocal(LEFT,IFAIL)
             CALL VecGetArray(LEFT,x_array,i_x,IFAIL)
             CALL FSCALE(x_array(i_x+1),N,NOFVAR)
             CALL VecRestoreArray(LEFT,x_array,i_x,IFAIL)
             CALL MatDiagonalScale(A,LEFT,PETSC_NULL,IFAIL)
             CALL VecPointwiseMult(RHS,LEFT,RHS,IFAIL)
             CALL VecDestroy(LEFT,IFAIL)
          endif
C
C     Check if we are using ml as precondiioner
C
          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-pc_type',
     &    string,flg,IFAIL)
          IF(flg.EQV.PETSC_TRUE)THEN
               IF(string(1:2).EQ."ml")THEN
                  flg = PETSC_TRUE
               ELSE
                  flg = PETSC_FALSE
               ENDIF
          ENDIF
          IF(flg.EQV.PETSC_TRUE)THEN
          CALL MatConvert(A,MATAIJ,MAT_INITIAL_MATRIX,AijMat,IFAIL)
caldo     CALL MatGetType(AijMat,mt,IFAIL)
caldo     WRITE(6,*)'Mat type is ',mt,' on PE # ',MY_PE
          ENDIF
C
C     Sets the matrix associated with the linear system
C
          CALL PetscGetTime(tbeg,IFAIL)
C
          IF(flg.EQV.PETSC_TRUE)THEN
             CALL KSPSetOperators(FlowSolver,AijMat,AijMat,
     +       SAME_NONZERO_PATTERN,IFAIL)
          ELSE
             CALL KSPSetOperators(FlowSolver,A,A,SAME_NONZERO_PATTERN,
     +       IFAIL)
          ENDIF
C
#if 0
C
C           Here we explicitly call KSPSetUp() and KSPSetUpOnBlocks() to
C           enable more precise profiling of setting up the preconditioner.
C           These calls are optional, since both will be called within
C           KSPSolve() if they haven't been called already.
C
          CALL KSPSetUp(FlowSolver,IFAIL)
          CALL KSPSetUpOnBlocks(FlowSolver,IFAIL)
          CALL KSPGetPC(FlowSolver,mypc,IFAIL)
          CALL PCASMGetSubKSP(mypc,PETSC_NULL_INTEGER,
     &    PETSC_NULL_INTEGER,my_sub_ksp,IFAIL)
          CALL KSPGetPC(my_sub_ksp,my_sub_pc,IFAIL)
          CALL PCGetFactoredMatrix(my_sub_pc,factored_mat,IFAIL)
!         CALL PCASMGetLocalSubmatrices(mypc,1,factored_mat,IFAIL)
          CALL MatValid(factored_mat,matflag,IFAIL)
!         WRITE(6,*)'Factored Matrix is valid? ',matflag,PETSC_TRUE,
!    +' on PE # ',MY_PE
          IF(matflag.EQ.PETSC_FALSE)
     &    CALL MPI_Abort(PETSC_COMM_WORLD,-1,IFAIL)
          LIWORK = 1
          LJWORK = 1
          LDWORK = 1
          CALL PrintMatCSR(factored_mat,RHS,X,ISTAK(LIWORK),
     & ISTAK(LJWORK),DSTAK(LDWORK),NR,NZR,ITER,0)
          LIWORK = ISTKGT(NR+1,KIND_INTEGER)
          LJWORK = ISTKGT(NZR,KIND_INTEGER)
          LDWORK = ISTKGT(NZR,KIND_REAL8)
          CALL PrintMatCSR(factored_mat,RHS,X,ISTAK(LIWORK),
     & ISTAK(LJWORK),DSTAK(LDWORK),NR,NZR,ITER,100)
          CALL ISTKRL(3)
C
#endif
#if 1
!         LIWORK = 1
!         LJWORK = 1
!         LDWORK = 1
!         CALL PrintMatCSR(A,RHS,X,ISTAK(LIWORK),
!    & ISTAK(LJWORK),DSTAK(LDWORK),NR,NZR,ITER,0)
!         LIWORK = ISTKGT(NR+1,KIND_INTEGER)
!         LJWORK = ISTKGT(NZR,KIND_INTEGER)
!         LDWORK = ISTKGT(NZR,KIND_REAL8)
!         CALL PrintMatCSR(A,RHS,X,ISTAK(LIWORK),
!    & ISTAK(LJWORK),DSTAK(LDWORK),NR,NZR,ITER,100)
!         CALL ISTKRL(3)
#endif
C
C
C     Solve THE linear system
C
          CALL KSPSolve(FlowSolver,RHS,X,IFAIL)
          CALL KSPGetIterationNumber(FlowSolver,ITS,IFAIL)
#ifdef PRINT_KSP
          LIWORK = 1
          LDWORK = 1
!         CALL PrintMat(A,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),NITEMS,ITER,
!    & 0)
          CALL PrintMat(factored_mat,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),
     & NITEMS,ITER,0)
          LIWORK = ISTKGT(NITEMS,KIND_INTEGER)
          LDWORK = ISTKGT(NITEMS,KIND_REAL8)
!         CALL PrintMat(A,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),NITEMS,ITER,
!    & 100)
          CALL PrintMat(factored_mat,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),
     & NITEMS,ITER,100)
          CALL ISTKRL(2)
#endif
          CALL PetscGetTime(tend,IFAIL)
C
          IF(flg.EQV.PETSC_TRUE)CALL MatDestroy(AijMat,IFAIL)
C
      ELSE
C
C     explicit time stepping
C
          X = RHS
C
      ENDIF
#else
C
      IF( TIMEIMPL )THEN
          IF( ITER .EQ. 1 )THEN
C
C     create a vector X to store the solution
C
              N = NPOIN*NOFVAR
#ifdef MPI
              CALL VecCreateMPI(PETSC_COMM_WORLD,N,PETSC_DECIDE,X,IFAIL)
#else
              CALL VecCreateSeq(PETSC_COMM_SELF,N,X,IFAIL)
#endif
              CALL VecSetBlockSize(X,NOFVAR,IFAIL)
C
          ENDIF
C
          call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-rowScale',
     &         flg,IFAIL)
          if (flg .EQV. PETSC_TRUE) then
             N = NPOIN*NOFVAR
#ifdef MPI
          CALL VecCreateMPI(PETSC_COMM_WORLD,N,PETSC_DECIDE,LEFT,IFAIL)
#else
          CALL VecCreateSeq(PETSC_COMM_SELF,N,LEFT,IFAIL)
#endif
          CALL VecSetBlockSize(LEFT,NOFVAR,IFAIL)
          call MatGetRowMax(A,LEFT,IFAIL)
          CALL VecReciprocal(LEFT,IFAIL)
          CALL MatDiagonalScale(A,LEFT,PETSC_NULL,IFAIL)
          CALL VecPointwiseMult(RHS,LEFT,RHS,IFAIL)
          CALL VecDestroy(LEFT,IFAIL)
          endif
C
C     Sets the matrix associated with the linear system
C
          CALL PetscGetTime(tbeg,IFAIL)
          CALL KSPSetOperators(FlowSolver,A,A,SAME_NONZERO_PATTERN,
     +    IFAIL)
          call KSPSetUp(FlowSolver,IFAIL)
          call KSPSetInitialGuessNonzero(FlowSolver,flag,IFAIL)
C
C
C     Solve THE linear system
C
          CALL KSPSolve(FlowSolver,RHS,X,IFAIL)
          CALL KSPGetIterationNumber(FlowSolver,ITS,IFAIL)
#ifdef PRINT_KSP
          LIWORK = 1
          LDWORK = 1
          CALL PrintMat(A,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),NITEMS,ITER,
     & 0)
          LIWORK = ISTKGT(NITEMS,KIND_INTEGER4)
          LDWORK = ISTKGT(NITEMS,KIND_REAL8)
          CALL PrintMat(A,RHS,X,ISTAK(LIWORK),DSTAK(LDWORK),NITEMS,ITER,
     & 100)
          CALL ISTKRL(2)
#endif
          CALL PetscGetTime(tend,IFAIL)
 
C
      ELSE
C
C     explicit time stepping
C
          X = RHS
C
      ENDIF
#endif
#ifdef DO_NOT_UPDATE
      write(6,*)'Not updating the mean flow eqns'
      GOTO 456
#endif
C
C     Update the nodal unknown vector by forming Z := X + Z
C     REM: observe that ghost values are NOT updated by VecAXPY
C     therefore the X array has incorrect values within the
C     ghost locations; these will be corrected by a subsequent call
C     to RHSFunction() but, if the program terminates, because
C     convergence or max its has been reached, the datafile will
C     contain incorrect ghost values
C
      CALL VecAXPY(ZROE,ONE,X,IFAIL)
  456 continue
C
C     Monitor the L2 and L infinity norms of the update ..
C
      CALL FindVecStrideAbsMinMax(X,DELL2(1,1),DELMAX(1,1),INDEL(1,1),
     &WKSP2,NDIM)
      CALL FindVecStrideMinMax(ZROE,AL2(1,5),AMAX(1,5),WHEREMAX(1,5),
     &WKSP5,NDIM)
C
C     Print out the convergence history ..
C
      IF(MY_PE.EQ.0)THEN
          IF ((ITER/ISTMP)*ISTMP.EQ.ITER) THEN
              WRITE (IWUNIT,FMT=200) ITER,ITS,CNST
              WRITE (IWUNIT,FMT=215)
              DO 20 IVAR = 1,NOFVAR
                  WRITE (IWUNIT,FMT=210) IVAR,
     &            DLOG10(MAX(EPS,RESL2(IVAR,1))),
     +            DLOG10(MAX(EPS,RESMAX(IVAR,1))),
     &            INMAX(IVAR,1), (WKSP1(I,IVAR),I=1,3)
   20         CONTINUE
              WRITE (IWUNIT,FMT=225)
              DO 32 IVAR = 1,NOFVAR
                  WRITE (IWUNIT,FMT=210) IVAR,
     &            DLOG10(MAX(EPS,DELL2(IVAR,1))),
     +            DLOG10(MAX(EPS,DELMAX(IVAR,1))),
     &            INDEL(IVAR,1), (WKSP2(I,IVAR),I=1,3)
   32         CONTINUE
              WRITE (IWUNIT,FMT=240)
              DO 34 IVAR = 1,NOFVAR
                  WRITE (IWUNIT,FMT=212) IVAR,
     &            AL2(IVAR,4), AMAX(IVAR,4),
     &            WHEREMAX(IVAR,4), (WKSP4(I,IVAR),I=1,3)
   34         CONTINUE
              WRITE (IWUNIT,FMT=230)
              DO 36 IVAR = 1,NOFVAR
                  WRITE (IWUNIT,FMT=212) IVAR,
     &            AL2(IVAR,5), AMAX(IVAR,5),
     &            WHEREMAX(IVAR,5), (WKSP5(I,IVAR),I=1,3)
   36         CONTINUE
          ENDIF
      ENDIF
C
C     Write convergence history to file ...
C
      IF(MY_PE.EQ.0)THEN
          CALL PetscGetTime(telapsed,IFAIL)
          telapsed=telapsed-tbegall
          WRITE (IHST1,FMT=235) NITER,ITS,tend-tbeg,telapsed, (RESL2
     +    (IVAR,1),IVAR=1,NOFVAR),CNST
          WRITE (IHST2,FMT=235) NITER,ITS,tend-tbeg,telapsed, (RESMAX
     +    (IVAR,1),IVAR=1,NOFVAR),CNST
      ENDIF
C
C     Check to make sure that no negative viscosities appear
C
#ifdef NOWHERE
      IF(TURBULENT.AND.COUPLED)THEN
          THRESH = -1.D-2
          THRESH = -.5D+0
          THRESH = ZERO
          IF(AL2(NOFVAR,5) .LT. THRESH )THEN
!     CALL solzne("res010.dat",x_array(i_x+1),nofvar,npoin,"w")
!             OMEGA = OMEGA*HALF
              IF(MY_PE.EQ.0)THEN
              write(IWUNIT,*)'min(Z) is now = ',AL2(NOFVAR,5)
              write(IWUNIT,*)'Omega is now set = ',omega
              ENDIF
!             CALL VecGetArray(ZROE,x_array,i_x,IFAIL)
!             NITEMS = (NPOIN)*NOFVAR
!             DO i = NOFVAR, NITEMS, NOFVAR
!                  x_array(i_x+i) = max(ZERO,x_array(i_x+i))
!             ENDDO
!             CALL VecRestoreArray(ZROE,x_array,i_x,IFAIL)
              CALL VecAXPY(ZROE,MONE,X,IFAIL)
          ELSE
              OMEGA = ONE
          ENDIF
      ENDIF
#endif
C
#ifndef INITIAL_GUESS_NONZERO
C
C     clear memory allocated for the solution vector
C
      IF(TIMEIMPL)CALL VecDestroy(X,IFAIL)
#endif
C
C     We do not need DT any more
C
      CALL VecDestroy(DT,IFAIL)
C
      RETURN
 
  200 FORMAT (5X, 70('-'),/,25X,'ITERATION # ',I4,' (',I4,') CFL = ',E10
     +.4,/,5X, 70('-'),/,5X, 70('-'),/,5X,'Var.',4X,'L2-norm',3X,
     +'L_infty',3X,'node #',3X,'(',4X,'x',7X,'y',7X,'z',4X,')',/, 5X, 70
     +('-'))
  210 FORMAT (5X,I1,5X,F10.5,1X,F10.5,1X,I6,3X,'(', 2(F8.4,','),F8.4,
     +')')
  212 FORMAT (5X,I1,5X,E10.3,1X,E10.3,1X,I6,3X,'(', 2(F8.4,','),F8.4,
     +')')
  215 FORMAT (5X,'Nodal Residual',/)
  225 FORMAT (5X,'Nodal Update',/)
  230 FORMAT (5X,'Nodal Values',/)
  235 FORMAT (I5,1X,I4,9 (1X,E10.4))
  240 FORMAT (5X,'Timestep ',/)
      END
 
 
 
      subroutine mysillyroutine(dt,NOFVAR,NOFEQN,NPOIN)
C
      IMPLICIT NONE
C
      integer nofvar,nofeqn,npoin
      double precision dt(nofvar,*)
      double precision v
      integer i,j

      do i = 1,npoin
         v = dt(1,i)
         do j = 2,nofeqn
         if(dt(j,i).ne.0.d0)then
            write(6,*)'silly',i,j,dt(j,i)
         endif
         dt(j,i) = v
         enddo
         if(nofvar.gt.nofeqn)then
             if(ABS( dt(nofvar,i) ) .LT. 1.E-12 )dt(nofvar,i)=v
         endif 
!        dt(nofvar,i)=v
!        write(40,*)'mesh point ',i
!        write(40,FMT="(4(E12.3,1X))")(dt(j,i),j=1,nofvar)
      enddo
      return
      end

      SUBROUTINE FindVecStrideAbsMinMax(X,ANL2,ANMAX,INMAX,WKSP,NDIM)
C
      IMPLICIT NONE
C
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
C
      INCLUDE 'nloc'
C
      INTEGER NDIM
      Vec X
      DOUBLE PRECISION ANMAX(*),ANL2(*),WKSP(3,*)
      INTEGER INMAX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION XMIN,XMAX
      INTEGER I,IVAR,IFAIL,MY_PE,ROWBGN,ROWEND,IPOIN,INDX,
     &IADDR,INDXMIN,INDXMAX,NOFVAR,IOFF
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON /MPICOM/MY_PE
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
      CALL VecGetOwnerShipRange(X,ROWBGN,ROWEND,IFAIL)
      CALL VecGetBlockSize(X,NOFVAR,IFAIL)
      DO 24 IVAR = 1,NOFVAR
         CALL VecStrideNorm(X,IVAR-1,NORM_2,ANL2(IVAR),IFAIL)
         CALL VecStrideMin(X,IVAR-1,INDXMIN,XMIN,IFAIL)
         CALL VecStrideMax(X,IVAR-1,INDXMAX,XMAX,IFAIL)
         IF( ABS(XMAX) .GT. ABS(XMIN) )THEN
             INDX = INDXMAX 
             ANMAX(IVAR) = ABS(XMAX)
         ELSE
             INDX = INDXMIN
             ANMAX(IVAR) = ABS(XMIN)
         ENDIF
         IF(INDX.NE.-1)THEN
            IPOIN = (INDX + 1 - IVAR)/NOFVAR + 1
         ELSE 
            IPOIN = INDX
         ENDIF
         INMAX(IVAR) = IPOIN
#ifndef MPI
         IF( ROWBGN .LE. INDX .AND. INDX .LE. ROWEND)THEN
             IADDR = LCORG + (IPOIN-1)*NDIM
             CALL DCOPY(NDIM,DSTAK(IADDR),1,WKSP(1,IVAR),1)
         ENDIF
#endif
   24 CONTINUE
      RETURN
      END
C
      SUBROUTINE FindVecStrideMinMax(X,ANMIN,ANMAX,INMAX,WKSP,NDIM)
C
C     returns min/max values of an array along with
C     the location where the max ABS value occurs
C
      IMPLICIT NONE
C
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
C
      INCLUDE 'nloc'
C
      INTEGER NDIM
      Vec X
      DOUBLE PRECISION ANMAX(*),ANMIN(*),WKSP(3,*)
      INTEGER INMAX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION XMIN,XMAX
      INTEGER I,IVAR,IFAIL,MY_PE,ROWBGN,ROWEND,IPOIN,INDX,
     &IADDR,INDXMIN,INDXMAX,NOFVAR
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON /MPICOM/MY_PE
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
      CALL VecGetOwnerShipRange(X,ROWBGN,ROWEND,IFAIL)
      CALL VecGetBlockSize(X,NOFVAR,IFAIL)
      DO 24 IVAR = 1,NOFVAR
         CALL VecStrideMin(X,IVAR-1,INDXMIN,ANMIN(IVAR),IFAIL)
         CALL VecStrideMax(X,IVAR-1,INDXMAX,ANMAX(IVAR),IFAIL)
C
C     get the nodenumber where the minimum occurs
C
         INDX = INDXMIN
         IF(INDX.NE.-1)THEN ! fixes a PETSc bug ?!?!?!
            IPOIN = (INDX + 1 - IVAR)/NOFVAR + 1
         ELSE 
            IPOIN = INDX
         ENDIF
         INMAX(IVAR) = IPOIN
#ifndef MPI
         IF( ROWBGN .LE. INDX .AND. INDX .LE. ROWEND)THEN
             IADDR = LCORG + (IPOIN-1)*NDIM
             CALL DCOPY(NDIM,DSTAK(IADDR),1,WKSP(1,IVAR),1)
         ENDIF
#endif
   24 CONTINUE
      RETURN
      END
      SUBROUTINE FSCALE(X,N,BS)
      IMPLICIT NONE
      INTEGER N,BS
      DOUBLE PRECISION X(N)
      INTEGER I
      DOUBLE PRECISION HELP
      INCLUDE "paramt.h" 
      INCLUDE "constants.h" 
      INCLUDE "stream.com" 
      CALL DINIT(N,ONE,X,1)
      HELP = M_INFTY*M_INFTY
      CALL DINIT(N/BS,HELP,X(BS),BS)
!     do i = 1,n
!        write(17,*)i,x(i)
!     enddo
!     call exit(3)
      RETURN
      END
