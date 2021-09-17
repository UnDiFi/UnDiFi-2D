!>
!> \brief
!> updates the nodal solution by implicit time integration when KAN = 2,4
!>
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVAR number of dofs and leading dimension of U
!> @param[in] NPOIN number of gridpoints and second dimension of U
!> @param[in] FlowSolver is the PETSc KSP used to solve the sparse linear system
!> @param[in,out] A (PETSc mat) is the jacobian matrix
!> @param[in,out] RHS (PETSc vec) nodal residual
!> @param[in,out] DT (PETSc vec) nodal timestep
!> @param[in] U (PETSc vec) the vector of conserved variables
!> @param[in] VMDCO the median dual cell area/volume at time level n
!> @param[in] VMDCN the median dual cell area/volume at time level n+1
!> @param[in] NodalBcs array of Index Sets addressing the nodal bcs
!> \author $Author: abonfi $
!> \version $Revision: 1.84 $
!> \date $Date: 2020/04/23 09:46:19 $
!> \todo Should rather destroy DT in the calling routine, where it is created
!> \todo Duplicate call to \c RHSBC4 in the \c ALE case should be changed/optimized
      SUBROUTINE UPDATE3(NDIM,NOFVAR,NPOIN,FlowSolver,A,RHS,DT,U,
     &                   VMDCO,VMDCN,NodalBcs)
C
C     $Id: update3.F,v 1.84 2020/04/23 09:46:19 abonfi Exp $
C
C#define INITIAL_GUESS_NONZERO
CC#define PRINT_KSP
CCC#define DEBUG
CC#define DO_NOT_UPDATE
CCC#define NOWHERE
CCC#define HASISNAN
C
C
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscvec.h"
      use petscvec
      use petscpc
      use petscksp
C
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'paramt.h'
      INTEGER NDNM
      PARAMETER (NDNM=3*MAXNOFVAR)
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.h'
      INCLUDE 'datatype.h'
      INCLUDE 'iset.h'
C
      INCLUDE 'time.com'
      INCLUDE 'conv.com'
      INCLUDE 'visco.com'
      INCLUDE 'nloc.com'
      INCLUDE 'verbose.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
C
C
      COMMON /TIMING/TBEGALL
C
C     Petsc stuff
C
      Mat A
      Vec RHS,DT,U,X,LEFT
      SAVE X
      KSP FlowSolver
      PetscLogDouble TBEGALL,telapsed,tbeg,tend
      DOUBLE PRECISION VMDCO(*),VMDCN(*) ! Volume of the Median Dual cell New=(n+1) Old=n
C
      Mat AijMat
      MatType mt
C
caldo
      PC mypc
      PC my_sub_pc
      KSP my_sub_ksp
      Mat factored_mat
      PetscBool matflag
      PetscBool flg
      IS, dimension(0:*) :: NodalBcs
      DOUBLE PRECISION INFO(MAT_INFO_SIZE)
caldo
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR,NPOIN
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
#ifdef NOWHERE
      DOUBLE PRECISION THRESH 
#endif
      DOUBLE PRECISION CNST,DTMAX,eps,S
      parameter(eps=1.d-20) 
      INTEGER I,IVAR,N,ITS,IFAIL,MY_PE,ROWBGN,ROWEND,IPOIN,NI,INDX,
     &IADDR,INDXMIN,INDXMAX,LIWORK(3),LDWORK,NITEMS,NR,NZR,BS,
     2IBGN,IEND,NofHangingNodes
      CHARACTER*23 MYFMT
      CHARACTER*11 fname
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP1(3,MAXNOFVAR),WKSP2(3,MAXNOFVAR),
     &  WKSP3(3,MAXNOFVAR),WKSP4(3,MAXNOFVAR),
     &  WKSP5(3,MAXNOFVAR), AL2(MAXNOFVAR,5),AMAX(MAXNOFVAR,5)
      INTEGER WHEREMAX(MAXNOFVAR,5)
      INTEGER ISTAK(1)
      PetscScalar x_array(1)
      PetscScalar rhs_array(1)
      PetscScalar z_array(1)
      INTEGER idx_v(1)
      PetscOffset idx_i,i_x,i_rhs,i_z
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
C 
C     Executable statements 
C
      MYFMT = "(I5,1X,I4,10(1X,E10.4))"
      WRITE(MYFMT(11:12),FMT="(I2.2)")NOFVAR+3
C
#ifdef DEBUG
      do ivar = 1,nofvar
      CALL VecStrideNorm(U,ivar-1,NORM_2,s,IFAIL)
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
             CNST = RESL20(1)/RESL2(IVCNVG,1) * OMEGA
             CNST = MIN(CFLMAX(1),CFL(1)*CNST)
!            write(6,*)'Omega is = ',omega,CFLMAX(1),CFL(1),cnst
!           write(6,*)'iter is = ',iter,RESL20(1),RESL2(IVCNVG,1),ivcnvg
          ELSEIF( CFL_RAMP .EQ. CFL_EXP )THEN
             CNST = REAL(ITER-1)/REAL(ITMAX-1)
             CNST = CFL(1)+((EXP(PE*CNST)-ONE)/(EXP(PE)-ONE))*
     &              (CFLMAX(1)-CFL(1))
          ENDIF
      ELSE
!         CNST=ONE/CFL(1)
          CNST=CFL(1)
      ENDIF
C
C     dump nodal residual when using -dump_nodal_residual
C     for debugging purposes
C
      IF(LDUMP(1))THEN
         fname = 'rhs0000.dat'
         WRITE(fname(4:7),FMT="(I4.4)")ITER
         CALL VecGetArray(rhs,x_array,i_x,IFAIL)
         CALL solzne(fname,x_array(i_x+1),nofvar,npoin,"w")
         CALL VecRestoreArray(rhs,x_array,i_x,IFAIL)
      ENDIF
C
C     For pure advection problems, the inverse of
C         the local timestep in the Dirichlet Nodes
C         is 0., so we set it to the maximum DT
C
C to be changed
C
      IF(KAN.EQ.-4.OR.KAN.EQ.-2.OR.KAN.EQ.-3)THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,x_array,i_x,IFAIL)
          CALL ISGetIndices(NodalBcs(SupersonicVariables),IDX_V,IDX_I,
     &                      IFAIL)
          CALL ISGetSize(NodalBcs(SupersonicVariables),NI,IFAIL)
          DO 10 I = 1, NI
C
C     supersonic nodes are 0-based indexed
C
              IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
C             write(6,*)x_array(i_x+IPOIN) , DTMAX
              x_array(i_x+IPOIN) = DTMAX
   10     CONTINUE
          CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          CALL ISRestoreIndices(NodalBcs(SupersonicVariables),IDX_V,
     &                          IDX_I,IFAIL)
      ENDIF
C
#ifdef DEBUG
      IF(TIMEIMPL)THEN
         call MatNorm(A,NORM_FROBENIUS,s,ifail)
         if(MY_PE.EQ.0)write(6,*)'||A|| before is ',s
      ENDIF
         call VecNorm(DT,NORM_2,s,ifail)
         if(MY_PE.EQ.0)write(6,*)'||dt|| before myroutine is ',s
#endif
C
C     Handles hanging nodes
C
      CALL ISGetSize(NodalBcs(HangingNodes),NofHangingNodes,IFAIL)
      IF( NofHangingNodes .GT. 0 .AND. IGLOB .EQ. 1 )THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,x_array,i_x,IFAIL)
          CALL ISGetIndices(NodalBcs(HangingNodes),IDX_V,IDX_I,IFAIL)
          DO 12 I = 1, NofHangingNodes !, NOFVAR
C     hanging nodes are 0-based indexed
              IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
              x_array(i_x+IPOIN) = DTMAX
   12     CONTINUE
          CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          CALL ISRestoreIndices(NodalBcs(HangingNodes),IDX_V,IDX_I,
     &                          IFAIL)
      ENDIF
C
C     time V/dt to file, if required
C
      IF(LDUMP(2))THEN
         CALL VecGetArray(DT,x_array,i_x,IFAIL)
         fname = 'dtv0000.dat'
         write(fname(4:7),FMT="(I4.4)")ITER
         CALL solzne(fname,x_array(i_x+1),nofvar,npoin,"w")
         CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
      ENDIF
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
          IF( LTIME )THEN ! Time-accurate
             CALL VecGetArray(DT,x_array,i_x,IFAIL)
             DO 14 IPOIN = 1, NPOIN ! May not be ok in mpi !
                IADDR = (IPOIN-1)*NOFVAR
                DO 14 IVAR = 1,NOFVAR
                   x_array(i_x+IVAR+IADDR) = VMDCO(IPOIN)/DELT
   14        CONTINUE
             IF(NofHangingNodes.NE.0)THEN
                CALL ISGetIndices(NodalBcs(HangingNodes),IDX_V,IDX_I,
     &                            IFAIL)
                DO 15 I = 1, NofHangingNodes !, NOFVAR
C     hanging nodes are 0-based indexed
                   IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
                   write(6,*)i,ipoin
                   x_array(i_x+IPOIN) = ONE
   15           CONTINUE
                CALL ISRestoreIndices(NodalBcs(HangingNodes),IDX_V,
     &                                IDX_I,IFAIL)
             ENDIF
             CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          ELSE ! using Local time-stepping
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
C     create a vector X to store the solution
C
          N = NPOIN*NOFVAR
          CALL VecCreate(PETSC_COMM_WORLD,X,IFAIL)
#ifdef MPI
          CALL VecSetType(X,VECMPI,IFAIL)
#else
          CALL VecSetType(X,VECSEQ,IFAIL)
#endif
          CALL VecSetBlockSize(X,NOFVAR,IFAIL)
          CALL VecSetSizes(X,N,PETSC_DECIDE,IFAIL)
C
C     Apply boundary conditions
C
#ifdef DEBUG
          call VecNorm(X,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||X|| before JacBCs is ',s
          call VecNorm(RHS,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||B|| before JacBCs is ',s
#endif
C
          CALL JacobianBoundaryConditions(-1,U,A,X,RHS,NodalBcs,
     &                                    NDIM,NOFVAR)
C
#ifdef DEBUG
          call VecNorm(X,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||X|| after JacBCs is ',s
          call VecNorm(RHS,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||B|| after JacBCs is ',s
#endif
C
#ifdef DEBUG
          call MatNorm(A,NORM_FROBENIUS,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||A|| after JacBCs is ',s
#endif
C
C     Sets the matrix associated with the linear system
C
          CALL PetscTime(tbeg,IFAIL)
C
          CALL KSPSetOperators(FlowSolver,A,A,IFAIL)
C
C   0 maxes the lines INvisible to the preprocessor
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
          CALL PCFactorGetMatrix(my_sub_pc,factored_mat,IFAIL)
          CALL PCASMGetLocalSubmatrices(mypc,1,factored_mat,IFAIL)
          CALL MatValid(factored_mat,matflag,IFAIL)
          WRITE(6,*)'Factored Matrix is valid? ',matflag,PETSC_TRUE,
     +' on PE # ',MY_PE
          IF(matflag.EQ.PETSC_FALSE)
     &    CALL MPI_Abort(PETSC_COMM_WORLD,-1,IFAIL)
          LIWORK(1) = 1
          LIWORK(2) = 1
          LDWORK = 1
          CALL PrintMatCSR(A,RHS,X,ISTAK(LIWORK(1)),
     & ISTAK(LIWORK(2)),DSTAK(LDWORK),NR,NZR,BS,ITER,0)
          LIWORK(1) = ISTKGT(NR+1,KIND_INTEGER)
          LIWORK(2) = ISTKGT(NZR,KIND_INTEGER)
          LDWORK = ISTKGT(NZR,KIND_REAL8)
          CALL PrintMatCSR(factored_mat,RHS,X,ISTAK(LIWORK(1)),
     & ISTAK(LIWORK(2)),DSTAK(LDWORK),NR,NZR,BS,ITER,100)
          CALL ISTKRL(3)
C
#endif
C
C     Solve THE linear system
C
          CALL KSPSolve(FlowSolver,RHS,X,IFAIL)
          CALL KSPGetIterationNumber(FlowSolver,ITS,IFAIL)
C
C     dump the jacobian matrix: -dump_jacobian_matrix
C
      IF(LDUMP(3))THEN
          LIWORK(1) = 1
          LIWORK(2) = 1
          LIWORK(3) = 1
          LDWORK = 1
          NR = NPOIN
          CALL PrintMatMM(A,RHS,X,
     3                    ISTAK(LIWORK(3)),NR,NZR,BS,ITER,0)
          LIWORK(3) = ISTKGT(NZR,KIND_INTEGER)! storage for ir
          CALL PrintMatMM(A,RHS,X,
     3                    ISTAK(LIWORK(3)),NR,NZR,BS,ITER,100)
          write(6,*)'Beyond PrintMatMM'
          CALL ISTKRL(1)
      ENDIF
C   0 maxes the lines INvisible to the preprocessor
#if 0 
      IF(LDUMP(3))THEN
          LIWORK(1) = 1
          LIWORK(2) = 1
          LDWORK = 1
          NR = NPOIN
          CALL PrintMatCSR(A,RHS,X,ISTAK(LIWORK(1)),
     & ISTAK(LIWORK(2)),DSTAK(LDWORK),NR,NZR,BS,ITER,0)
          LIWORK(1) = ISTKGT(NR+1,KIND_INTEGER)
          LIWORK(2) = ISTKGT(NZR,KIND_INTEGER)
          LDWORK = ISTKGT(NZR,KIND_REAL8)
          CALL PrintMatCSR(A,RHS,X,ISTAK(LIWORK(1)),
     & ISTAK(LIWORK(2)),DSTAK(LDWORK),NR,NZR,BS,ITER,100)
          CALL ISTKRL(3)
      ENDIF
#endif
C
#ifdef PRINT_KSP
          LIWORK(1) = 1
          LDWORK = 1
          CALL PrintMat(A,RHS,X,ISTAK(LIWORK(1)),DSTAK(LDWORK),
     &                  NITEMS,ITER,0)
caldo     CALL KSPGetPC(FlowSolver,mypc,IFAIL)
caldo     CALL PCFactorGetMatrix(mypc,factored_mat,IFAIL)
caldo     call MatGetInfo(factored_mat,MAT_LOCAL,info,ifail)
caldo     do i = 1,10
caldo        write(6,*)'item = ',i,info(i)
caldo     enddo
!         CALL PrintMat(factored_mat,RHS,X,ISTAK(LIWORK(1)),DSTAK(LDWORK),
!    & NITEMS,ITER,0)
          LIWORK(1) = ISTKGT(NITEMS,KIND_INTEGER)
          LDWORK = ISTKGT(NITEMS,KIND_REAL8)
          CALL PrintMat(A,RHS,X,ISTAK(LIWORK(1)),DSTAK(LDWORK),
     &                  NITEMS,ITER,100)
!         CALL PrintMat(factored_mat,RHS,X,ISTAK(LIWORK(1)),DSTAK(LDWORK),
!    & NITEMS,ITER,100)
          CALL ISTKRL(2)
#endif
          CALL PetscTime(tend,IFAIL)
C
C
      ELSE ! explicit time-stepping
C
C     explicit time stepping
C       
         IF(LALE)THEN
            CALL VecGetArray(rhs,rhs_array,i_rhs,IFAIL)
            CALL VecGetArray(U,z_array,i_z,IFAIL)
C
C   here we have to reset rhs=0 in the supersonic nodes
C
C     we should do this better, RHSBC4 has already been called within RHS_Function
            CALL VecGetOwnershipRange(rhs,ibgn,iend,IFAIL)
C
C     boundary conditions for the mean flow equations 
C
            CALL RHS4ALE(NOFVAR,NPOIN,rhs_array(i_rhs+1),
     &                   z_array(i_z+1),VMDCO,VMDCN)
            CALL RHSBC4(x_array(i_x+1),rhs_array(i_rhs+1),NodalBcs,
     &               ibgn,NDIM,DSTAK(LFREE),(ABS(KAN).EQ.4))
            CALL VecRestoreArray(U,z_array,i_z,IFAIL)
            CALL VecRestoreArray(rhs,rhs_array,i_rhs,IFAIL)
         END IF  ! ALE
         X=RHS
C
      ENDIF
C
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
      CALL VecAXPY(U,ONE,X,IFAIL)
  456 continue
C
C     Monitor the L2 and L infinity norms of the update ..
C
      CALL FindVecStrideAbsMinMax(X,DELL2(1,1),DELMAX(1,1),INDEL(1,1),
     &WKSP2,NDIM)
      CALL FindVecStrideMinMax(U,AL2(1,5),AMAX(1,5),WHEREMAX(1,5),
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
          CALL PetscTime(telapsed,IFAIL)
          telapsed=telapsed-tbegall
          WRITE (IHST1,FMT=MYFMT) NITER,ITS,tend-tbeg,telapsed, (RESL2
     +    (IVAR,1),IVAR=1,NOFVAR),CNST
          WRITE (IHST2,FMT=MYFMT) NITER,ITS,tend-tbeg,telapsed, (RESMAX
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
              CALL VecGetArray(U,x_array,i_x,IFAIL)
!             NITEMS = (NPOIN)*NOFVAR
!             DO i = NOFVAR, NITEMS, NOFVAR
!                  x_array(i_x+i) = max(ZERO,x_array(i_x+i))
!             ENDDO
              CALL LIMITVISCT(x_array(i_x+1),NOFVAR,1,NPOIN)
              CALL VecRestoreArray(U,x_array,i_x,IFAIL)
!             CALL VecAXPY(U,MONE,X,IFAIL)
!         ELSE
!             OMEGA = ONE
          ENDIF
      ENDIF
#endif
C
C
C     clear memory allocated for the solution vector
C
      IF(TIMEIMPL)CALL VecDestroy(X,IFAIL)
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
  240 FORMAT (5X,'Timestep ',/)
      END
!>
!>    \brief
!>    returns the L2-norm of an array along with
!>    the location where the max ABS value occurs
!> @param[in] X The Petsc vec
!> @param[out] ANL2 L2-norm of X
!> @param[out] ANMAX maximum value
!> @param[out] INMAX location where the maximum value occurs
!> @param[out] WKSP coordinates of the location where the maximum value occurs
!> @param[in] NDIM dimension of the space
!
!> \author $Author: abonfi $
!> \version $Revision: 1.84 $
!> \date $Date: 2020/04/23 09:46:19 $
!
!> \warning The coordinates are only returned when run on a single processor
      SUBROUTINE FindVecStrideAbsMinMax(X,ANL2,ANMAX,INMAX,WKSP,NDIM)
C
C
#include "petsc/finclude/petscvec.h"
      use petscvec
      IMPLICIT NONE
C
      INCLUDE 'nloc.com'
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
!>
!>    \brief
!>    returns min/max values of an array along with
!>    the location where the min/max values occur
!> @param[in] X The Petsc vec
!> @param[out] ANMIN minimum value
!> @param[out] ANMAX maximum value
!> @param[out] INMAX location where the maximum value occurs
!> @param[out] WKSP coordinates of the location where the maximum value occurs
!> @param[in] NDIM dimension of the space
!
!> \author $Author: abonfi $
!> \version $Revision: 1.84 $
!> \date $Date: 2020/04/23 09:46:19 $
!
!> \warning The coordinates are only returned when run on a single processor
!> \warning The location where the minimum value occurs is NOT returned
      SUBROUTINE FindVecStrideMinMax(X,ANMIN,ANMAX,INMAX,WKSP,NDIM)
C
C
#include "petsc/finclude/petscvec.h"
      use petscvec
      IMPLICIT NONE
C
C
      INCLUDE 'nloc.com'
C
      INTEGER NDIM
      Vec X
      DOUBLE PRECISION ANMAX(*),ANMIN(*),WKSP(3,*)
      INTEGER INMAX(*)
C     ..
C     .. Local Scalars ..
      INTEGER IVAR,IFAIL,MY_PE,ROWBGN,ROWEND,IPOIN,INDX,
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
!>
!>    \brief
!>    return the modified RHS if we are dealing with moving and/or deforming grids
!>
!>    the update for a two times levels scheme such as the Lax-Wendroff reads:
!>    \f$ K^{n+1} Z^{n+1}-K^{n} Z^{n} = R\left(Z^n\right) \f$, where \f$ K = V/ \Delta t \f$.
!>
!>    this routine updates \c RHS in such a way that, within the calling subroutine \c update3 we can later set
!>    \f$ U^{n+1} = U^n + \mbox{RHS} \f$
!>    note that upon entry \f$ \mbox{RHS} = \left(K^{n}\right)^{-1} R\left(Z^n\right) \f$
!>    rather than just being the nodal residual
!>
!>    This is achieved as follows:
!>    \f{eqnarray*}{
!>    RHS &=& \left(K^{n+1}\right)^{-1} \left[K^{n}U^{n}+R\left(U^n\right)\right] - U^n \\\
!>        &=&  \left\{ \left[ \left(K^{n+1}\right)^{-1}K^n-I\right]U^n+\left(K^{n+1}\right)^{-1}R\left(U^n\right)\right\} \\\
!>        &=&  \left\{ \left(K^{n+1}\right)^{-1}K^n \left[U^n+\left(K^n\right)^{-1}R\left(U^n\right)\right] - U^n \right\} \\\
!>        &=&  \left(K^{n+1}\right)^{-1}K^n \left(U^n + RHS\right) - U_n \\\
!>        &=&  \left[\left(K^{n+1}\right)^{-1}K^n -I\right]U^n + \left(K^{n+1}\right)^{-1} \,K^n\, RHS
!>    \f}
!>
!>    in practice, \f$ \left(K^{n+1}\right)^{-1} \,K^n = V^n/V^{n+1} \f$
!>
!>    when the grid is fixed, then \f$ (K^{n+1})^{-1}*K^n = I \f$
!>    and \f$ \mbox{RHS} \f$ reduces to: \f$ RHS = (K^{n})^{-1}*R(Z^n)  \f$
!>
!>    since on entry \f$ \mbox{DT} = K^n \f$ and \f$ \mbox{RHS} = (K^n)^{-1} R(Z^n) \f$
!>
!> @param[in] NOFVAR number of dofs and leading dimension of U
!> @param[in] NPOIN number of gridpoints and second dimension of U
!> @param[in,out] RHS on entry it is: \f$ \mbox{RHS} = \left(K^{n}\right)^{-1} R\left(Z^n\right) \f$
!> @param[in] U the vector of conserved variables
!> @param[in] VMDCO the median dual cell area/volume at time level n
!> @param[in] VMDCN the median dual cell area/volume at time level n+1
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.84 $
!> \date $Date: 2020/04/23 09:46:19 $
!>
!> \warning NPOIN should only include interior nodes only (no ghosts, no periodic ones)
!> \warning NOT sure it is going to work in parallel
!> \todo Could just pass the ratio between the Median Dual cells, rather than both
      SUBROUTINE RHS4ALE(NOFVAR,NPOIN,RHS,U,VMDCO,VMDCN)
C
C
      IMPLICIT NONE
C
      INTEGER IPOIN,NPOIN,NOFVAR
C
C
      DOUBLE PRECISION VMDCO(*),VMDCN(*),U(*),RHS(*)
      DOUBLE PRECISION HELP
      INTEGER IADDR,ROWBGN,IVAR,IFAIL
C     LOGICAL ISNAN
#define HASISNAN
#ifdef HASISNAN
      INTRINSIC ISNAN
#endif
C 
      INCLUDE "constants.h" !includes ONE
C
      IFAIL = 0
      DO IPOIN = 1,NPOIN ! loop over interior nodes only (no ghosts, no periodic ones)
         ROWBGN = (IPOIN-1)*NOFVAR
         HELP = VMDCO(IPOIN)/VMDCN(IPOIN) ! V^{n}/V^{n+1}
#ifdef HASISNAN
         IF(ISNAN(HELP))THEN
            HELP = ONE
            IFAIL = IFAIL + 1
         ENDIF
#endif
         DO IVAR = 1, NOFVAR
            IADDR = ROWBGN + IVAR 
!     on a unifrom flowfield, before updating, rhs = (V^{n+1}/V^{n}-1)*u
!           write(6,*)'b',ipoin,' rhs = ',rhs(iaddr),' should be = ',
!    & (one/help-one)*U(iaddr),' u = ',U(iaddr)
            RHS(IADDR) = HELP * RHS(IADDR) + (HELP-ONE) * U(IADDR) 
!     on a unifrom flowfield, after updating, rhs should be 0.d0
!           write(6,*)'a',ipoin,ivar,' rhs = ',rhs(iaddr),
!    &' U = ',U(iaddr)
         ENDDO
      ENDDO
      IF(IFAIL.NE.0)THEN
          WRITE(6,*)IFAIL,' NaNs detected in update3()'
      ENDIF
      RETURN
      END
