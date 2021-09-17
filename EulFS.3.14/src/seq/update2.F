      SUBROUTINE UPDATE2( NDIM, NOFVAR, NPOIN, solver, A, RHS, DT, Z,
     & VMDCO, VMDCN, NodalBcs )
C
C     $Id: update2.F,v 1.48 2020/04/23 09:46:19 abonfi Exp $
C
C     Implicit timestepping for scalar equations ..
C
CCCC#define DO_NOT_UPDATE
CCCC#define DEBUG
C
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscmat.h"
CCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscksp.h"
      use petscksp
C
      IMPLICIT NONE
C
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.h'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
      INTEGER NDNM
      PARAMETER(NDNM=3*NMAX)
C
      INCLUDE 'conv.com'
      INCLUDE 'nloc.com'
      INCLUDE 'verbose.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'time.com'
      COMMON /TIMING/TBEGALL
C
CCCC#define DEBUG
C
      INTEGER NDIM,NOFVAR,NPOIN,IFAIL,ITS,Number_of_DirichletNodes,
     &NITEMS
C
      INTEGER IVAR,ROWBGN,ROWEND,IPOIN,I,MY_PE,NI
      DOUBLE PRECISION CNST,S,DTMAX
C
      DOUBLE PRECISION VMDCO(*),VMDCN(*) ! Volume of the Median Dual cell New=(n+1) Old=n
      DOUBLE PRECISION WKSP1(3,NMAX),WKSP2(3,NMAX)
      INTEGER IDX_V(1)
      CHARACTER*11 fname
C
      DOUBLE PRECISION dnrm2
      EXTERNAL DNRM2
C
      DATA WKSP1,WKSP2,ITS/NDNM*ZERO,NDNM*ZERO,0/
C
C     Petsc stuff
C
      Mat A,B
      Vec RHS,DT,Z,X
      KSP solver
      PetscLogDouble TBEGALL,telapsed,tbeg,tend
      PetscScalar DT_V(1)
      PetscOffset DT_I,IDX_I
      PetscScalar x_array(1)
      PetscScalar rhs_array(1)
      PetscScalar z_array(1)
      PetscOffset i_x,i_rhs,i_z
      IS, dimension(0:*) :: NodalBcs
C     ..
C
      CALL MPI_Comm_rank(PETSC_COMM_WORLD,MY_PE,IFAIL)
C
      CALL VecNorm(RHS,NORM_2,RESL2(1,1),IFAIL)
      CALL VecNorm(RHS,NORM_INFINITY,RESMAX(1,1),IFAIL)
#ifdef DEBUG
      call VecNorm(RHS,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||RHS|| before  is ',s
      call VecNorm(Z,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||Z|| before everything else is ',s
#endif
#ifdef DEBUG
      IF(TIMEIMPL)THEN
         call MatNorm(A,NORM_FROBENIUS,s,ifail)
         if(MY_PE.EQ.0)write(6,*)'||A|| before is ',s
      ENDIF
#endif
C
      IF( ITER .EQ. 1 )THEN
          RESMAX0(1) = RESMAX(1,1)
          RESL20(1) = RESL2(1,1)
      ENDIF
C
C     For pure advection problems, the inverse of
C         the local timestep in the Dirichlet Nodes
C         is 0., so we set it to the maximum DT
C
      IF(KAN.EQ.-1)THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,DT_V,DT_I,IFAIL)
          CALL ISGetIndices(NodalBcs(SupersonicNodes),IDX_V,IDX_I,IFAIL)
          CALL ISGetSize(NodalBcs(SupersonicNodes),
     &                   Number_of_DirichletNodes,IFAIL)
          DO 10 I = 1, Number_of_DirichletNodes
C     supersonic nodes are 0-based indexed
              IPOIN = IDX_V(IDX_I+I)+ 1- ROWBGN
#ifdef DEBUG
              s = DT_V(DT_I+IPOIN)
#endif
              DT_V(DT_I+IPOIN) = DTMAX
#ifdef DEBUG
              write(6,*)'Supersonic node was ',s,
     &' now is ',dt_v(dt_i+IPOIN)
#endif
   10     CONTINUE
          CALL VecRestoreArray(DT,DT_V,DT_I,IFAIL)
          CALL ISRestoreIndices(NodalBcs(SupersonicNodes),IDX_V,IDX_I,
     &                          IFAIL)
      ENDIF
#ifdef DEBUG
      call VecNorm(DT,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||DT|| after supersonic nodes  is ',s
#endif
C
C     Handles hanging nodes
C
      CALL ISGetSize(NodalBcs(HangingNodes),NI,IFAIL)
      IF( NI .GT. 0 .AND. IGLOB .EQ. 1 )THEN
          CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
          CALL VecGetOwnerShipRange(DT,ROWBGN,ROWEND,IFAIL)
          CALL VecGetArray(DT,x_array,i_x,IFAIL)
          CALL ISGetIndices(NodalBcs(HangingNodes),IDX_V,IDX_I,IFAIL)
          DO 13 I = 1, NI, NOFVAR
C     hanging nodes are 0-based indexed
              IPOIN = (IDX_V(IDX_I+I))+ 1- ROWBGN
#ifdef DEBUG
              s = x_array(i_x+IPOIN)
#endif
              x_array(i_x+IPOIN) = DTMAX
#ifdef DEBUG
              write(6,*)'Hanging node was ',s,
     &' now is ',x_array(i_x+IPOIN)
#endif
   13     CONTINUE
          CALL VecRestoreArray(DT,x_array,i_x,IFAIL)
          CALL ISRestoreIndices(NodalBcs(HangingNodes),IDX_V,IDX_I,
     &                          IFAIL)
      ENDIF
#ifdef DEBUG
      call VecNorm(DT,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||DT|| after hanging nodes  is ',s
#endif
C
C     Global time-stepping (if required)
C
      IF( IGLOB .EQ. 0 )THEN
          IF( LTIME )THEN ! Time-accurate
             CALL VecGetArray(DT,DT_V,DT_I,IFAIL)
             DO 12 IPOIN = 1, NPOIN ! May not be ok in mpi !
                DT_V(DT_I+IPOIN) = VMDCO(IPOIN)/DELT
   12     CONTINUE
             CALL VecRestoreArray(DT,DT_V,DT_I,IFAIL)
          ELSE ! Pseudo-un-steady
             CALL VecMax(DT,PETSC_NULL_INTEGER,DTMAX,IFAIL)
             CALL VecSet(DT,DTMAX,IFAIL)
          ENDIF
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
C
      IF( TIMEIMPL )THEN
C
C     scale the CFL
C
          CNST = RESL20(1)/RESL2(1,1)
          CNST = MIN( CFLMAX(1), CFL(1)*CNST )
C
C     divide V_i/Dt by the CFL number ...
C
          CALL VecScale(DT,ONE/CNST,IFAIL)
C
C     Adds V_i/Dt to the diagonal elements of A ...
C
#ifdef DEBUG
          call VecNorm(DT,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||dt|| before MatDiagonalSet is ',s
#endif
          CALL MatDiagonalSet(A,DT,ADD_VALUES,IFAIL)
#ifdef DEBUG
          call MatNorm(A,NORM_FROBENIUS,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||A|| after MatDiagonalSet is ',s
#endif
C
      ELSE ! explicit timestepping
caldo
C
C   0 maxes the lines INvisible to the preprocessor
#if 0
          CALL VecGetArray(RHS,DT_V,DT_I,IFAIL)
             DO 16 IPOIN = 1, NPOIN ! May not be ok in mpi !
                write(6,*)ipoin,'dV/dt = ',
     &(vmdcn(ipoin)-VMDCO(IPOIN))/delt,' rhs = ',dt_v(dt_i+ipoin)
   16     CONTINUE
          CALL VecRestoreArray(RHS,DT_V,DT_I,IFAIL)
#endif
caldo
#ifdef DEBUG
      call VecNorm(DT,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||DT|| before VecPointwiseDivide is ',s
      call VecMin(DT,PETSC_NULL_INTEGER,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'min(||DT||) before is ',s
      call VecNorm(RHS,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||RHS|| before VecPointwiseDivide is ',s
#endif
          CALL VecPointwiseDivide(RHS,RHS,DT,IFAIL)
          IF( CFL(1) .NE. ONE )CALL VecScale(RHS,CFL(1),IFAIL)
          CNST = CFL(1)
      ENDIF
C
#ifdef DEBUG
      call VecNorm(RHS,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||RHS|| after VecPointwiseDivide is ',s
#endif
C
C     We do not need DT any more
C
      CALL VecDestroy(DT,IFAIL)
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
C     modify the stiffness matrix structure to account for
C     Dirichlet b.c. this is done only for the first iteration
C
      IF( TIMEIMPL )THEN
C
CCCCCCC#ifdef MPI
CCCCCCCCALL MatZeroRowsLocal(A,SupersonicNodes,ONE,IFAIL)
CCCCCCC#else
CCCCCCC#endif
C
          CALL MatZeroRowsIS(A,NodalBcs(SupersonicNodes),ONE,
     &                       PETSC_NULL_VEC,PETSC_NULL_VEC,IFAIL)
C
C     is it really needed to re-assemble the matrix???
C
          CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,IFAIL)
          CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,IFAIL)
C
C     This will remove the entries from the matrix; but when
C     you assemble the matrix for the next time-step it will not
C     insert those unneeded values (the values that would be zeroed
C     in the next call to MatZeroRows()). This is a very efficient
C     way to do multiple linear systems BUT assumes
C
C     1) the rows you want to zero remain the same for each successive
C        linear system in the series
C     2) you are not adding new nonzeros in other parts of the matrix
C        at later timesteps (because the above option will cause those
C        new values to be ignored).
C
caldo    CALL MatSetOption(A,MAT_NO_NEW_NONZERO_LOCATIONS,IFAIL)
#ifdef DEBUG
          call MatNorm(A,NORM_FROBENIUS,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||A|| after JacBCs is ',s
#endif
C
      ENDIF
C
      IF( TIMEIMPL )THEN
C
C     create a vector to store the solution
C
#ifdef MPI
          CALL VecCreateMPI(PETSC_COMM_WORLD,NPOIN,PETSC_DECIDE,X,IFAIL)
#else
          CALL VecCreateSeq(PETSC_COMM_SELF,NPOIN,X,IFAIL)
#endif
C
C     Sets the matrix (and preconditioner) associated with the linear system
C
          CALL KSPSetOperators(solver,A,A,IFAIL)
          CALL KSPSetUp(solver,IFAIL)
C
C     Solve THE linear system
C
          CALL PetscTime(tbeg,IFAIL)
          CALL KSPSolve(solver,RHS,X,IFAIL)
          CALL KSPGetIterationNumber(solver,ITS,IFAIL)
          CALL PetscTime(tend,IFAIL)
C
C     Update the nodal unknown vector by forming Z := X + Z
C
          CALL VecAXPY(Z,ONE,X,IFAIL)
C
      ELSE ! explicit timestepping
C
C     Update the nodal unknown vector by forming Z := RHS + Z
C
          IF(LALE)THEN
            CALL VecGetArray(rhs,rhs_array,i_rhs,IFAIL)
            CALL VecGetArray(z,z_array,i_z,IFAIL)
            CALL RHS4ALE(NOFVAR,NPOIN,rhs_array(i_rhs+1),
     &                   z_array(i_z+1),VMDCO,VMDCN)
            CALL VecRestoreArray(z,z_array,i_z,IFAIL)
            CALL VecRestoreArray(rhs,rhs_array,i_rhs,IFAIL)
C           we should do it better, RHSBC1 has already been called previously
            CALL RHSBC1( rhs, NodalBcs )
C   0 maxes the lines INvisible to the preprocessor
#if 0
            CALL VecGetArray(rhs,rhs_array,i_rhs,IFAIL)
            write(6,*)dnrm2(npoin,rhs_array(i_rhs+i),1)
            do i  = 1,npoin
               write(6,*)i,rhs_array(i_rhs+i)
            enddo 
            CALL VecRestoreArray(rhs,rhs_array,i_rhs,IFAIL)
#endif
          ENDIF
C
#ifdef DO_NOT_UPDATE
      write(6,*)'Not updating the mean flow eqns'
      GOTO 456
#endif
C
#ifdef DEBUG
          call VecNorm(Z,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||Z|| before VecAXPY is ',s
          call VecNorm(RHS,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||RHS|| before VecAXPY is ',s
#endif
          CALL VecAXPY(Z,ONE,RHS,IFAIL)
#ifdef DEBUG
          call VecNorm(RHS,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||RHS|| after VecAXPY is ',s
          call VecNorm(Z,NORM_2,s,ifail)
          if(MY_PE.EQ.0)write(6,*)'||Z|| after VecAXPY is ',s
#endif
C
#ifdef DO_NOT_UPDATE
  456 CONTINUE
#endif
          X = RHS
C
      ENDIF
#ifdef DEBUG
      call VecNorm(X,NORM_2,s,ifail)
      if(MY_PE.EQ.0)write(6,*)'||X|| after KSP solve  is ',s
#endif
C
C     Monitor the norms of the update ..
C
      CALL VecNorm(X,NORM_2,DELL2(1,1),IFAIL)
      CALL VecNorm(X,NORM_INFINITY,DELMAX(1,1),IFAIL)
C
      IF(MY_PE.EQ.0)THEN
C
          WRITE(IWUNIT,200)ITER,ITS,CNST
C
C     .. Print out the convergence history ..
C
          WRITE(IWUNIT,215)
          DO 20 IVAR =  1, NOFVAR
              WRITE(IWUNIT,210)IVAR,DLOG10(RESL2(IVAR,1)), DLOG10(RESMAX
     +        (IVAR,1)),INMAX(IVAR,1),(WKSP1(I,IVAR),I=1,3)
   20     CONTINUE
          WRITE(IWUNIT,225)
          DO 30 IVAR =  1, NOFVAR
              WRITE(IWUNIT,210)IVAR,DLOG10(DELL2(IVAR,1)), DLOG10(DELMAX
     +        (IVAR,1)),INDEL(IVAR,1),(WKSP2(I,IVAR),I=1,3)
   30     CONTINUE
C
C     Writing convergence history to file ...
C
          CALL PetscTime(telapsed,IFAIL)
          telapsed=telapsed-tbegall
          WRITE (IHST1,FMT=235) NITER,ITS,tend-tbeg,telapsed,RESL2(1,1),
     +    CNST
          WRITE (IHST2,FMT=235) NITER,ITS,tend-tbeg,telapsed,RESMAX
     +    (1,1),CNST
      ENDIF
C
C     clear memory allocated for the solution vector
C
      IF(TIMEIMPL)THEN
             CALL VecDestroy(X,IFAIL)
      ENDIF
C
      RETURN
C
C
  200 FORMAT(5X,70('-'),/,25X,'ITERATION # ',I4,' (',I4,') CFL = ',
     +E10.4/,5X,70('-'),/,
     +5X,70('-'),/,5X,'Var.',4X,'L2-norm',3X,'L_infty',3X,
     +'node #',3X,'(',4X,'x',7X,'y',7X,'z',4X,')',/,5X,70('-'))
  210 FORMAT(5X,I1,5X,F10.5,1X,F10.5,2X,I5,3X,
     +'(',2(F8.5,','),F8.5,')')
  215 FORMAT(5X,'Nodal Residual'/)
  225 FORMAT(5X,'Nodal Update'/)
  235 FORMAT (I5,1X,I4,4 (1X,E10.4))
C
      END
