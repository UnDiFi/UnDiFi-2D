      subroutine RHSFunction(ts,t,global_dt,global_in,global_out,
     &                       A,dummy,NodalBcs)
C
CCCC#define DEBUG
C
C     $Id: RHSFunction.F,v 1.26 2020/04/24 06:01:52 abonfi Exp $
C
C     This subroutine is "registered" in TSSetRHSFunction
C     the calling sequence CANNOT be changed
C     NodalBcs has been ADDED to cope with petsc-3.8.*
C
C     global_in  is the vector value where the function is evaluated
C     global_out is the function value, i.e.
C     global_out = RHSFunction(global_in)
C
C     by global/local we refer to the global and local representations
C     of the same vector, i.e. whether we address it using local
C     or global nodenumbering
C
C     one of the key ideas about the present routine is
C     to confine all message-passing stuff here
C
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscts.h"
      use petscvec
      use petscmat
      use petscts
C
      implicit none
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'datatype.h'
      INCLUDE 'conv.com'
      INCLUDE 'periodic.com'
      INCLUDE 'nloc.com'
      INCLUDE 'flags.com'
      INCLUDE 'bnd.h'
      INCLUDE 'iset.h'
C
      PetscLogDouble tbeg,tend
C
C  Input/output parameters:
      TS               ts
      IS, dimension(0:*) :: NodalBcs
      double precision t
      Vec              global_in
      Vec              global_out
      integer          dummy(*)
C
C  Local variables:
      Vec              local_in
      Vec              local_out
      Vec              local_dt
cxxx  Mat              A,B,C
      integer          ifail,NDIM,NOFVAR,NPOIN,NGHOST,
     &ibgn,iend,ivar,IADDR,LINDX,NPNOD
      integer i,ni
#ifdef MPI
      integer          bs
#endif
      double precision s 
      PetscOffset i_rhs,i_x
      PetscScalar rhs_array(1),x_array(1)
#ifdef DEBUG
#include "finclude/petscviewer.h"
      PetscBool flg
      PetscViewer MyOpenMindedViewer
      PetscScalar xx_v(1),yy_v(1)
      PetscErrorCode ierr
      PetscOffset xx_i,yy_i
      integer i,j,k,iunit 
      CHARACTER* 6 matfile,rhsfile,solfile
      double precision dnrm2
      DATA matfile,rhsfile,solfile/"bcsXXX","rhsXXX","solXXX"/
#endif
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER MY_PE
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON/MPICOM/MY_PE
      Mat              A
      Vec              global_dt
cxxx  COMMON/MYTSCOM/A,B,C,global_dt
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
      INTEGER  ISTKGT
      EXTERNAL ISTKGT
C
      NDIM      =dummy( 1)
      NOFVAR    =dummy( 3)
      NPOIN     =dummy( 5)
      NGHOST    =dummy( 6)
      NPNOD     =dummy(10)
C
!     write(6,*)(dummy(i),i=1,10)
C
!     i = NoSlipVelocities
!     CALL ISGetSize(NodalBcs(i),NI,IFAIL) 
!     write(6,*)'RHSF:Index set ',I,' has ',NI,' entries'
!     i = Isothermal
!     CALL ISGetSize(NodalBcs(i),NI,IFAIL) 
!     write(6,*)'RHSF:Index set ',I,' has ',NI,' entries'
!     i = SupersonicVariables
!     CALL ISGetSize(NodalBcs(i),NI,IFAIL) 
!     write(6,*)'RHSF:Index set ',I,' has ',NI,' entries'
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Get ready for local function computations
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
#ifdef MPI
C
C     update the ghost regions of the input vector
C     with correct values from the owning process
C
      CALL VecGhostUpdateBegin(global_in,INSERT_VALUES,
     &                         SCATTER_FORWARD,IFAIL)
      CALL VecGhostUpdateEnd(global_in,INSERT_VALUES,
     &                       SCATTER_FORWARD,IFAIL)
C
C     get local representation of the input vector
C
      CALL VecGhostGetLocalForm(global_in,local_in,IFAIL)
      CALL VecGetBlockSize(local_in,bs,IFAIL)
      IF(BS.NE.NOFVAR.AND.NOFVAR.GT.1)THEN
          WRITE(6,*)'bs of local_in is ',bs,' rather than ',NOFVAR
      ENDIF
C
C     get the local representation of the r.h.s. (output) vector
C
      CALL VecGhostGetLocalForm(global_out,local_out,IFAIL)
      CALL VecGetBlockSize(local_out,bs,IFAIL)
      IF(BS.NE.NOFVAR.AND.NOFVAR.GT.1)THEN
          WRITE(6,*)'bs of local_out is ',bs,' rather than ',NOFVAR
      ENDIF
C
C     get local representation of the time-step vector
C
      CALL VecGhostGetLocalForm(global_dt,local_dt,IFAIL)
      CALL VecGetBlockSize(local_dt,bs,IFAIL)
      IF(BS.NE.NOFVAR.AND.NOFVAR.GT.1)THEN
          WRITE(6,*)'bs of local_dt is ',bs,' rather than ',NOFVAR
      ENDIF
C
#else
C
C     for the uniprocessor case local and global numbering 
C     are the same, maybe the following assignement
C     is already done in VecGhostGetLocalForm(), controlla!
C
      local_out=global_out
      local_in=global_in
      local_dt=global_dt
C
#endif
#ifdef DEBUG
!     CALL VecNorm(global_in,NORM_2,s,IFAIL)
         do ivar = 1,nofvar
            CALL VecStrideNorm(global_in,ivar-1,NORM_2,s,IFAIL)
            if(my_pe.EQ.0)write(6,*)'Z norm inside RHSFunction =',ivar,s
         enddo
#endif
C
C     Update periodic nodes
C
      LINDX = ISTKGT(NOFVAR*NPNOD,KIND_INTEGER)
      IADDR = LZROE+NOFVAR*(NPOIN+NGHOST)
      CALL UPDTP(DSTAK(LZROE),DSTAK(IADDR),
     &NOFVAR,NPOIN,NGHOST,ISTAK(LPMAP),ISTAK(LINDX),NPNOD,ANNULAR)
      CALL ISTKRL(1)
C
C     **********************************************************
C
C     BUILD_RHS = .TRUE.
C     BUILD_JAC = .NOT. BUILD_RHS
C
C     ICOUNT=ICOUNT+1
C
      CALL PetscTime(tbeg,IFAIL)
      CALL VecGetArray(local_in,x_array,i_x,IFAIL) 
#ifdef DEBUG
      call VecValid(local_out,flg,iFAIL)
       if (flg .ne. PETSC_TRUE) then
          if (MY_PE .eq. 0) then
             write(6,*) 'Corrupted vector local_out!'
          endif
          STOP
          SETERRQ(1,' ',IFAIL)
       endif
#endif
      CALL ApplicationFunction(A,x_array(i_x+1),
     &                         local_out,local_dt,dummy,IFAIL)
      CALL VecRestoreArray(local_in,x_array,i_x,IFAIL) 
      CALL PetscTime(tend,IFAIL)
C
C     **********************************************************
C
C     Must assembly RHS now since we shall call
C     VecGetArray() in RHSBC4 
C
      CALL VecAssemblyBegin(global_out,IFAIL)
      CALL VecAssemblyEnd(global_out,IFAIL)
#ifdef DEBUG
!     do ivar = 1,nofvar
!     CALL VecStrideNorm(global_in,ivar-1,NORM_2,s,IFAIL)
!     if(my_pe.EQ.0)write(32,*)'zroe norm inside RHSF before GH',ivar,s
!     enddo
      do ivar = 1,nofvar
      CALL VecStrideNorm(global_out,ivar-1,NORM_2,s,IFAIL)
      if(my_pe.EQ.0)write(6,*)'rhs norm inside RHSF before GH',ivar,s
      enddo
#endif
#ifdef DEBUG
      CALL VecNorm(global_out,NORM_2,s,IFAIL)
      IF(MY_PE.EQ.0)THEN
         WRITE(6,*)'RHS norm in RHSFunction before GhostUpd is ',s
      ENDIF
C
C
C
#endif
C
C update ghost values in the r.h.s (output) vector
C
#ifdef MPI
C
      CALL VecGhostUpdateBegin(global_out,ADD_VALUES,
     &                         SCATTER_REVERSE,IFAIL) 
      CALL VecGhostUpdateEnd(global_out,ADD_VALUES,
     &                         SCATTER_REVERSE,IFAIL) 
caldo
#endif
#ifdef DEBUG
!     do ivar = 1,nofvar
!     CALL VecStrideNorm(global_in,ivar-1,NORM_2,s,IFAIL)
!     if(my_pe.EQ.0)write(32,*)'zroe norm inside RHSF after GH',ivar,s
!     enddo
      do ivar = 1,nofvar
      CALL VecStrideNorm(global_out,ivar-1,NORM_2,s,IFAIL)
      if(my_pe.EQ.0)write(6,*)'rhs norm inside RHSF after GH',ivar,s
      enddo
#endif
#ifdef DEBUG
      CALL VecNorm(global_out,NORM_2,s,IFAIL)
      IF(MY_PE.EQ.0)THEN
         WRITE(6,*)'RHS norm in RHSFunction AFTER GhostUpdate is ',s
      ENDIF
C
#endif
C
C     Do the assembly of the time-step vector ....
C
      CALL VecAssemblyBegin(global_dt,IFAIL)      
      CALL VecAssemblyEnd(global_dt,IFAIL)      
C
#ifdef MPI
C
C     .... Update ghost values ....
C
      CALL VecGhostUpdateBegin(global_dt,ADD_VALUES,SCATTER_REVERSE,
     &                         IFAIL) 
      CALL VecGhostUpdateEnd(global_dt,ADD_VALUES,SCATTER_REVERSE,
     &                       IFAIL) 
cnew
caldo CALL VecGhostUpdateBegin(global_dt,INSERT_VALUES,
caldo&                         SCATTER_FORWARD,IFAIL) 
caldo CALL VecGhostUpdateEnd(global_dt,INSERT_VALUES,
caldo&                         SCATTER_FORWARD,IFAIL) 
cnew
C
C     ... and restore global representation
C
      CALL VecGhostRestoreLocalForm(global_dt,local_dt,IFAIL)
C
#endif
C
C     do 44 ivar = 1,nofvar
C     CALL VecStrideNorm(global_dt,ivar-1,NORM_2,s,IFAIL)
C     if(my_pe.EQ.0)write(6,*)'dt norm inside RHSF',ivar,s
C  44 continue
C
C
C     Apply b.c. to the nodal residual
C     N.B. in the parallel case one can only access the local
C     portion of an array
C     nodal boundary conditions address nodes by global numbering,
C     thus we have to recover local numbering by substracting
C     the starting address (ibgn) of the vector on the current processor
C
      IF(NOFVAR.EQ.1)THEN
         CALL RHSBC1(global_out, NodalBcs)
      ELSE
         CALL VecGetOwnershipRange(global_in,ibgn,iend,IFAIL)
         CALL VecGetArray(global_in,x_array,i_x,IFAIL) 
         CALL VecGetArray(global_out,rhs_array,i_rhs,IFAIL) 
C
C     boundary conditions for the mean flow equations 
C
         CALL RHSBC4(x_array(i_x+1),rhs_array(i_rhs+1),NodalBcs,ibgn,
     &               NDIM,DSTAK(LFREE),(ABS(KAN).EQ.4))
C
C     boundary conditions for the turbulent flow quantities:
C     this are kept separate so that if a different turb. 
C     model is used a different routine will be called
C
         IF(TURBULENT.AND.COUPLED)
     &CALL RHSBC5c(rhs_array(i_rhs+1),ibgn,NodalBcs)
C
C
         CALL VecRestoreArray(global_in,x_array,i_x,IFAIL) 
         CALL VecRestoreArray(global_out,rhs_array,i_rhs,IFAIL) 
#ifdef MPI
caldo
caldo this is meaningful only if you want to have correct ghosted
caldo values in the r.h.s. for debugging purposes
caldo
!     CALL VecGhostUpdateBegin(global_out,INSERT_VALUES,
!    &                         SCATTER_FORWARD,IFAIL) 
!     CALL VecGhostUpdateEnd(global_out,INSERT_VALUES,
!    &                         SCATTER_FORWARD,IFAIL) 
caldo
#endif
      ENDIF
C
#ifdef DEBUG
      CALL VecNorm(global_out,NORM_2,s,IFAIL)
      do ivar = 1,nofvar
         CALL VecStrideNorm(global_out,ivar-1,NORM_2,s,IFAIL)
      IF(MY_PE.EQ.0)THEN
         WRITE(6,*)'RHS norm in RHSFunction after BCs is ',ivar,s
      ENDIF
      enddo
#endif
C
C
caldo CALL MPI_Barrier(PETSC_COMM_WORLD,IFAIL) 
caldo CALL MPI_Abort(PETSC_COMM_WORLD,-10,IFAIL) 
C
C
C     Before leaving this routine we restore the global
C     representation of the vectors
C
#ifdef MPI
C
      CALL VecGhostRestoreLocalForm(global_out,local_out,IFAIL)
      CALL VecGhostRestoreLocalForm(global_in,local_in,IFAIL)
C
#endif
C
C     CALL MPI_Barrier(PETSC_COMM_WORLD,IFAIL) 
C     CALL MPI_Abort(PETSC_COMM_WORLD,-10,IFAIL) 
#ifdef DEBUG
!     do ivar = 1,nofvar
!     CALL VecStrideNorm(global_dt,ivar-1,NORM_2,s,IFAIL)
!     if(my_pe.EQ.0)write(6,*)'dt norm inside RHSF',ivar,s
!     enddo
!     do ivar = 1,nofvar
!     CALL VecStrideNorm(global_in,ivar-1,NORM_2,s,IFAIL)
!     if(my_pe.EQ.0)write(32,*)'zroe norm inside RHSF after BCs',ivar,s
!     enddo
      do ivar = 1,nofvar
      CALL VecStrideNorm(global_out,ivar-1,NORM_2,s,IFAIL)
      if(my_pe.EQ.0)write(6,*)'rhs norm inside RHSF after BCs',ivar,s
      enddo
#endif
C
      RETURN
      END
