      subroutine RHSFunction_b(ts,t,global_in,global_out,dummy)
C
C     This subroutine is "registered" in TSSetRHSFunction
C     the calling sequence CANNOT be changed
C
C     global_in  is the vector value where the function is evaluated
C     global_out is the function value
C     global_out = RHSFunction(global_in)
C
      implicit none
C
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscts.h"
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'datatype.h'
      INCLUDE 'conv.com'
      INCLUDE 'nloc.com'
      INCLUDE 'flags.com'
C
      PetscLogDouble tbeg,tend
      INTEGER MY_PE
C
C  Input/output parameters:
      TS               ts
      double precision t
      Vec              global_in
      Vec              global_out
      integer          dummy(*)
C
C  Local variables:
      Mat              A,B
      Vec              local_in
      Vec              local_out
      Vec              local_dt
      integer          ifail,NDIM,NDIM,NPOIN,NGHOST,
     &                 LINDX,IADDR,NPNOD
#ifdef MPI
      integer          bs
#endif
      PetscOffset i_x
      PetscScalar x_array(1)
C
      Vec              global_dt
      COMMON/MYTSCOM/A,B,global_dt
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))

      INTEGER  ISTKGT
      EXTERNAL ISTKGT
C
      NDIM      =dummy( 1)
      NDIM    = dummy( 4)
      NPOIN    = dummy( 5)
      NGHOST   = dummy( 6)
      NPNOD    = dummy(10)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     Get ready for local function computations
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
#ifdef MPI
C
C     update the ghost regions with correct values from the owning process
C
      CALL VecGhostUpdateBegin(global_in,INSERT_VALUES,
     &                         SCATTER_FORWARD,IFAIL)
      CALL VecGhostUpdateEnd(global_in,INSERT_VALUES,
     &                       SCATTER_FORWARD,IFAIL)
C
      CALL VecGhostGetLocalForm(global_out,local_out,IFAIL)
!     CALL VecGetBlockSize(local_out,bs,IFAIL)
!     IF(BS.NE.NDIM.AND.BS.GT.1)
!    +CALL VecSetBlockSize(local_out,NDIM,IFAIL)
C
C     get local representation of the r.h.s. vector
C
      CALL VecGhostGetLocalForm(global_in,local_in,IFAIL)
!     CALL VecGetBlockSize(local_in,bs,IFAIL)
!     IF(BS.NE.NDIM.AND.BS.GT.1)
!    +CALL VecSetBlockSize(local_in,NDIM,IFAIL)
C
C     get local representation of the time-step vector
C
!     CALL VecGhostGetLocalForm(global_dt,local_dt,IFAIL)
C
#else
C
C     for the uniprocessor case local and global numbering 
C     are the same, maybe the following assignement
C     is already done in VecGhostGetLocalForm(), controlla!
C
      local_out=global_out
      local_in=global_in
!     local_dt=global_dt
C
#endif
C
C     Update periodic nodes
C
      LINDX = ISTKGT(NDIM*NPNOD,KIND_INTEGER)
      IADDR = LXYZDOT+NDIM*(NPOIN+NGHOST)
      CALL UPDTP(DSTAK(LXYZDOT),DSTAK(IADDR),
     &NDIM,NPOIN,NGHOST,ISTAK(LPMAP),ISTAK(LINDX),NPNOD,.FALSE.)
      CALL ISTKRL(1)
C
C     **********************************************************
C
C     BUILD_RHS = .TRUE.
C     BUILD_JAC = .NOT. BUILD_RHS
C
C
!     ICOUNT=ICOUNT+1
C
      CALL PetscTime(tbeg,IFAIL)
      CALL VecGetArray(local_in,x_array,i_x,IFAIL) 
      CALL ApplicationFunction_b(B,x_array(i_x+1),local_out,
     &                           dummy,IFAIL)
      CALL VecRestoreArray(local_in,x_array,i_x,IFAIL) 
      CALL PetscTime(tend,IFAIL)
C
C     **********************************************************
C
C     Must assembly RHS now since a call to MatSetValues(INSERT_VALUES)
C     might follow in RHSBC5

      CALL VecAssemblyBegin(global_out,IFAIL)
      CALL VecAssemblyEnd(global_out,IFAIL)
C
C update ghost values
C
#ifdef MPI
C
      CALL VecGhostUpdateBegin(global_out,ADD_VALUES,
     &                         SCATTER_REVERSE,IFAIL) 
      CALL VecGhostUpdateEnd(global_out,ADD_VALUES,
     &                         SCATTER_REVERSE,IFAIL) 
#endif
C
C
C     Apply b.c. to the nodal residual
C     this must be done AFTER the contributions from the other procs.
C     have been assembled
C
      CALL RHSBC5(global_out)
C
!     CALL VecAssemblyBegin(global_dt,IFAIL)      
!     CALL VecAssemblyEnd(global_dt,IFAIL)      
C
#ifdef MPI
C
      CALL VecGhostRestoreLocalForm(global_out,local_out,IFAIL)
!     CALL VecGhostRestoreLocalForm(global_dt,local_dt,IFAIL)
C
!     CALL VecGhostUpdateBegin(global_dt,ADD_VALUES,SCATTER_REVERSE,
!    &IFAIL) 
!     CALL VecGhostUpdateEnd(global_dt,ADD_VALUES,SCATTER_REVERSE,
!    &IFAIL) 
C
#endif
C
C
  100 FORMAT(4(F12.5,1X))
  235 FORMAT (E10.4,1X,1(I4,1X),6(E10.4,1X))
C
      RETURN
      END
