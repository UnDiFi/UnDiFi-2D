      SUBROUTINE myTS(ZROE,VISCT,A,B,IDUMMY,IWORK)
C
C     $Id: myTS.F,v 1.22 2008/02/22 11:04:15 abonfi Exp abonfi $
C
C     My own time-stepping routine
C
      IMPLICIT NONE
C
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscts.h"
C
#include "iset.com"
C
C#define DEBUG
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'constants'
      INCLUDE 'conv.com'
      INCLUDE 'stream.com'
      INCLUDE 'backup.com'
      INCLUDE 'flags.com'
      INCLUDE 'nloc'
      INCLUDE 'io.com'
      INCLUDE 'implicit.h'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION DSTAK(1)
      COMMON/CSTAK/DSTAK
C
      INTEGER ISTAK(1),IWORK(*)
C     ..
C     .. Petsc Vectors ..
      Vec rhs,dt,ZRoe,ViscT,rhs_t
C     ..
C     .. Petsc Matrices ..
      Mat A,B ,MatA,MatB
C     ..
C     .. Petsc KSP ..
      KSP FlowSolver,TurbulenceSolver
      MatStructure str
      TS ts_dummy
      PetscLogDouble SPACE,FRAGS,MAXS,MAXSIZ
      PetscLogDouble tbeg,tbegall
#ifdef DEBUG
      PetscScalar      x_array(1)
      PetscScalar      b_array(1),dt_v(1)
      PetscOffset i_x,dt_i,idx_i
      double precision s
      integer ivar,N,J,I,locy,locv,IFLAG
      character*8 fname
#endif
C
      COMMON/MYTSCOM/MatA,MatB,DT
C     ..
C     .. Local Scalars ..
C
      INTEGER NERR,IOPT
      CHARACTER*72 ERRMSG
      PARAMETER(NERR = 5, IOPT = 1)
C
C     NPOIN   number of nodes in the mesh
C     NELEM   number of elements in the mesh
C     NFACE   number of faces in the mesh
C     NBFAC   number of boundary faces in the mesh
C
C
C     DIM     is the space dimension
C     NOFVERT = DIM+1 is the number of vertices (ONLY linear elements
C                     are currently allowed)
C
C
      INTEGER ITIM,ITIMS,ITIME
      DOUBLE PRECISION TIME
      INTEGER NDIM,IFAIL,LWORK,NIT,NGHOST,NBFAC,NELEM,NFACE,
     +NOFVAR,NTURB,NOFVERT,NPOIN,NPNOD,IDUMMY,NITEMS
#ifdef MPI
      INTEGER LGHOST
#endif
C     ..
C     .. External Subroutines ..
C
C     ..
C     .. External Functions ..
C
      INTEGER ISTKGT,ISTKST,I1MACH
      EXTERNAL ISTKGT,ISTKST,I1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DLOG10,IABS,LOG10
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
      ts_dummy = -1
      MatA = A
      MatB = B
C
      NDIM   =IWORK(1)
      NOFVERT=IWORK(2)
      NOFVAR =IWORK(3)
      NTURB  =IWORK(4)
      NPOIN  =IWORK(5)
      NGHOST =IWORK(6)
      NELEM  =IWORK(7)
      NFACE  =IWORK(8)
      NBFAC  =IWORK(9)
      NPNOD  =IWORK(10)
C
caldo IF(MY_PE .EQ. 0)OPEN(32,FILE='fort.32')
C
C **************************************************************
C     Create linear solver context (for implicit timestepping only)
C **************************************************************
C
      IF( TIMEIMPL )THEN
          CALL KSPCreate(PETSC_COMM_WORLD,FlowSolver,IFAIL)
C
C         by keeping KSPSetFromOptions() outside the iteration loop 10
C         we make sure options are read once and forall
C
          CALL KSPSetFromOptions(FlowSolver,IFAIL)
C
          IF( TURBULENT .AND. SEGREGATED )THEN
C
              CALL KSPCreate(PETSC_COMM_WORLD,TurbulenceSolver,IFAIL)
              CALL KSPSetFromOptions(TurbulenceSolver,IFAIL)
C
          ENDIF
      ENDIF
C
C     Create vectors to keep the right hand side
C
      CALL VecDuplicate(ZROE,RHS,IFAIL)
      IF( TURBULENT .AND. SEGREGATED )
     &CALL VecDuplicate(VISCT,RHS_T,IFAIL)
C
C======================================================
C.....START TIME LOOP
C======================================================
C
      TIME = 0.d0
      ITIM = 0
      ITIMS=ITIM+1
      ITIME=ITIM+ITSTEP
C
      DO 400 ITIM=ITIMS,ITIME
         TIME=TIME+DELT
C
C     shift time levels
C     copy u^(n) into u^(n-1)
C
C
         NITEMS = (NPOIN+NGHOST+NPNOD)*NOFVAR
         CALL DCOPY(NITEMS,DSTAK(LZROE+NITEMS),1,
     &                     DSTAK(LZROE+2*NITEMS),1)
         CALL DCOPY(NITEMS,DSTAK(LZROE),1,DSTAK(LZROE+NITEMS),1)
C
         WRITE(6,*) '  '
         WRITE(6,*) '  TIME = ',TIME,' TIME-STEP = ',DELT
         WRITE(6,*) '  *********************************************'

C
C **************************************************************
C     Start iterating in pseudo-time (NONlinear iterations)
C **************************************************************
C
      CALL PetscGetTime(tbegall,IFAIL)
C
      DO 10 ITER = 1,ITMAX
C
          CALL PetscGetTime(tbeg,IFAIL)
C
          NITER = NITER + 1
caldo     IF(MY_PE.EQ.0)WRITE(32,*)ITER
C
C +------------------------------------------------------------+
C     create a vector to store the local timestep
C +------------------------------------------------------------+
C
#ifdef MPI
C
C -- set 0-based indexing 
C
      LGHOST = ISTKGT(NGHOST*NOFVAR,2)
      CALL GETIDX(NGHOST,NOFVAR,ISTAK(LTZX),ISTAK(LGHOST))
C
C -- create the time-step vector
C
          CALL VecCreateGhost(PETSC_COMM_WORLD,NPOIN*NOFVAR,
     +                        PETSC_DECIDE,NGHOST*NOFVAR,ISTAK(lghost),
     +                        DT,IFAIL)
C
C -- release memory allocated for LGHOST
C
          CALL ISTKRL(1) 
C
C
#else
C
          CALL VecCreateSeq(PETSC_COMM_SELF,NPOIN*NOFVAR,DT,IFAIL)
C
#endif
C
C -- set the block size equal to the number of variables NOFVAR
C
      IF(NOFVAR.GT.1) CALL VecSetBlockSize(DT,NOFVAR,IFAIL)
C
caldo
#ifdef DEBUG
      CALL VecNorm(ZROE,NORM_2,s,IFAIL)
      write(6,*)' ||Z|| upon entry = ',s
#endif
caldo
C
          CALL RHSFunction(ts_dummy,tbeg,ZROE,RHS,IWORK)
C
          IF( TIMEIMPL )
     1    CALL RHSJacobian(ts_dummy,tbeg,ZROE,A,A,str,IWORK)
C
#ifdef DEBUG
          CALL PetscOptionsHasName(PETSC_NULL_CHARACTER,'-test_lhs',
     +    IFLAG,IFAIL)
          IF(IFLAG.NE.0)THEN
              locy = istkgt(npoin*nofvar,4)
              locv = istkgt(npoin*nofvar,4)
              call testlhs(A,ZROE,RHS,DSTAK(LOCY),DSTAK(LOCV), ISTAK
     +        (LNODCOD),NPOIN,NOFVAR)
              stop
          ENDIF
#endif
C
C         Update the solution
C
          IF(ABS(KAN).EQ.1)THEN
              CALL UPDATE2(NDIM,NOFVAR,NPOIN,FlowSolver,A,RHS,DT, ZROE)
#ifdef DEBUG
      CALL VecNorm(ZROE,NORM_2,s,IFAIL)
      write(6,*)' ||Z|| before leaving = ',s
#endif
C
          ELSE
              IF (IABS(KAN).EQ.4) CALL PARM_TO_CONS(DSTAK(LZROE), NDIM,
     +        NOFVAR,NPOIN+NGHOST,LTSQR,IFAIL)
Caldo
CC#define DUMP
Caldo
#ifdef DUMP
      CALL VecGetArray(rhs,x_array,i_x,IFAIL)
      call dcopy(nofvar*npoin,x_array(i_x+1),1,dstak(lzroe),1)
      CALL VecRestoreArray(rhs,x_array,i_x,IFAIL)
      LWORK = ISTKGT(MAX(NOFVAR,NTURB)*NPNOD,2)
      CALL BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB,ISTAK(LPMAP),
     +            ISTAK(LWORK),"res010.dat","res015.dat")
      call ISTKRL(1)
      stop
#endif
Caldo
Caldo
C
              CALL UPDATE3(NDIM,NOFVAR,NPOIN,FlowSolver,A,RHS,DT, ZROE)
 
              IF (IABS(KAN).EQ.4) THEN
                 CALL CONS_TO_PARM(DSTAK(LZROE), NDIM,
     +        NOFVAR,NPOIN+NGHOST,LTSQR,IFAIL)
                IF(IFAIL.NE.0)THEN
                    WRITE(ERRMSG(1:64),FMT=105)'myTS______',MY_PE,"(0)"
                    CALL SETERR(ERRMSG(1:64),64,NERR,IOPT)
                ENDIF   ! ifail
              ENDIF ! KAN
          ENDIF
C
C    Turbulent equations
C
          IF( TURBULENT .AND. SEGREGATED )THEN
C
#ifdef MPI
C
C     update the ghost regions of the Z vector
C     with correct values from the owning process
C
      CALL VecGhostUpdateBegin(ZROE,INSERT_VALUES,
     &                         SCATTER_FORWARD,IFAIL)
      CALL VecGhostUpdateEnd(ZROE,INSERT_VALUES,
     &                       SCATTER_FORWARD,IFAIL)
C
#endif
C
              DO 7 NIT = 1,NSUBIT
C
C     create a vector to store the local timestep
C     note that the timestep vector (DT) had been destroyed
C     in the UPDATE_ routines
C
#ifdef MPI
C
              LGHOST = ISTKGT(NGHOST,2)
              CALL GetIdx(NGHOST,1,ISTAK(Ltzx),ISTAK(LGHOST))
              CALL VecCreateGhost(PETSC_COMM_WORLD,NPOIN,PETSC_DECIDE,
     +                            NGHOST,ISTAK(lghost),DT,IFAIL)
              CALL ISTKRL(1) 
C
#else
C
              CALL VecCreateSeq(PETSC_COMM_SELF,NPOIN,DT,IFAIL)
#endif
C
              CALL RHSFunction_t(ts_dummy,tbeg,VISCT,RHS_T,IWORK)
C
              CALL RHSJacobian_t(ts_dummy,tbeg,VISCT,B,B,str,IWORK)
c
              CALL UPDATE4(NDIM,NTURB,NPOIN,TurbulenceSolver,B, RHS_T,
     +                     ViscT,DT,NIT,IFAIL)
C
              IF(IFAIL.EQ.0)GOTO 9
C
    7         CONTINUE
    9         CONTINUE
#ifdef MPI
C
C     update the ghost regions of the VISCT vector
C     with correct values from the owning process
C
      CALL VecGhostUpdateBegin(VISCT,INSERT_VALUES,
     &                         SCATTER_FORWARD,IFAIL)
      CALL VecGhostUpdateEnd(VISCT,INSERT_VALUES,
     &                       SCATTER_FORWARD,IFAIL)
C
#endif
C
          ENDIF
C
C  =========================================================
C  Check NONlinear convergence
C  =========================================================
C
!         IF (LOG10(RESL2(IVCNVG,1)/RESL20(1)).LE.TOLER) GOTO 12
          IF (LOG10(RESL2(IVCNVG,1)).LE.TOLER) GOTO 12
C
C
C  =========================================================
C  Backing up (if required)
C  =========================================================
C
#if 0
      LWORK = ISTKGT(MAX(NOFVAR,NTURB)*NPNOD,2)
      fname(1:1) = "t"
      fname(5:8) = ".dat"
      write(fname(2:4),FMT="(I3.3)")iter
      CALL BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB,ISTAK(LPMAP),
     +            ISTAK(LWORK),fname,VISCTFILE)
      CALL ISTKRL(1)
#endif
C
C
          IF ((ITER/IBAK)*IBAK.EQ.ITER) THEN
C
caldo
caldo WRITE(BAKFILE(4:7),FMT="(I4.4)")iter+1000
caldo
      LWORK = ISTKGT(MAX(NOFVAR,NTURB)*NPNOD,2)
      CALL BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB,ISTAK(LPMAP),
     +            ISTAK(LWORK),BAKFILE,VISCTFILE)
      IF(IABS(KAN).NE.1)
     +CALL CHKFLX(ISTAK(LBNDFAC),NBFAC,ISTAK(LCELFAC),ISTAK(LCELNOD),
     +NOFVERT,NELEM,DSTAK(LFACNOR),DSTAK(LZROE),NOFVAR,NDIM,NFACE,
     +(IABS(KAN).EQ.4),ITER)
      CALL ISTKRL(1)
C
C
C     Write skin friction data to a file
C
#ifndef MPI
#if 0
Crem
Crem  core dumps
Crem
              LWORK = ISTKGT(2*NBFAC,2)
              CALL WSKIN(DSTAK(LSKINF),DSTAK(LHEAT),ISTAK(LWORK),NWFAC,
     +        ISTAK(LBNDFAC),NBODY6,'file016.dat')
              CALL ISTKRL(1)
#endif
#endif
C
C
C **************************************************************
C     Checking memory usage
C **************************************************************
C
              CALL PetscMemoryGetCurrentUsage(MAXS,IFAIL)
!             CALL PetscMemoryGetMaximumUsage(MAXSIZ,IFAIL)
              WRITE (NOUT,FMT=9998) ISTKST(3)/2, (ISTKST(3)*100.)/ISTKST
     +        (4)
              WRITE (NOUT,FMT=9999) MAXS
!             WRITE (NOUT,FMT=9994) MAXSIZ
C
          ENDIF
C
   10 CONTINUE ! end loop over ITER
   12 CONTINUE
  400 CONTINUE ! end time loop
C
caldo IF(MY_PE.EQ.0)CLOSE(32)
C
      CALL VecDestroy(RHS,IFAIL)
      IF(TIMEIMPL) CALL KSPDestroy(FlowSolver,IFAIL)
      IF(TURBULENT .AND. SEGREGATED )THEN
          CALL VecDestroy(RHS_T,IFAIL)
          CALL KSPDestroy(TurbulenceSolver,IFAIL)
      ENDIF
C
      RETURN
  105 FORMAT('CONS_TO_PARM failed when called from ',A10,' on PE # ',I4,
     &1X,A3)
 9994 FORMAT (/,5X,'MAXIMUM MEMORY USED BY THE PROGRAM ',F12.0, ' BYTES'
     +)
 9995 FORMAT (/,5X,'MEMORY ALLOCATED IN DSTAK ',I8,' REAL*8 WORDS')
 9996 FORMAT (5X,'Skin friction file written to file016.dat')
 9997 FORMAT (I6,5(1X,E12.6))
 9998 FORMAT (/,5X,'MAXIMUM MEMORY USED IN DSTAK',I8, ' REAL*8 WORDS (',
     +F5.2,' %)')
 9999 FORMAT (/,5X,'MAXIMUM MEMORY EVER ALLOCATED BY PETSc ',F12.0,
     +' BYTES')
 
      END
