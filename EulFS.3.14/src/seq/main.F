!> \par Purpose
!>                                                                   
!> \verbatim
!>       An unstructured triangular/tetrahedral code for             
!>         solving the compressible and incompressible               
!>           Euler and (RANS) NavierStokes equations                
!>           (as well as scalar advectiondiffusion)                 
!>            using multidimensional upwind residual                 
!>         distribution (Fluctuation Splitting) schemes.             
!> \endverbatim
!>                                                                   
CC   $Id: main.F,v 1.98 2020/04/24 06:02:45 abonfi Exp $
!>                                                                   
!>                                                                   
!>       Written by:                                                 
!>       Aldo Bonfiglioli                                            
!>       Dipartimento di Ingegneria e Fisica dell' Ambiente          
!>       Universita` della Basilicata                                
!>       Contrada Macchia Romana                                     
!>       85100 Potenza Italy                                         
!>       Tel: ++39.0971.205.203                                      
!>       Email:aldo.bonfiglioli@unibas.it                           
!>       Home page: http://www.unibas.it/utenti/bonfiglioli/www.html 
!>                                                                   
!> \author $Author: abonfi $
!> \version $Revision: 1.98 $
!> \date $Date: 2020/04/24 06:02:45 $
!> 
!> 
      PROGRAM EULFS
C
#define XPRINT
C
CCCC#include "petsc/finclude/petscsys.h"
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscmat.h"
CCCC#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscts.h"
      use petscts
C
      IMPLICIT NONE
C
#include "mapping.com"
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'datatype.h'
      INCLUDE 'implicit.h'
      INCLUDE 'time.h'
      INCLUDE 'constants.h'
C
      INCLUDE 'stack.com'
      INCLUDE 'conv.com'
      INCLUDE 'verbose.com'
      INCLUDE 'stream.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'backup.com'
      INCLUDE 'flags.com'
      INCLUDE 'nloc.com'
      INCLUDE 'io.com'
      INCLUDE 'periodic.com'
      INCLUDE 'newtoneps.com'
      INCLUDE 'time.com'
C
      INTEGER MY_PE,NPROCS
      COMMON/MPICOM/MY_PE,NPROCS
C
C
      INTEGER ISTAK(1)
C     ..
C     .. Petsc Vectors ..
      Vec ZRoe,ViscT
C     ..
C     .. Petsc Matrices ..
      Mat, dimension (3) :: PetscMatrices
C     ..
C     .. Petsc array of Isets ..
      IS, dimension (0:NCOLOR+9) :: NodalBcs
C     ..
C     .. Petsc Sles ..
      PetscLogDouble tbeg,tend,totime
C     ..
C     .. Local Scalars ..
C
C
C GRID RELATED VARIABLES
C
C
C     NPOIN   number of nodes in the mesh
C     NBPOIN  number of boundary points in the mesh
C     NELEM   number of elements in the mesh
C     NFACE   number of faces in the mesh
C     NBFAC   number of boundary faces in the mesh
C     NWFAC   number of viscous wall boundary faces in the mesh
C     NHOLE   number of bodies(holes) in the mesh (ONLY for 2D)
C
C
C     NDIM     is the space dimension
C     NOFVERT = NDIM+1 is the number of vertices (ONLY linear elements
C                     are currently allowed)
C
C
      INTEGER NDIM,I,J,IFAIL,LWORKD,NGHOST,LWORKA,LWORKB,
     +        NBFAC,NBODY4,NBPOIN,NELEM,NFACE,NHOLE,NOPT,
     +        NOFVAR,NTURB,NOFVERT,NPOIN,NPNOD,NITEMS,
     +        rowbgn,rowend,NBINT,IWORK(15),OPT(10),NPROBES
      integer ilev,k
      double precision s,dnrm2
      DOUBLE PRECISION HELP
      LOGICAL LFLAG
#ifdef MPI
      INTEGER LGHOST
      CHARACTER OUTFILE*10
#endif
      CHARACTER DATADIR*255,MESHFILE*255,NGHBFILE*255,PTOTFILE*255,
     +          PBCSFILE*255,STARFILE*255,WDISTFILE*255,PROBEOUT*16,
     +          EXT*3,TIMEFILE*10,CLFILE*255,SKINFILE*255,PROBEDATA*17,
     4          SOMEFILE*255
C     ..
C     .. External Subroutines ..
C
C     ..
C     .. External Functions ..
C
      INTEGER ISTKGT,ISTKST,I1MACH
      DOUBLE PRECISION D1MACH,STORE
      EXTERNAL ISTKGT,ISTKST,I1MACH,D1MACH,STORE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DLOG10,IABS,LOG10
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
      PetscLogDouble tbegall
      COMMON/TIMING/tbegall
C     ..
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                 Beginning of program
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C
      HELP = STORE(DSQRT(D1MACH(4)))
C
C **************************************************************
C     Set the length of the static stack ..
C **************************************************************
C
      CALL ISTKIN(NVA,4)
C
C **************************************************************
C     Initialize PETSc
C **************************************************************
C
      CALL PetscInitialize(PETSC_NULL_CHARACTER,IFAIL)
      CALL PetscTime(tbeg,IFAIL)
      totime = 0.d0
      tbegall = tbeg
!     CALL PetscMemorySetGetMaximumUsage(IFAIL)
C
C **************************************************************
C     Get the number (MY_PE) of the current processor; will be
C         0 for the uni-processor case
C **************************************************************
C
      CALL MPI_Comm_rank(PETSC_COMM_WORLD,MY_PE,IFAIL)
      CALL MPI_Comm_size(PETSC_COMM_WORLD,NPROCS,IFAIL)
C
C **************************************************************
C     Set OUTPUT devices
C **************************************************************
C
#ifdef MPI
      IWUNIT = 8
      NOUT = 10
      WRITE(EXT,FMT="(I3.3)")MY_PE+1
      OUTFILE = 'output.' // EXT
      OPEN(UNIT=NOUT,FILE=OUTFILE,FORM="FORMATTED",STATUS="UNKNOWN")
      IF(MY_PE.EQ.0)THEN
          OPEN(UNIT=IWUNIT,FILE='fspl.out',FORM="FORMATTED",
     +    STATUS="UNKNOWN")
      ENDIF
#else
      NOUT = I1MACH(2)
      IWUNIT = NOUT
#endif
C
      DO 778 I = 1,255
         MESHFILE(I:I)  = CHAR(0)
         NGHBFILE(I:I)  = CHAR(0)
         STARFILE(I:I)  = CHAR(0)
         FOLDFILE(I:I)  = CHAR(0)
         BAKFILE(I:I)   = CHAR(0)
         WDISTFILE(I:I) = CHAR(0)
         VISCTFILE(I:I) = CHAR(0)
         PBCSFILE(I:I) = CHAR(0)
         PTOTFILE(I:I) = CHAR(0)
         SKINFILE(I:I) = CHAR(0)
         CLFILE(I:I) = CHAR(0)
  778 CONTINUE
C
C
C **************************************************************
C     file for timing (written by all processors)
C **************************************************************
C
      WRITE(EXT,FMT="(I3.3)")MY_PE
      TIMEFILE = 'timing.' // EXT
      OPEN(UNIT=ITIM1,FILE=TIMEFILE,FORM="FORMATTED",STATUS="UNKNOWN")
C
C **************************************************************
C     Read info from the input datafile
C **************************************************************
C
      CALL READVAL(MESHFILE,NGHBFILE,STARFILE,FOLDFILE,BAKFILE,
     +             VISCTFILE,WDISTFILE,PBCSFILE,PTOTFILE,CLFILE,
     2             SKINFILE,DATADIR,NTURB)
C
C **************************************************************
C     Open files to dump the convergence history
C **************************************************************
C
      IF(MY_PE.EQ.0)THEN
      OPEN (IHST1,FILE='convhst.l2',FORM='formatted',STATUS='UNKNOWN')
      OPEN (IHST2,FILE='convhst.max',FORM='formatted',STATUS='UNKNOWN')
          IF(TURBULENT .AND. SEGREGATED)
     +OPEN(IHST3,FILE='convhst.turb',FORM='formatted',STATUS='UNKNOWN')
      ENDIF 
C
      IF( IABS(KAN) .EQ. 3 .AND. LAPLACE )THEN
      IWEFLX = 117
      OPEN(UNIT=IWEFLX,FILE='eflux.out',FORM="FORMATTED",
     +    STATUS="UNKNOWN")
      ENDIF
C
C **************************************************************
C     Reading datafiles
C **************************************************************
C
      CALL READAT(NELEM,NPOIN,NGHOST,NPNOD,NBPOIN,NFACE,NBFAC,NBINT,
     +            NHOLE,NOFVERT,NDIM,NOFVAR,NTURB,
     2            MESHFILE,NGHBFILE,PBCSFILE,CLFILE)
caldo
      IF(ABS(KAN).EQ.2.OR.ABS(KAN).EQ.4.OR.ABS(KAN).EQ.3)
     &CALL SETIBC(NDIM,NOFVAR)
      IF(IABS(KAN).EQ.3) THEN
         CALL readArgon()
      ENDIF
caldo
C
C **************************************************************
C     Open probe files, if any 
C **************************************************************
C
      PROBEDATA = "probedataXXX.dat" 
      WRITE(PROBEDATA(10:12),FMT="(I3.3)")MY_PE+1
      INQUIRE(FILE=PROBEDATA,EXIST=LFLAG)
      IF(LFLAG)THEN
         OPEN(IPROBE,FILE=PROBEDATA,STATUS="OLD")
         READ(IPROBE,*)NPROBES
         LPROBE(1) = ISTKGT(NPROBES,KIND_INTEGER)
         LPROBE(2) = ISTKGT(NPROBES*NOFVERT,KIND_REAL8)
         CALL PROBEIN(ISTAK(LPROBE(1)),DSTAK(LPROBE(2)),NPROBES,NOFVERT,
     &                IPROBE)
         CLOSE(IPROBE)
         PROBEOUT = "probeoutXXX.dat" 
         WRITE(PROBEOUT(9:11),FMT="(I3.3)")MY_PE+1
         OPEN(IPROBE,FILE=PROBEOUT,STATUS="UNKNOWN")
      ELSE
         NPROBES = 0
      ENDIF
      IWORK(11) = NPROBES
C
C **************************************************************
C     Open files to dump the integral of the conserved quantities
C **************************************************************
C
      IF(LDUMP(5))THEN
         SOMEFILE = 'integral.' // EXT
         OPEN(IHST4,FILE=SOMEFILE(1:12),FORM='formatted',
     &        STATUS='UNKNOWN')
         WRITE(IHST4,*)'# itim, iter, t, integrals'
      ENDIF
C
C **************************************************************
C     Computing geometrical quantities ..
C **************************************************************
C
C     REM: routine EXGEO() is called with the number of 
C          meshpoints being equal to NPOIN+NGHOST+NPNOD
C
      CALL EXGEO(NELEM,NOFVAR,NPOIN,NGHOST,NPNOD,NPOIN+NGHOST+NPNOD,
     1           NBPOIN,NFACE,NBFAC,NBINT,NWFAC,NBODY4,NBODY6,NHOLE,
     2           NOFVERT,NDIM)
C
      WRITE(NOUT,FMT=9995)NVA
C
C     changes the connectivity of those cells laying close
C     to one of the periodic boundaries
C     Note that all FORTRAN arrays are allocated for NPOIN+NGHOST+NPNOD
C     while PETSc vectors for NPOIN+NGHOST
C     this extra room is used in subr. BACKUP
C     (a little bit tricky)
C
C     observe that it may not be irrelevant to move
C     PERIODIC() BEFORE EXGEO() because of the addressing
C     of periodic gridpoints
C
      IF (PERIODIC_MESH) CALL PERIODIC(NPOIN,NGHOST,NPNOD,NELEM,
     +                   NOFVAR,NOFVERT,ISTAK(LCELNOD),ISTAK(LPMAP))
C
C **************************************************************
C     form PETSc Vectors to store the right hand side (RHS),
C                solution (ZROE) vectors and turbulent viscosity (ViscT)
C **************************************************************
C
#ifdef MPI
C
C
C -- set 0-based indexing; uses block indexing since petsc-3.3 
C
      LGHOST = ISTKGT(NGHOST,KIND_INTEGER)
      CALL GETIDX(NGHOST,1,ISTAK(LTZX),ISTAK(LGHOST))
C
C -- create the vector of the unknowns:
C           memory is provided by the FORTRAN static vector
C           whose pointer is LZROE
C
      CALL VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,NOFVAR,
     2                                  NOFVAR*NPOIN,PETSC_DECIDE,
     3                                  NGHOST,ISTAK(LGHOST),
     +                                  DSTAK(LZROE),ZROE,IFAIL)
C
C     WRITE(6,*)'VecCreateGhostBlockWithArray has ret IOK = ',IFAIL,
C    &' on array Z and PE = ',MY_PE
C
      CALL ISTKRL(1) ! release LGHOST
C
      IF(TURBULENT.AND.SEGREGATED)THEN
C
         LGHOST = ISTKGT(NGHOST,KIND_INTEGER)
         CALL GETIDX(NGHOST,1,ISTAK(LTZX),ISTAK(LGHOST))
C
         CALL VecCreateGhostWithArray(PETSC_COMM_WORLD,NPOIN,
     +                             PETSC_DECIDE,NGHOST,
     +                             ISTAK(LGHOST),DSTAK(LTURB),
     +                             ViscT,IFAIL)
C
         CALL ISTKRL(1) ! release LGHOST
C
      ENDIF
C
#else
C
C    create vectors for: flow variables, turbulent viscosity
C
      CALL VecCreateSeqWithArray(PETSC_COMM_SELF,NOFVAR,NOFVAR*NPOIN,
     +                           DSTAK(LZROE),ZROE,IFAIL)
      IF(TURBULENT.AND.SEGREGATED)
     +CALL VecCreateSeqWithArray(PETSC_COMM_SELF,NTURB,NPOIN,
     +                           DSTAK(LTURB),ViscT,IFAIL)
C
#endif
C
C
c
c  =========================================================
c               Set -restart_file filename to restart
c  =========================================================
c
      IF (ISTART.EQ.0) THEN
C
cxxxx    
C         IF(TURBULENT.AND.COUPLED)THEN
C         u_infty(nofvar) = 1.d0
C         ENDIF
cxxxx   
          J = LZROE - 1
          DO 1 I = 1,NOFVAR
              J = J + 1
              CALL DINIT(NPOIN+NGHOST,U_INFTY(I),DSTAK(J),NOFVAR) ! why not + NPNOD ?!?!?!
    1     CONTINUE
C
          WRITE (NOUT,FMT=115) (U_INFTY(I),I=1,NOFVAR)
C
C     if (TURBULENT.AND.COUPLED) the freestream viscosity is initialized
C         within TURBINI (see below)
C
      ELSE
C
C     in the un-steady case the initial solution is read into
C     U^{n+1} then copied onto U^{n} and U^{n-1}; see below
C
          NITEMS = NPOIN+NGHOST+NPNOD
          CALL SOLZNE(STARFILE,DSTAK(LZROE),NOFVAR,NITEMS,'r')
!     CALL R8Mat_Print('General',' ',Nofvar,(NPOIN+NGHOST+NPNOD)*NTIMLEVS,
!    +    DSTAK(LZROE),Nofvar,'Z array just after solzne ',IFAIL)
      ENDIF
C
C **************************************************************
C     Create the index sets that store informations on
C     the boundary conditions of the boundary nodes
C **************************************************************
C
      LWORKA = ISTKGT(NOFVAR*(NPOIN+NGHOST+NPNOD),KIND_INTEGER)
      LWORKB = ISTKGT(NPOIN+NGHOST+NPNOD,KIND_INTEGER)
      CALL BNDRY_ISET(NodalBcs,ISTAK(LWORKA),ISTAK(LBNDFAC),NBFAC,
     &          ISTAK(LCELNOD),NOFVERT,ISTAK(LWORKB),NDIM,NOFVAR,NPOIN,
     &          NGHOST,NPNOD,NELEM)
      CALL VecGetOwnerShipRange(ZROE,ROWBGN,ROWEND,IFAIL)
      CALL ISET(NodalBcs,ISTAK(LWORKA),ISTAK(LBNDFAC),NBFAC,
     &          ISTAK(LCELNOD),NOFVERT,ISTAK(LWORKB),NDIM,NOFVAR,NPOIN,
     &          NGHOST,NPNOD,NELEM,ROWBGN,(ABS(KAN).EQ.4))
      CALL ISTKRL(2)
C
C **************************************************************
C     Check if there are subsonic inflow boundaries and
C     the corresponding file with the total pressure profile. 
C     REM: NLIST = no. meshpoints on the inflow bnd 
C     *** NLIST is set in the previous call to ISET() ***
C **************************************************************
C
      IF(NLIST.NE.0)THEN
C
C        LKLIST,LVLIST is stored in ibc8.com
C
         LKLIST = ISTKGT(NLIST,KIND_INTEGER)
         LVLIST = ISTKGT(nVarsInlet*NLIST,KIND_REAL8)
         LWORKA = ISTKGT(NLIST,KIND_REAL8)
         LWORKB = ISTKGT(NLIST,KIND_INTEGER)
C
         CALL SETPTO(PTOTFILE,ISTAK(LKLIST),DSTAK(LVLIST),
     +               DSTAK(LWORKA),ISTAK(LWORKB),NDIM,NOFVAR)
C     clear only workarrays
         CALL ISTKRL(2)
      ENDIF
C
c     CALL MPI_Barrier(PETSC_COMM_WORLD,IFAIL) 
c     CALL MPI_Abort(PETSC_COMM_WORLD,-10,IFAIL) 
C
C **************************************************************
C     Set velocity to zero on no-slip surfaces and
C     temperature on isothermal walls; skip if a scalar problem
C **************************************************************
C
      IF(ABS(KAN).EQ.2.OR.ABS(KAN).EQ.4)
     &          CALL INIT(NodalBcs,ISTAK(LBNDFAC),NBFAC,ISTAK(LCELNOD),
     &          NOFVERT,DSTAK(LZROE),DSTAK(LCORG),ROWBGN,NDIM,NOFVAR,
     &          NPOIN,(ABS(KAN).EQ.4),(ISTART.NE.0))
C
C **************************************************************
C     Initialize turbulent stuff (only for KAN = +2,+4)
C **************************************************************
C
C
      IF (TURBULENT) THEN
          IF(SEGREGATED)THEN
              CALL VecGetOwnerShipRange(VISCT,ROWBGN,ROWEND,IFAIL)
          ELSE
              CALL VecGetOwnerShipRange(ZROE,ROWBGN,ROWEND,IFAIL)
          ENDIF
          LWORKA = ISTKGT(NOFVAR*NPOIN,2)
          LWORKB = ISTKGT(NPOIN+NGHOST,2)
          CALL ISET_T(NodalBcs,ISTAK(LWORKA),ISTAK(LBNDFAC),NBFAC,
     &          ISTAK(LCELNOD),NOFVERT,ISTAK(LWORKB),NDIM,NOFVAR,NPOIN,
     &          NGHOST,ROWBGN,(ABS(KAN).EQ.4))
          CALL ISTKRL(2)
          CALL TURBINI(NodalBcs,NOFVERT,NOFVAR,NDIM,NPOIN+NGHOST+NPNOD,
     +                 DSTAK(LTURB),DSTAK(LZROE),DSTAK(LTD),DSTAK(LTTD),
     +                 (ISTART.NE.0),ROWBGN,VISCTFILE,WDISTFILE)
C
      ENDIF
C
C **************************************************************
C     Build the stiffness matrix ...
C **************************************************************
C
      OPT(1) = MAT_KEEP_NONZERO_PATTERN
      OPT(2) = MAT_STRUCTURALLY_SYMMETRIC
      IF (TIMEIMPL ) THEN
c
c we do not use OPT(2) at present
c
          NOPT = 1
#ifdef MPI
          CALL MatAllocaMPI(PetscMatrices(1),NDIM,NOFVERT,NOFVAR,NELEM,
     +                      NPOIN,NGHOST,OPT,NOPT,MAPPING(1),
     +                      'Flow equation matrix')
          IF( TURBULENT .AND. SEGREGATED )
     +    CALL MatAllocaMPI(PetscMatrices(2),NDIM,NOFVERT,NTURB,NELEM,
     +                      NPOIN,NGHOST,OPT,NOPT,MAPPING(2),
     +                      'Turbulence equation matrix')
#else
          CALL MatAllocaSeq(PetscMatrices(1),NDIM,NOFVERT,NOFVAR,NELEM,
     +                      NPOIN,NGHOST,OPT,NOPT,
     3                      'Flow equation matrix')
          IF( TURBULENT .AND. SEGREGATED )
     +    CALL MatAllocaSeq(PetscMatrices(2),NDIM,NOFVERT,NTURB,NELEM,
     +                      NPOIN,NGHOST,OPT,NOPT,
     &                      'Turbulence equation matrix')
#endif
C
      ENDIF 
      IF ( LALE .AND. ALE_LAPLACE_GRID ) THEN
          NOPT = 1
#ifdef MPI
          CALL MatAllocaMPI(PetscMatrices(3),NDIM,NOFVERT,1,NELEM,NPOIN,
     +                      NGHOST,OPT,NOPT,MAPPING(3),
     +                      'Motion solver matrix')
#else
          CALL MatAllocaSeq(PetscMatrices(3),NDIM,NOFVERT,1,NELEM,NPOIN,
     +                      NGHOST,OPT,NOPT,'Motion solver matrix')
#endif
      ENDIF 
C
C **************************************************************
C     Checking memory usage
C **************************************************************
C
!     CALL PetscMemoryGetCurrentUsage(MAXSIZ,IFAIL)
!     WRITE (NOUT,FMT=9994) MAXSIZ
      WRITE (NOUT,FMT=9998) ISTKST(3)/2,
     +          (ISTKST(3)*100.)/ISTKST(4)
C
      IWORK(1) = NDIM
      IWORK(2) = NOFVERT
      IWORK(3) = NOFVAR
      IWORK(4) = NTURB
      IWORK(5) = NPOIN
      IWORK(6) = NGHOST
      IWORK(7) = NELEM
      IWORK(8) = NFACE
      IWORK(9) = NBFAC
      IWORK(10) = NPNOD
      IWORK(12) = NHOLE
      IWORK(13) = NBINT
C
      IF(LTIME)THEN 
C
C
C     Copy current time level into the previous ones
C     this should be done after all initializations
C     have been performed
C     Z(LZROE) contains u^{0} which is then copied
C     onto u^{n} and u^{n-1}
C     when re-starting reads time level n-1
C
         DO I = 2,MAXTIMLEVS
            NITEMS = (NPOIN+NGHOST+NPNOD)*NOFVAR
            J = LZROE + (I-1)*NITEMS
            CALL DCOPY(NITEMS,DSTAK(LZROE),1,DSTAK(J),1)
         ENDDO
         IF( LALE )THEN ! save a copy of the original mesh
            NITEMS = (NPOIN+NGHOST+NPNOD)*NDIM
            CALL DCOPY(NITEMS,DSTAK(LCORG),1,DSTAK(LCORG+2*NITEMS),1)
C
C     Copy the current volumes into those at time levels n-1 n-2
C     LVOL        |---- n   ----| |---- n ----| |---- n-1 ----| |---- n-2 ----|
C
            IF( (NTIMLEVS.EQ.3) )THEN
               DO I = 1,MAXTIMLEVS
                  J = LVOL + I*NELEM
                  CALL DCOPY(NELEM,DSTAK(LVOL),1,DSTAK(J),1)
               ENDDO
            ENDIF
C        Open log file
              OPEN(UNIT=IOALE,FILE='ale.log',FORM='formatted')
         ENDIF ! test on LALE
C
C
C     time level (n-1) is read from a datafile and then
C     copied into n; it will then
C     be shifted into (n-1) at the beginning of the time loop
C
         IF((ISTART.NE.0).AND.(NTIMLEVS.EQ.3))THEN
             INQUIRE(FILE=FOLDFILE,EXIST=LFOOLD)
             IF(LFOOLD)THEN
                NITEMS = NPOIN+NGHOST+NPNOD
                J  = LZROE+NOFVAR*NITEMS ! points to the solution at (n)
                CALL SOLZNE(FOLDFILE,DSTAK(J),NOFVAR,NITEMS,'r')
caldo
#ifdef XPRINT
            NITEMS = (NGHOST+NPNOD+NPOIN)
            DO ILEV = 1,MAXTIMLEVS
                K = LZROE+(ILEV-1)*NITEMS*NOFVAR
                DO J = 0,NOFVAR-1
                   S= DNRM2(NITEMS,DSTAK(K+J),NOFVAR)
                   write(6,*)' ILEV = ',ILEV ,' IVAR = ',J,' S = ',S
                ENDDO
            ENDDO
#endif
caldo
             ELSE ! there is no restart file for time level n-1
                WRITE(NOUT,FMT=105)
             ENDIF
         ENDIF ! on ISTART
      ENDIF ! on LTIME
C
C     Pseudo-time stepping
C
      CALL myTS(ZROE,VISCT,PetscMatrices,NodalBcs,IWORK)
C
      CALL PetscTime(tend,IFAIL)
      totime = totime + tend-tbeg
C
C  =========================================================
C
      LWORKD = ISTKGT(2*NBFAC,2)
      CALL WSKIN(DSTAK(LSKINF),DSTAK(LHEAT),ISTAK(LWORKD),
     +           NWFAC,ISTAK(LBNDFAC),NBODY6,SKINFILE)
      CALL ISTKRL(1)
C
      LWORKD = ISTKGT(MAX(NOFVAR,NTURB)*NPNOD,2)
      CALL BACKUP(NPOIN,NPNOD,NGHOST,NOFVAR,NTURB,ISTAK(LPMAP),
     +            ISTAK(LWORKD),BAKFILE,VISCTFILE,FOLDFILE)
      IF(IABS(KAN).NE.1)
     1       CALL CHKFLX(ISTAK(LBNDFAC),NBFAC,ISTAK(LCELFAC),
     2                   ISTAK(LCELNOD),NOFVERT,NELEM,DSTAK(LFACNOR),
     3                   DSTAK(LZROE),DSTAK(LXYZDOT),NOFVAR,NDIM,NFACE,
     4                   (IABS(KAN).EQ.4),LDUMP(4),ITER)
caldo        CALL AVGFLX(ISTAK(LBNDFAC),NBFAC,ISTAK(LCELFAC),
caldo2                   ISTAK(LCELNOD),NOFVERT,NELEM,DSTAK(LFACNOR),
caldo3                   DSTAK(LZROE),NOFVAR,NDIM,NFACE ,
caldo4                   (IABS(KAN).EQ.4),ITER)
      CALL ISTKRL(1)
C
!     CALL PetscMemoryGetMaximumUsage(MAXSIZ,IFAIL)
      WRITE (NOUT,FMT=9998) ISTKST(3)/2, (ISTKST(3)*100.)/ISTKST(4)
!     WRITE (NOUT,FMT=9994) MAXSIZ
C
caldo CLOSE(NOUT)
caldo CLOSE(IWUNIT)
C
      CALL clearmem( ZROE, VISCT, PetscMatrices, NodalBcs, NOFVAR )
C
      CALL PetscFinalize(IFAIL)
C
      CLOSE(UNIT=ITIM1)
      IF(LDUMP(5))CLOSE(IHST4)
      IF(MY_PE.EQ.0)THEN
          WRITE(IWUNIT,FMT=*)'time spend ',totime
          CLOSE(IHST1)
          CLOSE(IHST2)
          IF(TURBULENT .AND. SEGREGATED)CLOSE(IHST3)
      ENDIF 
C
C *************** I/O FORMATS *****************
C
  105 FORMAT (5X,'Solution at n=-1 NOT found; using current solution',/)
  115 FORMAT (5X,'Initial conditions are :',/,5X,8 (E12.6,2X),/)
 9995 FORMAT (/,5X,'MEMORY ALLOCATED IN DSTAK ',I9,' REAL*8 WORDS')
 9998 FORMAT (/,5X,'MAXIMUM MEMORY USED IN DSTAK ',I9,
     +       ' REAL*8 WORDS (',F5.2,' %)')

      END
