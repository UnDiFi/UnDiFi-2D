      SUBROUTINE MatAllocaMPI(A,NDIM,NOFVERT,NBLK,NELEM,NR,NGHOST,
     +                        OPT,NOPT,MAPPING,TITLE)
C
C Purpose:
C          the driver for allocating memory and building the
C          stiffness matrix in parallel
C
C
C     $Id: MatAllocaMPI.F,v 1.35 2021/01/29 09:01:29 abonfi Exp $
C
C
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
      use petscvec 
      use petscmat 
      use petscis 
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'bnd.h'
      INCLUDE 'nloc.com'
      INCLUDE 'verbose.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
C
C
C     .. Scalar Arguments ..
      INTEGER NBLK,NDIM,NELEM,NOFVERT,NR,NGHOST,NOPT
      CHARACTER*(*) TITLE
C
C     .. Petsc Stuff ..
C     ..
      Mat A
      MatOption OPT(NOPT)
      ISLocalToGlobalMapping MAPPING
      PetscBool IFLAG
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER MY_PE
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,N,kspace,NERR,IOPT
      INTEGER ROWBGN,ROWEND,I,J,LINDX,IXDRS,LD_NNZ,LO_NNZ
      LOGICAL VERBOSE,FNAME
      CHARACTER MATRFILE*256,DATADIR*256,EXT*3,ERRMSG*72
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
      DOUBLE PRECISION INFO(MAT_INFO_SIZE)
C     ..
C     .. External Functions ..
      INTEGER ISTKGT,ISTKST,initxdr,ixdrint,ixdrimat,IXDRCLOSE
      EXTERNAL ISTKGT,ISTKST,initxdr,ixdrint,ixdrimat,IXDRCLOSE
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT,ISTKRL
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON/MPICOM/MY_PE
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
C     ..
C
C
C
C     ------ Allocate space for the stiffness matrix -------
C
      VERBOSE = (IVERBOSE.GT.0)
      VERBOSE = .TRUE.
C
C     the sparsity pattern of the parallel matrix is read from
C     files matr$num_proc$.dat
C     If these are not found; program terminates
C
C
      CALL PetscOptionsGetString(PETSC_NULL_OPTIONS,
     &     PETSC_NULL_CHARACTER,'-data_dir',DATADIR,IFLAG,IFAIL)
      DO 70 I=1,80
          MATRFILE(I:I) = CHAR(0)
   70 CONTINUE
      WRITE(EXT,FMT='(I3.3)')MY_PE+1
      KSPACE = INDEX(DATADIR," ")
      MATRFILE(1:KSPACE+10) = DATADIR(1:KSPACE-1)//'matr'//ext//'.dat'
      KSPACE = KSPACE+10
C
      INQUIRE(FILE=MATRFILE(1:KSPACE),EXIST=FNAME)
      IF(FNAME)THEN
         IXDRS = INITXDR(MATRFILE(1:KSPACE),'r',.FALSE.) 
         LD_NNZ = ISTKGT(NR,2)
         LO_NNZ = ISTKGT(NR,2)
         IFAIL = IXDRINT(IXDRS,N)
         IF(N.NE.NR)THEN
            WRITE(ERRMSG(1:72),500)
  500 FORMAT('MatAllocaMPI Number of rows in matr$$$ file does NOT match
     >',2(1X,I6))
            NERR = 3
            IOPT = 1
            CALL SETERR(ERRMSG,NERR,IOPT)
         ENDIF
C
C     read diagonal and off-diagonal CSR structure
C
         IFAIL = IXDRIMAT(IXDRS,NR,ISTAK(LD_NNZ))
         IFAIL = IXDRIMAT(IXDRS,NR,ISTAK(LO_NNZ))
C     close file
         IFAIL = IXDRCLOSE(IXDRS)
      ELSE
         WRITE(6,520)MY_PE,MATRFILE(1:KSPACE)
         WRITE(ERRMSG(1:37),510)MY_PE
         WRITE(6,*)ERRMSG(1:37)
         WRITE(6,*)'Exist / PE ',FNAME,MY_PE
  520 FORMAT('PE # ',I4,' File : ',A120)
  510 FORMAT('MatAllocaMPI PE# ',I4,'CANNOT OPEN FILE')
            NERR = 9
            IOPT = 1
            CALL SETERR(ERRMSG,NERR,IOPT)
      ENDIF
C
C
      N=NR*NBLK
C
      IF( NBLK .EQ. 1 )THEN
         CALL MatCreateAIJ(PETSC_COMM_WORLD,NR,NR,PETSC_DETERMINE,
     +         PETSC_DETERMINE,PETSC_DEFAULT_INTEGER,ISTAK(LD_NNZ),
     +         PETSC_DEFAULT_INTEGER,ISTAK(LO_NNZ),A,IFAIL)
      ELSE
         CALL MatCreateBAIJ(PETSC_COMM_WORLD,NBLK,N,N,
     +                PETSC_DETERMINE,PETSC_DETERMINE,
     +                PETSC_NULL_INTEGER,ISTAK(LD_NNZ),
     +                PETSC_NULL_INTEGER,ISTAK(LO_NNZ),A,IFAIL)
      ENDIF
C
C     create Local to Global mapping
C     for the parallel matrix
C
      CALL MatGetOwnershipRange(A,ROWBGN,ROWEND,IFAIL)
C
C     note that N = ROWEND-ROWBGN
C
C     the indexing is blocked of block size NBLK
C
      LINDX = ISTKGT((NR+NGHOST),2)
C
C     mapping for the processor owned (or local)  rows
C
      J = 0
      DO 55 I = ROWBGN,ROWEND-1,NBLK
         ISTAK(LINDX+J) = I/NBLK
         J = J + 1
   55 CONTINUE
C
C     mapping for ghost rows: GetIdx copies the ghost locations 
C                             into the array INDX and sets the
C                             0-based indexing 
C
      CALL GETIDX(NGHOST,1,ISTAK(LTZX),ISTAK(LINDX+NR))
C
      CALL ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD,NBLK,NR+NGHOST,
     &                                  ISTAK(LINDX),PETSC_COPY_VALUES,
     3                                  mapping,IFAIL)
C
C     free LINDX
C
      CALL ISTKRL(1)
C
C     Sets local to global mapping:
C
          CALL MatSetLocalToGlobalMapping(A,mapping,mapping,IFAIL) 
C
C
C     Set options (if any)
C
      DO 12 IOPT = 1, NOPT
          CALL MatSetOption(A,OPT(IOPT),PETSC_TRUE,IFAIL) 
   12 CONTINUE
C
C     Blocks to be inserted in A will be COLUMN ORIENTED
C
      CALL MatSetOption(A,MAT_ROW_ORIENTED,PETSC_FALSE,IFAIL) 
      CALL MatZeroEntries(A,IFAIL)
C
C     ***************  Local Matrix  ***************
C
      CALL MatGetInfo(A,MAT_LOCAL,INFO,IFAIL)
      WRITE(NOUT,346)TITLE
C
      WRITE(NOUT,347)N,N,
     +               INT(INFO(MAT_INFO_BLOCK_SIZE)),
     +               INT(INFO(MAT_INFO_NZ_ALLOCATED)),
     +               INFO(MAT_INFO_MEMORY)
C
C     ***************  Global Matrix  ***************
C
      CALL MatGetInfo(A,MAT_GLOBAL_SUM,INFO,IFAIL)
      IF(MY_PE.EQ.0)THEN
          WRITE(IWUNIT,346)TITLE
C
          WRITE(IWUNIT,347)N,N,
     +                     INT(INFO(MAT_INFO_BLOCK_SIZE)),
     +                     INT(INFO(MAT_INFO_NZ_ALLOCATED)),
     +                     INFO(MAT_INFO_MEMORY)
      ENDIF
C
      RETURN
C
  346 FORMAT(10X,5("*"),A)
  347 FORMAT(10X,'Number of (block) rows and columns ',2(I6,1X),/,
     +       10X,'Block size                    ',I2,/,
     +       10X,'Number of allocated  nonzeros ',I12,/,
     +       10X,'Memory allocated              ',F16.0)
C
      END  
