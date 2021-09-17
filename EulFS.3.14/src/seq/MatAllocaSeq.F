      SUBROUTINE MatAllocaSeq(A,NDIM,NOFVERT,NBLK,NELEM,NR,NGHOST,
     +                        OPT,NOPT,TITLE)
C
C     $Id: MatAllocaSeq.F,v 1.19 2020/04/23 09:37:34 abonfi Exp $
C
CCCCC#define DEBUG
CCCCC#define CSR_PRINT
C
C
C     This routine creates the Compressed Sparse Row structure
C        of the stiffness matrix, i.e. the arrays:
C        ia(1:NR+1) and ja(1:NNZ)
C
C     $Id: MatAllocaSeq.F,v 1.19 2020/04/23 09:37:34 abonfi Exp $
C
C
#include "petsc/finclude/petscmat.h"
CCCC#include "petsc/finclude/petscvec.h"
CCCC#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscis.h"
              use petscis
              use petscmat
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'implicit.h'
      INCLUDE 'verbose.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'nloc.com'
C
C     .. Scalar Arguments ..
      INTEGER NBLK,NDIM,NELEM,NOFVERT,NR,NGHOST,NOPT
      CHARACTER*(*) TITLE
C
C     .. Petsc Arrays ..
C     ..
      Mat A
      MatOption OPT(NOPT)
C
C     Input:
C     -----
C     NR          numerb of (block) rows in the matrix
C     NBLK        block dimension of the matrix
C     NGHOST      number of ghost points on this processor
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,LDEGREE,NNZR,N,IOPT,IBGN,IEND,I
      DOUBLE PRECISION HELP
#ifdef DEBUG
      INTEGER LIWK,LJWK,LIA,LIAO,LJA,LJAO,JOB,IOUT,IPOS
#endif
      LOGICAL VERBOSE
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
      DOUBLE PRECISION INFO(MAT_INFO_SIZE)
#ifdef DEBUG
      DOUBLE PRECISION AMAT(1),BMAT(1)
#endif
C     ..
C     .. External Functions ..
      INTEGER ISTKGT,ISTKST
      EXTERNAL ISTKGT,ISTKST
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT,ISTKRL
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
C     ..
C
C     ------ Allocate space for the stiffness matrix -------
C
      VERBOSE = (IVERBOSE.GT.0)
      VERBOSE = .TRUE.
C
      LDEGREE = ISTKGT(NR+1,2)
C
C     MatGetSizeSeq finds the number of non zero entries (NNZR)
C     of the stiffness matrix and returns in DEGREE(i) the number
C     of non zero entries in the i-th row of the stiffness matrix ..
C
      CALL MatGetSizeSeq(ISTAK(LCELNOD),ISTAK(LDEGREE),ISTAK(LCELCEL),
     +              ISTAK(LNODCOD),NDIM,NOFVERT,NELEM,NR,NNZR,NOUT,
     +              VERBOSE,IFAIL)
caldo
!     IFAIL = 3
caldo
C
C     If MatGetSizeSeq fails, we use tet_mesh_order4_adj_count
C
      IF(IFAIL.NE.0)THEN
      WRITE(NOUT,FMT=300)IFAIL
         IF( NDIM .EQ. 3)THEN
            WRITE(NOUT,FMT=302)"tet_mesh_order4_adj_count"
            CALL tet_mesh_order4_adj_count ( NR, NELEM, ISTAK(LCELNOD),
     &        NNZR, ISTAK(LDEGREE) )
         ELSEIF( NDIM .EQ. 2 )THEN
            GOTO 44
#ifdef CSR_PRINT
            OPEN(87,FILE="csr_1.log")
            do I = 1,NR
               WRITE(87,FMT=*)I,ISTAK(LDEGREE+I-1)
            ENDDO
            CLOSE(87)
#endif
caldo
caldo looks like triangulation_order3_adj_count does NOT work properly
caldo tested on 26.01.2019
caldo
            CALL EXIT(IFAIL)
caldo
caldo triangulation_order3_adj_count has been modified in such a way that
caldo icelcel is set to a negative value whenever icelcel > NELEM
caldo
            WRITE(NOUT,FMT=302)"triangulation_order3_adj_count"
            CALL triangulation_order3_adj_count ( NR,
     &      NELEM, ISTAK(LCELNOD), ISTAK(LCELCEL), NNZR,
     &      ISTAK(LDEGREE) )
         ENDIF ! test on NDIM
C
#ifdef CSR_PRINT
         OPEN(87,FILE="csr_2.log")
#endif
         DO I = LDEGREE,LDEGREE+NR-1
            IBGN = ISTAK(I)
            IEND = ISTAK(I+1)
            N = IEND-IBGN+1            
            ISTAK(I) = N
         ENDDO
         IBGN = 1000000
         IEND = 0
         HELP = ZERO
         DO I = LDEGREE,LDEGREE+NR-1
            N = ISTAK(I)
#ifdef CSR_PRINT
            WRITE(87,FMT=*)I-LDEGREE+1,N
#endif
            HELP = HELP + REAL(N)/REAL(NR)
            IBGN = MIN(IBGN,N)
            IEND = MAX(IEND,N)
         ENDDO
#ifdef CSR_PRINT
         CLOSE(87)
#endif
         WRITE(NOUT,FMT=320)IBGN,IEND,HELP,NNZR
      ENDIF ! test on IFAIL
C
#ifdef DEBUG
C
C     BAIJ matrix
C
      lia = istkgt(nr+1,2)
      lja = istkgt(nnzr,2)
      liwk = istkgt(nr+1,2)
      ljwk = istkgt(nnzr,2)
c
c     copy degree into ia, since the routine csrfmt needs
c     needs the info stored in degree
c
      call icopy(nr+1,istak(ldegree),1,istak(lia),1)
c
      call CSRFMT(NR,NELEM,NOFVERT,istak(lcelnod),
     >istak(lia),istak(lja),istak(liwk),istak(ljwk))
c
c     release liwk ljwk
c
      call istkrl(2)
c
c     AIJ matrix
c
      liao = istkgt(nblk*(nr+1),2)
      ljao = istkgt((nblk**2)*nnzr,2)
c
c     values are NOT copied
c
      job=0
c
      IF( NBLK .EQ. 1 )THEN
         ipos = 1
         call copmat (nr,amat,ISTAK(LJA),ISTAK(LIA),bmat,
     &                ISTAK(Ljao),ISTAK(Liao),ipos,job)
      ELSE 
c
c     convert BAIJ into AIJ (not in place)
c
c
         call bsrcsr (job, nr, nblk, nblk, amat, istak(lja), istak(lia),
     >bmat, istak(ljao), istak(liao) )
      ENDIF
c

      iout=100
      write(6,*)'Writing sparsity pattern to matAIJ'
      open(unit=iout,file='matAIJ',status='unknown')
c
c     dump the sparsity pattern of the AIJ matrix
c
      call dump (1,nr*nblk,.false.,amat,istak(ljao),istak(liao),iout)
      write(6,*)'....done'
      close(iout)
C
      call istkrl(4)
!     stop
C
#endif
C
   44 CONTINUE
      N=NR*NBLK
C
      IF( NBLK .EQ. 1 )THEN
          CALL MatCreateSeqAIJ(PETSC_COMM_SELF,NR,NR,NNZR/NR,
     +                         ISTAK(LDEGREE),A,IFAIL)
      ELSE
          CALL MatCreateSeqBAIJ(PETSC_COMM_SELF,NBLK,N,N,NNZR/NR,
     +                          ISTAK(LDEGREE),A,IFAIL)
      ENDIF
C
C     Blocks to be inserted in A will be COLUMN ORIENTED
C
      WRITE(NOUT,348)
      WRITE(NOUT,349)0,MAT_ROW_ORIENTED
      CALL MatSetOption(A,MAT_ROW_ORIENTED,PETSC_FALSE,IFAIL) 
      DO 12 IOPT = 1, NOPT
          WRITE(NOUT,349)IOPT,OPT(IOPT)
          CALL MatSetOption(A,OPT(IOPT),PETSC_TRUE,IFAIL) 
   12 CONTINUE
      CALL MatZeroEntries(A,IFAIL)
C
C     free DEGREE
C
      CALL ISTKRL(1)
C
      CALL MatGetInfo(A,MAT_LOCAL,INFO,IFAIL)
      WRITE(NOUT,346)TITLE
C
      WRITE(NOUT,347)N,N,
     +            INT(INFO(MAT_INFO_BLOCK_SIZE)),
     +            INT(INFO(MAT_INFO_NZ_ALLOCATED)),
     +            INFO(MAT_INFO_MEMORY)
C
      RETURN
C
  300 FORMAT(/,/,5X,'SUBR. MatGetSizeSeq has RETURNED IFAIL = ',I6)
  302 FORMAT(/,/,5X,'Calling subr. ',A,/,/)
! 305 FORMAT(/,/,5X,'SUBR. MatGetSizeSeq has RETURNED IFAIL = ',I6,/,
!    &15X,'using triangulation_order3_adj_count',/,/)
  320 FORMAT (/,/,'MatGetSizeSeq; STIFFNESS MATRIX',/,' ',15 ('='),/,
     +       10X,'MIN/MAX/AVG VERTEX DEGREE : ',2 (I2,2X),F4.1,/,10X,
     +       'NNZR = ',I12,' NONZERO BLOCK ENTRIES',/)
  346 FORMAT(10X,5("*"),A)
  347 FORMAT(10X,'Number of rows and columns    ',2(I8,1X),/,
     +       10X,'Block size                    ',I2,/,
     +       10X,'Number of allocated  nonzeros ',I14,/,
     +       10X,'Memory allocated              ',F16.0)
  348 FORMAT(10X,'Setting matrix options :',/)
  349 FORMAT(10X,'Matrix option # ',I2,' is ',I2)
C
      END  
