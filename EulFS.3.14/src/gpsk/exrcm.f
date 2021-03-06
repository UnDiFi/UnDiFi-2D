      SUBROUTINE EXRCM(NROWS,IA,JA,PERMUT,DEGREE,RSTART,
     +CONNEC,WORK,WRKLEN)
C
C     Applies an RCM reordering to a mesh
C     RCM=Reverse Kuthill-McKee
C
C     The RCM algorithm has been taken from
C     ACM TOMS
C     Some routines are from Y.Saad's SPARSKIT
C     library
C     Other routines are from ITPACK
C
C     needs ain input file called inp
C     which specifies the mesh files 
C     the output is a new, reorderes
C
      IMPLICIT NONE 
C
C
C     .. External Subroutines ..
C
C     ..
C     .. External Functions ..
C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INT,LOG10,MAX,MIN,REAL
C     ..
C     .. Equivalences ..
C     ..
C     .. Local Arrays ..
C     ..
C     .. Parameters ..
C     ..
C     ..
C     .. Local Scalars ..
      REAL*8 A
      INTEGER LWORK,ISYM,NROWS,WRKLEN
      INTEGER CONNEC(*),DEGREE(*),I,IFAIL,
     +        LEN,JOB,PERMUT(*),WORK(WRKLEN),RSTART(*),
     +        IA(*),JA(*),NNZR
C     ..
C     .. Local Scalars ..
      INTEGER BANDWD,ERROR,IOFF,LEVEL,PROFIL,SPACE
      LOGICAL OPTPRO
C
      CALL IINIT(NROWS,0,PERMUT,1)
C
C     arrays WRKLEN,DEGREE,RSTART,CONNEC are described
C     in subroutine GPSKCA
C     RSTART,CONNEC are basically IA,JA with
C     the diagonal elements removed
C
      NNZR = IA(NROWS+1)-IA(1)
      write(6,*)'nrows, nnzr = ',nrows, nnzr
      write(6,*)'wrklen = ',wrklen
      CALL ICOPY(NROWS+1,IA,1,RSTART,1)
      CALL ICOPY(NNZR,JA,1,CONNEC,1)
      CALL IINIT(NROWS,0,DEGREE,1)
      CALL IINIT(WRKLEN,0,WORK,1)
C
C job   = integer. job indicator.  if job = 0 then
C         the matrix a, ja, ia, is not altered on return.
C         if job.ne.0  then getdia will remove the entries
C         collected in diag from the original matrix.
C         this is done in place.
C
      JOB = 100
      IOFF = 0
C
C  remove diagonal entries
C
      WRITE(6,*)'Removing diagonal entries '
C
C     N.B. gpskca requires that diagonal elements
C     are NOT in the connectivity array;
C     we therefore remove them from
C
      CALL GETDIA2(NROWS,NROWS,JOB,CONNEC,RSTART,LEN,
     +             WORK,IOFF)
C
      WRITE(6,*)'Done'
      WRITE(6,*)'Computing RCM reordering'
C
      DO 120 I= 1,NROWS
          PERMUT(I) = I
          DEGREE(I) = RSTART(I+1)-RSTART(I)
  120 CONTINUE
C
C
C
C     This section drives the GPSK algorithm from ACM
C
C
      OPTPRO = .FALSE.
C
      CALL GPSKCA(NROWS,DEGREE,RSTART,CONNEC,OPTPRO,WRKLEN,
     +            PERMUT,WORK,BANDWD,PROFIL,ERROR,SPACE)
C
      WRITE (6,FMT=*) 'GPSKCA has returned error/space ,',ERROR,SPACE
      WRITE (6,FMT=*) 'GPSKCA has returned BANDWD/PROFIL ,',
     +BANDWD,PROFIL
C
      CALL M01ZAF(PERMUT,1,NROWS,ERROR)
C
      WRITE(6,*)'Done'
C
      END
