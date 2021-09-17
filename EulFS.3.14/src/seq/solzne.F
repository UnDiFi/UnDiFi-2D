!>
!> \par Purpose
!>
!> Read/write a solution-type array from/to a XDR file using call to library FXDR
!>
!> @param[in] FILENAME name of the XDR file where data are read from or written to
!> @param[in] VARRAY the array to be read or written
!> @param[in] NOFVAR number of dofs
!> @param[in] NPOIN number of gridpoints
!> @param[in] MODE starts with 'R' or 'r' for reading and with 'w' or 'W' for writing
      SUBROUTINE SOLZNE(FILENAME,VARRAY,NOFVAR,NPOIN,MODE)
!> \author $Author: abonfi $
!> \version $Revision: 1.13 $
!> \date $Date: 2020/03/28 09:42:53 $
C
C Subroutine for reading(mode='r') and writing(mode='w') nodal values
C
      IMPLICIT NONE
      INCLUDE 'io.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NPOIN
      CHARACTER FILENAME* (*),MODE* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VARRAY(NOFVAR,NPOIN)
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,IXDRS,NPOLD,NVOLD,LENFN
      LOGICAL RETOERR
!     PARAMETER(RETOERR=.TRUE.) 
      PARAMETER(RETOERR=.FALSE.) 
c     RETOERR returnonerror: logical: if .TRUE., then
c                        routines return even if there is
c                        an I/O error.  If .FALSE., then
c                        routines halt on I/O error,
c                        printing out a useful message
c                        (this is like what fortran does, so set to .FALSE. for most fortran-
c                        like behavior)
      CHARACTER*72 ERRMSG
C     ..
C     .. External Functions ..
      INTEGER INITXDR,IXDRCLOSE,IXDRINT,I1MACH
      EXTERNAL INITXDR,IXDRCLOSE,IXDRINT,I1MACH
#ifdef PARCH_t3d
      INTEGER IXDRRMAT64
      EXTERNAL IXDRRMAT64
#else
      INTEGER IXDRDMAT
      EXTERNAL IXDRDMAT
#endif
C     ..
C     .. External Subroutines ..
      EXTERNAL SETERR
C     ..
C     .. Data statements ..
      DATA IFAIL/0/
C     ..
C
C
C     .. Reading or Backing up ..
C
      LENFN = INDEX( FILENAME, CHAR(0) ) -1
      IF(LENFN.LE.1)LENFN = 128
      IF (MODE(1:1).EQ.'r'.OR.MODE(1:1).EQ.'R') THEN
          WRITE (NOUT,FMT=110) FILENAME(1:LENFN)

      ELSEIF (MODE(1:1).EQ.'w'.OR.MODE(1:1).EQ.'W') THEN
          WRITE (NOUT,FMT=112) FILENAME(1:LENFN)
      ELSE
          STOP 'Unknown mode in SUBROUTINE solzne' 
      ENDIF
C
      IXDRS = INITXDR(FILENAME,MODE,RETOERR)
      IF( IXDRS.LT.0 )THEN
          WRITE(ERRMSG,FMT=2035)FILENAME(1:MIN(LENFN,42))
          CALL SETERR(ERRMSG,72,IXDRS,2)
!         WRITE(I1MACH(4),2020)'INITXDR',IXDRS,'SOLZNE',MY_PE
!         WRITE(I1MACH(4),2030)FILENAME
!         CALL EXIT(IXDRS)
      ENDIF
C
C
      NPOLD = NPOIN
      NVOLD = NOFVAR
      IFAIL = IXDRINT(IXDRS,NPOIN)
      IFAIL = IXDRINT(IXDRS,NOFVAR)
      IF (MODE(1:1).EQ.'r'.OR.MODE(1:1).EQ.'R') THEN
          IF (NPOIN.NE.NPOLD) THEN
             WRITE(ERRMSG,FMT=2040)FILENAME(1:MIN(LENFN,50))
             CALL SETERR (ERRMSG,72,1,2)
          ENDIF
          IF (NOFVAR.NE.NVOLD) THEN
          WRITE(6,*)'Expecting ',NVOLD,' found : ',NOFVAR
          CALL SETERR (31HINCONSISTENT NOFVAR IN DATAFILE,31,1,2)
          ENDIF
      ENDIF

#ifdef PARCH_t3d
      IFAIL = IXDRRMAT64(IXDRS,NOFVAR*NPOIN,VARRAY)
#else
      IFAIL = IXDRDMAT(IXDRS,NOFVAR*NPOIN,VARRAY)
#endif
C     CALL R8Mat_Print('General',' ',NOFVAR,NPOIN,VARRAY,NOFVAR,
C    +            'Nodal values',IFAIL)
      IFAIL = IXDRCLOSE(IXDRS)
C
      RETURN
C
  110 FORMAT (/,5X,'Reading solution from:'/,A,/)
  112 FORMAT (/,5X,'Writing solution to:'/,A,/)
  115 FORMAT (' done',/)
C
 2020 FORMAT (5X,'FUNCTION ',A31,/,
     +5X,'HAS RETURNED A NONZERO ERROR FLAG =',I6,/,
     +5X,'IN SUBROUTINE ',A10,' ON PE =',I4)
 2030 FORMAT (5X,'WHILE TRYING TO OPEN FILE:',/,A75)
!1234567890123456789012345678901234567890123456789012345678901234567890123456789012
 2035 FORMAT('INITXDR failed while reading ',A42)
 2040 FORMAT('INCONSISTENT NPOIN IN ',A50)
C
      END
