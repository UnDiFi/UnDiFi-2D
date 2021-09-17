      SUBROUTINE SOLZNE(FILENAME,VARRAY,NOFVAR,NPOIN,MODE)
C
C
C Subroutine for reading(mode='r') and writing(mode='w') nodal values
C
      IMPLICIT NONE
C
C
C
C
C
C
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NPOIN
      CHARACTER FILENAME* (*),MODE* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VARRAY(NOFVAR,NPOIN)
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,IXDRS,NPOLD,NVOLD
C     ..
C     .. External Functions ..
      INTEGER INITXDR
      INTEGER IXDRINT,IXDRIMAT,IXDRCLOSE,IXDRDMAT
      EXTERNAL INITXDR
C     ..
C     .. External Subroutines ..
      EXTERNAL SETERR,IXDRINT,IXDRIMAT,IXDRCLOSE,IXDRDMAT
C     ..
C     .. Data statements ..
      DATA IFAIL/0/
C     ..
C
C
C     .. Reading or Backing up ..
C
      IF (MODE.EQ.'r') THEN
          WRITE (6,FMT=110) FILENAME

      ELSE
          WRITE (6,FMT=112) FILENAME
      ENDIF
C
      IXDRS = INITXDR(FILENAME,MODE,.false.)
C
      NPOLD = NPOIN
      NVOLD = NOFVAR
      IFAIL = IXDRINT(IXDRS,NPOIN)
      IFAIL = IXDRINT(IXDRS,NOFVAR)
      IF (MODE.EQ.'r') THEN
          IF (NPOIN.NE.NPOLD) CALL SETERR
     +                             (30HINCONSISTENT NPOIN IN DATAFILE,
     +                             30,1,2)
          IF (NOFVAR.NE.NVOLD) CALL SETERR
     +                              (31HINCONSISTENT NOFVAR IN DATAFILE,
     +                              31,1,2)
      ENDIF

      IFAIL = IXDRDMAT(IXDRS,NOFVAR*NPOIN,VARRAY)
caldo CALL X04CAF('General',' ',NOFVAR,NPOIN,VARRAY,NOFVAR,
caldo+            'Nodal values',IFAIL)
      IFAIL = IXDRCLOSE(IXDRS)
C
      RETURN
C
  110 FORMAT (/,5X,'Reading solution from ',A40,/)
  112 FORMAT (/,5X,'Writing solution to ',A40,/)
  115 FORMAT (' done',/)

      END
