      SUBROUTINE WKSPCE(DIM,NPOIN,NPNOD,NELEM,NBFAC,NHOLE,NOFVAR,
     +NOFVERT)
C
      IMPLICIT NONE
C
C
      INTEGER*4 NIN,NOUT
      PARAMETER (NIN=5,NOUT=6)
C
C	NIN	is the OUTPUT device number
C	NOUT	is the INPUT device number
C
      INTEGER ixdrs(10)
      CHARACTER*80 filename(5)
      CHARACTER*72 ERRMSG
      COMMON/ES/ixdrs,filename
C
      INTEGER DIM,NPOIN,NPNOD,NELEM,NBFAC,NHOLE,NOFVAR,NOFVERT
      INTEGER IFAIL,NITEMS
      LOGICAL LFLAG
C
      INTEGER  INITXDR,IXDRINT,IXDRCLOSE
      EXTERNAL INITXDR,IXDRINT,IXDRCLOSE
C
C     .. Reading info from the mesh file ..
C
      write(6,*)' filename(1) = ',filename(1)
      ixdrs(1) = INITXDR( filename(1) , 'r' , .true.)
      write(6,*)' ixdrs(1) = ',ixdrs(1)
C
      IFAIL = IXDRINT( ixdrs(1) , DIM )
      IFAIL = IXDRINT( ixdrs(1) , NPOIN )
      IFAIL = IXDRINT( ixdrs(1) , NELEM )
      IFAIL = IXDRINT( ixdrs(1) , NBFAC )
      IFAIL = IXDRINT( ixdrs(1) , NHOLE )
C
      NOFVERT = DIM + 1
C
C     .. Reading info from the file with periodic bcs ..
C
      INQUIRE(FILE=FILENAME(4),EXIST=LFLAG)
      IF(LFLAG)THEN
      IXDRS(4) = INITXDR(FILENAME(4),'r',.TRUE.)
      IF( IXDRS(2).LT.0 )THEN
          WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(4)
          WRITE(ERRMSG(15:72),FMT="(A57)")FILENAME(4)(1:57)
          WRITE(6,*)ERRMSG
          CALL EXIT(1)
      ENDIF
      IFAIL = IXDRINT(IXDRS(4),NPNOD)
      write(6,*)'NPNOD = ',NPNOD
      ELSE
      write(6,*)'could not find file ',FILENAME(4)
          NPNOD=0
      ENDIF
C
C     .. Reading info from the solution file ..
C
      ixdrs(2) = INITXDR( filename(2) , 'r' , .true.)
C
      IFAIL = IXDRINT( ixdrs(2) , NITEMS )
      IF ((NPOIN+NPNOD).NE.NITEMS) THEN
         WRITE(6,*) 'INCONSISTENT NPOIN IN DATAFILE'
         WRITE(6,*) 'NITEMS = ',NITEMS
         WRITE(6,*) 'NPOIN = ',NPOIN
         WRITE(6,*) 'NPNOD = ',NPNOD
         CALL EXIT(1)
      ENDIF 
      IFAIL = IXDRINT( ixdrs(2) , NOFVAR )
      IFAIL = IXDRCLOSE( ixdrs(2) )
C
      RETURN
      END
