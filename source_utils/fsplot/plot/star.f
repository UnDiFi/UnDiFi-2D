      SUBROUTINE STAR(FILENAME,POINT,ICELNOD,Z,IBNDPTR,NBFAC)
C
      IMPLICIT NONE
C
C This routine writes a tecplot formatted file
C
      INCLUDE 'constants'
      INCLUDE 'mesh_i4'
      INCLUDE 'int_flags'
      INCLUDE 'dim_flags'
      INCLUDE 'IO'
C
C     .. Scalar Arguments ..
C
      CHARACTER*(*) FILENAME
C
C     .. Array Arguments ..
C
      INTEGER NBFAC
      DOUBLE PRECISION  POINT,Z
      INTEGER   ICELNOD,IBNDPTR(3,NBFAC)
      DIMENSION POINT(DIM,1),Z(NOFVAR,1),ICELNOD(NOFVERT,1)
C
C     .. Local Scalars ..
C
      INTEGER IVERT,IBC,K,I
      CHARACTER*11  Elementtype
      CHARACTER*256  VarString
      CHARACTER*10   STRNG(3)
      INTEGER*4 Debug,Vlength,IPOIN,IELEM,IFREQ
      INTEGER*4 J,LastChar,IVAR
      DATA STRNG/"wing      ","symmetry  ","far-field "/
C
C     .. Local Arrays ..
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
      INTEGER ICYCL
C
C     .. Intrinsic Functions ..
C
      INTRINSIC INDEX,MAX0
C
C     .. Executable Statements ..
C
C...Set defaults
C
      OPEN(40,FILE='sinus.vrt')
      OPEN(42,FILE='sinus.cel')
      OPEN(44,FILE='sinus.bnd')
C
      LastChar = INDEX(filename,CHAR(32))
C
C
C       writing the nodes
C
      IFREQ = MAX0( NPOIN/20 , 1 )
      DO 1 IPOIN = 1,NPOIN
        WRITE(40,FMT=125)IPOIN,(POINT(J,IPOIN),J=1,DIM)
    1 CONTINUE
      CLOSE(40)
C
C       writing the elements
C
      DO 2 IELEM = 1 , NELEM
        WRITE(42,FMT=135)IELEM,
     &  ICELNOD(1,IELEM),ICELNOD(2,IELEM),ICELNOD(3,IELEM),
     &  ICELNOD(3,IELEM),ICELNOD(4,IELEM),ICELNOD(4,IELEM),
     &  ICELNOD(4,IELEM),ICELNOD(4,IELEM),1,1
    2 CONTINUE
      CLOSE(42)
C
C
C       writing the elements
C
      DO 7 I= 1 , NBFAC
        IVERT=IBNDPTR(1,I)
        IELEM=IBNDPTR(2,I)
        IBC=IBNDPTR(3,I)
        WRITE(44,FMT=145)I,
     & (ICELNOD(ICYCL(IVERT+K,NOFVERT),IELEM),K=1,3),
     &  ICELNOD(ICYCL(IVERT+3,NOFVERT),IELEM),IBC,0,STRNG(IBC)
    7 CONTINUE
      CLOSE(44)
C
C
      RETURN
C
  111 FORMAT('.',$)
* 111 FORMAT('x')
  125 FORMAT(I9,6X,3G16.9)
  135 FORMAT(I9,6X,8I9,I9,2I5)
  140 FORMAT(5X,'IFMT = ',I2,' MUST be 1 or 2')
  145 FORMAT(I8,6X,4I9,2I7,6X,A10)
  155 FORMAT(5X,'I am writing the Tecplot file ... ',/)
  165 FORMAT(5X,'Writing nodes ',$)
  170 FORMAT(5X,'Writing elements ',$)
  175 FORMAT(' done !',/)
  180 FORMAT(5X,'Tecplot file WRITTEN',/)
  185 FORMAT(/,5X,'Tecplot file WRITTEN to ... ',A60/)
  999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END
