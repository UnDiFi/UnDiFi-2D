C
C ----------------------------------- + -----------------------------------
C
      SUBROUTINE CELNHB(ICELCEL,NOFVERT,NELEM,FILENAME)
C
C     .. This routine reads che cell neighbours ..
C
      IMPLICIT NONE
C
C     .. Parameters ..
C
C
C     .. Commons ..
C
C
C     .. Scalar Arguments ..
C
      INTEGER NOFVERT,NELEM 
      CHARACTER*(*) FILENAME 
C
C     .. Array Arguments ..
C
      INTEGER ICELCEL(NOFVERT,NELEM)
C
C ICELCEL -- Integer ICELCEL(1:NOFVERT,1:NELEM)
C            Cell to Cell pointer : ICELCEL(i,ielem) gives the
C            global number of the element sharing with ielem
C            the face opposite the i-th vertex of the ielem-th cell
C            If ICELCEL(i,ielem) = 0 or ICELCEL(i,ielem) > NELEM
C            the element ielem is a boundary element and the face
C            opposite its i-th vertex is a boundary face
C
C
C     .. Local Scalars ..
C
      INTEGER IXDRS,IFAIL,IDUMMY
C
C     .. Local Arrays ..
C
C
C     .. External Subroutines ..
C
C
C     .. External Functions ..
C
      INTEGER  INITXDR,IXDRINT,IXDRIMAT,IXDRCLOSE
      EXTERNAL INITXDR,IXDRINT,IXDRIMAT,IXDRCLOSE
C
C     .. Intrinsic Functions ..
C
C
      DATA IFAIL /0/
C
C     .. Executable Statements ..
C
      IXDRS = INITXDR( filename , 'r' , .FALSE. ) 
C
      IFAIL = IXDRINT( IXDRS , IDUMMY )
      IFAIL = IXDRINT( IXDRS , IDUMMY )
      IFAIL = IXDRIMAT( IXDRS , NOFVERT*NELEM , ICELCEL )
C     CALL I4Mat_Print('General',' ',NOFVERT,NELEM,ICELCEL,NOFVERT,
C    +            'Neighbouring elements',IFAIL)
      IFAIL = IXDRCLOSE( IXDRS )
C
      RETURN
      END
