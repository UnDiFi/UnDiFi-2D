      SUBROUTINE WKSPCE
C
      IMPLICIT NONE
C
C     .. Parameters ..
C
C
C     .. Commons ..
C
      INCLUDE 'mesh_i4'
      INCLUDE 'int_flags'
      INCLUDE 'dim_flags'
      INCLUDE 'IO'
      INCLUDE 'es.h'
C
C     .. Scalar Arguments ..
C
C
C     .. Array Arguments ..
C
C
C     .. Local Scalars ..
C
      INTEGER NBFAC
      LOGICAL LFLAG
      EQUIVALENCE(NBFAC,nBoundaryFaces)
C
C     .. Local Arrays ..
C
      integer ifail,ixdrint
C
C     .. External Subroutines ..
C
C
C     .. External Functions ..
C
      INTEGER  INITXDR
      EXTERNAL INITXDR
C
C     .. Intrinsic Functions ..
C
C     .. Executable Statements ..
C
C     .. Reading info from the mesh file ..
C
      ixdrs(1) = INITXDR( filename(1) , 'r' , .FALSE.)
C
      IFAIL = IXDRINT( ixdrs(1) , DIM )
      IFAIL = IXDRINT( ixdrs(1) , NPOIN )
      IFAIL = IXDRINT( ixdrs(1) , NELEM )
      IFAIL = IXDRINT( ixdrs(1) , NBFAC )
      IFAIL = IXDRINT( ixdrs(1) , NHOLE )
C
      NOFVERT = DIM + 1
C
C     .. Solution file ..
C
      ixdrs(2) = INITXDR( filename(2) , 'r', .false. )
      IFAIL = IXDRINT( ixdrs(2) , NPOIN )
      IFAIL = IXDRINT( ixdrs(2) , NOFVAR )
C
      INQUIRE(FILE=filename(4),EXIST=LFLAG)
      IF(LFLAG)THEN
          ixdrs(4) = INITXDR( filename(4) , 'r' , .FALSE. )
          IFAIL = IXDRINT( ixdrs(4) , NPNOD )
!         WRITE(6,*)NPNOD,' periodic gridpoints'
      ELSE
          NPNOD = 0
      ENDIF
C
C
      RETURN
      END
