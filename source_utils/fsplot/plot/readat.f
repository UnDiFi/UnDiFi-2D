      SUBROUTINE READAT(VCORG,ICELNOD,IBNDFAC,VZROE)
C
C     .. This routine reads the meshpoints, mesh connectivity
C        boundary structure AND solution ..
C
      IMPLICIT NONE
C
C     .. Parameters ..
C
      INCLUDE'es.h'
      INCLUDE'IO'
      INCLUDE'dim_flags'
      INCLUDE'int_flags'
      INCLUDE'mesh_i4'
C
C     .. Commons ..
C
C
C     .. Scalar Arguments ..
C
      integer ixdrclose,ixdrimat,ixdrdmat
      external ixdrclose,ixdrimat,ixdrdmat
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION VCORG(DIM,1),VZROE(NOFVAR,1)
      INTEGER ICELNOD(NOFVERT,1),IBNDFAC(3,1)
C
C     .. Local Scalars ..
C
      INTEGER NBFAC,IFAIL
      EQUIVALENCE(NBFAC,nBoundaryFaces)
C
      LOGICAL VERBOSE
C
      DATA IFAIL,VERBOSE /0,.FALSE./
C
C     .. Executable Statements ..
C
C     .. Reading nodal coordinates ..
C
      WRITE(NOUT,2000) NPOIN,NPNOD,NOFVAR,DIM,0
2000  FORMAT(//' INPUT OF NODES '/' ',15('=')/
     1  15X,'MAX. NUMBER OF NODES             (NPOIN)=',I7/
     1  15X,'MAX. NUMBER OF PERIODIC NODES    (NPNOD)=',I7/
     2  15X,'MAX. NUMBER OF D.O.F. PER NODE  (NOFVAR)=',I2/
     3  15X,'DIMENSIONS OF THE PROBLEM         (NDIM)=',I5/
     4  15X,'WORKSPACE IN REAL WORDS            (NVA)=',I12/)
C
      IFAIL = IXDRDMAT( ixdrs(1) , DIM*NPOIN , VCORG )
      IF(VERBOSE)
     &CALL X04CAF('General',' ',DIM,NPOIN,VCORG,DIM,
     +            'Nodal coordinates',IFAIL)
C
C     .. Reading mesh connectivity ..
C
      WRITE(NOUT,3000)NELEM,NOFVERT
3000  FORMAT(//' INPUT OF ELEMENTS '/' ',17('=')/
     1  15X,'MAX. NUMBER OF ELEMENTS            (NELEM)=',I8/
     2  15X,'MAX. NUMBER OF NODES PER ELEMENT (NOFVERT)=',6X,I1/)
C
      IFAIL = IXDRIMAT( ixdrs(1) , NOFVERT*NELEM , ICELNOD )
      IF(VERBOSE)
     &CALL X04EAF('General',' ',NOFVERT,NELEM,ICELNOD,NOFVERT,
     +            'Mesh connectivity',IFAIL)
C
C     .. Reading boundary data ..
C
      WRITE(NOUT,4000)NBFAC,NHOLE
4000  FORMAT(//' INPUT OF BOUNDARIES '/' ',19('=')/
     1  15X,'MAX. NUMBER OF BOUNDARY FACES   (NBFAC)=',I7/
     2  15X,'     NUMBER OF HOLES            (NHOLE)=',I4/)
C
      IFAIL = IXDRIMAT( ixdrs(1) , 3*NBFAC , IBNDFAC )
      IF(VERBOSE)
     &CALL X04EAF('General',' ',3,NBFAC,IBNDFAC,3,
     +            'Boundary info',IFAIL)
      IFAIL = IXDRCLOSE( ixdrs(1) )
C
C     .. Reading solution ..
C
      IFAIL = IXDRDMAT( ixdrs(2) , NOFVAR*NPOIN , VZROE )
      IFAIL = IXDRCLOSE( ixdrs(2) )
C
      RETURN
      END
