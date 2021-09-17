      SUBROUTINE READAT(VCORG,ICELNOD,ICELCEL,IBNDFAC,ZROE,NDIM,
     +NOFVAR,NOFVERT,NPOIN,NELEM,NBFAC)
C
C     This routine reads the meshpoints, mesh connectivity
C        boundary structure AND solution ..
C
      IMPLICIT NONE
C
      INTEGER*4 NIN,NOUT
      PARAMETER (NIN=5,NOUT=6)
C
C	NIN	is the OUTPUT device number
C	NOUT	is the INPUT device number
C
C
      INTEGER ixdrs(10)
      CHARACTER*80 filename(5)
      COMMON/ES/ixdrs,filename
C
      INTEGER NDIM,NOFVAR,NOFVERT,NPOIN,NELEM,NHOLE,IDUMMY
C
      DOUBLE PRECISION VCORG(NDIM,1),ZROE(NOFVAR,1)
      INTEGER ICELNOD(NOFVERT,1),ICELCEL(NOFVERT,NELEM),IBNDFAC(3,1)
C
      INTEGER INITXDR,IXDRINT,IXDRIMAT,IXDRCLOSE,IXDRDMAT
C
      INTEGER NBFAC,IFAIL
C
C
      DATA IFAIL /0/
C
C     Reading nodal coordinates ..
C
      WRITE(NOUT,2000) NPOIN,NOFVAR,NDIM,0
2000  FORMAT(//' INPUT OF NODES '/' ',15('=')/
     1  15X,'MAX. NUMBER OF NODES             (NPOIN)=',I5/
     2  15X,'MAX. NUMBER OF D.O.F. PER NODE  (NOFVAR)=',I5/
     3  15X,'DIMENSIONS OF THE PROBLEM         (NDIM)=',I5/
     4  15X,'WORKSPACE IN REAL WORDS            (NVA)=',I10/)
C
      IFAIL = IXDRDMAT( ixdrs(1) , NDIM*NPOIN , VCORG )
*     CALL X04CAF('General',' ',NDIM,NPOIN,VCORG,NDIM,
*    +            'Nodal coordinates',IFAIL)
C
C     Reading mesh connectivity ..
C
      WRITE(NOUT,3000)NELEM,NOFVERT
3000  FORMAT(//' INPUT OF ELEMENTS '/' ',17('=')/
     1  15X,'MAX. NUMBER OF ELEMENTS            (NELEM)=',I7/
     2  15X,'MAX. NUMBER OF NODES PER ELEMENT (NOFVERT)=',6X,I1/)
C
      IFAIL = IXDRIMAT( ixdrs(1) , NOFVERT*NELEM , ICELNOD )
*     CALL X04EAF('General',' ',NOFVERT,NELEM,ICELNOD,NOFVERT,
*    +            'Mesh connectivity',IFAIL)
C
C     Reading boundary data ..
C
      WRITE(NOUT,4000)NBFAC,NHOLE
4000  FORMAT(//' INPUT OF BOUNDARIES '/' ',19('=')/
     1  15X,'MAX. NUMBER OF BOUNDARY FACES   (NBFAC)=',I5/
     2  15X,'     NUMBER OF HOLES            (NHOLE)=',I5/)
C
      IFAIL = IXDRIMAT( ixdrs(1) , 3*NBFAC , IBNDFAC )
caldo CALL X04EAF('General',' ',3,NBFAC,IBNDFAC,3,
caldo+            'Boundary info',IFAIL)
      IFAIL = IXDRCLOSE( ixdrs(1) )
C
C     read neighbours from file FILE
C
caldo IXDRS(3) = INITXDR(filename(3),'r',.false.)
caldo IFAIL = IXDRINT(IXDRS(3),IDUMMY)
caldo IFAIL = IXDRINT(IXDRS(3),IDUMMY)
caldo IFAIL = IXDRIMAT(IXDRS(3),NOFVERT*NELEM,ICELCEL)
C
      CALL SOLZNE(FILENAME(2),ZROE,NOFVAR,NPOIN,"r")
C
      RETURN
      END
C
