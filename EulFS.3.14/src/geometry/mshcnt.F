      SUBROUTINE MSHCNT(VCORG,NP,NGHOST,GHOSTS,NDIM,ICELNOD,ICELCEL,
     +                  NOFVERT,NELEM,IBNDFAC,NBFAC,NBINT,NHOLE,NFACE,
     +                  MAP,NPNOD,A,JA,IA,DEGREE,RCLINE,LDA,NCL,NZR,
     +                  IXDRS)
C
C     $Id: mshcnt.F,v 1.17 2020/09/30 08:32:35 abonfi Exp $
C
C     This routine reads the meshpoints, mesh connectivity
C     and boundary structure ..
C     in the sequential case NP = NPOIN+NPNOD
C     in the parallel   case NP = NPOIN+NPNOD+NGHOST
C
C
      IMPLICIT NONE
      INCLUDE'verbose.com'
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
      INCLUDE'io.com'
C
C     .. Scalar Arguments ..
      INTEGER IXDRS(*),NBFAC,NDIM,NELEM,NFACE,NBINT,NHOLE,NOFVERT,
     +        NP,NGHOST,NPNOD,NCL,LDA,NZR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VCORG(NDIM,NP)
      INTEGER IBNDFAC(3,NBFAC),ICELNOD(NOFVERT,NELEM),
     +ICELCEL(NOFVERT,NELEM),GHOSTS(NGHOST),MAP(NPNOD),
     +IA(NCL+1),JA(LDA,NZR),DEGREE(NCL)
      DOUBLE PRECISION A(LDA,NZR),RCLINE(NCL)
C     ..
C     .. Local Scalars ..
      INTEGER IFAIL,NITEMS
C     ..
C     .. External Subroutines ..
      INTEGER  IXDRCLOSE,IXDRIMAT
      EXTERNAL IXDRCLOSE,IXDRIMAT
      EXTERNAL R8Mat_Print,I4Mat_Print
#ifdef PARCH_t3d
      INTEGER  IXDRRMAT64
      EXTERNAL IXDRRMAT64
#else
      INTEGER  IXDRDMAT
      EXTERNAL IXDRDMAT
#endif
C     ..
C     .. Data statements ..
      DATA IFAIL/0/
C     ..
C
C     Reading nodal coordinates ..
C
      WRITE (NOUT,FMT=2000) MY_PE,NP,NGHOST,
     +NGHOST*100./REAL(NP),NDIM
C
C---------- COMPUTE THE NUMBER OF NORMALS ACCORDING TO Euler's FORMULA
C
C     NBINT is the number of interprocessor edges/faces
C
      IF (NDIM.EQ.2) THEN
          NFACE = 3*NP - (NBFAC+NBINT) - 3 + 3*NHOLE
      ELSE
          NFACE = 2*NELEM + (NBFAC+NBINT)/2
      ENDIF
C
#ifdef PARCH_t3d
      IFAIL = IXDRRMAT64(IXDRS(1),NDIM*NP,VCORG)
#else
      IFAIL = IXDRDMAT(IXDRS(1),NDIM*NP,VCORG)
#endif
      IF (IVERBOSE.GT.1) CALL R8Mat_Print('General',' ',NDIM,NP,VCORG,
     +                               NDIM,'Nodal coordinates',IFAIL)
C
C     Reading mesh connectivity
C
      WRITE (NOUT,FMT=3000) MY_PE,NELEM,NOFVERT
C
      IFAIL = IXDRIMAT(IXDRS(1),NOFVERT*NELEM,ICELNOD)
      IF (IVERBOSE.GT.1) CALL I4Mat_Print('General',' ',NOFVERT,NELEM,
     +                               ICELNOD,NOFVERT,
     +                               'Mesh connectivity',IFAIL)
C
C     Reading boundary data
C
      WRITE (NOUT,FMT=4000) MY_PE,NBFAC,NBINT,NHOLE
C
      IFAIL = IXDRIMAT(IXDRS(1),3*NBFAC,IBNDFAC)
      IF (IVERBOSE.GT.1) CALL I4Mat_Print('General',' ',3,NBFAC,
     +                        IBNDFAC,3,'Boundary info',IFAIL)
      IF(NGHOST.GT.0)IFAIL = IXDRIMAT(IXDRS(1),NGHOST,GHOSTS)
      IFAIL = IXDRCLOSE(IXDRS(1))
C
C     Reading cell neighbours data
C
      IFAIL = IXDRIMAT( IXDRS(2) , NOFVERT*NELEM , ICELCEL )
      IF (IVERBOSE.GT.1) CALL I4Mat_Print('General',' ',NOFVERT,NELEM,
     +                   ICELCEL,NOFVERT,'Neighbouring elements',IFAIL)
      IFAIL = IXDRCLOSE( IXDRS(2) )
C
C     Reading periodic nodes ..
C
      WRITE (NOUT,FMT=2100) MY_PE,NPNOD
C
C
C MAP -- Integer MAP(1:NPNOD)
C
C     maps the NPNOD periodic nodes onto the interior and ghost
C     nodes; 1<=MAP(1:NPNOD)<=NPOIN+NGHOST
C
      IF(NPNOD.NE.0)THEN
        IFAIL = IXDRIMAT(IXDRS(4),NPNOD,MAP)
        IF (IVERBOSE.GT.1)CALL I4Mat_Print('General',' ',1,NPNOD,
     +            MAP,1,'Periodic nodes',IFAIL)
        IFAIL = IXDRCLOSE(IXDRS(4))
      ENDIF
C
C     .. Reading c-lines linked list from the mesh file ..
C
      IF(NCL.NE.0)THEN
C
C     Reading c-lines ..
C
          WRITE (NOUT,FMT=2300) MY_PE,NCL,NZR
          IFAIL = IXDRDMAT( ixdrs(5) , NCL ,  RCLINE )
          IFAIL = IXDRIMAT( ixdrs(5) , NCL+1 , IA )
          IFAIL = IXDRIMAT( ixdrs(5) , NCL   , DEGREE )
          NITEMS = IA(NCL+1)-IA(1)
          IF( NZR .NE. NITEMS )STOP' aborted while reading c-lines '
          IFAIL = IXDRIMAT( ixdrs(5) , LDA*NZR , JA )
          IFAIL = IXDRDMAT( ixdrs(5) , LDA*NZR ,  A )
          IFAIL = IXDRCLOSE( ixdrs(5) )
          IF (IVERBOSE.GT.1)THEN
               CALL I4Mat_Print('General',' ',1,NCL,DEGREE,1,
     &         'DEGREE pointer ',IFAIL)
               CALL I4Mat_Print('General',' ',1,NCL+1,IA,1,
     &         'IA pointer ',IFAIL)
               CALL I4Mat_Print('General',' ',LDA,NZR,JA,LDA,
     &         'Meshpoints neighbours ',IFAIL)
               CALL R8Mat_Print('General',' ',LDA,NZR,A,LDA,
     &         'Area coord + theta ',IFAIL)
          ENDIF
      ENDIF
C
C
      RETURN
 2000 FORMAT (/,/,' INPUT OF NODES',10X,'ON PROCESSOR ',I3,/,1X,
     +       15 ('='),/,10X,
     +       'MAX. NUMBER OF NODES (NPOIN+NGHOST+NPNOD)=',I7,/,10X,
     +       'MAX. NUMBER OF GHOST NODES       (NGHOST)=',I6,/,10X,
     +       'RATIO NGHOST/(NPOIN+NGHOST+NPNOD)        =',F6.2,/,10X,
     +       'DIMENSIONS OF THE PROBLEM          (NDIM)=',I1,/)
 2100 FORMAT (/,/,' PERIODIC NODES',10X,'ON PROCESSOR ',I3,/,1X,
     +       15 ('='),/,10X,
     +       'NUMBER OF PERIODIC NODES   (NPOIN)=',I6,/)
 2300 FORMAT (/,/,' INPUT OF C-LINES ',10X,'ON PROCESSOR ',I3,/,1X,
     +       15 ('='),/,10X,
     +       'NUMBER OF C-LINES         (NCL)=',I3,/,10X,
     +       'NUMBER OF NODES ON ALL C-LINES =',I6,/)
 3000 FORMAT (/,/,' INPUT OF ELEMENTS',10X,'ON PROCESSOR',I3,/,1X,
     +       17 ('='),/,10X,
     +       'NUMBER OF ELEMENTS            (NELEM)=',I8,/,10X,
     +       'NUMBER OF NODES PER ELEMENT (NOFVERT)=',6X,I1,/)
 4000 FORMAT (/,/,' INPUT OF BOUNDARIES',10X,'ON PROCESSOR ',I3,/,1X,
     +       19 ('='),/,10X,
     +       'NUMBER OF BOUNDARY FACES     (NBFAC)=',I5,/,10X,
     +       'NUMBER OF INTER-PROCS. FACES (NBINT)=',I5,/,10X,
     +       'NUMBER OF HOLES              (NHOLE)=',I5,/)
      END
