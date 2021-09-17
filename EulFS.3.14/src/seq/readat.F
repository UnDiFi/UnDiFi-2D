      SUBROUTINE READAT(NELEM,NPOIN,NGHOST,NPNOD,
     +                  NBPOIN,NFACE,NBFAC,NBINT,NHOLE,NOFVERT,
     +                  NDIM,NOFVAR,NTURB,
     &                  FNAME1,FNAME2,FNAME3,FNAME4)
C
C     This routine reads the mesh parameters from the meshfiles
C     (FNAME1) and allocates the space required for reading
C     the data structure
C
C     $Id: readat.F,v 1.32 2013/06/04 15:37:42 abonfi Exp $
C
C
C     for the UNI-processor case
C         NPOIN is the # of meshpoints 
C         NGHOST is 0
C
C     for the parallel case
C         NPOIN is the # of processor owned meshpoints 
C         NGHOST is the number of off-processor nodes or
C         inter-processor nodes
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'datatype.h'
      INCLUDE 'plasma.h'
      INCLUDE 'bnd.h'
      INCLUDE 'time.h'
      INCLUDE 'io.com'
      INCLUDE 'dofs.com'
      INCLUDE 'periodic.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'nloc.com'
      INCLUDE 'flags.com'
      INCLUDE 'time.com'
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NBPOIN,NDIM,NELEM,NFACE,NHOLE,NOFVAR,NOFVERT,NPOIN,
     +NTURB,NGHOST,NPNOD
      CHARACTER FNAME1* (*),FNAME2* (*), FNAME3* (*),
     &          FNAME4* (*)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER MY_PE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER IFAIL,IXDRS(5),NBINT,NERR,IOPT,NITEMS,NZR
      CHARACTER*72 ERRMSG
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER        INITXDR,ISTKGT,IXDRINT,I1MACH
      EXTERNAL DNRM2,INITXDR,ISTKGT,IXDRINT,I1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT,MSHCNT,NODCOD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON/MPICOM/MY_PE
C
      DATA ERRMSG(1:14) /'IXDRS=00;FILE='/
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
      INTEGER LGHOST
      EQUIVALENCE (LTZX,LGHOST)
C     ..
      WRITE (NOUT,FMT=4000)
C
C     .. Reading info from mesh file 001 ..
C
      IXDRS(1) = INITXDR(FNAME1,'r',.TRUE.)
      IF( IXDRS(1).LT.0 )THEN
          WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(1)
          WRITE(ERRMSG(15:72),FMT="(A57)")FNAME1(1:57)
          NERR = 3
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
      IFAIL = IXDRINT(IXDRS(1),NDIM)
      IFAIL = IXDRINT(IXDRS(1),NPOIN)
      IFAIL = IXDRINT(IXDRS(1),NELEM)
      IFAIL = IXDRINT(IXDRS(1),NBFAC)
#ifdef MPI
      IFAIL = IXDRint(ixdrs(1),NBINT)
#else
      NBINT = 0
#endif 
      IFAIL = IXDRINT(IXDRS(1),NHOLE)
#ifdef MPI
      IFAIL = IXDRINT(IXDRS(1),NGHOST)
#else
      NGHOST = 0
#endif 
C
C     .. Reading info from mesh file 002 ..
C
      IXDRS(2) = INITXDR( FNAME2 , 'r' , .FALSE. ) 
C
      IFAIL = IXDRINT( IXDRS(2) , NITEMS )
      IF( NITEMS .NE. NDIM )THEN
          WRITE(NOUT,*) 'Invalid no of vertices in file002',NITEMS 
          STOP
      ENDIF 
      IFAIL = IXDRINT( IXDRS(2) , NITEMS )
      IF( NITEMS .NE. NELEM )THEN
          STOP 'Invalid no of cells in file002' 
      ENDIF 
C
C     .. Reading info from mesh file 004 (periodic bcs) ...
C
      INQUIRE(FILE=FNAME3,EXIST=LFLAG)
      IF(LFLAG)THEN
          PERIODIC_MESH = .TRUE.
          IXDRS(4) = INITXDR(FNAME3,'r',.TRUE.)
          IF( IXDRS(4).LT.0 )THEN
              WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(4)
              WRITE(ERRMSG(15:72),FMT="(A57)")FNAME3(1:57)
              NERR = 3
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
          IFAIL = IXDRINT(IXDRS(4),NPNOD)
      ELSE
          NPNOD=0
          PERIODIC_MESH = .FALSE.
      ENDIF
C
C     .. Reading info from mesh file 006 (c-lines) ...
C
      INQUIRE(FILE=FNAME4,EXIST=LCLTHERE)
      IF(LCLTHERE)THEN
          IXDRS(5) = INITXDR(FNAME4,'r',.TRUE.)
          IF( IXDRS(5).LT.0 )THEN
              WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(5)
              WRITE(ERRMSG(15:72),FMT="(A57)")FNAME4(1:57)
              NERR = 3
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
          IFAIL = IXDRINT(IXDRS(5),NCL) ! nof c-lines
          IFAIL = IXDRINT(IXDRS(5),NZR) ! nof vertices on all c-lines
      ELSE
          NCL=0
          NZR=0
      ENDIF
C
C     Only triangles/tetrahedra are allowed ..
C
      NOFVERT = NDIM + 1
C
C     KAN stands for Kind of ANalysis
C
      IF (ABS(KAN).EQ.1 .OR. KAN.EQ.5) THEN ! scalar or potential (Laplace's eqn)
          NTURB = 0
          NOFVAR = 1
      ENDIF
      IF (ABS(KAN).EQ.2) THEN ! Incompressible
          NOFVAR = NDIM + 1
          IE = -23414
          IX = 2
          IY = 3
          IZ = 4
      ENDIF
      IF (ABS(KAN).EQ.4) THEN ! Compressible
          NOFVAR = NDIM + 2
          IE = 2
          IX = 3
          IY = 4
          IZ = 5
      ENDIF
      IF (ABS(KAN).EQ.3) THEN ! Argon
          NOFVAR = NDIM + 5
          IE = NSP+1
          IX = NSP+2
          IY = NSP+3
          IZ = NSP+4
          IF(LAPLACE)nofvar = nofvar + 1
      ENDIF
cxxxx
      IF(TURBULENT .AND. COUPLED) then
         NOFVAR = NOFVAR + 1
C
C     It is important to set NTURB=0 or subr. backup()
C     will try to save the array DSTAK(LTURB)
C
         NTURB = 0
      endif
cxxxx
cxxxx
caldo        NOFVAR = NOFVAR + 1
cxxxx
cxxxx
C
C
C ---------- ALLOCATE SPACE
C
C     Nodal coordinates ..
C
      NITEMS = (NPOIN+NGHOST+NPNOD)*NDIM 
      write(6,*)nitems
      IF(LALE) NITEMS = NITEMS*MAXTIMLEVS ! store the coordinates of all time levels
      write(6,*)nitems,lale
      LCORG = ISTKGT(NITEMS,KIND_REAL8)
      CALL DINIT(NITEMS,ZERO,DSTAK(LCORG),1)
C
C     Nodal grid velocities ..
C
      IF(LALE)THEN
         NITEMS = (NPOIN+NGHOST+NPNOD)*NDIM 
         LXYZDOT = ISTKGT(NITEMS,KIND_REAL8)
         CALL DINIT(NITEMS,ZERO,DSTAK(LXYZDOT),1)
      ENDIF
C
C     Nodal values ..
C
      NITEMS = (NPOIN+NGHOST+NPNOD)*NOFVAR 
      IF(LTIME) NITEMS = NITEMS*MAXTIMLEVS
C
      LZROE = ISTKGT(NITEMS,KIND_REAL8)
      CALL DINIT(NITEMS,ZERO,DSTAK(LZROE),1)
C
C     Turbulent stuff ..
C
      IF( TURBULENT )THEN
C     Turbulent variables ..
          IF(SEGREGATED)THEN
             NITEMS = (NPOIN+NGHOST+NPNOD)*NTURB 
             IF(LTIME) NITEMS = NITEMS*MAXTIMLEVS
             LTURB = ISTKGT(NITEMS,KIND_REAL8)
             CALL DINIT(NITEMS,ZERO,DSTAK(LTURB),1)
          ENDIF
C     Nearest wall distance ..
          LTD = ISTKGT(NPOIN+1+NGHOST+NPNOD,KIND_REAL8)
          LTTD = ISTKGT(NPOIN+1+NGHOST+NPNOD,KIND_REAL8)
          CALL DINIT((NPOIN+1+NGHOST),ZERO,DSTAK(LTD),1)
          CALL DINIT((NPOIN+1+NGHOST),ZERO,DSTAK(LTTD),1)
      ENDIF
C
C     Cell to node pointer ..
C
      NITEMS = 2*(NELEM*NOFVERT)
      LCELNOD = ISTKGT(NITEMS,KIND_INTEGER)
      CALL IINIT(NITEMS,0,ISTAK(LCELNOD),1)
C
C     Cell to cell pointer ..
C
      LCELCEL = ISTKGT(NELEM*NOFVERT,KIND_INTEGER)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELCEL),1)
C
C     Boundary structure ..
C
      LBNDFAC = ISTKGT(3*NBFAC,KIND_INTEGER)
      CALL IINIT(3*NBFAC,0,ISTAK(LBNDFAC),1)
C
C     Nodal code ..
C
      LNODCOD = ISTKGT(NPOIN+NGHOST+NPNOD,KIND_INTEGER)
      CALL IINIT(NPOIN+NGHOST+NPNOD,0,ISTAK(LNODCOD),1)
C
C     ghost code ..
C
#ifdef MPI
      LGHOST = ISTKGT(NGHOST,KIND_INTEGER)
      CALL IINIT(NGHOST,0,ISTAK(LGHOST),1)
#endif
C
C     periodic mapping ..
C
      IF(NPNOD.GT.0)THEN
          LPMAP = ISTKGT(NPNOD,KIND_INTEGER)
          CALL IINIT(NPNOD,0,ISTAK(LPMAP),1)
      ENDIF
C
C     c-lines linked list ..
C
      IF(NCL.NE.0)THEN
          LCLRC = ISTKGT(    NCL,KIND_REAL8)
          LCLDEG= ISTKGT(    NCL,KIND_INTEGER)
          LCLA  = ISTKGT(LCLDA*NZR,KIND_REAL8)
          LCLJA = ISTKGT(LCLDA*NZR,KIND_INTEGER)
          LCLIA = ISTKGT(  NCL+1,KIND_INTEGER)
          LCLZB = ISTKGT(NOFVAR*NCL,KIND_REAL8)
          CALL DINIT(    NCL,ZERO,DSTAK(LCLRC),1) ! radial position
          CALL DINIT(LCLDA*NZR,ZERO,DSTAK(LCLA ),1) ! a(LCLDA,*) area coords
          CALL IINIT(LCLDA*NZR,0,ISTAK(LCLJA),1)    ! ja(LCLDA,*)
          CALL IINIT(  NCL+1,0,ISTAK(LCLIA),1)    ! ia(1:NCL+1)
          CALL IINIT(  NCL  ,0,ISTAK(LCLDEG),1)   ! degree(1:NCL)
          CALL DINIT(NOFVAR*NCL,ZERO,DSTAK(LCLZB),1) ! zb(NOFVAR,*)
      ENDIF
C
C     Reads mesh points, connectivity and boundary structure ..
C     REM: routine MSHCNT() is called with the number of 
C          meshpoints being equal to NPOIN+NGHOST+NPNOD
C
C
      CALL MSHCNT(DSTAK(LCORG),NPOIN+NGHOST+NPNOD,NGHOST,ISTAK(LGHOST),
     +            NDIM,ISTAK(LCELNOD),ISTAK(LCELCEL),NOFVERT,NELEM,
     +            ISTAK(LBNDFAC),NBFAC,NBINT,NHOLE,NFACE,ISTAK(LPMAP),
     +            NPNOD,DSTAK(LCLA),ISTAK(LCLJA),ISTAK(LCLIA),
     &            ISTAK(LCLDEG),DSTAK(LCLRC),LCLDA,NCL,NZR,IXDRS)
C
C     store a copy of the cell to vertex pointer; this
C     is used for periodic grids
C
      CALL ICOPY(NELEM*NOFVERT,ISTAK(LCELNOD),1,
     &           ISTAK(LCELNOD+NELEM*NOFVERT),1)
C
C     Sets the node code: a flag for labeling boundary nodes ..
C     Does it REALLY make sense to call NODCOD() here ?!?! 
C
      CALL NODCOD(ISTAK(LNODCOD),NPOIN+NGHOST+NPNOD,NBPOIN,
     +            ISTAK(LCELNOD),NOFVERT,NELEM,ISTAK(LBNDFAC),NBFAC)
C
      RETURN
C
 2020 FORMAT (15X,'FUNCTION ',A31,/,
     +5X,'HAS RETURNED A NONZERO ERROR FLAG =',I6,/,
     +5X,'IN SUBROUTINE ',A10,' ON PE =',I4)
 2030 FORMAT (5X,'WHILE TRYING TO OPEN FILE:',/,A75)
 4000 FORMAT (/,/,' MEMORY ALLOCATION   ',/,' ',19 ('='),/)
      END
