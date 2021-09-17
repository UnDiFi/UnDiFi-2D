      SUBROUTINE EXGEO(NELEM,NOFVAR,NPOIN,NGHOST,NPNOD,NP,NBPOIN,NFACE,
     +                 NBFAC,NBINT,NWFAC,NBODY4,NBODY6,NHOLE,NOFVERT,
     3                 NDIM)
C
C     $Id: exgeo.F,v 1.19 2020/03/25 11:05:33 abonfi Exp $
C
C     This routine computes geometrical entities and 
C     performs some checks
C     in the sequential case NP = NPOIN+NPNOD (NGHOST=0)
C     in the parallel   case NP = NPOIN+NPNOD+NGHOST
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'time.h'
      INCLUDE 'datatype.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'flags.com'
      INCLUDE 'time.com'
      INCLUDE 'io.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'nloc.com'
C
C     .. Scalar Arguments ..
      INTEGER NOFVAR,NBFAC,NBODY4,NBODY6,NBPOIN,NDIM,NELEM,NFACE,NHOLE,
     +        NOFVERT,NP,NWFAC,NBINT,NPOIN,NGHOST,NPNOD
C
C     On entry:
C     --------
C     NELEM    no. of processor owned elements (triangles/tetrahedra);
C              global number of elements in the uni-processor case
C     NBFAC    no. of processor owned boundary faces/edges;
C              global number of boundary faces/edges 
C              in the uni-processor case
C
C     Upon return:
C     -----------
C
C
C
C
C
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER MY_PE 
C     ..
C     .. Local Scalars ..
      INTEGER I,IFAIL,IOFF,LINDX,LWORK,LWKSP,NITEMS
      LOGICAL LFLAG
      CHARACTER*12 FNAME
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1)
C     ..
C     .. External Functions ..
      INTEGER ISTKGT,ISTKQU
      EXTERNAL ISTKGT,ISTKQU
C     ..
C     .. External Subroutines ..
      EXTERNAL CHKBND,CMPVOL,DINIT,FF,ICOPY,IINIT,ISORTRX,ISTKRL,REORDR
C     ..
C     .. Common blocks ..
      COMMON /MPICOM/MY_PE
      COMMON /CSTAK/DSTAK
C
      LOGICAL VERBOSE
C     ..
C     .. Equivalences ..
      EQUIVALENCE (ISTAK(1),DSTAK(1))
C     ..
C
C     set VERBOSE = .TRUE. if you want to dump the boundary nodes and colours
C
      VERBOSE = .FALSE.
C
C---------- PRINT BLOCK PARAMETERS
C
      WRITE (NOUT,FMT=2000) MY_PE
      WRITE (NOUT,FMT=2005) NFACE
C
C     allocate an array to store the edges/faces of the mesh 
C
      LFACNOR = ISTKGT(NDIM* (NFACE+1),KIND_REAL8)
      CALL DINIT(NDIM* (NFACE+1),ZERO,DSTAK(LFACNOR),1)
C
C     allocate cell to face pointers
C
      LCELFAC = ISTKGT(NOFVERT*NELEM,KIND_INTEGER)
      CALL IINIT(NOFVERT*NELEM,0,ISTAK(LCELFAC),1)
C
C     Finds all faces of the mesh
C
      CALL FF(ISTAK(LCELNOD),ISTAK(LCELCEL),ISTAK(LCELFAC),NOFVERT,
     +        NELEM,DSTAK(LCORG),NDIM,NP,DSTAK(LFACNOR),NFACE,NBFAC,
     +        NBINT)
C
C     Allocate space to store the volumes ..
C
      NITEMS = NELEM
      IF(LALE)NITEMS=NITEMS*(MAXTIMLEVS+1) ! leaves room for time levels: n+a,n+1,n,n-1
      LVOL = ISTKGT(NITEMS,KIND_REAL8)
      CALL DINIT(NITEMS,ZERO,DSTAK(LVOL),1)
C
C     The volume of EACH element is computed
C
      CALL CMPVOL(ISTAK(LCELNOD),ISTAK(LCELFAC),NOFVERT,NELEM,
     +            DSTAK(LFACNOR),NFACE,DSTAK(LCORG),NDIM,NP,
     +            DSTAK(LVOL),ZERO)
C
      IF( LTIME )THEN ! if time-accurate
C
C     Allocate space to store the median dual cell volumes ..
C
         NITEMS = NP
         IF(LALE)NITEMS = NITEMS*MAXTIMLEVS
         LMEDIAN = ISTKGT(NITEMS,KIND_REAL8)
         CALL DINIT(NITEMS,ZERO,DSTAK(LMEDIAN),1)
C
         CALL MEDIAN(DSTAK(LMEDIAN),NPOIN,NGHOST,NPNOD,NP,DSTAK(LVOL),
     &   ISTAK(LCELNOD),ISTAK(LPMAP),NOFVERT,NELEM)
      ENDIF
C
C     Sort boundary info for ascending color
C
      LWORK = ISTKGT(NBFAC,KIND_INTEGER)
      LINDX = ISTKGT(NBFAC,KIND_INTEGER)
C
      IOFF = 2
C     copy colours into the work-array
      CALL ICOPY(NBFAC,ISTAK(LBNDFAC+IOFF),3,ISTAK(LWORK),1)
      CALL ISORTRX(NBFAC,ISTAK(LWORK),ISTAK(LINDX))
      CALL REORDR(NBFAC,ISTAK(LBNDFAC+IOFF),3,ISTAK(LINDX),ISTAK(LWORK))
C
C     now sort IBNDFAC(1:2,*) according to the permutation
C
      DO 10 IOFF = 0,1
          CALL ICOPY(NBFAC,ISTAK(LBNDFAC+IOFF),3,ISTAK(LWORK),1)
          CALL REORDR(NBFAC,ISTAK(LBNDFAC+IOFF),3,ISTAK(LINDX),
     +                ISTAK(LWORK))
   10 CONTINUE
      IF(NCL.NE.0)THEN
C
C     RNKIDX INVERTS A PERMUTATION, AND HENCE CONVERTS A RANK VECTOR TO
C     AN INDEX VECTOR, OR VICE VERSA.
C
cnag      CALL M01ZAF(ISTAK(LINDX),1,NBFAC,IFAIL)
          CALL RNKIDX(ISTAK(LINDX),1,NBFAC,IFAIL)
          IF (IFAIL.NE.0) THEN
              WRITE (6,FMT=*) 'M01ZAF RNKIDX within exgeo '
              CALL EXIT(IFAIL)
          ENDIF
          CALL RANKCL(ISTAK(LINDX),DSTAK(LCLA),ISTAK(LCLJA),
     &            ISTAK(LCLIA),LCLDA,NCL)
      ENDIF
C
C     Free INDX and IWORK
C
      CALL ISTKRL(2)
C
C     perform some checks on the boundary data structure
C
      IF(VERBOSE)THEN
         LWORK = ISTKGT(NOFVERT*NBFAC,KIND_INTEGER)
         CALL IINIT(NOFVERT*NBFAC,0,ISTAK(LWORK),1)
      ELSE
         LWORK = 1
      ENDIF
      CALL CHKBND(ISTAK(LBNDFAC),NBFAC,ISTAK(LCELFAC),ISTAK(LCELNOD),
     +            NOFVERT,NELEM,DSTAK(LFACNOR),DSTAK(LCORG),NDIM,NFACE,
     +            NWFAC,NBODY4,NBODY6,DSTAK(LCLA),ISTAK(LCLJA),
     &            ISTAK(LCLIA),LCLDA,NCL,VERBOSE,ISTAK(LWORK))
      IF(VERBOSE)CALL ISTKRL(1)
C
C     allocate space for storing the skin friction coefficient
C
      IF (NWFAC.NE.0) THEN
          LSKINF = ISTKGT(NWFAC,KIND_REAL8)
          LHEAT = ISTKGT(NWFAC,KIND_REAL8)
      ENDIF
C
C     prescribed flux boundary condition
C
      IF(.NOT.LREAD(2))GOTO 120
      IF(MY_PE.NE.0)THEN
         WRITE(6,*)'Reading bflux files un-implemented in MPI'
         CALL EXIT(128)
      ENDIF
      FNAME = "bflux000.dat"
      WRITE(FNAME(6:8),FMT="(I3.3)")MY_PE
      OPEN(120,FILE=FNAME,STATUS="OLD")
      READ(120,*)NBFLX(2) ! nof bndry faces with prescribed flux
      LBNDFLX(1) = ISTKGT(NOFVAR*NBFLX(2),KIND_REAL8) ! store fluxes
      LBNDFLX(2) = ISTKGT(NBFLX(2),KIND_INTEGER) ! store pointer
      LWORK = ISTKGT((NOFVERT-1)*NBFLX(2),KIND_INTEGER) ! temporary pointer storage
      LWKSP = ISTKGT(NBFLX(2),KIND_REAL8) ! temporary storage for sorting
      CALL SETBNDFLXPTR(DSTAK(LBNDFLX(1)),ISTAK(LWORK),
     2                  ISTAK(LBNDFLX(2)),NOFVAR,NDIM,NOFVERT,
     3                  ISTAK(LBNDFAC),NBFAC,ISTAK(LCELNOD),
     4                  DSTAK(LWKSP),120,1)
      CALL SETBNDFLXPTR(DSTAK(LBNDFLX(1)),ISTAK(LWORK),
     2                  ISTAK(LBNDFLX(2)),NOFVAR,NDIM,NOFVERT,
     3                  ISTAK(LBNDFAC),NBFAC,ISTAK(LCELNOD),
     4                  DSTAK(LWKSP),120,2)
      
      CALL ISTKRL(2) ! release both work arrays
  120 CONTINUE
C
      RETURN
C
C     I/O FORMATS
C
 2000 FORMAT (/,/,' COMPUTATION OF THE NORMALS ',10X,
     +'ON PROCESSOR ',I3,/,' ',26 ('='))
 2005 FORMAT (10X,'ESTIMATED NUMBER OF FACES(EDGES) : ',I8,/)

      END
