      SUBROUTINE GRIDPTS(NPOIN,NGHOST,NPNOD,NELEM,NBFAC,NBINT,NHOLE,
     +                   NOFVERT,NDIM,VCORG,ICELNOD,IBNDPTR,
     3                   FNAME1,FNAME3,JOB)
C
C     This routine reads or writes the mesh parameters from the meshfiles
C     (FNAME1)
C
C     $Id: rgrdpts.F,v 1.4 2020/03/28 09:42:41 abonfi Exp $
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
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'datatype.h'
      INCLUDE 'io.com'
      INCLUDE 'periodic.com'
      INCLUDE 'ibc2.com'
      INCLUDE 'flags.com'
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NHOLE,NOFVERT,NPOIN,
     +NGHOST,NPNOD
      DOUBLE PRECISION VCORG(NDIM,*)
      INTEGER ICELNOD(NOFVERT,NELEM),IBNDPTR(3,NBFAC)
      CHARACTER FNAME1* (*), FNAME3* (*), JOB*1
C     ..
C     .. Arrays in Common ..
!     DOUBLE PRECISION DSTAK(1)
      INTEGER MY_PE
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER IFAIL,IXDRS(5),NBINT,NERR,IOPT,NITEMS,IVERBOSE,NP
      CHARACTER*72 ERRMSG
      LOGICAL LFLAG
C     ..
C     .. Local Arrays ..
!     INTEGER ISTAK(1)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER        INITXDR,ISTKGT,IXDRINT,IXDRDMAT,IXDRIMAT,IXDRCLOSE
      EXTERNAL DNRM2,INITXDR,ISTKGT,IXDRINT,IXDRDMAT,IXDRIMAT,IXDRCLOSE
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
!     COMMON /CSTAK/DSTAK
      COMMON/MPICOM/MY_PE
C
      DATA ERRMSG(1:14) /'IXDRS=00;FILE='/
C     ..
C     .. Equivalences ..
C     ..
      IF( .NOT. ((JOB.EQ.'r').OR.(JOB.EQ.'w')) )THEN
          STOP "JOB should be either 'r' or 'w'"
      ENDIF
      IF(JOB.EQ.'r')THEN
         WRITE (NOUT,FMT=400)FNAME1
      ELSE
         WRITE (NOUT,FMT=405)FNAME1
      ENDIF
      IVERBOSE = 0
C
      IXDRS(1) = INITXDR(FNAME1,JOB,.TRUE.)
      IF( IXDRS(1).LT.0 )THEN
          WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(1)
          WRITE(ERRMSG(15:72),FMT="(A57)")FNAME1(1:57)
          NERR = 3
          IOPT = 1
          CALL SETERR(ERRMSG,72,NERR,IOPT)
      ENDIF
C
      IF(JOB.EQ.'w')NITEMS=NDIM
      IFAIL = IXDRINT(IXDRS(1),NITEMS)
      IF(NITEMS.NE.NDIM)THEN
          WRITE(6,*)'There is a mismatch in NDIM ',NITEMS
      ENDIF 
      IF(JOB.EQ.'w')NITEMS=NPOIN
      IFAIL = IXDRINT(IXDRS(1),NITEMS)
      IF(NITEMS.NE.NPOIN)THEN
          WRITE(6,*)'There is a mismatch in NPOIN ',NITEMS
      ENDIF 
      IF(JOB.EQ.'w')NITEMS=NELEM
      IFAIL = IXDRINT(IXDRS(1),NITEMS)
      IF(NITEMS.NE.NELEM)THEN
          WRITE(6,*)'There is a mismatch in NELEM ',NITEMS
      ENDIF 
      IF(JOB.EQ.'w')NITEMS=NBFAC
      write(6,*)'NBFAC = ',nbfac
      IFAIL = IXDRINT(IXDRS(1),NITEMS)
      IF(NITEMS.NE.NBFAC)THEN
          WRITE(6,*)'There is a mismatch in NBFAC ',NITEMS
      ENDIF 
#ifdef MPI
      IF(JOB.EQ.'w')NITEMS=NBINT
      IFAIL = IXDRint(ixdrs(1),NITEMS)
      IF(NITEMS.NE.NBINT)THEN
          WRITE(6,*)'There is a mismatch in NBINT ',NITEMS
      ENDIF 
#else
      NBINT = 0
#endif 
      IF(JOB.EQ.'w')NITEMS=NHOLE
      IFAIL = IXDRINT(IXDRS(1),NITEMS)
      IF(NITEMS.NE.NHOLE)THEN
          WRITE(6,*)'There is a mismatch in NHOLE ',NITEMS
      ENDIF 
#ifdef MPI
      IF(JOB.EQ.'w')NITEMS=NGHOST
      IFAIL = IXDRINT(IXDRS(1),NGHOST)
      IF(NITEMS.NE.NHOLE)THEN
          WRITE(6,*)'There is a mismatch in NGHOST ',NITEMS
      ENDIF 
#else
      NGHOST = 0
#endif 
C
C     .. Reading info from mesh file 004 (periodic bcs) ...
C
      IF(JOB.EQ."w")THEN
         IF(PERIODIC_MESH) THEN
             LFLAG = .TRUE.
         ELSE
             LFLAG = .FALSE.
         ENDIF
      ELSE ! read the mesh
         INQUIRE(FILE=FNAME3,EXIST=LFLAG)
         IF( (LFLAG.AND.(.NOT.PERIODIC_MESH)) .OR.
     &       (PERIODIC_MESH.AND.(.NOT.LFLAG)) )THEN
             WRITE(6,*)'There is a problem with periodicity ',
     &                 LFLAG,PERIODIC_MESH
             CALL EXIT(1)
         ENDIF
      ENDIF
C Not sure this is ok for periodic grids!!!
      IF(LFLAG)THEN
          PERIODIC_MESH = .TRUE.
          IXDRS(4) = INITXDR(FNAME3,JOB,.TRUE.)
          IF( IXDRS(4).LT.0 )THEN
              WRITE(ERRMSG(7:8),FMT="(I2.2)")IXDRS(4)
              WRITE(ERRMSG(15:72),FMT="(A57)")FNAME3(1:57)
              NERR = 3
              IOPT = 1
              CALL SETERR(ERRMSG,72,NERR,IOPT)
          ENDIF
          IF(JOB.EQ."w")NITEMS=NPNOD
          IFAIL = IXDRINT(IXDRS(4),NITEMS)
          IF(NITEMS.NE.NPNOD)THEN
              WRITE(6,*)'There is a mismatch in NPNOD ',NITEMS
          ENDIF 
          IFAIL = IXDRCLOSE(IXDRS(4))
      ELSE
          NPNOD=0
          PERIODIC_MESH = .FALSE.
      ENDIF
C
      NP = NPOIN+NGHOST+NPNOD
C
#ifdef PARCH_t3d
      IFAIL = IXDRRMAT64(IXDRS(1),NDIM*NP,VCORG)
#else
      IFAIL = IXDRDMAT(IXDRS(1),NDIM*NP,VCORG)
#endif
      IF (IVERBOSE.GE.1) CALL R8Mat_Print('General',' ',NDIM,NP,VCORG,
     +                   NDIM,'Nodal coordinates in vgrdpts',IFAIL)
      IF(JOB.EQ."w")THEN
         IFAIL = IXDRIMAT(IXDRS(1),NELEM*NOFVERT,ICELNOD)
         IFAIL = IXDRIMAT(IXDRS(1),NBFAC*3,IBNDPTR)
      ENDIF 
      IFAIL = IXDRCLOSE(IXDRS(1))
      RETURN
  400 FORMAT(10X,'Reading grid data from ',A)
  405 FORMAT(10X,'Writing grid data to ',A)
      END
