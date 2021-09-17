      SUBROUTINE SETPTO( FNAME, KLIST, VLIST, DWKSP, IRANK, 
     &                   NDIM, NOFVAR )
C
C     $Id: setpto.f,v 1.12 2020/03/25 15:19:11 abonfi Exp $
C
C     This routine reads the total pressure profile
C
      IMPLICIT NONE
C
      INCLUDE 'bnd.h'
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'ibc8.com'
      INCLUDE 'flags.com'
      INCLUDE 'stream.com'
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE 
C
C     Scalar arguments:
C
      INTEGER NDIM,NOFVAR,K
      CHARACTER*(*) FNAME
      CHARACTER*72 ERRMSG
C
C     Array arguments:
C
      INTEGER KLIST(NLIST),IRANK(NLIST)
      DOUBLE PRECISION VLIST(nVarsInlet,NLIST),DWKSP(NLIST)
C
C     Local scalars:
C
      INTEGER NN,IFAIL,IOPT,NERR,IUNIT,I
      INTEGER I1MACH
C
C     If there is no file (file005.dat/pbcs???.dat) with a
C     inflow profile, we total pressure and temperature
C     at the inflow bndry as set to 1.d0
C     (because of the non-dimensionalisation) and the
C     flow direction is taken equal to the one specified
C     with -flow_angles
C     filling VLIST is not actually needed, since
C     it is never used when LREAD(1) = .FALSE.
C     (see subr ghost2.F) 
C
      INQUIRE(FILE=FNAME,EXIST=LREAD(1))
      IF( .NOT. LREAD(1) )THEN
          IF(MY_PE.EQ.0)WRITE(IWUNIT,FMT=400)
  400 FORMAT(/,5X,'SETTING UNIFORM TOTAL PRESSURE ON THE INLET FACE',/,
     +       5X,48("="))
C
C         set total pressure to one (because of the non-dimensionalisation)
C
          DO 3 I = 1,NLIST
               VLIST(1,I) = ONE ! total pressure/ref pressure
               VLIST(2,I) = ONE ! total temperature/ref temperature
               VLIST(3,I) = 1.d38 ! unused (backward compatibility)
               VLIST(4,I) = FLOWDIR(1)
               VLIST(5,I) = FLOWDIR(2)
               VLIST(6,I) = FLOWDIR(3)
    3     CONTINUE
          RETURN
C
      ELSE
          IF(MY_PE.EQ.0)WRITE(NOUT,1000)MY_PE,FNAME
      ENDIF
C
      CALL IINIT(NLIST,0,IRANK,1) 
C
      IUNIT=77
      OPEN(IUNIT,FILE=FNAME)
      READ(IUNIT,*)NN
      IF(NN.NE.NLIST)THEN
          WRITE(I1MACH(4),FMT=500)NN,FNAME,NLIST
          ERRMSG(1:15) = 'SETPTO, PROC # '
          WRITE(ERRMSG(16:19),FMT="(I4.4)")MY_PE
          NERR = 12
          IOPT = 1
          CALL SETERR(ERRMSG,19,NERR,IOPT)
      ENDIF 
C
C     +-----------------------------------------------------+
C     read the global mesh number and corresponding value
C     of total pressure (relative)
C     +-----------------------------------------------------+
C
      DO 1 I = 1 , NLIST
         READ(IUNIT,*)KLIST(I),(VLIST(K,I),K=1,nVarsInlet)
    1 CONTINUE
      CLOSE(IUNIT) 
C
C     +-----------------------------------------------------+
C     sort the array for increasing nodenumber
C     +-----------------------------------------------------+
!     CALL M01DBF(KLIST,1,NLIST,'Ascending',IRANK,IFAIL)
      CALL QSORTI(IRANK,NLIST,KLIST)
      CALL RNKIDX(IRANK,1,NLIST,IFAIL) ! converti in ranking
      CALL I4RANK(KLIST,1,NLIST,IRANK,IFAIL)
      DO 12 K = 1,nVarsInlet
         CALL DCOPY(NLIST,VLIST(K,1),nVarsInlet,DWKSP,1)
         CALL R8RANK(DWKSP,1,NLIST,IRANK,IFAIL)
         CALL DCOPY(NLIST,DWKSP,1,VLIST(K,1),nVarsInlet)
   12 CONTINUE
C
C
      RETURN
  500 FORMAT(5X,'THE NUMBER OF MESHPOINTS ',I7,' DECLARED IN FILE ',/,
     +A,/,' DOES NOT MATCH THE NO. OF ENTRIES ',I7,
     +' IN THE Index Set')
  550 FORMAT(5X,'MESHPOINT ',I7,' IN THE Index Set',/, 
     +' IS NOT AMONG THOSE DECLARED IN FILE ',/,
     +A)
 1000 FORMAT(/,5X,'READING INFLOW BOUNDARY CONDITIONS ON PE # ',I4,/,
     &       5X,48("=")/,5X,'FROM FILE: ',A,/)
      END
