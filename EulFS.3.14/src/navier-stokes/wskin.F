      SUBROUTINE WSKIN(TAUW,HEAT,IWORK,NWFAC,IBNDFAC,NBODY6,FILENAME)
C
C     Write skin friction data to a file
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'stream.com'
C
C     .. Scalar Arguments ..
      INTEGER NBODY6,NWFAC
      CHARACTER FILENAME* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION TAUW(NWFAC),HEAT(NWFAC)
      INTEGER IWORK(2,NWFAC),IBNDFAC(3,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IC,IELEM,IFACE,IFRST,ILAST,IVERT,IXDRS,IFAIL
C     ..
C     .. External Functions ..
      INTEGER INITXDR
      EXTERNAL INITXDR
      INTEGER IXDRCLOSE,IXDRIMAT,IXDRINT
      EXTERNAL IXDRCLOSE,IXDRIMAT,IXDRINT
#ifdef PARCH_t3d
      INTEGER IXDRRMAT64
      EXTERNAL IXDRRMAT64
#else
      INTEGER IXDRDMAT
      EXTERNAL IXDRDMAT
#endif
C     ..
C
C     do not scale, left to the post-processor;
C     now TAUW is in fact the wall shear stress 
C
      IXDRS = INITXDR(FILENAME,'w',.FALSE.)
C
      IFAIL = IXDRINT(IXDRS,NBODY6)
      IFAIL = IXDRINT(IXDRS,NWFAC)
C
      IC = 0
      DO 10 I = 1,NBODY6
C
          IFRST = IBGN(I)
          ILAST = IEND(I)
C
C     write number of faces for the current body
C
          IFAIL = IXDRINT(IXDRS,ILAST-IFRST+1)
C
          DO 10 IFACE = IFRST,ILAST
              IELEM = IBNDFAC(1,IFACE)
              IVERT = IBNDFAC(2,IFACE)
              IC = IC + 1
              IWORK(1,IC) = IELEM
              IWORK(2,IC) = IVERT
   10 CONTINUE

      IF (IC.NE.NWFAC) THEN
          WRITE (6,FMT=*) 'Ic NFAC',IC,NWFAC,NBODY6
          STOP

      ENDIF
      IFAIL = IXDRIMAT(IXDRS,2*IC,IWORK)
#ifdef PARCH_t3d
      IFAIL = IXDRRMAT64(IXDRS,IC,TAUW)
      IFAIL = IXDRRMAT64(IXDRS,IC,HEAT)
#else
      IFAIL = IXDRDMAT(IXDRS,IC,TAUW)
      IFAIL = IXDRDMAT(IXDRS,IC,HEAT)
#endif
C
      IFAIL = IXDRCLOSE(IXDRS)
      WRITE(NOUT,9996)FILENAME
C
      RETURN
 9996 FORMAT (/,5X,'Skin friction file has been written to ',A)

      END
