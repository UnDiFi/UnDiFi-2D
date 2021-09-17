      SUBROUTINE MEDIAN(VMEDIAN,NPOIN,NGHOST,NPNOD,NP,VOLUME,ICELNOD,
     &MAP,NOFVERT,NELEM)
C
C     $Id: median.f,v 1.3 2008/06/10 10:13:12 abonfi Exp $
C
C     Routine for computing the median dual cell ..
C     note that ICELNOD has not yet been re-mapped
C     at the time this routine is being called
C     otherwise we could address ICELNOD(IVERT,NELEM+IELEM)
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'io.com'
C
C     .. Scalar Arguments ..
      INTEGER NELEM,NP,NOFVERT,NPOIN,NGHOST,NPNOD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VMEDIAN(*),VOLUME(NELEM)
      INTEGER ICELNOD(NOFVERT,NELEM),MAP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP,HELP,TVOLUME,VOLMIN,VOLMAX
      INTEGER IELEM,IPOIN,IVERT,I,JPOIN
C     ..
C
      WRITE(NOUT,2000)
2000  FORMAT(//' COMPUTATION OF THE MEDIAN DUAL CELL VOLUMES'/' ',31('='
     &)/)
      HELP = ONE/REAL(NOFVERT)
      CALL DINIT(NP,ZERO,VMEDIAN,1)
      DO 100 IELEM = 1,NELEM
          TEMP = VOLUME(IELEM)!*HELP
c what about periodic meshes ?!?!
          DO 90 IVERT = 1,NOFVERT
              IPOIN = ICELNOD(IVERT,IELEM)
              VMEDIAN(IPOIN) = VMEDIAN(IPOIN) + TEMP
   90     CONTINUE
  100 CONTINUE
      CALL DSCAL(NP,HELP,VMEDIAN,1)
C
C account for periodicity
C
      DO 300 I = 1, NPNOD
         IPOIN = (NPOIN+NGHOST)+I
         JPOIN = MAP(I)
         HELP = VMEDIAN(IPOIN)
         TEMP = VMEDIAN(JPOIN)
         VMEDIAN(IPOIN) = VMEDIAN(IPOIN) + TEMP
         VMEDIAN(JPOIN) = VMEDIAN(JPOIN) + HELP
  300 CONTINUE
C
      TVOLUME = ZERO
      VOLMIN = 1.D+38
      VOLMAX = - ONE
      DO 200 IPOIN = 1,NPOIN+NGHOST
          HELP = VMEDIAN(IPOIN)
          TVOLUME = TVOLUME + HELP
          VOLMIN = MIN(HELP,VOLMIN)
          VOLMAX = MAX(HELP,VOLMAX)
!         write(6,*)ipoin,help
  200 CONTINUE
!     DO IPOIN = 1,NPOIN+NGHOST+NPNOD
!        VMEDIAN(IPOIN) = ONE/1600.d0
!     enddo
      WRITE(NOUT,500)TVOLUME,VOLMIN,VOLMAX,VOLMAX/VOLMIN
C
C     I/O FORMATS
C
  500 FORMAT(/10X,'TOTAL VOLUME ......... ',E10.3/, 10X,
     +'MIN.  VOLUME ......... ',E10.3/, 10X,'MAX.  VOLUME ......... ',
     +E10.3/, 10X,'MAX/MIN RATIO ........ ',E10.3/)
C
      RETURN

      END
