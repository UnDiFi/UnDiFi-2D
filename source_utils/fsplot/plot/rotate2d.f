      SUBROUTINE rotate2d(FILENAME,POINT,ndim,ICELNOD,nofvert,nelem,
     +nhole,Z,nofvar,npoin,npnod,ibndfac,nbfac)
C
      IMPLICIT NONE
C
C This routine writes a tecplot formatted file
C
      INTEGER ndim,nofvert,nofvar,npoin,npnod,nelem,nbfac,nhole,
     +ifail,ixdrs
      DOUBLE PRECISION seno,coseno,x,y,ux,uy
      DOUBLE PRECISION xle,yledn,yleup,gamma,deg2rad,pi
      parameter(xle=-0.471088,yledn=-0.178009,yleup=yledn,gamma=20.5)
C
C     .. Scalar Arguments ..
C
      CHARACTER*(*) FILENAME
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION  POINT,Z
      INTEGER   ICELNOD,ibndfac
      DIMENSION
     +POINT(NDIM,NPOIN),
     +Z(NOFVAR,NPOIN),ICELNOD(NOFVERT,NPOIN),ibndfac(3,nbfac)
C
C     .. Local Scalars ..
C
      INTEGER*4 Debug,Vlength,IPOIN,IELEM,IFREQ
      INTEGER*4 J,LastChar
C
C     .. Local Arrays ..
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
      INTEGER INITXDR,IXDRINT,IXDRIMAT,IXDRCLOSE,IXDRDMAT
C
C     .. Intrinsic Functions ..
C
      INTRINSIC INDEX,MAX0
C
C     .. Executable Statements ..
C
C...Set defaults
C
      pi = acos(-1.d0)
      deg2rad = pi/180.d0
      seno = sin(gamma*deg2rad)
      coseno = cos(gamma*deg2rad)
C
      WRITE(6,*)'Performing a rotation'
C
C       tecplot formatted interface
C
C       Opening the tecplot file
C
C
C       writing the nodes
C
      DO 1 IPOIN = 1,NPOIN
        x = POINT(1,IPOIN)
        y = POINT(2,IPOIN)
        POINT(1,IPOIN) =  (x-xle)*coseno+(y-yledn)*seno
        POINT(2,IPOIN) = -(x-xle)*  seno+(y-yledn)*coseno
        ux = Z(2,IPOIN)
        uy = Z(3,IPOIN)
        Z(2,IPOIN) =  (ux)*coseno+(uy)*seno
        Z(3,IPOIN) = -(ux)*  seno+(uy)*coseno
    1 CONTINUE
        WRITE(6,*)'Done!! '
C
C       writing the elements
C
C
!     call SOLZNE(FILENAME,Z,NOFVAR,NPOIN,"w")
C
C
C    Writing the mesh & connectivity file
C
!     ixdrs = INITXDR( 'file101.dat' , 'w' ,.FALSE.)
!     IFAIL = IXDRINT( ixdrs , NDIM )
!     IFAIL = IXDRINT( ixdrs , NPOIN-NPNOD )
!     IFAIL = IXDRINT( ixdrs , NELEM )
!     IFAIL = IXDRINT( ixdrs , NBFAC )
!     IFAIL = IXDRINT( ixdrs , NHOLE )
!     IFAIL = IXDRDMAT( ixdrs , NDIM*NPOIN , POINT )
!     IFAIL = IXDRIMAT( ixdrs , NOFVERT*NELEM , ICELNOD )
!     IFAIL = IXDRIMAT( ixdrs , 3*NBFAC , IBNDFAC )
!     IFAIL = IXDRCLOSE( ixdrs )
C
*     WRITE(6,180)
C
      RETURN
C
  111 FORMAT('.',$)
* 111 FORMAT('x')
  135 FORMAT(5X,'Writing the coordinates ... ')
  140 FORMAT(5X,'IFMT = ',I2,' MUST be 1 or 2')
  145 FORMAT(5X,'Writing the variables ... ')
  155 FORMAT(5X,'I am writing the Tecplot file ... ',/)
  165 FORMAT(5X,'Writing nodes ',$)
  170 FORMAT(5X,'Writing elements ',$)
  175 FORMAT(' done !',/)
  180 FORMAT(5X,'Tecplot file WRITTEN',/)
  185 FORMAT(/,5X,'Tecplot file WRITTEN to ... ',A60/)
  999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END
