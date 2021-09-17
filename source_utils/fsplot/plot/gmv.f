      SUBROUTINE GMV(FILENAME,POINT,NDIM,NPOIN,ICELNOD,NOFVERT,NELEM,
     &Z,NOFVAR)
C
      IMPLICIT NONE
C
C This routine writes a tecplot formatted file
C
      INCLUDE 'constants'
      INCLUDE 'IO'
C
C     .. Scalar Arguments ..
C
      CHARACTER*(*) FILENAME
C
C     .. Array Arguments ..
C
      INTEGER NDIM,NPOIN,NELEM,NOFVERT,NOFVAR
      DOUBLE PRECISION  POINT,Z
      INTEGER   ICELNOD
      DIMENSION POINT(NDIM,*),Z(NOFVAR,*),ICELNOD(NOFVERT,*)
C
C     .. Local Scalars ..
C
      CHARACTER*5  Elementtype
      CHARACTER*256  VarString
      CHARACTER*256  VarFmt
      INTEGER*4 IPOIN,IELEM,IFREQ
      INTEGER*4 J,IVAR,IUNIT,data_type
      PARAMETER(data_type=1)
C
C     .. Local Arrays ..
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
      DOUBLE PRECISION PRESSC,MACHNO,PTOT
C
C     .. Intrinsic Functions ..
C
      INTRINSIC INDEX,MAX0
C
C     .. Executable Statements ..
C
C...Set defaults
C
      IF( NDIM .EQ. 2 )THEN
        Elementtype = 'tri 3'
      ELSE
        Elementtype = 'tet 4'
      ENDIF
C
C
      WRITE(6,155)
C
C       tecplot formatted interface
C
C       Opening the tecplot file
C
      IUNIT = 3
      OPEN(IUNIT,FILE=filename,STATUS='UNKNOWN')
      WRITE(IUNIT,FMT="(A14)")'gmvinput ascii'
C
C       writing the nodes
C
      IFREQ = MAX0( NPOIN/20 , 1 )
      WRITE(6,165)
      WRITE(IUNIT,FMT="(A5,1X,I6)")'nodev',NPOIN
      IF( NDIM .EQ. 2 )THEN
         DO 1 IPOIN = 1,NPOIN
           IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
           WRITE(IUNIT,*)(POINT(J,IPOIN),J=1,NDIM),0.d0
    1 CONTINUE
      ELSE
         DO 7 IPOIN = 1,NPOIN
           IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
           WRITE(IUNIT,*)(POINT(J,IPOIN),J=1,NDIM)
    7 CONTINUE
      ENDIF
      WRITE(6,175)
C
C       writing the elements
C
      IFREQ = MAX0( NELEM/20 , 1 )
      WRITE(IUNIT,FMT="(A5,1X,I7)")'cells',NELEM
      WRITE(6,170)
      DO 2 IELEM = 1 , NELEM
        IF((IELEM/IFREQ)*IFREQ .EQ. IELEM)WRITE(*,111)
        WRITE(IUNIT,FMT="(A5)")elementtype
        WRITE(IUNIT,*)(ICELNOD(J,IELEM),J=1,NOFVERT)
    2 CONTINUE
C
      WRITE(6,175)
C
C       writing the variables
C
      WRITE(IUNIT,FMT="(A8)")'variable'
C
      DO IVAR = 1,NOFVAR
        WRITE(IUNIT,FMT="(A2,I1,A1,1X,I1)")'Z(',IVAR,')',data_type
        WRITE(IUNIT,*)(Z(IVAR,IPOIN),IPOIN=1,NPOIN)
      END DO
C
      IF( NOFVAR .EQ. 1 )THEN
         WRITE(IUNIT,FMT="(A7)")'endvars'
         GOTO 100
      ENDIF
      IF( NOFVAR .EQ.(NDIM+2) )THEN ! compressible
        WRITE(IUNIT,FMT="(A1,1X,I1)")'P',data_type
        WRITE(IUNIT,*)(PRESSC(NDIM,Z(1,IPOIN)),IPOIN=1,NPOIN)
        WRITE(IUNIT,FMT="(A4,1X,I1)")'Mach',data_type
        WRITE(IUNIT,*)(MachNO(NDIM,Z(1,IPOIN)),IPOIN=1,NPOIN)
        WRITE(IUNIT,FMT="(A3,1X,I1)")'Rho',data_type
        WRITE(IUNIT,*)((Z(1,IPOIN)*Z(1,IPOIN)),IPOIN=1,NPOIN)
        WRITE(IUNIT,FMT="(A1,1X,I1)")'H',data_type
        WRITE(IUNIT,*)((Z(2,IPOIN)/Z(1,IPOIN)),IPOIN=1,NPOIN)
        WRITE(IUNIT,FMT="(A4,1X,I1)")'Ptot',data_type
        WRITE(IUNIT,*)(PTOT(NDIM,Z(1,IPOIN)),IPOIN=1,NPOIN)
      ENDIF
      WRITE(IUNIT,FMT="(A7)")'endvars'
      WRITE(IUNIT,FMT="(A8,1X,I1)")'velocity',data_type
C
      IF    ( NOFVAR.EQ.(NDIM+2) )THEN ! compressible
         DO IVAR = 3,2+NDIM
           WRITE(IUNIT,*)(Z(IVAR,IPOIN)/Z(1,IPOIN),IPOIN=1,NPOIN)
         END DO
         IF(NDIM.EQ.2)  WRITE(IUNIT,*)(0.d0,IPOIN=1,NPOIN)
      ELSEIF( NOFVAR.EQ.(NDIM+1) )THEN 
         DO IVAR = 2,1+NDIM
           WRITE(IUNIT,*)(Z(IVAR,IPOIN),IPOIN=1,NPOIN)
         END DO
         IF(NDIM.EQ.2)  WRITE(IUNIT,*)(0.d0,IPOIN=1,NPOIN)
      END IF
C
  100 WRITE(IUNIT,FMT="(A6)")'endgmv'

      CLOSE(IUNIT)
      WRITE(6,185)FILENAME
C
C
      RETURN
C
  111 FORMAT('.',$)
  135 FORMAT(5X,'Writing the coordinates ... ')
  140 FORMAT(5X,'IFMT = ',I2,' MUST be 1 or 2')
  145 FORMAT(5X,'Writing the variables ... ')
  155 FORMAT(5X,'I am writing the GMV file ... ',/)
  165 FORMAT(5X,'Writing nodes ',$)
  170 FORMAT(5X,'Writing elements ',$)
  175 FORMAT(' done !',/)
  180 FORMAT(5X,'GMV file WRITTEN',/)
  185 FORMAT(/,5X,'GMV file WRITTEN to ... ',A60/)
  999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END

      DOUBLE PRECISION FUNCTION PRESSC( NDIM , ZROE )
C
C    .. This function computes PRESSURE from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
C
      INTEGER NDIM
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION TEMP
C
      TEMP       = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF * TEMP
      PRESSC = GM1OG * ( ZROE(1)*ZROE(2) - TEMP )
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION MACHNO( NDIM , ZROE )
C
C    .. This function computes Mach number from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
C
      INTEGER NDIM
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION TEMP,ASQR,KINE
C
      TEMP       = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(5)*ZROE(5)
      KINE = HALF * TEMP / (ZROE(1)*ZROE(1))
      ASQR = GM1 * ( ZROE(2)/ZROE(1) - KINE )
      MACHNO = SQRT( 2.d0 * KINE/ASQR )
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION PTOT( NDIM , ZROE )
C
C    .. This function computes total pressure from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
C
      INTEGER NDIM
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION TEMP,ASQR,KINE,MACHSQR
      DOUBLE PRECISION PRESSC
C
      TEMP       = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(5)*ZROE(5)
      KINE = HALF * TEMP / (ZROE(1)*ZROE(1))
      ASQR = GM1 * ( ZROE(2)/ZROE(1) - KINE )
      MACHSQR = 2.d0 * KINE/ASQR 
      TEMP = (1.d0+0.5*GM1*MACHSQR)**(GAM/GM1)
      PTOT = PRESSC( NDIM, ZROE ) * TEMP
C
      RETURN
      END
C
