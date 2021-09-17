      SUBROUTINE AVS(FILENAME,POINT,ICELNOD,Z)
C
C     IMPLICIT NONE
C
C This routine writes a tecplot formatted file
C
      INCLUDE 'constants'
      INCLUDE 'mesh_i4'
      INCLUDE 'int_flags'
      INCLUDE 'dim_flags'
      INCLUDE 'IO'
C
C     .. Scalar Arguments ..
C
      CHARACTER*(*) FILENAME
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION  POINT,Z
      INTEGER   ICELNOD
      DIMENSION POINT(DIM,1),Z(NOFVAR,1),ICELNOD(NOFVERT,1)
C
C     .. Local Scalars ..
C
      CHARACTER*3      Elementtype
      CHARACTER*30      VarString
      INTEGER*4 Debug,Vlength,IPOIN,IELEM,IFREQ
      INTEGER*4 J,LastChar
C
C     .. Local Arrays ..
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
C
C     .. Intrinsic Functions ..
C
      INTRINSIC INDEX,MAX0
C
C     .. Executable Statements ..
C
C...Set defaults
C
      num_cdata = 0
      num_mdata = 0
      mat_id = 1
C
      IF( DIM .EQ. 2 )THEN
        VarString = 'X Y'
        Vlength = 3
        Elementtype = 'tri'
      ELSE
        VarString = 'X Y Z'
        Vlength = 5
        Elementtype = 'tet'
      ENDIF
      IF    ( NOFVAR .EQ. 1 )THEN
	VarString = VarString(1:Vlength) // ' Z(1)' 
        Vlength = Vlength + 5
      ELSEIF( NOFVAR .EQ. 2 )THEN
	VarString = VarString(1:Vlength) // ' Z(1) Z(2)' 
        Vlength = Vlength + 10
      ELSEIF( NOFVAR .EQ. 3 )THEN
	VarString = VarString(1:Vlength) // ' p u v' 
        Vlength = Vlength + 6
      ELSEIF( NOFVAR .EQ. 4 )THEN
        VarString = VarString(1:Vlength) // ' Z(1) Z(2) Z(3) Z(4)'
        Vlength = Vlength + 20
      ELSEIF( NOFVAR .EQ. 5 )THEN
        VarString = VarString(1:Vlength) // ' Z(1) Z(2) Z(3) Z(4) Z(5)'
        Vlength = Vlength + 25
      ELSE
	WRITE(6,*)"Uh Oh! I don't know what to do with NOFVAR = ",NOFVAR
	STOP
      ENDIF
C
C
      WRITE(6,155)
C
C       tecplot formatted interface
C
C       Opening the tecplot file
C
      OPEN(3,FILE=filename,STATUS='UNKNOWN')
C
C     writing nodes
C
      IFREQ = MAX0( NPOIN/20 , 1 )
      WRITE(6,165)
      WRITE(3,*)NPOIN,NELEM,NOFVAR,num_cdata,num_mdata
      DO 1 IPOIN = 1,NPOIN
        IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
        WRITE(3,*)IPOIN,(POINT(J,IPOIN),J=1,DIM)
    1 CONTINUE
      WRITE(6,175)
C
C     writing elements
C
      IFREQ = MAX0( NELEM/20 , 1 )
      WRITE(6,170)
      DO 2 IELEM = 1 , NELEM
        IF((IELEM/IFREQ)*IFREQ .EQ. IELEM)WRITE(*,111)
        WRITE(3,*)IELEM,mat_id,elementtype,
     >            (ICELNOD(J,IELEM),J=1,NOFVERT)
    2 CONTINUE
C
C     num data components, size of each component
C
C     doesn't work for scalar problems
C
      num_comp = nofvar - dim + 1
      WRITE(3,*)num_comp,(1,i=1,num_comp-1),dim
      WRITE(3,*)'density, dimensionless'
      IF(num_comp.EQ.3)WRITE(3,*)'enthalpy, dimensionless'
      WRITE(3,*)'velocity vector, dimensionless'
      DO 4 IPOIN = 1,NPOIN
        IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
        WRITE(3,*)IPOIN,(Z(J,IPOIN),J=1,NOFVAR)
    4 CONTINUE
C
      WRITE(6,175)
      WRITE(6,185)FILENAME
C
      CLOSE(3)
C
C     WRITE(6,180)
C
      RETURN
C
  111 FORMAT('.',$)
* 111 FORMAT('x')
  135 FORMAT(5X,'Writing the coordinates ... ')
  140 FORMAT(5X,'IFMT = ',I2,' MUST be 1 or 2')
  145 FORMAT(5X,'Writing the variables ... ')
  155 FORMAT(5X,'I am writing the AVS file ... ',/)
  165 FORMAT(5X,'Writing nodes ',$)
  170 FORMAT(5X,'Writing elements ',$)
  175 FORMAT(' done !',/)
  180 FORMAT(5X,'AVS file WRITTEN',/)
  185 FORMAT(/,5X,'AVS file WRITTEN to ... ',A60/)
  999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END
