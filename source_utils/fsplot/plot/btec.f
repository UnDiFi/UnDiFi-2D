      SUBROUTINE BTEC(FILENAME,POINT,NDIM,NPOIN,ICELNOD,NOFVERT,
     &                NELEM,Z,NOFVAR,MOVE,IWRK)
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
      INTEGER NDIM,NPOIN,NELEM,NOFVERT,NOFVAR,IWRK
C
C     .. Array Arguments ..
C
      INTEGER*4   ICELNOD(NOFVERT,*),MOVE(IWRK)
      DOUBLE PRECISION POINT(NDIM,*),Z(NOFVAR,*)
C
C     .. Local Scalars ..
C
      CHARACTER*11  Elementtype
      CHARACTER*256  VarString
      CHARACTER*6   STRNG
      INTEGER*4 Vlength,IPOIN,IELEM,IFREQ,IFAIL
      INTEGER*4 J,LastChar,IVAR
      INTEGER*4 FileType
      INTEGER*4 Debug
      INTEGER*4 VIsDouble
      INTEGER*4 IsDouble
!     CHARACTER*(*) ZoneTitle
      INTEGER*4 ZoneType
      INTEGER*4 IMxOrNumPts
      INTEGER*4 JMxOrNumElements
      INTEGER*4 KMxOrNumFaces
      INTEGER*4 ICellMax
      INTEGER*4 JCellMax
      INTEGER*4 KCellMax
      DOUBLE PRECISION Solution_Time
      INTEGER*4 StrandID
      INTEGER*4 ParentZone
      INTEGER*4 IsBlock
      INTEGER*4 NumFaceConnections
      INTEGER*4 FaceNeighborMode
      INTEGER*4 TotalNumFaceNodes
      INTEGER*4 NumConnectedBoundaryFaces
      INTEGER*4 TotalNumBoundaryConnections
      INTEGER*4 PassiveVarList(15)
      INTEGER*4 ValueLocation(15)
      INTEGER*4 ShareVarFromZone(15)
      INTEGER*4 ShareConnectivityFromZone
      INTEGER*4 GridOnly,SolutionOnly,Full
      PARAMETER(GridOnly=1,SolutionOnly=2,Full=0)
C
C     .. Local Arrays ..
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
      INTEGER*4 TECINI111,TECZNE111,TECDAT111,TECNOD111,
     &TECEND111
C
C     .. Intrinsic Functions ..
C
      INTRINSIC INDEX,MAX0
C
C     .. Executable Statements ..
C
C...Set defaults
C
      FileType = Full
      FileType = GridOnly
      FileType = SolutionOnly
      Debug = 1
      VIsDouble = 1
      DO IVAR = 1, NOFVAR+NDIM
         PassiveVarList(Ivar) = 0
         ShareVarFromZone(Ivar) = 0
         ValueLocation(Ivar) = 1
      ENDDO
C
      IF( NDIM .EQ. 2 )THEN
           ZoneType=2
      ELSE
           ZoneType=4
      ENDIF
      IF(FileType.EQ.0.OR.FileType.EQ.GridOnly)THEN ! solution only
         IF( NDIM .EQ. 2 )THEN
           VarString = 'X Y'
           Vlength = 3
         ELSE
           VarString = 'X Y Z'
           Vlength = 5
         ENDIF
      ELSE
         Vlength = 0
      ENDIF
      IF(FileType.EQ.Full.OR.FileType.EQ.SolutionOnly)THEN
         STRNG(1:6) = " Z(xx)"
         IF    ( NOFVAR .GT. 5 )THEN
           DO IVAR = 1, NOFVAR
              WRITE(STRNG(4:5),FMT="(I2.2)")IVAR
              VarString = VarString(1:Vlength) // STRNG
              Vlength = Vlength + 6
           ENDDO
         ELSEIF( NOFVAR .EQ. 1 )THEN
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
           VarString = VarString(1:Vlength) // 
     &' Z(1) Z(2) Z(3) Z(4) Z(5)'
           Vlength = Vlength + 25
         ELSE
           WRITE(6,*)"Uh Oh! I don't know what to do with NOFVAR = ",
     &               NOFVAR
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      LastChar = INDEX(filename,CHAR(32))
C
   10 CONTINUE
      WRITE(6,155)
C
C     tecplot unformatted interface
C
      IFAIL =  TECINI111("No Title"//char(0), 
     & VarString(1:Vlength)//char(0),
     & FILENAME//char(0),
     & "."//char(0),
     & FileType,
     & Debug,
     & VIsDouble)
      WRITE(6,FMT=345)'TECINI111',IFAIL
      IF(IFAIL.NE.0) CALL EXIT(IFAIL)
      KMxOrNumFaces=0 ! does not apply to triangles tets
      ICellMax=0
      JCellMax=0
      KCellMax=0
      Solution_Time=1.d0
      StrandID=0
      StrandID=1
      ParentZone=0
      IsBlock=0 ! point format
      IsBlock=1
      NumFaceConnections=0
      FaceNeighborMode=0
      TotalNumFaceNodes=0 ! does not apply to triangles tets
      NumConnectedBoundaryFaces=0 ! does not apply to triangles tets
      TotalNumBoundaryConnections=0 ! does not apply to triangles tets
      ValueLocation=1
      ShareConnectivityFromZone=0
C
      IFAIL = TECZNE111("No ZoneTitle yet"//char(0),
     & ZoneType,
     & NPOIN,
     & NELEM,
     & KMxOrNumFaces,
     & ICellMax,
     & JCellMax,
     & KCellMax,
     & Solution_Time,
     & StrandID,
     & ParentZone,
     & IsBlock,
     & NumFaceConnections,
     & FaceNeighborMode,
     & TotalNumFaceNodes,
     & NumConnectedBoundaryFaces,
     & TotalNumBoundaryConnections,
     & PassiveVarList,
     & ValueLocation,
     & ShareVarFromZone,
     & ShareConnectivityFromZone)
      WRITE(6,FMT=345)'TECZNE111',IFAIL
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
C
      IsDouble = 1
C
      IF(FileType.EQ.Full.OR.FileType.EQ.GridOnly)THEN
         CALL TRANS(POINT,NDIM,NPOIN,NDIM*NPOIN,MOVE,IWRK,IFAIL) 
         WRITE(6,FMT=345)'TRANS',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         IFAIL = TECDAT111( NDIM*NPOIN, POINT, IsDouble)
         WRITE(6,FMT=345)'TECDAT111',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         CALL TRANS(POINT,NPOIN,NDIM,NDIM*NPOIN,MOVE,IWRK,IFAIL) 
         WRITE(6,FMT=345)'TRANS',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
      ENDIF
      IF(FileType.EQ.Full.OR.FileType.EQ.SolutionOnly)THEN
         CALL TRANS(Z,NOFVAR,NPOIN,NOFVAR*NPOIN,MOVE,IWRK,IFAIL) 
         WRITE(6,FMT=345)'TRANS',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         IFAIL = TECDAT111( NOFVAR*NPOIN, Z, IsDouble)
         WRITE(6,FMT=345)'TECDAT111',IFAIL
         IF(IFAIL.NE.0) CALL EXIT(IFAIL)
         CALL TRANS(Z,NPOIN,NOFVAR,NOFVAR*NPOIN,MOVE,IWRK,IFAIL) 
         WRITE(6,FMT=345)'TRANS',IFAIL
      ENDIF
C
      IF(FileType.EQ.Full.OR.FileType.EQ.GridOnly)THEN
         IFAIL = TECNOD111(ICELNOD)
         WRITE(6,FMT=345)'TECNOD111',IFAIL
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
      ENDIF
      IFAIL = TECEND111()
      WRITE(6,FMT=345)'TECEND111',IFAIL
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
C
C       Opening the tecplot file
C
      WRITE(6,180)
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
  345 FORMAT(5X,'Subroutine ',A9,' has returned IFAIL = ',I4)
  999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END
