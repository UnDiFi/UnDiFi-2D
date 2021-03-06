      SUBROUTINE TECPLOT(FILENAME,POINT,NDIM,NPOIN,ICELNOD,NOFVERT,
     &                   NELEM,Z,NOFVAR)
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
      INTEGER NDIM,NPOIN,NELEM,NOFVERT,NOFVAR
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION  POINT,Z
      INTEGER   ICELNOD
      DIMENSION POINT(NDIM,1),Z(NOFVAR,1),ICELNOD(NOFVERT,1)
C
C     .. Local Scalars ..
C
      CHARACTER*13  Elementtype
      CHARACTER*256  VarString
      CHARACTER*6   STRNG
      INTEGER*4 Debug,Vlength,IPOIN,IELEM,IFREQ
      INTEGER*4 J,LastChar,IVAR,NDOF
      INTEGER*4 GridOnly,SolutionOnly,Full,FileType
      PARAMETER(GridOnly=1,SolutionOnly=2,Full=0)
      CHARACTER*8 FTYPE(0:2)
      DATA FTYPE/'FULL    ','GRID    ','SOLUTION'/
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
      FileType = SolutionOnly 
      FileType = GridOnly 
      FileType = FULL 
      IF(FileType.EQ.Full)THEN
         NDOF = NOFVAR+NDIM
      ELSEIF(FileType.EQ.SolutionOnly)THEN
         NDOF = NOFVAR
      ELSEIF(FileType.EQ.GridOnly)THEN
         NDOF = NDIM
      ELSE
         WRITE(6,*)'Unknown FileType ',FileType
         CALL EXIT(1)
      ENDIF
      IF( NDIM .EQ. 2 )THEN
        Elementtype = 'FETRIANGLE'
      ELSE
        Elementtype = 'FETETRAHEDRON'
      ENDIF
      STRNG(1:6) = " Z(xx)"
      Vlength = 0
      IF(FileType.EQ.Full.OR.FileType.EQ.GridOnly)THEN ! solution only
         IF( NDIM .EQ. 2 )THEN
           VarString = 'X Y'
           Vlength = 3
         ELSE
           VarString = 'X Y Z'
           Vlength = 5
         ENDIF
      ENDIF
      IF(FileType.EQ.Full.OR.FileType.EQ.SolutionOnly)THEN
         IF    ( NOFVAR .GT. 5 )THEN
           DO IVAR = 1, NOFVAR
              WRITE(STRNG(4:5),FMT="(I2.2)")IVAR
              VarString = VarString(1:Vlength) // STRNG
              Vlength = Vlength + 6
           ENDDO
         ELSEIF( NOFVAR .EQ. 1 )THEN
           VarString = VarString(1:Vlength) // ' Z1' 
           Vlength = Vlength + 3
         ELSEIF( NOFVAR .EQ. 2 )THEN
            VarString = VarString(1:Vlength) // ' Z1 Z2' 
           Vlength = Vlength + 6
         ELSEIF( NOFVAR .EQ. 3 )THEN
            VarString = VarString(1:Vlength) // ' p u v' 
           Vlength = Vlength + 6
         ELSEIF( NOFVAR .EQ. 4 )THEN
           VarString = VarString(1:Vlength) // ' Z1 Z2 Z3 Z4'
           Vlength = Vlength + 12
         ELSEIF( NOFVAR .EQ. 5 )THEN
           VarString = VarString(1:Vlength) // 
     &' Z1 Z2 Z3 Z4 Z5'
           Vlength = Vlength + 15
         ELSE
           WRITE(6,*)"Uh Oh! I don't know what to do with NOFVAR = ",
     &NOFVAR
           CALL EXIT(1)
         ENDIF
      ENDIF
C
      LastChar = INDEX(filename,CHAR(32))
C
      WRITE(6,155)
C
C       tecplot formatted interface
C
C       Opening the tecplot file
C
      OPEN(3,FILE=filename,STATUS='UNKNOWN')
C
C       writing the nodes
C
      IFREQ = MAX0( NPOIN/20 , 1 )
      WRITE(6,165)
      WRITE(3,*)'TITLE = "Full grid"'
      WRITE(3,*)'VARIABLES = ',VarString(1:Vlength)
!     WRITE(3,*)'FILETYPE = ',FTYPE(FileType)
      WRITE(3,120)'"sampletext"',NPOIN,NELEM,elementtype
 120  FORMAT("ZONE T=",A,", DATAPACKING=POINT, NODES=",I7,
     2", ELEMENTS=",I8,", ZONETYPE=",A)
C
      WRITE(3,FMT="(A4,$)")'DT=('
      DO IVAR = 1,NDOF
         WRITE(3,FMT="(A7,$)")'DOUBLE '
      ENDDO
      WRITE(3,*)')'
C
      IF(FileType.EQ.Full)THEN
         DO IPOIN = 1,NPOIN
           IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
           WRITE(3,*)(POINT(J,IPOIN),J=1,NDIM),(Z(J,IPOIN),J=1,NOFVAR)
         ENDDO
      ELSEIF(FileType.EQ.SolutionOnly)THEN
         DO IPOIN = 1,NPOIN
           IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
           WRITE(3,*)(Z(J,IPOIN),J=1,NOFVAR)
         ENDDO
      ELSEIF(FileType.EQ.GridOnly)THEN
         DO IPOIN = 1,NPOIN
           IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)
           WRITE(3,*)(POINT(J,IPOIN),J=1,NDIM)
         ENDDO
      ENDIF
C
      WRITE(6,175)
C
C       writing the elements
C
      IFREQ = MAX0( NELEM/20 , 1 )
      IF(FileType.EQ.Full.OR.FileType.EQ.GridOnly)THEN ! solution only
         WRITE(6,170)
         DO 2 IELEM = 1 , NELEM
           IF((IELEM/IFREQ)*IFREQ .EQ. IELEM)WRITE(*,111)
           WRITE(3,*)(ICELNOD(J,IELEM),J=1,NOFVERT)
    2    CONTINUE
         WRITE(6,175)
      ENDIF
      WRITE(6,185)FILENAME
C
      CLOSE(3)
C
*     WRITE(6,180)
C
      RETURN
C
  111 FORMAT('.',$)
* 111 FORMAT('x')
! 135 FORMAT(5X,'Writing the coordinates ... ')
! 140 FORMAT(5X,'IFMT = ',I2,' MUST be 1 or 2')
! 145 FORMAT(5X,'Writing the variables ... ')
  155 FORMAT(5X,'I am writing the Tecplot file ... ',/)
  165 FORMAT(5X,'Writing nodes ',$)
  170 FORMAT(5X,'Writing elements ',$)
  175 FORMAT(' done !',/)
! 180 FORMAT(5X,'Tecplot file WRITTEN',/)
  185 FORMAT(/,5X,'Tecplot file WRITTEN to ... ',A60/)
! 999 FORMAT(5X,A6,' has returned an error message, IFAIL = ',I2)
      END
