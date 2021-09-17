C
C -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
C
      SUBROUTINE bndplot(FILENAME,POINT,VSCL,ICELNOD,IBNDFAC,
     +IWKSP,IRANK,LIWK)
C
      IMPLICIT NONE
C
C This routine writes a tecplot formatted file
C
      INCLUDE 'constants'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd'
      INCLUDE 'mesh_i4'
      INCLUDE 'int_flags'
      INCLUDE 'dim_flags'
      INCLUDE 'IO'
C
C     .. Scalar Arguments ..
C
      CHARACTER*(*) FILENAME
      INTEGER LIWK
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION POINT,VSCL
      INTEGER ICELNOD,IBNDFAC,IWKSP,IRANK
      DIMENSION POINT(DIM,1),VSCL(NOFVAR,1),ICELNOD(NOFVERT,1),
     &IBNDFAC(3,1),IWKSP(LIWK),IRANK(LIWK)
C
C     .. Local Scalars ..
C
      CHARACTER*11 Elementtype
      CHARACTER*80 VarString
      INTEGER Vlength
      INTEGER LastChar,ISUM,IBC,IBDFACE
      INTEGER I,J,IGAOK,IPOIN,INODE,IELEM,IFREQ,IVERT,IFAIL,IFACE,
     +JCOLOR,N,IC,Lo,Me,Up,NSIZE,IVAR
C
C     .. Local Arrays ..
C
      INTEGER IFRST(0:NCOLOR),ILAST(0:NCOLOR),NBNODE(0:NCOLOR),
     1ICNVT(3)
      CHARACTER*6   STRNG
C
C     .. Pointers ..
C
C#ifdef Linux
CINTEGER IRANK(1000),IWKSP(1000)
C#else
CINTEGER IRANK(1),IWKSP(1)
CPOINTER(PRANK,IRANK)
CPOINTER(PWKSP,IWKSP)
C#endif
C
C     .. External Subroutines ..
C
C     .. External Functions ..
C
      INTEGER  ICYCL,Binsrch
      EXTERNAL ICYCL,Binsrch
C#ifdef DEC
C      INTEGER*8 MALLOC
C      EXTERNAL  MALLOC
C#endif
C#ifdef HP-UX
C
C     .. Intrinsic Functions ..
C
C      INTRINSIC MALLOC
C#endif
C
      INTRINSIC INDEX,MAX0,FLOAT
C
      INTEGER NBFAC
      EQUIVALENCE(nBoundaryFaces,NBFAC)
C
      DATA IFAIL,NSIZE /2*0/
C
C     IF(NOFVAR.GT.5) THEN
C         WRITE(6,*)'Not writing bndry tecplot file; NOFVAR too big'
C         RETURN
C     ENDIF
C
C     .. Executable Statements ..
C
      STRNG(1:4) = " Zxx"
      DO 16 IBDFACE = 1 , NBFAC
         JCOLOR = IBNDFAC(3,IBDFACE)
         IF( JCOLOR .LT. 0 .OR. JCOLOR .GT. NCOLOR )THEN
            WRITE(NOUT,500)3,IBDFACE,JCOLOR
            CALL EXIT(1)
         ENDIF
         IBC = ICOLOR(JCOLOR)
         IF( IBC .LT. 0 .OR. IBC .GT. NBTYPE )THEN
            WRITE(NOUT,510)IBC,JCOLOR
            CALL EXIT(1)
         ENDIF
         MCOLOR(JCOLOR) = MCOLOR(JCOLOR) + 1 
   16 CONTINUE
C
      NSIZE = NBFAC
C
C     .. Allocates memory for the two work arrays: 
C        note that NSIZE = NBFAC might not be enough ....
C
C#ifdef Linux
C#else
C      PWKSP = MALLOC(4*NSIZE)
C      PRANK = MALLOC(4*NSIZE)
C#endif
C
C     .. Colors are copied into the work array ..
C
      IF( NBFAC .GT. LIWK )THEN
          WRITE(NOUT,888)
          CALL EXIT(1)
      ENDIF
      DO 10 IFACE = 1 , NBFAC
         IWKSP(IFACE) = IBNDFAC(3,IFACE)
   10 CONTINUE
C
C     .. which is ranked according to the different colors ..
C
      CALL QSORTI(IRANK,NBFAC,IWKSP)
c
!     do iface = 1 , NBFAC
!        igaok = irank(iface)
!        WRITE(16,*)iface,(ibndfac(i,igaok),i=1,3)
!     enddo
C
      IF( DIM .EQ. 2 )THEN
         VarString = 'X Y '
         Vlength = 4
      ELSE
         VarString = 'X Y Z '
         Vlength = 6
      ENDIF
        DO IVAR = 1, NOFVAR
           WRITE(STRNG(3:4),FMT="(I2.2)")IVAR
           VarString = VarString(1:Vlength) // STRNG
           Vlength = Vlength + 4
        ENDDO
        VarString = VarString(1:Vlength) // ' COLOR'
        Vlength = Vlength + 6
      Elementtype = 'TRIANGLE'
C
      LastChar = INDEX(filename,CHAR(32))
C
      WRITE(6,155)
C
C	tecplot formatted interface
C
C	Opening the tecplot file
C
      OPEN(3,FILE=filename,STATUS='UNKNOWN')
C
C	writing the nodes
C
      IFREQ = MAX0( NPOIN/20 , 1 )
ccc   WRITE(6,165)
C
C     .. Write header ..
C
      WRITE(3,*)'TITLE = "Full grid"'
      WRITE(3,*)'VARIABLES = ',VarString(1:Vlength)
C
C     .. Sets pointers for the first and last boundary
C        face of a given boundary color ..
C
      IC = 1
      ISUM = 0
      DO 1 JCOLOR = 0 , NCOLOR
         N = MCOLOR(JCOLOR)
         ISUM = ISUM + N
         IF( N .EQ. 0 )THEN
            IFRST(JCOLOR) = -1
            ILAST(JCOLOR) = -1
         ELSE
            IFRST(JCOLOR) = IC
            IC = IC + N
            ILAST(JCOLOR) = IC - 1
            write(6,*)IFRST(JCOLOR) , ILAST(JCOLOR)
         ENDIF
    1 CONTINUE
      write(6,*)'Pointers have been set'
      write(6,*)'Maximum number of faces is ',ISUM,' NBFAC = ',NBFAC,
     &' LIWK = ',LIWK
C
C     .. Outermost loop over the boundary colors ..
C
      DO 5 JCOLOR = 0 , NCOLOR
C
C     .. N is the number of bnd. faces colored JCOLOR ..
C
         N = MCOLOR(JCOLOR)
         NBNODE(JCOLOR) = 0
         IF( N .EQ. 0)GOTO 5
C
C     .. loop over the faces colored JCOLOR ..
C
         DO 7 IFACE = IFRST(JCOLOR) , ILAST(JCOLOR)
!        write(6,*)'Pointers are :',jcolor,IFRST(JCOLOR),ILAST(JCOLOR)
            IGAOK = IRANK(IFACE)
            IELEM = IBNDFAC(1,IGAOK)
            IVERT = IBNDFAC(2,IGAOK)
!           write(6,*)'check :',igaok,ielem,ivert
C
C     .. loop over the vertices of the boundary face ..
C
            DO 8 I =  1, DIM
               INODE = ICELNOD( ICYCL( I + IVERT , NOFVERT ) , IELEM )
C
C     .. checks whether node INODE is already in the list (IWKSP) of
C        boundary nodes, if not it is add to the list ..
C
C              write(6,*)'candidate node is :',inode
               CALL Add_to_List( IWKSP , NBNODE(JCOLOR) , INODE )
               IF( NBNODE(JCOLOR) .EQ. LIWK )THEN
                  WRITE(NOUT,888)
                  CALL EXIT(1)
               ENDIF
    8       CONTINUE
C
    7    CONTINUE
C
         WRITE(NOUT,200)JCOLOR,N,NBNODE(JCOLOR)
C
         WRITE(3,*)'ZONE F=FEPOINT,ET=',elementtype,',N=',NBNODE
     +   (JCOLOR),      ',E=',N
C
C     .. Write nodes ..
C
         DO 21 IPOIN = 1,NBNODE(JCOLOR)
ccc         IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(NOUT,111)
               INODE = IWKSP(IPOIN)
               WRITE(3,*)(POINT(J,INODE),J=1,DIM),(VSCL(J,INODE),J=1,
     +         NOFVAR),FLOAT(JCOLOR)
   21    CONTINUE
         WRITE(3,*)
ccc      WRITE(6,175)
C
C	writing the elements
C
         IFREQ = MAX0( N/ 20,  1)
ccc      WRITE(6,170)
         DO 25 IFACE = IFRST(JCOLOR) , ILAST(JCOLOR)
            IGAOK = IRANK(IFACE)
            IELEM = IBNDFAC(1,IGAOK)
            IVERT = IBNDFAC(2,IGAOK)
            DO 28 I =  1, DIM
               INODE = ICELNOD( ICYCL( I + IVERT , NOFVERT ) , IELEM )
               ICNVT(I) = Binsrch( IWKSP , NBNODE(JCOLOR) , INODE ,
     +         Lo , Me , Up )
   28       CONTINUE
ccc         IF((IELEM/IFREQ)*IFREQ .EQ. IELEM)WRITE(NOUT,111)
            IF(DIM.EQ.3)THEN
               WRITE(3,*)(ICNVT(J),J=1,DIM)
            ELSE
               WRITE(3,*)(ICNVT(J),J=1,DIM),ICNVT(1)
            ENDIF
   25    CONTINUE
ccc      WRITE(6,175)
C
    5 CONTINUE ! End of the outermost loop on colors
      CLOSE(3)
      WRITE(6,185)
C
C
*     WRITE(6,180)
C
      RETURN
 
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
  185 FORMAT(/,5X,'Tecplot file WRITTEN .......',/)
  200 FORMAT(5X,'BOUNDARY COLORED ',I2,' BND. FACES : ',I6,
     &' BND. POINTS ',I6)
  500 FORMAT(/10X,'Inconsistent data : IBNDFAC(',I1,',',I6,') = ',I6)
  510 FORMAT(/10X,'Incorrect boundary type ',I1,
     1' has been assigned to colour ',I2)
  888 FORMAT(10X,'Insufficient memory allocated ... sorry')
 9996 FORMAT(10X,A6,' PASSED   IFAIL = ',I2)
 9999 FORMAT(10X,A6,' RETURNED IFAIL = ',I2)
*
      END
