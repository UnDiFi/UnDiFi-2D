      PROGRAM DATRIANGLE
C
C     $Id: main.f,v 1.1 2013/01/25 16:01:33 abonfi Exp abonfi $
C
C     this program reads files in EulFS format
C     then writes two input files for triangle: *.node *.poly
C
C     28.11.2012: now asks for the colour of the hole (if any)
C     and tries to compute the coordinates of a point
C     inside the hole; it will NOT (or may not) work if the hole
C     is NOT a convex polygon
C
      IMPLICIT NONE
C
C
C
C       NIN     is the OUTPUT device number
C       NOUT    is the INPUT device number
C
C
C
C       NPOIN   number of nodes in the mesh
C       NBFAC:  number of boundary faces in the mesh
C       NELEM:  number of elements in the mesh
C       NHOLE:  number of bodies(holes) in the mesh
C
C
C       NDIM     is the space dimension
C       NOFVERT = NDIM+1 is the number of vertices
C
C
C
C
C
C
C
C
C
C     .. Parameters ..
      REAL*8 ZERO,HALF,ONE,TWO
      PARAMETER (ZERO=0.00d0,HALF=0.5d0,ONE=1.00d0,TWO=2.00d0)
      INTEGER*4 NIN,NOUT,NMAXBODIES
      PARAMETER (NIN=5,NOUT=6,NMAXBODIES=20)
      INTEGER NVA
      PARAMETER (NVA=33990000)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(NVA)
      INTEGER IXDRS(10)
      CHARACTER FILENAME(5)*80,FNAME*5,EXECMD*80
      CHARACTER*255 DATADIR
C     ..
C     .. Local Scalars ..
      INTEGER NDIM,I,LBNDFAC,LCELCEL,LCELNOD,LCORG,NPNOD,
     +        LNODCOD,LZROE,NBFAC,NBKGRD,LWKSP(3),
     +        NBPOIN,NELEM,NHOLE,NOFVAR,NOFVERT,NPOIN
C     ..
C     .. Local Arrays ..
      INTEGER ISTAK(1),ICOLOR(NMAXBODIES),ICHOLES(NMAXBODIES)
C
      INTEGER kspace,lenstr
C     ..
C     .. External Functions ..
      INTEGER INITXDR,ISTKGT
      EXTERNAL INITXDR,ISTKGT
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DINIT,
     +         IINIT,ISTKIN,ISTKRL,NODCOD,READAT,WKSPCE
C     ..
C     .. Common blocks ..
      COMMON /CSTAK/DSTAK
      COMMON /ES/IXDRS,FILENAME
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
      write(6,*)
      write(6,*)
     &'This code converts the EulFS code files into triangle files'
      write(6,*)
      NOFVAR = 4
C     OPEN (1,FILE='inp',STATUS='UNKNOWN')
C     READ (1,'(12X,A)') DATADIR
      DATADIR = "./"
      kspace = lenstr(datadir)
!     write(6,*)datadir(1:kspace)
      filename(2) = "file010.dat"
C     READ (1,'(12X,A)') filename(2)
C
C     filename(5) is the triangle node file to be overwritten
C
      write(6,*)'Enter triangle file'
      READ (5,'(A)') filename(5)
      write(6,*)filename(2)
C
      filename(1) = datadir(1:KSPACE) // 'file001.dat'
      filename(3) = datadir(1:KSPACE) // 'file002.dat'
      filename(4) = datadir(1:KSPACE) // 'file004.dat'
C
C     write(6,*)filename(1)
C     write(6,*)filename(2)
C     write(6,*)filename(3)
C     write(6,*)filename(4)
C
!     write(6,*)'before wkspce'
      CALL WKSPCE(NDIM,NPOIN,NPNOD,NELEM,NBFAC,NHOLE,NOFVAR,NOFVERT)
!     NBKGRD = 6399-141
!     write(6,*)'NBKGRD = ',NBKGRD
!     write(6,*)'beyond wkspce'
!     IF(NHOLE.GT.0)THEN
!        DO I = 1,NHOLE
!           WRITE(6,FMT="('Enter the colour of hole # ',I2,' ---> ',$)")
!    &I
!           READ(5,*)ICHOLES(I)
!        ENDDO
!     ENDIF
C
      CALL ISTKIN(NVA,4)
C
C ---------- ALLOCATE SPACE
C
      WRITE (NOUT,FMT=4000)

 4000 FORMAT (/,/,' MEMORY ALLOCATION   ',/,' ',19 ('='),/)
C
C     Nodal coordinates
C
      LCORG = ISTKGT((NPOIN+NPNOD)*NDIM,4)
      CALL DINIT((NPOIN+NPNOD)*NDIM,ZERO,DSTAK(LCORG),1)
C
C     Roe variables
C
      LZROE = ISTKGT(NOFVAR*(NPOIN+NPNOD),4)
      CALL DINIT(NOFVAR*(NPOIN+NPNOD),ZERO,DSTAK(LZROE),1)

C     Cell to node pointer
C
      LCELNOD = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELNOD),1)
C     Cell to node pointer
C
      LNODCOD = ISTKGT(NPOIN+NPNOD,2)
      CALL IINIT(NPOIN+NPNOD,0,ISTAK(LNODCOD),1)
C
C     Boundary structure
C
      LBNDFAC = ISTKGT(3*NBFAC,2)
      CALL IINIT(3*NBFAC,0,ISTAK(LBNDFAC),1)
C
C     Cell to cell pointer
C
      LCELCEL = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELCEL),1)
C
C      Reading data
C
!     write(6,*)'before readat'
      CALL READAT(DSTAK(LCORG),ISTAK(LCELNOD),ISTAK(LCELCEL),
     +            ISTAK(LBNDFAC),DSTAK(LZROE),NDIM,NOFVAR,NOFVERT,
     +            NPOIN+NPNOD,NELEM,NBFAC)
C
!     write(6,*)'beyond readat'
      CALL NODCOD(ISTAK(LNODCOD),NPOIN+NPNOD,NBPOIN,ISTAK(LCELNOD),
     +                  NOFVERT,NELEM,ISTAK(LBNDFAC),NBFAC)
C
      write(6,*)'NPNOD =',NPNOD
      IF(NPNOD.GT.0)THEN
         LWKSP(1) = ISTKGT(NOFVAR*(NPOIN+NPNOD),4)
         LWKSP(2) = ISTKGT(NOFVAR*(NPOIN+NPNOD),2) ! integer
         CALL DINIT(NOFVAR*(NPOIN+NPNOD),ZERO,DSTAK(LWKSP(1)),1)
         CALL IINIT(NOFVAR*(NPOIN+NPNOD),0,ISTAK(LWKSP(2)),1)
         CALL renumber(DSTAK(LCORG),ndim,npoin+npnod,DSTAK(LZROE),
     &              DSTAK(Lwksp(1)),ISTAK(LWKSP(2)),nOFVAR,
     &              ISTAK(LCELNOD),ISTAK(LCELNOD),NOFVERT,NELEM)
         CALL ISTKRL(2)
      ENDIF
C
      kspace = lenstr(filename(5))
      EXECMD = "mv "//filename(5)(1:KSPACE)//".node "
     &//filename(5)(1:KSPACE)//".node.BAK"
      WRITE(6,*)EXECMD
      CALL SYSTEM(EXECMD) 
C
      FNAME(1:KSPACE) = filename(5)(1:KSPACE)
      CALL WTRI(ISTAK(LBNDFAC),NBFAC,ISTAK(LCELNOD),NOFVERT,
     +DSTAK(LCORG),NDIM,DSTAK(LZROE),NOFVAR,ISTAK(LNODCOD),NPOIN+NPNOD,
     &FILENAME(5),ICHOLES,NHOLE)
C
      CALL EXIT(0)
C
  100 FORMAT (12X,A)
  120 FORMAT (12X,I6)
      END
