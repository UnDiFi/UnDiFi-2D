      PROGRAM TRIANGLEDAT
C
C     $Id:$
C
C     converts a triangle MESH to a dat file fmt
C     the input file for triangle can be created
C     using dat2triangle 
C     this program needs edge neigh and poly files
C     hence run
C     triangle -n -e -p filename.poly 
C
C     28.11.2012: now deals with periodic meshes
C
C     after the triangle filename
C     the user is now asked to enter new data:
C     2) enter y/Y if the mesh is periodic; otherwise enter n/N 
C     3-4) enter the two (different) colours of the two periodic surfaces
C     5) enter x/X if the corresponding nodes on the two periodic patches have
C        the same x coordinate; otherwise enter y/Y (have the same y coordinate) 
C     6) enter the number of holes
C     it is no longer required to give the number of holes;
C     this information is taken from the poly file
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NELEM,NPOIN,ns,nfr,nt,K,NHOLE,IXY,mdim,nbc,i
C
      INTEGER NVA
      PARAMETER (NVA=39800000)
C
      DOUBLE PRECISION DSTAK(NVA)
      INTEGER ISTAK(1)
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      COMMON /CSTAK/DSTAK
C
      INTEGER LCORG,LCELCEL,LCELNOD,LBNDPTR,LOCIA,LOCJA,LOCZ,
     &LOCW,NITEMS
C
      INTEGER   NOFVERT,NDIM,NDOF
      PARAMETER(NDIM=2,NOFVERT=(NDIM+1))
C
      INTEGER ifail,IXDRS,N,LWORK,ILEN(10),NPNOD,ID,ICLR(2)
C
      INTEGER  ISTKGT,INITXDR,IXDRINT,IXDRDMAT,IXDRIMAT,IXDRCLOSE,LENSTR
      EXTERNAL ISTKGT,INITXDR,IXDRINT,IXDRDMAT,IXDRIMAT,IXDRCLOSE,LENSTR
C
      DOUBLE PRECISION EPS
      COMMON /EPSCOM/EPS
      CHARACTER FNAME(10)*255
      CHARACTER ANSW*1
C
      equivalence(nt,nelem)
      equivalence(ns,npoin)
      equivalence(nfr,nbfac)
C
C     .. Executable Statements ..
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'This code converts triangle files into EulFS files '
      WRITE(6,*)
      WRITE(6,*)
      FNAME(10) = "file.1"
      WRITE(6,*)'Enter fname'
      READ(5,*)FNAME(10)
      K = LENSTR(FNAME(10))
      FNAME(1)(1:K+5) = FNAME(10)(1:K)//".node"
      FNAME(2)(1:K+4) = FNAME(10)(1:K)//".ele"
      FNAME(3)(1:K+6) = FNAME(10)(1:K)//".neigh"
      FNAME(4)(1:K+5) = FNAME(10)(1:K)//".edge"
      FNAME(5)(1:K+5) = FNAME(10)(1:K)//".poly"
      ILEN(1) = K+5
      ILEN(2) = K+4
      ILEN(3) = K+6
      ILEN(4) = K+5
      ILEN(5) = K+5
      WRITE(6,*)'Trying to open ',FNAME(1)(1:K+5)
      open(1,file=FNAME(1)(1:K+5),status="old")
      rewind(1)
      read(1,*) ns,id,ndof,id
!     if(ndof.NE.4)then
!        write(6,*)'Nof attributes should be 4 in ',FNAME(1)(1:K+5),
!    &'while it is ',ndof
!        CALL EXIT(1)
!     endif
      if(id.NE.1)then
         write(6,*)'Nof bndry markers should be 1 in ',FNAME(1)(1:K+5),
     &'while it is ',id
         CALL EXIT(1)
      endif
      close(1)
      WRITE(6,*)'Trying to open ',FNAME(2)(1:K+4)
      open(1,file=FNAME(2)(1:K+4),status="old")
      rewind(1)
      read(1,*) nt
      close(1)
      WRITE(6,*)'Trying to open ',FNAME(4)(1:K+5)
      open(1,file=FNAME(4)(1:K+5),status="old")
      rewind(1)
      read(1,*) nfr ! this is actually the nof edges in the whole grid
      close(1)
C     reading poly file
      WRITE(6,*)'Trying to open ',FNAME(5)(1:K+5)
      open(1,file=FNAME(5)(1:K+5),status="old")
      rewind(1)
      read(1,*)nitems,mdim,ndof,nbc
      do i = 1,nitems
         read(1,*)
      enddo
      read(1,*)nitems,nbc
      do i = 1,nitems
         read(1,*)
      enddo
      read(1,*)nhole
      close(1)
C
      write(6,50)ns,nt,nfr,nhole
   50 format(5X,'NPOIN = ',I8,5X,'NELEM ',I8,5X,'NFACE = ',I8,' NHOLE = 
     &',I2)
C
C     Set the length of the stack
C
      CALL ISTKIN(NVA,4)
C
C ---------- ALLOCATE SPACE
C
      WRITE (6,4000)

 4000 FORMAT (//' MEMORY ALLOCATION   '/' ',19 ('=')/)
C
C     Nodal coordinates ..
C
      LCORG = ISTKGT(NPOIN*NDIM,4)
      CALL DINIT(NPOIN*NDIM,0.d0,DSTAK(LCORG),1)
C
      LOCZ = ISTKGT(NDOF*NPOIN,4)
      CALL DINIT(NPOIN*NDOF,0.d0,DSTAK(LOCZ),1)
C
C     Cell to node pointer ..
C
      LCELNOD = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELNOD),1)
C
C     Cell to cell pointer ..
C
      LCELCEL = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELCEL),1)
C
C     Boundary structure ..
C
      LBNDPTR = ISTKGT(3*NBFAC,2)
      CALL IINIT(3*NBFAC,0,ISTAK(LBNDPTR),1)
C
C     Reading data from the MESH file
C
      WRITE(6,*)'Reading mesh..........'
      CALL RTRI(DSTAK(LCORG),NPOIN,DSTAK(LOCZ),NDOF,ISTAK(LCELNOD),
     &ISTAK(LCELCEL),ISTAK(LBNDPTR),NELEM,NBFAC,FNAME,ILEN)
      WRITE(6,*)'Done'
C
C 
C
      WRITE(6,FMT="('Are there any periodic surfaces? y/n ',$)")
      READ(5,*)ANSW
C
      IF(ANSW.EQ.'y'.OR.ANSW.EQ.'Y')THEN
         NITEMS = 2*NPOIN
         LOCW = ISTKGT(NITEMS,2)
         WRITE(6,FMT=120)
         READ(5,*)ICLR(1)
         WRITE(6,FMT=125)
         READ(5,*)ICLR(2)
    8    WRITE(6,FMT=135)
         READ(5,*)ANSW
         IF(ANSW.EQ.'x'.OR.ANSW.EQ.'X')THEN
            IXY = 1
         ELSEIF(ANSW.EQ.'y'.OR.ANSW.EQ.'Y')THEN
            IXY = 2
         ELSE
             GOTO 8
         ENDIF
  120 FORMAT('Enter the colour of the 1st periodic surface : ',$)
  125 FORMAT('Enter the colour of the 2nd periodic surface : ',$)
  135 FORMAT('The periodic surfaces have the same x/y ? : ',$)
         CALL xa112d(istak(lcelnod),istak(lcelcel),nelem,istak(lbndptr),
     +nbfac,dstak(lcorg),istak(locw),nITEMS,ndim,npoin,iclr,IXY)
         write(6,*)nitems
         LOCIA = ISTKGT(NITEMS,2)
         LWORK = ISTKGT(NITEMS,4)
      CALL periodic(dstak(lcorg),DSTAK(LWORK),istak(lcelnod),
     &istak(lcelcel),istak(locia),ndim,nofvert,
     +npoin,nelem,nbfac,istak(locw),npnod,dstak(locz),NDOF,IXY)
         CALL ISTKRL(3)
      ELSE
         NPNOD = 0
      ENDIF
      write(6,*)'There are ',npnod,' periodic nodes' 
      write(6,*)'NPOIN is currently set to ',npoin
!     write(6,*)'Enter the number of holes '
!     read(5,*)NHOLE
!     NHOLE = 0
C
CCCCCCC#endif
C
C    Writing the mesh & connectivity file
C
      ixdrs = INITXDR( 'file001.dat' , 'w' ,.FALSE.)
      IFAIL = IXDRINT( ixdrs , NDIM )
      IFAIL = IXDRINT( ixdrs , NPOIN-NPNOD )
      IFAIL = IXDRINT( ixdrs , NELEM )
      IFAIL = IXDRINT( ixdrs , NBFAC )
      IFAIL = IXDRINT( ixdrs , NHOLE )
      IFAIL = IXDRDMAT( ixdrs , NDIM*NPOIN , DSTAK(LCORG) )
      IFAIL = IXDRIMAT( ixdrs , NOFVERT*NELEM , ISTAK(LCELNOD) )
      IFAIL = IXDRIMAT( ixdrs , 3*NBFAC , ISTAK(LBNDPTR) )
      IFAIL = IXDRCLOSE( ixdrs )
C
C
C    Writing the neighbouring elements file
C
      ixdrs = INITXDR( 'file002.dat' , 'w' ,.FALSE.)
      IFAIL = IXDRINT( ixdrs , NDIM )
      IFAIL = IXDRINT( ixdrs , NELEM )
      IFAIL = IXDRIMAT( ixdrs , NOFVERT*NELEM , ISTAK(LCELCEL) )
      IFAIL = IXDRCLOSE( ixdrs )
C
      ixdrs = INITXDR('file003.dat','w',.FALSE.)
      IFAIL = IXDRINT(ixdrs,NPOIN)
      IFAIL = IXDRINT(ixdrs,NDOF)
      IFAIL = IXDRDMAT(ixdrs,NDOF*NPOIN,DSTAK(LOCZ))
      IFAIL = IXDRCLOSE(ixdrs)
C
C
C
      CALL EXIT(0)
C
  111 FORMAT ('.',$)
  112 FORMAT (/,5X,A11,1X,'BACKUP FILE WRITTEN .......',/)
  113 FORMAT (/,5X,A11,1X,'BACKUP FILE OPEN ..........',/)
  114 FORMAT (/,5X,A11,1X,'BACKUP FILE READ ..........',/)
  115 FORMAT (' done')
  300 FORMAT (/,10X,'IFORM = ',I2,' while it should be 1 or 2')

      END
