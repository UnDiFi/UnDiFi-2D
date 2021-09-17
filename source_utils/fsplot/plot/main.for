      PROGRAM PLOT
C
C     $Id: main.f,v 1.9 2020/02/05 10:46:08 abonfi Exp abonfi $
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
C
C
      INCLUDE 'IO'
      INCLUDE 'es.h'
      INCLUDE 'mesh_i4'
      INCLUDE 'int_flags'
      INCLUDE 'dim_flags'
C
      INTEGER MAX_NUM_RAKES
      PARAMETER(MAX_NUM_RAKES=15)
      INTEGER NVA
      PARAMETER (NVA=267000000)
!     PARAMETER (NVA=9900000)
C
      DOUBLE PRECISION DSTAK(NVA)
      INTEGER ISTAK(1)
      REAL RSTAK(1)
      EQUIVALENCE (DSTAK(1),ISTAK(1),RSTAK(1))
      COMMON /CSTAK/DSTAK
C
      INTEGER      LCORG,LZROE,LCELNOD,LBNDFAC,LCELCEL,LHEAT,LQ2,
     &          LNODPTR,LBNDPTR,LNODCOD
      COMMON /NLOC/LCORG,LZROE,LCELNOD,LBNDFAC,LCELCEL,LHEAT,LQ2,
     &          LNODPTR,LBNDPTR,LNODCOD
C
      INTEGER ADIM
      DOUBLE PRECISION C1,C2,C3,C4,fren,ADIM
      common /freestream/ fren
      common /visclaw/ C1,C2,C3,C4,ADIM
C
      INTEGER NBFAC,ifail,ICLR,LRANK,LENW,LSKIN,LWORK(10),LTINDX,
     >I,iunit,IRAKE,NRAKES,it,nlen,LENW2,IDUMMY,LWKSP(3),
     +NELEMTYPES,NNODESFR,NDOF,NITEMS,IWRK
      EQUIVALENCE(nBoundaryFaces,NBFAC)
      INTEGER IV(3,5) 
      DOUBLE PRECISION LINEBGN(3,MAX_NUM_RAKES),LINEEND(3,MAX_NUM_RAKES)
      DOUBLE PRECISION AA,BB,CC,REYNO,M_INFTY,S,T,T_REF,T0_REF
      DOUBLE PRECISION angle_min, angle_max, retval
      DOUBLE PRECISION alpha_min, alpha_ave, alpha_area, 
     &triangulation_area
      INTEGER angle_max_triangle , angle_min_triangle
      LOGICAL DPLOT_,TECPLOT_,NORMAL_TO_THE_WALL,LFLAG,KAN(3)
C
C     KAN Kind of ANalysis
C
C     KAN(1) = .TRUE. KAN(2:3) = .FALSE. -----> compressible
C     KAN(2) = .TRUE. KAN(1,3) = .FALSE. -----> Incompressible
C     KAN(3) = .TRUE. KAN(1:2) = .FALSE. -----> scalar
C
C     ADIM non-dimensionalization for External flows = 0
C     ADIM non-dimensionalization for Internal flows = 1
C
      CHARACTER*2 answ
      CHARACTER*10 fname
      CHARACTER*80 title,probefname
      CHARACTER*255 datadir
      LOGICAL lprobe
C
      INTEGER  ISTKGT,ISTKST,ISTKQU,INITXDR,IXDRCLOSE
      EXTERNAL ISTKGT,ISTKST,ISTKQU,INITXDR,IXDRCLOSE
C     ..
C     .. External Subroutines ..
      INTEGER  IXDRIMAT,IXDRINT
      EXTERNAL IXDRIMAT,IXDRINT
      DOUBLE PRECISION DNRM2,DASUM
      EXTERNAL DNRM2,DASUM
C
      INTEGER kspace,lenstr
C
C     .. Executable Statements ..
C
      KAN(1) = .FALSE.
      KAN(2) = .FALSE.
      KAN(3) = .FALSE.
C
      OPEN (1,file='inp',STATUS='OLD')
      READ (1,'(12X,A)') DATADIR
      kspace = lenstr(datadir)
      READ (1,'(12X,A)') filename(2)
      READ (1,'(12X,A2)')answ
      ADIM = -1
      IF    (answ(1:1).EQ."C".OR.answ.EQ."c")THEN
           KAN(1) = .TRUE.
           IF(answ(2:2).EQ."E".OR.answ(2:2).EQ."e")THEN
              ADIM = 0 ! External flows
           ELSEIF(answ(2:2).EQ."I".OR.answ(2:2).EQ."i")THEN
              ADIM = 1 ! Internal flows
           ELSE
              WRITE(6,*)'Need to specify either CI(ci) or CE(ce)' 
              CALL EXIT(1)
           ENDIF
      ELSEIF(answ(1:1).EQ."I".OR.answ.EQ."i")THEN
           KAN(2) = .TRUE.
      ELSEIF(answ(1:1).EQ."S".OR.answ.EQ."s")THEN
           KAN(3) = .TRUE.
      ELSE
           WRITE(6,*)'Need to specify whether I/C/S !'
           CALL EXIT(1)
      ENDIF
      READ (1,'(12X,F12.8)') T_REF
      READ (1,'(12X,F12.8)') T0_REF
      DPLOT_ = .TRUE.
      TECPLOT_ = .TRUE.
      READ (1,'(12X,F12.8)') REYNO
      READ (1,'(12X,F12.8)') M_INFTY
      READ (1,'(12X,I3)') ICLR
      READ (1,'(12X,F12.8)') AA
      READ (1,'(12X,F12.8)') BB
      READ (1,'(12X,F12.8)') CC
      READ (1,'(12X,I3)') I
      NORMAL_TO_THE_WALL = (I .NE. 0)
C
C     Check the nondimensionalization when dealing with compressible
C     flows
C
      IF(KAN(1))THEN
         IF(ADIM.EQ.0)THEN ! External flows
            C2 = M_INFTY*M_INFTY 
            C1 = M_INFTY*C2
            C3 = S0/T_REF
            C4 = ONE+C3
         ELSEIF(ADIM.EQ.1)THEN ! Internal flows
            C1 = ONE/GAM**1.5d0
            C2 = ONE/GAM
            C3 = S0/T0_REF
            C4 = ONE+C3
         ELSE  
         ENDIF
      ENDIF
c
c     read 3D rakes
c
      READ (1,'(12X,I3)') nrakes
      IF(nrakes.GT.MAX_NUM_RAKES)THEN
         stop 'far too many rakes! increase MAX_NUM_RAKES'
      ENDIF
      DO 178 irake = 1,  nrakes
         READ (1,'(11X,6(1X,F7.3))') 
     +   (linebgn(i,irake),i=1,3),(lineend(i,irake),i=1,3)
!        write (6,'(11X,6(1X,F7.3))') 
!    +   (linebgn(i,irake),i=1,3),(lineend(i,irake),i=1,3)
  178 CONTINUE
      READ (1,'(12X,L)') lprobe
      IF(LPROBE)THEN
         READ (1,'(12X,A)') probefname
      ENDIF
      CLOSE (1)
C
      filename(1) = datadir(1:KSPACE) // 'file001.dat'
      filename(3) = datadir(1:KSPACE) // 'file002.dat'
      filename(4) = datadir(1:KSPACE) // 'file004.dat'
C
      fren = (2.d0 + 0.4d0*M_INFTY*M_INFTY*1.4d0)/
     +       (1.4d0*0.8d0*M_INFTY*M_INFTY)
C
C
      CALL ISTKIN(NVA,4)
C
C
      CALL WKSPCE
C
C ---------- ALLOCATE SPACE
C
      WRITE (NOUT,4000)

 4000 FORMAT (//' MEMORY ALLOCATION   '/' ',19 ('=')/)
C
C     .. Nodal coordinates ..
C
      LCORG = ISTKGT((NPOIN+NPNOD)*DIM,4)
      CALL DINIT((NPOIN+NPNOD)*DIM,ZERO,DSTAK(LCORG),1)
C
C     .. Nodal values ..
C
      LZROE = ISTKGT((NPOIN+NPNOD)*NOFVAR,4)
      CALL DINIT((NPOIN+NPNOD)*NOFVAR,ZERO,DSTAK(LZROE),1)
C
C     .. Cell to node pointer ..
C
      LCELNOD = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELNOD),1)
C
C     .. Boundary structure ..
C
      LBNDFAC = ISTKGT(3*NBFAC,2)
      CALL IINIT(3*NBFAC,0,ISTAK(LBNDFAC),1)
C
C     ... Reading data ...
C
      CALL READAT(DSTAK(LCORG),ISTAK(LCELNOD),ISTAK(LBNDFAC),
     +            DSTAK(LZROE))
C
      WRITE(6,*)' Computing L2 norm with stride '
      DO I = 0,NOFVAR-1
         IDUMMY = LZROE + I
         S=DNRM2((NPOIN+NPNOD),DSTAK(IDUMMY),NOFVAR)
         T=DASUM((NPOIN+NPNOD),DSTAK(IDUMMY),NOFVAR)/(NPOIN+NPNOD)
         WRITE(6,*)' L2 norm of dof ',I,' is = ',S,' <Q> = ',T
      ENDDO
C
C rotate the mesh about a pole
C
!     CALL rotate2d("file110.dat",DSTAK(LCORG),dim,ISTAK(LCELNOD),
!    +nofvert,nelem,nhole,DSTAK(LZROE),nofvar,npoin,npnod,
!    +ISTAK(LBNDFAC),nbfac)
C
C     ... Dplot ...
C
      IF (DIM.EQ.3) GOTO 100
      IF (.FALSE.) THEN

          IF (NOFVAR.EQ.4) CALL PARM_TO_CONS(DSTAK(LZROE))

          IF (NOFVAR.EQ.3) then
           nofvar = 1
           Lwksp = ISTKGT((NPOIN+NPNOD),4)
           CALL Dcopy((NPOIN+NPNOD),DSTAK(LZROE),3,dstak(lwksp),1)
          CALL DPLOT('p.dpl',DSTAK(LCORG),ISTAK(LCELNOD),
     +               DSTAK(lwksp),ISTAK(LBNDFAC))
           CALL Dcopy((NPOIN+NPNOD),DSTAK(LZROE+1),3,dstak(lwksp),1)
          CALL DPLOT('u.dpl',DSTAK(LCORG),ISTAK(LCELNOD),
     +               DSTAK(lwksp),ISTAK(LBNDFAC))
           CALL Dcopy((NPOIN+NPNOD),DSTAK(LZROE+2),3,dstak(lwksp),1)
          CALL DPLOT('v.dpl',DSTAK(LCORG),ISTAK(LCELNOD),
     +               DSTAK(Lwksp),ISTAK(LBNDFAC))
           nofvar = 3
          else
          CALL DPLOT('file012.dpl',DSTAK(LCORG),ISTAK(LCELNOD),
     +               DSTAK(LZROE),ISTAK(LBNDFAC))
          endif
          IF (NOFVAR.EQ.4) CALL CONS_TO_PARM(DSTAK(LZROE))
      ENDIF
C
  100 CONTINUE
C     ... Tecplot ...
C
caldo     CALL GMV('file012.inp',DSTAK(LCORG),DIM,NPOIN,
caldo+             ISTAK(LCELNOD),NOFVERT,NELEM,DSTAK(LZROE),NOFVAR)
C     ... Tecplot ...
C
      IF (TECPLOT_) THEN
          fname = "tecXXX.dat"
          write(fname(4:6),FMT="(I3.3)")it
          IF(NPOIN.GE.1000000)THEN
             WRITE(6,*)'This is a large dataset with ',NPOIN,
     &       ' gridpoints; do u want to produce a volume file?'
             READ(5,*)answ
          ELSE
             answ = 'Y'
          ENDIF
             IF(answ.EQ.'y'.OR.answ.EQ.'Y')tHEN
             IWRK=(MAX(DIM,NOFVAR)+NPOIN+NPNOD)/2
caldo        LWKSP(1) = ISTKGT(IWRK,2)
caldo     CALL BTEC('file012.plt',DSTAK(LCORG),DIM,NPOIN,
caldo&                 ISTAK(LCELNOD),NOFVERT,NELEM,
caldo+                 DSTAK(LZROE),NOFVAR,ISTAK(LWKSP(1)),IWRK)
caldo        CALL ISTKRL(1)
          CALL TECPLOT('file012.dat',DSTAK(LCORG),DIM,NPOIN,
     &                 ISTAK(LCELNOD),NOFVERT,NELEM,
     +                 DSTAK(LZROE),NOFVAR)
             ENDIF
          NLEN = ISTKQU(2)
          nlen = 0.45*nlen
          LWKSP(1) = ISTKGT(Nlen,2)
          LRANK = ISTKGT(Nlen,2)
          CALL BNDPLOT('file020.dat',DSTAK(LCORG),DSTAK(LZROE),
     +                 ISTAK(LCELNOD),ISTAK(LBNDFAC),ISTAK(LWKSP(1)),
     +                 ISTAK(LRANK),Nlen)
          CALL ISTKRL(2)
      ENDIF
C
      INQUIRE(FILE="file016.dat",EXIST=LFLAG)
      IF(.NOT.LFLAG)GOTO 77
      IF(DIM.EQ.2 .AND. NOFVAR .GT. 1)THEN
C
C     .. Cell to node pointer ..
C
!         LNODCOD = ISTKGT(NPOIN+NPNOD,2)
!         CALL IINIT(NPOIN+NPNOD,0,ISTAK(LNODCOD),1)
!         CALL NODCOD(ISTAK(LNODCOD),NPOIN+NPNOD,NBPOIN,
!    &                ISTAK(LCELNOD),NOFVERT,NELEM,ISTAK(LBNDFAC),NBFAC)
!         CALL SetBndryNodePtr(LBNDFAC,LNODCOD,LNODPTR,LCELNOD,NBFAC,
!    &                         NPOIN+NPNOD,NBPOIN)
!         LBNDPTR = ISTKGT(6*NBFAC,2)
!         CALL IINIT(6*NBFAC,0,ISTAK(LBNDPTR),1)
!         CALL SetBndryPntr(ISTAK(LCELNOD),ISTAK(LCELCEL),NOFVERT,NELEM,
!    &    ISTAK(LBNDFAC),NBFAC,ISTAK(LBNDPTR),ISTAK(LNODPTR),NBPOIN)
C
          LENW = NBFAC
          LSKIN = ISTKGT(LENW,4)
          LHEAT = ISTKGT(LENW,4)
          LWORK(1) = ISTKGT(2*LENW,2)
          LTINDX = ISTKGT(LENW,4)
          LWKSP(1) = ISTKGT(2*LENW,2)
          LRANK = ISTKGT(LENW,2)
          CALL CF(DSTAK(LSKIN),DSTAK(LHEAT),DSTAK(LTINDX),
     +            ISTAK(LWKSP(1)),ISTAK(LWORK(1)),ISTAK(LRANK),LENW,
     +            ISTAK(LCELNOD),DSTAK(LCORG),DSTAK(LZROE),DIM,NOFVERT,
     +            NOFVAR,REYNO,M_INFTY,'file016.dat')
          CALL ISTKRL(3)
      ELSEIF(DIM.EQ.3 .AND. NOFVAR .GT. 1)THEN
          LENW = NBFAC
          LSKIN = ISTKGT(LENW,4)
          LHEAT = ISTKGT(LENW,4)
          LWORK(1) = ISTKGT(2*LENW,2)
          LTINDX = ISTKGT(LENW,4)
          LWKSP(1) = ISTKGT((NOFVERT-1)*LENW,2)
          LWORK(5) = ISTKGT((NOFVERT-1)*LENW,2)
          NDOF = 5
          LWORK(2) = ISTKGT(NDOF*LENW,4)
          LWORK(3) = ISTKGT(LENW,4)
          LWORK(4) = ISTKGT(LENW,4)
          LWORK(6) = ISTKGT(LENW,4)
          LWORK(7) = ISTKGT(LENW,4)
          CALL CH3D(DSTAK(LSKIN),DSTAK(LHEAT),DSTAK(LTINDX),
     &         ISTAK(LWKSP(1)),ISTAK(LWORK(1)),ISTAK(LWORK(5)),
     &         DSTAK(LWORK(2)),NDOF,DSTAK(LWORK(3)),DSTAK(LWORK(4)),
     &         DSTAK(LWORK(6)),DSTAK(LWORK(7)),LENW,ISTAK(LCELNOD),
     +         DSTAK(LCORG),DIM,NOFVERT,NOFVAR,
     &    REYNO,M_INFTY,'file016.dat')
          CALL ISTKRL(11)
      ENDIF
   77 CONTINUE 
C
C     Make room for the cell 2 cell pointer 
C
      IF (ICLR.GE.0.OR.LPROBE) THEN
C
         LCELCEL = ISTKGT(NELEM*NOFVERT,2)
         CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELCEL),1)
C
         IXDRS(3) = INITXDR(FILENAME(3),'r',.FALSE.)
         IFAIL = IXDRINT(IXDRS(3),IDUMMY)
         IFAIL = IXDRINT(IXDRS(3),IDUMMY)
         IFAIL = IXDRIMAT(IXDRS(3),NOFVERT*NELEM,ISTAK(LCELCEL))
         IFAIL = IXDRCLOSE(IXDRS(3))
C
      ENDIF
C
      IF (ICLR.GE.0) THEN
C
C     .. Cell to cell pointer ..
C
 
          IF(DIM.EQ.2)THEN
          CALL CUT(ICLR,AA,BB,CC,DSTAK(LZROE),NOFVAR,DSTAK(LCORG),DIM,
     +             NPOIN+NPNOD,ISTAK(LCELNOD),ISTAK(LCELCEL),NOFVERT,
     +             NELEM,ISTAK(LBNDFAC),NBFAC,ISTAK(LWORK(1)),
     +             DSTAK(LSKIN),LENW,REYNO,M_INFTY,NORMAL_TO_THE_WALL,
     +             FILENAME(3),KAN)
!         CALL BL(ICLR,AA,BB,CC,DSTAK(LZROE),NOFVAR,DSTAK(LCORG),DIM,
!    +             NPOIN+NPNOD,ISTAK(LCELNOD),ISTAK(LCELCEL),NOFVERT,
!    +             NELEM,ISTAK(LBNDFAC),NBFAC,ISTAK(LWORK),DSTAK(LSKIN),
!    +             LENW,REYNO,M_INFTY,NORMAL_TO_THE_WALL,FILENAME(3))
          ELSEIF(DIM.EQ.3)THEN
             IUNIT = 20
             OPEN(UNIT=IUNIT,FILE='rake3d.dat',FORM="FORMATTED")
             title = "this is a test"
             CALL WTecplotHeader(IUNIT,title,NOFVAR)
caldo
caldo
caldo
             DO 145 IRAKE = 1,NRAKES
c
c LWKSP(1) will keep the rake(s)
c leave room for the distance (hence NOFVAR+1)
c
!            LENW2 = 0.9*ISTKQU(2)/(NOFVAR+1)
             LENW2 = 1000
             LWKSP(1) = ISTKGT(LENW2*(NOFVAR+1),4)
             LWKSP(2)= ISTKGT(LENW2*DIM,4)
             CALL DINIT(LENW2*(NOFVAR+1),ZERO,DSTAK(LWKSP(1)),1)
             CALL DINIT(LENW2*DIM,ZERO,DSTAK(LWKSP(2)),1)
c
             CALL CUT3D(ICLR,IRAKE,LINEBGN(1,IRAKE),LINEEND(1,IRAKE),
     +       DSTAK(LZROE),DSTAK(LWKSP(1)),DSTAK(LWKSP(2)),LENW2,NOFVAR,
     +       DSTAK(LCORG),DIM,NPOIN+NPNOD,ISTAK(LCELNOD),ISTAK(LCELCEL),
     +       NOFVERT,NELEM,ISTAK(LBNDFAC),NBFAC,ISTAK(LWORK(1)),
     +       DSTAK(LSKIN),LENW,REYNO,NORMAL_TO_THE_WALL)
             CALL ISTKRL(1) ! can we free both?
  145        CONTINUE
             CLOSE(IUNIT)
c
          ENDIF ! ndim
 
      ENDIF
C
      IF (LPROBE) THEN
         IUNIT = 12
         OPEN(IUNIT,FILE=probefname,FORM="FORMATTED")
         READ(IUNIT,*)NITEMS
         WRITE(6,*)'Reading ',NITEMS,' probes from file ',probefname
         LWKSP(1) = ISTKGT(NITEMS*DIM,4) 
         LWKSP(2) = ISTKGT(NITEMS*NOFVAR,4) 
         LWKSP(3) = ISTKGT((NPOIN+NPNOD),2) 
         CALL RPROBE(IUNIT,DSTAK(LWKSP(1)),DIM,NITEMS)
         CLOSE(IUNIT)
         IV(2,1) = NPOIN
         IV(3,1) = NITEMS
         IV(2,2) = NPNOD
         IV(3,2) = 0
         IV(2,3) = NELEM
         IV(3,3) = 0
         CALL BURKARDT(NOFVAR,DIM,IV(1,1),IV(1,2),NOFVAR,
     &                 NOFVERT,IV(1,3),ISTAK(LCELNOD),ISTAK(LCELCEL),
     &                 DSTAK(LCORG),DSTAK(LWKSP(1)),DSTAK(LZROE),
     &                 DSTAK(LWKSP(2)),ISTAK(LWKSP(3)))
         CALL r8mat_print(NOFVAR,NITEMS,DSTAK(LWKSP(2)),'values')
         IUNIT = 12
         WRITE(6,*)'Writing probe data in probeout.dat'
         OPEN(IUNIT,FILE="probeout.dat",FORM="FORMATTED")
         CALL WPROBE(IUNIT,DSTAK(LWKSP(1)),DIM,NITEMS,DSTAK(LWKSP(2)),
     &               NOFVAR)
         CLOSE(IUNIT)
         CALL ISTKRL(3)
      ENDIF
caldo
      IF(DIM.EQ.2)THEN
         LWKSP(3) = ISTKGT(NELEM,4)
         call triangulation_areas ( NPOIN+NPNOD, DSTAK(LCORG), NOFVERT,
     &       NELEM, ISTAK(LCELNOD), DSTAK(LWKSP(3)), 
     &       triangulation_area )
         WRITE(6,*)'Triangulation area is ',triangulation_area
!        OPEN(77,FILE='vol.log',FORM="FORMATTED")
         IT = 0
         DO I = LWKSP(3),LWKSP(3)+NELEM-1
            IT = IT+1
            IF( DSTAK(I). LE. ZERO )THEN
                WRITE(6,*)'Triangle ',IT,' HAS NEGATIVE OR ZERO AREA = '
     &          ,DSTAK(I)
            ELSE
!               WRITE(77,*)'Triangle ',IT,' HAS AREA = ',DSTAK(I)
            ENDIF
         ENDDO
         CALL ISTKRL(1) ! release LWKSP(3)
!        CLOSE(77)
         CALL triangulation_order3_check ( NPOIN+NPNOD, NELEM,
     &         ISTAK(LCELNOD), IFAIL )
         IF( (IFAIL .EQ. 0) .OR. (IFAIL.EQ.8.AND.NHOLE.NE.0) )THEN
         WRITE(6,*)'triangulation_order3_check has returned ',IFAIL
         ELSE
             CALL EXIT(IFAIL)
         ENDIF
         CALL triangulation_order3_edge_check ( NELEM, ISTAK(LCELNOD),
     &        NNODESFR, IFAIL )
         WRITE(6,*)'triangulation_order3_edge_check has returned ',IFAIL
         IF( IFAIL .NE. 0 )CALL EXIT(IFAIL)
         call alpha_measure ( NPOIN+NPNOD, DSTAK(LCORG), NOFVERT, 
     &        NELEM, ISTAK(LCELNOD), alpha_min, alpha_ave, alpha_area )
         WRITE(6,*)'Alpha criterion: min, avg, weighted = ',alpha_min, 
     &             alpha_ave, alpha_area
         LQ2 = ISTKGT(NOFVERT*NELEM,2)
         call triangulation_neighbor_elements ( NOFVERT, NELEM,
     &        ISTAK(LCELNOD), ISTAK(LQ2) )
         CALL triangulation_delaunay_discrepancy_compute ( NPOIN+NPNOD,
     &        DSTAK(LCORG), NOFVERT, NELEM, ISTAK(LCELNOD), ISTAK(LQ2),
     &        angle_min, angle_min_triangle, angle_max,
     &        angle_max_triangle, retval )
        CALL ISTKRL(1) ! release LQ2
!
      write(6,*)
      write(6,*)'The Delaunay check has returned ',retval
      write(6,*)'The Delaunay check has returned min angle ',angle_min
      write(6,*)'in triangle ',angle_min_triangle
      write(6,*)'The Delaunay check has returned max angle ',angle_max
      write(6,*)'in triangle ',angle_max_triangle
      write(6,*)
      ENDIF ! on NDIM
C
      call mshsiz ( ISTAK(LCELNOD), NOFVERT, NELEM, DSTAK(LCORG),
     &       DIM)
C
      CALL ISTKRL(ISTKST(1))
C
      CALL EXIT(0)
C

      END
