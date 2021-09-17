      SUBROUTINE CHKBND(IBNDFAC,NBFAC,ICELFAC,ICELNOD,NOFVERT,NELEM,
     &FACENORM,COORD,NDIM,NFACE,NWFAC,NBODY4,NBODY6,A,JA,IA,LDA,NCL,
     2LFLAG,IWORK)
C
C     $Id: bndchk.f,v 1.17 2020/03/28 09:45:52 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'bnd.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.com'
      INCLUDE 'io.com'
      INCLUDE 'flags.com'
C
      INTEGER NBFAC,NOFVERT,NELEM,NDIM,NFACE,NWFAC,NBODY4,NBODY6,
     &LDA,NCL
      INTEGER NN
      PARAMETER(NN=NBTYPE+1)
      CHARACTER EXT1*2,EXT2*3,FILENAME*16,ERRMSG*72
C
C
      INTEGER IBNDFAC(3,*),ICELFAC(NOFVERT,*),ICELNOD(NOFVERT,*),
     &JA(LDA,*),IA(NCL+1),IWORK(NOFVERT,*)
      DOUBLE PRECISION FACENORM(NDIM,*),COORD(NDIM,*),A(LDA,*)
C
      INTEGER IBC,IL,I,J,K,JCOLOR,NERR,IFRST,N1,N2,IOPT
      INTEGER IELEM,IFREQ,IVERT,IBDFACE,JBGN,JEND,IFACE
      DOUBLE PRECISION TEMP,TEMP1,
     +AR,ARMIN,ARMAX,ARAVG,H,HMIN,HMAX,HAVG
C
C
      INTEGER          IWKSP(0:NBTYPE)
      DOUBLE PRECISION WKSP(3),BNDSURF(0:NBTYPE)
C
C
      INTEGER ICYCL
      EXTERNAL ICYCL,DAXPY
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
C
      LOGICAL LFLAG
C
      INTRINSIC MAX0
C
C
      DATA BNDSURF,TEMP,ARMIN,ARMAX,ARAVG,HMIN,HMAX,HAVG / NN*ZERO,
     +ZERO,1.E+38,2*ZERO,1.E+38,2*ZERO/
      DATA (IWKSP(J),J=0,NBTYPE),NERR/ NN*0,0 /
      DATA ERRMSG(1:7)/'BNDCHK '/
C
      NWFAC = 0
      IOPT = 1
C
      WRITE(NOUT,2000)MY_PE
2000  FORMAT(//' CHECKING BOUNDARIES ON PE #',I4,/' ',19('=')/)
C
      DO 14 J = 1 , NDIM
   14 WKSP(J) = ZERO
C
C     each boundary face has a color, which is stored in IBNDFAC(3,*);
C     the option "-color" assigns a boundary type (see bnd.h) to each color
C
C
      IFREQ = MAX0( 1 , NBFAC / 20 )
      WRITE(NOUT,114)
      DO 16 IBDFACE = 1 , NBFAC
         IF((IBDFACE/IFREQ)*IFREQ .EQ. IBDFACE)WRITE(NOUT,111)
         JCOLOR = IBNDFAC(3,IBDFACE)
C
C     colors can range between 0 and NCOLOR; if a color outside range
C     is found in the datafile, issue an error
C
         IF( JCOLOR .LT. 0 .OR. JCOLOR .GT. NCOLOR )THEN
            WRITE(ERRMSG(8:72),500)3,IBDFACE,JCOLOR
            NERR = 6 
            CALL SETERR(ERRMSG,72,NERR,IOPT)
         ENDIF
c
c     assign the boundary type to a certain color
c
         IBC = ICOLOR(JCOLOR,1)
c
c     check if that boundary type exists
c
         IF( IBC .LT. 0 .OR. IBC .GT. NBTYPE )THEN
            WRITE(ERRMSG(8:72),510)IBC,JCOLOR
            NERR = 7
            CALL SETERR(ERRMSG,72,NERR,IOPT)
         ENDIF
         MCOLOR(JCOLOR) = MCOLOR(JCOLOR) + 1
         IWKSP(IBC) = IWKSP(IBC) + 1
C
         IELEM = IBNDFAC(1,IBDFACE)
         IVERT = IBNDFAC(2,IBDFACE)
C
         IF    ( IELEM .LT. 1 .OR. IELEM .GT. NELEM   )THEN
            WRITE(ERRMSG(8:72),500)1,IBDFACE,IELEM
            NERR = 6 
            CALL SETERR(ERRMSG,72,NERR,IOPT)
         ELSEIF( IVERT .LT. 1 .OR. IVERT .GT. NOFVERT )THEN
            WRITE(ERRMSG(8:72),500)2,IBDFACE,IVERT
            NERR = 6 
            CALL SETERR(ERRMSG,72,NERR,IOPT)
         ENDIF
C
         IF (LFLAG) THEN
            DO I = 1, (NOFVERT-1) 
               IWORK(I,IBDFACE) = ICELNOD(ICYCL(IVERT+I,NOFVERT),IELEM)
            ENDDO
            IWORK(NOFVERT,IBDFACE) = JCOLOR
         ENDIF 
C
         IFACE = ICELFAC(IVERT,IELEM)
         CALL DAXPY( NDIM , ONE , FACENORM(1,IFACE) ,  1, WKSP ,  1)
C
         TEMP1 = DNRM2(NDIM,FACENORM(1,IFACE),1)
         BNDSURF(IBC) = BNDSURF(IBC) + TEMP1
         SCOLOR(JCOLOR) = SCOLOR(JCOLOR) + TEMP1
         TEMP = TEMP + TEMP1
C
C     .. If the current face is a solid wall with no slip condition
C        applied to it, compute cell aspect ratio ad minimum distance
C        from the wall
C
         IF(IBC.EQ.BC_TYPE_NO_SLIP)THEN
            CALL ARATIO(IVERT,ICELNOD(1,IELEM),COORD,
     +                  FACENORM(1,IFACE),TEMP1,NDIM,NOFVERT,
     +                  AR,H)
C
            NWFAC = NWFAC+1
            ARAVG = ARAVG+AR
            ARMIN = MIN(AR,ARMIN)
            ARMAX = MAX(AR,ARMAX)
            HAVG  = HAVG+H
            HMIN  = MIN(H,HMIN)
            HMAX  = MAX(H,HMAX)
C
         ENDIF
C
   16 CONTINUE ! End loop over boundary faces
C
      IF( NWFAC .NE. 0 )THEN
          HAVG  = HAVG/NWFAC
          ARAVG = ARAVG/NWFAC
      ENDIF
C
      WRITE(NOUT,115)NBFAC
C
C     Print out the number of faces of a given colour
C
      NBODY4=0
      NBODY6=0
C
      DO 12 J = 0 , NCOLOR
      IF(MCOLOR(J).NE.0)THEN
         WRITE(NOUT,355)MCOLOR(J),J
         WRITE(EXT1,FMT="(I2.2)")J
         WRITE(EXT2,FMT="(I3.3)")MY_PE
         FILENAME = 'bndflux.'//EXT1//'.'//EXT2
         OPEN(UNIT=IFUNIT(J),FILE=FILENAME,STATUS='UNKNOWN')
         IF    (ICOLOR(J,1).EQ.BC_TYPE_NO_SLIP)THEN
            NBODY6 = NBODY6+1
            FILENAME = 'no-slip.'//EXT1//'.'//EXT2
            OPEN(UNIT=IMUNIT(J),FILE=FILENAME,STATUS='UNKNOWN')
C           WRITE(NOUT,720)'No-slip  ',J,IMUNIT(J) 
         ELSEIF(ICOLOR(J,1).EQ.BC_TYPE_SLIP_FREE)THEN
            NBODY4 = NBODY4+1
            FILENAME = 'slip-free.'//EXT1//'.'//EXT2
            OPEN(UNIT=IMUNIT(J),FILE=FILENAME,STATUS='UNKNOWN')
C           WRITE(NOUT,720)'Slip-free',J,IMUNIT(J) 
         ENDIF
      ENDIF
   12 CONTINUE
C
C     set pointers to the beginning of the viscous wall boundary
C     faces in IBNDFAC; this is used for referencing the array SKINF
C
C        viscous                          viscous       
C         body#1                           body#2              Array
C     +-------------+-----------+----------------------------+ IBNDFAC
C      ^           ^             ^                          ^
C      |           |             |                          |
C IBGN(1)     IEND(1)       IBGN(2)                    IEND(2)
C
      I = 0
      IFRST = 1
      DO 18 J = 0 , NCOLOR
         IF( ICOLOR(J,1) .EQ. 6 .AND. MCOLOR(J) .GT. 0 )THEN
            I = I+1
            IBGN(I) = IFRST
            IEND(I) = IFRST + MCOLOR(J) -1
         ENDIF
         IFRST = IFRST+MCOLOR(J)
   18 CONTINUE
C
      IF(NBODY4.NE.0)WRITE(NOUT,700)NBODY4
      IF(NBODY6.NE.0)WRITE(NOUT,710)NBODY6
C
C     ... Print out the number of faces and total surface of a given
C         boundary type ...
C
      DO 10 J = 0 , NBTYPE
   10 IF(IWKSP(J).NE.0)WRITE(NOUT,350)CBTYPE(J),IWKSP(J),BNDSURF(J)
      WRITE(NOUT,454)TEMP
C
      LREAD(2) = (IWKSP(BC_TYPE_PRESCRIBED_FLUX).NE.0)
C
      IF(NWFAC.NE.0)
     +WRITE(NOUT,455)NWFAC,HMIN,HMAX,HAVG,ARMIN,ARMAX,ARAVG
C
C Check the surface integral of the scaled normals
C (should be zero)
C
      WRITE(NOUT,330)DNRM2( NDIM , WKSP , 1 )
C
      IF(LFLAG)THEN
         CALL I4Mat_Print('General',' ',NOFVERT,NBFAC,IWORK,
     +            NOFVERT,'Boundary vertices & colour',NERR)
         IF(NERR.NE.0)CALL EXIT(NERR)
      ENDIF
C
C Perform a check on c-lines 
C
      IF(NCL.EQ.0)RETURN
      DO IL = 1, NCL
         JBGN = IA(IL) 
         JEND = IA(IL+1)-1
         DO J = JBGN,JEND
            IBDFACE = JA(4,J)
            IELEM = IBNDFAC(1,IBDFACE)
            IVERT = IBNDFAC(2,IBDFACE)
            NERR = 0
            DO K = 1,NOFVERT-1
               N1 = ICELNOD(ICYCL(IVERT+K,NOFVERT),IELEM)
               N2 = JA(K,J)
               IF(N1.NE.N2)NERR = NERR+1
            END DO
            IF(NERR.NE.0)THEN
               WRITE(6,*)'Inconsistency detected on c-lines'
               WRITE(6,*)'c-line # ',IL,' vertex # ',J
               WRITE(6,*)'jbgn = ',JBGN,' jend = ',JEND
               WRITE(6,*)'bdry f # ',IBDFACE,' tetra # ',IELEM
               WRITE(6,*)(ICELNOD(ICYCL(IVERT+K,NOFVERT),IELEM),K=1,3)
               WRITE(6,*)(JA(K,J),K=1,4)
               STOP
            ENDIF
         END DO ! loop over vertices of the c-line
      END DO ! loop over c-lines
      RETURN
C
C     I/O FORMATS
C
  111 FORMAT('.',$)
  114 FORMAT(10X,'CHECKING BOUNDARY FACES ',$)
  115 FORMAT(I6,/)
  330 FORMAT(/,10X,'BOUNDARY SCALED NORMALS SUM UP TO ',E9.4,/30X,
     .'(should be 0.E+00)',/)
  350 FORMAT(10X,A20,' BOUNDARY FACES : ',I6,2X,'(',D9.4,')')
  355 FORMAT(10X,'THERE ARE ',I5,' BOUNDARY FACES COLOURED ',I2)
  454 FORMAT(57X,9("-"),/57X,D9.4)
  455 FORMAT(10X,'THERE ARE ',I6,' VISCOUS WALL FACES',//,
     +15X,'Minimum wall distance ',E10.4,/, 
     +15X,'Maximum wall distance ',E10.4,/, 
     +15X,'Average wall distance ',E10.4,/, 
     +15X,'Minimum aspect ratio  ',E10.4,/, 
     +15X,'Maximum aspect ratio  ',E10.4,/, 
     +15X,'Average aspect ratio  ',E10.4)
  500 FORMAT('Inconsistent data : IBNDFAC(',I1,',',I6,') = ',I6)
  510 FORMAT('Incorrect boundary type ',I2,
     1' has been assigned to colour ',I2)
  700 FORMAT(/10X,'THERE ARE ',I2,' "Slip-free" BODIES') 
  710 FORMAT(10X,'THERE ARE ',I2,' "No-slip"  BODIES'/) 
  720 FORMAT(10X,'Aerodynamic coefficients for ',A9,' body color ',I2,/
     +20X,'will be written to unit ',I2) 
C
      END
