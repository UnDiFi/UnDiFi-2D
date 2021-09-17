      SUBROUTINE CH3D(SKINF,QFLUX,TINDX,IFACNOD,IBNDPTR,KLIST,RWORK,
     +              NDOF,TAREA,VAREA,YPLUS,YY,LWORK,ICELNOD,COOR,
     +              NDIM,NOFVERT,NOFVAR,REYNO,M_INFTY,FILENAME)
C
C      $Id: ch3d.f,v 1.1 2006/04/03 09:20:49 abonfi Exp abonfi $
C
C      Writes heat flux data in 3D using a Tecplot fmt.
C
       IMPLICIT NONE
C
C     .. Parameters ..
      INTEGER MBODIES
      PARAMETER (MBODIES=10)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION M_INFTY,REYNO
      INTEGER LWORK,NDIM,NOFVAR,NOFVERT,NDOF
      CHARACTER FILENAME* (*)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,*),SKINF(LWORK),TINDX(*),YY(*),
     +QFLUX(*),YPLUS(*),RWORK(NDOF,LWORK),TAREA(LWORK),VAREA(LWORK)
      INTEGER IFACNOD(3,LWORK),ICELNOD(NOFVERT,*),KLIST(LWORK),
     +        IBNDPTR(2,LWORK)
C
C     TAREA is the surface traingle's area
C     VAREA is the median dual area
C
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DNX,DNY,DX,DY,SURF,UTAU,VISCL,DENS
      INTEGER I,IBODY,IC,NPOIN,IELEM,IFLAG,IFRST,ILAST,IPOIN,IUNIT,
     +        IVERT,IXDRS,IXDRS2,J,N1,N2,NBODY6,NEDGE,NWFAC,IFAIL,
     2        IPOS,LAST,INODE,IFACE,K
      LOGICAL COMPRESSIBLE,TURBO
C     ..
C     .. Local Arrays ..
      INTEGER IBGN(MBODIES+1),NFACE(MBODIES)
      DOUBLE PRECISION VARR(3,3)
C     ..
C     .. External Functions ..
      INTEGER  INITXDR,ICYCL
      EXTERNAL INITXDR,ICYCL
C     ..
C     .. External Subroutines ..
      INTEGER  IXDRDMAT,IXDRIMAT,IXDRINT,IXDRCLOSE
      EXTERNAL IXDRDMAT,IXDRIMAT,IXDRINT,IXDRCLOSE
C     ..
C     .. Intrinsic Functions ..
      DOUBLE PRECISION DNRM2,YDIST
C     INTRINSIC DNRM2
C     ..
      IF (NDIM.EQ.2) THEN
         WRITE(6,*)'SUBROUTINE CH only works with NDIM=3'
         CALL EXIT(1)
      ENDIF
C
      COMPRESSIBLE = (NOFVAR.EQ.5)
      IF( COMPRESSIBLE )THEN
          VISCL = 1.d0
          DENS  = 1.d0
      WRITE (6,FMT=*) 'Warning!Warning!Warning!Warning!Warning!Warning!'
      WRITE (6,FMT=*) 'Warning!Warning!Warning!Warning!Warning!Warning!'
      WRITE (6,FMT=*) 'Surface density is set to ',DENS
      WRITE (6,FMT=*) 'Surface viscosity is set to ',VISCL
      WRITE (6,FMT=*) 'Warning!Warning!Warning!Warning!Warning!Warning!'
      WRITE (6,FMT=*) 'Warning!Warning!Warning!Warning!Warning!Warning!'
      ELSEIF( .NOT. COMPRESSIBLE )THEN
          VISCL = 1.d0
          DENS  = 1.d0
      ENDIF
C
      IUNIT = 30
      WRITE (6,FMT=*) 'Writing heat transfer data to unit = ',IUNIT
      WRITE (6,FMT=*) 'Reynolds number is ',REYNO
      WRITE (6,FMT=*) 'Freestream Mach is ',M_INFTY
C
C
      IXDRS = INITXDR(FILENAME,'r',.FALSE.)
C
C     file017.dat is the file containing the turbulent index
C
      INQUIRE(FILE='file017.dat',EXIST=TURBO)
      IF(TURBO)THEN
         IXDRS2 = INITXDR('file017.dat','r',.FALSE.)
      ELSE
         WRITE(6,*)'No turbulent data'
      ENDIF
C
      IFAIL = IXDRINT(IXDRS,NBODY6)
      IFAIL = IXDRINT(IXDRS,NWFAC)
C
      IBGN(1) = 1
      NWFAC = 0
C
      WRITE (6,FMT=*) 'There are ',NBODY6,' No-slip boundaries'
      DO 10 IBODY = 1,NBODY6
C
C     number of faces for the current body
C
         IFAIL = IXDRINT(IXDRS,NFACE(IBODY))
C
         IBGN(IBODY+1) = IBGN(IBODY) + NFACE(IBODY)
         WRITE (6,FMT=*) 'No-slip body # ',IBODY,NFACE(IBODY),' faces'

         NWFAC = NWFAC + NFACE(IBODY)
   10 CONTINUE

      IFAIL = IXDRIMAT(IXDRS,2*NWFAC,IBNDPTR)
      IFAIL = IXDRDMAT(IXDRS,NWFAC,SKINF)
      IFAIL = IXDRDMAT(IXDRS,NWFAC,QFLUX)
      IF(TURBO)THEN
          IFAIL = IXDRDMAT(IXDRS2,NWFAC,TINDX)
      ELSE
          CALL DINIT(NWFAC,ZERO,TINDX,1)
      ENDIF
C        
      OPEN(UNIT=IUNIT,FILE='skin.dat')
      WRITE(IUNIT,*)'TITLE = "skin data"'
      WRITE(IUNIT,*)'VARIABLES = "X", "Y", "Z", "Ch" "Cf" "Tindex" "yn" 
     &"y+"'
C
C     loop over no-slip bodies
C
      WRITE (6,FMT=*) 'Sorting out boundary faces .......' 
      DO 15 IBODY = 1,NBODY6
      WRITE (6,FMT=*) 'No-slip body # ',IBODY,NFACE(IBODY),' faces'
        IFRST = IBGN(IBODY)
        ILAST = IBGN(IBODY+1) - 1
        IFLAG = 0
        DO 17 J = IFRST,ILAST
          IELEM = IBNDPTR(1,J)
          IVERT = IBNDPTR(2,J)
C
C     here we build the face to vertex pointer  
C     and copy the edge lenghts into a temp array
C     nodal indices are in global numbering
C
          DO 18 I = 1, (NOFVERT-1) 
            IC = ICELNOD(ICYCL(IVERT+I,NOFVERT),IELEM) 
            IFACNOD(I,J) = IC
            CALL DCOPY(NDIM,COOR(1,IC),1,VARR(1,I),1)
   18     CONTINUE
C
C         compute (P(3)-P(1)) and (P(3)-P(2))
C
            CALL DAXPY(NDIM,-1.d0,VARR(1,1),1,VARR(1,2),1)
            CALL DAXPY(NDIM,-1.d0,VARR(1,1),1,VARR(1,3),1)
C
C         compute the area of the triangle... 
C
            CALL CROSS_PROD( VARR(1,2), VARR(1,3), VARR(1,1) )
            TAREA(J) = DNRM2(NDIM,VARR(1,1),1)
C
            YY(J) = YDIST(ICELNOD,NOFVERT,COOR,NDIM,IELEM,IVERT)
c non-dimensional viscosity and density
            UTAU = SQRT(SKINF(J)/DENS)
            YPLUS(J) = REYNO*UTAU/VISCL*YY(J)
C
   17   CONTINUE
C
        N1 = NFACE(IBODY)*(NOFVERT-1)
        CALL ICOPY(N1,IFACNOD(1,IFRST),1,KLIST,1)
C
C       sort the list of nodes and get rid of duplicates
C
        WRITE(6,*)'Removing duplicated entries .....'
        CALL SORTSP(N1,KLIST,N2)
        NPOIN = N2
        WRITE(6,*)'Done ! '
C
C       now KLIST keeps the global nodenumbers of the meshpoints
C       belonging to the bondary faces
C
C       initialize median dual area and nodal values to 0.d0
C
        CALL DINIT(NPOIN,0.d0,VAREA,1)
        CALL DINIT(NDOF*NPOIN,0.d0,RWORK,1)
C
        WRITE (6,FMT=*)'No-slip body # ',IBODY,' found ',N2,' bnd nodes'
C
        DO 19 J = IFRST,ILAST
          DO 21 I = 1, (NOFVERT-1) 
            SURF = TAREA(J)
            N1 = IFACNOD(I,J)
            CALL BINSRC(N1,KLIST,NPOIN,IPOS,LAST)
C
C       we now update the connectivity so that it addresses
C       local nodenumbers 
C
            IF( IPOS .EQ. 0 )THEN

                WRITE(6,*)'IFRST, J, ILAST = ',IFRST,J,ILAST
                WRITE(6,*)'Cannot find ',N1,' in KLIST'
                WRITE(6,*)'ifacnod = ',(IFACNOD(K,J),K=1,3)
                WRITE(6,*)'klist = ',(KLIST(K),K=1,NPOIN)
                CALL X04EAF('General',' ',NOFVERT-1,NFACE(IBODY),
     +          IFACNOD,NOFVERT-1,'Face connectivity',IFAIL)
                CALL EXIT(1)
            ELSE
                IFACNOD(I,J) = IPOS
            ENDIF 
   21     CONTINUE
C
C    compute median dual
C
          DO 23 I = 1, (NOFVERT-1) 
             IPOIN = IFACNOD(I,J)
             VAREA(IPOIN) = VAREA(IPOIN) + SURF
             RWORK(1,IPOIN) = RWORK(1,IPOIN) + SURF * QFLUX(J)
             RWORK(2,IPOIN) = RWORK(2,IPOIN) + SURF * SKINF(J) * 2.d0
             RWORK(3,IPOIN) = RWORK(3,IPOIN) + SURF * TINDX(J)
             RWORK(4,IPOIN) = RWORK(4,IPOIN) + SURF * YY(J)
             RWORK(5,IPOIN) = RWORK(5,IPOIN) + SURF * YPLUS(J)
   23     CONTINUE
   19   CONTINUE ! end loop over faces
C
        DO 33 IPOIN = 1, NPOIN
           DO 33 IPOS = 1,NDOF
           RWORK(IPOS,IPOIN) = RWORK(IPOS,IPOIN) / VAREA(IPOIN)
   33   CONTINUE
C
      WRITE(IUNIT,*)'ZONE T=" Boundary # ',IBODY,
     &'", F=FEPOINT,ET=TRIANGLE,N=',NPOIN,',E=',NFACE(IBODY)
C
C     .. Write nodes ..
C
         DO 31 IPOIN = 1,NPOIN
               INODE = KLIST(IPOIN)
               WRITE(IUNIT,*)(COOR(J,INODE),J=1,NDIM),
     &         (RWORK(J,IPOIN),J=1,NDOF)
   31    CONTINUE
         WRITE(IUNIT,*)
C
C	writing the elements
C
         DO 25 IFACE = IFRST , ILAST
            WRITE(IUNIT,*)(IFACNOD(J,IFACE),J=1,(NOFVERT-1))
   25    CONTINUE
C
    5 CONTINUE
C
   15 CONTINUE ! End of the outermost loop on colors
C
      CLOSE (IUNIT)
      CLOSE(3)

      LWORK = NWFAC
      IFAIL = IXDRCLOSE( IXDRS )
      IF(TURBO)IFAIL = IXDRCLOSE( IXDRS2 )
      WRITE(*,*)'Done writing tecplot file skin.dat'

  100 FORMAT (F14.9,2X,5 (E12.6,1X))

      RETURN

      END
