      SUBROUTINE CUT3D(ICLR,IRAKE,RBGN,REND,ZROE,ZRAKE,XYZRAKE,LRAKE,
     + NOFVAR,COOR,NDIM,NPOIN,ICELNOD,ICELCEL,NOFVERT,NELEM,IBNDFAC,
     + NBFAC,IPNTR,SKINF,NWFAC,REYNO,NORMAL_TO_THE_WALL)
C
C     SId:$
C
C
      IMPLICIT NONE
C
C     This routine finds the values of the
C     dependent variables along the segment of endpoints:
C     A and B
C     The line is traced starting from its intersection
C     with the boundary face coloured ICLR
C
C
C     .. Scalar Arguments ..
      REAL*8 REYNO
      INTEGER IRAKE,ICLR,NBFAC,NDIM,NELEM,NOFVAR,NOFVERT,NPOIN,
     +NWFAC,LRAKE
      CHARACTER ZHEAD*80
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,NPOIN),SKINF(NWFAC),ZROE(NOFVAR,NPOIN),
     +ZRAKE(NOFVAR+1,LRAKE),RBGN(*),REND(*),XYZRAKE(NDIM,LRAKE)
      INTEGER IBNDFAC(3,NBFAC),ICELCEL(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM),IPNTR(2,NWFAC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DIST,HELP
      INTEGER I,JVERT,IELEM,IFACE,NINTSC,IUNIT,IVERT,IVAR,IXDRS,J,JV,
     +        IDUMMY,IP,JELEM,KUNIT,IFAIL,INFO,N,IDIM
      LOGICAL COMPRESSIBLE,NORMAL_TO_THE_WALL,INSIDE
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION X(3),Y(3),Z(3),POUT(3),AC(3),XLAST(3),tri(3,3),
     &vvt(3,3)
      INTEGER IDX(4)
C     ..
C     .. External Functions ..
      INTEGER  ICYCL
      EXTERNAL ICYCL
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
      DOUBLE PRECISION DNRM2
C     ..
C     .. Data statements ..

C
      DATA AC/3*0.d0/
C
C
      write(6,*)'Log file is cut3d.log'
      open(16,FILE='cut3d.log')
      write(16,*)'Rake starts at ',(rbgn(i),i=1,3)
      write(16,*)'Rake ends   at ',(rend(i),i=1,3)
      write(16,*)
      write(16,*)'Searching the boundary face that intersects the rake'
      write(16,*)
C
C     first of all we find which of the bndry faces
C     is intersected by the rake
C
C     no. of intersections
C
      NINTSC = 0
C
      DO 1 I = 1,NBFAC
C ... find the intersection of the rake with the boundary
          IF (IBNDFAC(3,I).NE.ICLR) GOTO 1
          IELEM = IBNDFAC(1,I)
          IVERT = IBNDFAC(2,I)
C ... nodes on the boundary
          DO 2 N = 1,3
             IP = ICELNOD(ICYCL(IVERT+N,NOFVERT),IELEM)
             IDX(N) = IP
             X(N) = COOR(1,IP)
             Y(N) = COOR(2,IP)
             Z(N) = COOR(3,IP)
             do idim = 1, ndim 
                tri(idim,n) = coor(idim,ip)
             enddo
    2     CONTINUE
          write(16,*)'Scanning bndry face ',i,' of cell ',ielem,
     &    'facing vertex ', ivert
C
C     check for the intersection
C
          call triangle_contains_line_exp_3d( tri, rbgn, rend, 
     &                                        inside, pout )
!         CALL XAA23S(x,y,z,rbgn,rend,Pout,ac,info)
          if( inside )then
C
             call triangle_area_3d ( tri, help )
             help = 1.d0/help 
             do j = 1, (nofvert-1)
                call dcopy(ndim*ndim,tri,1,vvt,1)
!               call dcopy(ndim,tri(1,j),1,vvt(1,j),1)
                call dcopy(ndim,pout,1,vvt(1,j),1)
                call triangle_area_3d ( vvt, ac(j) )
                ac(j) = ac(j)*help
             enddo
               write(6,*)'Intersection found face is ',i
               write(6,*)(pout(j),j=1,3)
               write(16,*)
               write(16,*)'Intersection found face is ',i
               write(16,*)(pout(j),j=1,3)
               write(16,*)'weigths are ',(ac(idim),idim=1,ndim),' sum up
     & to ',ac(1)+ac(2)+ac(3)
               NINTSC = NINTSC + 1
               call intp3d(idx,ac,zroe,zrake(1,NINTSC),nofvar,nofvert-1)
               goto 4
          else
!              CALL X04CAF('General',' ',NDIM,NOFVERT-1,tri,NDIM,
!    +            'Triangle vertices ',IFAIL)
          endif

    1 CONTINUE ! end loop over all boundary faces
      WRITE(16,*) 'smthg went wrong while searching the bndry face colou
     &red ',iclr
      WRITE(6,*) 'smthg went wrong while searching the bndry face colour
     &ed ',iclr
      CLOSE(16)
      CALL EXIT(1)
    4 CONTINUE
C
C     POUT is where the rake starts
C
      CALL DCOPY(NDIM,POUT,1,XYZRAKE(1,1),1)
C
C     the distance in the 1st point is 0.d0
C
      ZRAKE(NOFVAR+1,NINTSC) = 0.d0
C
      WRITE(ZHEAD(1:51),FMT=450)IRAKE,(pout(j),j=1,3)
  450 FORMAT('Rake # ',I2.2,' starting at',3(2X,F8.4)) 
      WRITE(6,*)ZHEAD(1:51)
      WRITE(6,*)'starting bndry face is ',i
      WRITE(6,*)'starting bndry elmt is ',ielem
      WRITE(16,*)ZHEAD(1:51)
      WRITE(16,*)'starting bndry face is ',i
      WRITE(16,*)'starting bndry elmt is ',ielem
C
C     IELEM is the element where the rake starts
C     IVERT is the vertex opposite that face
C
   12 CONTINUE
      write(16,*)
      WRITE(16,*)'current tet is ',ielem
C
C     scan the other faces of the current elmt (IELEM)
C
      DO 5 I = 1,(NOFVERT-1)
          JVERT = ICYCL(IVERT+I,NOFVERT)
c
c     store the coords of the face opposite vertex JVERT 
c
          DO 6 N = 1,3
             IP = ICELNOD(ICYCL(JVERT+N,NOFVERT),IELEM)
             IDX(N) = IP
             X(N) = COOR(1,IP)
             Y(N) = COOR(2,IP)
             Z(N) = COOR(3,IP)
             do idim = 1, ndim 
                tri(idim,n) = coor(idim,ip)
             enddo
    6     CONTINUE
C
C     check for the intersection
C
!         CALL XAA23S(x,y,z,rbgn,rend,XLAST,ac,info)
          call triangle_contains_line_exp_3d( tri, rbgn, rend, 
     &                                        inside, XLAST )
          write(16,*)'checking face opposite vertex ',icelnod(jvert,iele
     &m),' of tet ',ielem,' neighb tet is ',icelcel(jvert,ielem)
          write(16,*)'face is bounded by vertices: ',
     &(ICELNOD(ICYCL(JVERT+N,NOFVERT),IELEM),n=1,nofvert-1)
          do n = 1,nofvert-1
          write(16,*)(tri(idim,n),idim=1,ndim)
          enddo
          if( inside )then
!         if( info .EQ. 0 )then
               write(16,*)'Intersection found; face is  ',i,' x,y,z ',
     >         (xlast(j),j=1,3),zrake(NOFVAR+1,NINTSC)
               call triangle_area_3d ( tri, help )
               help = 1.d0/help 
               do j = 1, (nofvert-1)
                  call dcopy(ndim*(nofvert-1),tri,1,vvt,1)
                  call dcopy(ndim,xlast,1,vvt(1,j),1)
                  call triangle_area_3d ( vvt, ac(j) )
                  ac(j) = ac(j)*help
               enddo
               write(16,*)'weigths are ',(ac(idim),idim=1,ndim),' sum up
     & to ',ac(1)+ac(2)+ac(3)
               NINTSC = NINTSC + 1
               call intp3d(idx,ac,zroe,zrake(1,NINTSC),nofvar,nofvert-1)
C     compute the distance and store
C     since we move along a straight line it is sufficient
C     to compute the norm of (XLAST-POUT)
C
               zrake(NOFVAR+1,NINTSC) = zrake(NOFVAR+1,NINTSC) + sqrt(
     >         (POUT(1)-XLAST(1))**2+
     >         (POUT(2)-XLAST(2))**2+
     >         (POUT(3)-XLAST(3))**2)
               CALL DCOPY(NDIM,xlast,1,XYZRAKE(1,NINTSC),1)
C              write(6,*)nintsc,(XYZRAKE(j,NINTSC),j=1,ndim)
C              call DCOPY(3,POUT,1,XLAST,1)
C
C         JELEM is the element that shares the face
C               intersected by the rake
C
                  JELEM = ICELCEL(JVERT,IELEM)
                  write(16,*)'nearby elmt is ',jelem
C
C         if a bndry element skip
C
                  IF(JELEM.LE.0.OR.JELEM.GT.NELEM)GOTO 9
C
C         we look for the vertex facing the intersected face
C         in element JELEM
C
                  DO 7 JV = 1,NOFVERT
                     IF( ICELCEL(JV,JELEM) .EQ. IELEM )THEN
                         IELEM = JELEM
                         IVERT = JV
                         GOTO 12 
                     ENDIF
    7             CONTINUE
                  write(16,*)(ICELCEL(J,JELEM),J=1,4)
                  write(16,*)(ICELCEL(J,IELEM),J=1,4)
                  WRITE(6,*) 'smthg else is wrong'
                  CALL EXIT(1)
          else
              write(16,*)'triangle_contains_line_exp_3d returned ',
     &inside
          endif ! info

          
    5 CONTINUE ! loop over the faces of the current tetrahedron
          WRITE(16,*) ' should not get here'
          WRITE(16,*) ' current tet is ',ielem
          WRITE(6,*) ' should not get here'
          WRITE(6,*) ' current tet is ',ielem
          CALL EXIT(1)
    9 CONTINUE
      write(16,*)'Found a boundary face '
      write(16,*)'Last neighbour is ',jelem
      write(6,*)'Found ',nintsc,' intersections for rake # ',irake
C
      call WriteXYZone(ZHEAD(1:51),zrake,nofvar,nintsc)
      call WriteMacro(irake,xyzrake,ndim,nintsc)
C
      RETURN

   80 FORMAT (8 (E16.9,1X))
  100 FORMAT (/10("*"),'Probing normal to the wall',/)
      END
