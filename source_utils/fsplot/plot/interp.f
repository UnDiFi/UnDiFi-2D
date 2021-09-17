      subroutine interp(icelnod,coor,ndim,ZROE,NDOF,xyIN,
     &ZOUT,ielem,info,a)
C
C     input:
C            xyIN nodal coords (belonging to the background grid)
C                 to be located inside the current grid 
C            icelnod cell to node pointer
C            coor nodal coordinates
C            ndim space dimension
C            ndof nof degrees of freedom
C            ielem is the cell node xyIN falls inside
C     output:
C            info = 0 node found !=0 search failed
C            ZOUT(*) is filled with the interpolated value
C
      implicit none
      integer ielem,ndim,NDOF,info
      integer icelnod(3,*)
      double precision coor(ndim,*),ZROE(NDOF,*)
      double precision xyIN(ndim),ZOUT(*)
      double precision x0,y0,xp(4),yp(4),a(3),t,s,help
      integer idxs(3),iv,ipoin,ivar,ilog,j
      double precision EPS
      parameter(EPS=1.e-9,ilog=1)
      logical LTEST
      double precision area
      integer icycl
c
      info = 0
c
c     x0,y0 are the coordinates of the point to be located
c
      x0 = xyIN(1)
      y0 = xyIN(2)
c
      do 50 iv = 1,3
         ipoin = icelnod(iv,ielem) ! global nodenumber
         xp(iv) = coor(1,ipoin)
         yp(iv) = coor(2,ipoin)
   50 continue
         xp(4) = x0 ! node to be located
         yp(4) = y0 ! node to be located
c
         idxs(1) = 1
         idxs(2) = 2
         idxs(3) = 3
         help = 1.d0/area(xp,yp,3,idxs)
         idxs(3) = 4 ! node to be located
c
c     compute area coordinates (in the x-y plane)
c
         do iv = 1,3
            idxs(1) = icycl(1+iv,3)
            idxs(2) = icycl(2+iv,3)
            a(iv) = area(xp,yp,3,idxs)*help
!           if(ABS(a(iv)).LE.EPS)a(iv)=0.d0
         enddo
c
         s = min( a(1), a(2), a(3) )
         t = max( a(1), a(2), a(3) )
         LTEST = ( s .GE. 0.d0 .AND. s .LE. 1.d0 ) .AND.
     +           ( t .GE. 0.d0 .AND. t .LE. 1.d0 )
!        LTEST = ( ABS(s) .LE. EPS .AND. ABS(s-1.d0) .LE. EPS ) .AND.
!    +           ( ABS(t) .LE. EPS .AND. ABS(t-1.d0) .LE. EPS )
!        IF( .NOT. LTEST )THEN
!            info = -1
!            RETURN
!        ENDIF
!     write(6,*)'Area coords ',(a(j),j=1,3),' cell ',ielem,help
!     write(6,*)'y/z coords ',(xp(j),j=1,4),(yp(j),j=1,4)
!        if( ( s .GE. 0.d0 .AND. s .LE. 1.d0 ) .AND.
!    +       ( t .GE. 0.d0 .AND. t .LE. 1.d0 ) )then
!        if(ilog.EQ.0)write(6,FMT=*)(icelnod(iv,ielem),iv=1,3)
          do ivar = 1, NDOF 
             ZOUT(ivar) =  0.d0
          enddo
          do iv = 1, 3 
             ipoin = icelnod(iv,ielem)
             help = a(iv)
             do ivar = 1, NDOF 
                ZOUT(ivar) = ZOUT(ivar) + help* ZROE(ivar,ipoin)
             enddo
          enddo
          if(ilog.EQ.0)then
             write(6,*)ielem
             do ivar = 1, NDIM 
             write(6,FMT=400)ivar,xyIN(ivar),
     3         (coor(ivar,icelnod(iv,ielem)),iv=1,3),
     3         (a(iv),iv=1,3)
             enddo
             do ivar = 1, NDOF 
             write(6,FMT=300)ivar,ZOUT(ivar),
     &         (ZROE(ivar,icelnod(iv,ielem)),iv=1,3),
     3         (a(iv),iv=1,3)
             enddo
          endif ! ilog
          info = 0
          return
!     write(6,*)'Search failed for Vertex coords ',x0,y0
!     write(6,FMT=1100)(a(iv),iv=1,3),a(1)+a(2)+a(3),s,t
      return
 1100 format('a(i),S,s,t ',6(E12.4,1X))
  300 format(I1,7(F10.5,1X))
  400 format(I1,9(F10.5,1X))
      end
