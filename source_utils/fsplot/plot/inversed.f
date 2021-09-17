      subroutine inversed(ord,coor,ndim,npoin,ZROE,NDOF,xyIN,
     &                    ZOUT,info)
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
      integer npoin,ndim,NDOF,info
      integer ord(npoin)
      double precision coor(ndim,npoin),ZROE(NDOF,npoin)
      double precision xyIN(ndim),ZOUT(*),d(npoin)
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
      do ipoin = 1, npoin
c
c     x0,y0 are the coordinates of the point to be located
c
         x0 = xyIN(1)
         y0 = xyIN(2)
         xp(1) = coor(1,ipoin)
         yp(1) = coor(2,ipoin)
c
         d(ipoin) = sqrt((x0-xp(1))**2+(y0-yp(1))**2)
c
      enddo
c
      call qsortd(ord,npoin,d)
c
c     write(6,*)(d(ord(j)),j=1,npoin) 
c
         help = 1.d0/(d(ord(1))+d(ord(2)))
c
         a(1) = d(ord(2))*help
         a(2) = d(ord(1))*help
c
         write(12,*)a(1),a(2),d(ord(1)),d(ord(2))
         do ivar = 1, NDOF 
             ZOUT(ivar) =  0.d0
         enddo
         do iv = 1, 2 
             ipoin = ord(iv)
             help = a(iv)
             do ivar = 1, NDOF 
                ZOUT(ivar) = ZOUT(ivar) + help* ZROE(ivar,ipoin)
             enddo
          enddo
          info = 0
!     write(6,*)'Search failed for Vertex coords ',x0,y0
!     write(6,FMT=1100)(a(iv),iv=1,3),a(1)+a(2)+a(3),s,t
      return
 1100 format('a(i),S,s,t ',6(E12.4,1X))
  300 format(I1,7(F10.5,1X))
  400 format(I1,9(F10.5,1X))
      end
