      subroutine gcoor(xyc,xy,nelem)

      implicit none

      double precision xyc(2,*),xy(2,*)
      double precision xg,yg
      integer nelem
      integer icelnod(3,nelem)
      integer ielem,ivert,i 


      do ielem = 1,nelem
         xg=0.d0 
         yg=0.d0 
         do ivert=1,3
            i=icelnod(ivert,ielem)
            xg=xg+xy(1,i)
            yg=yg+xy(2,i)
         enddo
         xg=xg/3.
         yg=yg/3.
      enddo

      return
      end  
        
