      subroutine areacoords( t, p, xsi )
      implicit none
      double precision t(2,3),p(3),xsi(3) 
      double precision absarea 
      double precision xA,xB,xC,help
      double precision yA,yB,yC
      integer i,j
      absarea(xA,xB,xC,yA,yB,yC) = 0.5d0* abs(
     &(xA-xC)*(yB-yA)-(xA-xB)*(yC-yA) )
c     
      help = absarea(t(1,1),t(1,2),t(1,3),
     &t(2,1),t(2,2),t(2,3))
      if(help.LT.1.d-12)then
         write(6,*)'Triangle area may be too small'
         write(6,*)((t(j,i),j=1,2),i=1,3)
         call exit(1)
      endif
      help = 1.d0/help
      xsi(1) = absarea(p(1),t(1,2),t(1,3),
     &p(2),t(2,2),t(2,3)) * help
      xsi(2) = absarea(p(1),t(1,3),t(1,1),
     &p(2),t(2,3),t(2,1)) * help
      xsi(3) = 1.d0-xsi(1)-xsi(2)
c
      return
      end
