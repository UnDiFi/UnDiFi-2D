      subroutine subyz(icn,nofvert,coor,ndim,pcn,lda,rcline,zb,ncl)
      implicit none
c
c     $Id: subyz.f,v 1.1 2005/12/23 09:50:39 abonfi Exp $     
c
c     interpolates pressure on the outflow plane
c     rcline(1:ncl) radial coord of the c-lines
c     zb(1:ncl,nofvar) Roe variables, circumferentially averaged
c                      along the c-lines
c
      include 'paramt.h'
      include 'bnd.h'
      include 'stream.com'
      integer ndim,nofvert,lda,ncl
      integer icn(nofvert)
      double precision coor(ndim,*),pcn(lda,nofvert),
     &rcline(ncl),zb(ncl,*)
      integer ipoin,ivert
      double precision rr,yp,zp
      double precision value
c
      do 1 ivert= 1, nofvert-1
c
c     here we assume that the axis of the machine is the x axis
c
         ipoin = icn(ivert)
         yp = coor(2,ipoin)
         zp = coor(3,ipoin)
         rr = sqrt(yp*yp+zp*zp)
         pcn(3,ivert) = value(rr,rcline,zb,ncl-1)
    1 continue
      return
  300 FORMAT('SUBYZ CANNOT LOCATE NODE ',I6.6,
     +       ' AMONG THOSE IN THE LIST')
      end

      double precision function value(x0,x,y,n)

      implicit none
      integer n
      integer i,im1
      double precision x(0:n),y(0:n),x0

      if(x0.LT.x(0))THEN
!        write(6,*)'x0 ',x0,' out of interval  [',x(0),',',x(n),']'
         value = y(0)+(x0-x(0))/(x(1)-x(0))*(y(1)-y(0)) 
         return
      endif 
      if(x0.GT.x(n))THEN
!        write(6,*)'x0 ',x0,' out of interval  [',x(0),',',x(n),']'
         value = y(n)+(x0-x(n))/(x(n)-x(n-1))*(y(n)-y(n-1)) 
         return
      endif 

      do 3 i = 1, n

         im1 = i-1
         if( (x0-x(im1))*(x0-x(i)) .LE. 0 )THEN
              value = y(im1)+(x0-x(im1))/(x(i)-x(im1))*(y(i)-y(im1)) 
caldo         write(6,*)i,x(im1),x0,x(i),value
         return
         endif
    3 continue
      write(6,*)x0
      write(6,*)'xmin = ',x(0)
      write(6,*)'xmax = ',x(n)
      stop 'problem in function value'
      end
