! Change the shock topology for the regular reflection of a shock on a wedge

      subroutine ch_sh_topology(
     +           xysh,
     +           zroeshu,
     +           zroeshd,
     +           wsh,
     +           xywedge,
     +           nshocks,
     +           typesh,
     +           nshockpoints,
     +           nspecpoints,
     +           typespecpoints,
     +           isppnts,
     +           shinspps,
     +           ispclr,
     +           nclr,
     +           iclr,
     +           ia,
     +           ja,
     +           zroe,
     +           corg)

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer nshocks,nspecpoints,nshockpoints(*),
     &        shinspps(2,5,*),
     &        ispclr(5,*)
      character*5 typespecpoints(*)
      character*1 typesh(*)

!     .. array arguments ..
      double precision
     +                 xysh(ndim, npshmax, *),
     +                 zroeshu(ndof, npshmax, *),
     +                 zroeshd(ndof, npshmax, *),
     +                 wsh(ndim, npshmax, *)
      integer nclr
      integer ia(*),ja(*),iclr(nclr)
      double precision corg(ndim,*),zroe(ndof,*)

!     .. local scalars ..
      integer isppnts,ip,ip1,ish1,i,j
      integer clr,bbgn,bend, bgnwedge
      double precision x1wedge,y1wedge,
     +                 x2wedge,y2wedge,
     +                 x1wall,y1wall,
     +                 x2wall,y2wall,
     +                 x1sh,y1sh,x2sh,y2sh
      double precision xywedge(2)
      double precision dum,dum1,dum2,dumx1,dumy1
      double precision dx,dy,rad,wedgeangle
      double precision ui,vi
      double precision tau,taux,tauy,nx,ny,xi,yi
      integer nn,ish
      parameter (nn=2)
      double precision a(nn,nn), b(nn),x(nn)

      open(8,file="log/ch_sh_topology.log")

      write(8,*)'the shock topology is changed'

      ish1 = shinspps(1,1,isppnts)
      i    = shinspps(2,1,isppnts)-1
      ip1  = 1+i*(nshockpoints(ish1)-1)

!     update the number of shocks in the domain
      nshocks=nshocks+1
      typesh(nshocks)='s'
      write(8,'(a,i3)')'the number of shocks is now',nshocks

!     update the number of special points in the domain
      nspecpoints=nspecpoints+1
      write(8,'(a,i4)')'the number of special points is now',
     &          nspecpoints

!     change the type of isppnts-special point from FWP to rr
!     and fix the shinspps elements for the second shock
      typespecpoints(isppnts)='rr'

caldo
c     shinspps(1,2,isppnts)=shinspps(1,1,isppnts)
c     shinspps(2,2,isppnts)=shinspps(2,1,isppnts)
caldo

      shinspps(1,2,isppnts)=nshocks
      shinspps(2,2,isppnts)=1
      ispclr(2,isppnts)=ispclr(1,isppnts)

      write(8,'(a,i2)')'added rr special point on
     &           the coloured boundary',ispclr(2,isppnts)

!     add wpnrx special point as last special point
!     and fix the shinspps elements

      typespecpoints(nspecpoints)='FWP'
!     typespecpoints(nspecpoints)='WPNRX'
      shinspps(1,1,nspecpoints)=nshocks
      shinspps(2,1,nspecpoints)=2
      ispclr(1,nspecpoints)=ispclr(1,isppnts)

      write(8,'(a,i2)')'added FWP special point on
     +           the coloured boundary',ispclr(1,nspecpoints)

!     change the topology of the first shock: the point
!     crossing the wedge is aligned with the other shock points

!     set the colour of the boundary on which the shock is located
      clr=5
      do clr=1,nclr
      if(iclr(clr).eq.ispclr(1,isppnts))exit
      enddo

!     find the first and the last points of the wedge
      bbgn=ia(clr)
      bend=ia(clr+1)-1
      do j = bbgn, bend-1
       if((xywedge(1)-corg(1,ja(j))).lt.1e-13.and.
     &    (xywedge(2)-corg(2,ja(j))).lt.1e-13) then
        bgnwedge=ja(j)
        exit
       endif
      enddo

      x1wedge=xywedge(1)
      y1wedge=xywedge(2)
      x2wedge=corg(1,ja(bend))
      y2wedge=corg(2,ja(bend))

      x1sh=xysh(1,nshockpoints(ish1),ish1)
      y1sh=xysh(2,nshockpoints(ish1),ish1)
      x2sh=xysh(1,ip1+1,ish1)
      y2sh=xysh(2,ip1+1,ish1)

      a(1,1)=(y1sh-y2sh)
      a(1,2)=(x2sh-x1sh)
      b(1)=y1sh*x1sh-y2sh*x1sh

      a(2,1)=y1wedge-y2wedge
      a(2,2)=x2wedge-x1wedge
      b(2)=y1wedge*x2wedge-y2wedge*x1wedge

      call solg(nn,nn,a,b,x)
      xi = x(1)
      yi = x(2)

!     set new position of the shock point
      xysh(1,ip1,ish1)=xi
      xysh(2,ip1,ish1)=yi

      write(8,*)'shock no.',ish1
      write(8,'(a,1x,i2,1x,a)')'position of point',ip1,'is changed'
      write(8,'(a,i2,a,1x,f16.8)')'x(',ip1,'):',xysh(1,ip1,ish1)
      write(8,'(a,i2,a,1x,f16.8)')'y(',ip1,'):',xysh(2,ip1,ish1)
      write(8,'(a,1x,i2,1x,a)')'position of point',ip1+1,':'
      write(8,'(a,i2,a,1x,f16.8)')'x(',ip1+1,'):',xysh(1,ip1+1,ish1)
      write(8,'(a,i2,a,1x,f16.8)')'y(',ip1+1,'):',xysh(2,ip1+1,ish1)

!     set the topology of the second shock:
!     the first point is the RR point
!     the second point is located at dxcell from the wedge,
!     perpendicular to the wedge
!     the last point is located on a
!     circumference with radius set equal to the distance
!     between the second point and the beginning of the wedge

      nshockpoints(nshocks)=3

      write(8,'(a,i3)')'shock no.',nshocks

!     first shock point
      xysh(1,1,nshocks)=xysh(1,ip1,ish1)
      xysh(2,1,nshocks)=xysh(2,ip1,ish1)

      write(8,'(a,1x,i1,1x,a)')'position of point',1,'is'
      write(8,'(2(1x,a,1x,f16.8))')'x(1):',xysh(1,1,nshocks),
     &                             'y(1):',xysh(2,1,nshocks)

!     second shock point

!     compute the tangent and the normal vectors to the wedge
      taux=x2wedge-x1wedge
      tauy=y2wedge-y1wedge
      tau=sqrt(taux**2+tauy**2)
      taux=taux/tau
      tauy=tauy/tau
      nx=-taux
      ny=tauy

!     compute the angle of the wedge
      dx=abs(xysh(1,1,nshocks)-x1wedge)
      dy=abs(xysh(2,1,nshocks)-y1wedge)
      wedgeangle=atan(dy/dx)

!     compute the coordinates of projection of the second point on the wedge
      dumx1=xysh(1,1,nshocks)-taux*dxcell*sin(wedgeangle)**2
      dumy1=xysh(2,1,nshocks)-
     +      tauy*dxcell*sin(wedgeangle)*cos(wedgeangle)

!     compute the coordinates of the second point
      xysh(1,2,nshocks)=dumx1+nx*dxcell*cos(wedgeangle)
      xysh(2,2,nshocks)=dumy1+ny*dxcell*cos(wedgeangle)

      write(8,'(a,1x,i1,1x,a)')'position of point',2,'is'
      write(8,'(2(1x,a,1x,f16.8))')'x(2):',xysh(1,2,nshocks),
     &                             'y(2):',xysh(2,2,nshocks)

      dum=sqrt((xysh(1,2,nshocks)-xysh(1,1,nshocks))**2+
     +         (xysh(2,2,nshocks)-xysh(2,1,nshocks))**2)

!     third shock point

!     find the first and the last points of the wall
      x1wall=corg(1,ja(bbgn))
      y1wall=corg(2,ja(bbgn))
      x2wall=xywedge(1)
      y2wall=xywedge(2)

!     compute the tangential and normal vectors to the wall
      taux = x2wall-x1wall
      tauy = 0d0
      tau=sqrt(taux**2+tauy**2)
      taux=taux/tau
      tauy=tauy/tau
      nx = 0d0
      ny = taux

!     compute the distance between the beginning of the wedge and the second point
      rad = sqrt((xysh(1,2,nshocks)-x1wedge)**2+
     &           (xysh(2,2,nshocks)-y1wedge)**2)

!     find the coordinates of the fifth shock point
      xysh(1,3,nshocks)=x1wedge-taux*rad
      xysh(2,3,nshocks)=y1wedge

!     find the coordinates of the fourth point
      write(8,'(a,1x,i1,1x,a)')'position of point',3,'is'
      write(8,'(2(1x,a,1x,f16.8))')'x(3):',xysh(1,3,nshocks),
     &                             'y(3):',xysh(2,3,nshocks)

!     set the upstream and downstream state for the new shock points
      write(8,*)'upstream and downstream state of the new'
      write(8,*)'shock points'

      do ip=1,nshockpoints(nshocks)
       do i=1,ndof
        zroeshu(i,ip,nshocks)=zroeshd(i,ip1+1,ish1)
        write(8,'(a,i1,a,i1,a,i2,a,f16.8)')
     +   'zroeshu(',i,',',ip,',',nshocks,
     +   '):',zroeshu(i,ip,nshocks)
       enddo
      enddo

      do ip=1,nshockpoints(nshocks)
       do i=1,ndof
        zroeshd(i,ip,nshocks)=zroe(i,bgnwedge)
        write(8,'(a,i1,a,i1,a,i2,a,f16.8)')
     +   'zroeshd(',i,',',ip,',',nshocks,
     +   '):',zroeshd(i,ip,nshocks)
       enddo
      enddo

!     set the shock speed of the new shock points
      write(8,*)'shock speed of the new shock points'

      do ip=1,nshockpoints(nshocks)
       do i=1,ndim
        wsh(i,ip,nshocks)=0d0
        if(ip.eq.1) wsh(i,ip,nshocks)=wsh(i,ip1,ish1)
        write(8,'(a,i1,a,i1,a,i2,a,f16.8)')
     +   'wsh(',i,',',ip,',',nshocks,
     +   '):',wsh(i,ip,nshocks)
       enddo
      enddo

      close(8)

      return
      end
