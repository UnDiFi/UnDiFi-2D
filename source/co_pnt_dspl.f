! Compute displaced positions of the shock/discontinuity points

      subroutine co_pnt_dspl(                   ! not used
     +                       xysh,
     +                       xyshu,
     +                       xyshd,
     +                       zroeshu,
     +                       nodcodsh,
     +                       vshnor,
     +                       nshocks,
     +                       nshockpoints,
     +                       nshockedges,       ! not used
     +                       typesh,
     +                       nspecpoints,
     +                       typespecpoints,
     +                       shinspps,
     +                       ispclr)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer nspecpoints,shinspps(2,5,*),ispclr(5,*)
      integer nshocks,nshockedges(*),nshockpoints(*)
      character*1 typesh(*)
      character*5 typespecpoints(*)

!     .. array arguments ..
      double precision xysh(ndim,npshmax,*),
     &                 xyshu(ndim,npshmax,*),
     &                 xyshd(ndim,npshmax,*),
     &                 vshnor(ndim,npshmax,*),
     &                 zroeshu(ndof,npshmax,*)

      integer          nodcodsh(npshmax,*)

!     .. array arguments ..
      double precision xc(2), yc(2), xs(2), ys(2)

!     character*(*) fname

!     .. local scalars ..
      double precision xi, yi
      double precision dx,dy,dum,tx,ty
      double precision alpha,dds,f1,f2,f3
      integer i,k,n1,n2
      integer nnn,np,isppnts
      integer ip,ip1,ip2,ip3,ip4,ip5
      integer ish,ish1,ish2,ish3,ish4,ish5

      do ish = 1,nshmax
        do i = 1,npshmax
          nodcodsh(i,ish)=-99
         enddo
      enddo

!     add a second layer of shock nodes
!     place extra shock nodes on the shock normal
      do ish = 1,nshocks
        do i = 1,nshockpoints(ish)
          dx = vshnor(1,i,ish)
          dy = vshnor(2,i,ish)
          xyshu(1,i,ish)=xysh(1,i,ish)+0.5d0*eps*dx
          xyshu(2,i,ish)=xysh(2,i,ish)+0.5d0*eps*dy
          xyshd(1,i,ish)=xysh(1,i,ish)-0.5d0*eps*dx
          xyshd(2,i,ish)=xysh(2,i,ish)-0.5d0*eps*dy
          nodcodsh(i,ish)=10
        enddo
      enddo

!     correct the displacement in the boundary special points
      do isppnts = 1,nspecpoints
!       write(*,*)isppnts,typespecpoints(isppnts)
!       write(*,*)shinspps(1,1,isppnts),shinspps(2,1,isppnts)
        if(typespecpoints(isppnts).eq.'IPX'.or.
     +    typespecpoints(isppnts).eq.'OPX'.or.
     +    typespecpoints(isppnts).eq.'WPNRX')then

          ish = shinspps(1,1,isppnts)
          i   = shinspps(2,1,isppnts)-1
          ip  = 1+i*(nshockpoints(ish)-1)
          dx  = vshnor(1,ip,ish)
          dy  = vshnor(2,ip,ish)
          xyshu(1,ip,ish)=xysh(1,ip,ish)+0.5d+0*eps/dx
          xyshu(2,ip,ish)=xysh(2,ip,ish)
          xyshd(1,ip,ish)=xysh(1,ip,ish)-0.5d+0*eps/dx
          xyshd(2,ip,ish)=xysh(2,ip,ish)

        elseif(typespecpoints(isppnts).eq.'IPY'.or.
     +         typespecpoints(isppnts).eq.'OPY'.or.
     +         typespecpoints(isppnts).eq.'WPNRY')then

          ish = shinspps(1,1,isppnts)
          i   = shinspps(2,1,isppnts)-1
          ip  = 1+i*(nshockpoints(ish)-1)
          dx  = vshnor(1,ip,ish)
          dy  = vshnor(2,ip,ish)
          xyshu(1,ip,ish)=xysh(1,ip,ish)
          xyshu(2,ip,ish)=xysh(2,ip,ish)+0.5d+0*eps/dy
          xyshd(1,ip,ish)=xysh(1,ip,ish)
          xyshd(2,ip,ish)=xysh(2,ip,ish)-0.5d+0*eps/dy

        elseif(typespecpoints(isppnts).eq.'TP')then

!         incident shock
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         reflected shock
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

!         mach stem
          ish3 = shinspps(1,3,isppnts)
          i    = shinspps(2,3,isppnts)-1
          ip3  = 1+i*(nshockpoints(ish3)-1)

!         contact discontinuity
          ish4 = shinspps(1,4,isppnts)
          i    = shinspps(2,4,isppnts)-1
          ip4  = 1+i*(nshockpoints(ish4)-1)

caldo
!         determine the family of shock 1
          f1=vshnor(1, ip1, ish1)*zroeshu(4,ip1,ish1)-
     .       vshnor(2, ip1, ish1)*zroeshu(3,ip1,ish1)

          f1=-sign(1.d0,f1)
!         write(*,*)'f1:',f1

!         determine the family of shock 2
          f2=vshnor(1, ip2, ish2)*zroeshu(4,ip2,ish2)-
     .       vshnor(2, ip2, ish2)*zroeshu(3,ip2,ish2)

          f2=-sign(1.d0,f2)
!         write(*,*)'f2:',f2

!         pause
caldo

!         move incident shock point
          if(ip1.eq.1)then
            tx=xysh(1,ip1,ish1)-xysh(1,ip1-1,ish1)
            ty=xysh(2,ip1,ish1)-xysh(2,ip1-1,ish1)
          else
            tx=xysh(1,ip1-1,ish1)-xysh(1,ip1,ish1)
            ty=xysh(2,ip1-1,ish1)-xysh(2,ip1,ish1)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip1,ish1)=xyshu(1,ip1,ish1)+eps*tx
          xyshu(2,ip1,ish1)=xyshu(2,ip1,ish1)+eps*ty
          xyshd(1,ip1,ish1)=xyshd(1,ip1,ish1)+eps*tx
          xyshd(2,ip1,ish1)=xyshd(2,ip1,ish1)+eps*ty

!         move reflected shock point
          if(ip2.eq.1)then
            tx=xysh(1,ip2,ish2)-xysh(1,ip2-1,ish2)
            ty=xysh(2,ip2,ish2)-xysh(2,ip2-1,ish2)
          else
            tx=xysh(1,ip2-1,ish2)-xysh(1,ip2,ish2)
            ty=xysh(2,ip2-1,ish2)-xysh(2,ip2,ish2)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip2,ish2)=xyshu(1,ip2,ish2)+eps*tx
          xyshu(2,ip2,ish2)=xyshu(2,ip2,ish2)+eps*ty
          xyshd(1,ip2,ish2)=xyshd(1,ip2,ish2)+eps*tx
          xyshd(2,ip2,ish2)=xyshd(2,ip2,ish2)+eps*ty

!         move mach stem point
          if(ip3.eq.1)then
            tx=xysh(1,ip3,ish3)-xysh(1,ip3-1,ish3)
            ty=xysh(2,ip3,ish3)-xysh(2,ip3-1,ish3)
          else
            tx=xysh(1,ip3-1,ish3)-xysh(1,ip3,ish3)
            ty=xysh(2,ip3-1,ish3)-xysh(2,ip3,ish3)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip3,ish3)=xyshu(1,ip3,ish3)+eps*tx
          xyshu(2,ip3,ish3)=xyshu(2,ip3,ish3)+eps*ty
          xyshd(1,ip3,ish3)=xyshd(1,ip3,ish3)+eps*tx
          xyshd(2,ip3,ish3)=xyshd(2,ip3,ish3)+eps*ty

!         move contact discontinuity point
          if(ip4.eq.1)then
            tx=xysh(1,ip4,ish4)-xysh(1,ip4-1,ish4)
            ty=xysh(2,ip4,ish4)-xysh(2,ip4-1,ish4)
          else
            tx=xysh(1,ip4-1,ish4)-xysh(1,ip4,ish4)
            ty=xysh(2,ip4-1,ish4)-xysh(2,ip4,ish4)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip4,ish4)=xyshu(1,ip4,ish4)+eps*tx
          xyshu(2,ip4,ish4)=xyshu(2,ip4,ish4)+eps*ty
          xyshd(1,ip4,ish4)=xyshd(1,ip4,ish4)+eps*tx
          xyshd(2,ip4,ish4)=xyshd(2,ip4,ish4)+eps*ty

caldo

!     if incident shock and reflected shock belong to opposite family
!     overlap incident shock upstream point with mach stem upstream point
!     write(*,*)'f1,f2:',f1,f2
!     pause
      if(f1*f2.lt.0.)then
          do k=1,ndim
           dum=0.5*(xyshu(k,ip1,ish1)+xyshu(k,ip3,ish3))
           xyshu(k,ip1,ish1)=dum
           xyshu(k,ip3,ish3)=dum
          enddo

!         overlap incident shock downstream point with reflected shock upstream point
          do k=1,ndim
           dum=0.5*(xyshd(k,ip1,ish1)+xyshu(k,ip2,ish2))
           xyshd(k,ip1,ish1)=dum
           xyshu(k,ip2,ish2)=dum
          enddo
      else
!     if otherwise incident shock and reflected shock belong to the same family
!     overlap incident shock downstream point with mach stem upstream point
         do k=1,ndim
           dum=0.5*(xyshd(k,ip1,ish1)+xyshu(k,ip3,ish3))
           xyshd(k,ip1,ish1)=dum
           xyshu(k,ip3,ish3)=dum
          enddo

!         overlap incident shock uptream point with reflected shock upstream point
          do k=1,ndim
           dum=0.5*(xyshu(k,ip1,ish1)+xyshu(k,ip2,ish2))
           xyshu(k,ip1,ish1)=dum
           xyshu(k,ip2,ish2)=dum
          enddo
       endif

caldo

!         overlap mach stem downstream point with downstream point
!         of the contact discontinuity (be careful to the orientation
!         of the contact discontinuity)
          do k=1,ndim
           dum=0.5*(xyshd(k,ip3,ish3)+xyshd(k,ip4,ish4))
           xyshd(k,ip3,ish3)=dum
           xyshd(k,ip4,ish4)=dum
          enddo

!         overlap reflected shock downstream point with upstream point
!         of the contact discontinuity (be careful to the orientation
!         of the contact discontinuity)
          do k=1,ndim
           dum=0.5*(xyshd(k,ip2,ish2)+xyshu(k,ip4,ish4))
           xyshd(k,ip2,ish2)=dum
           xyshu(k,ip4,ish4)=dum
          enddo

          elseif(typespecpoints(isppnts).eq.'RRX')then

!         incident shock
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         reflected shock
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

          xyshu(1,ip1,ish1)=xysh(1,ip1,ish1)+0.5d+0*eps/dx
          xyshu(2,ip1,ish1)=xysh(2,ip1,ish1)
          xyshd(1,ip1,ish1)=xysh(1,ip1,ish1)-0.5d+0*eps/dx
          xyshd(2,ip1,ish1)=xysh(2,ip1,ish1)

          xyshu(1,ip2,ish2)=xysh(1,ip2,ish2)+0.5d+0*eps/dx
          xyshu(2,ip2,ish2)=xysh(2,ip2,ish2)
          xyshd(1,ip2,ish2)=xysh(1,ip2,ish2)-0.5d+0*eps/dx
          xyshd(2,ip2,ish2)=xysh(2,ip2,ish2)

!         overlap incident shock downstream point with incident shock upstream point
!         xyshd(1,ip1,ish1)=xysh(1,ip1,ish1)
!         xyshd(2,ip1,ish1)=eps
!         xyshu(1,ip2,ish2)=xysh(1,ip2,ish2)
!         xyshu(2,ip2,ish2)=eps

          xc(1)=xyshd(1,ip1,ish1)
          yc(1)=xyshd(2,ip1,ish1)

          if(ip1.eq.1)then
           xc(2)=xyshd(1,ip1+1,ish1)
           yc(2)=xyshd(2,ip1+1,ish1)
          else
           xc(2)=xyshd(1,ip1-1,ish1)
           yc(2)=xyshd(2,ip1-1,ish1)
          endif

          xs(1)=xyshu(1,ip2,ish2)
          ys(1)=xyshu(2,ip2,ish2)

          if(ip2.eq.1)then
           xs(2)=xyshu(1,ip2+1,ish2)
           ys(2)=xyshu(2,ip2+1,ish2)
          else
           xs(2)=xyshu(1,ip2-1,ish2)
           ys(2)=xyshu(2,ip2-1,ish2)
          endif

!         find interpolation point
          call co_intr_pnt(xi,yi,xc,yc,xs,ys)

          xyshd(1,ip1,ish1)=xi
          xyshd(2,ip1,ish1)=yi

          xyshu(1,ip2,ish2)=xi
          xyshu(2,ip2,ish2)=yi

!         xyshd(1,ip1,ish1)=xysh(1,ip1,ish1)
!         xyshd(2,ip1,ish1)=eps
!         xyshu(1,ip2,ish2)=xysh(1,ip2,ish2)
!         xyshu(2,ip2,ish2)=eps

          elseif(typespecpoints(isppnts).eq.'RR')then

!         incident shock
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         reflected shock
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

!         overlap incident shock downstream point with reflected shock upstream point
          xc(1)=xyshd(1,ip1,ish1)
          yc(1)=xyshd(2,ip1,ish1)

          if(ip1.eq.1)then
           xc(2)=xyshd(1,ip1+1,ish1)
           yc(2)=xyshd(2,ip1+1,ish1)
          else
           xc(2)=xyshd(1,ip1-1,ish1)
           yc(2)=xyshd(2,ip1-1,ish1)
          endif

          xs(1)=xyshu(1,ip2,ish2)
          ys(1)=xyshu(2,ip2,ish2)

          if(ip2.eq.1)then
           xs(2)=xyshu(1,ip2+1,ish2)
           ys(2)=xyshu(2,ip2+1,ish2)
          else
           xs(2)=xyshu(1,ip2-1,ish2)
           ys(2)=xyshu(2,ip2-1,ish2)
          endif

!         find interpolation point
          call co_intr_pnt(xi,yi,xc,yc,xs,ys)

          xyshd(1,ip1,ish1)=xi
          xyshd(2,ip1,ish1)=yi

          xyshu(1,ip2,ish2)=xi
          xyshu(2,ip2,ish2)=yi

         elseif(typespecpoints(isppnts).eq.'QP')then

!         incident shock 1
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         reflected shock 1
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

!         incident shock 2
          ish3 = shinspps(1,3,isppnts)
          i    = shinspps(2,3,isppnts)-1
          ip3  = 1+i*(nshockpoints(ish3)-1)

!         reflected shock 2
          ish4 = shinspps(1,4,isppnts)
          i    = shinspps(2,4,isppnts)-1
          ip4  = 1+i*(nshockpoints(ish4)-1)

!         contact discontinuity
          ish5 = shinspps(1,5,isppnts)
          i    = shinspps(2,5,isppnts)-1
          ip5  = 1+i*(nshockpoints(ish5)-1)

!         determine family of shock 1
          f1=vshnor(1, ip1, ish1)*zroeshu(4,ip1,ish1)-
     .       vshnor(2, ip1, ish1)*zroeshu(3,ip1,ish1)

          f1=-sign(1.d0,f1)

!         determine family of shock 3
          f3=vshnor(1, ip3, ish3)*zroeshu(4,ip3,ish3)-
     .       vshnor(2, ip3, ish3)*zroeshu(3,ip3,ish3)

          f3=-sign(1.d0,f3)

!         move incident shock 1 point
          if(ip1.eq.1)then
            tx=xysh(1,ip1,ish1)-xysh(1,ip1-1,ish1)
            ty=xysh(2,ip1,ish1)-xysh(2,ip1-1,ish1)
          else
            tx=xysh(1,ip1-1,ish1)-xysh(1,ip1,ish1)
            ty=xysh(2,ip1-1,ish1)-xysh(2,ip1,ish1)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip1,ish1)=xyshu(1,ip1,ish1)+eps*tx
          xyshu(2,ip1,ish1)=xyshu(2,ip1,ish1)+eps*ty

          xyshd(1,ip1,ish1)=xyshd(1,ip1,ish1)+eps*tx
          xyshd(2,ip1,ish1)=xyshd(2,ip1,ish1)+eps*ty

!         move reflected shock 1 point
          if(ip2.eq.1)then
            tx=xysh(1,ip2,ish2)-xysh(1,ip2-1,ish2)
            ty=xysh(2,ip2,ish2)-xysh(2,ip2-1,ish2)
          else
            tx=xysh(1,ip2-1,ish2)-xysh(1,ip2,ish2)
            ty=xysh(2,ip2-1,ish2)-xysh(2,ip2,ish2)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip2,ish2)=xyshu(1,ip2,ish2)+eps*tx
          xyshu(2,ip2,ish2)=xyshu(2,ip2,ish2)+eps*ty

          xyshd(1,ip2,ish2)=xyshd(1,ip2,ish2)+eps*tx
          xyshd(2,ip2,ish2)=xyshd(2,ip2,ish2)+eps*ty

!         move incident shock 2 point
          if(ip3.eq.1)then
            tx=xysh(1,ip3,ish3)-xysh(1,ip3-1,ish3)
            ty=xysh(2,ip3,ish3)-xysh(2,ip3-1,ish3)
          else
            tx=xysh(1,ip3-1,ish3)-xysh(1,ip3,ish3)
            ty=xysh(2,ip3-1,ish3)-xysh(2,ip3,ish3)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip3,ish3)=xyshu(1,ip3,ish3)+eps*tx
          xyshu(2,ip3,ish3)=xyshu(2,ip3,ish3)+eps*ty

          xyshd(1,ip3,ish3)=xyshd(1,ip3,ish3)+eps*tx
          xyshd(2,ip3,ish3)=xyshd(2,ip3,ish3)+eps*ty

!         move reflected shock 2 point
          if(ip4.eq.1)then
            tx=xysh(1,ip4,ish4)-xysh(1,ip4-1,ish4)
            ty=xysh(2,ip4,ish4)-xysh(2,ip4-1,ish4)
          else
            tx=xysh(1,ip4-1,ish4)-xysh(1,ip4,ish4)
            ty=xysh(2,ip4-1,ish4)-xysh(2,ip4,ish4)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip4,ish4)=xyshu(1,ip4,ish4)+eps*tx
          xyshu(2,ip4,ish4)=xyshu(2,ip4,ish4)+eps*ty

          xyshd(1,ip4,ish4)=xyshd(1,ip4,ish4)+eps*tx
          xyshd(2,ip4,ish4)=xyshd(2,ip4,ish4)+eps*ty

!         move contact discontinuity point
          if(ip5.eq.1)then
            tx=xysh(1,ip5,ish5)-xysh(1,ip5-1,ish5)
            ty=xysh(2,ip5,ish5)-xysh(2,ip5-1,ish5)
          else
            tx=xysh(1,ip5-1,ish5)-xysh(1,ip5,ish5)
            ty=xysh(2,ip5-1,ish5)-xysh(2,ip5,ish5)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip5,ish5)=xyshu(1,ip5,ish5)+eps*tx
          xyshu(2,ip5,ish5)=xyshu(2,ip5,ish5)+eps*ty

          xyshd(1,ip5,ish5)=xyshd(1,ip5,ish5)+eps*tx
          xyshd(2,ip5,ish5)=xyshd(2,ip5,ish5)+eps*ty

!       if incident shocks belong to opposite family
!       overlap incident shock 1 upstream point with
!       incident shock 2 upstream point

!       if shock 1

        if(f1.gt.0.d0)then
         xc(1)=xyshu(1,ip1,ish1)
         yc(1)=xyshu(2,ip1,ish1)

         if(ip1.eq.1)then
          xc(2)=xyshu(1,ip1+1,ish1)
          yc(2)=xyshu(2,ip1+1,ish1)
         else
          xc(2)=xyshu(1,ip1-1,ish1)
          yc(2)=xyshu(2,ip1-1,ish1)
         endif
        else
         xc(1)=xyshd(1,ip1,ish1)
         yc(1)=xyshd(2,ip1,ish1)

         if(ip1.eq.1)then
          xc(2)=xyshd(1,ip1+1,ish1)
          yc(2)=xyshd(2,ip1+1,ish1)
         else
          xc(2)=xyshd(1,ip1-1,ish1)
          yc(2)=xyshd(2,ip1-1,ish1)
         endif
        endif

        if(f3.lt.0.d0)then
         xs(1)=xyshu(1,ip3,ish3)
         ys(1)=xyshu(2,ip3,ish3)

         if(ip3.eq.1)then
          xs(2)=xyshu(1,ip3+1,ish3)
          ys(2)=xyshu(2,ip3+1,ish3)
         else
          xs(2)=xyshu(1,ip3-1,ish3)
          ys(2)=xyshu(2,ip3-1,ish3)
         endif

        else
         xs(1)=xyshd(1,ip3,ish3)
         ys(1)=xyshd(2,ip3,ish3)

         if(ip3.eq.1)then
          xs(2)=xyshd(1,ip3+1,ish3)
          ys(2)=xyshd(2,ip3+1,ish3)
         else
          xs(2)=xyshd(1,ip3-1,ish3)
          ys(2)=xyshd(2,ip3-1,ish3)
         endif
       endif

!       find intersection point
        call co_intr_pnt(xi,yi,xc,yc,xs,ys)

        if(f1.gt.0.d0)then
         xyshu(1,ip1,ish1)=xi
         xyshu(2,ip1,ish1)=yi
        else
         xyshd(1,ip1,ish1)=xi
         xyshd(2,ip1,ish1)=yi
        endif

        if(f3.lt.0.d0)then
         xyshu(1,ip3,ish3)=xi
         xyshu(2,ip3,ish3)=yi
        else
         xyshd(1,ip3,ish3)=xi
         xyshd(2,ip3,ish3)=yi
        endif

!       overlap incident shock 1 downstream point with
!       reflected shock 1 upstream point
!       do k=1,ndim
!        dum=0.5*(xyshd(k,ip1,ish1)+xyshu(k,ip2,ish2))
!        xyshd(k,ip1,ish1)=dum
!        xyshu(k,ip2,ish2)=dum
!       enddo

        if(f1.gt.0.d0)then
          xc(1)=xyshd(1,ip1,ish1)
          yc(1)=xyshd(2,ip1,ish1)
          if(ip1.eq.1)then
           xc(2)=xyshd(1,ip1+1,ish1)
           yc(2)=xyshd(2,ip1+1,ish1)
          else
           xc(2)=xyshd(1,ip1-1,ish1)
           yc(2)=xyshd(2,ip1-1,ish1)
          endif
        else
          xc(1)=xyshu(1,ip1,ish1)
          yc(1)=xyshu(2,ip1,ish1)
          if(ip1.eq.1)then
           xc(2)=xyshu(1,ip1+1,ish1)
           yc(2)=xyshu(2,ip1+1,ish1)
          else
           xc(2)=xyshu(1,ip1-1,ish1)
           yc(2)=xyshu(2,ip1-1,ish1)
          endif
        endif

        xs(1)=xyshu(1,ip2,ish2)
        ys(1)=xyshu(2,ip2,ish2)

        if(ip2.eq.1)then
         xs(2)=xyshu(1,ip2+1,ish2)
         ys(2)=xyshu(2,ip2+1,ish2)
        else
         xs(2)=xyshu(1,ip2-1,ish2)
         ys(2)=xyshu(2,ip2-1,ish2)
        endif

!       find intersection point
        call co_intr_pnt(xi,yi,xc,yc,xs,ys)

        if(f1.gt.0.d0)then
         xyshd(1,ip1,ish1)=xi
         xyshd(2,ip1,ish1)=yi
        else
         xyshu(1,ip1,ish1)=xi
         xyshu(2,ip1,ish1)=yi
        endif

        xyshu(1,ip2,ish2)=xi
        xyshu(2,ip2,ish2)=yi

!       overlap incident shock 2 downstream point with
!       reflected shock 2 upstream point
!       do k=1,ndim
!        dum=0.5*(xyshd(k,ip3,ish3)+xyshu(k,ip4,ish4))
!        xyshd(k,ip3,ish3)=dum
!        xyshu(k,ip4,ish4)=dum
!       enddo

        if(f3.gt.0.d+0)then
         xc(1)=xyshu(1,ip3,ish3)
         yc(1)=xyshu(2,ip3,ish3)

         if(ip3.eq.1)then
          xc(2)=xyshu(1,ip3+1,ish3)
          yc(2)=xyshu(2,ip3+1,ish3)
         else
          xc(2)=xyshu(1,ip3-1,ish3)
          yc(2)=xyshu(2,ip3-1,ish3)
         endif
        else
         xc(1)=xyshd(1,ip3,ish3)
         yc(1)=xyshd(2,ip3,ish3)

         if(ip3.eq.1)then
          xc(2)=xyshd(1,ip3+1,ish3)
          yc(2)=xyshd(2,ip3+1,ish3)
         else
          xc(2)=xyshd(1,ip3-1,ish3)
          yc(2)=xyshd(2,ip3-1,ish3)
         endif
        endif

        xs(1)=xyshu(1,ip4,ish4)
        ys(1)=xyshu(2,ip4,ish4)

        if(ip4.eq.1)then
         xs(2)=xyshu(1,ip4+1,ish4)
         ys(2)=xyshu(2,ip4+1,ish4)
        else
         xs(2)=xyshu(1,ip4-1,ish4)
         ys(2)=xyshu(2,ip4-1,ish4)
        endif

!       find intersection point
        call co_intr_pnt(xi,yi,xc,yc,xs,ys)

        if(f3.lt.0.d0)then
         xyshd(1,ip3,ish3)=xi
         xyshd(2,ip3,ish3)=yi
        else
         xyshu(1,ip3,ish3)=xi
         xyshu(2,ip3,ish3)=yi
        endif

        xyshu(1,ip4,ish4)=xi
        xyshu(2,ip4,ish4)=yi

!       overlap reflected shock 1 downstream point with
!       contact discontinuity upstream point
!       do k=1,ndim
!        dum=0.5*(xyshd(k,ip2,ish2)+xyshu(k,ip5,ish5))
!        xyshd(k,ip2,ish2)=dum
!        xyshu(k,ip5,ish5)=dum
!       enddo

        xc(1)=xyshd(1,ip2,ish2)
        yc(1)=xyshd(2,ip2,ish2)

        if(ip2.eq.1)then
         xc(2)=xyshd(1,ip2+1,ish2)
         yc(2)=xyshd(2,ip2+1,ish2)
        else
         xc(2)=xyshd(1,ip2-1,ish2)
         yc(2)=xyshd(2,ip2-1,ish2)
        endif

        xs(1)=xyshu(1,ip5,ish5)
        ys(1)=xyshu(2,ip5,ish5)

        if(ip5.eq.1)then
         xs(2)=xyshu(1,ip5+1,ish5)
         ys(2)=xyshu(2,ip5+1,ish5)
        else
         xs(2)=xyshu(1,ip5-1,ish5)
         ys(2)=xyshu(2,ip5-1,ish5)
        endif

!       find intersection point
        call co_intr_pnt(xi,yi,xc,yc,xs,ys)

        xyshd(1,ip2,ish2)=xi
        xyshd(2,ip2,ish2)=yi

        xyshu(1,ip5,ish5)=xi
        xyshu(2,ip5,ish5)=yi

!       overlap reflected shock 2 downstream point with
!       contact discontinuity downstream point
!       do k=1,ndim
!        dum=0.5*(xyshd(k,ip4,ish4)+xyshd(k,ip5,ish5))
!        xyshd(k,ip4,ish4)=dum
!       xyshd(k,ip5,ish5)=dum
!       enddo

        xc(1)=xyshd(1,ip4,ish4)
        yc(1)=xyshd(2,ip4,ish4)

        if(ip4.eq.1)then
         xc(2)=xyshd(1,ip4+1,ish4)
         yc(2)=xyshd(2,ip4+1,ish4)
        else
         xc(2)=xyshd(1,ip4-1,ish4)
         yc(2)=xyshd(2,ip4-1,ish4)
        endif

        xs(1)=xyshd(1,ip5,ish5)
        ys(1)=xyshd(2,ip5,ish5)

        if(ip5.eq.1)then
         xs(2)=xyshd(1,ip5+1,ish5)
         ys(2)=xyshd(2,ip5+1,ish5)
        else
         xs(2)=xyshd(1,ip5-1,ish5)
         ys(2)=xyshd(2,ip5-1,ish5)
        endif

!       find intersection point
        call co_intr_pnt(xi,yi,xc,yc,xs,ys)

        xyshd(1,ip4,ish4)=xi
        xyshd(2,ip4,ish4)=yi

        xyshd(1,ip5,ish5)=xi
        xyshd(2,ip5,ish5)=yi

         elseif(typespecpoints(isppnts).eq.'TE')then

!         shock 1
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         contact discontinuity
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)

!         shock 2
          ish3 = shinspps(1,3,isppnts)
          i    = shinspps(2,3,isppnts)-1
          ip3  = 1+i*(nshockpoints(ish3)-1)

!         move contact discontinuity point
          if(ip1.eq.1)then
            tx=xysh(1,ip2,ish2)-xysh(1,ip2-1,ish2)
            ty=xysh(2,ip2,ish2)-xysh(2,ip2-1,ish2)
          else
            tx=xysh(1,ip2-1,ish2)-xysh(1,ip2,ish2)
            ty=xysh(2,ip2-1,ish2)-xysh(2,ip2,ish2)
          endif
          dum=sqrt(tx**2+ty**2)
          tx=tx/dum
          ty=ty/dum
          xyshu(1,ip2,ish2)=xyshu(1,ip2,ish2)+eps*tx
          xyshu(2,ip2,ish2)=xyshu(2,ip2,ish2)+eps*ty

          xyshd(1,ip2,ish2)=xyshd(1,ip2,ish2)+eps*tx
          xyshd(2,ip2,ish2)=xyshd(2,ip2,ish2)+eps*ty

!         move upstream and downstream shock 1 point
!         xyshu(1,ip1,ish1)=xysh(1,ip2,ish2)
!         xyshu(2,ip1,ish1)=xysh(2,ip2,ish2)

          xyshd(1,ip1,ish1)=xyshd(1,ip2,ish2)
          xyshd(2,ip1,ish1)=xyshd(2,ip2,ish2)

!         move upstream and downstream shock 2 point
!         xyshu(1,ip3,ish3)=xysh(1,ip2,ish2)
!         xyshu(2,ip3,ish3)=xysh(2,ip2,ish2)

          xyshd(1,ip3,ish3)=xyshu(1,ip2,ish2)
          xyshd(2,ip3,ish3)=xyshu(2,ip2,ish2)

         elseif(typespecpoints(isppnts).eq.'EP')then

           ish = shinspps(1,1,isppnts)
           i    = shinspps(2,1,isppnts)-1
           ip  = 1+i*(nshockpoints(ish)-1)

!          overlap downstream with upstream point
           xyshd(1,ip,ish)=xysh(1,ip,ish)
           xyshd(2,ip,ish)=xysh(2,ip,ish)
           xyshu(1,ip,ish)=xysh(1,ip,ish)
           xyshu(2,ip,ish)=xysh(2,ip,ish)

         elseif(typespecpoints(isppnts).eq.'C')then

!         incident shock 1
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

!         reflected shock 1
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2 = 1+i*(nshockpoints(ish2)-1)

!         overlap downstream and upstream points with
!         those of the connection point
          xyshd(1,ip2,ish2)=xyshd(1,ip1,ish1)
          xyshd(2,ip2,ish2)=xyshd(2,ip1,ish1)
          xyshu(1,ip2,ish2)=xyshu(1,ip1,ish1)
          xyshu(2,ip2,ish2)=xyshu(2,ip1,ish1)

        elseif(typespecpoints(isppnts).eq.'SP')then

           ish = shinspps(1,1,isppnts)
           i    = shinspps(2,1,isppnts)-1
           ip  = 1+i*(nshockpoints(ish)-1)

!         overlap downstream point with upstream point
          xyshd(1,ip,ish)=xysh(1,ip,ish)
          xyshd(2,ip,ish)=xysh(2,ip,ish)
          xyshu(1,ip,ish)=xysh(1,ip,ish)
          xyshu(2,ip,ish)=xysh(2,ip,ish)

         elseif(typespecpoints(isppnts).eq.'FWP')then

         elseif(typespecpoints(isppnts).eq.'PC')then

         else

           write(8,*)'condition not implemented'
           write(*,*)'condition not implemented'
           stop

         endif
      enddo

      return
      end

!     find interpolation point
      subroutine co_intr_pnt(xi,yi,xc,yc,xs,ys)

      integer nn
      parameter (nn=2)
      double precision a(nn,nn), b(nn),x(nn)
      double precision xi, yi
      double precision xc(2), yc(2), xs(2), ys(2)

      rdshp=-1.0

      a(1,1)=(ys(2)-ys(1))
      a(1,2)=(xs(1)-xs(2))
      b(1)=xs(2)*(ys(2)-ys(1))+ys(2)*(xs(1)-xs(2))
!     b(1)=xs(1)*(ys(2)-ys(1))+ys(1)*(xs(1)-xs(2))
      a(2,1)=(yc(2)-yc(1))
      a(2,2)=(xc(1)-xc(2))
      b(2)=xc(2)*(yc(2)-yc(1))+yc(2)*(xc(1)-xc(2))
!     b(2)=xc(1)*(yc(2)-yc(1))+yc(1)*(xc(1)-xc(2))

      call solg(nn,nn,a,b,x)

      xi=x(1)
      yi=x(2)

      return
      end
