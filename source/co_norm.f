! Compute the normal unit vectors to the shocks and discontinuities in shock/discontinuity points

      subroutine co_norm (xysh,
     +                    zroeshu, !upstream
     +                    zroesh,  !downstream
     +                    vshnor,
     +                    nshocks,
     +                    nshockpoints,
     +                    typesh,
     +                    nspecpoints,
     +                    typespecpoints,
     +                    shinspps,
     +                    ispclr,
     +                    ia,
     +                    ja,
     +                    iclr,
     +                    nclr,
     +                    corg)

      implicit none
      include 'paramt.h'

      integer nshocks,nshockpoints(*),nspecpoints,shinspps(2,5,*)
      integer ispclr(5,*)
      integer i,j,j2,k,kp1,depip1,depim1,ii

      integer nclr
      integer ia(*),ja(*),iclr(nclr)
      double precision corg(ndim,*)

      integer clr,bbgn,bend

      double precision xysh  (ndim,  npshmax,*)
      double precision vshnor(ndim,npshmax,*)
      double precision zroesh(ndof,npshmax,*)

! vale
      double precision zroeshu(ndof,npshmax,*)
      double precision pu, pd
! vale

      double precision xi,yi,xj,yj,ush,vsh,tau,dum
      double precision x1,y1,x2,y2
      double precision dumx1,dumy1,dumx2,dumy2
      double precision um,vm,rom,am,pm,help,mm,alpham,thetam
      double precision uj,vj,roj,aj,pj,     mj,alphaj,thetaj
      double precision uv,vv,rov,av,pv,     mv,alphav,thetav
      double precision tauxim1,tauyim1,tauxip1,tauyip1,taux,tauy
      double precision xj2,yj2,tauxip2,tauyip2,tauxim2,tauyim2
      double precision tauxjp2,tauyjp2,tauxjm2,tauyjm2
      double precision lp12,lm12,ui,vi,nx2,ny2,nx3,ny3,nx4,ny4,nx1,ny1
      double precision dum1,dum2,lp1,lp2,lm1,lm2,lp22,lm22
      double precision a,b,c, nx,ny,dist
      integer shp_dpndnc,dcp_dpndnc,ish,isppnts
      integer ip,ip1,ip2,ip3,ip4,ish1,ish2,ish3,ish4
      character*1 typesh(*)
      character*5 typespecpoints(*)

!     .. arrays in common ..
      double precision dstak(1)
      integer istak(1)
      common/cstak/dstak

!     .. equivalences ..
      equivalence (dstak(1),istak(1))

!     input:
!     -----
!     xysh   x,y coords of the shock points
!     zroesh status variables
!     nshockpoints number of shock points
!
!     output:
!     ------
!     vshnor x,y components of the unit normal to the shock

!     open log file
      open(8,file='log/co_norm.log')

      write(8,*)'n. of shocks:',nshocks
      do ish=1,nshocks
      write(8,*)'n. of shock points:',nshockpoints(ish)
      write(8,*)'shock-point coordinates, downtream'
       do i = 1, nshockpoints(ish)
         write(8,*)i,xysh(1,i,ish),xysh(2,i,ish),
     &               zroesh(1,i,ish),zroesh(2,i,ish),
     &               zroesh(3,i,ish),zroesh(4,i,ish)
       enddo
      enddo

!     normals computation for each shock
      do ish=1,nshocks
       do i = 1, nshockpoints(ish)

        ush=0.d0
        vsh=0.d0

        xi  = xysh(1,i,ish)
        yi  = xysh(2,i,ish)

!       upward points
        j=i+1

!       recover coordinates
        xj=xysh(1,j,ish)
        yj=xysh(2,j,ish)

        j2=i+2
        xj2=xysh(1,j2,ish)
        yj2=xysh(2,j2,ish)

!       tangent vectors computation
        tauxip1=xj-xi
        tauyip1=yj-yi

        tauxip2=xj2-xi
        tauyip2=yj2-yi
        tauxjp2=xj2-xj
        tauyjp2=yj2-yj

!       recover state
        uj=zroesh(3,j,ish)/zroesh(1,j,ish)
        vj=zroesh(4,j,ish)/zroesh(1,j,ish)
        roj=zroesh(1,j,ish)*zroesh(1,j,ish)
        help=zroesh(3,j,ish)**2+zroesh(4,j,ish)**2
        pj = gm1/ga*( zroesh(1,j,ish)*zroesh(2,j,ish)-0.5d0*help)
        aj=sqrt(ga*pj/roj)

!       dependency evaluation
!       if it is the last points, then compute backward
        if(i.ne.1.and.i.ne.nshockpoints(ish))then
         if(typesh(ish).eq.'S')
     +      depip1=shp_dpndnc(xi,yi,ush,vsh,xj,yj,uj,vj,aj)
         if(typesh(ish).eq.'D')
     +      depip1=dcp_dpndnc(xi,yi,ush,vsh,xj,yj,uj,vj,aj)
        elseif(i.eq.nshockpoints(ish))then
         depip1=0
         depim1=1
!        goto 999
        endif

!       backward points
        j=i-1

!       recover coordinates
        xj=xysh(1,j,ish)
        yj=xysh(2,j,ish)

        j2=i-2
        xj2=xysh(1,j2,ish)
        yj2=xysh(2,j2,ish)

!       tangent vectors computation
        tauxim1=xi-xj
        tauyim1=yi-yj

        tauxim2=xi-xj2
        tauyim2=yi-yj2
        tauxjm2=xj-xj2
        tauyjm2=yj-yj2

!       recover state
        uj=zroesh(3,j,ish)/zroesh(1,j,ish)
        vj=zroesh(4,j,ish)/zroesh(1,j,ish)
        roj=zroesh(1,j,ish)*zroesh(1,j,ish)
        help=zroesh(3,j,ish)**2+zroesh(4,j,ish)**2
        pj = gm1/ga*( zroesh(1,j,ish)*zroesh(2,j,ish)-0.5d0*help)
        aj=sqrt(ga*pj/roj)

!       dependency evaluation
!       if it is the first point, then compute upward
        if(i.ne.1.and.i.ne.nshockpoints(ish))then
         if(typesh(ish).eq.'S')
     +      depim1=shp_dpndnc(xi,yi,ush,vsh,xj,yj,uj,vj,aj)
         if(typesh(ish).eq.'D')
     +      depim1=dcp_dpndnc(xi,yi,ush,vsh,xj,yj,uj,vj,aj)
        elseif(I.eq.1)then
         depim1=0
         depip1=1
!        goto 999
        endif

999    continue

!       tangent versors computation
        lp12=tauxip1**2+tauyip1**2
        lm12=tauxim1**2+tauyim1**2
        lp1=sqrt(lp12)
        lm1=sqrt(lm12)

        lm22=tauxjm2**2+tauyjm2**2
        lp22=tauxjp2**2+tauyjp2**2

        lp2=sqrt(lp22)
        lm2=sqrt(lm22)

        if(i.eq.1)then
          depim1=0
          depip1=1
          lm12=1.0
        endif

        if(i.eq.nshockpoints(ish))then
          depim1=1
          depip1=0
          lp12=1.0
        endif

!       TODO: verify the case when the point has no dependence
        if(depim1.eq.0.and.depip1.eq.0)then
            depim1=1
            depip1=1
        endif

        taux=(tauxim1*lp12+tauxip1*lm12)
        tauy=(tauyim1*lp12+tauyip1*lm12)
        if(depim1*depip1.eq.0.0)then
          if(depim1.eq.1)then
           taux=tauxim1*(lm1+lm2)**2-tauxim2*lm12
           tauy=tauyim1*(lm1+lm2)**2-tauyim2*lm12
          endif
          if(depip1.eq.1)then
           taux=tauxip1*(lp1+lp2)**2-tauxip2*lp12
           tauy=tauyip1*(lp1+lp2)**2-tauyip2*lp12
          endif
        endif

        write(8,*)I,'tau',taux,tauy,depim1,depip1

        tau=sqrt(taux*taux+tauy*tauy)
        taux=taux/tau
        tauy=tauy/tau

!       calculate and assign of the normal versor
!       TODO: check orientation downstream->upstream of the versor
!       if(ish.eq.2) write(*,*)'taux,tauy:',i,taux,tauy,
!    &                          xysh(1,i,ish),xysh(2,i,ish)

         vshnor(1,i,ish)=tauy
         vshnor(2,i,ish)=-taux

      enddo

      enddo

!     compute and/or assign the normal in the union point of two shocks
!     in this case, the normal of the first point of the second shock
!     coincides ith the normal of the last point of the first shock
!
!          vshnor(1,nshockpoints(1),1)=vshnor(1,1,2)
!          vshnor(2,nshockpoints(1),1)=vshnor(2,1,2)
!
!     check the correct orientation of the normals
!     the orientation is downstream toward upstream
!     verify if u * n <= 0 otherwise change direction of n
!     do ish=1,nshocks
!      if(typesh(ish).eq.'s')then
!      do i=1,nshockpoints(ish)
!       ui=zroesh(3,i,ish)/zroesh(1,i,ish)
!       vi=zroesh(4,i,ish)/zroesh(1,i,ish)
!       dum=ui*vshnor(1,i,ish)+vi*vshnor(2,i,ish)
!       if(dum.gt.0.d0)then
!         vshnor(1,i,ish)= -vshnor(1,i,ish)
!         vshnor(2,i,ish)= -vshnor(2,i,ish)
!       endif
!      enddo
!      endif
!     enddo

      do ish=1,nshocks
       if(typesh(ish).eq.'S')then
       ii=0
       do i=1,nshockpoints(ish)
        ui=zroesh(3,i,ish)/zroesh(1,i,ish)
        vi=zroesh(4,i,ish)/zroesh(1,i,ish)
        dum=ui*vshnor(1,i,ish)+vi*vshnor(2,i,ish)
        if(dum.gt.0.d0)then
          ii=ii+1
        endif
       enddo
       if(ii.lt.nshockpoints(ish)/2.)goto 234
       do i=1,nshockpoints(ish)
! vale   vshnor(1,i,ish)= -vshnor(1,i,ish)
! vale   vshnor(2,i,ish)= -vshnor(2,i,ish)
         vshnor(1,i,ish)= -vshnor(1,i,ish)
         vshnor(2,i,ish)= -vshnor(2,i,ish)
       enddo
234    continue

       endif
      enddo

!     in the case of triple points, it forces the orientation of the
!     normal to the contact discontinuity such that it forms an
!     angle > 90° with the normal to the mach stem
      do isppnts=1,nspecpoints

!     if the shock point is floating on the wall directed along x-axis

       if(typespecpoints(isppnts).eq.'WPNRX')then
          write(8,*)'correction for WPNRX point '

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          vshnor(1,ip1,ish1)=vshnor(1,ip1,ish1)/abs(vshnor(1,ip1,ish1))
          vshnor(2,ip1,ish1)=0.

!      if the shock point is floating on the boundary with a specific colour

       elseif(typespecpoints(isppnts).eq.'FWP') then
          write(8,*)'correction for FWP point '

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          xi= xysh(1,ip1,ish1)
          yi= xysh(2,ip1,ish1)

!         set colour of the coundary on which the shocs point moves
          clr=5

          do clr=1,nclr
           if(iclr(clr).eq.ispclr(1,isppnts))exit
          enddo

          bbgn=ia(clr)
          bend=ia(clr+1)-1
          do j = bbgn, bend-1
             k   = ja(j)
             kp1 = ja(j+1)
             x1=corg(1,k)
             y1=corg(2,k)
             x2=corg(1,kp1)
             y2=corg(2,kp1)

             dumx1=x1
             dumy1=y1
             dumx2=x2
             dumy2=y2

             if(dumx2.lt.dumx1)then
              dum=dumx1
              dumx1=dumx2
              dumx2=dum
             endif

             if(dumy2.lt.dumy1)then
              dum=dumy1
              dumy1=dumy2
              dumy2=dum
             endif

             if(xi.le.dumx2.and.
     +          xi.ge.dumx1.and.
     +          yi.le.dumy2.and.
     +          yi.ge.dumy1)then

              taux=x2-x1
              tauy=y2-y1
              dum= sqrt(taux**2+tauy**2)
              taux=taux/dum
              tauy=tauy/dum

              ui=zroesh(3,ip1,ish1)/zroesh(1,ip1,ish1)
              vi=zroesh(4,ip1,ish1)/zroesh(1,ip1,ish1)

              dum=ui*taux+vi*tauy

!             vshnor(1,ip1,ish1)=taux
!             vshnor(2,ip1,ish1)=tauy
              dum=taux*vshnor(1,ip1,ish1)+tauy*vshnor(2,ip1,ish1)
              if(dum.lt.0.)then
                   taux=-taux
                   tauy=-tauy
              endif
caldo
!             taux=1.0
!             tauy=0.0

              vshnor(1,ip1,ish1)=taux
              vshnor(2,ip1,ish1)=tauy

! vale
!            help=zroesh(3,ip1,ish1)**2+zroesh(4,ip1,ish1)**2
!            pd=gm1/ga*( zroesh(1,ip1,ish1)*zroesh(2,ip1,ish1)
!    &            -0.5d0*help)
!
!            help=zroeshu(3,ip1,ish1)**2+zroeshu(4,ip1,ish1)**2
!            pu=gm1/ga*(zroeshu(1,ip1,ish1)*zroeshu(2,ip1,ish1)
!    &            -0.5d0*help)
! vale

              if(dum.gt.0.)then
!vale           vshnor(1,ip1,ish1)=-taux
!vale           vshnor(2,ip1,ish1)=-tauy
              endif

              endif
          enddo
caldo

!         compute normal of the first internal point
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          ip2  = 2  ! tmp

          x1= xysh(1,ip1,ish1)
          y1= xysh(2,ip1,ish1)
          x2= xysh(1,ip2,ish1)
          y2= xysh(2,ip2,ish1)
          nx=vshnor(1,ip1,ish1)
          ny=vshnor(2,ip1,ish1)

          b=(-x1+x2-ny/(2.0*y1*nx)*(y1**2-y2**2))/ ! tmp
     &    (-y1+y2+1./(2.0*y1   )*(y1**2-y2**2))    ! tmp
          a=-(b*nx+ny)/(2*y1*nx)                   ! tmp
          c=x1-a*y1**2-b*y1                        ! tmp

!         ip2  = 2             ! tmp
!         x2= xysh(1,ip2,ish1) ! tmp
!         y2= xysh(2,ip2,ish1) ! tmp

          tauy=1.0
          taux=2*a*y2+b
          dum=sqrt(taux**2+tauy**2)

!         vshnor(1,ip2,ish1)=-tauy/dum
!         vshnor(2,ip2,ish1)=taux/dum

caldo

!      if the shock point is a connection or periodic connection

       elseif(typespecpoints(isppnts).eq.'C'
     +        .or.typespecpoints(isppnts).eq.'PC')then
          write(8,*)'correction for C and PC point '

!       determine shock and indices of extrema
!       shock 1
        ish1 = shinspps(1,1,isppnts)
        i    = shinspps(2,1,isppnts)-1
        ip1  = 1+i*(nshockpoints(ish1)-1)
!       shock 2
        ish2 = shinspps(1,2,isppnts)
        i    = shinspps(2,2,isppnts)-1
        ip2  = 1+i*(nshockpoints(ish2)-1)

!       correction of the normals
        nx1=vshnor(1,ip1,ish1)
        ny1=vshnor(2,ip1,ish1)
        nx2=vshnor(1,ip2,ish2)
        ny2=vshnor(2,ip2,ish2)

!       nx1=nx1+ nx2
!       ny1=ny1+ ny2

!       dum=sqrt(nx1*nx1+ny1*ny1)
!       nx1=nx1/dum
!       ny1=ny1/dum

!        vshnor(1,ip1,ish1)=nx1
!        vshnor(2,ip1,ish1)=ny1
!
!        vshnor(1,ip2,ish2)=nx1
!        vshnor(2,ip2,ish2)=ny1

! Attention: section of the code not general!

        if(ip1.eq.1)then
         nx1=nx2
         ny1=ny2
        else
         nx2=nx1
         ny2=ny1
        endif

        vshnor(1,ip1,ish1)=nx1
        vshnor(2,ip1,ish1)=ny1

        vshnor(1,ip2,ish2)=nx2
        vshnor(2,ip2,ish2)=ny2

! Attention: section of the code not general!
!
!      if the shock point in on the inlet section
!
       elseif(typespecpoints(isppnts).eq.'TP')then
          write(8,*)'correction for TP point '

!         determine shocks and indices of extrema
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

          nx2=vshnor(1,ip2,ish2)
          ny2=vshnor(2,ip2,ish2)

          nx4=vshnor(1,ip4,ish4)
          ny4=vshnor(2,ip4,ish4)

          dum=nx2*nx4+ny2*ny4
          if(dum.lt.0.0d0)then
            do i=1,nshockpoints(ish4)
!            vshnor(1,i,ish4)= -vshnor(1,i,ish4)
!            vshnor(2,i,ish4)= -vshnor(2,i,ish4)
            enddo
          endif

!      if Start Point

       elseif(typespecpoints(isppnts).eq.'SP')then
          write(8,*)'correction for sp point '

!         recover index of the end point (ip) and of the internal point (ip1)
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish1)-1)
          ip1=2+i*(nshockpoints(ish1)-3)

!       recover state upstream in the internal point
        um=zroeshu(3,ip1,ish1)/zroeshu(1,ip1,ish1)
        vm=zroeshu(4,ip1,ish1)/zroeshu(1,ip1,ish1)
        thetam=atan(vm/um)
        rom=zroeshu(1,ip1,ish1)*zroeshu(1,ip1,ish1)
        help=zroeshu(3,ip1,ish1)**2+zroeshu(4,ip1,ish1)**2
        pm =gm1/ga*( zroeshu(1,ip1,ish1)*zroeshu(2,ip1,ish1)-0.5d0*help)
!       compute upstream mach
        am=sqrt(ga*pm/rom)
        mm=sqrt(um**2+vm**2)/am
        if(mm.lt.1.0000)then
         write(*,*)'upstream mach number negative'
         write(*,*)'at shock point', ip1
         write(*,*)'shock n.',ish1
         stop
        endif
        alpham=asin(1./mm)
!       write(*,*)'monte'
!       write(*,*)'mach, thetaj, alphaj'
!       write(*,*)mm, thetam, alpham
!       write(*,*)cos(thetam-alpham),sin(thetam-alpham)

!       compute normals to a characteristic direction and
!       assign the correct slope to the point ip1
        nx1=-sin(thetam-alpham)
        ny1= cos(thetam-alpham)

        nx2=-sin(thetam+alpham)
        ny2= cos(thetam+alpham)

        dum1=nx1*vshnor(1,ip1,ish1)+ny1*vshnor(2,ip1,ish1)
        dum2=nx2*vshnor(1,ip1,ish1)+ny2*vshnor(2,ip1,ish1)

        vshnor(1,ip1,ish1)=nx2
        vshnor(2,ip1,ish1)=ny2

        if(abs(dum1).gt.abs(dum2))then
           vshnor(1,ip1,ish1)=nx1
           vshnor(2,ip1,ish1)=ny1
        endif

        dum=(um*vshnor(1,ip1,ish1)+vm*vshnor(2,ip1,ish1))/am
        if(dum.gt.0.)then
         vshnor(1,ip1,ish1)=-vshnor(1,ip1,ish1)
         vshnor(2,ip1,ish1)=-vshnor(2,ip1,ish1)
        endif

!       write(*,*)'normal mach:',
!    +             (um*vshnor(1,ip1,ish1)+vm*vshnor(2,ip1,ish1))/am

!       assignment of IP point slope
        vshnor(1,ip,ish1)= vshnor(1,ip1,ish1)
        vshnor(2,ip,ish1)= vshnor(2,ip1,ish1)

!       correction end point position
        nx=xysh(1,ip,ish1)-xysh(1,ip1,ish1)
        ny=xysh(2,ip,ish1)-xysh(2,ip1,ish1)
        dist=sqrt((xysh(1,ip1,ish1)-xysh(1,ip,ish1))**2+
     +            (xysh(2,ip1,ish1)-xysh(2,ip,ish1))**2)
        nx=nx/dist
        ny=ny/dist
        dum1= vshnor(2,ip,ish1)
        dum2=-vshnor(1,ip,ish1)
        if(dum1*nx+dum2*ny.lt.0.) then
           dum1=-dum1
           dum2=-dum2
        end if
        nx=dum1
        ny=dum2

        xysh(1,ip,ish1)=xysh(1,ip1,ish1)+nx*dist
        xysh(2,ip,ish1)=xysh(2,ip1,ish1)+ny*dist

       elseif(typespecpoints(isppnts).eq.'OPX')then

       elseif(typespecpoints(isppnts).eq.'OPY')then

       elseif(typespecpoints(isppnts).eq.'IPX')then

       elseif(typespecpoints(isppnts).eq.'IPY')then

       elseif(typespecpoints(isppnts).eq.'RRX')then

       elseif(typespecpoints(isppnts).eq.'RR')then

       elseif(typespecpoints(isppnts).eq.'QP')then

       elseif(typespecpoints(isppnts).eq.'EP')then

       elseif(typespecpoints(isppnts).eq.'TE')then

       else
          write(*,*) typespecpoints(isppnts)
          write(*,*)'condition not defined'
          write(8,*)'condition not defined'
          stop
       endif

      enddo

!     read normals of the Mach stem and of the reflected shock in the triple point

!     ish=3
!     open(66,file='norm_tp1',err=661)
!     read(66,*)vshnor(1,nshockpoints(ish),ish),
!    .          vshnor(2,nshockpoints(ish),ish)
!
!     ish=2
!     read(66,*)vshnor(1,nshockpoints(ish),ish),
!    .          vshnor(2,nshockpoints(ish),ish)
!661   close(66)

!     impose on the wall the direction of the shock normal parallel to that of the wall

!     vshnor(1,1,3)= 1.0d0
!     vshnor(2,1,3)= 0.0d0

!     write the tecplot file with computed normals

      write(8,*)'start writing file shocknor.dat'
      open(12,file='shocknor.dat')

      do ish=1,nshocks
      WRITE(12,*)'TITLE = Shock normals'
      WRITE(12,*)'VARIABLES = X Y Z(1) Z(2) NX NY'
      WRITE(12,*)'ZONE T="sampletext",F=FEPOINT,ET=TRIANGLE, N=',
     &nShockPoints(ISH),', E=',nShockPoints(ISH)-1
      do i = 1, nshockpoints(ish)
         write(12,*)(xysh(k,i,ish),k=1,ndim),1.d0,1.d0,
     .              (vshnor(k,i,ish),k=1,ndim)
      enddo
      do i = 1, nshockpoints(ish)-1
         write(12,*)i,i+1,i
      enddo
      end do
      close(12)
      write(8,*)'end writing file '
      close(8)

!     call exit(-1)

      return
      end

      integer function shp_dpndnc_old(x,y,ush,vsh,xi,yi,ui,vi,ai)
      double precision x,y,ush,vsh,xi,yi,ui,vi,ai
      double precision  taux,tauy,tau,utau

!     determine tangent versor
      taux=xi-x
      tauy=yi-y
      tau=sqrt(taux*taux+tauy*tauy)
      taux=taux/tau
      tauy=tauy/tau

!     determine the velocity component along the verson
      utau=ui*taux+vi*tauy

!     dependence evaluation
      shp_dpndnc=0
      if((utau-ai).lt.0.0d0)shp_dpndnc=1

      return
      end

      integer function shp_dpndnc_old2(x,y,ush,vsh,xi,yi,ui,vi,ai)
      double precision x,y,ush,vsh,xi,yi,ui,vi,ai,umod
      double precision xx,yy,xxi,yyi, dsh1, dsh2,dt
      double precision taux1,tauy1,tau1,taux2,tauy2,tau2,dum,dx1,dx2

!     calculate dt for the test
      rl=sqrt((x-xi)**2+(y-yi)**2)
      umod=sqrt(ui**2+vi**2)
      dt=0.5*rl/(ai+umod)

!     compute position of shock point after dt=1
      xx=x+ush*dt
      yy=y+vsh*dt

!     compute position of point i after dt=1
!     amplification factor of the shock velocity
!     to partially extend to the subsonic zone the computation
!     of the upwind formula for the calculation of shock normals
!
!     xxi=xi+3.0*ui*dt
!     yyi=yi+3.0*vi*dt

      xxi=xi+ui*dt
      yyi=yi+vi*dt

!     calculate original distance between shock point and point i
      dsh1=(x-xi)**2+(y-yi)**2
      dsh1=sqrt(dsh1)
      dx1=dsh1

!     determine versor
      taux1=xi-x
      tauy1=yi-y
      tau1=sqrt(taux1*taux1+tauy1*tauy1)
      taux1=taux1/tau1
      tauy1=tauy1/tau1

!     calculate the perturbation distance of point i and shock point after dt=1
      dsh2=(xx-xxi)**2+(yy-yyi)**2
      dx2=sqrt(dsh2)
!     dsh2=sqrt(dsh2)-ai*dt
      dsh2=sqrt(dsh2)-ai*dt
!     dsh1=dsh2-umod*dt

!     determine versor
      taux2=xxi-xx
      tauy2=yyi-yy
      tau2=sqrt(taux2*taux2+tauy2*tauy2)
      taux2=taux2/tau2
      tauy2=tauy2/tau2

      dum=taux1*taux2+tauy1*tauy2

!     distance evaluation
      shp_dpndnc=0
      if(dsh2.lt.dsh1-(umod*dt-abs(dx2-dx1)))shp_dpndnc=1
!     write(*,*)'dsh1,dsh2:',dsh1-(umod*dt-abs(dx2-dx1)),dsh2
!     if(dsh2.lt.0.d0)shp_dpndnc=1

      return
      end

      integer function shp_dpndnc(x,y,ush,vsh,xi,yi,ui,vi,ai)
      double precision x,y,ush,vsh,xi,yi,ui,vi,ai
      double precision xx,yy,uu,vv,aa,dum,dt1,dt2
      double precision  taux,tauy,tau,utau

      xx=x-xi
      yy=y-yi

      uu=ui
      vv=vi

      aa=ai*1.00

!     dum=2.0*xx*yy*uu*vv - xx*xx*vv*vv + xx*xx*aa*aa
!     dum=dum             - yy*yy*uu*uu + yy*yy*aa*aa

!     dum=(xx*uu+yy*vv)**2-(xx**2+yy**2)*(uu**2+vv**2-aa**2)
      dum= aa**2*(xx**2+yy**2)-(uu*yy-vv*xx)**2

      shp_dpndnc=0
      dt1=0.0
      dt2=0.0
      if(dum.gt.0.0)then
        dt1=((xx*uu+yy*vv)+sqrt(dum))/(uu**2+vv**2-aa**2)
        dt2=((xx*uu+yy*vv)-sqrt(dum))/(uu**2+vv**2-aa**2)
        if(dt1.gt.0.0.or.dt2.gt.0.0)shp_dpndnc=1
      endif
!     write(*,*) dum,dt1,dt2

      return
      end

      integer function dcp_dpndnc(x,y,ush,vsh,xi,yi,ui,vi,ai)
      double precision x,y,ush,vsh,xi,yi,ui,vi,ai
      double precision xx,yy,uu,vv,aa,dum,dt1,dt2
      double precision  taux,tauy,tau,utau

      xx=x-xi
      yy=y-yi

      uu=ui
      vv=vi

      aa=ai*1.00

      dum= (ui*xx+vi*yy)

      dcp_dpndnc=0
      if(dum.gt.0.0)dcp_dpndnc=1
!     write(*,*) dum,dt1,dt2

      return
      end
