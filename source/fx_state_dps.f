! Fix the state of the special points (ending points) of discontinuities

      subroutine fx_state_dps(
     +           xysh,
     +           xyshu,
     +           xyshd,
     +           zroeshu,
     +           zroeshd,
     +           zroeshuold,
     +           zroeshdold,
     +           vshnor,
     +           wsh,
     +           nshocks,
     +           nshockpoints,
     +           nshockedges,
     +           typesh,
     +           iter,
     +           nspecpoints,
     +           typespecpoints,
     +           shinspps,
     +           ispclr,
     +           ia,
     +           ja,
     +           iclr,
     +           nclr,
     +           corg)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer iter,nshocks,nshockpoints(nshmax),nshockedges(nshmax),
     +        nspecpoints,shinspps(2,5,*),ispclr(*)
      character*1 typesh(*)
      character*5 typespecpoints(*)

!     .. array arguments ..
      double precision
     +                 xysh       (ndim, npshmax, *),
     +                 xyshu      (ndim, npshmax, *),
     +                 xyshd      (ndim, npshmax, *),
     +                 zroeshu    (ndof, npshmax, *),
     +                 zroeshd    (ndof, npshmax, *),
     +                 zroeshuold (ndof, npshmax, *),
     +                 zroeshdold (ndof, npshmax, *),
     +                 vshnor     (ndim, npshmax, *),
     +                 wsh        (ndim, npshmax, *)
      integer nclr
      integer ia(*),ja(*),iclr(nclr)
      double precision corg(ndim,*)

!     .. local scalars ..
      double precision dx,dy,kine,ws,hh,wws,taux,tauy,dum,wrr,wss
      double precision wqpx,wqpy,w,cs,dum1,f1,f2,f3,dxr14,dyr14
      double precision help,x1(ndof),x2(ndof),taux12,tauy12,nx,ny
      double precision r2(npshmax,nshmax),xtpi(24),xtp(24),r23,r14
      double precision r14old,varz(ndof),avarz(ndof)
      double precision dumm2,dumx,dumy,rmach,gam,p0,rho0,astar
      double precision dumz1v,dumz2v,dumz3v,dumz4v
      double precision chxpu,chxmu,chypu,chymu,uu,vv,a,wshmod
      double precision chxpd,chxmd,chypd,chymd,rmachip1d,dump,dumm
      double precision pu,pd
      double precision dumx1,dumy1,dumx2,dumy2,xi,yi
      integer i,im,iv,k,totnshockpoints
      integer j,kp1,bbgn,bend,clr
      integer isppnts
      integer ip,ish,ip1,ip2,ip3,ip4,ip5,ish1,ish2,ish3,ish4,ish5
      logical flag1,ifail

!     assign constants
      gam=ga

!     open log file
      open(8,file='log/fx_state_sps.log')

      do isppnts=1,nspecpoints

!     if the shock point is on the inlet section
      if(typespecpoints(isppnts).eq.'IPX'.or.
     +   typespecpoints(isppnts).eq.'IPY')then

!       determine shock and index of the extreme
        ish = shinspps(1,1,isppnts)
        i   = shinspps(2,1,isppnts)-1
        ip  = 1+i*(nshockpoints(ish)-1)

!       restore previous upstream and downstream state
        do k=1, ndof
         zroeshu(k,ip,ish)=zroeshuold(k,ip,ish)
         zroeshd(k,ip,ish)=zroeshdold(k,ip,ish)
        end do

!       nullify shock velocity
        do k=1, ndim
         wsh(k,ip,ish)=0.0d+0
        enddo

!       if the shock point is on the exit section
        elseif(typespecpoints(isppnts).eq.'OPX'.or.
     +         typespecpoints(isppnts).eq.'OPY')then

!         determine shock and indices of the extremal points
          ish=shinspps(1,1,isppnts)
          i=  shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish)-1)

!         modify velocity for shock velocity
          ws=sqrt(wsh(1,ip,ish)**2+wsh(2,ip,ish)**2)
          dx=wsh(1,ip,ish)/ws
          dy=wsh(2,ip,ish)/ws

          if(typespecpoints(isppnts).eq.'OPX')then
            wsh(1,ip,ish)=ws/dx
            wsh(2,ip,ish)=0.0d+0
          else
            wsh(1,ip,ish)=0.0d+0
            wsh(2,ip,ish)=ws/dy
          endif

!       if wall shock point
        elseif(typespecpoints(isppnts).eq.'WPNRX'.or.
     +         typespecpoints(isppnts).eq.'WPNRY')then

!         determine shock and indices of the extremal points
          ish=shinspps(1,1,isppnts)
          i=  shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish)-1)

          ws=sqrt(wsh(1,ip,ish)**2+wsh(2,ip,ish)**2)
          dx=wsh(1,ip,ish)/(ws)
          dy=wsh(2,ip,ish)/(ws)

          if(typespecpoints(isppnts).eq.'WPNRX')then
            wsh(1,ip,ish)=ws/dx
            wsh(2,ip,ish)=0.0d+0
          else
            wsh(1,ip,ish)=0.0d+0
            wsh(2,ip,ish)=ws/dy
          endif

       elseif(typespecpoints(isppnts).eq.'FWP') then

!         determine shock and indices of the extremal points
          ish=shinspps(1,1,isppnts)
          i=  shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish)-1)

          dx=vshnor(1,ip,ish)
          dy=vshnor(2,ip,ish)

          dum=zroeshu(3,ip,ish)*dx+zroeshu(4,ip,ish)*dy

          dumx=zroeshu(3,ip,ish)-dum*dx
          dumy=zroeshu(4,ip,ish)-dum*dy

          zroeshu(3,ip,ish)=zroeshu(3,ip,ish)-dumx
          zroeshu(4,ip,ish)=zroeshu(4,ip,ish)-dumy

          dum=zroeshd(3,ip,ish)*dx+zroeshd(4,ip,ish)*dy

          dumx=zroeshd(3,ip,ish)-dum*dx
          dumy=zroeshd(4,ip,ish)-dum*dy
          zroeshd(3,ip,ish)=zroeshd(3,ip,ish)-dumx
          zroeshd(4,ip,ish)=zroeshd(4,ip,ish)-dumy

caldo
          if(vshnor(2,ip,ish).gt.0.3)then
             wsh(1,ip,ish)=wsh(1,ip,ish)/dx
             wsh(2,ip,ish)=wsh(1,ip,ish)/dx
          endif
caldo

!       if triple point
        elseif(typespecpoints(isppnts).eq.'TP')then

! determine shocks and indices of extremal points

! incident shock
          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)
! reflected shock
          ish2 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip2  = 1+i*(nshockpoints(ish2)-1)
! mach stem
          ish3 = shinspps(1,3,isppnts)
          i    = shinspps(2,3,isppnts)-1
          ip3  = 1+i*(nshockpoints(ish3)-1)
! contact discontinuity
          ish4 = shinspps(1,4,isppnts)
          i    = shinspps(2,4,isppnts)-1
          ip4  = 1+i*(nshockpoints(ish4)-1)

caldo
! determine family of shock  1
           f1=vshnor(1, ip1, ish1)*zroeshu(4,ip1,ish1)-
     .        vshnor(2, ip1, ish1)*zroeshu(3,ip1,ish1)

           f1=-sign(1.d0,f1)

! determine family of shock 2
           f2=vshnor(1, ip2, ish2)*zroeshu(4,ip2,ish2)-
     .        vshnor(2, ip2, ish2)*zroeshu(3,ip2,ish2)

           f2=-sign(1.d0,f2)

! retrieve states

! state 1

!          if the incident and reflected shocks belong to the opposite families, then ...
!          state 1 = upstream state of incident shock
           if(f1*f2.lt.0.)then
           xtpi(4) = zroeshu(4,ip1,ish1)/zroeshu(1,ip1,ish1) ! y-component
           xtpi(3) = zroeshu(3,ip1,ish1)/zroeshu(1,ip1,ish1) ! x-component
           xtpi(1) = zroeshu(1,ip1,ish1)*zroeshu(1,ip1,ish1) ! density
           help    = zroeshu(3,ip1,ish1)**2+zroeshu(4,ip1,ish1)**2
           xtpi(2) = gm1/ga*( zroeshu(1,ip1,ish1)*zroeshu(2,ip1,ish1)
     .              -0.5d0*help) ! pressure
           else
!          if the incident and reflected shocks belong to the same family, then ...
!          state 1 = downstream state of incident shock
           xtpi(4) = zroeshd(4,ip1,ish1)/zroeshd(1,ip1,ish1) ! y-component
           xtpi(3) = zroeshd(3,ip1,ish1)/zroeshd(1,ip1,ish1) ! x-component
           xtpi(1) = zroeshd(1,ip1,ish1)*zroeshd(1,ip1,ish1) ! density
           help    = zroeshd(3,ip1,ish1)**2+zroeshd(4,ip1,ish1)**2
           xtpi(2) = gm1/ga*( zroeshd(1,ip1,ish1)*zroeshd(2,ip1,ish1)
     .              -0.5d0*help) ! pressure
           endif
caldo

! state 2
           xtpi(8) = ZROESHu(4,IP2,ish2)/ZROESHu(1,IP2,ish2) ! y-component
           xtpi(7) = ZROESHu(3,IP2,ish2)/ZROESHu(1,IP2,ish2) ! x-component
           xtpi(5) = ZROESHu(1,IP2,ish2)*ZROESHu(1,IP2,ish2) ! density
           HELP    = ZROESHu(3,IP2,ish2)**2+ZROESHu(4,IP2,ish2)**2
           xtpi(6) = GM1/GA*( ZROESHu(1,IP2,ish2)*ZROESHu(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure
! state 3
           xtpi(12)= ZROESHd(4,IP2,ish2)/ZROESHd(1,IP2,ish2) ! y-component
           xtpi(11)= ZROESHd(3,IP2,ish2)/ZROESHd(1,IP2,ish2) ! x-component
           xtpi(9) = ZROESHd(1,IP2,ish2)*ZROESHd(1,IP2,ish2) ! density
           HELP    = ZROESHd(3,IP2,ish2)**2+ZROESHd(4,IP2,ish2)**2
           xtpi(10)= GM1/GA*( ZROESHd(1,IP2,ish2)*ZROESHd(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

           DX=VSHNOR(1, IP2, ISH2)
           DY=VSHNOR(2, IP2, ISH2)

      R23=sqrt(GA*xtpi(10)/xtpi( 9))+
     &        GM1*0.5d0*(xtpi(11)*DX+xtpi(12)*DY)

! state 4
           xtpi(16)= ZROESHd(4,IP3,ish3)/ZROESHd(1,IP3,ish3) ! y-component
           xtpi(15)= ZROESHd(3,IP3,ish3)/ZROESHd(1,IP3,ish3) ! x-component
           xtpi(13)= ZROESHd(1,IP3,ish3)*ZROESHd(1,IP3,ish3) ! density
           HELP    = ZROESHd(3,IP3,ish3)**2+ZROESHd(4,IP3,ish3)**2
           xtpi(14)= GM1/GA*( ZROESHd(1,IP3,ish3)*ZROESHd(2,IP3,ish3)
     .              -0.5d0*HELP) ! pressure

           DX=VSHNOR(1, IP3, ISH3)
           DY=VSHNOR(2, IP3, ISH3)

           DXR14=DX
           DYR14=DY

      R14=sqrt(GA*xtpi(14)/xtpi(13))+
     &        GM1*0.5d0*(xtpi(15)*DXR14+xtpi(16)*DYR14)
!          write(*,*)'DXR14,DYR14:',DXR14,DYR14
!          write(*,*)
!          pause

caldo

! state 3
          xtpi(12)=ZROESHdOLD(4,IP2,ish2)/ZROESHdOLD(1,IP2,ish2) ! y-component
          xtpi(11)=ZROESHdOLD(3,IP2,ish2)/ZROESHdOLD(1,IP2,ish2) ! x-component
          xtpi(9) =ZROESHdOLD(1,IP2,ish2)*ZROESHdOLD(1,IP2,ish2) ! density
          HELP    =ZROESHdOLD(3,IP2,ish2)**2+ZROESHdOLD(4,IP2,ish2)**2
          xtpi(10)=GM1/GA*(ZROESHdOLD(1,IP2,ish2)*ZROESHdOLD(2,IP2,ish2)
     .            -0.5d0*HELP)           ! pressure
! state 4
          xtpi(16)=ZROESHdOLD(4,IP3,ish3)/ZROESHdOLD(1,IP3,ish3) ! y-component
          xtpi(15)=ZROESHdOLD(3,IP3,ish3)/ZROESHdOLD(1,IP3,ish3) ! x-component
          xtpi(13)=ZROESHdOLD(1,IP3,ish3)*ZROESHdOLD(1,IP3,ish3) ! density
          HELP    =ZROESHdOLD(3,IP3,ish3)**2+ZROESHdOLD(4,IP3,ish3)**2
          xtpi(14)=GM1/GA*(ZROESHdOLD(1,IP3,ish3)*ZROESHdOLD(2,IP3,ish3)
     .            -0.5d0*HELP)           ! pressure

       R14OLD=sqrt(GA*xtpi(14)/xtpi(13))+
     &        GM1*0.5d0*(xtpi(15)*DXR14+xtpi(16)*DYR14)

caldo

! slopes of the shocks

!         shock 1 (incident shock)
          dx=vshnor(1, ip1, ish1)
          dy=vshnor(2, ip1, ish1)
          xtpi(17)=atan2( -dx,dy)

!         shock 2 (reflected shock)
          dx=vshnor(1, ip2, ish2)
          dy=vshnor(2, ip2, ish2)
          xtpi(18)=atan2( -dx,dy)

!         shock 3 (mach stem)
          dx=vshnor(1, ip3, ish3)
          dy=vshnor(2, ip3, ish3)
          xtpi(19)=atan2( -DX,DY)

!         set velocity of triple point
!         xtpi(20)=0.
          xtpi(20)=-0.01

!         recover normal velocity of the incident shock
          wws=sqrt(wsh(1,ip1,ish1)**2+wsh(2,ip1,ish1)**2)

          do i=1,20
            xtp(i)=xtpi(i)
          enddo

100       FLAG1=.false.
!         if(ITER.GT.0)FLAG1=.true.
!         IF(ITER.GT.01)THEN
          IFAIL=.false.
          call co_utp(xtpi,R14,DXR14,DYR14,R23,WWS,xtp,FLAG1,IFAIL)
!         xtp(16)  ! y-component
!         xtp(15)  ! x-component
!         xtp(13)  ! density
!         xtp(14)) ! pressure

!        taux12=cos(xtp(17))
!        tauy12=sin(xtp(17))
!        dumx=xtp(15)-taux12*xtp(20)
!        dumy=xtp(16)-tauy12*xtp(20)
!        dumm2=DXR14*xtp(15)+DYR14*xtp(16)
!        dumm2=dumm2-xtp(20)*taux12*DXR14
!        dumm2=dumm2-xtp(20)*tauy12*DYR14
!        dumm2=dumm2*dumm2/(GA*xtp(14)/xtp(13))
!        write(*,*)'Mach a valle:',sqrt(dumm2)
!        if(dumm2.GT.1.0)then
!            write(*,*)'Ricalculate'
!            R14=R14OLD
!            GOTO 100
!        endIF
!
!        write(*,*)
!        write(*,*)'Downstream Mach:',sqrt(dumm2)
!        pause
!        continue

          IF(IFAIL)THEN
           IFAIL=.false.
           xtpi(20)=-xtpi(20)
           call co_utp(xtpi,R14,DXR14,DYR14,R23,WWS,xtp,FLAG1,IFAIL)

           IF(IFAIL)THEN
             stop 'tp does not converge'
           endif
          endif

!         ENDIF

!         do K=1,20
!          write(*,*)k,xtpi(k),xtp(k)
!         enddo

!  assign values in zone 1

         z1v = sqrt(xtp(1))
         kine = 0.5d0*(xtp( 3)*xtp( 3)+xtp( 4)*xtp( 4))
         hh = GA/GM1*xtp( 2)/xtp(1)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 3)
         z4v = z1v*xtp( 4)

! if the incident and reflected shocks belong to opposite families, then  ...
! .. to upstream incident shock 
!         if(f1*f2.lt.0.)then
!
!          ZROESHu(1,IP1,ISH1)=z1v
!          ZROESHu(2,IP1,ISH1)=z2v
!          ZROESHu(3,IP1,ISH1)=z3v
!          ZROESHu(4,IP1,ISH1)=z4v
!         else
! if the incident and reflected shocks belong to the same family, then  ...
! .. to downstream incident shock
!          ZROESHd(1,IP1,ISH1)=z1v
!          ZROESHd(2,IP1,ISH1)=z2v
!          ZROESHd(3,IP1,ISH1)=z3v
!          ZROESHd(4,IP1,ISH1)=z4v
!         endif

! .. to upstream mach stem
          ZROESHu(1,IP3,ISH3)=z1v
          ZROESHu(2,IP3,ISH3)=z2v
          ZROESHu(3,IP3,ISH3)=z3v
          ZROESHu(4,IP3,ISH3)=z4v

! assign values in zone 2 ..
         z1v = sqrt(xtp(5))
         kine = 0.5d0*(xtp( 7)*xtp( 7)+xtp( 8)*xtp( 8))
         hh = GA/GM1*xtp( 6)/xtp(5)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 7)
         z4v = z1v*xtp( 8)

!         if(f1*f2.lt.0.)then
! if the incident and reflected shocks belong to opposite family, then  ...
! .. to downstream incident shock
!          ZROESHd(1,IP1,ISH1)=z1v
!          ZROESHd(2,IP1,ISH1)=z2v
!          ZROESHd(3,IP1,ISH1)=z3v
!          ZROESHd(4,IP1,ISH1)=z4v
!         else
! if the incident and reflected shocks belong to same family, then  ...
! .. to upstream incident shock
!          ZROESHu(1,IP1,ISH1)=z1v
!          ZROESHu(2,IP1,ISH1)=z2v
!          ZROESHu(3,IP1,ISH1)=z3v
!          ZROESHu(4,IP1,ISH1)=z4v
!         endif

! .. to upstream reflected shock
         ZROESHu(1,IP2,ISH2)=z1v
         ZROESHu(2,IP2,ISH2)=z2v
         ZROESHu(3,IP2,ISH2)=z3v
         ZROESHu(4,IP2,ISH2)=z4v

! assign values in zone 3
         z1v = sqrt(xtp(9))
         kine = 0.5d0*(xtp(11)*xtp(11)+xtp(12)*xtp(12))
         hh = GA/GM1*xtp(10)/xtp(9)+kine
         z2v = z1v*hh
         z3v = z1v*xtp(11)
         z4v = z1v*xtp(12)

caldo
!        dumz1v=z1v-ZROESHd(1,IP2,ISH2)
!        dumz2v=z2v-ZROESHd(2,IP2,ISH2)
!        dumz3v=z3v-ZROESHd(3,IP2,ISH2)
!        dumz4v=z4v-ZROESHd(4,IP2,ISH2)
!
!        write(*,*)'dumz1v:',dumz1v
!        write(*,*)'dumz2v:',dumz2v
!        write(*,*)'dumz3v:',dumz3v
!        write(*,*)'dumz4v:',dumz4v
!        write(*,*)
!
!        ZROESHd(1,IP2+1,ISH2)=ZROESHd(1,IP2+1,ISH2)+dumz1v*0.5
!        ZROESHd(2,IP2+1,ISH2)=ZROESHd(2,IP2+1,ISH2)+dumz2v*0.5
!        ZROESHd(3,IP2+1,ISH2)=ZROESHd(3,IP2+1,ISH2)+dumz3v*0.5
!        ZROESHd(4,IP2+1,ISH2)=ZROESHd(4,IP2+1,ISH2)+dumz4v*0.5
!
!        ZROESHu(1,IP4+1,ISH4)=ZROESHu(1,IP4+1,ISH4)+dumz1v*0.5
!        ZROESHu(2,IP4+1,ISH4)=ZROESHu(2,IP4+1,ISH4)+dumz2v*0.5
!        ZROESHu(3,IP4+1,ISH4)=ZROESHu(3,IP4+1,ISH4)+dumz3v*0.5
!        ZROESHu(4,IP4+1,ISH4)=ZROESHu(4,IP4+1,ISH4)+dumz4v*0.5

! .. to downstream reflected shock
         ZROESHd(1,IP2,ISH2)=z1v
         ZROESHd(2,IP2,ISH2)=z2v
         ZROESHd(3,IP2,ISH2)=z3v
         ZROESHd(4,IP2,ISH2)=z4v

! .. to upstream of the contact discontinuity (attention to the orientation of the contact discontinuity)
         ZROESHu(1,IP4,ISH4)=z1v
         ZROESHu(2,IP4,ISH4)=z2v
         ZROESHu(3,IP4,ISH4)=z3v
         ZROESHu(4,IP4,ISH4)=z4v

! assign values in zone 4
         z1v = sqrt(xtp(13))
         kine = 0.5d0*(xtp(15)*xtp(15)+xtp(16)*xtp(16))
         hh = GA/GM1*xtp(14)/xtp(13)+kine
         z2v = z1v*hh
         z3v = z1v*xtp(15)
         z4v = z1v*xtp(16)

! .. to downstream of the mach stem
         ZROESHd(1,IP3,ISH3)=z1v
         ZROESHd(2,IP3,ISH3)=z2v
         ZROESHd(3,IP3,ISH3)=z3v
         ZROESHd(4,IP3,ISH3)=z4v

! .. to downstream of the contact discontinuity (attention to the orientation of the contact discontinuity)
         ZROESHd(1,IP4,ISH4)=z1v
         ZROESHd(2,IP4,ISH4)=z2v
         ZROESHd(3,IP4,ISH4)=z3v
         ZROESHd(4,IP4,ISH4)=z4v
caldo

! sum the contribution along incident shock to the velocity of triple point
         taux12=cos(xtp(17))
         tauy12=sin(xtp(17))

         WSH(1,IP1,ISH1)= WSH(1,IP1,ISH1)+ taux12*xtp(20)*1.0
         WSH(2,IP1,ISH1)= WSH(2,IP1,ISH1)+ tauy12*xtp(20)*1.0

!        write(*,*)'taux12,tauy12:',taux12,tauy12
!        write(*,*)'xtp(20)',xtp(20)
!        write(*,*)'WSH:',WSH(1,IP1,ISH1),WSH(2,IP1,ISH1)
!        write(*,*)

! propagate to the other shocks the new velocity
         WSH(1,IP2,ISH2)= WSH(1,IP1,ISH1)
         WSH(2,IP2,ISH2)= WSH(2,IP1,ISH1)

         WSH(1,IP3,ISH3)= WSH(1,IP1,ISH1)
         WSH(2,IP3,ISH3)= WSH(2,IP1,ISH1)

         WSH(1,IP4,ISH4)= WSH(1,IP1,ISH1)
         WSH(2,IP4,ISH4)= WSH(2,IP1,ISH1)


! if quadruple point
        ELSEIF(TypeSpecPoints(ISPPNTS).EQ.'QP')THEN

! determine shocks and indices of extremal points

! incident shock 1
          ISH1 = SHinSPPs(1,1,ISPPNTS)
          I    = SHinSPPs(2,1,ISPPNTS)-1
          IP1  = 1+I*(nShockPoints(ISH1)-1)
! reflected shock 1
          ISH2 = SHinSPPs(1,2,ISPPNTS)
          I    = SHinSPPs(2,2,ISPPNTS)-1
          IP2  = 1+I*(nShockPoints(ISH2)-1)
! incident shock 2
          ISH3 = SHinSPPs(1,3,ISPPNTS)
          I    = SHinSPPs(2,3,ISPPNTS)-1
          IP3  = 1+I*(nShockPoints(ISH3)-1)
! reflected shock 2
          ISH4 = SHinSPPs(1,4,ISPPNTS)
          I    = SHinSPPs(2,4,ISPPNTS)-1
          IP4  = 1+I*(nShockPoints(ISH4)-1)
! conctact discontinuity
          ISH5 = SHinSPPs(1,5,ISPPNTS)
          I    = SHinSPPs(2,5,ISPPNTS)-1
          IP5  = 1+I*(nShockPoints(ISH5)-1)

! determine family of shock 1
           f1=VSHNOR(1, IP1, ISH1)*ZROESHu(4,IP1,ish1)-
     .        VSHNOR(2, IP1, ISH1)*ZROESHu(3,IP1,ish1)

           f1=-SIGN(1.d0,f1)

! determine family of shock 3
           f3=VSHNOR(1, IP3, ISH3)*ZROESHu(4,IP3,ish3)-
     .        VSHNOR(2, IP3, ISH3)*ZROESHu(3,IP3,ish3)

           f3=-SIGN(1.d0,f3)

! recover states

! state 1
           if(f1.gt.0.d0)then
           xtpi(4) = ZROESHu(4,IP1,ish1)/ZROESHu(1,IP1,ish1) ! y-component
           xtpi(3) = ZROESHu(3,IP1,ish1)/ZROESHu(1,IP1,ish1) ! x-component
           xtpi(1) = ZROESHu(1,IP1,ish1)*ZROESHu(1,IP1,ish1) ! density
           HELP    = ZROESHu(3,IP1,ish1)**2+ZROESHu(4,IP1,ish1)**2
           xtpi(2) = GM1/GA*( ZROESHu(1,IP1,ish1)*ZROESHu(2,IP1,ish1)
     .              -0.5d0*HELP) ! pressure
           else
           xtpi(4) = ZROESHd(4,IP1,ish1)/ZROESHd(1,IP1,ish1) ! y-component
           xtpi(3) = ZROESHd(3,IP1,ish1)/ZROESHd(1,IP1,ish1) ! x-component
           xtpi(1) = ZROESHd(1,IP1,ish1)*ZROESHd(1,IP1,ish1) ! density
           HELP    = ZROESHd(3,IP1,ish1)**2+ZROESHd(4,IP1,ish1)**2
           xtpi(2) = GM1/GA*( ZROESHd(1,IP1,ish1)*ZROESHd(2,IP1,ish1)
     .              -0.5d0*HELP) ! pressure
           endif

! state 2
           xtpi(8) = ZROESHu(4,IP2,ish2)/ZROESHu(1,IP2,ish2) ! y-component
           xtpi(7) = ZROESHu(3,IP2,ish2)/ZROESHu(1,IP2,ish2) ! x-component
           xtpi(5) = ZROESHu(1,IP2,ish2)*ZROESHu(1,IP2,ish2) ! density
           HELP    = ZROESHu(3,IP2,ish2)**2+ZROESHu(4,IP2,ish2)**2
           xtpi(6) = GM1/GA*( ZROESHu(1,IP2,ish2)*ZROESHu(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

! state 3
           xtpi(12)= ZROESHu(4,IP4,ish4)/ZROESHu(1,IP4,ish4) ! y-component
           xtpi(11)= ZROESHu(3,IP4,ish4)/ZROESHu(1,IP4,ish4) ! x-component
           xtpi(9) = ZROESHu(1,IP4,ish4)*ZROESHu(1,IP4,ish4) ! density
           HELP    = ZROESHu(3,IP4,ish4)**2+ZROESHu(4,IP4,ish4)**2
           xtpi(10)= GM1/GA*( ZROESHu(1,IP4,ish4)*ZROESHu(2,IP4,ish4)
     .              -0.5d0*HELP) ! pressure

! state 4
           xtpi(16)= ZROESHd(4,IP2,ish2)/ZROESHd(1,IP2,ish2) ! y-component
           xtpi(15)= ZROESHd(3,IP2,ish2)/ZROESHd(1,IP2,ish2) ! x-component
           xtpi(13)= ZROESHd(1,IP2,ish2)*ZROESHd(1,IP2,ish2) ! density
           HELP    = ZROESHd(3,IP2,ish2)**2+ZROESHd(4,IP2,ish2)**2
           xtpi(14)= GM1/GA*( ZROESHd(1,IP2,ish2)*ZROESHd(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

! state 5
           xtpi(20)= ZROESHd(4,IP4,ish4)/ZROESHd(1,IP4,ish4) ! y-component
           xtpi(19)= ZROESHd(3,IP4,ish4)/ZROESHd(1,IP4,ish4) ! x-component
           xtpi(17)= ZROESHd(1,IP4,ish4)*ZROESHd(1,IP4,ish4) ! density
           HELP    = ZROESHd(3,IP4,ish4)**2+ZROESHd(4,IP4,ish4)**2
           xtpi(18)= GM1/GA*( ZROESHd(1,IP4,ish4)*ZROESHd(2,IP4,ish4)
     .              -0.5d0*HELP) ! pressure

! slopes of the shocks
 
! shock 1 (incident shock 1)
          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)
          xtpi(21)=atan2( -DX,DY)

! shock 2 (reflected shock 1)
          DX=VSHNOR(1, IP2, ISH2)
          DY=VSHNOR(2, IP2, ISH2)
          xtpi(22)=atan2( -DX,DY)

! shock 3 (incident shock 2)
          DX=VSHNOR(1, IP3, ISH3)
          DY=VSHNOR(2, IP3, ISH3)
          xtpi(23)=atan2( -DX,DY)

! shock 4 (reflected shock 2)
          DX=VSHNOR(1, IP4, ISH4)
          DY=VSHNOR(2, IP4, ISH4)
          xtpi(24)=atan2( -DX,DY)

! recover velocity of quadruple point
!         WQPX=WSH(1,IP1,ISH1)+WSH(1,IP3,ISH3)
!         WQPY=WSH(2,IP1,ISH1)+WSH(2,IP3,ISH3)

          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)
          taux=-DY
          tauy=+DX
          DX=VSHNOR(1, IP3, ISH3)
          DY=VSHNOR(2, IP3, ISH3)

          cs=DX*taux+DY*tauy
          dum1=DX*WSH(1,IP3,ISH3)+DY*WSH(2,IP3,ISH3)

          w=sign(1.d0,dum1)*sqrt(WSH(1,IP3,ISH3)**2+WSH(2,IP3,ISH3)**2)

          WQPX=w/cs*taux
          WQPY=w/cs*tauy

          DX=VSHNOR(1, IP3, ISH3)
          DY=VSHNOR(2, IP3, ISH3)
          taux=-DY
          tauy=+DX
          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)

          cs=DX*taux+DY*tauy

          dum1=DX*WSH(1,IP1,ISH1)+DY*WSH(2,IP1,ISH1)

          w=sign(1.d0,dum1)*sqrt(WSH(1,IP1,ISH1)**2+WSH(2,IP1,ISH1)**2)

          WQPX=WQPX+w/cs*taux
          WQPY=WQPY+w/cs*tauy

          DO I=1,24
            xtp(i)=xtpi(i)
!           write(*,*)i,xtp(i)
          ENDDO

!         if(iter.gt.01)then
!           calculate non-stationary quadruple point
            call co_uqp(xtpi,wqpx,wqpy,xtp)
!         endif

! assign values in zone 1
         z1v = sqrt(xtp(1))
         kine = 0.5d0*(xtp( 3)*xtp( 3)+xtp( 4)*xtp( 4))
         hh = GA/GM1*xtp( 2)/xtp(1)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 3)
         z4v = z1v*xtp( 4)

! .. to upstream of the incident shock 1
         if(f1.gt.0.d0)then
          ZROESHu(1,IP1,ISH1)=z1v
          ZROESHu(2,IP1,ISH1)=z2v
          ZROESHu(3,IP1,ISH1)=z3v
          ZROESHu(4,IP1,ISH1)=z4v
         else
          ZROESHd(1,IP1,ISH1)=z1v
          ZROESHd(2,IP1,ISH1)=z2v
          ZROESHd(3,IP1,ISH1)=z3v
          ZROESHd(4,IP1,ISH1)=z4v
         endif

! .. to upstream of the incident shock 2
        if(f3.lt.0.d0)then
         ZROESHu(1,IP3,ISH3)=z1v
         ZROESHu(2,IP3,ISH3)=z2v
         ZROESHu(3,IP3,ISH3)=z3v
         ZROESHu(4,IP3,ISH3)=z4v
        else

! .. to downstream of the incident shock 2
         ZROESHd(1,IP3,ISH3)=z1v
         ZROESHd(2,IP3,ISH3)=z2v
         ZROESHd(3,IP3,ISH3)=z3v
         ZROESHd(4,IP3,ISH3)=z4v
        endif

! assign values in zone 2
         z1v = sqrt(xtp(5))
         kine = 0.5d0*(xtp( 7)*xtp( 7)+xtp( 8)*xtp( 8))
         hh = GA/GM1*xtp( 6)/xtp(5)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 7)
         z4v = z1v*xtp( 8)

! .. to downstream of the incident shock 1
         if(f1.gt.0.d0)then
          ZROESHd(1,IP1,ISH1)=z1v
          ZROESHd(2,IP1,ISH1)=z2v
          ZROESHd(3,IP1,ISH1)=z3v
          ZROESHd(4,IP1,ISH1)=z4v
         else
          ZROESHu(1,IP1,ISH1)=z1v
          ZROESHu(2,IP1,ISH1)=z2v
          ZROESHu(3,IP1,ISH1)=z3v
          ZROESHu(4,IP1,ISH1)=z4v
         endif

! .. to upstream of the reflected shock 1
         ZROESHu(1,IP2,ISH2)=z1v
         ZROESHu(2,IP2,ISH2)=z2v
         ZROESHu(3,IP2,ISH2)=z3v
         ZROESHu(4,IP2,ISH2)=z4v

! assign values in zone 3
         z1v = sqrt(xtp(9))
         kine = 0.5d0*(xtp( 11)*xtp( 11)+xtp( 12)*xtp( 12))
         hh = GA/GM1*xtp( 10)/xtp(9)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 11)
         z4v = z1v*xtp( 12)

! .. to downstream of the incident shock 2
         if(f3.lt.0.d0)then
          ZROESHd(1,IP3,ISH3)=z1v
          ZROESHd(2,IP3,ISH3)=z2v
          ZROESHd(3,IP3,ISH3)=z3v
          ZROESHd(4,IP3,ISH3)=z4v
         else

! .. to upstream of the incident shock 2
          ZROESHu(1,IP3,ISH3)=z1v
          ZROESHu(2,IP3,ISH3)=z2v
          ZROESHu(3,IP3,ISH3)=z3v
          ZROESHu(4,IP3,ISH3)=z4v
         endif

! .. to upstream of the reflected shock 2
         ZROESHu(1,IP4,ISH4)=z1v
         ZROESHu(2,IP4,ISH4)=z2v
         ZROESHu(3,IP4,ISH4)=z3v
         ZROESHu(4,IP4,ISH4)=z4v

! assign values in zone 4
         z1v = sqrt(xtp(13))
         kine = 0.5d0*(xtp( 15)*xtp( 15)+xtp( 16)*xtp( 16))
         hh = GA/GM1*xtp( 14)/xtp(13)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 15)
         z4v = z1v*xtp( 16)

! .. to downstream of the reflected shock 1
         ZROESHd(1,IP2,ISH2)=z1v
         ZROESHd(2,IP2,ISH2)=z2v
         ZROESHd(3,IP2,ISH2)=z3v
         ZROESHd(4,IP2,ISH2)=z4v

! .. to downstream of the contact discontinuity
         ZROESHu(1,IP5,ISH5)=z1v
         ZROESHu(2,IP5,ISH5)=z2v
         ZROESHu(3,IP5,ISH5)=z3v
         ZROESHu(4,IP5,ISH5)=z4v

! assign values in zone 5
         z1v = sqrt(xtp(17))
         kine = 0.5d0*(xtp( 19)*xtp( 19)+xtp( 20)*xtp( 20))
         hh = GA/GM1*xtp( 18)/xtp(17)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 19)
         z4v = z1v*xtp( 20)

! .. to downstream of the reflected shock 
         ZROESHd(1,IP4,ISH4)=z1v
         ZROESHd(2,IP4,ISH4)=z2v
         ZROESHd(3,IP4,ISH4)=z3v
         ZROESHd(4,IP4,ISH4)=z4v

! .. to downstream of the contact discontinuity
         ZROESHd(1,IP5,ISH5)=z1v
         ZROESHd(2,IP5,ISH5)=z2v
         ZROESHd(3,IP5,ISH5)=z3v
         ZROESHd(4,IP5,ISH5)=z4v

! assign velocity of the triple point to ... 

! .... incident shock 1
         WSH(1,IP1,ISH1)= WQPX
         WSH(2,IP1,ISH1)= WQPY

! .... reflected shock 1
         WSH(1,IP2,ISH2)= WSH(1,IP1,ISH1)
         WSH(2,IP2,ISH2)= WSH(2,IP1,ISH1)

! .... incident shock 2
         WSH(1,IP3,ISH3)= WSH(1,IP1,ISH1)
         WSH(2,IP3,ISH3)= WSH(2,IP1,ISH1)

! .... reflected shock 1
         WSH(1,IP4,ISH4)= WSH(1,IP1,ISH1)
         WSH(2,IP4,ISH4)= WSH(2,IP1,ISH1)

! .... contact discontinuity
         WSH(1,IP5,ISH5)= WSH(1,IP1,ISH1)
         WSH(2,IP5,ISH5)= WSH(2,IP1,ISH1)

! if the point is trailing edge
        ELSEIF(TypeSpecPoints(ISPPNTS).EQ.'TE')THEN

! determine shocks and indices of extremal points

! incident shock 1
          ISH1 = 0
          I    = 0
          IP1  = 0
! reflected shock 1
          ISH2 = SHinSPPs(1,1,ISPPNTS)
          I    = SHinSPPs(2,1,ISPPNTS)-1
          IP2  = 1+I*(nShockPoints(ISH2)-1)
! incident shock 2
          ISH3 = 0
          I    = 0
          IP3  = 0
! incident shock 2
          ISH4 = SHinSPPs(1,3,ISPPNTS)
          I    = SHinSPPs(2,3,ISPPNTS)-1
          IP4  = 1+I*(nShockPoints(ISH4)-1)
! contact discontinuity
          ISH5 = SHinSPPs(1,2,ISPPNTS)
          I    = SHinSPPs(2,2,ISPPNTS)-1
          IP5  = 1+I*(nShockPoints(ISH5)-1)

! recover states

! state 1
           xtpi(4) = 1.0
           xtpi(3) = 1.0
           xtpi(1) = 1.0
           xtpi(2) = 1.0

! state 2
           xtpi(8) = ZROESHu(4,IP2,ish2)/ZROESHu(1,IP2,ish2) ! y-component
           xtpi(7) = ZROESHu(3,IP2,ish2)/ZROESHu(1,IP2,ish2) ! x-component
           xtpi(5) = ZROESHu(1,IP2,ish2)*ZROESHu(1,IP2,ish2) ! density
           HELP    = ZROESHu(3,IP2,ish2)**2+ZROESHu(4,IP2,ish2)**2
           xtpi(6) = GM1/GA*( ZROESHu(1,IP2,ish2)*ZROESHu(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

! state 3
           xtpi(12)= ZROESHu(4,IP4,ish4)/ZROESHu(1,IP4,ish4) ! y-component
           xtpi(11)= ZROESHu(3,IP4,ish4)/ZROESHu(1,IP4,ish4) ! x-component
           xtpi(9) = ZROESHu(1,IP4,ish4)*ZROESHu(1,IP4,ish4) ! density
           HELP    = ZROESHu(3,IP4,ish4)**2+ZROESHu(4,IP4,ish4)**2
           xtpi(10)= GM1/GA*( ZROESHu(1,IP4,ish4)*ZROESHu(2,IP4,ish4)
     .              -0.5d0*HELP) ! pressure

! state 4
           xtpi(16)= ZROESHd(4,IP2,ish2)/ZROESHd(1,IP2,ish2) ! y-component
           xtpi(15)= ZROESHd(3,IP2,ish2)/ZROESHd(1,IP2,ish2) ! x-component
           xtpi(13)= ZROESHd(1,IP2,ish2)*ZROESHd(1,IP2,ish2) ! density
           HELP    = ZROESHd(3,IP2,ish2)**2+ZROESHd(4,IP2,ish2)**2
           xtpi(14)= GM1/GA*( ZROESHd(1,IP2,ish2)*ZROESHd(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

! state 5
           xtpi(20)= ZROESHd(4,IP4,ish4)/ZROESHd(1,IP4,ish4) ! y-component
           xtpi(19)= ZROESHd(3,IP4,ish4)/ZROESHd(1,IP4,ish4) ! x-component
           xtpi(17)= ZROESHd(1,IP4,ish4)*ZROESHd(1,IP4,ish4) ! density
           HELP    = ZROESHd(3,IP4,ish4)**2+ZROESHd(4,IP4,ish4)**2
           xtpi(18)= GM1/GA*( ZROESHd(1,IP4,ish4)*ZROESHd(2,IP4,ish4)
     .              -0.5d0*HELP) ! pressure

! slopes of the shocks

! shock 1 (incident shock 1)
          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)
          xtpi(21)=1.0

! shock 2 (reflected shock 1)
          DX=VSHNOR(1, IP2, ISH2)
          DY=VSHNOR(2, IP2, ISH2)
          xtpi(22)=atan2( -DX,DY)

! shock 3 (incident shock 2)
          DX=VSHNOR(1, IP3, ISH3)
          DY=VSHNOR(2, IP3, ISH3)
          xtpi(23)=1.0

! shock 4 (reflected shock 2)
          DX=VSHNOR(1, IP4, ISH4)
          DY=VSHNOR(2, IP4, ISH4)
          xtpi(24)=atan2( -DX,DY)

!         do i=1,24
!           xtp(i)=xtpi(i)
!           write(*,*)i,xtp(i)
!         enddo

          call co_uqp(xtpi,WQPX,WQPY,xtp)

!  assign values zone 1 ..
         z1v = sqrt(xtp(1))
         kine = 0.5d0*(xtp( 3)*xtp( 3)+xtp( 4)*xtp( 4))
         hh = GA/GM1*xtp( 2)/xtp(1)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 3)
         z4v = z1v*xtp( 4)

!  assign values zone 2 ..
         z1v = sqrt(xtp(5))
         kine = 0.5d0*(xtp( 7)*xtp( 7)+xtp( 8)*xtp( 8))
         hh = GA/GM1*xtp( 6)/xtp(5)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 7)
         z4v = z1v*xtp( 8)

! .. to upstream of the reflected shock 1
         ZROESHu(1,IP2,ISH2)=z1v
         ZROESHu(2,IP2,ISH2)=z2v
         ZROESHu(3,IP2,ISH2)=z3v
         ZROESHu(4,IP2,ISH2)=z4v

!  assign values zone 3 ..
         z1v = sqrt(xtp(9))
         kine = 0.5d0*(xtp( 11)*xtp( 11)+xtp( 12)*xtp( 12))
         hh = GA/GM1*xtp( 10)/xtp(9)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 11)
         z4v = z1v*xtp( 12)

! .. to upstream of the reflected shock 2
         ZROESHu(1,IP4,ISH4)=z1v
         ZROESHu(2,IP4,ISH4)=z2v
         ZROESHu(3,IP4,ISH4)=z3v
         ZROESHu(4,IP4,ISH4)=z4v

!  assign values zone 4 ..
         z1v = sqrt(xtp(13))
         kine = 0.5d0*(xtp( 15)*xtp( 15)+xtp( 16)*xtp( 16))
         hh = GA/GM1*xtp( 14)/xtp(13)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 15)
         z4v = z1v*xtp( 16)

! .. to downstream of the reflected shock1
         ZROESHd(1,IP2,ISH2)=z1v
         ZROESHd(2,IP2,ISH2)=z2v
         ZROESHd(3,IP2,ISH2)=z3v
         ZROESHd(4,IP2,ISH2)=z4v

! .. to upstream of the contact discontinuity
         ZROESHd(1,IP5,ISH5)=z1v
         ZROESHd(2,IP5,ISH5)=z2v
         ZROESHd(3,IP5,ISH5)=z3v
         ZROESHd(4,IP5,ISH5)=z4v

!  assign values zone 5 ..
         z1v = sqrt(xtp(17))
         kine = 0.5d0*(xtp( 19)*xtp( 19)+xtp( 20)*xtp( 20))
         hh = GA/GM1*xtp( 18)/xtp(17)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 19)
         z4v = z1v*xtp( 20)

! .. to downstream of the reflected shock 2
         ZROESHd(1,IP4,ISH4)=z1v
         ZROESHd(2,IP4,ISH4)=z2v
         ZROESHd(3,IP4,ISH4)=z3v
         ZROESHd(4,IP4,ISH4)=z4v

! .. to downstream of the contact discontinuity
         ZROESHu(1,IP5,ISH5)=z1v
         ZROESHu(2,IP5,ISH5)=z2v
         ZROESHu(3,IP5,ISH5)=z3v
         ZROESHu(4,IP5,ISH5)=z4v

!  assign velocity of triple point to ..

! .... incident shock 1
!        WSH(1,IP1,ISH1)= 0.0
!        WSH(2,IP1,ISH1)= 0.0

! .... reflected shock 1
         WSH(1,IP2,ISH2)= 0.0
         WSH(2,IP2,ISH2)= 0.0

! .... incident shock 2
c        WSH(1,IP3,ISH3)= 0.0
c        WSH(2,IP3,ISH3)= 0.0

! .... reflected shock 1
         WSH(1,IP4,ISH4)= 0.0
         WSH(2,IP4,ISH4)= 0.0

! .... to the contact discontinuity
         WSH(1,IP5,ISH5)= 0.0
         WSH(2,IP5,ISH5)= 0.0

! if regular reflection point along the line Y=cost
        ELSEIF(TypeSpecPoints(ISPPNTS).EQ.'RRX'.OR.
     &         TypeSpecPoints(ISPPNTS).EQ.'RR')THEN

! determine shock and indices of extremal points

! incident shock
          ISH1 = SHinSPPs(1,1,ISPPNTS)
          I    = SHinSPPs(2,1,ISPPNTS)-1
          IP1  = 1+I*(nShockPoints(ISH1)-1)

! reflected shock
          ISH2 = SHinSPPs(1,2,ISPPNTS)
          I    = SHinSPPs(2,2,ISPPNTS)-1
          IP2  = 1+I*(nShockPoints(ISH2)-1)

! recover states

! state 1
           xtpi(4) = ZROESHu(4,IP1,ish1)/ZROESHu(1,IP1,ish1) ! y-component
           xtpi(3) = ZROESHu(3,IP1,ish1)/ZROESHu(1,IP1,ish1) ! x-component
           xtpi(1) = ZROESHu(1,IP1,ish1)*ZROESHu(1,IP1,ish1) ! density
           HELP    = ZROESHu(3,IP1,ish1)**2+ZROESHu(4,IP1,ish1)**2
           xtpi(2) = GM1/GA*( ZROESHu(1,IP1,ish1)*ZROESHu(2,IP1,ish1)
     .              -0.5d0*HELP) ! pressure

! state 2
!          xtpi(8) = ZROESHu(4,IP2,ish2)/ZROESHu(1,IP2,ish2) ! y-component
!          xtpi(7) = ZROESHu(3,IP2,ish2)/ZROESHu(1,IP2,ish2) ! x-component
!          xtpi(5) = ZROESHu(1,IP2,ish2)*ZROESHu(1,IP2,ish2) ! density
!          HELP    = ZROESHu(3,IP2,ish2)**2+ZROESHu(4,IP2,ish2)**2
!          xtpi(6) = GM1/GA*( ZROESHu(1,IP2,ish2)*ZROESHu(2,IP2,ish2)
!    .              -0.5d0*HELP) ! pressure

           xtpi(8) = zroeshd(4,ip1,ish1)/zroeshd(1,ip1,ish1) ! y-component
           xtpi(7) = zroeshd(3,ip1,ish1)/zroeshd(1,ip1,ish1) ! x-component
           xtpi(5) = zroeshd(1,ip1,ish1)*zroeshd(1,ip1,ish1) ! density
           help    = zroeshd(3,ip1,ish1)**2+zroeshd(4,ip1,ish1)**2
           xtpi(6) = gm1/ga*( zroeshd(1,ip1,ish1)*zroeshd(2,ip1,ish1)
     .              -0.5d0*help) ! pressure

!          write(*,*)'stato 1'
!          write(*,*)'z1:',zroeshu(1,ip1,ish1)
!          write(*,*)'z2:',zroeshu(2,ip1,ish1)
!          write(*,*)'z3:',zroeshu(3,ip1,ish1)
!          write(*,*)'z4:',zroeshu(4,ip1,ish1)

!          write(*,*)'stato 2'
!          write(*,*)'z1:',zroeshu(1,ip2,ish2),zroeshd(1,ip1,ish1)
!          write(*,*)'z2:',zroeshu(2,ip2,ish2),zroeshd(2,ip1,ish1)
!          write(*,*)'z3:',zroeshu(3,ip2,ish2),zroeshd(3,ip1,ish1)
!          write(*,*)'z4:',zroeshu(4,ip2,ish2),zroeshd(4,ip1,ish1)

!          zroeshd(1,ip2,ish2)=2.1
!          zroeshd(2,ip2,ish2)=8.95
!          zroeshd(3,ip2,ish2)=0.1
!          zroeshd(4,ip2,ish2)=0.1

!          write(*,*)"state 3"
!          write(*,*)'z1:',zroeshd(1,ip2,ish2)
!          write(*,*)'z2:',zroeshd(2,ip2,ish2)
!          write(*,*)'z3:',zroeshd(3,ip2,ish2)
!          write(*,*)'z4:',zroeshd(4,ip2,ish2)
caldo
!          zroeshd(4,ip2,ish2)=
!    +      (zroeshd(3,ip2,ish2)+zroeshd(4,ip2,ish2))*0.5
!          zroeshd(3,ip2,ish2)=zroeshd(4,ip2,ish2)
!          zroeshd(1,ip2,ish2)=2.1
!          zroeshd(2,ip2,ish2)=8.95
!          zroeshd(3,ip2,ish2)=0.1
!          zroeshd(4,ip2,ish2)=0.1
caldo

! state 3
           xtpi(12)= ZROESHd(4,IP2,ish2)/ZROESHd(1,IP2,ish2) ! y-component
           xtpi(11)= ZROESHd(3,IP2,ish2)/ZROESHd(1,IP2,ish2) ! x-component
           xtpi(9) = ZROESHd(1,IP2,ish2)*ZROESHd(1,IP2,ish2) ! density
           HELP    = ZROESHd(3,IP2,ish2)**2+ZROESHd(4,IP2,ish2)**2
           xtpi(10)= GM1/GA*( ZROESHd(1,IP2,ish2)*ZROESHd(2,IP2,ish2)
     .              -0.5d0*HELP) ! pressure

! slopes of shocks

! shock 1 (incident shock)
          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)
          xtpi(13)=atan2( -DX,DY)

! shock 2 (reflected shock)
          DX=VSHNOR(1, IP2, ISH2)
          DY=VSHNOR(2, IP2, ISH2)
          xtpi(14)=atan2( -DX,DY)

! recover normal velocity of the incident shock
          WWS=SQRT(WSH(1,IP1,ISH1)**2+WSH(2,IP1,ISH1)**2)

! normal of shock 1 (incident shock)
          DX=VSHNOR(1, IP1, ISH1)
          DY=VSHNOR(2, IP1, ISH1)

! versor tangent to the wall
! if RRX
          taux=1.0d0
          tauy=0.0d0

!??????????????????????????????????
! in the general case RR
          IF(TypeSpecPoints(ISPPNTS).EQ.'RR')THEN

           xi= XYSH(1,IP1,ISH1)
           yi= XYSH(2,IP1,ISH1)

! find boundary segment where the shock point moves
           clr=5
           DO clr=1,NCLR
            IF(ICLR(clr).eq.ispclr(ISPPNTS))exit
           enddo
          BBGN=IA(clr)
          BEND=IA(clr+1)-1
          DO J = BBGN, BEND-1
             K   = JA(J)
             KP1 = JA(J+1)
             dumX1=CORG(1,K)
             dumY1=CORG(2,K)
             dumX2=CORG(1,KP1)
             dumY2=CORG(2,KP1)

!            if(dumx2.LT.dumx1)THEN
!             dum=dumx1
!             dumX1=dumX2
!             dumX2=dum
!            ENDIF
!
!            if(dumy2.LT.dumy1)THEN
!             dum=dumy1
!             dumY1=dumY2
!             dumY2=dum
!            ENDIF

!            if(xi.le.dumX2.and.
!    +          xi.ge.dumX1.and.
!    +          yi.le.dumY2.and.
!    +          yi.ge.dumY1)     THEN

               dum=sqrt((dumX1-dumX2)**2+(dumY1-dumY2)**2)
               dum1=sqrt((xi-dumX2)**2+(yi-dumY2)**2)+
     +              sqrt((xi-dumX1)**2+(yi-dumY1)**2)
              if(dum1-dum.le.1.0e-5)then

                taux=-(dumX1-dumX2)
                tauy=-(dumY1-dumY2)
                dum=sqrt(taux**2+tauy**2)
                taux=taux/dum
                tauy=tauy/dum
             endif

          ENDDO

          ENDIF
!         write(*,*)
!         write(*,*)'taux,tauy:',taux,tauy
!         write(*,*)
!????????????????????????????????

          dum=DX*taux+DY*tauy

! calculate the velocity of the reflected point
          dum1=DX*WSH(1,IP1,ISH1)+DY*WSH(2,IP1,ISH1)
          WRR=sign(1.d0,dum1)*WWS/dum

! set velocity of the reflected point
          xtpi(15)=WRR
!         xtpi(15)=2.5

!         write(*,*)
!         write(*,*)'WRR:',WRR
!         do I=1,15
!          write(*,*)I,xtpi(i)
!         enddo

!         call subroutine to compute the regular reflection
          call co_urr(xtpi,taux,tauy,xtp)

!         do I=1,15
!          write(*,*)I,xtp(i)
!         enddo

! assign values zone 1 ..
         z1v = sqrt(xtp(1))
         kine = 0.5d0*(xtp( 3)*xtp( 3)+xtp( 4)*xtp( 4))
         hh = GA/GM1*xtp( 2)/xtp(1)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 3)
         z4v = z1v*xtp( 4)

! .. to upstream of the incident shock
         zroeshu(1,ip1,ish1)=z1v
         zroeshu(2,ip1,ish1)=z2v
         zroeshu(3,ip1,ish1)=z3v
         zroeshu(4,ip1,ish1)=z4v

! assign values of zone 2 ..
         z1v = sqrt(xtp(5))
         kine = 0.5d0*(xtp( 7)*xtp( 7)+xtp( 8)*xtp( 8))
         hh = GA/GM1*xtp( 6)/xtp(5)+kine
         z2v = z1v*hh
         z3v = z1v*xtp( 7)
         z4v = z1v*xtp( 8)

! .. to downstream of the incident shock
         zroeshd(1,ip1,ish1)=z1v
         zroeshd(2,ip1,ish1)=z2v
         zroeshd(3,ip1,ish1)=z3v
         zroeshd(4,ip1,ish1)=z4v

! .. to upstream of the reflected shock
         zroeshu(1,ip2,ish2)=z1v
         zroeshu(2,ip2,ish2)=z2v
         zroeshu(3,ip2,ish2)=z3v
         zroeshu(4,ip2,ish2)=z4v

! assign values of zone 3 ..
         z1v = sqrt(xtp(9))
         kine = 0.5d0*(xtp(11)*xtp(11)+xtp(12)*xtp(12))
         hh = ga/gm1*xtp(10)/xtp(9)+kine
         z2v = z1v*hh
         z3v = z1v*xtp(11)
         z4v = z1v*xtp(12)

! .. to downstream of the reflected shock
         zroeshd(1,ip2,ish2)=z1v
         zroeshd(2,ip2,ish2)=z2v
         zroeshd(3,ip2,ish2)=z3v
         zroeshd(4,ip2,ish2)=z4v

!        write(*,*)"z(1)_new:",zroeshu(1,ip1,ish1)
!        write(*,*)"z(2)_new:",zroeshu(2,ip1,ish1)
!        write(*,*)"z(3)_new:",zroeshu(3,ip1,ish1)
!        write(*,*)"z(4)_new:",zroeshu(4,ip1,ish1)

!        write(*,*)"z(1)_new:",zroeshd(1,ip1,ish1)
!        write(*,*)"z(2)_new:",zroeshd(2,ip1,ish1)
!        write(*,*)"z(3)_new:",zroeshd(3,ip1,ish1)
!        write(*,*)"z(4)_new:",zroeshd(4,ip1,ish1)

!        write(*,*)"z(1)_new:",zroeshu(1,ip2,ish2)
!        write(*,*)"z(2)_new:",zroeshu(2,ip2,ish2)
!        write(*,*)"z(3)_new:",zroeshu(3,ip2,ish2)
!        write(*,*)"z(4)_new:",zroeshu(4,ip2,ish2)

!        write(*,*)"z(1)_new:",zroeshd(1,ip2,ish2)
!        write(*,*)"z(2)_new:",zroeshd(2,ip2,ish2)
!        write(*,*)"z(3)_new:",zroeshd(3,ip2,ish2)
!        write(*,*)"z(4)_new:",zroeshd(4,ip2,ish2)

! sum the contribution along incident shock to the velocity of triple point

! Vale
         WRR=xtp(15)
! Vale

         wsh(1,ip1,ish1)= wrr*taux
         wsh(2,ip1,ish1)= wrr*tauy

!       write(*,*)"WRR",WRR*taux,WRR*tauy
!       write(*,*)

! propagate to the other shock points the new velocity
         wsh(1,ip2,ish2)= wsh(1,ip1,ish1)
         wsh(2,ip2,ish2)= wsh(2,ip1,ish1)

! if termination point
        ELSEIF(TypeSpecPoints(ISPPNTS).EQ.'EP')THEN

! determine shock and indices of the extremum and internal point
          ish=shinspps(1,1,isppnts)
          i=  shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish)-1)
          ip1=2+i*(nshockpoints(ish)-3)

! recover normal and tangential versor of the internal point
          nx=vshnor(1, ip1, ish)
          ny=vshnor(2, ip1, ish)
          taux=ny
          tauy=-nx

! calculate geometric tangential vector
          dx=xysh(1,ip,ish)-xysh(1,ip1,ish)
          dy=xysh(2,ip,ish)-xysh(2,ip1,ish)

! check direction of tangential versors and if necessary rectify the
! direction of the versor determined with the normal to the internal point
          dum=taux*dx+tauy*dy
          if(dum.lt.0.d+0)then
            taux=-taux
            tauy=-tauy
          endif

! determine the new position of the termination point
          xysh(1,ip,ish)=xysh(1,ip1,ish)+taux*dxcell
          xysh(2,ip,ish)=xysh(2,ip1,ish)+tauy*dxcell

! assign the discontinuity velocity equal to the velocity of the
! internal point
         wsh(1,ip,ish)= wsh(1,ip1,ish)
         wsh(2,ip,ish)= wsh(2,ip1,ish)

! fix the upstream and downstream state by assigning the average of the
! upstream and downstream state of the internal point
         zroeshd(1,ip,ish)=0.5d0*(zroeshd(1,ip1,ish)+zroeshu(1,ip1,ish))
         zroeshd(2,ip,ish)=0.5d0*(zroeshd(2,ip1,ish)+zroeshu(2,ip1,ish))
         zroeshd(3,ip,ish)=0.5d0*(zroeshd(3,ip1,ish)+zroeshu(3,ip1,ish))
         zroeshd(4,ip,ish)=0.5d0*(zroeshd(4,ip1,ish)+zroeshu(4,ip1,ish))

         zroeshu(1,ip,ish)=0.5d0*(zroeshd(1,ip1,ish)+zroeshu(1,ip1,ish))
         zroeshu(2,ip,ish)=0.5d0*(zroeshd(2,ip1,ish)+zroeshu(2,ip1,ish))
         zroeshu(3,ip,ish)=0.5d0*(zroeshd(3,ip1,ish)+zroeshu(3,ip1,ish))
         zroeshu(4,ip,ish)=0.5d0*(zroeshd(4,ip1,ish)+zroeshu(4,ip1,ish))

! if shock point is a start point
        elseif(typespecpoints(isppnts).eq.'SP')then

! determine shock and index of the extreme and internal point
          ish=shinspps(1,1,isppnts)
          i=  shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish)-1)
          ip1=2+i*(nshockpoints(ish)-3)

          wsh(1,ip,ish)= 0.0
          wsh(2,ip,ish)= 0.0

! the point is moved in a different way

! if the point if a connection or periodic connection
         elseif(typespecpoints(isppnts).eq.'C'
     +      .or.typespecpoints(isppnts).eq.'PC')then

! determine shock and indices of extremal points
! point  1
         ish1 = shinspps(1,1,isppnts)
         i    = shinspps(2,1,isppnts)-1
         ip1  = 1+i*(nshockpoints(ish1)-1)

! point 2
         ish2 = shinspps(1,2,isppnts)
         i = shinspps(2,2,isppnts)-1
         ip2 = 1+i*(nshockpoints(ish2)-1)

         if(typespecpoints(isppnts).eq.'C')then

! set coordinates of point 2 with point 1
          xysh(1,ip2,ish2)=xysh(1,ip1,ish1)
          xysh(2,ip2,ish2)=xysh(2,ip1,ish1)
         endif

! set shock velocity of point 2 with point 1
         wsh(1,ip2,ish2)= wsh(1,ip1,ish1)
         wsh(2,ip2,ish2)= wsh(2,ip1,ish1)

! set state in point 2 with point 1
         zroeshd(1,ip2,ish2)=zroeshd(1,ip1,ish1)
         zroeshd(2,ip2,ish2)=zroeshd(2,ip1,ish1)
         zroeshd(3,ip2,ish2)=zroeshd(3,ip1,ish1)
         zroeshd(4,ip2,ish2)=zroeshd(4,ip1,ish1)

         zroeshu(1,ip2,ish2)=zroeshu(1,ip1,ish1)
         zroeshu(2,ip2,ish2)=zroeshu(2,ip1,ish1)
         zroeshu(3,ip2,ish2)=zroeshu(3,ip1,ish1)
         zroeshu(4,ip2,ish2)=zroeshu(4,ip1,ish1)

         else
          write(*,*)'condition not defined'
          write(8,*)'condition not defined'
          stop
         endif

      end do

!     calculate variation of the shock points state
      wws=0.0
      do iv=1, ndof
       varz(iv)=0.0
       avarz(iv)=0.0
      end do
      do ish=1,nshocks
      write(8,*)'variations of z downstream of the shock'
       write(8,*)'shock n.',ish
       do i=1, nshockpoints(ish)
        do iv=1, ndof
         varz(iv)=zroeshd(iv,i,ish)-zroeshdold(iv,i,ish)
         avarz(iv)=avarz(iv)+sqrt(varz(iv)**2)
        end do
        wws=wws+sqrt(wsh(1,i,ish)**2+wsh(2,i,ish)**2)

        write(8,'(1x,i3,5(1x,f15.7))')i,(varz(iv),iv=1,ndof)
       enddo
      enddo

      write(8,'(1x,a19,5(1x,f15.7))')'average values var z:',
     .                           (avarz(iv),iv=1,ndof)
      write(87,'(1x,i5,5(1x,f15.7))')iter,(avarz(iv),iv=1,ndof),wws

!      save the old upstream state and assign
!      the recomputed shock upstream state
       do ish=1,nshocks
        do ip=1, nshockpoints(ish)
         do iv=1,ndof
          zroeshdold(iv,ip,ish)=zroeshd(iv,ip,ish)
          zroeshuold(iv,ip,ish)=zroeshu(iv,ip,ish)
         enddo
        enddo
       enddo

      close(8)

      return
      end
