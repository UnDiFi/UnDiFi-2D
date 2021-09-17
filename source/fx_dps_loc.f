! Fix and correct the nodal position in all special points

      subroutine fx_dps_loc(
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
     +           zroe,          !vale
     +           corg,
     +           shtopolchanged)!vale

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
     +                 xysh    (ndim, npshmax, *),
     +                 xyshu   (ndim, npshmax, *),
     +                 xyshd   (ndim, npshmax, *),
     +                 zroeshu    (ndof, npshmax, *),
     +                 zroeshd    (ndof, npshmax, *),
     +                 zroeshuold (ndof, npshmax, *),
     +                 zroeshdold (ndof, npshmax, *),
     +                 vshnor     (ndim, npshmax, *),
     +                 wsh        (ndim, npshmax, *)
      integer nclr
      integer ia(*),ja(*),iclr(nclr)
      double precision corg(ndim,*), zroe(ndof,*)

! vale
      double precision xywedge(2), dum1, dum2
      logical shtopolchanged
! vale

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
      integer isppnts,ip1,ish1,ip2,ish2,i,j,k,kp1
      integer clr,bbgn,bend
      double precision xi,yi,x1,y1,x2,y2,dumx1,dumy1,dumx2,dumy2
      double precision dum,xi1,yi1
      integer nn,ish
      parameter (nn=2)
      double precision a(nn,nn), b(nn),x(nn)

!     double precision dx,dy,kine,ws,hh,wws,taux,tauy,dum,wrr,wss
!     double precision wqpx,wqpy,w,cs,dum1,f1,f3,dxr14,dyr14
!     double precision help,x1(ndof),x2(ndof),taux12,tauy12,nx,ny
!     double precision r2(npshmax,nshmax),xtpi(24),xtp(24),r23,r14
!     double precision r14old
!     double precision dumm2,dumx,dumy
!     integer i,im,iv,k,totnshockpoints
!     integer isppnts
!     integer ip,ish,ip1,ip2,ip3,ip4,ip5,ish1,ish2,ish3,ish4,ish5
!     logical flag1,ifail

!     open log file
      open(8,file='log/fx_dps_loc.log')

      do isppnts=1,nspecpoints

!     if shock point is in inlet section
      if(typespecpoints(isppnts).eq.'IPX'.or.
     +   typespecpoints(isppnts).eq.'IPY')then

!     if shock point is in outlet section
      elseif(typespecpoints(isppnts).eq.'OPX'.or.
     +       typespecpoints(isppnts).eq.'OPY')then

!     if shock point on X or Y wall
        elseif(typespecpoints(isppnts).eq.'WPNRX'.or.
     +         typespecpoints(isppnts).eq.'WPNRY')then

!     if shock point is on curved wall
        elseif(typespecpoints(isppnts).eq.'FWP'.or.
     +         typespecpoints(isppnts).eq.'RR')then

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          xi = xysh(1,ip1,ish1)
          yi = xysh(2,ip1,ish1)

!         set colour of the boundary on which the shock point moves
          clr=5
          do clr=1,nclr
          if(iclr(clr).eq.ispclr(isppnts))exit
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

!            if(xi.le.dumX2.and.
!    +          xi.ge.dumX1.and.
!    +          yi.le.dumY2.and.
!    +          yi.ge.dumY1)THEN
!               write(*,*)'trovato edge'
!               write(*,*)'xi:',xi
!               write(*,*)'x1:',x1
!               write(*,*)'x2:',x2

! write the equation of the straight line passing for the two boundary points
                a(1,1)=(y2-y1)
                a(1,2)=(x1-x2)
!               b(1)=x1*(y2-y1)+y1*(x1-x2)
                b(1)=x2*(y2-y1)+y2*(x1-x2)

! write the equation of the straight line perpendicular to previous boundary line
! and passing for the vertex node
                a(2,1)=a(1,2)
                a(2,2)=-a(1,1)
                b(2)= a(2,1)*xi+a(2,2)*yi

! solve the linear system and find the intersection point coordinates
               call solg(nn,nn,a,b,x)
               xi1=x(1)
               yi1=x(2)
!              write(*,*)'xi:',xi
!              write(*,*)'x1:',x1
!              write(*,*)'x2:',x2
               dum=sqrt((xi1-xi)**2+ (yi1-yi)**2)

              if(dum.lt.dxcell*0.5.and.
     +          xi.le.dumx2.and.
     +          xi.ge.dumx1.and.
     +          yi.le.dumy2.and.
     +          yi.ge.dumy1)then

              xysh(1,ip1,ish1)=xi1
              xysh(2,ip1,ish1)=yi1

             if(typespecpoints(isppnts).eq.'RR')then

                ish1 = shinspps(1,2,isppnts)
                i    = shinspps(2,2,isppnts)-1
                ip1  = 1+i*(nshockpoints(ish1)-1)

                xysh(1,ip1,ish1)=xi1
                xysh(2,ip1,ish1)=yi1
              endif

             endif
           ENDDO
! vale
! check if the shock is after the beginning of the wedge

!           xywedge(1)=0.5d0
!           xywedge(2)=0d0

!           dum1=sqrt(xi1**2+yi1**2)
!           dum2=sqrt(xywedge(1)**2+xywedge(2)**2)

! include the velocity dependence to retrieve the verse
!           if((dum1-dum2).gt.1e-16) then

!             if((dum1-dum2).lt.dxcell*1.0) then
!             if((dum1-dum2).lt.dxcell*0.8) then
!
!             if((dum1-dum2).lt.dxcell) then
!             write(8,*)"shock after the wedge"
!             write(8,*) "topology not changed"
!             write(*,*)"shock after the wedge"
!             write(*,*) "topology not changed"
!             write(*,*)

!            else

!             if(.not.shtopolchanged) then

! change the shock topology

!             write(8,*)"Shock after the wedge"
!             write(8,*) "Topology changed"
!             write(*,*)"Shock after the wedge"
!             write(*,*) "Topology changed"
!             write(*,*)

!              call ch_sh_topology(xysh,zroeshu,zroeshd,wsh,
!     &                            xywedge,nshocks,typesh,
!     &                            nshockpoints,nspecpoints,
!     &                            typespecpoints,isppnts,
!     &                            shinspps,ispclr,
!     &                            nclr,iclr,ia,ja,
!     &                            zroe,corg)

!              shtopolchanged=.true.
!             else
!              write(8,*)"shock after the wedge"
!              write(8,*) "topology not changed"
!             endif
!            endif

!           elseif((dum1-dum2).lt.1e-13) then
!            write(8,*)"shock before the wedge"
!           endif
! vale

! if point of periodic conenction
        elseif(typespecpoints(isppnts).eq.'PC')then

! fix point 1

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          xi= xysh(1,ip1,ish1)
          yi= xysh(2,ip1,ish1)

! set colour of the boundary on which the shock point moves
          do clr=1,nclr
           if(iclr(clr).eq.ispclr(isppnts))exit
          enddo

!         write(*,*)'ispclr:',ispclr(isppnts+1)
!         write(*,*)'iclr:',iclr(clr)

          bbgn=ia(clr)
          bend=ia(clr+1)-1
!         write(*,*)'bbgn,bend:',ja(bbgn),ja(bend)
          do j = bbgn, bend-1
             k   = ja(j)
             kp1 = ja(j+1)
             x1=corg(1,k)
             y1=corg(2,k)
             x2=corg(1,kp1)
             y2=corg(2,kp1)
!            write(*,*)'x 1',j, x1

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

! write the equation of the straight line passing for the two boundary points
                a(1,1)=(y2-y1)
                a(1,2)=(x1-x2)
!               b(1)=x1*(y2-y1)+y1*(x1-x2)
                b(1)=x2*(y2-y1)+y2*(x1-x2)

! write the equation of the straight line perpendicular to previous boundary line
! and passing for the vertex node
                a(2,1)=a(1,2)
                a(2,2)=-a(1,1)
                b(2)= a(2,1)*xi+a(2,2)*yi

! solve the linear system and find the intersection point coordinates
               call solg(nn,nn,a,b,x)
               xi1=x(1)
               yi1=x(2)
               dum=sqrt((xi1-xi)**2+ (yi1-yi)**2)

              if(dum.lt.dxcell*0.5.and.
     +          xi1.le.dumx2.and.
     +          xi1.ge.dumx1.and.
     +          yi1.le.dumy2.and.
     +          yi1.ge.dumy1)then

!              write(*,*)'xysh(1,ip1,ish1):',xysh(1,ip1,ish1)
!              write(*,*)'dumx2:',dumx2
!              write(*,*)'dumx1:',dumx1
!              write(*,*)'dumy2:',dumy2
!              write(*,*)'dumy1:',dumy1
!              write(*,*)'dum:',dum
!              write(*,*)'xi,yi:',xi,yi
!              write(*,*)'xi1,yi1:',xi1,yi1
!              pause
!              continue

              xysh(1,ip1,ish1)=xi1
              xysh(2,ip1,ish1)=yi1

             endif

        enddo

!       write(*,*)
!       write(*,*)'1',xysh(1,ip1,ish1)
!       write(*,*)

! fix point 2

          ish1 = shinspps(1,2,isppnts)
          i    = shinspps(2,2,isppnts)-1
          ip1  = 1+i*(nshockpoints(ish1)-1)

          xi= xysh(1,ip1,ish1)
          yi= xysh(2,ip1,ish1)

! set colour of the boundary on which the shock point moves
          clr=5
          do clr=1,nclr
           if(iclr(clr).eq.ispclr(isppnts+1))exit    ! FIXME: correct this part of code, could be dangerous
          enddo
!         write(*,*)'ispclr:',ispclr(isppnts+1)
!         write(*,*)'iclr:',iclr(clr)

          bbgn=ia(clr)
          bend=ia(clr+1)-1
!         write(*,*)'bbgn,bend:',ja(bbgn),ja(bend)
          do j = bbgn, bend-1
             k   = ja(j)
             kp1 = ja(j+1)
             x1=corg(1,k)
             y1=corg(2,k)
             x2=corg(1,kp1)
             y2=corg(2,kp1)
!            write(*,*)'x 2',j, x1

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

!            if(xi.le.dumx2.and.
!    +          xi.ge.dumx1.and.
!    +          yi.le.dumy2.and.
!    +          yi.ge.dumy1)then
!               write(*,*)'edge found'
!               write(*,*)'xi:',xi
!               write(*,*)'x1:',x1
!               write(*,*)'x2:',x2

! write the equation of the straight line passing for the two boundary points
                a(1,1)=(y2-y1)
                a(1,2)=(x1-x2)
!               b(1)=x1*(y2-y1)+y1*(x1-x2)
                b(1)=x2*(y2-y1)+y2*(x1-x2)

! write the equation of the straight line perpendicular to previous boundary line
! and passing for the vertex node
                a(2,1)=a(1,2)
                a(2,2)=-a(1,1)
                b(2)= a(2,1)*xi+a(2,2)*yi

! solve the linear system and find the intersection point coordinates
               call solg(nn,nn,a,b,x)
               xi1=x(1)
               yi1=x(2)
!              write(*,*)'xi:',xi
!              write(*,*)'x1:',x1
!              write(*,*)'x2:',x2
               dum=sqrt((xi1-xi)**2+ (yi1-yi)**2)

              if(dum.lt.dxcell*0.5.and.
     +          xi1.le.dumx2.and.
     +          xi1.ge.dumx1.and.
     +          yi1.le.dumy2.and.
     +          yi1.ge.dumy1)then

!              write(*,*)'xysh(1,ip1,ish1):',xysh(1,ip1,ish1)
!              write(*,*)'dumx2:',dumx2
!              write(*,*)'dumx1:',dumx1
!              write(*,*)'dumy2:',dumy2
!              write(*,*)'dumy1:',dumy1
!              write(*,*)'dum:',dum
!              write(*,*)'xi,yi:',xi,yi
!              write(*,*)'xi1,yi1:',xi1,yi1
!              pause
!              continue

!             write(*,*)'xysh(1,ip2,ish2):',xysh(1,ip1,ish1)

              xysh(1,ip1,ish1)=xi1
              xysh(2,ip1,ish1)=yi1

             endif

        enddo

!       write(*,*)
!       write(*,*)'2',xysh(1,ip1,ish1)
!       write(*,*)

caldo
          ish1 = shinspps(1,1,isppnts)                    ! temporary code
          i    = shinspps(2,1,isppnts)-1                  ! temporary code
          ip1  = 1+i*(nshockpoints(ish1)-1)               ! temporary code

          ish2 = shinspps(1,2,isppnts)                    ! temporary code
          i    = shinspps(2,2,isppnts)-1                  ! temporary code
          ip2  = 1+i*(nshockpoints(ish1)-1)               ! temporary code

          xi= 0.5*(xysh(1,ip1,ish1)+xysh(1,ip2,ish2))     ! temporary code

!         xysh(1,ip1,ish1)=xi                             ! temporary code
!         xysh(1,ip2,ish2)=xi                             ! temporary code

caldo

! if triple point
        elseif(typespecpoints(isppnts).eq.'TP')then

! if quadruple point
        elseif(typespecpoints(isppnts).eq.'QP')then

! if regular reflection point along line y=cost
        elseif(typespecpoints(isppnts).eq.'RRX')then

! if termination point
        elseif(typespecpoints(isppnts).eq.'EP')then

! if sonic point moving along x direction
        elseif(typespecpoints(isppnts).eq.'SP')then

! if connection point
        elseif(typespecpoints(isppnts).eq.'C')then

! if trailing edge point
        elseif(typespecpoints(isppnts).eq.'TE')then

         else
          write(*,*)'condition not defined'
          write(8,*)'condition not defined'
          stop
         endif

      end do

!    Save old upstream state and assign recomputed shock upstream state

!       do ish=1,nshocks
!        do ip=1, nshockpoints(ish)
!         do iv=1,ndof
!          zroeshdold(iv,ip,ish)=zroeshd(iv,ip,ish)
!          zroeshuold(iv,ip,ish)=zroeshu(iv,ip,ish)
!         enddo
!        enddo
!       enddo

!      do ish=1,nshocks
!       if(typesh(ish).eq.'s'.and.iter.le.1001.and.mod(iter,50).eq.0)then
!
!        j=nshockpoints(ish)
!        xysh(1,j+1,ish)=2*xysh(1,j,ish)-xysh(1,j-1,ish)
!        xysh(2,j+1,ish)=2*xysh(2,j,ish)-xysh(2,j-1,ish)
!        zroeshu(1,j+1,ish)=zroeshu(1,j,ish)
!        zroeshu(2,j+1,ish)=zroeshu(2,j,ish)
!        zroeshu(3,j+1,ish)=zroeshu(3,j,ish)
!        zroeshu(4,j+1,ish)=zroeshu(4,j,ish)
!
!        zroeshd(1,j+1,ish)=zroeshd(1,j,ish)
!        zroeshd(2,j+1,ish)=zroeshd(2,j,ish)
!        zroeshd(3,j+1,ish)=zroeshd(3,j,ish)
!        zroeshd(4,j+1,ish)=zroeshd(4,j,ish)

!        zroeshu(1,j,ish)=zroeshu(1,j-1,ish)
!        zroeshu(2,j,ish)=zroeshu(2,j-1,ish)
!        zroeshu(3,j,ish)=zroeshu(3,j-1,ish)
!        zroeshu(4,j,ish)=zroeshu(4,j-1,ish)
!
!        zroeshd(1,j,ish)=zroeshd(1,j-1,ish)
!        zroeshd(2,j,ish)=zroeshd(2,j-1,ish)
!        zroeshd(3,j,ish)=zroeshd(3,j-1,ish)
!        zroeshd(4,j,ish)=zroeshd(4,j-1,ish)
!
!        nshockpoints(ish)=j+1
!       endif
!      enddo

      return
      end
