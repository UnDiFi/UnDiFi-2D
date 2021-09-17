! Updates values in the phantom nodes of the background mesh(0)
! interpolating values in the shocked mesh (1)

      subroutine interp(
     .           ibndfac,        !not used
     .           nbfac,          !not used
     .           icelnod,
     .           nvt,
     .           nelem,
     .           xy,
     .           zroe,
     .           xysh,
     .           xyshu,
     .           xyshd,
!    .           work,
     .           nphpoin,        !not used
     .           xybkg,
     .           zbkg,
     .           nodcod,
     .           npoin,
     .           nshocks,
     .           nshockpoints,
     +           ia,
     +           ja,
     +           iclr,
     +           nclr)

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer nelem,nvt,nbfac
      integer nshocks,nshockpoints(nshmax),nphpoin

!     .. array arguments ..
      double precision xysh(ndim,npshmax,*),
     +                 xyshu(ndim,npshmax,*),
     +                 xyshd(ndim,npshmax,*),
     +                 xy(ndim,*),
!    &                 work(ndim,npshmax,*),
     &                 zroe(ndof,*),
     +                 xybkg(ndim,*),
     +                 zbkg(ndof,*)

      integer          ibndfac(3,*),
     +                 icelnod(nvt,nelem),
     +                 nodcod(*),
     +                 npoin(0:*)

      integer nclr
      integer ia(*),ja(*),iclr(nclr)

!     .. local scalars ..
      integer ipoin,ielem,i,ii,k,n,ifail,ish,clr,bbgn,bend,j,kp1,ibc
      double precision x0,y0,x1,y1,x2,y2,dum,dum1,dum2

!     open log file
      open(8,file='log/interp.log')

      write(8,*)'enter in interp'

!     make upstream node coordinates coincide with downstream node coordinates
!     these are the coordinates of the shocked mesh (1)
!     Note: the nof of shock points is that on the shocked mesh
!     not the one on the background mesh since this one might have
!     been updated in the shock redistribution routine called by shockmov
      do ish = 1, nshocks
        do ipoin = 1, nshockpoints(ish)
          xyshu(1,ipoin,ish) = xysh(1,ipoin,ish)
          xyshu(2,ipoin,ish) = xysh(2,ipoin,ish)
          xyshd(1,ipoin,ish) = xysh(1,ipoin,ish)
          xyshd(2,ipoin,ish) = xysh(2,ipoin,ish)
        enddo
      enddo

!     interpolate in the background mesh nodes, using the connectivity of the shocked mesh;
!     the interpolation is necessary only for the ghost nodes
      do ipoin = 1, npoin(0)
!       if((nodcod(ipoin).eq.-1) or.(nodcod(ipoin).eq.-2)) then
        if((nodcod(ipoin).eq.-1)) then
          write(8,*)'trying to locate ',ipoin,
     &'(',xybkg(1,ipoin),',',xybkg(2,ipoin),')'
             ifail=0
!            if(ipoin.eq.1449)ifail=99
             call finder(icelnod,nelem,xy,ndim,zroe,ndof,xybkg(1,ipoin),
     &                   zbkg(1,ipoin),ielem,ifail)
             write(8,*)'found in cell ',ielem,ifail
             if(ifail.ne.0)then
                write(8,*)'search failed for vertex ',ipoin
                write(8,*)(xybkg(i,ipoin),i=1,ndim)
                write(8,*)'cell no is ',ielem
                stop
             endif
         else
!           call dcopy(ndof,zroe(1,ipoin),1,zbkg(1,ipoin),1)
             do i = 1, ndof
                zbkg(i,ipoin) = zroe(i,ipoin)
             enddo
         endif
      enddo

!     interpolate phantom node on the boundary
      do ipoin = 1, npoin(0)
        if((nodcod(ipoin).eq.-2))then
          write(8,*)'trying to locate bounday point',ipoin,
     &'(',xybkg(1,ipoin),',',xybkg(2,ipoin),')'

         x0=xybkg(1,ipoin)
         y0=xybkg(2,ipoin)
!        write(*,*)'x0,y0:',x0,y0

         do i=1,nbfac

           k   = ibndfac(1,i)
           kp1 = ibndfac(2,i)
           ibc = ibndfac(3,i)
           if(ibc.lt.10.and.ibc.gt.0)then

             x1=xy(1,k)
             y1=xy(2,k)
             x2=xy(1,kp1)
             y2=xy(2,kp1)

             dum  = ((x2-x1)**2+(y2-y1)**2)
             dum1 = ((x0-x1)**2+(y0-y1)**2)
             dum2 = ((x0-x2)**2+(y0-y2)**2)

              if((dum1+dum2).le.dum)then
                write(8,*)'search succesfully for vertex ',ipoin
                write(8,*)(xybkg(ii,ipoin),ii=1,ndim)

               dum=sqrt(dum1)+sqrt(dum2)
               dum1=sqrt(dum1)/dum
               dum2=sqrt(dum2)/dum

               write(8,*)'dum:',dum,'dum1:',dum1,'dum2:',dum2
               do j = 1, ndof
                 zbkg(j,ipoin) = dum2*zroe(j,k)+dum1*zroe(j,kp1)
                 write(8,*)'k:',k,'kp1:',kp1
                 write(8,*)'z',j,':',zbkg(j,ipoin),zroe(j,k),zroe(j,kp1),

               enddo

              endif
             endif

         enddo

       endif

      enddo

!    goto 65
!    do ish=1,nshocks
!     do ipoin = 1, nshockpoints(ish)
!        k = ipoin ! shock point
!        n = k + nshockpoints(ish) ! duplicated shock point
!        xysh(1,n,ish) = work(1,ipoin,ish)
!        xysh(2,n,ish) = work(2,ipoin,ish)
!     enddo
!     enddo
   65 continue

      write(8,*)'exit interp'

      close(8)

      return
      end

      subroutine finder(icelnod,nelem,coor,ndim,zroe,ndof,xyin,zout,
     &ielem,info)

!     input:
!            xyin nodal coords (belonging to the background grid)
!                 to be located inside the current grid
!            icelnod cell to node pointer
!            coor nodal coordinates
!            ndim space dimension
!            ndof nof degrees of freedom
!     output:
!            ielem is the cell node xyin falls inside
!            info = 0 node found !=0 search failed
!            zout(*) is filled with the interpolated value

      implicit none
      integer ielem,ndim,ndof,info,nelem,info1
      integer icelnod(3,*)
      double precision coor(ndim,*),zroe(ndof,*),aa
      double precision xyin(ndim),zout(ndof)
      double precision x0,y0,xp(4),yp(4),a(3),t,s,help
      integer idxs(3),iv,ipoin,ivar,ilog
      double precision eps
      parameter(eps=1.e-08,ilog=1)
      double precision area
      integer icycl

      info1 = info
      info = 0

!     x0,y0 are the coordinates of the point to be located
      x0 = xyin(1)
      y0 = xyin(2)

      do 1000 ielem = 1, nelem
         do 50 iv = 1,3
            ipoin = icelnod(iv,ielem) ! global node number
            xp(iv) = coor(1,ipoin)
            yp(iv) = coor(2,ipoin)
   50    continue
         xp(4) = x0 ! node to be located
         yp(4) = y0 ! node to be located

         idxs(1) = 1
         idxs(2) = 2
         idxs(3) = 3
         help = 1.d0/area(xp,yp,3,idxs)
         idxs(3) = 4 ! node to be located

!        compute area coordinates (in the x-y plane)
         do iv = 1,3
            idxs(1) = icycl(1+iv,3)
            idxs(2) = icycl(2+iv,3)
            a(iv) = area(xp,yp,3,idxs)*help
!           if(abs(a(iv)).le.eps)a(iv)=0.d0
         enddo

         s = min( a(1), a(2), a(3) )
         t = max( a(1), a(2), a(3) )
         aa=abs(a(1))+ abs(a(2))+ abs(a(3))
         if(info1.eq.99)then
!          if(t.le.2.0.and.s.ge.-1.0)then
!          if(aa.le.2.5.and.aa.gt.0.5)then
           write(*,*)'x0,y0:',x0,y0
           write(*,*)'ielem',ielem
           write(*,*)t,s,aa
           write(*,*)xp(1), yp(1)
           write(*,*)xp(2), yp(2)
           write(*,*)xp(3), yp(3)

           write(*,*)
!          endif
         endif
!        write(6,*)'area coords ',(a(j),j=1,3),' cell ',ielem,help
!        write(6,*)'y/z coords ',(xp(j),j=1,4),(yp(j),j=1,4)
         if( ( s .ge. 0.d0 .and. s .le. 1.d0 ) .and.
     +       ( t .ge. 0.d0 .and. t .le. 1.d0 ) )then
         if(ilog.eq.0)write(6,fmt=*)(icelnod(iv,ielem),iv=1,3)
          do ivar = 1, ndof
             zout(ivar) =  0.d0
          enddo
          do iv = 1, 3
             ipoin = icelnod(iv,ielem)
             help = a(iv)
             do ivar = 1, ndof
                zout(ivar) = zout(ivar) + help* zroe(ivar,ipoin)
             enddo
          enddo
          if(ilog.eq.0)then
             do ivar = 1, ndim
             write(6,fmt=400)xyin(ivar),
     &(coor(ivar,icelnod(iv,ielem)),iv=1,3),(a(iv),iv=1,3),s,t
             enddo
             do ivar = 1, ndof
             write(6,fmt=300)zout(ivar),
     &(zroe(ivar,icelnod(iv,ielem)),iv=1,3),(a(iv),iv=1,3)
             enddo
          endif
          info = 0
          return
          endif
 1000 continue
      info = 1
      write(6,*)'search failed for vertex coords ',x0,y0
      write(6,fmt=1100)(a(iv),iv=1,3),a(1)+a(2)+a(3),s,t
      return
 1100 format('a(i),s,s,t ',6(e12.4,1x))
  300 format(7(f10.5,1x))
  400 format(9(f10.5,1x))
      end

! *************************************************************
! Updates values in the phantom nodes of the background mesh(0)
! interpolating values in the shocked mesh (1)
! *************************************************************

      subroutine interp_sp(
     .             icelnod,
     .             nvt,
     .             nelem,
     .             xy,
     .             zroe,
     .             xysh,
     +             zroeshu, !upstream
     +             zroesh,  !downstream
     +             vshnor,
     .             npoin,
     .             nshocks,
     .             nshockpoints,
     +             typesh,
     +             nspecpoints,
     +             typespecpoints,
     +             shinspps)

      implicit none
      include 'paramt.h'

!     .. scalar arguments ..
      integer nelem,nvt,nbfac
      integer nshocks,nshockpoints(nshmax),nphpoin
      integer nspecpoints,shinspps(2,5,*)

!     .. array arguments ..
      double precision xysh(ndim,npshmax,*),
     +                 xy(ndim,*),
     &                 zroesh(ndof,npshmax,*),
     &                 zroeshu(ndof,npshmax,*),
     +                 zroe(ndim,*)

      integer          icelnod(nvt,nelem),
     +                 npoin(0:*)
      double precision vshnor(ndim,npshmax,*)

!     .. array arguments ..
      double precision xybkg(ndim),
     +                 zbkg(ndof)

!     .. character array arguments
      character*1 typesh(*)
      character*5 typespecpoints(*)

!     .. local scalars ..
      integer ipoin,ielem,i,ii,k,n,ifail,ish,clr,bbgn,bend,j,kp1,ibc
      integer ip,ip1,ish1,isppnts
      double precision x0,y0,x1,y1,x2,y2,dum,dum1,dum2
      double precision uv,vv,av,thetav,rov,help,pv,mv,alphav

!     open log file
      open(8,file='log/interp_sp.log')

      write(8,*)'enter in interp_sp'

!     make upstream node coordinates coincide with downstream node coordinates
!     these are the coordinates of the shocked mesh (1)
!     Note: the nof of shock points is that on the shocked mesh
!     not the one on the background mesh since this one might have
!     been updated in the shock redistribution routine called by shockmov
      do isppnts=1,nspecpoints

       if(typespecpoints(isppnts).eq.'SP')then

          ish1 = shinspps(1,1,isppnts)
          i    = shinspps(2,1,isppnts)-1
          ip=1+i*(nshockpoints(ish1)-1)
          ip1=2+i*(nshockpoints(ish1)-3)

          xybkg(1)=xysh(1,ip,ish1)
          xybkg(2)=xysh(2,ip,ish1)
          write(8,*)'vertecx coordinates to find ',xybkg(1),xybkg(2)
     +             ,ip,ip1,ish1,nshockpoints(ish1)

          ifail=0
          call finder(icelnod,nelem,xy,ndim,zroe,ndof,xybkg,
     &                  zbkg,ielem,ifail)
          if(ifail.ne.0)then
               write(8,*)'cell not found '
               stop
          endif
          write(8,*)'found in cell ',ielem,ifail

          zroesh(1,ip,ish1)=zbkg(1)
          zroesh(2,ip,ish1)=zbkg(2)
          zroesh(3,ip,ish1)=zbkg(3)
          zroesh(4,ip,ish1)=zbkg(4)

          zroeshu(1,ip,ish1)=zbkg(1)
          zroeshu(2,ip,ish1)=zbkg(2)
          zroeshu(3,ip,ish1)=zbkg(3)
          zroeshu(4,ip,ish1)=zbkg(4)

        end if
      end do

      close(8)
      return
      end
