! Filter the speed of points of the discontinuity lines (optional) and move them

      subroutine mv_dps(xysh,
     +                  zroesh,
     +                  wsh,
     +                  iter,
     +                  nshocks,
     +                  nshockpoints,
     +                  nshockedges,
     +                  typesh)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer iter,nshocks,nshockpoints(nshmax),nshockedges(nshmax)
      character*1 typesh(*)

!     .. array arguments ..
      double precision
     +                 xysh   (ndim, npshmax, *),
     +                 zroesh (ndof, npshmax, *),
     +                 wsh    (ndim, npshmax, *)

!     .. local scalar  ..
      double precision dum,p,a,help,ro,u,v

!     .. local scalars ..
      double precision dt,wshmod
      double precision xyshnew(ndim, npshmax, npshmax),
     +                 wshnew (ndim, npshmax, npshmax)

      integer i,im,iv,ish,k

!     open log file
      open(8,file='log/mv_dps.log')

!     compute the dt max
!     dxcell=0.01*0.4

      dt=1.0d+39
      do ish = 1,nshocks
        write(8,*)'shock n.', ish

        do iv = 1,nshockpoints(ish)

!         u    = zroesh(3,iv,ish)/zroesh(1,iv,ish)
!         v    = zroesh(4,iv,ish)/zroesh(1,iv,ish)
          ro   = zroesh(1,iv,ish)*zroesh(1,iv,ish)
          help = zroesh(3,iv,ish)**2+zroesh(4,iv,ish)**2
          p    = gm1/ga*(zroesh(1,iv,ish)*zroesh(2,iv,ish)-0.5d0*help)
          a    = sqrt(ga*p/ro)
          dum  = 0.0d+0
          do k = 1,2
            dum=dum+wsh(k,iv,ish)**2
          end do
          wshmod=sqrt(dum)
          write(8,*) 'shock pnt n.', iv,' speed:', wshmod
!         dum=dxcell*sndmin/a
!         cfl=0.010
!         if(iter.lt.801)cfl=0.010
          dum=shrelax*dxcell*sndmin/(a+wshmod)
          if(dt.gt.dum) dt=dum
        end do
      end do
      write(8,*) 'dt max:', dt

!     filter shock speeds
      do ish = 1,nshocks
        write(8,*)'shock/disc. n.', ish
        do iv = 2, nshockpoints(ish)-1, 1
          do k = 1,2
            wshnew(k,iv,ish)=(1.0-flt_dspeed*0.5) * wsh(k,iv,ish) +
     +             flt_dspeed*0.25*(wsh(k,iv+1,ish)+wsh(k,iv-1,ish))
          end do
        end do

        do iv = 2, nshockpoints(ish)-1,1
          do k = 1,2
            wsh(k,iv,ish) = wshnew(k,iv,ish)
          end do
        end do
      end do

!     if(iter.le.10)dt=0.0d0

      do ish=1,nshocks
        write(8,*)'shock/disc. n.', ish
        do iv = 1,nshockpoints(ish)
          do k = 1,2
            xysh(k,iv,ish)=xysh(k,iv,ish)+wsh(k,iv,ish)*dt
          end do
          write(8,*) xysh(1,iv,ish), xysh(2,iv,ish)
        end do
      end do

      close(8)

      return
      end
