! Filter the position of the points of the discontinuities

      subroutine fltr_dls(
     +           xysh,
     +           zroesh,
     +           wsh,
     +           iter,
     +           nshocks,
     +           nshockpoints,
     +           nshockedges,
     +           typesh)

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

      double precision fltrsh(nshmax)

!     .. local scalar  ..
      double precision dum,p,a,help,ro,u,v

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
      double precision dt,wshmod
      double precision xyshnew(ndim, npshmax, npshmax)

      integer i,im,iv,ish,k

!     open log file
      open(8,file='log/fltr_dls.log')

      write(8,*)'Warning: active filter!'

      fltrsh(1)=0.0d+0
      fltrsh(2)=0.0d+0
      fltrsh(3)=0.0d+0
      fltrsh(4)=0.0d+0
      fltrsh(5)=0.0d+0
      fltrsh(6)=0.0d+0
      fltrsh(7)=0.0d+0
      fltrsh(8)=0.0d+0

      do ish=1,nshocks
        if(fltrsh(ish).eq.0.0d+0)then
          write(8,*)'shock/disc. n.', ish,'filter inactive'
        else
          write(8,*)'shock/disc. n.', ish,'filter active -->',fltrsh(5)
        endif

        do iv = 2, nshockpoints(ish)-1
          do k=1,2
            xyshnew(k,iv,ish)=(1.0d+0-fltrsh(ish))*xysh(k,iv,ish)+
     +                                fltrsh(ish) *(xysh(k,iv+1,ish)+
     +                                              xysh(k,iv-1,ish))*.5
          end do
        end do
        do iv = 2, nshockpoints(ish)-1
          do k=1,2
            xysh(k,iv,ish)=xyshnew(k,iv,ish)
!           if(iter.le.10) xysh(k,iv,ish)=xyshnew(k,iv,ish)
          end do
        end do

      end do

      close(8)
      return
      end
