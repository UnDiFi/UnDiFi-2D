
      subroutine fx_sh_state(zroesh,
     +                       nshocks,
     +                       nshockpoints,
     +                       nshockedges)

      implicit none
      include 'paramt.h'
      include 'shock.com'

!     .. scalar arguments ..
      integer nshocks,nspecpoints,nshockedges(*),nshockpoints(*),
     +        isppnts,idummy,nshe

!     .. array arguments ..
      double precision zroesh(ndof,npshmax,*)

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
      integer i,k,ish

!     open log file
      open(8,file='log/re_sdw_info.log')

      do ish = 1,nshocks
        do k = 1,nshockpoints(ish)
          zroesh(1,k,ish)=1.0d+0
          zroesh(2,k,ish)=0.50625d+0
          zroesh(3,k,ish)=-1.0d+0
          zroesh(4,k,ish)=0.0d+0
        enddo
      enddo

      close(8)
      return
      end
