! Read the file containing the information concerning the shocks and discontinuities

      subroutine pr_sh_state(
     +                       zroesh,
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

      do ish=1,nshocks
        do  k=1, nshockpoints(ish)
          write(*,*)k,zroesh(1,k,ish),zroesh(2,k,ish)
        enddo
      enddo
      write(*,*)
      pause
      continue

      close(8)
      return
      end
