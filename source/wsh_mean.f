! Compute the mean shock velocity      

      subroutine wsh_mean(
     +           wsh_n,
     +           wsh_n1,
     +           wshmean)

      implicit none
      include 'paramt.h'

!     ..array definition..
      double precision
     +                 wsh_n   (ndim,npshmax,nshmax),
     +                 wsh_n1  (ndim,npshmax,nshmax),
     +                 wshmean (ndim,npshmax,nshmax)

!     ..integer definition..
      integer i1,i2,i3

      open(8,file='log/ush_mean.log')
      open(9,file='log/vsh_mean.log')
!      write(8,*) 'n        ush_n                  ',
!     +           'ush_n+1                ','ush_mean'

!      write(9,*) 'n        vsh_n                  ',
!     +           'vsh_n1                 ', 'vshmean'

      do i1=1,nshmax
        do i2=1,npshmax
          do i3=1,2
             wshmean(i3,i2,i1)=(wsh_n(i3,i2,i1)+wsh_n1(i3,i2,i1))*half
          enddo
         write(8,*) wsh_n(1,i2,i1),wsh_n1(1,i2,i1),wshmean(1,i2,i1)
         write(9,*) wsh_n(2,i2,i1),wsh_n1(2,i2,i1),wshmean(2,i2,i1)
        end do
      end do 

      close(8)
      close(9)

      return
      end
