! Compute the mean shock velocity

subroutine wsh_mean(&
&wsh_n,&
&wsh_n1,&
&wshmean)

  use mod_kinds, only: wp, i4
  use mod_constants, only: half, ndim, npshmax, nshmax
  use mod_log, only: log_line
  implicit none(type, external)

!     ..array definition..
  real(wp)&
  &wsh_n(ndim, npshmax, nshmax),&
  &wsh_n1(ndim, npshmax, nshmax),&
  &wshmean(ndim, npshmax, nshmax)

!     ..integer definition..
  integer(i4) i1, i2, i3
  character(len=200) :: log_buf

  open (8, file='log/ush_mean.log')
  open (9, file='log/vsh_mean.log')
!      write(8,*) 'n        ush_n                  ',
!     +           'ush_n+1                ','ush_mean'

!      write(9,*) 'n        vsh_n                  ',
!     +           'vsh_n1                 ', 'vshmean'

! Phase 4.5 (ROADMAP.md #16): this loop is nshmax*npshmax-bounded, the
! same shape as the already-OMP'd kernels -- writes go through mod_log
! so it stays safe if this loop is ever parallelized.
  do i1 = 1, nshmax
    do i2 = 1, npshmax
      do i3 = 1, 2
        wshmean(i3, i2, i1) = (wsh_n(i3, i2, i1) + wsh_n1(i3, i2, i1))*half
      end do
      write (log_buf, *) wsh_n(1, i2, i1), wsh_n1(1, i2, i1), wshmean(1, i2, i1)
      call log_line(8, log_buf)
      write (log_buf, *) wsh_n(2, i2, i1), wsh_n1(2, i2, i1), wshmean(2, i2, i1)
      call log_line(9, log_buf)
    end do
  end do

  close (8)
  close (9)

  return
end subroutine wsh_mean
