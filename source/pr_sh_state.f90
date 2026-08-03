! Read the file containing the information concerning the shocks and discontinuities

subroutine pr_sh_state(&
&zroesh,&
&nshocks,&
&nshockpoints,&
&nshockedges)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndof, npshmax
  implicit none(type, external)
  include 'shock.com'

!     .. scalar arguments ..
  integer(i4) nshocks, nspecpoints, nshockedges(*), nshockpoints(*),&
  &isppnts, idummy, nshe

!     .. array arguments ..
  real(wp) zroesh(ndof, npshmax, *)

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
  integer(i4) i, k, ish

!     open log file
  open (8, file='log/re_sdw_info.log')

  do ish = 1, nshocks
    do k = 1, nshockpoints(ish)
      write (*, *) k, zroesh(1, k, ish), zroesh(2, k, ish)
    end do
  end do
  write (*, *)
  pause
  continue

  close (8)
  return
end subroutine pr_sh_state
