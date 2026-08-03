
subroutine fx_sh_state(zroesh,&
&nshocks,&
&nshockpoints,&
&nshockedges)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndof, npshmax
  implicit none(type, external)

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
      zroesh(1, k, ish) = 1.0d+0
      zroesh(2, k, ish) = 0.50625d+0
      zroesh(3, k, ish) = -1.0d+0
      zroesh(4, k, ish) = 0.0d+0
    end do
  end do

  close (8)
  return
end subroutine fx_sh_state
