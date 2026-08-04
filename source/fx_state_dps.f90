! Fix the state of the special points (ending points) of discontinuities

subroutine fx_state_dps(&
&xysh,&
&xyshu,&
&xyshd,&
&zroeshu,&
&zroeshd,&
&zroeshuold,&
&zroeshdold,&
&vshnor,&
&wsh,&
&nshocks,&
&nshockpoints,&
&nshockedges,&
&typesh,&
&iter,&
&nspecpoints,&
&typespecpoints,&
&shinspps,&
&ispclr,&
&ia,&
&ja,&
&iclr,&
&nclr,&
&corg)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax, nshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_unpack
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) iter, nshocks, nshockpoints(nshmax), nshockedges(nshmax),&
  &nspecpoints, shinspps(2, 5, *), ispclr(*)
  character*1 typesh(*)
  character*5 typespecpoints(*)

!     .. array arguments ..
  real(wp)&
  &xysh(ndim, npshmax, *),&
  &xyshu(ndim, npshmax, *),&
  &xyshd(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *),&
  &zroeshuold(ndof, npshmax, *),&
  &zroeshdold(ndof, npshmax, *),&
  &vshnor(ndim, npshmax, *),&
  &wsh(ndim, npshmax, *)
  integer(i4) nclr
  integer(i4) ia(*), ja(*), iclr(nclr)
  real(wp) corg(ndim, *)

!     .. local scalars ..
  real(wp) wws, varz(ndof), avarz(ndof)
  integer(i4) i, iv, ish, ip
  integer(i4) isppnts
  class(special_point_t), allocatable :: sp

!     open log file
  open (8, file='log/fx_state_sps.log')

  do isppnts = 1, nspecpoints
    call make_special_point(typespecpoints(isppnts), sp)
    call sp_unpack(sp, shinspps(:, :, isppnts))
    call sp%solve_state(xysh(:, :, 1:nshmax), zroeshu(:, :, 1:nshmax),&
    &zroeshd(:, :, 1:nshmax), zroeshuold(:, :, 1:nshmax),&
    &zroeshdold(:, :, 1:nshmax), vshnor(:, :, 1:nshmax),&
    &wsh(:, :, 1:nshmax), nshockpoints, ia, ja, iclr, nclr, corg,&
    &ispclr(isppnts))
  end do

!     calculate variation of the shock points state
  wws = 0.0
  do iv = 1, ndof
    varz(iv) = 0.0
    avarz(iv) = 0.0
  end do
  do ish = 1, nshocks
    write (8, *) 'variations of z downstream of the shock'
    write (8, *) 'shock n.', ish
    do i = 1, nshockpoints(ish)
      do iv = 1, ndof
        varz(iv) = zroeshd(iv, i, ish) - zroeshdold(iv, i, ish)
        avarz(iv) = avarz(iv) + sqrt(varz(iv)**2)
      end do
      wws = wws + sqrt(wsh(1, i, ish)**2 + wsh(2, i, ish)**2)

      write (8, '(1x,i3,5(1x,f15.7))') i, (varz(iv), iv=1, ndof)
    end do
  end do

  write (8, '(1x,a19,5(1x,f15.7))') 'average values var z:',&
  &(avarz(iv), iv=1, ndof)
  write (87, '(1x,i5,5(1x,f15.7))') iter, (avarz(iv), iv=1, ndof), wws

!      save the old upstream state and assign
!      the recomputed shock upstream state
  do ish = 1, nshocks
    do ip = 1, nshockpoints(ish)
      do iv = 1, ndof
        zroeshdold(iv, ip, ish) = zroeshd(iv, ip, ish)
        zroeshuold(iv, ip, ish) = zroeshu(iv, ip, ish)
      end do
    end do
  end do

  close (8)

  return
end subroutine fx_state_dps
