! Write the file containing the information concerning shocks discontinuities

subroutine wrt_sdw_info(xysh,&
&zroeshu,&
&zroeshd,&
&nodcodsh,&
&nshocks,&
&nshockpoints,&
&nshockedges,&
&typesh,&
&nspecpoints,&
&typespecpoints,&
&shinspps,&
&ispclr)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_write
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) nshocks, nspecpoints, nshockedges(*), nshockpoints(*),&
  &isppnts
  character*1 typesh(*)
  character*5 typespecpoints(*)

!     .. array arguments ..
  real(wp) xysh(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *)

  integer(i4) nodcodsh(npshmax, *),&
  &shinspps(2, 5, *),&
  &ispclr(5, *)

!     .. local scalars ..
  integer(i4) i, k, ish
  class(special_point_t), allocatable :: sp

!     open log file
  open (8, file='log/wrt_sdw_info.log')

!     open file sh99 containing the informations concerning shocks and discontinuities
  open (12, file='sh99.dat')
  write (8, *) ' open file sh99.dat'
  write (12, *) nshocks
  write (8, *) ' n.', nshocks, 'shocks/discontinuities'
  do ish = 1, nshocks
    write (8, *) ' shock/discontinuity n.', ish
    write (12, *) nshockpoints(ish), typesh(ish)
    write (8, *) ' kind of discontinuity:', typesh(ish)
    write (8, *) ' n. of points', nshockpoints(ish)

    do k = 1, nshockpoints(ish)
      write (12, *) (xysh(i, k, ish), i=1, ndim),&
      &(zroeshd(i, k, ish), i=1, ndof),&
      &(zroeshu(i, k, ish), i=1, ndof)

    end do

  end do

  write (12, *) nspecpoints
  write (8, *) nspecpoints
  do isppnts = 1, nspecpoints
    write (12, *) typespecpoints(isppnts)
    write (8, *) typespecpoints(isppnts)

    call make_special_point(typespecpoints(isppnts), sp)
    call sp_write(12, 8, sp, shinspps(:, :, isppnts), ispclr(:, isppnts))
  end do

  close (8)

  return
end subroutine wrt_sdw_info
