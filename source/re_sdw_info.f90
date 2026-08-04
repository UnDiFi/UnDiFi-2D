! Read the file containing the information concerning shocks and discontinuities

subroutine re_sdw_info(xysh,&
&zroeshu,&
&zroeshd,&
&zroeshuold,&
&zroeshdold,&
&nodcodsh,&
&coor,&      !vale
&ibndfac,&   !vale
&nbfac,&     !vale
&npoin,&     !vale
&nshocks,&
&nshockpoints,&
&nshockedges,&
&typesh,&
&nspecpoints,&
&typespecpoints,&
&shinspps,&
&ispclr)

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax, nshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_read
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) nshocks, nspecpoints, nshockedges(*), nshockpoints(*),&
  &isppnts, idummy
  character*1 typesh(*)
  character*5 typespecpoints(*)

!     .. array arguments ..
  real(wp) xysh(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *),&
  &zroeshuold(ndof, npshmax, *),&
  &zroeshdold(ndof, npshmax, *)
! vale
  integer(i4) npoin, nbfac, ibndfac(3, *)
  real(wp) coor(ndim, *)
! vale
  integer(i4) nodcodsh(npshmax, *),&
  &shinspps(2, 5, *),&
  &ispclr(5, *)

!     .. local scalars ..
  integer(i4) i, k, ish
  class(special_point_t), allocatable :: sp

!     open log file
  open (8, file='log/re_sdw_info.log')

!     initialize nodcodsh which is part of nodcod
!     if the code -99 is used this means no shock point
!     if the code 10 is used this means shock point
  do ish = 1, nshmax
    do k = 1, npshmax
      nodcodsh(k, ish) = -99
    end do
  end do

!     open file outnew containing the information concerning shocks and discontinuities
  open (12, file='sh00.dat', status='old')
  write (8, *) ' open file sh00.dat'
  read (12, *) nshocks
  write (8, *) ' found n.', nshocks, 'shocks/discontinuities'
  do ish = 1, nshocks
    write (8, *) ' shock/discontinuity n.', ish
    read (12, *) nshockpoints(ish), typesh(ish)
    write (8, *) ' kind of discontinuity:', typesh(ish)
    write (8, *) ' n. of points', nshockpoints(ish)
    write (8, *) 'shock-point coordinates, downstream, upstream states'

    nshockedges(ish) = nshockpoints(ish) - 1
    do k = 1, nshockpoints(ish)
      read (12, *) (xysh(i, k, ish), i=1, ndim),&
      &(zroeshd(i, k, ish), i=1, ndof),&
      &(zroeshu(i, k, ish), i=1, ndof)
      write (8, *) (xysh(i, k, ish), i=1, ndim),&
      &(zroeshd(i, k, ish), i=1, ndof),&
      &(zroeshu(i, k, ish), i=1, ndof)

      nodcodsh(k, ish) = 10
    end do

    do k = 1, nshockpoints(ish)
      do i = 1, ndof
        zroeshdold(i, k, ish) = zroeshd(i, k, ish)
        zroeshuold(i, k, ish) = zroeshu(i, k, ish)
      end do
    end do

  end do

  idummy = 0

  read (12, *) nspecpoints
  write (8, *) nspecpoints
  do isppnts = 1, nspecpoints
    read (12, *) typespecpoints(isppnts)
    write (8, *) typespecpoints(isppnts)

    call make_special_point(typespecpoints(isppnts), sp)
    idummy = idummy + sp%nshe
    call sp_read(12, 8, sp, shinspps(:, :, isppnts), ispclr(:, isppnts))
  end do

!     check condition on special points
  if (idummy .ne. 2*nshocks) then
    write (8, *) 'wrong n. of conditions on special points '
    write (8, *) 'n. imposed condition on special points:', nspecpoints
    write (8, *) 'n. required conditions:', 2*nshocks
    write (8, *) 'n. imposed conditions:', idummy
    write (*, *) 'wrong n. of conditions on special points '
    stop
  end if

  close (8)

  return
end subroutine re_sdw_info
