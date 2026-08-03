! Write the current shock/discontinuity state to a caller-given file name.
!
! This is a read-only instrumentation clone of wrt_sdw_info: it writes the
! exact same sh99.dat-compatible record layout (nshocks; per shock/front
! nshockpoints+type then one row per point of x,y,zroeshd(1:ndof),
! zroeshu(1:ndof); nspecpoints then per special point its tag and
! connections), but to an arbitrary file name instead of the hardcoded
! 'sh99.dat', so it can be called at multiple points within one iteration
! (after fnd_phps, co_norm, co_pnt_dspl, co_state_dps, ...) without
! clobbering the "official" end-of-iteration sh99.dat. Added for R-020
! (per-subroutine golden dumps for step-level Phase-3 port validation) -
! it does not read or mutate any solver state, only dumps it, so it
! cannot change the numerical results of a run.

subroutine dump_sdw_info(dumpfname,&
&xysh,&
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
  implicit none(type, external)

!     .. scalar arguments ..
  character*(*) dumpfname
  integer(i4) nshocks, nspecpoints, nshockedges(*), nshockpoints(*),&
  &isppnts, nshe
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

!     open the dump file (fresh each call: same semantics as sh99.dat)
  open (14, file=dumpfname)
  write (14, *) nshocks
  do ish = 1, nshocks
    write (14, *) nshockpoints(ish), typesh(ish)
    do k = 1, nshockpoints(ish)
      write (14, *) (xysh(i, k, ish), i=1, ndim),&
      &(zroeshd(i, k, ish), i=1, ndof),&
      &(zroeshu(i, k, ish), i=1, ndof)
    end do
  end do

  write (14, *) nspecpoints
  do isppnts = 1, nspecpoints
    write (14, *) typespecpoints(isppnts)
    if (typespecpoints(isppnts) .eq. 'TP') then
      nshe = 4
    elseif (typespecpoints(isppnts) .eq. 'QP') then
      nshe = 5
    elseif (typespecpoints(isppnts) .eq. 'TE') then
      nshe = 3
    elseif (typespecpoints(isppnts) .eq. 'RRX') then
      nshe = 2
    elseif (typespecpoints(isppnts) .eq. 'RR') then
      nshe = 2
    elseif (typespecpoints(isppnts) .eq. 'WPNRX') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'WPNRY') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'FWP') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'IPX') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'IPY') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'OPX') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'OPY') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'EP') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'SP') then
      nshe = 1
    elseif (typespecpoints(isppnts) .eq. 'C') then
      nshe = 2
    elseif (typespecpoints(isppnts) .eq. 'PC') then
      nshe = 2
    else
      write (*, *) 'dump_sdw_info: special point type not implemented!'
      stop
    end if

    if (typespecpoints(isppnts) .eq. 'RR' .or.&
    &typespecpoints(isppnts) .eq. 'FWP' .or.&
    &typespecpoints(isppnts) .eq. 'PC') then
      do k = 1, nshe
        write (14, *) shinspps(1, k, isppnts), shinspps(2, k, isppnts),&
        &ispclr(k, isppnts)
      end do
    else
      do k = 1, nshe
        write (14, *) shinspps(1, k, isppnts), shinspps(2, k, isppnts)
      end do
    end if
  end do

  close (14)

  return
end subroutine dump_sdw_info
