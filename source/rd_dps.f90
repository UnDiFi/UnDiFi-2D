! Remove discontinuity points where they are too close and insert new
! discontinuity points where they are too far.

subroutine rd_dps(&
&xysh,&
&zroeshu,&
&zroeshd,&
&nshocks,&
&nshockpoints,&
&nshockedges)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use mod_log, only: log_line
  implicit none(type, external)
  include 'paramt.h'

!     .. scalar arguments ..
  integer(i4) iter, nshocks, nshockpoints(nshmax), nshockedges(nshmax)

!     .. array arguments ..
  real(wp)&
  &xysh(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *)

!     .. local scalar  ..
  real(wp) dum, length_rel_min, length_rel_max

!     .. array arguments ..
!     character*(*) fname
  real(wp) sh_edge_lgth(npshmax)

!     .. local scalars ..
  real(wp) dt
  integer(i4) i, im, iv, ish, k, ile_min, ile_max, np_c, np_i
  character(len=200) :: log_buf

!     open log file
  open (8, file='log/rd_dps.log')

  do ish = 1, nshocks
    ile_min = 0
    length_rel_min = 1.0
    ile_max = 0
    length_rel_max = 1.0

!       compute length of shock edge
    do iv = 1, nshockedges(ish)
      sh_edge_lgth(iv) = 0.0d0
      do k = 1, ndim
        sh_edge_lgth(iv) = sh_edge_lgth(iv) +&
        &(xysh(k, iv, ish) - xysh(k, iv + 1, ish))**2
      end do
      sh_edge_lgth(iv) = sqrt(sh_edge_lgth(iv))
      dum = sh_edge_lgth(iv)/dxcell
      if (dum .lt. 0.5) then
        if (length_rel_min .gt. dum) then
          length_rel_min = dum
          ile_min = iv
        end if
      end if
      if (dum .gt. 1.5) then
        if (length_rel_max .lt. dum) then
          length_rel_max = dum
          ile_max = iv
        end if
      end if
    end do

    if (nshockpoints(ish) .le. 3) then
      ile_min = 0.0d0
      ile_max = 0.0d0
    end if

! Phase 4.5 (ROADMAP.md #16): all four write loops below are
! nshockpoints(ish)-bounded, the same shape as the already-OMP'd
! kernels -- writes go through mod_log so they stay safe if this outer
! do-ish loop is ever parallelized across shocks.
    if (ile_min .ne. 0) then
      call log_line(8, 'before')
      do iv = 1, nshockpoints(ish)
        write (log_buf, *) iv, zroeshd(1, iv, ish), zroeshd(2, iv, ish)
        call log_line(8, log_buf)
      end do

      np_c = ile_min
      if (sh_edge_lgth(ile_min - 1) .gt. sh_edge_lgth(ile_min + 1))&
      &np_c = ile_min + 1
      if (ile_min .eq. 1) np_c = 2
      if (ile_min .eq. nshockedges(ish)) np_c = ile_min

      do iv = np_c + 1, nshockpoints(ish)
        do k = 1, ndim
          xysh(k, iv - 1, ish) = xysh(k, iv, ish)
        end do
        do k = 1, ndof
          zroeshu(k, iv - 1, ish) = zroeshu(k, iv, ish)
          zroeshd(k, iv - 1, ish) = zroeshd(k, iv, ish)
        end do
      end do
      do k = 1, ndim
        xysh(k, iv, ish) = 0.0d+0
      end do
      do k = 1, ndof
        zroeshu(k, iv, ish) = 0.0d+0
        zroeshd(k, iv, ish) = 0.0d+0
      end do

      nshockpoints(ish) = nshockpoints(ish) - 1
      nshockedges(ish) = nshockedges(ish) - 1

      call log_line(8, 'after')
      do iv = 1, nshockpoints(ish)
        write (log_buf, *) iv, zroeshd(1, iv, ish), zroeshd(2, iv, ish)
        call log_line(8, log_buf)
      end do

    end if

    if (ile_max .ne. 0) then
      call log_line(8, 'before')
      do iv = 1, nshockpoints(ish)
        write (log_buf, *) iv, zroeshd(1, iv, ish), zroeshd(2, iv, ish)
        call log_line(8, log_buf)
      end do

      np_i = ile_max

      do iv = nshockpoints(ish), np_i + 1, -1
        do k = 1, ndim
          xysh(k, iv + 1, ish) = xysh(k, iv, ish)
        end do
        do k = 1, ndof
          zroeshu(k, iv + 1, ish) = zroeshu(k, iv, ish)
          zroeshd(k, iv + 1, ish) = zroeshd(k, iv, ish)
        end do
      end do
      do k = 1, ndim
        xysh(k, np_i + 1, ish) = 0.5*(xysh(k, np_i, ish) + xysh(k, np_i + 2, ish))
      end do
      do k = 1, ndof
        zroeshu(k, np_i + 1, ish) = 0.5*(zroeshu(k, np_i, ish) +&
        &zroeshu(k, np_i + 2, ish))
        zroeshd(k, np_i + 1, ish) = 0.5*(zroeshd(k, np_i, ish) +&
        &zroeshd(k, np_i + 2, ish))
      end do

      nshockpoints(ish) = nshockpoints(ish) + 1
      nshockedges(ish) = nshockedges(ish) + 1

      call log_line(8, 'after')
      do iv = 1, nshockpoints(ish)
        write (log_buf, *) iv, zroeshd(1, iv, ish), zroeshd(2, iv, ish)
        call log_line(8, log_buf)
      end do

    end if

  end do

  close (8)
  return
end subroutine rd_dps
