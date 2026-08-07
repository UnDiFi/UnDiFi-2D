! Compute all discontinuity points, shocks and contact discontinuities

!    +           xyshu
!    +           xyshd,
subroutine co_state_dps(&
&xysh,&
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
&iter)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  implicit none(type, external)
  external co_dc, co_shock
  include 'paramt.h'

!     .. scalar arguments ..
  integer(i4) iter, nshocks, nshockpoints(nshmax), nshockedges(nshmax)
  character*1 typesh(*)

!     .. array arguments ..
!    +                 xyshu      (ndim, npshmax, *),
!    +                 xyshd      (ndim, npshmax, *),
  real(wp)&
  &xysh(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *),&
  &zroeshuold(ndof, npshmax, *),&
  &zroeshdold(ndof, npshmax, *),&
  &vshnor(ndim, npshmax, *),&
  &wsh(ndim, npshmax, *)

!     .. local scalars ..
  real(wp) dx, dy, kine, ws, hh, mmn
  real(wp) help, x1(ndof), x2(ndof)
  real(wp) r2(npshmax, nshmax)
  integer(i4) i, im, iv, ish, k, totnshockpoints

! z1m/z1v/z2m/z2v/z3m/z3v/z4m/z4v were mod_freestream module scratch in
! the original -- write-only-then-immediately-consumed within one loop
! iteration, never read by anything else (confirmed: no other file
! `use`s them). Made local, same fix already applied to the identical
! pattern in mod_special_point.f90's tp_solve_state; a per-iteration
! write to a module global is a data race under `!$omp parallel do`
! (Phase 4.2 target).
  real(wp) z1m, z2m, z3m, z4m
  real(wp) z1v, z2v, z3v, z4v

!     open log file
  open (8, file='log/co_dps_state.log')

  totnshockpoints = 0
  do ish = 1, nshocks
    write (8, *) 'shock/disc. n.', ish

    write (20, *) 'Zone'

    do iv = 1, nshockpoints(ish)
      totnshockpoints = totnshockpoints + 1

      i = iv
      dx = vshnor(1, i, ish)
      dy = vshnor(2, i, ish)

!        load downstream state (for the shock case)
      x1(4) = (-zroeshd(3, iv, ish)*dy + zroeshd(4, iv, ish)*dx)
      x1(4) = x1(4)/zroeshd(1, iv, ish)                                ! tangential
      x1(3) = (zroeshd(3, iv, ish)*dx + zroeshd(4, iv, ish)*dy)
      x1(3) = x1(3)/zroeshd(1, iv, ish)                                ! normal
      x1(1) = zroeshd(1, iv, ish)*zroeshd(1, iv, ish)                    ! density
      help = zroeshd(3, iv, ish)**2 + zroeshd(4, iv, ish)**2
      x1(2) = gm1/ga*(zroeshd(1, iv, ish)*zroeshd(2, iv, ish) - 0.5d0*help)! pressure
      r2(iv, ish) = sqrt(ga*x1(2)/x1(1)) + 0.5*gm1*x1(3)

      write (8, *) 'zd(1)', zroeshd(1, iv, ish)
      write (8, *) 'zd(2)', zroeshd(2, iv, ish)
      write (8, *) 'zd(3)', zroeshd(3, iv, ish)
      write (8, *) 'zd(4)', zroeshd(4, iv, ish)
      write (8, *) 'r2d  ', r2(iv, ish)

!        im is the node corresponding to the shock point iv
      im = iv

!        load upstream state (for the shock case)
      x2(4) = (-zroeshu(3, im, ish)*dy + zroeshu(4, im, ish)*dx)
      x2(4) = x2(4)/zroeshu(1, im, ish)
      x2(3) = (zroeshu(3, im, ish)*dx + zroeshu(4, im, ish)*dy)
      x2(3) = x2(3)/zroeshu(1, im, ish)
      x2(1) = zroeshu(1, im, ish)*zroeshu(1, im, ish)
      help = zroeshu(3, im, ish)**2 + zroeshu(4, im, ish)**2
      x2(2) = gm1/ga*(zroeshu(1, im, ish)*zroeshu(2, im, ish) - 0.5d0*help)

      write (8, *) 'zu(1)', zroeshu(1, im, ish)
      write (8, *) 'zu(2)', zroeshu(2, im, ish)
      write (8, *) 'zu(3)', zroeshu(3, im, ish)
      write (8, *) 'zu(4)', zroeshu(4, im, ish)

!        initialize discontinuity velocity
      ws = 0.0d0

!        compute shock or discontinuity points
      if (typesh(ish) .eq. 'S') call co_shock(x1, x2, ws, r2(iv, ish))

      if (typesh(ish) .eq. 'D') call co_dc(x1, x2, ws)

      mmn = abs(ws - x2(3))/sqrt(ga*x2(2)/x2(1))

      write (20, fmt=400) (xysh(k, iv, ish), k=1, ndim),&
      &(x1(k), k=1, 4), (x2(k), k=1, 4), mmn, ws, dx, dy, i

!        wsavg = wsavg + abs(ws)
!        wsmax = max(abs(ws),wsmax)
400   format(1x, 14(d12.4, 1x), i3)
!400     format(a1,1x,i3,1x,11(d12.4,1x))

!        impose equality of tangential components for the shock
      if (typesh(ish) .eq. 'S') x1(4) = x2(4)

!        assign computed values

!        compute downstream primitive variables (for the shock case)
      z1v = sqrt(x1(1))
      kine = 0.5d0*(x1(3)*x1(3) + x1(4)*x1(4))
      hh = ga/gm1*x1(2)/x1(1) + kine
      z2v = z1v*hh
      z3v = z1v*(x1(3)*dx - x1(4)*dy)
      z4v = z1v*(x1(3)*dy + x1(4)*dx)

!        compute upstream primitive variables (for the shock case)
      z1m = sqrt(x2(1))
      kine = 0.5d0*(x2(3)*x2(3) + x2(4)*x2(4))
      hh = ga/gm1*x2(2)/x2(1) + kine
      z2m = z1m*hh
      z3m = z1m*(x2(3)*dx - x2(4)*dy)
      z4m = z1m*(x2(3)*dy + x2(4)*dx)

!        save old downstream state and assign shock downstream recomputed state
      do k = 1, ndof
        zroeshdold(k, iv, ish) = zroeshd(k, iv, ish)
      end do

      zroeshd(1, iv, ish) = z1v
      zroeshd(2, iv, ish) = z2v
      zroeshd(3, iv, ish) = z3v
      zroeshd(4, iv, ish) = z4v

      write (8, *) 'zd(1)_new', zroeshd(1, iv, ish)
      write (8, *) 'zd(2)_new', zroeshd(2, iv, ish)
      write (8, *) 'zd(3)_new', zroeshd(3, iv, ish)
      write (8, *) 'zd(4)_new', zroeshd(4, iv, ish)

!        save old upstream state and assign shock upstream recomputed state
!        do k = 1,ndof
!          zroeshuold(k,iv,ish)=zroeshu(k,iv,ish)
!        end do

      zroeshu(1, im, ish) = z1m
      zroeshu(2, im, ish) = z2m
      zroeshu(3, im, ish) = z3m
      zroeshu(4, im, ish) = z4m

      write (8, *) 'zu(1)_new', zroeshu(1, iv, ish)
      write (8, *) 'zu(2)_new', zroeshu(2, iv, ish)
      write (8, *) 'zu(3)_new', zroeshu(3, iv, ish)
      write (8, *) 'zu(4)_new', zroeshu(4, iv, ish)

!        assign discontinuity velocity
      wsh(1, iv, ish) = ws*dx
      wsh(2, iv, ish) = ws*dy

      write (8, *) 's/d pnt n.', iv, 'speed:', wsh(1, iv, ish), wsh(2, iv, ish)

    end do
  end do
  close (8)

  return
end subroutine co_state_dps
