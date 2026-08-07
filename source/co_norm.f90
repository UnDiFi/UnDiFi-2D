! Compute the normal unit vectors to the shocks and discontinuities in shock/discontinuity points

subroutine co_norm(xysh,&
&zroeshu,& !upstream
&zroesh,&  !downstream
&vshnor,&
&nshocks,&
&nshockpoints,&
&typesh,&
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
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_unpack
  implicit none(type, external)
  include 'paramt.h'

  integer(i4) nshocks, nshockpoints(*), nspecpoints, shinspps(2, 5, *)
  integer(i4) ispclr(5, *)
  integer(i4) i, j, j2, k, kp1, depip1, depim1, ii

  integer(i4) nclr
  integer(i4) ia(*), ja(*), iclr(nclr)
  real(wp) corg(ndim, *)

  real(wp) xysh(ndim, npshmax, *)
  real(wp) vshnor(ndim, npshmax, *)
  real(wp) zroesh(ndof, npshmax, *)

! vale
  real(wp) zroeshu(ndof, npshmax, *)
! vale

  real(wp) xi, yi, xj, yj, ush, vsh, tau, dum
  real(wp) help
  real(wp) uj, vj, roj, aj, pj, mj, alphaj, thetaj
  real(wp) uv, vv, rov, av, pv, mv, alphav, thetav
  real(wp) tauxim1, tauyim1, tauxip1, tauyip1, taux, tauy
  real(wp) xj2, yj2, tauxip2, tauyip2, tauxim2, tauyim2
  real(wp) tauxjp2, tauyjp2, tauxjm2, tauyjm2
  real(wp) lp12, lm12, ui, vi
  real(wp) lp1, lp2, lm1, lm2, lp22, lm22
  integer(i4) shp_dpndnc, dcp_dpndnc, ish, isppnts
  external shp_dpndnc, dcp_dpndnc
  character*1 typesh(*)
  character*5 typespecpoints(*)

  class(special_point_t), allocatable :: sp

!     input:
!     -----
!     xysh   x,y coords of the shock points
!     zroesh status variables
!     nshockpoints number of shock points
!
!     output:
!     ------
!     vshnor x,y components of the unit normal to the shock

!     open log file
  open (8, file='log/co_norm.log')

  write (8, *) 'n. of shocks:', nshocks
  do ish = 1, nshocks
    write (8, *) 'n. of shock points:', nshockpoints(ish)
    write (8, *) 'shock-point coordinates, downtream'
    do i = 1, nshockpoints(ish)
      write (8, *) i, xysh(1, i, ish), xysh(2, i, ish),&
      &zroesh(1, i, ish), zroesh(2, i, ish),&
      &zroesh(3, i, ish), zroesh(4, i, ish)
    end do
  end do

!     normals computation for each shock
  do ish = 1, nshocks
    do i = 1, nshockpoints(ish)

      ush = 0.d0
      vsh = 0.d0

      xi = xysh(1, i, ish)
      yi = xysh(2, i, ish)

!       upward points
      j = i + 1

!       recover coordinates
      xj = xysh(1, j, ish)
      yj = xysh(2, j, ish)

      j2 = i + 2
      xj2 = xysh(1, j2, ish)
      yj2 = xysh(2, j2, ish)

!       tangent vectors computation
      tauxip1 = xj - xi
      tauyip1 = yj - yi

      tauxip2 = xj2 - xi
      tauyip2 = yj2 - yi
      tauxjp2 = xj2 - xj
      tauyjp2 = yj2 - yj

!       recover state
      uj = zroesh(3, j, ish)/zroesh(1, j, ish)
      vj = zroesh(4, j, ish)/zroesh(1, j, ish)
      roj = zroesh(1, j, ish)*zroesh(1, j, ish)
      help = zroesh(3, j, ish)**2 + zroesh(4, j, ish)**2
      pj = gm1/ga*(zroesh(1, j, ish)*zroesh(2, j, ish) - 0.5d0*help)
      aj = sqrt(ga*pj/roj)

!       dependency evaluation
!       if it is the last points, then compute backward
      if (i .ne. 1 .and. i .ne. nshockpoints(ish)) then
        if (typesh(ish) .eq. 'S')&
        &depip1 = shp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
        if (typesh(ish) .eq. 'D')&
        &depip1 = dcp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
      elseif (i .eq. nshockpoints(ish)) then
        depip1 = 0
        depim1 = 1
!        goto 999
      end if

!       backward points
      j = i - 1

!       recover coordinates
      xj = xysh(1, j, ish)
      yj = xysh(2, j, ish)

      j2 = i - 2
      xj2 = xysh(1, j2, ish)
      yj2 = xysh(2, j2, ish)

!       tangent vectors computation
      tauxim1 = xi - xj
      tauyim1 = yi - yj

      tauxim2 = xi - xj2
      tauyim2 = yi - yj2
      tauxjm2 = xj - xj2
      tauyjm2 = yj - yj2

!       recover state
      uj = zroesh(3, j, ish)/zroesh(1, j, ish)
      vj = zroesh(4, j, ish)/zroesh(1, j, ish)
      roj = zroesh(1, j, ish)*zroesh(1, j, ish)
      help = zroesh(3, j, ish)**2 + zroesh(4, j, ish)**2
      pj = gm1/ga*(zroesh(1, j, ish)*zroesh(2, j, ish) - 0.5d0*help)
      aj = sqrt(ga*pj/roj)

!       dependency evaluation
!       if it is the first point, then compute upward
      if (i .ne. 1 .and. i .ne. nshockpoints(ish)) then
        if (typesh(ish) .eq. 'S')&
        &depim1 = shp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
        if (typesh(ish) .eq. 'D')&
        &depim1 = dcp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
      elseif (I .eq. 1) then
        depim1 = 0
        depip1 = 1
!        goto 999
      end if

999   continue

!       tangent versors computation
      lp12 = tauxip1**2 + tauyip1**2
      lm12 = tauxim1**2 + tauyim1**2
      lp1 = sqrt(lp12)
      lm1 = sqrt(lm12)

      lm22 = tauxjm2**2 + tauyjm2**2
      lp22 = tauxjp2**2 + tauyjp2**2

      lp2 = sqrt(lp22)
      lm2 = sqrt(lm22)

      if (i .eq. 1) then
        depim1 = 0
        depip1 = 1
        lm12 = 1.0
      end if

      if (i .eq. nshockpoints(ish)) then
        depim1 = 1
        depip1 = 0
        lp12 = 1.0
      end if

!       TODO: verify the case when the point has no dependence
      if (depim1 .eq. 0 .and. depip1 .eq. 0) then
        depim1 = 1
        depip1 = 1
      end if

      taux = (tauxim1*lp12 + tauxip1*lm12)
      tauy = (tauyim1*lp12 + tauyip1*lm12)
      if (depim1*depip1 .eq. 0.0) then
        if (depim1 .eq. 1) then
          taux = tauxim1*(lm1 + lm2)**2 - tauxim2*lm12
          tauy = tauyim1*(lm1 + lm2)**2 - tauyim2*lm12
        end if
        if (depip1 .eq. 1) then
          taux = tauxip1*(lp1 + lp2)**2 - tauxip2*lp12
          tauy = tauyip1*(lp1 + lp2)**2 - tauyip2*lp12
        end if
      end if

      write (8, *) I, 'tau', taux, tauy, depim1, depip1

      tau = sqrt(taux*taux + tauy*tauy)
      taux = taux/tau
      tauy = tauy/tau

!       calculate and assign of the normal versor
!       TODO: check orientation downstream->upstream of the versor
!       if(ish.eq.2) write(*,*)'taux,tauy:',i,taux,tauy,
!    &                          xysh(1,i,ish),xysh(2,i,ish)

      vshnor(1, i, ish) = tauy
      vshnor(2, i, ish) = -taux

    end do

  end do

!     compute and/or assign the normal in the union point of two shocks
!     in this case, the normal of the first point of the second shock
!     coincides ith the normal of the last point of the first shock
!
!          vshnor(1,nshockpoints(1),1)=vshnor(1,1,2)
!          vshnor(2,nshockpoints(1),1)=vshnor(2,1,2)
!
!     check the correct orientation of the normals
!     the orientation is downstream toward upstream
!     verify if u * n <= 0 otherwise change direction of n
!     do ish=1,nshocks
!      if(typesh(ish).eq.'s')then
!      do i=1,nshockpoints(ish)
!       ui=zroesh(3,i,ish)/zroesh(1,i,ish)
!       vi=zroesh(4,i,ish)/zroesh(1,i,ish)
!       dum=ui*vshnor(1,i,ish)+vi*vshnor(2,i,ish)
!       if(dum.gt.0.d0)then
!         vshnor(1,i,ish)= -vshnor(1,i,ish)
!         vshnor(2,i,ish)= -vshnor(2,i,ish)
!       endif
!      enddo
!      endif
!     enddo

! Phase 4.2 note: this count-then-bulk-flip pair is a reduction, not an
! independent per-point loop -- `ii` is accumulated across every `i` before
! the flip loop reads it. A naive `!$omp parallel do` on the count loop
! races on `ii`; it needs `omp reduction(+:ii)` (which also supplies the
! implicit end-of-loop barrier the flip loop's read of `ii` depends on).
! The bulk-flip loop itself has no hazard once `ii` is final -- each
! iteration only touches its own vshnor(:,i,ish) slot.
  do ish = 1, nshocks
    if (typesh(ish) .eq. 'S') then
      ii = 0
      do i = 1, nshockpoints(ish)
        ui = zroesh(3, i, ish)/zroesh(1, i, ish)
        vi = zroesh(4, i, ish)/zroesh(1, i, ish)
        dum = ui*vshnor(1, i, ish) + vi*vshnor(2, i, ish)
        if (dum .gt. 0.d0) then
          ii = ii + 1
        end if
      end do
      if (ii .ge. nshockpoints(ish)/2.) then
        do i = 1, nshockpoints(ish)
! vale   vshnor(1,i,ish)= -vshnor(1,i,ish)
! vale   vshnor(2,i,ish)= -vshnor(2,i,ish)
          vshnor(1, i, ish) = -vshnor(1, i, ish)
          vshnor(2, i, ish) = -vshnor(2, i, ish)
        end do
      end if

    end if
  end do

!     in the case of triple points, it forces the orientation of the
!     normal to the contact discontinuity such that it forms an
!     angle > 90 degrees with the normal to the mach stem
!
! Phase 3.2 increment 8: each isppnts' correct_normal now dispatches
! through the special_point_t registry instead of a 16-way if/elseif.
! make_special_point error-stops on any code outside the 16 handled
! here; WPNRY (which has no branch at all in the original) still hits
! its own fatal stop inside wf_correct_normal, not this error-stop --
! see mod_special_point.f90's wf_correct_normal header for why.
  do isppnts = 1, nspecpoints
    call make_special_point(typespecpoints(isppnts), sp)
    call sp_unpack(sp, shinspps(:, :, isppnts), ispclr(:, isppnts))
    call sp%correct_normal(xysh(:, :, 1:nshmax), zroeshu(:, :, 1:nshmax),&
    &vshnor(:, :, 1:nshmax), nshockpoints(1:nshmax), ia, ja, iclr, nclr, corg)
  end do

!     read normals of the Mach stem and of the reflected shock in the triple point

!     ish=3
!     open(66,file='norm_tp1',err=661)
!     read(66,*)vshnor(1,nshockpoints(ish),ish),
!    .          vshnor(2,nshockpoints(ish),ish)
!
!     ish=2
!     read(66,*)vshnor(1,nshockpoints(ish),ish),
!    .          vshnor(2,nshockpoints(ish),ish)
!661   close(66)

!     impose on the wall the direction of the shock normal parallel to that of the wall

!     vshnor(1,1,3)= 1.0d0
!     vshnor(2,1,3)= 0.0d0

!     write the tecplot file with computed normals

  write (8, *) 'start writing file shocknor.dat'
  open (12, file='shocknor.dat')

  do ish = 1, nshocks
    write (12, *) 'TITLE = Shock normals'
    write (12, *) 'VARIABLES = X Y Z(1) Z(2) NX NY'
    write (12, *) 'ZONE T="sampletext",F=FEPOINT,ET=TRIANGLE, N=',&
    &nShockPoints(ISH), ', E=', nShockPoints(ISH) - 1
    do i = 1, nshockpoints(ish)
      write (12, *) (xysh(k, i, ish), k=1, ndim), 1.d0, 1.d0,&
      &(vshnor(k, i, ish), k=1, ndim)
    end do
    do i = 1, nshockpoints(ish) - 1
      write (12, *) i, i + 1, i
    end do
  end do
  close (12)
  write (8, *) 'end writing file '
  close (8)

!     call exit(-1)

  return
end subroutine co_norm

integer(i4) function shp_dpndnc_old(x, y, ush, vsh, xi, yi, ui, vi, ai)
  use mod_kinds, only: wp, i4
  real(wp) x, y, ush, vsh, xi, yi, ui, vi, ai
  real(wp) taux, tauy, tau, utau

!     determine tangent versor
  taux = xi - x
  tauy = yi - y
  tau = sqrt(taux*taux + tauy*tauy)
  taux = taux/tau
  tauy = tauy/tau

!     determine the velocity component along the verson
  utau = ui*taux + vi*tauy

!     dependence evaluation
  shp_dpndnc = 0
  if ((utau - ai) .lt. 0.0d0) shp_dpndnc = 1

  return
end function shp_dpndnc_old

integer(i4) function shp_dpndnc_old2(x, y, ush, vsh, xi, yi, ui, vi, ai)
  use mod_kinds, only: wp, i4
  real(wp) x, y, ush, vsh, xi, yi, ui, vi, ai, umod
  real(wp) xx, yy, xxi, yyi, dsh1, dsh2, dt
  real(wp) taux1, tauy1, tau1, taux2, tauy2, tau2, dum, dx1, dx2

!     calculate dt for the test
  rl = sqrt((x - xi)**2 + (y - yi)**2)
  umod = sqrt(ui**2 + vi**2)
  dt = 0.5*rl/(ai + umod)

!     compute position of shock point after dt=1
  xx = x + ush*dt
  yy = y + vsh*dt

!     compute position of point i after dt=1
!     amplification factor of the shock velocity
!     to partially extend to the subsonic zone the computation
!     of the upwind formula for the calculation of shock normals
!
!     xxi=xi+3.0*ui*dt
!     yyi=yi+3.0*vi*dt

  xxi = xi + ui*dt
  yyi = yi + vi*dt

!     calculate original distance between shock point and point i
  dsh1 = (x - xi)**2 + (y - yi)**2
  dsh1 = sqrt(dsh1)
  dx1 = dsh1

!     determine versor
  taux1 = xi - x
  tauy1 = yi - y
  tau1 = sqrt(taux1*taux1 + tauy1*tauy1)
  taux1 = taux1/tau1
  tauy1 = tauy1/tau1

!     calculate the perturbation distance of point i and shock point after dt=1
  dsh2 = (xx - xxi)**2 + (yy - yyi)**2
  dx2 = sqrt(dsh2)
!     dsh2=sqrt(dsh2)-ai*dt
  dsh2 = sqrt(dsh2) - ai*dt
!     dsh1=dsh2-umod*dt

!     determine versor
  taux2 = xxi - xx
  tauy2 = yyi - yy
  tau2 = sqrt(taux2*taux2 + tauy2*tauy2)
  taux2 = taux2/tau2
  tauy2 = tauy2/tau2

  dum = taux1*taux2 + tauy1*tauy2

!     distance evaluation
  shp_dpndnc = 0
  if (dsh2 .lt. dsh1 - (umod*dt - abs(dx2 - dx1))) shp_dpndnc = 1
!     write(*,*)'dsh1,dsh2:',dsh1-(umod*dt-abs(dx2-dx1)),dsh2
!     if(dsh2.lt.0.d0)shp_dpndnc=1

  return
end function shp_dpndnc_old2

integer(i4) function shp_dpndnc(x, y, ush, vsh, xi, yi, ui, vi, ai)
  use mod_kinds, only: wp, i4
  real(wp) x, y, ush, vsh, xi, yi, ui, vi, ai
  real(wp) xx, yy, uu, vv, aa, dum, dt1, dt2
  real(wp) taux, tauy, tau, utau

  xx = x - xi
  yy = y - yi

  uu = ui
  vv = vi

  aa = ai*1.00

!     dum=2.0*xx*yy*uu*vv - xx*xx*vv*vv + xx*xx*aa*aa
!     dum=dum             - yy*yy*uu*uu + yy*yy*aa*aa

!     dum=(xx*uu+yy*vv)**2-(xx**2+yy**2)*(uu**2+vv**2-aa**2)
  dum = aa**2*(xx**2 + yy**2) - (uu*yy - vv*xx)**2

  shp_dpndnc = 0
  dt1 = 0.0
  dt2 = 0.0
  if (dum .gt. 0.0) then
    dt1 = ((xx*uu + yy*vv) + sqrt(dum))/(uu**2 + vv**2 - aa**2)
    dt2 = ((xx*uu + yy*vv) - sqrt(dum))/(uu**2 + vv**2 - aa**2)
    if (dt1 .gt. 0.0 .or. dt2 .gt. 0.0) shp_dpndnc = 1
  end if
!     write(*,*) dum,dt1,dt2

  return
end function shp_dpndnc

integer(i4) function dcp_dpndnc(x, y, ush, vsh, xi, yi, ui, vi, ai)
  use mod_kinds, only: wp, i4
  real(wp) x, y, ush, vsh, xi, yi, ui, vi, ai
  real(wp) xx, yy, uu, vv, aa, dum, dt1, dt2
  real(wp) taux, tauy, tau, utau

  xx = x - xi
  yy = y - yi

  uu = ui
  vv = vi

  aa = ai*1.00

  dum = (ui*xx + vi*yy)

  dcp_dpndnc = 0
  if (dum .gt. 0.0) dcp_dpndnc = 1
!     write(*,*) dum,dt1,dt2

  return
end function dcp_dpndnc
