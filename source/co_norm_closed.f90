! Compute the normal unit vectors along a CLOSED discontinuity curve
! (Issue #22/E2: closed discontinuity curves, prerequisite for shock-
! bubble interaction -- a gas bubble's interface has no ends).
!
! A direct cyclic-indexing adaptation of co_norm.f90's own per-point
! tangent/normal computation for a SINGLE curve. Deliberately a
! separate, standalone kernel rather than a branch inside co_norm.f90
! itself -- same reasoning as co_interface.f90/co_shock_2gas.f90 (issue
! #23/#21): co_norm.f90's real call site (co_norm's own caller in
! main.f90) has no closed-curve-aware plumbing (no `closed` flag
! anywhere in discontinuity_t's consumers, no hole-topology support in
! wtri.f90/fnd_phps.f90/co_pnt_dspl.f90 for the "annular hole with an
! interior sub-domain" issue #22 itself calls out as a separate,
! larger task) -- wiring this into the live multi-shock/multi-special-
! point dispatch co_norm.f90 handles would be premature before that
! mesh-topology work exists to actually USE a closed curve for
! anything. This file is the normal-computation kernel only, verified
! standalone (tests/unit/test_co_norm_closed.f90).
!
! Per issue #22's own task wording ("cyclic indexing (modulo npoints,
! no endpoint special points)"), a closed curve has no first/last point
! at all -- every point is topologically identical, so unlike
! co_norm.f90 there is no i==1/i==nshockpoints special-casing anywhere
! in this file. The per-point tangent/dependency logic itself
! (shp_dpndnc/dcp_dpndnc, the weighted central-vs-one-sided tangent
! formula) is otherwise UNCHANGED from co_norm.f90 -- only the neighbor-
! index arithmetic differs.
subroutine co_norm_closed(xysh, zroesh, vshnor, npoints, typesh)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax
  implicit none(type, external)
  include 'paramt.h'

  integer(i4), intent(in) :: npoints
  real(wp), intent(in) :: xysh(ndim, npoints)
  real(wp), intent(in) :: zroesh(ndof, npoints)
  real(wp), intent(out) :: vshnor(ndim, npoints)
  character(len=1), intent(in) :: typesh

  real(wp) xi, yi, xj, yj, ush, vsh, tau
  real(wp) help
  real(wp) uj, vj, roj, aj, pj
  real(wp) tauxim1, tauyim1, tauxip1, tauyip1, taux, tauy
  real(wp) xj2, yj2, tauxip2, tauyip2, tauxim2, tauyim2
  real(wp) tauxjp2, tauyjp2, tauxjm2, tauyjm2
  real(wp) lp12, lm12, lp1, lm1, lp2, lm2, lm22, lp22
  integer(i4) shp_dpndnc, dcp_dpndnc, i, j, j2, depip1, depim1
  external shp_dpndnc, dcp_dpndnc

  do i = 1, npoints

    ush = 0.0_wp
    vsh = 0.0_wp

    xi = xysh(1, i)
    yi = xysh(2, i)

! forward neighbor, cyclic
    j = cyc(i + 1, npoints)
    xj = xysh(1, j)
    yj = xysh(2, j)

    j2 = cyc(i + 2, npoints)
    xj2 = xysh(1, j2)
    yj2 = xysh(2, j2)

    tauxip1 = xj - xi
    tauyip1 = yj - yi
    tauxip2 = xj2 - xi
    tauyip2 = yj2 - yi
    tauxjp2 = xj2 - xj
    tauyjp2 = yj2 - yj

    uj = zroesh(3, j)/zroesh(1, j)
    vj = zroesh(4, j)/zroesh(1, j)
    roj = zroesh(1, j)*zroesh(1, j)
    help = zroesh(3, j)**2 + zroesh(4, j)**2
    pj = gm1/ga*(zroesh(1, j)*zroesh(2, j) - 0.5_wp*help)
    aj = sqrt(ga*pj/roj)

    if (typesh .eq. 'S') then
      depip1 = shp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
    else
      depip1 = dcp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
    end if

! backward neighbor, cyclic
    j = cyc(i - 1, npoints)
    xj = xysh(1, j)
    yj = xysh(2, j)

    j2 = cyc(i - 2, npoints)
    xj2 = xysh(1, j2)
    yj2 = xysh(2, j2)

    tauxim1 = xi - xj
    tauyim1 = yi - yj
    tauxim2 = xi - xj2
    tauyim2 = yi - yj2
    tauxjm2 = xj - xj2
    tauyjm2 = yj - yj2

    uj = zroesh(3, j)/zroesh(1, j)
    vj = zroesh(4, j)/zroesh(1, j)
    roj = zroesh(1, j)*zroesh(1, j)
    help = zroesh(3, j)**2 + zroesh(4, j)**2
    pj = gm1/ga*(zroesh(1, j)*zroesh(2, j) - 0.5_wp*help)
    aj = sqrt(ga*pj/roj)

    if (typesh .eq. 'S') then
      depim1 = shp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
    else
      depim1 = dcp_dpndnc(xi, yi, ush, vsh, xj, yj, uj, vj, aj)
    end if

    lp12 = tauxip1**2 + tauyip1**2
    lm12 = tauxim1**2 + tauyim1**2
    lp1 = sqrt(lp12)
    lm1 = sqrt(lm12)
    lm22 = tauxjm2**2 + tauyjm2**2
    lp22 = tauxjp2**2 + tauyjp2**2
    lp2 = sqrt(lp22)
    lm2 = sqrt(lm22)

    if (depim1 .eq. 0 .and. depip1 .eq. 0) then
      depim1 = 1
      depip1 = 1
    end if

    taux = (tauxim1*lp12 + tauxip1*lm12)
    tauy = (tauyim1*lp12 + tauyip1*lm12)
    if (depim1*depip1 .eq. 0) then
      if (depim1 .eq. 1) then
        taux = tauxim1*(lm1 + lm2)**2 - tauxim2*lm12
        tauy = tauyim1*(lm1 + lm2)**2 - tauyim2*lm12
      end if
      if (depip1 .eq. 1) then
        taux = tauxip1*(lp1 + lp2)**2 - tauxip2*lp12
        tauy = tauyip1*(lp1 + lp2)**2 - tauyip2*lp12
      end if
    end if

    tau = sqrt(taux*taux + tauy*tauy)
    taux = taux/tau
    tauy = tauy/tau

    vshnor(1, i) = tauy
    vshnor(2, i) = -taux

  end do

contains

! Wraps a 1-based index into [1,n] cyclically -- Fortran's MOD gives a
! result in [-(n-1), n-1], so shift to [0,n-1] first (idx-1+n), take
! mod n, shift back to 1-based. The "+n" before the first mod absorbs
! any idx as low as 1-2*n (this file only ever calls with idx=i+-1/i+-2
! and i in [1,n], so idx in [-1,n+2] -- one +n is always enough).
  integer(i4) function cyc(idx, n) result(w)
    integer(i4), intent(in) :: idx, n
    w = mod(idx - 1 + n, n) + 1
  end function cyc

end subroutine co_norm_closed
