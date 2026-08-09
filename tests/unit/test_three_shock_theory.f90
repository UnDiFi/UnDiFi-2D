module test_three_shock_theory
! Issue #25 (WP-A/M5), task "3-shock theory for the triple point":
! verifies co_utp against the classical theta-beta-M oblique-shock
! relations, applied per shock leg, rather than the abstract shared-
! angle construction used in Phase 7.5's own co_utp kernel test
! (ROADMAP.md #20, tests/unit/test_co_utp.f90).
!
! Classical 3-shock theory: two independent oblique shocks (the
! reflected shock 23 and the Mach stem 14) each turn their own upstream
! state (2, 1 respectively -- both fixed inputs to co_utp) toward a
! COMMON downstream flow direction phi_common, and must additionally
! agree on downstream pressure across the slip line between states 3
! and 4. Built here the same way as test_two_shock_theory.f90's RR case,
! applied independently to both legs: pick (M,phi,rho,p) for state1 and
! (M,phi,rho) for state2 freely, solve each leg's shock angle via
! theta-beta-M for the deflection needed to reach phi_common, then solve
! state2's base pressure in closed form so its resulting downstream
! pressure exactly matches state1's leg -- same "solve the free
! parameter in closed form" idea used for co_uqp's manufactured test in
! ROADMAP.md #20, just with theta-beta-M-derived deflections instead of
! an abstract Mach-number pick. Downstream velocities end up as
! different-magnitude, same-direction vectors (both along phi_common) --
! automatically parallel, matching co_utp's eq11.
!
! Each shock's orientation is derived from its own velocity difference
! (tangential-velocity continuity), same convention-independent
! technique as test_two_shock_theory.f90. Construction verified
! numerically (all 11 residuals ~1e-15) via a standalone Python script
! before writing the Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_three_shock_theory

  real(wp), parameter :: pi = 3.14159265358979323846_wp

contains

  subroutine collect_three_shock_theory(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("triple_point_matches_3shock_theory", test_tp)&
    &]
  end subroutine collect_three_shock_theory

  real(wp) function theta_from_beta(beta, M, gam) result(theta)
    real(wp), intent(in) :: beta, M, gam
    theta = atan(2.0_wp/tan(beta)*(M**2*sin(beta)**2 - 1.0_wp)/&
    &(M**2*(gam + cos(2.0_wp*beta)) + 2.0_wp))
  end function theta_from_beta

  real(wp) function solve_beta_weak(M, theta, gam) result(beta)
    real(wp), intent(in) :: M, theta, gam
    real(wp) :: lo, hi, mid, flo, fmid
    integer(i4), parameter :: n_bisect = 100
    integer(i4) :: i

    lo = asin(1.0_wp/M) + 1.0e-6_wp
    hi = pi/2.0_wp - 1.0e-6_wp
    flo = theta_from_beta(lo, M, gam) - theta
    do i = 1, n_bisect
      mid = 0.5_wp*(lo + hi)
      fmid = theta_from_beta(mid, M, gam) - theta
      if (flo*fmid <= 0.0_wp) then
        hi = mid
      else
        lo = mid; flo = fmid
      end if
    end do
    beta = 0.5_wp*(lo + hi)
  end function solve_beta_weak

  subroutine oblique_jump(M1, beta, gam, p_ratio, rho_ratio, M2)
    real(wp), intent(in) :: M1, beta, gam
    real(wp), intent(out) :: p_ratio, rho_ratio, M2
    real(wp) :: Mn1, delta, Mn2sq, theta

    Mn1 = M1*sin(beta)
    delta = (gam - 1.0_wp)/2.0_wp
    p_ratio = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(Mn1**2 - 1.0_wp)
    rho_ratio = (gam + 1.0_wp)*Mn1**2/((gam - 1.0_wp)*Mn1**2 + 2.0_wp)
    Mn2sq = (1.0_wp + delta*Mn1**2)/(gam*Mn1**2 - delta)
    theta = theta_from_beta(beta, M1, gam)
    M2 = sqrt(Mn2sq)/sin(beta - theta)
  end subroutine oblique_jump

! One shock leg: upstream (M_up, phi_up, rom, pm) deflected to
! phi_common. Returns the downstream state and the shock's own
! orientation (derived from the velocity difference, not phi/beta
! directly -- convention-independent, see module header).
  subroutine shock_leg(M_up, phi_up, rom, pm, phi_common, gam,&
  &u_up, v_up, ro_dn, p_dn, u_dn, v_dn, sx)
    real(wp), intent(in) :: M_up, phi_up, rom, pm, phi_common, gam
    real(wp), intent(out) :: u_up, v_up, ro_dn, p_dn, u_dn, v_dn, sx
    real(wp) :: theta, beta, p_ratio, rho_ratio, M_dn, a_up, a_dn
    real(wp) :: dx, dy, nlen, nx, ny, taux, tauy

    theta = phi_up - phi_common
    beta = solve_beta_weak(M_up, abs(theta), gam)
    call oblique_jump(M_up, beta, gam, p_ratio, rho_ratio, M_dn)

    p_dn = pm*p_ratio
    ro_dn = rom*rho_ratio
    a_up = sqrt(gam*pm/rom)
    a_dn = sqrt(gam*p_dn/ro_dn)

    u_up = M_up*a_up*cos(phi_up); v_up = M_up*a_up*sin(phi_up)
    u_dn = M_dn*a_dn*cos(phi_common); v_dn = M_dn*a_dn*sin(phi_common)

    dx = u_up - u_dn; dy = v_up - v_dn
    nlen = sqrt(dx**2 + dy**2)
    nx = dx/nlen; ny = dy/nlen
    taux = ny; tauy = -nx
    sx = atan2(tauy, taux)
  end subroutine shock_leg

  subroutine test_tp(error)
    type(error_type), allocatable, intent(out) :: error

    include 'paramt.h'
    external :: co_utp

    real(wp), parameter :: gam = 1.4_wp
    real(wp), parameter :: phi_common = 0.0_wp
    real(wp), parameter :: dxr14 = 1.0_wp, dyr14 = 0.0_wp
    real(wp), parameter :: M1 = 3.0_wp, phi1 = 20.0_wp*pi/180.0_wp
    real(wp), parameter :: rom1 = 1.0_wp, pm1 = 1.0_wp
    real(wp), parameter :: M2 = 2.2_wp, phi2 = 12.0_wp*pi/180.0_wp
    real(wp), parameter :: rom2 = 1.3_wp

    real(wp) :: u1, v1, ro4, p4, u4, v4, sx14
    real(wp) :: u2, v2, ro3, p3, u3, v3, sx23
    real(wp) :: pm2, p_ratio23, rho_ratio_prep, beta23_prep, theta23_prep, M_prep_unused
    real(wp) :: r14, y(20), yn1(20)
    logical :: ifail

    ga = gam
    gm1 = gam - 1.0_wp

    ! block 2 (shock14, Mach stem): state1 -> state4, fully determined
    call shock_leg(M1, phi1, rom1, pm1, phi_common, gam,&
    &u1, v1, ro4, p4, u4, v4, sx14)

    ! block 1 (shock23, reflected): state2 -> state3; solve pm2 in
    ! closed form so p3 comes out exactly equal to p4 (eq10). shock_leg
    ! itself needs a base pressure to scale p_ratio by, so compute the
    ! ratio first with a throwaway pm=1, then rescale.
    theta23_prep = phi2 - phi_common
    beta23_prep = solve_beta_weak(M2, abs(theta23_prep), gam)
    call oblique_jump(M2, beta23_prep, gam, p_ratio23, rho_ratio_prep, M_prep_unused)
    pm2 = p4/p_ratio23

    call shock_leg(M2, phi2, rom2, pm2, phi_common, gam,&
    &u2, v2, ro3, p3, u3, v3, sx23)

    r14 = sqrt(gam*p4/ro4) + (gam - 1.0_wp)/2.0_wp*(u4*dxr14 + v4*dyr14)

    y(1:4) = [rom1, pm1, u1, v1]        ! state1: exact fixed input
    y(5:8) = [rom2, pm2, u2, v2]        ! state2: exact fixed input
    y(9:12) = [1.05_wp*ro3, 0.95_wp*p3, 1.05_wp*u3, 0.05_wp]   ! perturbed
    y(13:16) = [0.95_wp*ro4, 1.05_wp*p4, 0.95_wp*u4, -0.05_wp] ! perturbed
    y(17) = 0.4_wp                      ! sx12: unused (wsh_true=0)
    y(18) = sx23 - 0.03_wp              ! perturbed shock-23-angle guess
    y(19) = sx14 + 0.02_wp              ! perturbed shock-14-angle guess
    y(20) = 0.01_wp                     ! perturbed wsh guess (true value 0)

    call co_utp(y, r14, dxr14, dyr14, 0.0_wp, 0.0_wp, yn1, .false., ifail)

    call check(error, .not. ifail, message="co_utp reported ifail")
    if (allocated(error)) return
    call check(error, yn1(9), ro3, rel=.true., thr=1.0e-6_wp,&
    &message="state3 density doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(10), p3, rel=.true., thr=1.0e-6_wp,&
    &message="state3 pressure doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(11), u3, rel=.true., thr=1.0e-6_wp,&
    &message="state3 u-velocity doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(12), v3, thr=1.0e-6_wp,&
    &message="state3 v-velocity doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(13), ro4, rel=.true., thr=1.0e-6_wp,&
    &message="state4 density doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(14), p4, rel=.true., thr=1.0e-6_wp,&
    &message="state4 pressure doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(15), u4, rel=.true., thr=1.0e-6_wp,&
    &message="state4 u-velocity doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(16), v4, thr=1.0e-6_wp,&
    &message="state4 v-velocity doesn't match 3-shock theory")
    if (allocated(error)) return
    call check(error, yn1(10) - yn1(14), 0.0_wp, thr=1.0e-6_wp,&
    &message="converged state3/state4 pressures don't match (slip line)")
    if (allocated(error)) return
  end subroutine test_tp

end module test_three_shock_theory
