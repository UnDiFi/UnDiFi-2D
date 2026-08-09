module test_two_shock_theory
! Issue #25 (WP-A/M5), task "add analytic unit tests: 2-shock theory for
! RR": verifies co_urr against the classical theta-beta-M oblique-shock
! relations (Ben-Dor, "Shock Wave Reflection Phenomena"), NOT the
! abstract manufactured construction used in Phase 7.5's own co_urr
! kernel test (ROADMAP.md #20, tests/unit/test_co_urr.f90) -- that test
! only checks co_urr's Newton solve is internally self-consistent (any
! state3 satisfying the R-H + wall-tangency equations); this one checks
! the PHYSICS is right by comparing against an independent, literature-
! standard computation of what a real regular reflection should produce.
!
! Classical 2-shock-theory RR: an incident shock turns the freestream by
! deflection angle theta1 into region 2 (flow no longer parallel to the
! wall); the reflected shock must turn region 2's flow by theta1 in the
! opposite sense so region 3's flow is parallel to the wall again. Given
! region 2's Mach number M2 and its flow angle theta1 relative to the
! wall, the reflected-shock angle beta2 (measured from region 2's OWN
! flow direction, the standard theta-beta-M convention) is found by
! solving the same relation used for the incident shock, then region 3's
! state follows from the standard oblique-shock jump formulas. This test
! only needs region 2 (co_urr's actual upstream input) -- it does not
! need to also model the incident shock/region 1, since co_urr's own
! contract takes region 2 as a given.
!
! The shock's orientation (sx23, co_urr's y(14)) is derived from the
! velocity vectors themselves rather than guessed from theta-beta-M's
! own angle convention: tangential velocity is conserved across a shock,
! so (velocity2 - velocity3) is necessarily purely along the shock
! NORMAL direction -- a convention-independent way to get the shock
! tangent right by construction. Construction verified numerically (all
! 5 residuals ~1e-15, and the theta-beta-M solver's own round-trip
! theta check) via a standalone Python script before writing the
! Fortran below.
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_two_shock_theory

  real(wp), parameter :: pi = 3.14159265358979323846_wp

contains

  subroutine collect_two_shock_theory(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("regular_reflection_matches_2shock_theory", test_rr)&
    &]
  end subroutine collect_two_shock_theory

! Classical theta-beta-M relation: tan(theta) = 2*cot(beta) *
! (M**2*sin(beta)**2 - 1) / (M**2*(gam+cos(2*beta)) + 2)
  real(wp) function theta_from_beta(beta, M, gam) result(theta)
    real(wp), intent(in) :: beta, M, gam
    theta = atan(2.0_wp/tan(beta)*(M**2*sin(beta)**2 - 1.0_wp)/&
    &(M**2*(gam + cos(2.0_wp*beta)) + 2.0_wp))
  end function theta_from_beta

! Weak-shock root: bisection over beta in (asin(1/M), pi/2), which
! brackets the physically relevant (attached, weak) solution for any
! theta below the detachment angle.
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

! Standard oblique-shock jump across a shock of angle beta (measured
! from the upstream flow direction) in an upstream flow of Mach M1.
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

  subroutine test_rr(error)
    type(error_type), allocatable, intent(out) :: error

    include 'paramt.h'
    external :: co_urr

    real(wp), parameter :: gam = 1.4_wp
    real(wp), parameter :: rom = 1.0_wp, pm = 1.0_wp, M2 = 2.0_wp
    real(wp), parameter :: theta1 = 15.0_wp*pi/180.0_wp

    real(wp) :: beta2, p_ratio, rho_ratio, M3
    real(wp) :: p3, ro3, a2, a3
    real(wp) :: u2, v2, u3, v3
    real(wp) :: dx, dy, nlen, nx, ny, taux, tauy, sx23_true
    real(wp) :: y(15), yn1(15)

    ga = gam
    gm1 = gam - 1.0_wp

    beta2 = solve_beta_weak(M2, theta1, gam)
    call oblique_jump(M2, beta2, gam, p_ratio, rho_ratio, M3)

    p3 = pm*p_ratio
    ro3 = rom*rho_ratio
    a2 = sqrt(gam*pm/rom)
    a3 = sqrt(gam*p3/ro3)

    u2 = M2*a2*cos(theta1); v2 = M2*a2*sin(theta1)
    u3 = M3*a3; v3 = 0.0_wp   ! region 3: parallel to the wall (1,0)

    dx = u2 - u3; dy = v2 - v3
    nlen = sqrt(dx**2 + dy**2)
    nx = dx/nlen; ny = dy/nlen
    taux = ny; tauy = -nx
    sx23_true = atan2(tauy, taux)

    y(1:4) = [rom, pm, u2, 0.0_wp]  ! state1: physically unused, any value
    y(5:8) = [rom, pm, u2, v2]      ! state2: exact fixed input
    y(9:12) = [1.05_wp*ro3, 0.95_wp*p3, 1.05_wp*u3, 0.02_wp]  ! perturbed guess
    y(13) = 0.0_wp
    y(14) = sx23_true - 0.03_wp     ! perturbed shock-angle guess
    y(15) = 0.0_wp                  ! wrr: stationary reflection point

    call co_urr(y, 1.0_wp, 0.0_wp, yn1)  ! wall tangent (tauwx,tauwy)=(1,0)

    call check(error, yn1(9), ro3, rel=.true., thr=1.0e-6_wp,&
    &message="region-3 density doesn't match 2-shock theory")
    if (allocated(error)) return
    call check(error, yn1(10), p3, rel=.true., thr=1.0e-6_wp,&
    &message="region-3 pressure doesn't match 2-shock theory")
    if (allocated(error)) return
    call check(error, yn1(11), u3, rel=.true., thr=1.0e-6_wp,&
    &message="region-3 u-velocity doesn't match 2-shock theory")
    if (allocated(error)) return
    call check(error, yn1(12), v3, thr=1.0e-6_wp,&
    &message="region-3 v-velocity (wall tangency) doesn't match 2-shock theory")
    if (allocated(error)) return
  end subroutine test_rr

end module test_two_shock_theory
