module test_co_shock
! Phase 7.5 increment 2 (ROADMAP.md #20): co_shock -- the 4x4 Newton solve
! for a moving normal shock (Rankine-Hugoniot mass/momentum/energy across
! the shock, closed by the downstream C+ characteristic invariant R14
! carried in from the interior solution) -- against the classical analytic
! normal-shock jump relations, across a shock-Mach-number and gamma sweep
! including the near-sonic and strong-shock limits.
!
! co_shock itself is a bare (non-module) subroutine reading GA from the
! COMMON/PARAMT/ block via `include 'paramt.h'`, same convention as every
! other production call site (see re_inp_data.f90); this test sets GA the
! same way any real driver run would, just without going through
! mod_config/input.dat.
  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, naddholesmax, nprdbndmax
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_co_shock

contains

  subroutine collect_co_shock(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("normal_shock_moderate", test_moderate),&
    &new_unittest("normal_shock_near_sonic", test_near_sonic),&
    &new_unittest("normal_shock_strong", test_strong),&
    &new_unittest("normal_shock_gamma_sweep", test_gamma_sweep)&
    &]
  end subroutine collect_co_shock

! Manufactured-solution check: for a given upstream state (rom, pm, um)
! and shock Mach number M1 (>1, relative to the shock), a shock held
! STATIONARY in this frame (w=0) gives an exact downstream density,
! pressure and velocity via the classical normal-shock relations. R14 is
! then computed from that exact state via the same characteristic formula
! co_shock's 4th residual uses, so feeding it back in -- together with a
! perturbed initial guess -- isolates whether the Newton solve recovers
! the *known* root rather than merely converging to *some*
! self-consistent-looking answer.
!
! w=0 is deliberate, not a simplification of convenience: co_shock's own
! Newton seed hardcodes the shock-speed unknown to a small constant
! (w = -0.001_wp, see co_shock.f90) on every call, regardless of x1/x2 --
! i.e. every real call site already assumes the shock is near-stationary
! in whatever frame co_shock is invoked in (consistent with the
! shock-fitting scheme's small time steps). A manufactured case with a
! large |w| starts the Newton solve far outside the basin that seed was
! designed for and converges to a spurious algebraic root of the same
! R-H equations instead (verified while developing this test) -- not a
! bug in co_shock, just not the regime it's built to be seeded for.
  subroutine check_normal_shock(error, gam, rom, pm, M1)
    type(error_type), allocatable, intent(out) :: error
    real(wp), intent(in) :: gam, rom, pm, M1

    include 'paramt.h'
    external :: co_shock

    real(wp) :: am, um, rho_ratio, p_ratio, delta, R14
    real(wp) :: rov_exact, pv_exact, uv_exact
    real(wp) :: x1(4), x2(4), wshk

    ga = gam
    gm1 = gam - 1.0_wp

    am = sqrt(gam*pm/rom)
    um = M1*am                 ! upstream velocity into a stationary shock (w=0)

    rho_ratio = (gam + 1.0_wp)*M1**2/((gam - 1.0_wp)*M1**2 + 2.0_wp)
    p_ratio = 1.0_wp + 2.0_wp*gam/(gam + 1.0_wp)*(M1**2 - 1.0_wp)

    rov_exact = rom*rho_ratio
    pv_exact = pm*p_ratio
    uv_exact = um/rho_ratio    ! continuity: rom*um = rov*uv (w=0 on both sides)

    delta = (gam - 1.0_wp)/2.0_wp
    R14 = sqrt(gam*pv_exact/rov_exact) + delta*uv_exact

    x2 = [rom, pm, um, 0.0_wp]
    ! Deliberately off the exact root so the Newton solve has to do real
    ! work to get back to it, not just confirm a lucky initial guess.
    x1 = [1.05_wp*rov_exact, 0.95_wp*pv_exact, 1.05_wp*uv_exact, 0.0_wp]

    call co_shock(x1, x2, wshk, R14)

    call check(error, x1(1), rov_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, x1(2), pv_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, x1(3), uv_exact, rel=.true., thr=1.0e-6_wp)
    if (allocated(error)) return
    call check(error, wshk, 0.0_wp, thr=1.0e-6_wp)
    if (allocated(error)) return
  end subroutine check_normal_shock

  subroutine test_moderate(error)
    type(error_type), allocatable, intent(out) :: error
    call check_normal_shock(error, 1.4_wp, 1.0_wp, 1.0_wp, 2.0_wp)
  end subroutine test_moderate

! Near-sonic limit (M1 -> 1+): the FD-Jacobian path elsewhere in
! mod_newton_solve is documented as fragile here; co_shock always uses its
! analytic Jacobian (jf), so this exercises that analytic path right at
! the edge of its validity, not the FD fallback.
  subroutine test_near_sonic(error)
    type(error_type), allocatable, intent(out) :: error
    call check_normal_shock(error, 1.4_wp, 1.0_wp, 1.0_wp, 1.02_wp)
  end subroutine test_near_sonic

  subroutine test_strong(error)
    type(error_type), allocatable, intent(out) :: error
    call check_normal_shock(error, 1.4_wp, 1.0_wp, 1.0_wp, 20.0_wp)
  end subroutine test_strong

! Different ratios of specific heats (1.667: monatomic gas; 1.2: a heavy
! polyatomic gas), each at a different Mach number, so this doubles as a
! sanity check that gam isn't hardcoded anywhere in the residual/Jacobian.
  subroutine test_gamma_sweep(error)
    type(error_type), allocatable, intent(out) :: error
    call check_normal_shock(error, 1.667_wp, 1.0_wp, 1.0_wp, 3.0_wp)
    if (allocated(error)) return
    call check_normal_shock(error, 1.2_wp, 1.0_wp, 1.0_wp, 5.0_wp)
  end subroutine test_gamma_sweep

end module test_co_shock
