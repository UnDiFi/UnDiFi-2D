! Shock genesis: create a fitted front from the singularity-formation and
! shock-development theory instead of from a tuned gradient sensor.
! ROADMAP.md item 192 / issue #185.  Manuscript and derivations:
! research/R-192-shock-genesis/.
!
! WHY THIS EXISTS.  A fitting code cannot create a shock: the front set is read
! from sh00.dat (re_sdw_info.f90) and thereafter only redistributed along its own
! arclength (rd_dps.f90).  The usual workaround -- watch a gradient sensor, plant
! a front, call co_shock -- fails for a structural reason, not a robustness one.
! At zero strength co_shock's Jacobian has
!
!     det J = 2*[p]/rho + O([p]^2)                                    (Thm 2.1)
!
! because a zero-strength discontinuity satisfies mass, momentum and energy at
! ANY front speed, so the speed column dies with the jump; and the 3x3
! conservation minor in (rho,p,u) is exactly (a^2 - v^2)/(gamma-1), which
! vanishes at sonic -- and a newborn front is sonic on both sides.  Both were
! confirmed against the compiled co_shock (research/.../undifi_probe/).
!
! The consequence is NOT that the solve becomes data-sensitive: it does not, the
! map is well posed in the data.  It is that the signal disappears.  What has to
! be resolved is the front's excess over the local characteristic speed, which is
! itself O([p]), while co_shock stops on an absolute criterion (the L1 norm of
! the damped Newton step, tol = 1e-7, damping = 0.2).  Measured on the production
! kernel, that excess is resolved to better than 1 % only for [p] > 5.4e-5 and is
! wrong by more than 100 % below [p] = 5.4e-7.  Waiting is the only remedy, and
! this module computes how long to wait rather than tuning it.
!
! WHAT IT PROVIDES.
!   genesis_fit_profile   -- spatial estimator: fits the local quartic of the
!                            characteristic speed and returns the modulation
!                            constants A, B, E of the self-similar attractor,
!                            hence t*, with an rms residual as its error bar.
!   genesis_fit_history   -- temporal estimator: -1/min(dlambda/dx) decays
!                            LINEARLY to zero at t*, so a straight-line fit gives
!                            t* independently, with no threshold anywhere.  Its
!                            slope must come out at -1; that is a second check.
!   genesis_emit          -- the genesis operator: birth position, speed, both
!                            states and the entropy jump at t* + dt, with every
!                            leading constant known in closed form.
!   genesis_handover_time -- the computable crossover at which co_shock may take
!                            over from the expansion.
!   genesis_selfsim       -- the attractor y = -U^3 - U, inverted exactly.
!
! The two t* estimators are genuinely independent -- one uses the spatial
! structure at one instant, the other the time history of one scalar -- so
! requiring them to AGREE WITHIN THEIR OWN RESIDUALS is a falsification test, not
! a consistency check.  genesis_fit_agree performs it.
!
! DEPENDENCIES.  mod_kinds only.  No COMMON, no include, no I/O, no solver
! coupling: this is meant to sit above mod_solver_iface as a CouplingData
! consumer and stay solver-agnostic across EulFS, NEO and SU2, so it carries its
! own 5x5 least-squares solve rather than reaching for solg (fnd_phps.f90),
! which would drag in the special-point closure.
!
! SCOPE, stated up front.  One-dimensional, and shocks born from smooth
! steepening only -- not fronts present in the initial or boundary data, which
! the code already handles.  Non-generic preshocks forced by an exact symmetry
! have different exponents and must be REFUSED via the fit residual rather than
! mishandled; that is what the ok flags are for.  Nothing can be extrapolated
! from an unresolved preshock: the estimator inherits the discretisation error of
! the field it is fitted to.
module mod_genesis
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private

  public :: genesis_fit_t, genesis_data_t
  public :: genesis_fit_profile, genesis_fit_history, genesis_fit_agree
  public :: genesis_emit, genesis_handover_time, genesis_selfsim

! Default half-width of the least-squares window, in nodes.  Not a tuned
! threshold: the fit residual reports whether the window was a good choice, and
! genesis_fit_profile refuses rather than guesses when it was not.
  integer(i4), parameter, public :: genesis_nwin_default = 6_i4

! Relative rms residual above which a profile fit is refused.  This is the one
! number in the module that is a judgement rather than a derivation; it is a
! goodness-of-fit bound, not a shock-detection threshold, and the distinction
! matters -- a sensor threshold decides whether a shock exists, this decides
! whether the local field is quartic enough for the expansion to mean anything.
  real(wp), parameter, public :: genesis_resid_max = 1.0e-3_wp

! Solver-floor constant c_tau of Cor 2.4, c_tau = (tol/damping)*4*gamma/((gamma+1)*a).
! tol = 1e-7 and damping = 0.2 are co_shock.f90's own hardcoded values; the
! measured value at gamma = 7/5, rho = p = 1 is 6.5e-7.
  real(wp), parameter, public :: genesis_ctau_ref = 6.5e-7_wp

! Round-off floor for the agreement test below.  The profile fit goes through
! normal equations, which square the condition number, so the best achievable
! relative accuracy is ~sqrt(epsilon) rather than epsilon -- roughly 1.5e-8 in
! double precision.  Without this floor the test fires on cases where BOTH
! estimators are exact to eleven digits and disagree in the twelfth, which is a
! false alarm and was one: it was found by genesis_validate.f90 on the datum
! lambda_0 = -sin(x), where the residuals are 1e-11 and 9e-16 while the two t*
! estimates differ by 1.3e-10.
  real(wp), parameter, public :: genesis_agree_floor = 1.0e-8_wp

! ---------------------------------------------------------------------------
! Output of either t* estimator.
  type :: genesis_fit_t
    real(wp) :: x0    = 0.0_wp   ! location of the steepening point
    real(wp) :: lam0  = 0.0_wp   ! characteristic speed there (= birth front speed)
    real(wp) :: amod  = 0.0_wp   ! A = -lambda'(x0)          > 0
    real(wp) :: bmod  = 0.0_wp   ! B = lambda'''(x0)/6       > 0
    real(wp) :: emod  = 0.0_wp   ! E = lambda''''(x0)/24     (sets the offset)
    real(wp) :: kmod  = 0.0_wp   ! K = A^2/sqrt(B), the single modulation constant
    real(wp) :: tstar = 0.0_wp   ! absolute singularity time
    real(wp) :: resid = 0.0_wp   ! relative rms fit residual -- the error bar
    real(wp) :: slope = 0.0_wp   ! history estimator only: must come out at -1
    logical  :: ok    = .false.
  end type genesis_fit_t

! Birth data emitted by the operator.
  type :: genesis_data_t
    real(wp) :: t     = 0.0_wp   ! absolute time = tstar + dt
    real(wp) :: xs    = 0.0_wp   ! front position
    real(wp) :: w     = 0.0_wp   ! front speed
    real(wp) :: dlam  = 0.0_wp   ! jump in the characteristic speed u + a
    real(wp) :: du    = 0.0_wp
    real(wp) :: dp    = 0.0_wp
    real(wp) :: ds_cv = 0.0_wp   ! entropy jump / c_v
    real(wp) :: rho_u = 0.0_wp, p_u = 0.0_wp, u_u = 0.0_wp   ! ahead
    real(wp) :: rho_d = 0.0_wp, p_d = 0.0_wp, u_d = 0.0_wp   ! behind
    logical  :: ok    = .false.
  end type genesis_data_t

contains

! ---------------------------------------------------------------------------
! The stable self-similar Burgers profile, inverted exactly.
!
!   -U/2 + (3y/2 + U) U' = 0    <=>    y = -U^3 - U
!
! so U is the real root of U^3 + U + y = 0.  The depressed cubic has p = 1 > 0,
! hence discriminant (y/2)^2 + 1/27 > 0 always: exactly one real root, given by
! Cardano in closed form with no branch choice to make.  The far field
! U ~ -y^(1/3) is the C^(1/3) preshock cusp.
  pure function genesis_selfsim(y) result(u)
    real(wp), intent(in) :: y
    real(wp) :: u, q, disc, s1, s2

    q = -0.5_wp*y
    disc = sqrt(q*q + 1.0_wp/27.0_wp)
    s1 = q + disc
    s2 = q - disc
    u = sign(abs(s1)**(1.0_wp/3.0_wp), s1) + sign(abs(s2)**(1.0_wp/3.0_wp), s2)
  end function genesis_selfsim

! ---------------------------------------------------------------------------
! SPATIAL ESTIMATOR.
!
! Fits lambda(x) = u + a, the genuinely nonlinear characteristic speed, by a
! quartic least squares in a window around the node where lambda' is most
! negative, then moves the origin to the exact inflection (lambda'' = 0) and
! reads off the modulation constants.  For Burgers -- and lambda obeys Burgers
! exactly in a simple wave -- the first singularity is at t = t_now + 1/A.
!
! x must be increasing.  Nothing here assumes a uniform spacing.
  subroutine genesis_fit_profile(x, lam, n, t_now, fit, nwin)
    real(wp), intent(in) :: x(:), lam(:)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: t_now
    type(genesis_fit_t), intent(out) :: fit
    integer(i4), intent(in), optional :: nwin

    integer(i4) :: i, imin, ilo, ihi, nw, m
    real(wp) :: dlmin, dl, h, xc, xi, xi0, c(5), rms, spread
    real(wp) :: c0, c1, c2, c3, c4, aa, bb, cc, r1, r2, lp, lpp3, lam0

    fit = genesis_fit_t()
    nw = genesis_nwin_default
    if (present(nwin)) nw = nwin
    if (n < 2*nw + 1) return

!   node where the characteristic speed steepens fastest (most negative slope)
    imin = 0
    dlmin = 0.0_wp
    do i = 2, n - 1
      dl = (lam(i + 1) - lam(i - 1))/(x(i + 1) - x(i - 1))
      if (dl < dlmin) then
        dlmin = dl
        imin = i
      end if
    end do
!   no compression anywhere: there is no shock forming, and saying so is the
!   correct answer, not a failure
    if (imin == 0) return

    ilo = max(1_i4, imin - nw)
    ihi = min(n, imin + nw)
    m = ihi - ilo + 1
    if (m < 5) return

!   scale to xi = (x - xc)/h so the normal equations stay conditioned
    xc = x(imin)
    h = 0.5_wp*(x(ihi) - x(ilo))
    if (h <= 0.0_wp) return

    call lsq_quartic(x(ilo:ihi), lam(ilo:ihi), m, xc, h, c, rms, fit%ok)
    if (.not. fit%ok) return

    spread = maxval(lam(ilo:ihi)) - minval(lam(ilo:ihi))
    if (spread > 0.0_wp) then
      fit%resid = rms/spread
    else
      fit%resid = rms
    end if

    c0 = c(1); c1 = c(2); c2 = c(3); c3 = c(4); c4 = c(5)

!   inflection: lambda'' = 0  <=>  12 c4 xi^2 + 6 c3 xi + 2 c2 = 0.
!   Take the root inside the window and nearest the centre; fall back to the
!   linear root when the quartic term is negligible.
    aa = 12.0_wp*c4
    bb = 6.0_wp*c3
    cc = 2.0_wp*c2
    if (abs(aa) <= epsilon(1.0_wp)*max(abs(bb), abs(cc))) then
      if (abs(bb) <= tiny(1.0_wp)) then
        fit%ok = .false.
        return
      end if
      xi0 = -cc/bb
    else
      dl = bb*bb - 4.0_wp*aa*cc
      if (dl < 0.0_wp) then
        fit%ok = .false.
        return
      end if
      r1 = (-bb + sqrt(dl))/(2.0_wp*aa)
      r2 = (-bb - sqrt(dl))/(2.0_wp*aa)
      xi0 = merge(r1, r2, abs(r1) <= abs(r2))
    end if
    if (abs(xi0) > 1.0_wp) then
!     the inflection is outside the fitted window: the window is in the wrong
!     place, and extrapolating out of it is exactly what must not be done
      fit%ok = .false.
      return
    end if

    xi = xi0
    lam0 = c0 + c1*xi + c2*xi**2 + c3*xi**3 + c4*xi**4
    lp = (c1 + 2.0_wp*c2*xi + 3.0_wp*c3*xi**2 + 4.0_wp*c4*xi**3)/h
    lpp3 = (6.0_wp*c3 + 24.0_wp*c4*xi)/h**3

    fit%x0 = xc + h*xi
    fit%lam0 = lam0
    fit%amod = -lp
    fit%bmod = lpp3/6.0_wp
    fit%emod = c4/h**4

!   A > 0 is compression; B > 0 is the generic (non-degenerate) preshock.  A
!   non-generic preshock forced by symmetry has B = 0 and different exponents,
!   and is refused here rather than mishandled.
    if (fit%amod <= 0.0_wp .or. fit%bmod <= 0.0_wp) then
      fit%ok = .false.
      return
    end if

    fit%kmod = fit%amod**2/sqrt(fit%bmod)
    fit%tstar = t_now + 1.0_wp/fit%amod
    fit%ok = (fit%resid <= genesis_resid_max)
  end subroutine genesis_fit_profile

! ---------------------------------------------------------------------------
! TEMPORAL ESTIMATOR, independent of the spatial one.
!
! For Burgers, lambda_x(t) = lambda_0'/(1 + t lambda_0'), so
!
!     mu(t) := -1/min_x lambda_x(t)   =   -1/lambda_0'  -  t
!
! decays LINEARLY to zero at t = t*, with slope exactly -1.  This is the inverse
! foliation density of the acoustical framework in its one-dimensional form: the
! Jacobian of the characteristic map, the quantity Moretti's own birth criterion
! (coalescence of characteristics) was reaching for.  Two things come free: t* as
! the zero crossing of a straight-line fit, predicted IN ADVANCE, and the fit
! residual as its own error bar.  There is no threshold to tune.
!
! Pass the history of min_x(dlambda/dx) -- negative -- not mu itself.
  subroutine genesis_fit_history(t, lamx_min, n, fit)
    real(wp), intent(in) :: t(:), lamx_min(:)
    integer(i4), intent(in) :: n
    type(genesis_fit_t), intent(out) :: fit

    integer(i4) :: k, m
    real(wp) :: mu, st, smu, stt, stmu, det, slope, inter, r, ss, mubar

    fit = genesis_fit_t()
    if (n < 3) return

    st = 0.0_wp; smu = 0.0_wp; stt = 0.0_wp; stmu = 0.0_wp; m = 0_i4
    do k = 1, n
      if (lamx_min(k) >= 0.0_wp) cycle            ! no compression at this time
      mu = -1.0_wp/lamx_min(k)
      st = st + t(k); smu = smu + mu
      stt = stt + t(k)*t(k); stmu = stmu + t(k)*mu
      m = m + 1_i4
    end do
    if (m < 3) return

    det = real(m, wp)*stt - st*st
    if (abs(det) <= tiny(1.0_wp)) return
    slope = (real(m, wp)*stmu - st*smu)/det
    inter = (smu*stt - st*stmu)/det
    if (slope >= 0.0_wp) return                   ! mu is not decaying: no birth

    fit%slope = slope
    fit%tstar = -inter/slope
    fit%amod = -1.0_wp/inter                      ! mu(0) = 1/A when t is measured from 0

!   relative rms residual of the straight line -- the error bar
    r = 0.0_wp; ss = 0.0_wp; mubar = smu/real(m, wp)
    do k = 1, n
      if (lamx_min(k) >= 0.0_wp) cycle
      mu = -1.0_wp/lamx_min(k)
      r = r + (mu - (slope*t(k) + inter))**2
      ss = ss + (mu - mubar)**2
    end do
    r = sqrt(r/real(m, wp))
    if (ss > 0.0_wp) then
      fit%resid = r/sqrt(ss/real(m, wp))
    else
      fit%resid = r
    end if

!   the slope is a prediction, not a fitted parameter: it must be -1.  Its
!   departure from -1 is a second, independent diagnostic of the fit window.
    fit%ok = (fit%resid <= genesis_resid_max) .and. (abs(slope + 1.0_wp) <= 0.05_wp)
  end subroutine genesis_fit_history

! ---------------------------------------------------------------------------
! FALSIFICATION TEST.  The two estimators use disjoint information -- spatial
! structure at one instant versus the time history of one scalar -- so requiring
! them to agree WITHIN THEIR OWN RESIDUALS is a real test, and disagreement
! beyond that means the preshock is not generic or is not resolved.  Returns
! .true. only if both fits are ok and their t* estimates are consistent.
  logical function genesis_fit_agree(fp, fh, tol_rel) result(agree)
    type(genesis_fit_t), intent(in) :: fp, fh
    real(wp), intent(in), optional :: tol_rel

    real(wp) :: bar, gap, allow

    agree = .false.
    if (.not. (fp%ok .and. fh%ok)) return
    bar = 0.5_wp*abs(fp%tstar + fh%tstar)
    if (bar <= tiny(1.0_wp)) return
    gap = abs(fp%tstar - fh%tstar)/bar
!   the allowance is the estimators' OWN residuals, not a tuned tolerance --
!   floored at the round-off level of the normal equations, below which a
!   difference carries no information
    allow = max(fp%resid + fh%resid, genesis_agree_floor)
    if (present(tol_rel)) allow = allow + tol_rel
    agree = (gap <= allow)
  end function genesis_fit_agree

! ---------------------------------------------------------------------------
! THE GENESIS OPERATOR.  Emits a complete front at t = t* + dt.  Every leading
! constant below is computed offline in closed form (manuscript Thm 5.2); none
! is fitted.
!
!   strength   |[lambda]| = 2 K dt^(1/2),        K = A^2/sqrt(B)
!   speed      w = lambda(x0), the characteristic speed at the birth point
!   position   xs = x0 + w dt + (1/5) E A^3 B^(-2) dt^2
!              -- the dt^2 term is the whole displacement: at leading order the
!              odd symmetry Ubar(-y) = -Ubar(y) puts the front exactly at the
!              cusp, and that degeneracy is why the constants are closed-form.
!   states     [u] = 2 [lambda]/(gamma+1)  (simple wave),  [p] = rho a [u]
!   entropy    [s]/c_v = ((gamma^2-1)/(12 gamma^2)) ([p]/p)^3  ~  dt^(3/2)
!
! rho0, p0 are the state at the birth point at t*; u0 follows from lam0.
  subroutine genesis_emit(fit, gam, rho0, p0, dt, out)
    type(genesis_fit_t), intent(in) :: fit
    real(wp), intent(in) :: gam, rho0, p0, dt
    type(genesis_data_t), intent(out) :: out

    real(wp) :: a0, half, du_half, dp_half, ce

    out = genesis_data_t()
    if (.not. fit%ok) return
    if (dt <= 0.0_wp .or. rho0 <= 0.0_wp .or. p0 <= 0.0_wp) return

    a0 = sqrt(gam*p0/rho0)

    out%t = fit%tstar + dt
    out%w = fit%lam0
    out%xs = fit%x0 + fit%lam0*dt&
    &+ 0.2_wp*fit%emod*fit%amod**3/fit%bmod**2*dt*dt

    out%dlam = 2.0_wp*fit%kmod*sqrt(dt)
    out%du = 2.0_wp*out%dlam/(gam + 1.0_wp)
    out%dp = rho0*a0*out%du

!   split symmetrically about the birth state: the two branches are at
!   zeta = -1 and +1, i.e. lambda0 -/+ K sqrt(dt)
    half = 0.5_wp
    du_half = half*out%du
    dp_half = half*out%dp

!   ahead of the front (weaker branch)
    out%p_u = p0 - dp_half
    out%rho_u = rho0 - dp_half/(a0*a0)
    out%u_u = fit%lam0 - a0 - du_half
!   behind the front (stronger branch)
    out%p_d = p0 + dp_half
    out%rho_d = rho0 + dp_half/(a0*a0)
    out%u_d = fit%lam0 - a0 + du_half

    if (out%rho_u <= 0.0_wp .or. out%p_u <= 0.0_wp) return

!   Landau-Lifshitz weak-shock entropy law, rederived from the RH relations in
!   the manuscript rather than quoted: (gamma^2-1)/(12 gamma^2) = 2/49 at 7/5
    ce = (gam*gam - 1.0_wp)/(12.0_wp*gam*gam)
    out%ds_cv = ce*(out%dp/p0)**3

    out%ok = .true.
  end subroutine genesis_emit

! ---------------------------------------------------------------------------
! THE CROSSOVER.  Two errors run in opposite directions as the handover time
! grows: the expansion's own truncation error is relative O(A dt) and GROWS,
! while the error co_shock would leave behind is c_tau/(c_p sqrt(dt)) and
! SHRINKS -- the latter because its stopping rule is absolute while the signal
! it has to resolve is O([p]) = O(sqrt(dt)).  Equating them,
!
!     dt* = ( c_tau / (A c_p) )^(2/3) ,      c_p = 4 rho a K/(gamma+1)
!
! Every symbol is known before the run.  Note c_tau is a property of the
! SOLVER'S STOPPING RULE, not of the field: dt* is therefore unchanged by mesh
! refinement.  (An earlier form of this argument put the field truncation error
! here instead, on the assumption that data error is amplified by the
! ill-conditioning.  It is not -- the map is well posed in the data -- and that
! version would wrongly have predicted dt* to shrink with the mesh.)
  pure function genesis_handover_time(fit, gam, rho0, p0, ctau) result(dt)
    type(genesis_fit_t), intent(in) :: fit
    real(wp), intent(in) :: gam, rho0, p0
    real(wp), intent(in), optional :: ctau

    real(wp) :: dt, a0, cp, ct

    dt = 0.0_wp
    if (.not. fit%ok) return
    if (rho0 <= 0.0_wp .or. p0 <= 0.0_wp .or. fit%amod <= 0.0_wp) return

    a0 = sqrt(gam*p0/rho0)
    cp = 4.0_wp*rho0*a0*fit%kmod/(gam + 1.0_wp)
    ct = genesis_ctau_ref
    if (present(ctau)) ct = ctau
    if (cp <= 0.0_wp) return

    dt = (ct/(fit%amod*cp))**(2.0_wp/3.0_wp)
  end function genesis_handover_time

! ---------------------------------------------------------------------------
! Private: least-squares quartic through (x, y), origin at xc, scaled by h.
! Normal equations on a 5x5 symmetric system, solved by Gaussian elimination
! with partial pivoting.  Kept local so the module depends on mod_kinds alone.
  subroutine lsq_quartic(x, y, m, xc, h, c, rms, ok)
    real(wp), intent(in) :: x(:), y(:)
    integer(i4), intent(in) :: m
    real(wp), intent(in) :: xc, h
    real(wp), intent(out) :: c(5), rms
    logical, intent(out) :: ok

    integer(i4) :: i, j, k, ip
    real(wp) :: a(5, 6), xi, pw(5), piv, fac, f, r

    ok = .false.
    c = 0.0_wp
    rms = 0.0_wp
    a = 0.0_wp

    do k = 1, m
      xi = (x(k) - xc)/h
      pw(1) = 1.0_wp
      do j = 2, 5
        pw(j) = pw(j - 1)*xi
      end do
      do i = 1, 5
        do j = 1, 5
          a(i, j) = a(i, j) + pw(i)*pw(j)
        end do
        a(i, 6) = a(i, 6) + pw(i)*y(k)
      end do
    end do

    do i = 1, 5
      ip = i
      piv = abs(a(i, i))
      do k = i + 1, 5
        if (abs(a(k, i)) > piv) then
          piv = abs(a(k, i))
          ip = k
        end if
      end do
      if (piv <= tiny(1.0_wp)) return
      if (ip /= i) then
        do j = i, 6
          f = a(i, j); a(i, j) = a(ip, j); a(ip, j) = f
        end do
      end if
      do k = i + 1, 5
        fac = a(k, i)/a(i, i)
        do j = i, 6
          a(k, j) = a(k, j) - fac*a(i, j)
        end do
      end do
    end do

    do i = 5, 1, -1
      f = a(i, 6)
      do j = i + 1, 5
        f = f - a(i, j)*c(j)
      end do
      c(i) = f/a(i, i)
    end do

    do k = 1, m
      xi = (x(k) - xc)/h
      f = c(1) + xi*(c(2) + xi*(c(3) + xi*(c(4) + xi*c(5))))
      r = y(k) - f
      rms = rms + r*r
    end do
    rms = sqrt(rms/real(m, wp))
    ok = .true.
  end subroutine lsq_quartic

end module mod_genesis
