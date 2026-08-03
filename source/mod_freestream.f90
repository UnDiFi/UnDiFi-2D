module mod_freestream
! Phase 2.5 (ROADMAP.md #13): replaces shock.com's COMMON/SHOCKCOM/.
! Shared scratch state for the Roe-averaged upstream(v)/downstream(m) flow
! quantities either side of a shock/discontinuity point, threaded
! implicitly between co_state_dps and fx_state_dps. Initialized in
! shockini (see co_state_dps.f90).

  use mod_kinds, only: wp
  implicit none(type, external)
  public

  real(wp) :: uv, av, pv, rv
  real(wp) :: um, am, pm, rm
  real(wp) :: z1m, z2m, z3m, z4m
  real(wp) :: z1v, z2v, z3v, z4v
  real(wp) :: wsavg, wsmax

end module mod_freestream
