module mod_timer
! Per-call wall-clock timing for the console 'LABEL -->  ok' trace
! printed throughout main.f90's timestep loop, mod_shock_advance.f90's
! advance(), and the flow_solver_t backends (mod_neo_solver.f90,
! mod_eulfs_solver.f90, mod_su2_solver.f90). Added to let a run
! attribute wall time per named stage -- motivated by the OpenMP
! thread-count sweep (ROADMAP.md #16, Phase 4) showing no speedup at
! 1-8 threads even on a shock-point-refined case, so it's unclear
! without per-stage numbers whether the external system() calls
! (triangle/NEO/EulFS/SU2) or the parallelized fitting kernels dominate
! wall time.
!
! timer_tic/timer_toc share a single module-level start time rather
! than a caller-supplied handle: every call site in this codebase is
! strictly sequential (tic immediately followed by the timed work, then
! toc, never nested or concurrent across labels), so no caller-local
! timer variable is needed. timer_toc returns an already-formatted,
! fixed-length string so call sites can concatenate it onto their
! existing ' ok' write without touching any FORMAT statement.

  use, intrinsic :: iso_fortran_env, only: int64
  use mod_kinds, only: wp
  implicit none(type, external)
  private
  public :: timer_tic, timer_toc

  integer(int64), save :: t_start = 0_int64
  integer(int64), save :: t_rate = 1_int64

contains

  subroutine timer_tic()
    call system_clock(count=t_start, count_rate=t_rate)
  end subroutine timer_tic

  function timer_toc() result(str)
    character(len=16) :: str
    integer(int64) :: t_end
    real(wp) :: elapsed

    call system_clock(count=t_end)
    elapsed = real(t_end - t_start, wp)/real(t_rate, wp)
    write (str, '("  ",f11.4," s")') elapsed
  end function timer_toc

end module mod_timer
