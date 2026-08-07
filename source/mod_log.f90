module mod_log
! Phase 4.5 (ROADMAP.md #16), first increment: a thread-safe log-line
! writer. The ~26 files with `open(unit, file='log/*.log')` each build
! their message with plain `write`, which races if the enclosing loop
! is ever put under `!$omp parallel do` (Phase 4.2/4.3). log_line takes
! an already-formatted line and funnels the actual I/O through a single
! named critical section, so concurrent callers serialize on the write
! itself instead of corrupting interleaved output. Callers build their
! line into a local (thread-private) buffer via an internal write first
! -- only the final `write(unit,...)` to the shared file needs the lock.
! Converting the full ~26-file audit is out of scope here; this lands
! just the primitive plus the two files Phase 4.2 actually parallelizes
! (co_norm.f90, co_state_dps.f90).

  use mod_kinds, only: i4
  implicit none(type, external)
  private
  public :: log_line

contains

  subroutine log_line(unit, line)
    integer(i4), intent(in) :: unit
    character(len=*), intent(in) :: line

    !$omp critical (mod_log_io)
    write (unit, '(a)') trim(line)
    !$omp end critical (mod_log_io)
  end subroutine log_line

end module mod_log
