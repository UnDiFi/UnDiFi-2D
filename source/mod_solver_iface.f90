module mod_solver_iface
! Phase 3.7 (ROADMAP.md #14, issue #15): abstract flow_solver_t, replacing
! main.f90's three-way if(EULFS)/elseif(SU2)/elseif(NEO) dispatch -- ~700
! lines spread across 7 call sites, not just the two obvious ~250-line
! blocks (see the increment commits for the full site-by-site map). Each
! solver is a separate external process (a binary in bindir, invoked via
! mod_run_external's run_external) -- "abstracting" this is a source-level
! relocation of command-string-building + run_external calls into
! type-bound procedures, not an FFI or build-system change; nothing here
! is gated by any CMake conditional (see ROADMAP.md status notes).
!
! Deferred bindings map onto ROADMAP.md's own "convert-in/run/convert-out"
! framing:
!   prepare  -- convert the just-regenerated Triangle mesh into the
!               solver's native input format
!   run      -- invoke the external binary for one outer-iteration step
!   harvest  -- convert the solver's output back into Triangle format
! Each takes one solver_run_ctx_t bundling the invocation parameters that
! were previously local variables in main.f90's iteration loop (Triangle
! basename, bindir/hostype/testcase, iteration counters, and a `corrector`
! flag selecting which of the two structurally-different call sites --
! predictor/always-run vs UNSTEADY-corrector-only -- is asking). One
! context type avoids a long deferred-interface argument list and matches
! this codebase's existing ctx-object idiom (mod_newton_solve's
! `class(*) :: ctx` Newton contexts). EulFS/SU2/NEO-specific inputs that
! come from paramt.h's COMMON block (nprdbnd, prdbndclr, ibak, imtf, ...)
! are NOT duplicated into the context -- each concrete type's own module
! just `include`s paramt.h itself, same as co_utp.f90/co_uqp.f90/etc.
! already do.
!
! Four more bindings are non-deferred with a shared no-op default,
! overridden only by the solver(s) that actually do something at that
! site -- this mirrors mod_special_point.f90's pattern of a shared
! no-op for behaviors most concrete types don't need:
!   setup            -- once, before the outer loop starts (EulFS-only:
!                        reset convergenza.dat, prime .petscrc)
!   pre_mesh_setup   -- once per iteration, BEFORE Triangle regenerates
!                        the mesh (NEO-only: first-iteration neogrid0.grd
!                        creation for the ShockExpansion testcase -- this
!                        genuinely has to run before Triangle, unlike
!                        prepare/run/harvest which all run after)
!   archive          -- at loop end, backup-or-cleanup file archiving
!                        (EulFS/NEO only; SU2 has never had an arm here)
!   log_convergence  -- at loop end, append this iteration's residual to
!                        convergenza.dat (EulFS only)
! supports_ale() is a capability query (default .false., true only for
! eulfs_t) replacing the `if (EULFS)` checks around the ALE `solzne` call
! -- both the standalone one in main.f90 and the one inside
! mod_shock_advance.f90's advance() (whose own `eulfs` logical argument
! is unchanged; main.f90 just passes `solver%supports_ale()` instead of
! the raw `eulfs` variable there).
!
! **SU2's known pre-existing gaps are deliberately NOT fixed here** (the
! user's explicit scoping decision, see the plan/commit message): there
! is no elseif(SU2) arm at 3 of the 7 sites in the legacy code (the
! UNSTEADY corrector step, and both halves of the loop-end archive/
! cleanup logic), so SU2 + UNSTEADY silently skips its corrector step
! entirely today, and gets no backup/cleanup archiving either. su2_t
! reproduces this as explicit no-ops at those sites (see mod_su2_solver.f90),
! not a mechanical port of nonexistent code. This is also, as of this
! writing, entirely untestable in this environment: SU2_CFD is not
! installed here, and zero regression-baseline coverage for any su2_*
! case exists in tests/baseline/ (confirmed via tests/baseline/README.md).

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private

  public :: flow_solver_t, solver_run_ctx_t

  ! Bundles the per-call invocation parameters that used to be plain
  ! main.f90 locals. fname/backdir/fnameback/bindir/hostype are fixed-
  ! length slices matching exactly how main.f90 already truncates them
  ! (e.g. `fname(1:7)`, `bindir(1:10)`) -- not full character*255 buffers.
  type :: solver_run_ctx_t
    character(len=7)  :: fname = ''      ! Triangle basename, e.g. "na00161"
    character(len=9)  :: backdir = ''    ! backup dir, e.g. "step00161" (archive() only)
    character(len=4)  :: fnameback = ''  ! "na99" (archive() only)
    character(len=10) :: bindir = ''     ! "../../bin/"
    character(len=6)  :: hostype = ''    ! "x86_64"
    character(len=20) :: testcase = ''
    integer(i4)       :: iter = 0
    integer(i4)       :: nbegin = 0
    logical           :: unsteady = .false.
    logical           :: corrector = .false.  ! .false.=predictor/always-run site, .true.=UNSTEADY-corrector-only site
    logical           :: backing_up = .false. ! archive() only: mod(iter-1,ibak)==0
  end type solver_run_ctx_t

  type, abstract :: flow_solver_t
  contains
    procedure(prepare_i), deferred :: prepare
    procedure(run_i), deferred :: run
    procedure(harvest_i), deferred :: harvest
    procedure :: setup => solver_default_setup
    procedure :: pre_mesh_setup => solver_default_pre_mesh_setup
    procedure :: archive => solver_default_archive
    procedure :: log_convergence => solver_default_log_convergence
    procedure :: supports_ale => solver_default_supports_ale
  end type flow_solver_t

  abstract interface
    subroutine prepare_i(self, ctx)
      import :: flow_solver_t, solver_run_ctx_t
      class(flow_solver_t), intent(inout) :: self
      type(solver_run_ctx_t), intent(in) :: ctx
    end subroutine prepare_i

    subroutine run_i(self, ctx)
      import :: flow_solver_t, solver_run_ctx_t
      class(flow_solver_t), intent(inout) :: self
      type(solver_run_ctx_t), intent(in) :: ctx
    end subroutine run_i

    subroutine harvest_i(self, ctx)
      import :: flow_solver_t, solver_run_ctx_t
      class(flow_solver_t), intent(inout) :: self
      type(solver_run_ctx_t), intent(in) :: ctx
    end subroutine harvest_i
  end interface

contains

  subroutine solver_default_setup(self, unsteady)
    class(flow_solver_t), intent(inout) :: self
    logical, intent(in) :: unsteady
  end subroutine solver_default_setup

  subroutine solver_default_pre_mesh_setup(self, ctx)
    class(flow_solver_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
  end subroutine solver_default_pre_mesh_setup

  subroutine solver_default_archive(self, ctx)
    class(flow_solver_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
  end subroutine solver_default_archive

  subroutine solver_default_log_convergence(self)
    class(flow_solver_t), intent(inout) :: self
  end subroutine solver_default_log_convergence

  logical function solver_default_supports_ale(self) result(ale)
    class(flow_solver_t), intent(in) :: self
    ale = .false.
  end function solver_default_supports_ale

end module mod_solver_iface
