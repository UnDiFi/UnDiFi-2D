module mod_neo_solver
! Phase 3.7 increment 3 (ROADMAP.md #14, issue #15): neo_t, a mechanical
! port of every NEO-specific branch of main.f90's former
! if(EULFS)/elseif(SU2)/elseif(NEO) dispatch. Every literal command
! string, flag, and comment is preserved exactly, INCLUDING the small
! cosmetic differences (progress-label spacing, 'ok' vs ' ok', a
! run_external name arg's case) between the predictor and corrector
! call sites -- these only affect stdout text, never a checksummed
! result, but are kept faithfully rather than "cleaned up" into a
! single shared string, same rigor as mod_shock_advance.f90's port.
!
! Phase 5a (ROADMAP.md #17), 2026-08-08: na2vvvv/triangle2grd/
! neo2triangle's neo_prepare/neo_harvest call sites now call the
! in-process ports (na2vvvv.f90, triangle2grd.f90, neo2triangle.f90)
! directly instead of spawning the source_utils/ standalone executables
! via run_external -- see na2vvvv.f90 for the full rationale (the
! issue's own benchmark task showed this round-trip, not the NEO solve
! itself, dominates wall time). NEO's own solve (neo_run) and the
! neogrid0/archive run_external calls are untouched -- out of this
! slice's scope, still spawned as before.
!
! KNOWN BEHAVIOR CHANGE, not yet resolved: the corrector-side na2vvvv
! call used to pass fatal_on_error=.false. to run_external, so a
! nonzero na2vvvv exit status (e.g. a missing .node file) would log and
! continue rather than abort. The in-process na2vvvv has no equivalent
! -- an open()/read() failure now hard-stops the whole UnDiFi-2D
! process (loses the subprocess-level fault isolation the old
! architecture gave for free). Never observed to fire in the regression
! suite; flagged here rather than silently dropped in case it mattered
! for some untested case.
!
! Phase 5a increment 2, 2026-08-08: neo_prepare's na2vvvv and
! triangle2grd calls now run in `!$omp parallel sections` -- both read
! the same just-regenerated Triangle mesh read-only and write to
! disjoint files (vvvv_input.dat vs neogrid.grd), confirmed by tracing
! every open() in both, so there is no shared-state hazard the way
! there was for co_state_dps/co_norm in Phase 4.2. This is real
! measured wall time (~8.6% of one outer iteration, architecture-review
! finding 2026-08-08), unlike Phase 4's fitting-kernel OpenMP which
! parallelizes ~1% of runtime -- but it only takes effect under the
! opt-in UNDIFI_ENABLE_OPENMP build, same as every other Phase 4/5a
! `!$omp` in this codebase; the default build stays fully sequential.
! mod_timer's timer_tic/timer_toc are NOT used inside the parallel
! region (their single module-level start time would race across the
! two section threads) -- each section takes its own system_clock
! reading into a local variable instead, printed after the join.
! Original main.f90 line numbers (pre-3.7, commit c5b9cfb) cited below.
!
! pre_mesh_setup()  <- lines 730-737 (once per outer iteration, called
!                       from main.f90 at that exact position -- BEFORE
!                       Triangle regenerates the mesh, unlike
!                       prepare/run/harvest which all run after; no
!                       corrector-side equivalent exists, so this is
!                       only ever invoked from the predictor position)
! prepare()         <- lines 940-963 (predictor: na2vvvv gated by a
!                       STEADY/first-iteration/ShockVortex guard, then
!                       triangle2grd) / 1200-1217 (corrector: na2vvvv
!                       unconditional, then triangle2grd) -- the
!                       dropped guard on the corrector side is a
!                       verified pre-existing asymmetry (it only ever
!                       matters on i==1+nbegin, and the corrector call
!                       site never fires until the loop's second half
!                       of an already-not-first iteration in practice),
!                       not something to "fix" into symmetry.
! run()             <- lines 965-1019 (predictor: an UNSTEADY-and-
!                       first-iteration-only CRD_euler bootstrap +
!                       inputfile-exp.txt swap, then an unconditional
!                       CRD_euler call that also fires on iteration 1 --
!                       both really do run on iter 1, not a double-
!                       counting bug introduced here) / 1223-1227
!                       (corrector: the single CRD_euler call only)
! harvest()         <- lines 1021-1033 / 1229-1237 (NEO2triangle)
!
! archive() and log_convergence() are NOT overridden here for the
! non-backing-up ("rm") half of NEO's own logic -- wait, they ARE:
! archive() covers both backing-up (lines 1575-1616) and non-backing-up
! (lines 1625-1627) NEO branches. log_convergence() has no NEO override
! (the legacy `elseif (NEO)` arm at lines 1635-1637 is empty, just a
! comment "can we do something similar with NEO?") -- neo_t inherits
! flow_solver_t's no-op default there, which is exactly equivalent.
! supports_ale() also inherits the .false. default (NEO has no ALE step).

  use, intrinsic :: iso_fortran_env, only: int64
  use mod_kinds, only: wp, i4
  use mod_solver_iface, only: flow_solver_t, solver_run_ctx_t
  use mod_timer, only: timer_tic, timer_toc
  implicit none(type, external)
  private
  external :: na2vvvv, triangle2grd, neo2triangle

  public :: neo_t

  type, extends(flow_solver_t) :: neo_t
  contains
    procedure :: pre_mesh_setup => neo_pre_mesh_setup
    procedure :: prepare => neo_prepare
    procedure :: run => neo_run
    procedure :: harvest => neo_harvest
    procedure :: archive => neo_archive
  end type neo_t

contains

  subroutine neo_pre_mesh_setup(self, ctx)
    use mod_run_external, only: run_external
    class(neo_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%iter == 1 + ctx%nbegin .and. trim(ctx%testcase) == "ShockExpansion") then
      write (*, '(a)', advance='no') 'neogrid0               -->  '
      call timer_tic()
      execmd = ctx%bindir//'neogrid0'
      ifail = run_external(execmd, 'neogrid0')
      write (*, '(a)') ' ok'//timer_toc()
    end if
  end subroutine neo_pre_mesh_setup

  subroutine neo_prepare(self, ctx)
    class(neo_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    integer(int64) :: t0_na, t1_na, t0_tri, t1_tri, clock_rate
    real(wp) :: dt_na, dt_tri

    if (ctx%corrector) then

      call system_clock(count_rate=clock_rate)
      !$omp parallel sections
      !$omp section
      call system_clock(count=t0_na)
      call na2vvvv(trim(ctx%fname)//".1")
      call system_clock(count=t1_na)
      !$omp section
      call system_clock(count=t0_tri)
      call triangle2grd(trim(ctx%fname)//".1")
      call system_clock(count=t1_tri)
      !$omp end parallel sections
      dt_na = real(t1_na - t0_na, wp)/real(clock_rate, wp)
      dt_tri = real(t1_tri - t0_tri, wp)/real(clock_rate, wp)

      write (*, '(a,"  ",f11.4," s")') 'na2vvvv                -->   ok', dt_na
      write (*, '(a,"  ",f11.4," s")') 'triangle2grd           -->   ok', dt_tri

    else

!      neogrid0 works only for the 1st iteration of ShockVortex
!      but in all other cases we need this conversion
      if (.not. ctx%unsteady .or. ctx%iter /= 1 + ctx%nbegin .or.&
      &trim(ctx%testcase) == "ShockVortex") then

        call system_clock(count_rate=clock_rate)
        !$omp parallel sections
        !$omp section
        call system_clock(count=t0_na)
        call na2vvvv(trim(ctx%fname)//".1")
        call system_clock(count=t1_na)
        !$omp section
        call system_clock(count=t0_tri)
        call triangle2grd(trim(ctx%fname)//".1")
        call system_clock(count=t1_tri)
        !$omp end parallel sections
        dt_na = real(t1_na - t0_na, wp)/real(clock_rate, wp)
        dt_tri = real(t1_tri - t0_tri, wp)/real(clock_rate, wp)

        write (*, '(a,"  ",f11.4," s")') 'na00xTovvvv            -->   ok', dt_na
        write (*, '(a,"  ",f11.4," s")') 'triangle2grd           -->   ok', dt_tri

      else

        write (*, '(a)', advance='no') 'triangle2grd           -->   '
        call timer_tic()
        call triangle2grd(trim(ctx%fname)//".1")
        write (*, '(a)') 'ok'//timer_toc()

      end if

    end if
  end subroutine neo_prepare

  subroutine neo_run(self, ctx)
    use mod_run_external, only: run_external
    class(neo_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%corrector) then

      write (*, '(a)', advance='no') 'NEO                    -->  '
      call timer_tic()
      execmd = ctx%bindir//"CRD_euler"&
      &//" > log/neo.log"
      ifail = run_external(execmd, 'NEO')
      write (*, '(a)') ' ok'//timer_toc()

    else

!     ******************
      if (ctx%unsteady) then
!     ******************

!       if it is the 1st iteration
!       **************************
        if (ctx%iter == 1 + ctx%nbegin) then

          write (*, '(a)', advance='no') 'NEO 1st iteration      -->  '
          call timer_tic()

          execmd = ctx%bindir//"CRD_euler"//"> log/neo.log"
          ifail = run_external(execmd, 'NEO (1st iteration)')

          execmd = "cp ./NEO_data/output/vvvv.dat "//&
          &"./NEO_data/output/vvvv0.dat "
          ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

          execmd = "mv ./NEO_data/output/vvvv.dat "//&
          &"./NEO_data/output/vvvv_input.dat "
          ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

!         Here the following happens (for unsteady cases):
!         - the 1st iteration uses NEO_data/textinput/inputfile-exp.txt
!         - After the 1st NEO call, inputfile-exp.txt is moved in BAK
!         - Then, the inputfile-exp.txt in the testcase folder is
!           copied in NEO_data/textinput/inputfile-exp.txt
!         This is done because the two inputfile-exp.txt files differ
!         for the "Initial state" value. In the first case, it is 14
!         which means that the NEO function initial_solution()
!         writes the initial solution for the centered expansion test,
!         while after the 1st iteration should be 0, since we don't need
!         to initialize it again but instead just read_solution() which
!         happens if the variable intial_solution = 0.

          execmd = "mv ./NEO_data/textinput/inputfile-exp.txt "//&
          &"./NEO_data/textinput/inputfile-exp.txt.BAK "
          ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

          execmd = "cp inputfile-exp.txt "//"./NEO_data/textinput/"
          ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

          write (*, '(a)') ' ok'//timer_toc()

        end if ! 1ST ITERATION

      end if ! UNSTEADY

!     for all the other iterations
!     ****************************
      write (*, '(a)', advance='no') 'NEO                    -->   '
      call timer_tic()
      execmd = ctx%bindir//"CRD_euler"&
      &//"> log/neo.log"
      ifail = run_external(execmd, 'neo')
      write (*, '(a)') 'ok'//timer_toc()

    end if
  end subroutine neo_run

  subroutine neo_harvest(self, ctx)
    class(neo_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx

    if (ctx%corrector) then
      write (*, '(a)', advance='no') 'NEO2triangle           -->  '
      call timer_tic()
      call neo2triangle(trim(ctx%fname)//".1")
      write (*, '(a)') ' ok'//timer_toc()
    else
      write (*, '(a)', advance='no') 'NEO2triangle           -->   '
      call timer_tic()
      call neo2triangle(trim(ctx%fname)//".1")
      write (*, '(a)') 'ok'//timer_toc()
    end if
  end subroutine neo_harvest

  subroutine neo_archive(self, ctx)
    use mod_run_external, only: run_external
    class(neo_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%backing_up) then
      execmd = "mv -v shocknor.dat sh99.dat&
      &                                        "//ctx%fname//".* "//ctx%fnameback//".node "//&
      &ctx%backdir
      ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

      execmd =&
      &"cp -vp ./NEO_data/input/neogrid.grd&
      &                              ./NEO_data/input/vel.dat "//ctx%backdir
      ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

      execmd =&
      &"cp -v ./NEO_data/output/vvvv.dat "//ctx%backdir
      ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

      execmd =&
      &"mv -v ./NEO_data/output/vvvv_input.dat "//ctx%backdir
      ifail = run_external(execmd, 'mv', fatal_on_error=.false.)
    else
      execmd = "rm shocknor.dat "//ctx%fname//".* "//ctx%fnameback//&
      &".node "//"sh99.dat "
      ifail = run_external(execmd, 'rm', fatal_on_error=.false.)
    end if
  end subroutine neo_archive

end module mod_neo_solver
