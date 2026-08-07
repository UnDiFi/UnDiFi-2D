module mod_su2_solver
! Phase 3.7 increment 4 (ROADMAP.md #14, issue #15): su2_t, a mechanical
! port of the SU2-specific branch of main.f90's former
! if(EULFS)/elseif(SU2)/elseif(NEO) dispatch. Original main.f90 line
! numbers (pre-3.7, commit c5b9cfb): prepare() <- lines 897-904
! (triangle2su2), run() <- lines 912-918 (SU2_CFD), harvest() <- lines
! 927-934 (su22triangle) -- all three only for the predictor/always-run
! call site (site D in the plan). No periodic-BC support (see
! source_utils/triangle2su2/main.f), preserved as-is.
!
! **Real, pre-existing, user-confirmed-out-of-scope gaps, not covered by
! this port**: the legacy dispatch never had an `elseif (SU2)` arm at
! the UNSTEADY-corrector call site, nor at either half of the loop-end
! backup/cleanup archiving, nor at the convergenza.dat append -- SU2 +
! UNSTEADY has always silently skipped its entire corrector step, and
! SU2 runs have never gotten backup/cleanup archiving. su2_t reproduces
! this as explicit no-ops (prepare/run/harvest return immediately when
! ctx%corrector is true; archive/log_convergence/supports_ale are left
! at flow_solver_t's shared no-op/.false. defaults) rather than writing
! new SU2 corrector/archive logic -- the user's explicit scoping
! decision for this increment, since (a) this project's convention is
! to preserve pre-existing latent gaps during a mechanical port rather
! than opportunistically fix what the port happens to touch, and (b)
! SU2_CFD is not installed in this environment and zero regression-
! baseline coverage exists for any su2_* case (confirmed via
! tests/baseline/README.md), so a "fix" here would ship unverified.

  use mod_kinds, only: wp, i4
  use mod_solver_iface, only: flow_solver_t, solver_run_ctx_t
  use mod_timer, only: timer_tic, timer_toc
  implicit none(type, external)
  private

  public :: su2_t

  type, extends(flow_solver_t) :: su2_t
  contains
    procedure :: prepare => su2_prepare
    procedure :: run => su2_run
    procedure :: harvest => su2_harvest
  end type su2_t

contains

  subroutine su2_prepare(self, ctx)
    use mod_run_external, only: run_external
    class(su2_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%corrector) return ! see module header: no SU2 corrector support, preserved

! **********************************************************************
!  Convert the triangle files into SU2's native mesh + restart state:
!  echo na0x.1 / su2case | triangle2su2
!  "su2case" is a fixed basename (not fname) so su2case.cfg's
!  MESH_FILENAME/SOLUTION_FILENAME/RESTART_FILENAME never have to
!  change across outer iterations even though the Triangle basename
!  (fname) does. No periodic-BC support yet (see source_utils/
!  triangle2su2/main.f) -- fine for now, CircularCylinder has none.
! **********************************************************************

    write (*, '(a)', advance='no') 'triangle2su2           -->  '
    call timer_tic()
    execmd = "printf '"//ctx%fname&
    &//".1\nsu2case'|"&
    &//ctx%bindir//"triangle2su2-"//ctx%hostype&
    &//" > log/triangle2su2.log"
    ifail = run_external(execmd, 'triangle2su2')

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine su2_prepare

  subroutine su2_run(self, ctx)
    use mod_run_external, only: run_external
    class(su2_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%corrector) return ! see module header: no SU2 corrector support, preserved

! **************************
!  Run one step of SU2 code
! **************************
!  su2case.cfg sets ITER=1 with RESTART_SOL=YES: one implicit step
!  per outer UNDIFI iteration, the SU2 analogue of EulFS's -itmax 1.

    write (*, '(a)', advance='no') 'su2                    -->  '
    call timer_tic()
    execmd = ctx%bindir//"SU2_CFD"&
    &//" su2case.cfg > log/su2.log"

    ifail = run_external(execmd, 'su2')

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine su2_run

  subroutine su2_harvest(self, ctx)
    use mod_run_external, only: run_external
    class(su2_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%corrector) return ! see module header: no SU2 corrector support, preserved

! **********************************************************************
!  Convert su2case's restart state back into triangle fmt:
!  echo na0x.1 / su2case | su22triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.BAK
! **********************************************************************

    write (*, '(a)', advance='no') 'su22triangle           -->  '
    call timer_tic()
    execmd = "printf '"//ctx%fname&
    &//".1\nsu2case'|"&
    &//ctx%bindir//"su22triangle-"//ctx%hostype&
    &//" > log/su22triangle.log"
    ifail = run_external(execmd, 'su22triangle')

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine su2_harvest

end module mod_su2_solver
