module mod_eulfs_solver
! Phase 3.7 increment 2 (ROADMAP.md #14, issue #15): eulfs_t, a
! mechanical port of every EulFS-specific branch of main.f90's former
! if(EULFS)/elseif(SU2)/elseif(NEO) dispatch. Every literal command
! string, flag, and comment is preserved exactly -- this is a
! relocation into type-bound procedures, not a rewrite. Original
! main.f90 line numbers (pre-3.7, commit c5b9cfb) cited below so the
! port can be checked directly against git history.
!
! setup()           <- lines 209-223 (once, before the outer loop:
!                       reset convergenza.dat, prime .petscrc for
!                       STEADY runs). Preserves a real pre-existing bug,
!                       confirmed via `git log -L209,222:source/main.f90`
!                       to predate even this repo's Phase 1.7
!                       run_external migration: the `.petscrc`-priming
!                       cp command is built into execmd but its
!                       `ifail = system(execmd)` call was ALREADY
!                       commented out in the legacy source, so it is
!                       built but never executed for STEADY EulFS runs.
!                       The `ifail` checked afterward is therefore still
!                       the preceding `rm`'s status. Not fixed here.
! prepare()         <- lines 789-833 / 1100-1144 (triangle2dat; the two
!                       original call sites are byte-identical, so this
!                       is a single procedure, not one per stage)
! run()             <- lines 839-863 (predictor) / 1150-1174 (corrector)
!                       merged via ctx%corrector: petscrc file selection
!                       differs (_predictor vs _corrector) and the
!                       predictor arm gates its petscrc-copy and
!                       step000001->file001 copy behind `if (unsteady)`
!                       while the corrector arm (only ever reached when
!                       unsteady is already true) does both
!                       unconditionally -- preserved via `ctx%corrector
!                       .or. ctx%unsteady`, equivalent to the legacy
!                       code's own commented-out-but-redundant
!                       `if (UNSTEADY)` wrapper there.
! harvest()         <- lines 872-881 / 1183-1192 (dat2triangle; also
!                       byte-identical between the two call sites)
! archive()         <- lines 1570-1574 (backing-up branch) + 1621-1624
!                       (non-backing-up/cleanup branch)
! log_convergence() <- lines 1632-1634 (convergenza.dat append)
! supports_ale()    <- .true., replacing the `if (EULFS)` guards around
!                       the ALE solzne call at lines 683-692 and inside
!                       mod_shock_advance.f90's advance()

  use mod_kinds, only: wp, i4
  use mod_solver_iface, only: flow_solver_t, solver_run_ctx_t
  use mod_timer, only: timer_tic, timer_toc
  implicit none(type, external)
  private

  public :: eulfs_t

  type, extends(flow_solver_t) :: eulfs_t
  contains
    procedure :: setup => eulfs_setup
    procedure :: prepare => eulfs_prepare
    procedure :: run => eulfs_run
    procedure :: harvest => eulfs_harvest
    procedure :: archive => eulfs_archive
    procedure :: log_convergence => eulfs_log_convergence
    procedure :: supports_ale => eulfs_supports_ale
  end type eulfs_t

contains

  subroutine eulfs_setup(self, unsteady)
    use mod_run_external, only: run_external
    use mod_error, only: fatal
    class(eulfs_t), intent(inout) :: self
    logical, intent(in) :: unsteady
    character(len=255) :: execmd
    integer(i4) :: ifail

    execmd = "rm -fv convergenza.dat"
    ifail = run_external(execmd, 'rm')

!        copy file .petsrc in home
!        for UNSTEADY EulFS simulations this file
!        will be overwritten with .petsrc_predictor and
!        .petsrc_corrector in their respective steps
!        They only differ in the dt value
    if (.not. unsteady) then
      execmd = "cp -fv .petscrc .petscrc"
    end if
!        ifail = system(execmd)
    if (ifail .ne. 0) call fatal('system command failed', ifail)
  end subroutine eulfs_setup

  subroutine eulfs_prepare(self, ctx)
    use mod_run_external, only: run_external
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(eulfs_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    include 'paramt.h'
    character(len=255) :: execmd
    character(len=2) :: color1, color2
    integer(i4) :: ifail

    write (*, '(a)', advance='no') 'triangle2dat           -->  '
    call timer_tic()
    if (nprdbnd .eq. 0) then                            ! for the cases without periodic BCs
      execmd = "printf '"//ctx%fname&
      &//".1\nn'|"&
      &//ctx%bindir//"triangle2dat-NEW-"//ctx%hostype&
      &//" > log/triangle2dat.log"
    elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 1) then    ! for the cases with only one periodic boundary
      write (color1, fmt="(i2.2)") prdbndclr(1, 1)          ! with points having the same x
      write (color2, fmt="(i2.2)") prdbndclr(2, 1)
      execmd = "printf '"//ctx%fname&
      &//".1\ny\n"&
      &//color1//"\n"&
      &//color2//"\nx'|"&
      &//ctx%bindir//"triangle2dat-NEW-"//ctx%hostype&
      &//" > log/triangle2dat.log"
    elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 2) then ! for the cases with only one periodic boundary
      write (color1, fmt="(i2.2)") prdbndclr(1, 1)           ! with points having the same y
      write (color2, fmt="(i2.2)") prdbndclr(2, 1)
      execmd = "printf '"//ctx%fname&
      &//".1\ny\n"&
      &//color1//"\n"&
      &//color2//"\ny'|"&
      &//ctx%bindir//"triangle2dat-NEW-"//ctx%hostype&
      &//" > log/triangle2dat.log"
    else ! for cases with more thatn one periodic boundary
      write (*, *) ' case not implemented!'
    end if

    ifail = run_external(execmd, 'triangle2dat')

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine eulfs_prepare

  subroutine eulfs_run(self, ctx)
    use mod_run_external, only: run_external
    class(eulfs_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%corrector) then
!          It runs the corrector step of the EulFS code (now we use the full dt)
      execmd = "cp -f .petscrc_corrector .petscrc"
      ifail = run_external(execmd, 'cp')
    elseif (ctx%unsteady) then
!          It runs the predictor step of the EulFS code (we need to use dt/2)
      execmd = "cp -f .petscrc_predictor .petscrc"
      ifail = run_external(execmd, 'cp')
    end if

    write (*, '(a)', advance='no') 'eulfs                  -->  '
    call timer_tic()
    execmd = ctx%bindir//"EulFS_"//ctx%hostype&
    &//" -itmax 1 > log/eulfs.log"

    ifail = run_external(execmd, 'eulfs')

    if (ctx%corrector .or. ctx%unsteady) then
      execmd = "cp step000001.dat file001.dat"
      ifail = run_external(execmd, 'cp')
    end if

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine eulfs_run

  subroutine eulfs_harvest(self, ctx)
    use mod_run_external, only: run_external
    class(eulfs_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    write (*, '(a)', advance='no') 'dat2triangle           -->  '
    call timer_tic()
    execmd = "printf '"//ctx%fname//".1' | "//ctx%bindir&
    &//"dat2triangle-NEW-"//ctx%hostype&
    &//">log/dat2triangle.log"
    ifail = run_external(execmd, 'dat2triangle')

    write (*, '(a)') ' ok'//timer_toc()
  end subroutine eulfs_harvest

  subroutine eulfs_archive(self, ctx)
    use mod_run_external, only: run_external
    class(eulfs_t), intent(inout) :: self
    type(solver_run_ctx_t), intent(in) :: ctx
    character(len=255) :: execmd
    integer(i4) :: ifail

    if (ctx%backing_up) then
      execmd = "mv -v shocknor.dat file00[1-4].dat file010.dat sh&
      &99.dat  "//ctx%fname//".* "//ctx%fnameback//".node&
      &                               "//ctx%backdir
      ifail = run_external(execmd, 'mv', fatal_on_error=.false.)
    else
      execmd = "rm shocknor.dat file00[1-4].dat file010.dat&
      &                               "//ctx%fname//".* "//ctx%fnameback//".node "//&
      &"sh99.dat "
      ifail = run_external(execmd, 'rm', fatal_on_error=.false.)
    end if
  end subroutine eulfs_archive

  subroutine eulfs_log_convergence(self)
    use mod_run_external, only: run_external
    class(eulfs_t), intent(inout) :: self
    character(len=255) :: execmd
    integer(i4) :: ifail

    execmd = "cut -c34- convhst.l2 >> convergenza.dat"
    ifail = run_external(execmd, 'cut', fatal_on_error=.false.)
  end subroutine eulfs_log_convergence

  logical function eulfs_supports_ale(self) result(ale)
    class(eulfs_t), intent(in) :: self
    ale = .true.
  end function eulfs_supports_ale

end module mod_eulfs_solver
