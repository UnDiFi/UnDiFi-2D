module mod_shock_advance
! Phase 3.6 (ROADMAP.md #14): collapses the two identically-shaped
! "absorb the just-solved mesh into the shock/discontinuity state"
! sequences that used to be duplicated inline in main.f90's time loop --
! once (mid-loop, UNSTEADY only) targeting the *new* predictor/corrector
! shadow arrays owned by mod_shock_system, and once (loop-end, always)
! targeting the *old* ones. The solver dispatch (if EULFS/SU2/NEO) is
! deliberately NOT folded in here -- ROADMAP.md assigns that three-way
! if/elseif to 3.7 (flow_solver_t, issue #15) separately, and it isn't
! duplicated in the same shape anyway (the corrector-only call site is
! missing SU2 support and several STEADY-only branches, so it was never
! a clean copy of the predictor call site to begin with).

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, nshmax, npshmax, nspmax
  use mod_mesh, only: mesh_t
  use mod_shock_system, only: xysh, zroeshuold, zroeshdold, norsh, wsh,&
  &xyshnew, norshnew, wshnew, zroeshuoldnew, zroeshdoldnew
  use mod_timer, only: timer_tic, timer_toc
  implicit none(type, external)
  private

  public :: advance

contains

  ! Reads the triangle mesh just produced (fname, already suffixed by the
  ! caller), absorbs its nodal state into bkg, re-fixes the R-H/special-
  ! point state at the shock points via co_state_dps/fx_state_dps, and
  ! copies the fixed shock state back into fit. new_shadow selects which
  ! of mod_shock_system's two shadow-array sets ("new"/corrector-feed vs
  ! "old"/next-iteration) this call targets -- see mod_shock_system.f90's
  ! header for why only some fields have a real independent "new" copy.
  ! When dt is present, also runs calc_vel (+ solzne if eulfs) to prepare
  ! the grid velocity for the next solve -- only the mid-loop/new_shadow
  ! call site needs this; the loop-end/old_shadow call site omits dt.
  subroutine advance(bkg, fit, fname, iter, nshocks, nshockpoints, nshocksegs,&
  &typeshocks, nspecpoints, typespecpoints, shinspps, ispclr,&
  &new_shadow, eulfs, varray, velfile, mode, testcase, dt, velflag)
    type(mesh_t), intent(inout) :: bkg
    type(mesh_t), intent(inout) :: fit
    character(len=*), intent(in) :: fname
    integer(i4), intent(in) :: iter, nshocks
    integer(i4), intent(in) :: nshockpoints(nshmax), nshocksegs(nshmax)
    character(len=1), intent(in) :: typeshocks(nshmax)
    integer(i4), intent(in) :: nspecpoints
    character(len=5), intent(in) :: typespecpoints(nspmax)
    integer(i4), intent(in) :: shinspps(2, 5, nspmax), ispclr(5, nspmax)
    logical, intent(in) :: new_shadow, eulfs
    real(wp), intent(inout) :: varray(ndim, 30000)
    character(len=*), intent(in) :: velfile, testcase
    character(len=1), intent(in) :: mode
    real(wp), intent(in), optional :: dt
    character(len=1), intent(in), optional :: velflag

    external readmesh, co_state_dps, fx_state_dps, dcopy, calc_vel, solzne

    integer(i4) :: totshockpoints
    logical :: fndbnds
    real(wp) :: nowtime

    write (*, '(a)', advance='no') 'readmesh               -->  '
    call timer_tic()
    fndbnds = .false.
    call readmesh(fit, fname, fndbnds)
    write (*, '(a)') ' ok'//timer_toc()

    totshockpoints = 2*nshmax*npshmax

    write (*, '(a)', advance='no') 'zroe(1)->zroe(0)       -->  '
    call timer_tic()
    if (fit%npoin .eq. (bkg%npoin + totshockpoints)) then
      call dcopy(ndof*fit%npoin, fit%zroe, 1, bkg%zroe, 1)
      write (*, '(a)') ' ok'//timer_toc()
    else
!         the nof gridpoints in grid(1) must equal the number of
!         gridpoints on the background mesh + 2 * nshockpoints
      write (6, *) 'there is a mismatch in the nof gridpoints'
      write (6, *) 'btw grid(0) and grid(1)'
      write (*, *) bkg%npoin, totshockpoints
      write (*, *) fit%npoin, totshockpoints
      error stop 1
    end if

    if (new_shadow) then

      write (*, '(a)', advance='no') 'co_state_dps           -->  '
      call timer_tic()
      call co_state_dps(&
      &xyshnew,&
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuoldnew,&
      &zroeshdoldnew,&
      &norshnew,&
      &wshnew,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &iter)
      write (*, '(a)') ' ok'//timer_toc()

      write (*, '(a)', advance='no') 'fx_state_dps           -->  '
      call timer_tic()
      call fx_state_dps(&
      &xyshnew,&                                           ! not used
      &bkg%xy(1, bkg%npoin + 1),&                     ! upstream   coord.
      &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream   state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuoldnew,&
      &zroeshdoldnew,&
      &norshnew,&
      &wshnew,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &iter,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps,&
      &ispclr,&
      &bkg%ia,&
      &bkg%ja,&
      &bkg%iclr,&
      &bkg%nclr,&
      &bkg%xy)
      write (*, '(a)') ' ok'//timer_toc()

    else

      write (*, '(a)', advance='no') 'co_state_dps           -->  '
      call timer_tic()
      call co_state_dps(&
      &xysh,&
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuold,&
      &zroeshdold,&
      &norsh,&
      &wsh,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &iter)
      write (*, '(a)') ' ok'//timer_toc()

      write (*, '(a)', advance='no') 'fx_state_dps           -->  '
      call timer_tic()
      call fx_state_dps(&
      &xysh,&
      &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
      &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuold,&
      &zroeshdold,&
      &norsh,&
      &wsh,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &iter,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps,&
      &ispclr,&
      &bkg%ia,&
      &bkg%ja,&
      &bkg%iclr,&
      &bkg%nclr,&
      &bkg%xy)
      write (*, '(a)') ' ok'//timer_toc()

    end if

    write (*, '(a)', advance='no') 'zroesh(0)->zroesh(1)   -->  '
    call timer_tic()
    call dcopy(ndof*totshockpoints,&
    &bkg%zroe(1, bkg%npoin + 1), 1,&
    &fit%zroe(1, bkg%npoin + 1), 1)
    write (*, '(a)') ' ok'//timer_toc()

    if (present(dt)) then

      write (*, '(a)', advance='no') 'calc_vel               -->  '
      call timer_tic()
      call calc_vel(&
      &bkg%npoin,&
      &varray,&
      &dt,&
      &bkg%xy,&
      &wsh,& ! matches the original call site's own unresolved "!WSHnew?" question -- not fixed here
      &iter,&
      &velflag,&
      &nowtime,&
      &testcase)
      write (*, '(a)') ' ok'//timer_toc()

      if (eulfs) then
        write (*, '(a)', advance='no') 'solzne                 -->   '
        call timer_tic()
        call solzne(&
        &velfile,&
        &varray,&
        &ndim,&
        &bkg%npoin + 2*npshmax*nshmax,&
        &mode)
        write (*, '(a)') ' ok'//timer_toc()
      end if

    end if

  end subroutine advance

end module mod_shock_advance
