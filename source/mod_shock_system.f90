module mod_shock_system
! Phase 3.1 (ROADMAP.md #14): owns the NSHMAX-sized flat shock/
! discontinuity arrays that used to be main.f90 locals, plus the
! disc(:)/disc_new(:) discontinuity_t arrays that alias them (see
! mod_discontinuity.f90 for why they're pointer-based, not independent
! allocatables).
!
! The legacy dispatch routines (co_norm, fx_state_dps, fx_msh_sps,
! co_pnt_dspl, fx_dps_loc, co_state_dps, rd_dps, rd_dps_eq, mv_dps,
! wtri, wtri0, ...) are UNCHANGED in this phase and still address these
! arrays directly by name with an internal do ish=1,nshocks loop --
! main.f90's call sites therefore need no edits beyond a `use` statement,
! since these public module variables keep the exact names main.f90 used
! before.
!
! disc(:)/disc_new(:) mirror the bkg/fit/bak mesh_t precedent (two named
! instances of one type for the predictor/corrector shadow state), but
! faithfully reflect that the legacy code does NOT actually duplicate
! everything: there is no separate "new" R-H state (both predictor and
! corrector steps read/write the same bkg%zroe tail), so disc_new's
! zu/zd/nodcod simply alias disc's -- only xy/nor/w/zu_old/zd_old have
! real, independently-updated predictor/corrector shadow copies today
! (xyshnew/norshnew/wshnew/zroeshuoldnew/zroeshdoldnew).

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax, nshmax
  use mod_discontinuity, only: discontinuity_t
  implicit none(type, external)
  private

  public :: xysh, zroeshuold, zroeshdold, norsh, wsh, nodcodsh
  public :: xyshnew, norshnew, wshnew, wshmean, zroeshuoldnew, zroeshdoldnew
  public :: disc, disc_new
  public :: shock_system_init, shock_system_refresh

  real(wp), target, save :: xysh(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: zroeshuold(ndof, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: zroeshdold(ndof, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: norsh(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: wsh(ndim, npshmax, nshmax) = 0.0_wp

! arrays for unsteady predictor-corrector time accurate integration
  real(wp), target, save :: xyshnew(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: norshnew(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: wshnew(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: wshmean(ndim, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: zroeshuoldnew(ndof, npshmax, nshmax) = 0.0_wp
  real(wp), target, save :: zroeshdoldnew(ndof, npshmax, nshmax) = 0.0_wp

  integer(i4), target, save :: nodcodsh(npshmax, nshmax) = 0_i4

  type(discontinuity_t), allocatable :: disc(:), disc_new(:)

contains

  ! Allocates disc(:)/disc_new(:) to the run's actual shock count and
  ! points every field at its backing storage. zroe must be the caller's
  ! TARGET mesh_t%zroe component (main.f90 passes bkg%zroe) -- disc(ish)%zu/
  ! %zd alias its padded tail exactly as the legacy sequence-association
  ! call sites do (npoin+1 upstream, npoin+1+nshmax*npshmax downstream),
  ! per shock ish's npshmax-wide block within that tail.
  subroutine shock_system_init(nshocks, npoin, zroe)
    integer(i4), intent(in) :: nshocks, npoin
    real(wp), target, intent(in) :: zroe(:, :)
    integer(i4) :: ish, ustart, dstart

    allocate (disc(nshocks), disc_new(nshocks))
    do ish = 1, nshocks
      ustart = npoin + (ish - 1)*npshmax + 1
      dstart = npoin + nshmax*npshmax + (ish - 1)*npshmax + 1

      disc(ish)%xy => xysh(:, :, ish)
      disc(ish)%zu => zroe(:, ustart:ustart + npshmax - 1)
      disc(ish)%zd => zroe(:, dstart:dstart + npshmax - 1)
      disc(ish)%zu_old => zroeshuold(:, :, ish)
      disc(ish)%zd_old => zroeshdold(:, :, ish)
      disc(ish)%nor => norsh(:, :, ish)
      disc(ish)%w => wsh(:, :, ish)
      disc(ish)%nodcod => nodcodsh(:, ish)

      disc_new(ish)%xy => xyshnew(:, :, ish)
      disc_new(ish)%nor => norshnew(:, :, ish)
      disc_new(ish)%w => wshnew(:, :, ish)
      disc_new(ish)%zu_old => zroeshuoldnew(:, :, ish)
      disc_new(ish)%zd_old => zroeshdoldnew(:, :, ish)
      ! no separate "new" R-H state or nodcod exists in the legacy code
      disc_new(ish)%zu => disc(ish)%zu
      disc_new(ish)%zd => disc(ish)%zd
      disc_new(ish)%nodcod => disc(ish)%nodcod
    end do
  end subroutine shock_system_init

  ! disc(:)/disc_new(:)'s scalar bookkeeping (npoints/nedges/kind) has no
  ! storage of its own -- it's copied in from the legacy nshockpoints/
  ! nshocksegs/typeshocks arrays, which rd_dps/rd_dps_eq mutate directly.
  ! Call after any point where those may have changed.
  subroutine shock_system_refresh(nshockpoints, nshocksegs, typeshocks)
    integer(i4), intent(in) :: nshockpoints(:), nshocksegs(:)
    character(len=1), intent(in) :: typeshocks(:)
    integer(i4) :: ish

    do ish = 1, size(disc)
      disc(ish)%npoints = nshockpoints(ish)
      disc(ish)%nedges = nshocksegs(ish)
      disc(ish)%kind = typeshocks(ish)
      disc_new(ish)%npoints = nshockpoints(ish)
      disc_new(ish)%nedges = nshocksegs(ish)
      disc_new(ish)%kind = typeshocks(ish)
    end do
  end subroutine shock_system_refresh

end module mod_shock_system
