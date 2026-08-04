module mod_discontinuity
! Phase 3.1 (ROADMAP.md #14): discontinuity_t, the named replacement for
! main.f90's NSHMAX-sized shock/discontinuity arrays. Array components
! are POINTER, not allocatable (unlike mod_mesh.f90's mesh_t): the legacy
! dispatch routines (co_norm, fx_state_dps, rd_dps, mv_dps, ...) still
! expect one contiguous (dim,npshmax,nshmax)-shaped block per array, with
! an internal do ish=1,nshocks loop -- independently heap-allocated
! per-shock arrays would break that sequence-association calling
! convention. Each disc(ish)/disc_new(ish) instead aliases a slice of
! mod_shock_system's flat backing arrays (see there), so every legacy
! call site is unchanged in this phase.
!
! Consequently each field's usable extent is still capped at npshmax
! (the pointer's declared shape is the *capacity*, not the *count*) --
! %npoints/%nedges record how many of those columns are meaningful,
! exactly like today's nshockpoints/nshocksegs. True unbounded growth
! (dropping the npshmax/nshmax ceilings) needs the legacy routines
! rewritten to operate on one disc(ish) at a time -- Phase 3.2/3.3.
!
! disc(:) has no consumer yet in this phase -- nothing in 3.1 reads or
! writes through it. It exists so Phase 3.2 has real, correctly-aliased
! storage to start binding type-bound procedures to.

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: discontinuity_t

  type :: discontinuity_t
    character(len=:), allocatable :: kind ! 's'=shock, 'd'=contact today (typeshocks); Part III widens this to 'shock'/'contact'/'slip'/'interface'
    logical     :: closed = .false. ! unused until Part III (F11, closed/bubble interfaces)
    integer(i4) :: npoints = 0, nedges = 0
    integer(i4) :: gas_up = 1, gas_dn = 1 ! unused until Part III (F11, multi-species)
    real(wp), pointer :: xy(:, :) => null()     ! (ndim, npshmax)  point coords, was xysh
    real(wp), pointer :: zu(:, :) => null()     ! (ndof, npshmax)  upstream R-H state, was the bkg%zroe tail slice
    real(wp), pointer :: zd(:, :) => null()     ! (ndof, npshmax)  downstream R-H state, was the bkg%zroe tail slice
    real(wp), pointer :: zu_old(:, :) => null() ! (ndof, npshmax)  upstream state, previous iteration, was zroeshuold
    real(wp), pointer :: zd_old(:, :) => null() ! (ndof, npshmax)  downstream state, previous iteration, was zroeshdold
    real(wp), pointer :: nor(:, :) => null()    ! (ndim, npshmax)  unit normal, was norsh/vshnor
    real(wp), pointer :: w(:, :) => null()      ! (ndim, npshmax)  discontinuity speed, was wsh
    integer(i4), pointer :: nodcod(:) => null() ! (npshmax)        boundary/interior point code, was nodcodsh
  end type discontinuity_t

end module mod_discontinuity
