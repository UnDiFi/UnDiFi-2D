module mod_special_point
! Phase 3.2 (merged 3.2+3.3+3.4, ROADMAP.md #14): type hierarchy for the
! special-point object model that replaces all six typespecpoints-keyed
! dispatch chains (fx_state_dps, co_pnt_dspl, fx_msh_sps, fx_dps_loc,
! co_norm, and the re_sdw_info/wrt_sdw_info serialization pair). See
! /home/lk/.claude/plans/proud-jingling-donut.md's "Merged 3.2" section
! for the full design (type-hierarchy rationale, per-code nshe table,
! which of the five per-chain behaviors each type overrides vs. leaves
! at the shared no-op default, and the six known latent bugs that must
! be preserved bit-for-bit when each chain is ported).
!
! Increment 4 (this revision) ports fx_state_dps.f90's real logic into
! every type's solve_state, replacing the increment-2 stubs. The other
! four behaviors (remesh_boundary/displace/relocate/correct_normal)
! still default to the shared no-op; increments 5-8 add real overrides
! there only where the corresponding legacy chain does real work.
!
! special_point_t's array components are plain allocatable, NOT
! pointer-into-shared-storage like mod_discontinuity.f90's
! discontinuity_t. discontinuity_t needs pointers because the legacy
! shock-curve routines (co_norm, fx_state_dps, rd_dps, ...) still expect
! one contiguous (dim,npshmax,nshmax) block with an internal
! do ish=1,nshocks loop. special_point_t is the opposite: every one of
! the six dispatch chains already loops do isppnts=1,nspecpoints and
! builds small, transient, per-point locals (xtpi(20), xtp(24), ...)
! fresh from shinspps/ispclr before doing anything -- so a
! special_point_t built fresh via the registry, populated via
! mod_special_point_registry's sp_unpack, used, and discarded, matches
! the existing call pattern with zero pointer-bridge complexity.
!
! solve_state's abstract interface is deliberately broad (it takes
! essentially every array fx_state_dps.f90 itself receives as a dummy
! argument) because the nine concrete types' real needs are genuinely
! heterogeneous -- regular_reflection_t's curved variant alone needs the
! boundary-mesh arrays (ia/ja/iclr/nclr/corg) that the other eight never
! touch. Fortran requires every override of a deferred binding to share
! one exact signature, so the interface is sized for the union of all
! nine types' needs; each type's implementation simply ignores the
! arguments it doesn't need (matching how e.g. sp_noop already ignores
! its sole argument for the sixteen no-op codes elsewhere).

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim
  implicit none(type, external)
  private
  public :: special_point_t
  public :: triple_point_t, quad_point_t, trailing_edge_t
  public :: regular_reflection_t, wall_float_t, floating_wall_point_t
  public :: end_point_t, start_point_t, connection_t
  public :: new_triple_point_t, new_quad_point_t, new_trailing_edge_t
  public :: new_regular_reflection_t, new_wall_float_t
  public :: new_floating_wall_point_t, new_end_point_t, new_start_point_t
  public :: new_connection_t

  type, abstract :: special_point_t
    character(len=5) :: code = '' ! 'TP','QP',... for logging/serialization round-trip
    integer(i4) :: nshe = 0 ! number of (shock,leg) records this point references
    integer(i4), allocatable :: ish(:) ! (nshe) shock indices, from shinspps(1,:,isppnts)
    integer(i4), allocatable :: leg(:) ! (nshe) endpoint flags,  from shinspps(2,:,isppnts)
  contains
    procedure(sp_solve_if), deferred :: solve_state ! fx_state_dps.f90 -- every type overrides (increment 4)
    procedure :: remesh_boundary => sp_noop ! fx_msh_sps.f90  -- default no-op (increment 6)
    procedure :: displace => sp_displace_noop ! co_pnt_dspl.f90 -- most types override (this increment)
    procedure :: relocate => sp_noop ! fx_dps_loc.f90  -- default no-op (increment 7)
    procedure :: correct_normal => sp_noop ! co_norm.f90     -- default no-op (increment 8)
  end type special_point_t

  abstract interface
    subroutine sp_solve_if(this, xysh, zroeshu, zroeshd, zroeshuold,&
    &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
    &ispclr_flat)
      import :: special_point_t, wp, i4, ndim
      class(special_point_t), intent(inout) :: this
      real(wp), intent(inout) :: xysh(:, :, :)
      real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
      real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
      real(wp), intent(in) :: vshnor(:, :, :)
      real(wp), intent(inout) :: wsh(:, :, :)
      integer(i4), intent(in) :: nshockpoints(:)
      integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
      real(wp), intent(in) :: corg(ndim, *)
      integer(i4), intent(in) :: ispclr_flat ! = ispclr(isppnts) via the legacy rank-1 mis-indexed read (RR only, see rr_solve_state)
    end subroutine sp_solve_if

! dx_carry/dy_carry reproduce a confirmed pre-existing latent bug in
! co_pnt_dspl.f90 (bug #2 in the merged-3.2 plan): its RRX branch reads
! a `dx` local that the ENCLOSING dispatch subroutine never resets per
! point -- it holds whatever wall_float_t's x/y branch (or, before any
! such branch has run, the generic per-point offset loop's very last
! iteration) last left it as. That shared-mutable-local dependency has
! no equivalent once each point's displace is an independent call, so
! the caller (co_pnt_dspl.f90) threads its own persistent dx/dy locals
! through every displace call via these two arguments -- wf_displace
! updates them (mirroring the legacy code's own assignment into the
! same shared local), rr_displace's RRX path reads (but does not
! update) them, matching the original bit-for-bit. Every other type
! ignores them.
    subroutine sp_displace_if(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
    &nshockpoints, dx_carry, dy_carry)
      import :: special_point_t, wp, i4
      class(special_point_t), intent(inout) :: this
      real(wp), intent(in) :: xysh(:, :, :)
      real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
      real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
      integer(i4), intent(in) :: nshockpoints(:)
      real(wp), intent(inout) :: dx_carry, dy_carry
    end subroutine sp_displace_if
  end interface

! TP -- triple point (internal), nshe=4. Newton-solved via co_utp today.
  type, extends(special_point_t) :: triple_point_t
  contains
    procedure :: solve_state => tp_solve_state
    procedure :: displace => tp_displace
  end type triple_point_t

! QP -- quadruple point (internal), nshe=5. Newton-solved via co_uqp today.
  type, extends(special_point_t) :: quad_point_t
  contains
    procedure :: solve_state => qp_solve_state
    procedure :: displace => qp_displace
  end type quad_point_t

! TE -- trailing edge point, nshe=3. Also Newton-solved via co_uqp today
! (fx_state_dps.f90 reuses the quadruple-point solver for TE, with two
! of the five states defaulted to placeholder values).
  type, extends(special_point_t) :: trailing_edge_t
  contains
    procedure :: solve_state => te_solve_state
    procedure :: displace => te_displace
  end type trailing_edge_t

! RRX (curved=.false.) / RR (curved=.true.) -- regular reflection off a
! wall, nshe=2. Newton-solved via co_urr today. iclr (boundary color per
! leg) is only meaningful when curved -- RR recomputes its wall tangent
! from boundary geometry (reading one ispclr(k,isppnts) per of its 2
! legs, confirmed in re_sdw_info.f90/wrt_sdw_info.f90), RRX assumes an
! axis-aligned wall and carries no color at all.
  type, extends(special_point_t) :: regular_reflection_t
    logical :: curved = .false.
    integer(i4) :: iclr(2) = 0
  contains
    procedure :: solve_state => rr_solve_state
    procedure :: displace => rr_displace
  end type regular_reflection_t

! WPNRX/WPNRY/IPX/IPY/OPX/OPY -- floating points constrained to slide
! along an axis-aligned boundary, nshe=1. restore_only distinguishes
! IP*'s "restore old state + zero wsh" behavior from OP*/WPN*'s
! "re-project onto axis" behavior (byte-identical between OP* and WPN*
! in fx_state_dps.f90, confirmed during design) -- both behaviors are
! axis-generic, keyed only by `axis`.
  type, extends(special_point_t) :: wall_float_t
    character(len=1) :: axis = 'X'
    logical :: restore_only = .false.
  contains
    procedure :: solve_state => wf_solve_state
    procedure :: displace => wf_displace
  end type wall_float_t

! FWP -- floating wall point on a coloured (possibly curved) boundary,
! nshe=1. Kept separate from wall_float_t: it alone carries a boundary
! color and has real (non-no-op) work in 3 of the 5 behavior chains.
  type, extends(special_point_t) :: floating_wall_point_t
    integer(i4) :: iclr = 0
  contains
    procedure :: solve_state => fwp_solve_state
  end type floating_wall_point_t

! EP -- end point, nshe=1. Kept distinct from start_point_t despite
! byte-identical co_pnt_dspl.f90 displace bodies today -- EP and SP
! diverge for real in fx_state_dps.f90 and co_norm.f90 (see plan).
  type, extends(special_point_t) :: end_point_t
  contains
    procedure :: solve_state => ep_solve_state
    procedure :: displace => ep_displace
  end type end_point_t

! SP -- start/sonic point (characteristic coalescence), nshe=1.
  type, extends(special_point_t) :: start_point_t
  contains
    procedure :: solve_state => sonic_solve_state
    procedure :: displace => sonic_displace
  end type start_point_t

! C (periodic=.false.) / PC (periodic=.true.) -- connection between two
! shocks, nshe=2. iclr is only meaningful when periodic.
  type, extends(special_point_t) :: connection_t
    logical :: periodic = .false.
    integer(i4) :: iclr(2) = 0
  contains
    procedure :: solve_state => conn_solve_state
    procedure :: displace => conn_displace
  end type connection_t

contains

! TP: triple point. Ports fx_state_dps.f90's 'TP' branch verbatim
! (lines 176-488 of the pre-increment-4 file): recovers the 4 zones'
! states from zroeshu/zroeshd(old), the 3 incident shock-normal slopes,
! seeds the triple-point speed at -0.01, calls co_utp with its
! sign-flip-on-IFAIL retry (stop on double failure, matching the
! original exactly), then scatters the converged 4-zone state back into
! zroeshu/zroeshd and propagates the triple-point velocity to all 4
! legs. z1v/z2v/z3v/z4v were mod_freestream module scratch in the
! original (shared across every branch of the old if/elseif) -- made
! local here since they're pure write-then-immediately-consumed scratch
! within this one branch, never read by anything else while this runs;
! this changes no computed value, only removes an unnecessary shared-
! mutable-global dependency.
  subroutine tp_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(triple_point_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat
    external co_utp
    include 'paramt.h'

    real(wp) :: xtpi(20), xtp(20)
    real(wp) :: f1, f2, help, dx, dy, dxr14, dyr14, r14, r23, r14old, wws
    real(wp) :: kine, hh, taux12, tauy12, z1v, z2v, z3v, z4v
    integer(i4) :: i, ish1, ish2, ish3, ish4, ip1, ip2, ip3, ip4
    logical :: flag1, ifail

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = this%ish(3); i = this%leg(3) - 1; ip3 = 1 + i*(nshockpoints(ish3) - 1)
    ish4 = this%ish(4); i = this%leg(4) - 1; ip4 = 1 + i*(nshockpoints(ish4) - 1)

! determine family of shock 1
    f1 = vshnor(1, ip1, ish1)*zroeshu(4, ip1, ish1) -&
    &vshnor(2, ip1, ish1)*zroeshu(3, ip1, ish1)
    f1 = -sign(1.d0, f1)

! determine family of shock 2
    f2 = vshnor(1, ip2, ish2)*zroeshu(4, ip2, ish2) -&
    &vshnor(2, ip2, ish2)*zroeshu(3, ip2, ish2)
    f2 = -sign(1.d0, f2)

! state 1: upstream of incident shock if opposite families, else downstream
    if (f1*f2 .lt. 0.) then
      xtpi(4) = zroeshu(4, ip1, ish1)/zroeshu(1, ip1, ish1)
      xtpi(3) = zroeshu(3, ip1, ish1)/zroeshu(1, ip1, ish1)
      xtpi(1) = zroeshu(1, ip1, ish1)*zroeshu(1, ip1, ish1)
      help = zroeshu(3, ip1, ish1)**2 + zroeshu(4, ip1, ish1)**2
      xtpi(2) = gm1/ga*(zroeshu(1, ip1, ish1)*zroeshu(2, ip1, ish1)&
      &- 0.5d0*help)
    else
      xtpi(4) = zroeshd(4, ip1, ish1)/zroeshd(1, ip1, ish1)
      xtpi(3) = zroeshd(3, ip1, ish1)/zroeshd(1, ip1, ish1)
      xtpi(1) = zroeshd(1, ip1, ish1)*zroeshd(1, ip1, ish1)
      help = zroeshd(3, ip1, ish1)**2 + zroeshd(4, ip1, ish1)**2
      xtpi(2) = gm1/ga*(zroeshd(1, ip1, ish1)*zroeshd(2, ip1, ish1)&
      &- 0.5d0*help)
    end if

! state 2
    xtpi(8) = zroeshu(4, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(7) = zroeshu(3, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(5) = zroeshu(1, ip2, ish2)*zroeshu(1, ip2, ish2)
    help = zroeshu(3, ip2, ish2)**2 + zroeshu(4, ip2, ish2)**2
    xtpi(6) = gm1/ga*(zroeshu(1, ip2, ish2)*zroeshu(2, ip2, ish2)&
    &- 0.5d0*help)
! state 3
    xtpi(12) = zroeshd(4, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(11) = zroeshd(3, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(9) = zroeshd(1, ip2, ish2)*zroeshd(1, ip2, ish2)
    help = zroeshd(3, ip2, ish2)**2 + zroeshd(4, ip2, ish2)**2
    xtpi(10) = gm1/ga*(zroeshd(1, ip2, ish2)*zroeshd(2, ip2, ish2)&
    &- 0.5d0*help)

    dx = vshnor(1, ip2, ish2)
    dy = vshnor(2, ip2, ish2)
    r23 = sqrt(ga*xtpi(10)/xtpi(9)) +&
    &gm1*0.5d0*(xtpi(11)*dx + xtpi(12)*dy)

! state 4
    xtpi(16) = zroeshd(4, ip3, ish3)/zroeshd(1, ip3, ish3)
    xtpi(15) = zroeshd(3, ip3, ish3)/zroeshd(1, ip3, ish3)
    xtpi(13) = zroeshd(1, ip3, ish3)*zroeshd(1, ip3, ish3)
    help = zroeshd(3, ip3, ish3)**2 + zroeshd(4, ip3, ish3)**2
    xtpi(14) = gm1/ga*(zroeshd(1, ip3, ish3)*zroeshd(2, ip3, ish3)&
    &- 0.5d0*help)

    dx = vshnor(1, ip3, ish3)
    dy = vshnor(2, ip3, ish3)
    dxr14 = dx
    dyr14 = dy

    r14 = sqrt(ga*xtpi(14)/xtpi(13)) +&
    &gm1*0.5d0*(xtpi(15)*dxr14 + xtpi(16)*dyr14)

! state 3 (again, from the OLD downstream state this time)
    xtpi(12) = zroeshdold(4, ip2, ish2)/zroeshdold(1, ip2, ish2)
    xtpi(11) = zroeshdold(3, ip2, ish2)/zroeshdold(1, ip2, ish2)
    xtpi(9) = zroeshdold(1, ip2, ish2)*zroeshdold(1, ip2, ish2)
    help = zroeshdold(3, ip2, ish2)**2 + zroeshdold(4, ip2, ish2)**2
    xtpi(10) = gm1/ga*(zroeshdold(1, ip2, ish2)*zroeshdold(2, ip2, ish2)&
    &- 0.5d0*help)
! state 4 (again, from the OLD downstream state)
    xtpi(16) = zroeshdold(4, ip3, ish3)/zroeshdold(1, ip3, ish3)
    xtpi(15) = zroeshdold(3, ip3, ish3)/zroeshdold(1, ip3, ish3)
    xtpi(13) = zroeshdold(1, ip3, ish3)*zroeshdold(1, ip3, ish3)
    help = zroeshdold(3, ip3, ish3)**2 + zroeshdold(4, ip3, ish3)**2
    xtpi(14) = gm1/ga*(zroeshdold(1, ip3, ish3)*zroeshdold(2, ip3, ish3)&
    &- 0.5d0*help)

    r14old = sqrt(ga*xtpi(14)/xtpi(13)) +&
    &gm1*0.5d0*(xtpi(15)*dxr14 + xtpi(16)*dyr14)

! slopes of the shocks
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    xtpi(17) = atan2(-dx, dy)
    dx = vshnor(1, ip2, ish2); dy = vshnor(2, ip2, ish2)
    xtpi(18) = atan2(-dx, dy)
    dx = vshnor(1, ip3, ish3); dy = vshnor(2, ip3, ish3)
    xtpi(19) = atan2(-dx, dy)

    xtpi(20) = -0.01

    wws = sqrt(wsh(1, ip1, ish1)**2 + wsh(2, ip1, ish1)**2)

    xtp = xtpi

    flag1 = .false.
    ifail = .false.
    call co_utp(xtpi, r14, dxr14, dyr14, r23, wws, xtp, flag1, ifail)

    if (ifail) then
      ifail = .false.
      xtpi(20) = -xtpi(20)
      call co_utp(xtpi, r14, dxr14, dyr14, r23, wws, xtp, flag1, ifail)
      if (ifail) then
        stop 'tp does not converge'
      end if
    end if

! zone 1 -> upstream mach stem
    z1v = sqrt(xtp(1))
    kine = 0.5d0*(xtp(3)*xtp(3) + xtp(4)*xtp(4))
    hh = ga/gm1*xtp(2)/xtp(1) + kine
    z2v = z1v*hh; z3v = z1v*xtp(3); z4v = z1v*xtp(4)
    zroeshu(1, ip3, ish3) = z1v; zroeshu(2, ip3, ish3) = z2v
    zroeshu(3, ip3, ish3) = z3v; zroeshu(4, ip3, ish3) = z4v

! zone 2 -> upstream reflected shock
    z1v = sqrt(xtp(5))
    kine = 0.5d0*(xtp(7)*xtp(7) + xtp(8)*xtp(8))
    hh = ga/gm1*xtp(6)/xtp(5) + kine
    z2v = z1v*hh; z3v = z1v*xtp(7); z4v = z1v*xtp(8)
    zroeshu(1, ip2, ish2) = z1v; zroeshu(2, ip2, ish2) = z2v
    zroeshu(3, ip2, ish2) = z3v; zroeshu(4, ip2, ish2) = z4v

! zone 3 -> downstream reflected shock + upstream contact discontinuity
    z1v = sqrt(xtp(9))
    kine = 0.5d0*(xtp(11)*xtp(11) + xtp(12)*xtp(12))
    hh = ga/gm1*xtp(10)/xtp(9) + kine
    z2v = z1v*hh; z3v = z1v*xtp(11); z4v = z1v*xtp(12)
    zroeshd(1, ip2, ish2) = z1v; zroeshd(2, ip2, ish2) = z2v
    zroeshd(3, ip2, ish2) = z3v; zroeshd(4, ip2, ish2) = z4v
    zroeshu(1, ip4, ish4) = z1v; zroeshu(2, ip4, ish4) = z2v
    zroeshu(3, ip4, ish4) = z3v; zroeshu(4, ip4, ish4) = z4v

! zone 4 -> downstream mach stem + downstream contact discontinuity
    z1v = sqrt(xtp(13))
    kine = 0.5d0*(xtp(15)*xtp(15) + xtp(16)*xtp(16))
    hh = ga/gm1*xtp(14)/xtp(13) + kine
    z2v = z1v*hh; z3v = z1v*xtp(15); z4v = z1v*xtp(16)
    zroeshd(1, ip3, ish3) = z1v; zroeshd(2, ip3, ish3) = z2v
    zroeshd(3, ip3, ish3) = z3v; zroeshd(4, ip3, ish3) = z4v
    zroeshd(1, ip4, ish4) = z1v; zroeshd(2, ip4, ish4) = z2v
    zroeshd(3, ip4, ish4) = z3v; zroeshd(4, ip4, ish4) = z4v

! triple-point velocity: contribution along incident shock, propagated to all 4 legs
    taux12 = cos(xtp(17)); tauy12 = sin(xtp(17))
    wsh(1, ip1, ish1) = wsh(1, ip1, ish1) + taux12*xtp(20)*1.0
    wsh(2, ip1, ish1) = wsh(2, ip1, ish1) + tauy12*xtp(20)*1.0
    wsh(1, ip2, ish2) = wsh(1, ip1, ish1); wsh(2, ip2, ish2) = wsh(2, ip1, ish1)
    wsh(1, ip3, ish3) = wsh(1, ip1, ish1); wsh(2, ip3, ish3) = wsh(2, ip1, ish1)
    wsh(1, ip4, ish4) = wsh(1, ip1, ish1); wsh(2, ip4, ish4) = wsh(2, ip1, ish1)
  end subroutine tp_solve_state

! QP: quadruple point. Ports fx_state_dps.f90's 'QP' branch verbatim
! (lines 491-798): recovers 5 zones' states, 4 shock-normal slopes,
! computes the quadruple-point velocity from both incident shocks'
! projected wsh, calls co_uqp, scatters the converged 5-zone state back.
  subroutine qp_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(quad_point_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat
    external co_uqp
    include 'paramt.h'

    real(wp) :: xtpi(24), xtp(24)
    real(wp) :: f1, f3, help, dx, dy, taux, tauy, cs, dum1, w
    real(wp) :: wqpx, wqpy, kine, hh, z1v, z2v, z3v, z4v
    integer(i4) :: i, ish1, ish2, ish3, ish4, ish5, ip1, ip2, ip3, ip4, ip5

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = this%ish(3); i = this%leg(3) - 1; ip3 = 1 + i*(nshockpoints(ish3) - 1)
    ish4 = this%ish(4); i = this%leg(4) - 1; ip4 = 1 + i*(nshockpoints(ish4) - 1)
    ish5 = this%ish(5); i = this%leg(5) - 1; ip5 = 1 + i*(nshockpoints(ish5) - 1)

    f1 = vshnor(1, ip1, ish1)*zroeshu(4, ip1, ish1) -&
    &vshnor(2, ip1, ish1)*zroeshu(3, ip1, ish1)
    f1 = -sign(1.d0, f1)

    f3 = vshnor(1, ip3, ish3)*zroeshu(4, ip3, ish3) -&
    &vshnor(2, ip3, ish3)*zroeshu(3, ip3, ish3)
    f3 = -sign(1.d0, f3)

! state 1
    if (f1 .gt. 0.d0) then
      xtpi(4) = zroeshu(4, ip1, ish1)/zroeshu(1, ip1, ish1)
      xtpi(3) = zroeshu(3, ip1, ish1)/zroeshu(1, ip1, ish1)
      xtpi(1) = zroeshu(1, ip1, ish1)*zroeshu(1, ip1, ish1)
      help = zroeshu(3, ip1, ish1)**2 + zroeshu(4, ip1, ish1)**2
      xtpi(2) = gm1/ga*(zroeshu(1, ip1, ish1)*zroeshu(2, ip1, ish1)&
      &- 0.5d0*help)
    else
      xtpi(4) = zroeshd(4, ip1, ish1)/zroeshd(1, ip1, ish1)
      xtpi(3) = zroeshd(3, ip1, ish1)/zroeshd(1, ip1, ish1)
      xtpi(1) = zroeshd(1, ip1, ish1)*zroeshd(1, ip1, ish1)
      help = zroeshd(3, ip1, ish1)**2 + zroeshd(4, ip1, ish1)**2
      xtpi(2) = gm1/ga*(zroeshd(1, ip1, ish1)*zroeshd(2, ip1, ish1)&
      &- 0.5d0*help)
    end if

! state 2
    xtpi(8) = zroeshu(4, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(7) = zroeshu(3, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(5) = zroeshu(1, ip2, ish2)*zroeshu(1, ip2, ish2)
    help = zroeshu(3, ip2, ish2)**2 + zroeshu(4, ip2, ish2)**2
    xtpi(6) = gm1/ga*(zroeshu(1, ip2, ish2)*zroeshu(2, ip2, ish2)&
    &- 0.5d0*help)

! state 3
    xtpi(12) = zroeshu(4, ip4, ish4)/zroeshu(1, ip4, ish4)
    xtpi(11) = zroeshu(3, ip4, ish4)/zroeshu(1, ip4, ish4)
    xtpi(9) = zroeshu(1, ip4, ish4)*zroeshu(1, ip4, ish4)
    help = zroeshu(3, ip4, ish4)**2 + zroeshu(4, ip4, ish4)**2
    xtpi(10) = gm1/ga*(zroeshu(1, ip4, ish4)*zroeshu(2, ip4, ish4)&
    &- 0.5d0*help)

! state 4
    xtpi(16) = zroeshd(4, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(15) = zroeshd(3, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(13) = zroeshd(1, ip2, ish2)*zroeshd(1, ip2, ish2)
    help = zroeshd(3, ip2, ish2)**2 + zroeshd(4, ip2, ish2)**2
    xtpi(14) = gm1/ga*(zroeshd(1, ip2, ish2)*zroeshd(2, ip2, ish2)&
    &- 0.5d0*help)

! state 5
    xtpi(20) = zroeshd(4, ip4, ish4)/zroeshd(1, ip4, ish4)
    xtpi(19) = zroeshd(3, ip4, ish4)/zroeshd(1, ip4, ish4)
    xtpi(17) = zroeshd(1, ip4, ish4)*zroeshd(1, ip4, ish4)
    help = zroeshd(3, ip4, ish4)**2 + zroeshd(4, ip4, ish4)**2
    xtpi(18) = gm1/ga*(zroeshd(1, ip4, ish4)*zroeshd(2, ip4, ish4)&
    &- 0.5d0*help)

! slopes of the shocks
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    xtpi(21) = atan2(-dx, dy)
    dx = vshnor(1, ip2, ish2); dy = vshnor(2, ip2, ish2)
    xtpi(22) = atan2(-dx, dy)
    dx = vshnor(1, ip3, ish3); dy = vshnor(2, ip3, ish3)
    xtpi(23) = atan2(-dx, dy)
    dx = vshnor(1, ip4, ish4); dy = vshnor(2, ip4, ish4)
    xtpi(24) = atan2(-dx, dy)

! velocity of the quadruple point
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    taux = -dy; tauy = +dx
    dx = vshnor(1, ip3, ish3); dy = vshnor(2, ip3, ish3)
    cs = dx*taux + dy*tauy
    dum1 = dx*wsh(1, ip3, ish3) + dy*wsh(2, ip3, ish3)
    w = sign(1.d0, dum1)*sqrt(wsh(1, ip3, ish3)**2 + wsh(2, ip3, ish3)**2)
    wqpx = w/cs*taux; wqpy = w/cs*tauy

    dx = vshnor(1, ip3, ish3); dy = vshnor(2, ip3, ish3)
    taux = -dy; tauy = +dx
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    cs = dx*taux + dy*tauy
    dum1 = dx*wsh(1, ip1, ish1) + dy*wsh(2, ip1, ish1)
    w = sign(1.d0, dum1)*sqrt(wsh(1, ip1, ish1)**2 + wsh(2, ip1, ish1)**2)
    wqpx = wqpx + w/cs*taux; wqpy = wqpy + w/cs*tauy

    xtp = xtpi

    call co_uqp(xtpi, wqpx, wqpy, xtp)

! zone 1
    z1v = sqrt(xtp(1))
    kine = 0.5d0*(xtp(3)*xtp(3) + xtp(4)*xtp(4))
    hh = ga/gm1*xtp(2)/xtp(1) + kine
    z2v = z1v*hh; z3v = z1v*xtp(3); z4v = z1v*xtp(4)
    if (f1 .gt. 0.d0) then
      zroeshu(1, ip1, ish1) = z1v; zroeshu(2, ip1, ish1) = z2v
      zroeshu(3, ip1, ish1) = z3v; zroeshu(4, ip1, ish1) = z4v
    else
      zroeshd(1, ip1, ish1) = z1v; zroeshd(2, ip1, ish1) = z2v
      zroeshd(3, ip1, ish1) = z3v; zroeshd(4, ip1, ish1) = z4v
    end if
    if (f3 .lt. 0.d0) then
      zroeshu(1, ip3, ish3) = z1v; zroeshu(2, ip3, ish3) = z2v
      zroeshu(3, ip3, ish3) = z3v; zroeshu(4, ip3, ish3) = z4v
    else
      zroeshd(1, ip3, ish3) = z1v; zroeshd(2, ip3, ish3) = z2v
      zroeshd(3, ip3, ish3) = z3v; zroeshd(4, ip3, ish3) = z4v
    end if

! zone 2
    z1v = sqrt(xtp(5))
    kine = 0.5d0*(xtp(7)*xtp(7) + xtp(8)*xtp(8))
    hh = ga/gm1*xtp(6)/xtp(5) + kine
    z2v = z1v*hh; z3v = z1v*xtp(7); z4v = z1v*xtp(8)
    if (f1 .gt. 0.d0) then
      zroeshd(1, ip1, ish1) = z1v; zroeshd(2, ip1, ish1) = z2v
      zroeshd(3, ip1, ish1) = z3v; zroeshd(4, ip1, ish1) = z4v
    else
      zroeshu(1, ip1, ish1) = z1v; zroeshu(2, ip1, ish1) = z2v
      zroeshu(3, ip1, ish1) = z3v; zroeshu(4, ip1, ish1) = z4v
    end if
    zroeshu(1, ip2, ish2) = z1v; zroeshu(2, ip2, ish2) = z2v
    zroeshu(3, ip2, ish2) = z3v; zroeshu(4, ip2, ish2) = z4v

! zone 3
    z1v = sqrt(xtp(9))
    kine = 0.5d0*(xtp(11)*xtp(11) + xtp(12)*xtp(12))
    hh = ga/gm1*xtp(10)/xtp(9) + kine
    z2v = z1v*hh; z3v = z1v*xtp(11); z4v = z1v*xtp(12)
    if (f3 .lt. 0.d0) then
      zroeshd(1, ip3, ish3) = z1v; zroeshd(2, ip3, ish3) = z2v
      zroeshd(3, ip3, ish3) = z3v; zroeshd(4, ip3, ish3) = z4v
    else
      zroeshu(1, ip3, ish3) = z1v; zroeshu(2, ip3, ish3) = z2v
      zroeshu(3, ip3, ish3) = z3v; zroeshu(4, ip3, ish3) = z4v
    end if
    zroeshu(1, ip4, ish4) = z1v; zroeshu(2, ip4, ish4) = z2v
    zroeshu(3, ip4, ish4) = z3v; zroeshu(4, ip4, ish4) = z4v

! zone 4
    z1v = sqrt(xtp(13))
    kine = 0.5d0*(xtp(15)*xtp(15) + xtp(16)*xtp(16))
    hh = ga/gm1*xtp(14)/xtp(13) + kine
    z2v = z1v*hh; z3v = z1v*xtp(15); z4v = z1v*xtp(16)
    zroeshd(1, ip2, ish2) = z1v; zroeshd(2, ip2, ish2) = z2v
    zroeshd(3, ip2, ish2) = z3v; zroeshd(4, ip2, ish2) = z4v
    zroeshu(1, ip5, ish5) = z1v; zroeshu(2, ip5, ish5) = z2v
    zroeshu(3, ip5, ish5) = z3v; zroeshu(4, ip5, ish5) = z4v

! zone 5
    z1v = sqrt(xtp(17))
    kine = 0.5d0*(xtp(19)*xtp(19) + xtp(20)*xtp(20))
    hh = ga/gm1*xtp(18)/xtp(17) + kine
    z2v = z1v*hh; z3v = z1v*xtp(19); z4v = z1v*xtp(20)
    zroeshd(1, ip4, ish4) = z1v; zroeshd(2, ip4, ish4) = z2v
    zroeshd(3, ip4, ish4) = z3v; zroeshd(4, ip4, ish4) = z4v
    zroeshd(1, ip5, ish5) = z1v; zroeshd(2, ip5, ish5) = z2v
    zroeshd(3, ip5, ish5) = z3v; zroeshd(4, ip5, ish5) = z4v

! quadruple-point velocity, propagated to all 5 legs
    wsh(1, ip1, ish1) = wqpx; wsh(2, ip1, ish1) = wqpy
    wsh(1, ip2, ish2) = wsh(1, ip1, ish1); wsh(2, ip2, ish2) = wsh(2, ip1, ish1)
    wsh(1, ip3, ish3) = wsh(1, ip1, ish1); wsh(2, ip3, ish3) = wsh(2, ip1, ish1)
    wsh(1, ip4, ish4) = wsh(1, ip1, ish1); wsh(2, ip4, ish4) = wsh(2, ip1, ish1)
    wsh(1, ip5, ish5) = wsh(1, ip1, ish1); wsh(2, ip5, ish5) = wsh(2, ip1, ish1)
  end subroutine qp_solve_state

! TE: trailing edge. Ports fx_state_dps.f90's 'TE' branch verbatim
! (lines 801-991), including its dead-but-present ish1=0/ip1=0 and
! ish3=0/ip3=0 assignments and the two VSHNOR(:,0,0) reads that follow
! them -- both are confirmed unreachable in every regression fixture
! (TE appears in zero tests/*/sh00.dat), so this is preserved bug-for-
! bug rather than "cleaned up" per the merged-3.2 plan's bug policy.
! TE reuses co_uqp (the quadruple-point solver) with states 1 and 3
! defaulted to placeholder values and no incident-shock legs at all.
  subroutine te_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(trailing_edge_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat
    external co_uqp
    include 'paramt.h'

    real(wp) :: xtpi(24), xtp(24)
    real(wp) :: help, dx, dy, wqpx, wqpy, kine, hh, z1v, z2v, z3v, z4v
    integer(i4) :: i, ish1, ish2, ish3, ish4, ish5, ip1, ip2, ip3, ip4, ip5

    ish1 = 0; i = 0; ip1 = 0
    ish2 = this%ish(1); i = this%leg(1) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = 0; i = 0; ip3 = 0
    ish4 = this%ish(3); i = this%leg(3) - 1; ip4 = 1 + i*(nshockpoints(ish4) - 1)
    ish5 = this%ish(2); i = this%leg(2) - 1; ip5 = 1 + i*(nshockpoints(ish5) - 1)

! state 1: placeholder (no real incident shock 1 at a trailing edge)
    xtpi(4) = 1.0; xtpi(3) = 1.0; xtpi(1) = 1.0; xtpi(2) = 1.0

! state 2
    xtpi(8) = zroeshu(4, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(7) = zroeshu(3, ip2, ish2)/zroeshu(1, ip2, ish2)
    xtpi(5) = zroeshu(1, ip2, ish2)*zroeshu(1, ip2, ish2)
    help = zroeshu(3, ip2, ish2)**2 + zroeshu(4, ip2, ish2)**2
    xtpi(6) = gm1/ga*(zroeshu(1, ip2, ish2)*zroeshu(2, ip2, ish2)&
    &- 0.5d0*help)

! state 3
    xtpi(12) = zroeshu(4, ip4, ish4)/zroeshu(1, ip4, ish4)
    xtpi(11) = zroeshu(3, ip4, ish4)/zroeshu(1, ip4, ish4)
    xtpi(9) = zroeshu(1, ip4, ish4)*zroeshu(1, ip4, ish4)
    help = zroeshu(3, ip4, ish4)**2 + zroeshu(4, ip4, ish4)**2
    xtpi(10) = gm1/ga*(zroeshu(1, ip4, ish4)*zroeshu(2, ip4, ish4)&
    &- 0.5d0*help)

! state 4
    xtpi(16) = zroeshd(4, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(15) = zroeshd(3, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(13) = zroeshd(1, ip2, ish2)*zroeshd(1, ip2, ish2)
    help = zroeshd(3, ip2, ish2)**2 + zroeshd(4, ip2, ish2)**2
    xtpi(14) = gm1/ga*(zroeshd(1, ip2, ish2)*zroeshd(2, ip2, ish2)&
    &- 0.5d0*help)

! state 5
    xtpi(20) = zroeshd(4, ip4, ish4)/zroeshd(1, ip4, ish4)
    xtpi(19) = zroeshd(3, ip4, ish4)/zroeshd(1, ip4, ish4)
    xtpi(17) = zroeshd(1, ip4, ish4)*zroeshd(1, ip4, ish4)
    help = zroeshd(3, ip4, ish4)**2 + zroeshd(4, ip4, ish4)**2
    xtpi(18) = gm1/ga*(zroeshd(1, ip4, ish4)*zroeshd(2, ip4, ish4)&
    &- 0.5d0*help)

! slopes: shock 1 and shock 3 are placeholders (1.0); the DX/DY reads
! below feed nothing (dead, unreachable -- see subroutine header note)
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    xtpi(21) = 1.0
    dx = vshnor(1, ip2, ish2); dy = vshnor(2, ip2, ish2)
    xtpi(22) = atan2(-dx, dy)
    dx = vshnor(1, ip3, ish3); dy = vshnor(2, ip3, ish3)
    xtpi(23) = 1.0
    dx = vshnor(1, ip4, ish4); dy = vshnor(2, ip4, ish4)
    xtpi(24) = atan2(-dx, dy)

    xtp = xtpi
    wqpx = 0.0_wp; wqpy = 0.0_wp
    call co_uqp(xtpi, wqpx, wqpy, xtp)

! zone 1 (computed but not scattered anywhere, matching the original)
    z1v = sqrt(xtp(1))
    kine = 0.5d0*(xtp(3)*xtp(3) + xtp(4)*xtp(4))
    hh = ga/gm1*xtp(2)/xtp(1) + kine
    z2v = z1v*hh; z3v = z1v*xtp(3); z4v = z1v*xtp(4)

! zone 2 -> upstream reflected shock
    z1v = sqrt(xtp(5))
    kine = 0.5d0*(xtp(7)*xtp(7) + xtp(8)*xtp(8))
    hh = ga/gm1*xtp(6)/xtp(5) + kine
    z2v = z1v*hh; z3v = z1v*xtp(7); z4v = z1v*xtp(8)
    zroeshu(1, ip2, ish2) = z1v; zroeshu(2, ip2, ish2) = z2v
    zroeshu(3, ip2, ish2) = z3v; zroeshu(4, ip2, ish2) = z4v

! zone 3 -> upstream reflected shock 2
    z1v = sqrt(xtp(9))
    kine = 0.5d0*(xtp(11)*xtp(11) + xtp(12)*xtp(12))
    hh = ga/gm1*xtp(10)/xtp(9) + kine
    z2v = z1v*hh; z3v = z1v*xtp(11); z4v = z1v*xtp(12)
    zroeshu(1, ip4, ish4) = z1v; zroeshu(2, ip4, ish4) = z2v
    zroeshu(3, ip4, ish4) = z3v; zroeshu(4, ip4, ish4) = z4v

! zone 4 -> downstream reflected shock 1 + upstream contact discontinuity
    z1v = sqrt(xtp(13))
    kine = 0.5d0*(xtp(15)*xtp(15) + xtp(16)*xtp(16))
    hh = ga/gm1*xtp(14)/xtp(13) + kine
    z2v = z1v*hh; z3v = z1v*xtp(15); z4v = z1v*xtp(16)
    zroeshd(1, ip2, ish2) = z1v; zroeshd(2, ip2, ish2) = z2v
    zroeshd(3, ip2, ish2) = z3v; zroeshd(4, ip2, ish2) = z4v
    zroeshd(1, ip5, ish5) = z1v; zroeshd(2, ip5, ish5) = z2v
    zroeshd(3, ip5, ish5) = z3v; zroeshd(4, ip5, ish5) = z4v

! zone 5 -> downstream reflected shock 2 + downstream contact discontinuity
    z1v = sqrt(xtp(17))
    kine = 0.5d0*(xtp(19)*xtp(19) + xtp(20)*xtp(20))
    hh = ga/gm1*xtp(18)/xtp(17) + kine
    z2v = z1v*hh; z3v = z1v*xtp(19); z4v = z1v*xtp(20)
    zroeshd(1, ip4, ish4) = z1v; zroeshd(2, ip4, ish4) = z2v
    zroeshd(3, ip4, ish4) = z3v; zroeshd(4, ip4, ish4) = z4v
    zroeshu(1, ip5, ish5) = z1v; zroeshu(2, ip5, ish5) = z2v
    zroeshu(3, ip5, ish5) = z3v; zroeshu(4, ip5, ish5) = z4v

! velocity: zero at both real legs
    wsh(1, ip2, ish2) = 0.0; wsh(2, ip2, ish2) = 0.0
    wsh(1, ip4, ish4) = 0.0; wsh(2, ip4, ish4) = 0.0
    wsh(1, ip5, ish5) = 0.0; wsh(2, ip5, ish5) = 0.0
  end subroutine te_solve_state

! RRX/RR: regular reflection. Ports fx_state_dps.f90's 'RRX'/'RR' branch
! verbatim (lines 994-1230). RRX assumes an axis-aligned wall
! (taux=1,tauy=0); RR recomputes the wall tangent from boundary geometry
! using ia/ja/iclr/nclr/corg and the legacy rank-1-mis-indexed
! ispclr_flat (= ispclr(isppnts) in the original, where ispclr is
! actually a (5,*) array elsewhere -- a real latent bug, confirmed
! reachable only via 'RR' which no regression fixture exercises;
! preserved bit-for-bit rather than fixed, per the merged-3.2 plan's bug
! policy). Newton-solved via co_urr.
  subroutine rr_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(regular_reflection_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat
    external co_urr
    include 'paramt.h'

    real(wp) :: xtpi(15), xtp(15)
    real(wp) :: help, dx, dy, taux, tauy, wws, wrr, xi, yi
    real(wp) :: dum, dum1, dumx1, dumy1, dumx2, dumy2
    real(wp) :: kine, hh, z1v, z2v, z3v, z4v
    integer(i4) :: i, ish1, ish2, ip1, ip2, clr, j, k, kp1, bbgn, bend

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)

! state 1
    xtpi(4) = zroeshu(4, ip1, ish1)/zroeshu(1, ip1, ish1)
    xtpi(3) = zroeshu(3, ip1, ish1)/zroeshu(1, ip1, ish1)
    xtpi(1) = zroeshu(1, ip1, ish1)*zroeshu(1, ip1, ish1)
    help = zroeshu(3, ip1, ish1)**2 + zroeshu(4, ip1, ish1)**2
    xtpi(2) = gm1/ga*(zroeshu(1, ip1, ish1)*zroeshu(2, ip1, ish1)&
    &- 0.5d0*help)

! state 2 (downstream of shock 1, not shock 2 -- matches original comment)
    xtpi(8) = zroeshd(4, ip1, ish1)/zroeshd(1, ip1, ish1)
    xtpi(7) = zroeshd(3, ip1, ish1)/zroeshd(1, ip1, ish1)
    xtpi(5) = zroeshd(1, ip1, ish1)*zroeshd(1, ip1, ish1)
    help = zroeshd(3, ip1, ish1)**2 + zroeshd(4, ip1, ish1)**2
    xtpi(6) = gm1/ga*(zroeshd(1, ip1, ish1)*zroeshd(2, ip1, ish1)&
    &- 0.5d0*help)

! state 3
    xtpi(12) = zroeshd(4, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(11) = zroeshd(3, ip2, ish2)/zroeshd(1, ip2, ish2)
    xtpi(9) = zroeshd(1, ip2, ish2)*zroeshd(1, ip2, ish2)
    help = zroeshd(3, ip2, ish2)**2 + zroeshd(4, ip2, ish2)**2
    xtpi(10) = gm1/ga*(zroeshd(1, ip2, ish2)*zroeshd(2, ip2, ish2)&
    &- 0.5d0*help)

! slopes of shocks
    dx = vshnor(1, ip1, ish1); dy = vshnor(2, ip1, ish1)
    xtpi(13) = atan2(-dx, dy)
    dx = vshnor(1, ip2, ish2); dy = vshnor(2, ip2, ish2)
    xtpi(14) = atan2(-dx, dy)

    wws = sqrt(wsh(1, ip1, ish1)**2 + wsh(2, ip1, ish1)**2)

    dx = vshnor(1, ip1, ish1)
    dy = vshnor(2, ip1, ish1)

! versor tangent to the wall -- axis-aligned default (RRX)
    taux = 1.0d0
    tauy = 0.0d0

! curved wall (RR): recompute tangent from the boundary segment the
! shock point sits on
    if (this%curved) then
      xi = xysh(1, ip1, ish1)
      yi = xysh(2, ip1, ish1)

      clr = 5
      do clr = 1, nclr
        if (iclr(clr) .eq. ispclr_flat) exit
      end do
      bbgn = ia(clr)
      bend = ia(clr + 1) - 1
      do j = bbgn, bend - 1
        k = ja(j)
        kp1 = ja(j + 1)
        dumx1 = corg(1, k)
        dumy1 = corg(2, k)
        dumx2 = corg(1, kp1)
        dumy2 = corg(2, kp1)

        dum = sqrt((dumx1 - dumx2)**2 + (dumy1 - dumy2)**2)
        dum1 = sqrt((xi - dumx2)**2 + (yi - dumy2)**2) +&
        &sqrt((xi - dumx1)**2 + (yi - dumy1)**2)
        if (dum1 - dum .le. 1.0e-5) then
          taux = -(dumx1 - dumx2)
          tauy = -(dumy1 - dumy2)
          dum = sqrt(taux**2 + tauy**2)
          taux = taux/dum
          tauy = tauy/dum
        end if
      end do
    end if

    dum = dx*taux + dy*tauy
    dum1 = dx*wsh(1, ip1, ish1) + dy*wsh(2, ip1, ish1)
    wrr = sign(1.d0, dum1)*wws/dum

    xtpi(15) = wrr

    xtp = xtpi
    call co_urr(xtpi, taux, tauy, xtp)

! zone 1 -> upstream incident shock
    z1v = sqrt(xtp(1))
    kine = 0.5d0*(xtp(3)*xtp(3) + xtp(4)*xtp(4))
    hh = ga/gm1*xtp(2)/xtp(1) + kine
    z2v = z1v*hh; z3v = z1v*xtp(3); z4v = z1v*xtp(4)
    zroeshu(1, ip1, ish1) = z1v; zroeshu(2, ip1, ish1) = z2v
    zroeshu(3, ip1, ish1) = z3v; zroeshu(4, ip1, ish1) = z4v

! zone 2 -> downstream incident shock + upstream reflected shock
    z1v = sqrt(xtp(5))
    kine = 0.5d0*(xtp(7)*xtp(7) + xtp(8)*xtp(8))
    hh = ga/gm1*xtp(6)/xtp(5) + kine
    z2v = z1v*hh; z3v = z1v*xtp(7); z4v = z1v*xtp(8)
    zroeshd(1, ip1, ish1) = z1v; zroeshd(2, ip1, ish1) = z2v
    zroeshd(3, ip1, ish1) = z3v; zroeshd(4, ip1, ish1) = z4v
    zroeshu(1, ip2, ish2) = z1v; zroeshu(2, ip2, ish2) = z2v
    zroeshu(3, ip2, ish2) = z3v; zroeshu(4, ip2, ish2) = z4v

! zone 3 -> downstream reflected shock
    z1v = sqrt(xtp(9))
    kine = 0.5d0*(xtp(11)*xtp(11) + xtp(12)*xtp(12))
    hh = ga/gm1*xtp(10)/xtp(9) + kine
    z2v = z1v*hh; z3v = z1v*xtp(11); z4v = z1v*xtp(12)
    zroeshd(1, ip2, ish2) = z1v; zroeshd(2, ip2, ish2) = z2v
    zroeshd(3, ip2, ish2) = z3v; zroeshd(4, ip2, ish2) = z4v

! velocity of the reflection point, propagated to both legs
    wrr = xtp(15)
    wsh(1, ip1, ish1) = wrr*taux
    wsh(2, ip1, ish1) = wrr*tauy
    wsh(1, ip2, ish2) = wsh(1, ip1, ish1)
    wsh(2, ip2, ish2) = wsh(2, ip1, ish1)
  end subroutine rr_solve_state

! WPNRX/WPNRY/IPX/IPY/OPX/OPY. Ports fx_state_dps.f90's IPX/IPY branch
! (lines 82-100) when restore_only, else its OPX/OPY/WPNRX/WPNRY branch
! (lines 103-143, confirmed byte-identical between OPX/OPY and
! WPNRX/WPNRY apart from variable names -- both fold into one
! axis-keyed implementation here).
  subroutine wf_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: ndof, ndim
    class(wall_float_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat

    real(wp) :: ws, dx, dy
    integer(i4) :: i, ish, ip, k

    ish = this%ish(1); i = this%leg(1) - 1; ip = 1 + i*(nshockpoints(ish) - 1)

    if (this%restore_only) then
      do k = 1, ndof
        zroeshu(k, ip, ish) = zroeshuold(k, ip, ish)
        zroeshd(k, ip, ish) = zroeshdold(k, ip, ish)
      end do
      do k = 1, ndim
        wsh(k, ip, ish) = 0.0d+0
      end do
    else
      ws = sqrt(wsh(1, ip, ish)**2 + wsh(2, ip, ish)**2)
      dx = wsh(1, ip, ish)/ws
      dy = wsh(2, ip, ish)/ws
      if (this%axis .eq. 'X') then
        wsh(1, ip, ish) = ws/dx
        wsh(2, ip, ish) = 0.0d+0
      else
        wsh(1, ip, ish) = 0.0d+0
        wsh(2, ip, ish) = ws/dy
      end if
    end if
  end subroutine wf_solve_state

! FWP: floating wall point. Ports fx_state_dps.f90's 'FWP' branch
! verbatim (lines 145-173), including its confirmed pre-existing bug --
! wsh(2,ip,ish) is set from the just-overwritten wsh(1,ip,ish)/dx, not
! from the original wsh(2,ip,ish) -- preserved bit-for-bit (not fixed)
! since NACA0012_M080_A0 (the only fixture exercising FWP) checksums
! against this exact behavior.
  subroutine fwp_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    class(floating_wall_point_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat

    real(wp) :: dx, dy, dum, dumx, dumy
    integer(i4) :: i, ish, ip

    ish = this%ish(1); i = this%leg(1) - 1; ip = 1 + i*(nshockpoints(ish) - 1)

    dx = vshnor(1, ip, ish)
    dy = vshnor(2, ip, ish)

    dum = zroeshu(3, ip, ish)*dx + zroeshu(4, ip, ish)*dy
    dumx = zroeshu(3, ip, ish) - dum*dx
    dumy = zroeshu(4, ip, ish) - dum*dy
    zroeshu(3, ip, ish) = zroeshu(3, ip, ish) - dumx
    zroeshu(4, ip, ish) = zroeshu(4, ip, ish) - dumy

    dum = zroeshd(3, ip, ish)*dx + zroeshd(4, ip, ish)*dy
    dumx = zroeshd(3, ip, ish) - dum*dx
    dumy = zroeshd(4, ip, ish) - dum*dy
    zroeshd(3, ip, ish) = zroeshd(3, ip, ish) - dumx
    zroeshd(4, ip, ish) = zroeshd(4, ip, ish) - dumy

    if (vshnor(2, ip, ish) .gt. 0.3) then
      wsh(1, ip, ish) = wsh(1, ip, ish)/dx
      wsh(2, ip, ish) = wsh(1, ip, ish)/dx
    end if
  end subroutine fwp_solve_state

! EP: end point. Ports fx_state_dps.f90's 'EP' branch verbatim (lines
! 1233-1278): relocates the endpoint along the internal point's
! corrected tangent by dxcell, copies the internal point's velocity,
! and averages the internal point's up/downstream state into both the
! endpoint's up- and downstream slots.
  subroutine ep_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    use mod_constants, only: naddholesmax, ndim, nprdbndmax
    class(end_point_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat
    include 'paramt.h'

    real(wp) :: nx, ny, taux, tauy, dx, dy, dum
    integer(i4) :: i, ish, ip, ip1

    ish = this%ish(1); i = this%leg(1) - 1
    ip = 1 + i*(nshockpoints(ish) - 1)
    ip1 = 2 + i*(nshockpoints(ish) - 3)

    nx = vshnor(1, ip1, ish)
    ny = vshnor(2, ip1, ish)
    taux = ny
    tauy = -nx

    dx = xysh(1, ip, ish) - xysh(1, ip1, ish)
    dy = xysh(2, ip, ish) - xysh(2, ip1, ish)

    dum = taux*dx + tauy*dy
    if (dum .lt. 0.d+0) then
      taux = -taux
      tauy = -tauy
    end if

    xysh(1, ip, ish) = xysh(1, ip1, ish) + taux*dxcell
    xysh(2, ip, ish) = xysh(2, ip1, ish) + tauy*dxcell

    wsh(1, ip, ish) = wsh(1, ip1, ish)
    wsh(2, ip, ish) = wsh(2, ip1, ish)

    zroeshd(1, ip, ish) = 0.5d0*(zroeshd(1, ip1, ish) + zroeshu(1, ip1, ish))
    zroeshd(2, ip, ish) = 0.5d0*(zroeshd(2, ip1, ish) + zroeshu(2, ip1, ish))
    zroeshd(3, ip, ish) = 0.5d0*(zroeshd(3, ip1, ish) + zroeshu(3, ip1, ish))
    zroeshd(4, ip, ish) = 0.5d0*(zroeshd(4, ip1, ish) + zroeshu(4, ip1, ish))

    zroeshu(1, ip, ish) = 0.5d0*(zroeshd(1, ip1, ish) + zroeshu(1, ip1, ish))
    zroeshu(2, ip, ish) = 0.5d0*(zroeshd(2, ip1, ish) + zroeshu(2, ip1, ish))
    zroeshu(3, ip, ish) = 0.5d0*(zroeshd(3, ip1, ish) + zroeshu(3, ip1, ish))
    zroeshu(4, ip, ish) = 0.5d0*(zroeshd(4, ip1, ish) + zroeshu(4, ip1, ish))
  end subroutine ep_solve_state

! SP: start/sonic point. Ports fx_state_dps.f90's 'SP' branch verbatim
! (lines 1281-1292) -- just zeroes the shock velocity; the point's
! actual relocation happens elsewhere (co_norm.f90's real SP branch,
! increment 8), matching the original's own comment "the point is moved
! in a different way".
  subroutine sonic_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    class(start_point_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat

    integer(i4) :: i, ish, ip

    ish = this%ish(1); i = this%leg(1) - 1
    ip = 1 + i*(nshockpoints(ish) - 1)

    wsh(1, ip, ish) = 0.0
    wsh(2, ip, ish) = 0.0
  end subroutine sonic_solve_state

! C/PC: connection / periodic connection. Ports fx_state_dps.f90's
! 'C'/'PC' branch verbatim (lines 1295-1329): copies point 2's state
! from point 1's, and (C only) point 2's coordinates too.
  subroutine conn_solve_state(this, xysh, zroeshu, zroeshd, zroeshuold,&
  &zroeshdold, vshnor, wsh, nshockpoints, ia, ja, iclr, nclr, corg,&
  &ispclr_flat)
    class(connection_t), intent(inout) :: this
    real(wp), intent(inout) :: xysh(:, :, :)
    real(wp), intent(inout) :: zroeshu(:, :, :), zroeshd(:, :, :)
    real(wp), intent(in) :: zroeshuold(:, :, :), zroeshdold(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :)
    real(wp), intent(inout) :: wsh(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    integer(i4), intent(in) :: nclr, ia(*), ja(*), iclr(nclr)
    real(wp), intent(in) :: corg(ndim, *)
    integer(i4), intent(in) :: ispclr_flat

    integer(i4) :: i, ish1, ish2, ip1, ip2

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)

    if (.not. this%periodic) then
      xysh(1, ip2, ish2) = xysh(1, ip1, ish1)
      xysh(2, ip2, ish2) = xysh(2, ip1, ish1)
    end if

    wsh(1, ip2, ish2) = wsh(1, ip1, ish1)
    wsh(2, ip2, ish2) = wsh(2, ip1, ish1)

    zroeshd(1, ip2, ish2) = zroeshd(1, ip1, ish1)
    zroeshd(2, ip2, ish2) = zroeshd(2, ip1, ish1)
    zroeshd(3, ip2, ish2) = zroeshd(3, ip1, ish1)
    zroeshd(4, ip2, ish2) = zroeshd(4, ip1, ish1)

    zroeshu(1, ip2, ish2) = zroeshu(1, ip1, ish1)
    zroeshu(2, ip2, ish2) = zroeshu(2, ip1, ish1)
    zroeshu(3, ip2, ish2) = zroeshu(3, ip1, ish1)
    zroeshu(4, ip2, ish2) = zroeshu(4, ip1, ish1)
  end subroutine conn_solve_state

! Default displace: does nothing to xyshu/xyshd, matching co_pnt_dspl.f90's
! empty FWP and PC branches (the only two codes that reach it -- every
! other type overrides below).
  subroutine sp_displace_noop(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    class(special_point_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
  end subroutine sp_displace_noop

! wall_float_t: ports co_pnt_dspl.f90's {IPX,OPX,WPNRX}/{IPY,OPY,WPNRY}
! branches verbatim (lines 78-104) -- purely axis-keyed here (unlike
! solve_state, this chain treats IP*/OP*/WPN* identically regardless of
! restore_only). Updates dx_carry/dy_carry from vshnor, mirroring the
! legacy code's own assignment into its shared dx/dy locals -- see
! sp_displace_if's header for why that matters (rr_displace's RRX path
! depends on reading exactly this stale value).
  subroutine wf_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    use mod_constants, only: naddholesmax, nprdbndmax
    class(wall_float_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
    include 'paramt.h'

    integer(i4) :: i, ish, ip

    ish = this%ish(1); i = this%leg(1) - 1
    ip = 1 + i*(nshockpoints(ish) - 1)
    dx_carry = vshnor(1, ip, ish)
    dy_carry = vshnor(2, ip, ish)

    if (this%axis .eq. 'X') then
      xyshu(1, ip, ish) = xysh(1, ip, ish) + 0.5d+0*eps/dx_carry
      xyshu(2, ip, ish) = xysh(2, ip, ish)
      xyshd(1, ip, ish) = xysh(1, ip, ish) - 0.5d+0*eps/dx_carry
      xyshd(2, ip, ish) = xysh(2, ip, ish)
    else
      xyshu(1, ip, ish) = xysh(1, ip, ish)
      xyshu(2, ip, ish) = xysh(2, ip, ish) + 0.5d+0*eps/dy_carry
      xyshd(1, ip, ish) = xysh(1, ip, ish)
      xyshd(2, ip, ish) = xysh(2, ip, ish) - 0.5d+0*eps/dy_carry
    end if
  end subroutine wf_displace

! TP: ports co_pnt_dspl.f90's 'TP' branch verbatim (lines 106-258) --
! moves all 4 legs along their local shock tangent by eps, then overlaps
! (averages) the appropriate up/downstream point pairs depending on
! whether the incident and reflected shocks belong to opposite or the
! same family.
  subroutine tp_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    use mod_constants, only: naddholesmax, nprdbndmax
    class(triple_point_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
    include 'paramt.h'

    real(wp) :: tx, ty, dum, f1, f2
    integer(i4) :: i, k, ish1, ish2, ish3, ish4, ip1, ip2, ip3, ip4

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = this%ish(3); i = this%leg(3) - 1; ip3 = 1 + i*(nshockpoints(ish3) - 1)
    ish4 = this%ish(4); i = this%leg(4) - 1; ip4 = 1 + i*(nshockpoints(ish4) - 1)

    f1 = vshnor(1, ip1, ish1)*zroeshu(4, ip1, ish1) -&
    &vshnor(2, ip1, ish1)*zroeshu(3, ip1, ish1)
    f1 = -sign(1.d0, f1)

    f2 = vshnor(1, ip2, ish2)*zroeshu(4, ip2, ish2) -&
    &vshnor(2, ip2, ish2)*zroeshu(3, ip2, ish2)
    f2 = -sign(1.d0, f2)

! move incident shock point
    if (ip1 .eq. 1) then
      tx = xysh(1, ip1, ish1) - xysh(1, ip1 - 1, ish1)
      ty = xysh(2, ip1, ish1) - xysh(2, ip1 - 1, ish1)
    else
      tx = xysh(1, ip1 - 1, ish1) - xysh(1, ip1, ish1)
      ty = xysh(2, ip1 - 1, ish1) - xysh(2, ip1, ish1)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip1, ish1) = xyshu(1, ip1, ish1) + eps*tx
    xyshu(2, ip1, ish1) = xyshu(2, ip1, ish1) + eps*ty
    xyshd(1, ip1, ish1) = xyshd(1, ip1, ish1) + eps*tx
    xyshd(2, ip1, ish1) = xyshd(2, ip1, ish1) + eps*ty

! move reflected shock point
    if (ip2 .eq. 1) then
      tx = xysh(1, ip2, ish2) - xysh(1, ip2 - 1, ish2)
      ty = xysh(2, ip2, ish2) - xysh(2, ip2 - 1, ish2)
    else
      tx = xysh(1, ip2 - 1, ish2) - xysh(1, ip2, ish2)
      ty = xysh(2, ip2 - 1, ish2) - xysh(2, ip2, ish2)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip2, ish2) = xyshu(1, ip2, ish2) + eps*tx
    xyshu(2, ip2, ish2) = xyshu(2, ip2, ish2) + eps*ty
    xyshd(1, ip2, ish2) = xyshd(1, ip2, ish2) + eps*tx
    xyshd(2, ip2, ish2) = xyshd(2, ip2, ish2) + eps*ty

! move mach stem point
    if (ip3 .eq. 1) then
      tx = xysh(1, ip3, ish3) - xysh(1, ip3 - 1, ish3)
      ty = xysh(2, ip3, ish3) - xysh(2, ip3 - 1, ish3)
    else
      tx = xysh(1, ip3 - 1, ish3) - xysh(1, ip3, ish3)
      ty = xysh(2, ip3 - 1, ish3) - xysh(2, ip3, ish3)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip3, ish3) = xyshu(1, ip3, ish3) + eps*tx
    xyshu(2, ip3, ish3) = xyshu(2, ip3, ish3) + eps*ty
    xyshd(1, ip3, ish3) = xyshd(1, ip3, ish3) + eps*tx
    xyshd(2, ip3, ish3) = xyshd(2, ip3, ish3) + eps*ty

! move contact discontinuity point
    if (ip4 .eq. 1) then
      tx = xysh(1, ip4, ish4) - xysh(1, ip4 - 1, ish4)
      ty = xysh(2, ip4, ish4) - xysh(2, ip4 - 1, ish4)
    else
      tx = xysh(1, ip4 - 1, ish4) - xysh(1, ip4, ish4)
      ty = xysh(2, ip4 - 1, ish4) - xysh(2, ip4, ish4)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip4, ish4) = xyshu(1, ip4, ish4) + eps*tx
    xyshu(2, ip4, ish4) = xyshu(2, ip4, ish4) + eps*ty
    xyshd(1, ip4, ish4) = xyshd(1, ip4, ish4) + eps*tx
    xyshd(2, ip4, ish4) = xyshd(2, ip4, ish4) + eps*ty

    if (f1*f2 .lt. 0.) then
      do k = 1, ndim
        dum = 0.5*(xyshu(k, ip1, ish1) + xyshu(k, ip3, ish3))
        xyshu(k, ip1, ish1) = dum; xyshu(k, ip3, ish3) = dum
      end do
      do k = 1, ndim
        dum = 0.5*(xyshd(k, ip1, ish1) + xyshu(k, ip2, ish2))
        xyshd(k, ip1, ish1) = dum; xyshu(k, ip2, ish2) = dum
      end do
    else
      do k = 1, ndim
        dum = 0.5*(xyshd(k, ip1, ish1) + xyshu(k, ip3, ish3))
        xyshd(k, ip1, ish1) = dum; xyshu(k, ip3, ish3) = dum
      end do
      do k = 1, ndim
        dum = 0.5*(xyshu(k, ip1, ish1) + xyshu(k, ip2, ish2))
        xyshu(k, ip1, ish1) = dum; xyshu(k, ip2, ish2) = dum
      end do
    end if

    do k = 1, ndim
      dum = 0.5*(xyshd(k, ip3, ish3) + xyshd(k, ip4, ish4))
      xyshd(k, ip3, ish3) = dum; xyshd(k, ip4, ish4) = dum
    end do

    do k = 1, ndim
      dum = 0.5*(xyshd(k, ip2, ish2) + xyshu(k, ip4, ish4))
      xyshd(k, ip2, ish2) = dum; xyshu(k, ip4, ish4) = dum
    end do
  end subroutine tp_displace

! RRX (curved=.false.) / RR (curved=.true.): ports co_pnt_dspl.f90's
! 'RRX'/'RR' branches verbatim (lines 260-366). RRX additionally does an
! axis-constrained pre-step using dx_carry -- see sp_displace_if's
! header comment for the latent-bug preservation this requires. Both
! variants then find the intersection point of the two legs' local
! tangent segments via co_intr_pnt and overlap both legs onto it.
  subroutine rr_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    use mod_constants, only: naddholesmax, nprdbndmax
    class(regular_reflection_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
    include 'paramt.h'
    external co_intr_pnt

    real(wp) :: xc(2), yc(2), xs(2), ys(2), xi, yi
    integer(i4) :: i, ish1, ish2, ip1, ip2

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)

    if (.not. this%curved) then
      xyshu(1, ip1, ish1) = xysh(1, ip1, ish1) + 0.5d+0*eps/dx_carry
      xyshu(2, ip1, ish1) = xysh(2, ip1, ish1)
      xyshd(1, ip1, ish1) = xysh(1, ip1, ish1) - 0.5d+0*eps/dx_carry
      xyshd(2, ip1, ish1) = xysh(2, ip1, ish1)

      xyshu(1, ip2, ish2) = xysh(1, ip2, ish2) + 0.5d+0*eps/dx_carry
      xyshu(2, ip2, ish2) = xysh(2, ip2, ish2)
      xyshd(1, ip2, ish2) = xysh(1, ip2, ish2) - 0.5d+0*eps/dx_carry
      xyshd(2, ip2, ish2) = xysh(2, ip2, ish2)
    end if

    xc(1) = xyshd(1, ip1, ish1)
    yc(1) = xyshd(2, ip1, ish1)
    if (ip1 .eq. 1) then
      xc(2) = xyshd(1, ip1 + 1, ish1)
      yc(2) = xyshd(2, ip1 + 1, ish1)
    else
      xc(2) = xyshd(1, ip1 - 1, ish1)
      yc(2) = xyshd(2, ip1 - 1, ish1)
    end if

    xs(1) = xyshu(1, ip2, ish2)
    ys(1) = xyshu(2, ip2, ish2)
    if (ip2 .eq. 1) then
      xs(2) = xyshu(1, ip2 + 1, ish2)
      ys(2) = xyshu(2, ip2 + 1, ish2)
    else
      xs(2) = xyshu(1, ip2 - 1, ish2)
      ys(2) = xyshu(2, ip2 - 1, ish2)
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    xyshd(1, ip1, ish1) = xi
    xyshd(2, ip1, ish1) = yi
    xyshu(1, ip2, ish2) = xi
    xyshu(2, ip2, ish2) = yi
  end subroutine rr_displace

! QP: ports co_pnt_dspl.f90's 'QP' branch verbatim (lines 368-753) --
! moves all 5 legs along their local tangent by eps, then does 5
! successive co_intr_pnt calls to overlap: incident-1/incident-2 (side
! picked by family), incident-1(other side)/reflected-1-upstream,
! incident-2(other side)/reflected-2-upstream, reflected-1-downstream/
! contact-upstream, reflected-2-downstream/contact-downstream.
  subroutine qp_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    use mod_constants, only: naddholesmax, nprdbndmax
    class(quad_point_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
    include 'paramt.h'
    external co_intr_pnt

    real(wp) :: xc(2), yc(2), xs(2), ys(2), xi, yi
    real(wp) :: tx, ty, dum, f1, f3
    integer(i4) :: i, ish1, ish2, ish3, ish4, ish5, ip1, ip2, ip3, ip4, ip5

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = this%ish(3); i = this%leg(3) - 1; ip3 = 1 + i*(nshockpoints(ish3) - 1)
    ish4 = this%ish(4); i = this%leg(4) - 1; ip4 = 1 + i*(nshockpoints(ish4) - 1)
    ish5 = this%ish(5); i = this%leg(5) - 1; ip5 = 1 + i*(nshockpoints(ish5) - 1)

    f1 = vshnor(1, ip1, ish1)*zroeshu(4, ip1, ish1) -&
    &vshnor(2, ip1, ish1)*zroeshu(3, ip1, ish1)
    f1 = -sign(1.d0, f1)

    f3 = vshnor(1, ip3, ish3)*zroeshu(4, ip3, ish3) -&
    &vshnor(2, ip3, ish3)*zroeshu(3, ip3, ish3)
    f3 = -sign(1.d0, f3)

! move incident shock 1 point
    if (ip1 .eq. 1) then
      tx = xysh(1, ip1, ish1) - xysh(1, ip1 - 1, ish1)
      ty = xysh(2, ip1, ish1) - xysh(2, ip1 - 1, ish1)
    else
      tx = xysh(1, ip1 - 1, ish1) - xysh(1, ip1, ish1)
      ty = xysh(2, ip1 - 1, ish1) - xysh(2, ip1, ish1)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip1, ish1) = xyshu(1, ip1, ish1) + eps*tx
    xyshu(2, ip1, ish1) = xyshu(2, ip1, ish1) + eps*ty
    xyshd(1, ip1, ish1) = xyshd(1, ip1, ish1) + eps*tx
    xyshd(2, ip1, ish1) = xyshd(2, ip1, ish1) + eps*ty

! move reflected shock 1 point
    if (ip2 .eq. 1) then
      tx = xysh(1, ip2, ish2) - xysh(1, ip2 - 1, ish2)
      ty = xysh(2, ip2, ish2) - xysh(2, ip2 - 1, ish2)
    else
      tx = xysh(1, ip2 - 1, ish2) - xysh(1, ip2, ish2)
      ty = xysh(2, ip2 - 1, ish2) - xysh(2, ip2, ish2)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip2, ish2) = xyshu(1, ip2, ish2) + eps*tx
    xyshu(2, ip2, ish2) = xyshu(2, ip2, ish2) + eps*ty
    xyshd(1, ip2, ish2) = xyshd(1, ip2, ish2) + eps*tx
    xyshd(2, ip2, ish2) = xyshd(2, ip2, ish2) + eps*ty

! move incident shock 2 point
    if (ip3 .eq. 1) then
      tx = xysh(1, ip3, ish3) - xysh(1, ip3 - 1, ish3)
      ty = xysh(2, ip3, ish3) - xysh(2, ip3 - 1, ish3)
    else
      tx = xysh(1, ip3 - 1, ish3) - xysh(1, ip3, ish3)
      ty = xysh(2, ip3 - 1, ish3) - xysh(2, ip3, ish3)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip3, ish3) = xyshu(1, ip3, ish3) + eps*tx
    xyshu(2, ip3, ish3) = xyshu(2, ip3, ish3) + eps*ty
    xyshd(1, ip3, ish3) = xyshd(1, ip3, ish3) + eps*tx
    xyshd(2, ip3, ish3) = xyshd(2, ip3, ish3) + eps*ty

! move reflected shock 2 point
    if (ip4 .eq. 1) then
      tx = xysh(1, ip4, ish4) - xysh(1, ip4 - 1, ish4)
      ty = xysh(2, ip4, ish4) - xysh(2, ip4 - 1, ish4)
    else
      tx = xysh(1, ip4 - 1, ish4) - xysh(1, ip4, ish4)
      ty = xysh(2, ip4 - 1, ish4) - xysh(2, ip4, ish4)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip4, ish4) = xyshu(1, ip4, ish4) + eps*tx
    xyshu(2, ip4, ish4) = xyshu(2, ip4, ish4) + eps*ty
    xyshd(1, ip4, ish4) = xyshd(1, ip4, ish4) + eps*tx
    xyshd(2, ip4, ish4) = xyshd(2, ip4, ish4) + eps*ty

! move contact discontinuity point
    if (ip5 .eq. 1) then
      tx = xysh(1, ip5, ish5) - xysh(1, ip5 - 1, ish5)
      ty = xysh(2, ip5, ish5) - xysh(2, ip5 - 1, ish5)
    else
      tx = xysh(1, ip5 - 1, ish5) - xysh(1, ip5, ish5)
      ty = xysh(2, ip5 - 1, ish5) - xysh(2, ip5, ish5)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip5, ish5) = xyshu(1, ip5, ish5) + eps*tx
    xyshu(2, ip5, ish5) = xyshu(2, ip5, ish5) + eps*ty
    xyshd(1, ip5, ish5) = xyshd(1, ip5, ish5) + eps*tx
    xyshd(2, ip5, ish5) = xyshd(2, ip5, ish5) + eps*ty

! overlap incident shock 1 with incident shock 2 (side picked by family)
    if (f1 .gt. 0.d0) then
      xc(1) = xyshu(1, ip1, ish1); yc(1) = xyshu(2, ip1, ish1)
      if (ip1 .eq. 1) then
        xc(2) = xyshu(1, ip1 + 1, ish1); yc(2) = xyshu(2, ip1 + 1, ish1)
      else
        xc(2) = xyshu(1, ip1 - 1, ish1); yc(2) = xyshu(2, ip1 - 1, ish1)
      end if
    else
      xc(1) = xyshd(1, ip1, ish1); yc(1) = xyshd(2, ip1, ish1)
      if (ip1 .eq. 1) then
        xc(2) = xyshd(1, ip1 + 1, ish1); yc(2) = xyshd(2, ip1 + 1, ish1)
      else
        xc(2) = xyshd(1, ip1 - 1, ish1); yc(2) = xyshd(2, ip1 - 1, ish1)
      end if
    end if

    if (f3 .lt. 0.d0) then
      xs(1) = xyshu(1, ip3, ish3); ys(1) = xyshu(2, ip3, ish3)
      if (ip3 .eq. 1) then
        xs(2) = xyshu(1, ip3 + 1, ish3); ys(2) = xyshu(2, ip3 + 1, ish3)
      else
        xs(2) = xyshu(1, ip3 - 1, ish3); ys(2) = xyshu(2, ip3 - 1, ish3)
      end if
    else
      xs(1) = xyshd(1, ip3, ish3); ys(1) = xyshd(2, ip3, ish3)
      if (ip3 .eq. 1) then
        xs(2) = xyshd(1, ip3 + 1, ish3); ys(2) = xyshd(2, ip3 + 1, ish3)
      else
        xs(2) = xyshd(1, ip3 - 1, ish3); ys(2) = xyshd(2, ip3 - 1, ish3)
      end if
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    if (f1 .gt. 0.d0) then
      xyshu(1, ip1, ish1) = xi; xyshu(2, ip1, ish1) = yi
    else
      xyshd(1, ip1, ish1) = xi; xyshd(2, ip1, ish1) = yi
    end if
    if (f3 .lt. 0.d0) then
      xyshu(1, ip3, ish3) = xi; xyshu(2, ip3, ish3) = yi
    else
      xyshd(1, ip3, ish3) = xi; xyshd(2, ip3, ish3) = yi
    end if

! overlap incident shock 1 (other side) with reflected shock 1 upstream
    if (f1 .gt. 0.d0) then
      xc(1) = xyshd(1, ip1, ish1); yc(1) = xyshd(2, ip1, ish1)
      if (ip1 .eq. 1) then
        xc(2) = xyshd(1, ip1 + 1, ish1); yc(2) = xyshd(2, ip1 + 1, ish1)
      else
        xc(2) = xyshd(1, ip1 - 1, ish1); yc(2) = xyshd(2, ip1 - 1, ish1)
      end if
    else
      xc(1) = xyshu(1, ip1, ish1); yc(1) = xyshu(2, ip1, ish1)
      if (ip1 .eq. 1) then
        xc(2) = xyshu(1, ip1 + 1, ish1); yc(2) = xyshu(2, ip1 + 1, ish1)
      else
        xc(2) = xyshu(1, ip1 - 1, ish1); yc(2) = xyshu(2, ip1 - 1, ish1)
      end if
    end if

    xs(1) = xyshu(1, ip2, ish2); ys(1) = xyshu(2, ip2, ish2)
    if (ip2 .eq. 1) then
      xs(2) = xyshu(1, ip2 + 1, ish2); ys(2) = xyshu(2, ip2 + 1, ish2)
    else
      xs(2) = xyshu(1, ip2 - 1, ish2); ys(2) = xyshu(2, ip2 - 1, ish2)
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    if (f1 .gt. 0.d0) then
      xyshd(1, ip1, ish1) = xi; xyshd(2, ip1, ish1) = yi
    else
      xyshu(1, ip1, ish1) = xi; xyshu(2, ip1, ish1) = yi
    end if
    xyshu(1, ip2, ish2) = xi; xyshu(2, ip2, ish2) = yi

! overlap incident shock 2 (other side) with reflected shock 2 upstream
    if (f3 .gt. 0.d+0) then
      xc(1) = xyshu(1, ip3, ish3); yc(1) = xyshu(2, ip3, ish3)
      if (ip3 .eq. 1) then
        xc(2) = xyshu(1, ip3 + 1, ish3); yc(2) = xyshu(2, ip3 + 1, ish3)
      else
        xc(2) = xyshu(1, ip3 - 1, ish3); yc(2) = xyshu(2, ip3 - 1, ish3)
      end if
    else
      xc(1) = xyshd(1, ip3, ish3); yc(1) = xyshd(2, ip3, ish3)
      if (ip3 .eq. 1) then
        xc(2) = xyshd(1, ip3 + 1, ish3); yc(2) = xyshd(2, ip3 + 1, ish3)
      else
        xc(2) = xyshd(1, ip3 - 1, ish3); yc(2) = xyshd(2, ip3 - 1, ish3)
      end if
    end if

    xs(1) = xyshu(1, ip4, ish4); ys(1) = xyshu(2, ip4, ish4)
    if (ip4 .eq. 1) then
      xs(2) = xyshu(1, ip4 + 1, ish4); ys(2) = xyshu(2, ip4 + 1, ish4)
    else
      xs(2) = xyshu(1, ip4 - 1, ish4); ys(2) = xyshu(2, ip4 - 1, ish4)
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    if (f3 .lt. 0.d0) then
      xyshd(1, ip3, ish3) = xi; xyshd(2, ip3, ish3) = yi
    else
      xyshu(1, ip3, ish3) = xi; xyshu(2, ip3, ish3) = yi
    end if
    xyshu(1, ip4, ish4) = xi; xyshu(2, ip4, ish4) = yi

! overlap reflected shock 1 downstream with contact discontinuity upstream
    xc(1) = xyshd(1, ip2, ish2); yc(1) = xyshd(2, ip2, ish2)
    if (ip2 .eq. 1) then
      xc(2) = xyshd(1, ip2 + 1, ish2); yc(2) = xyshd(2, ip2 + 1, ish2)
    else
      xc(2) = xyshd(1, ip2 - 1, ish2); yc(2) = xyshd(2, ip2 - 1, ish2)
    end if

    xs(1) = xyshu(1, ip5, ish5); ys(1) = xyshu(2, ip5, ish5)
    if (ip5 .eq. 1) then
      xs(2) = xyshu(1, ip5 + 1, ish5); ys(2) = xyshu(2, ip5 + 1, ish5)
    else
      xs(2) = xyshu(1, ip5 - 1, ish5); ys(2) = xyshu(2, ip5 - 1, ish5)
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    xyshd(1, ip2, ish2) = xi; xyshd(2, ip2, ish2) = yi
    xyshu(1, ip5, ish5) = xi; xyshu(2, ip5, ish5) = yi

! overlap reflected shock 2 downstream with contact discontinuity downstream
    xc(1) = xyshd(1, ip4, ish4); yc(1) = xyshd(2, ip4, ish4)
    if (ip4 .eq. 1) then
      xc(2) = xyshd(1, ip4 + 1, ish4); yc(2) = xyshd(2, ip4 + 1, ish4)
    else
      xc(2) = xyshd(1, ip4 - 1, ish4); yc(2) = xyshd(2, ip4 - 1, ish4)
    end if

    xs(1) = xyshd(1, ip5, ish5); ys(1) = xyshd(2, ip5, ish5)
    if (ip5 .eq. 1) then
      xs(2) = xyshd(1, ip5 + 1, ish5); ys(2) = xyshd(2, ip5 + 1, ish5)
    else
      xs(2) = xyshd(1, ip5 - 1, ish5); ys(2) = xyshd(2, ip5 - 1, ish5)
    end if

    call co_intr_pnt(xi, yi, xc, yc, xs, ys)

    xyshd(1, ip4, ish4) = xi; xyshd(2, ip4, ish4) = yi
    xyshd(1, ip5, ish5) = xi; xyshd(2, ip5, ish5) = yi
  end subroutine qp_displace

! TE: ports co_pnt_dspl.f90's 'TE' branch verbatim (lines 756-802). Uses
! its OWN leg-index roles (leg1=shock1, leg2=contact discontinuity,
! leg3=shock2) -- distinct from te_solve_state's leg roles, since this
! is independent legacy code that happens to share the same nshe=3
! shinspps records. The `if (ip1 .eq. 1)` guard on a tangent computed
! from ip2/ish2 is exactly what the original does (not a transcription
! slip here) -- preserved as-is.
  subroutine te_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    use mod_constants, only: naddholesmax, nprdbndmax
    class(trailing_edge_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry
    include 'paramt.h'

    real(wp) :: tx, ty, dum
    integer(i4) :: i, ish1, ish2, ish3, ip1, ip2, ip3

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)
    ish3 = this%ish(3); i = this%leg(3) - 1; ip3 = 1 + i*(nshockpoints(ish3) - 1)

! move contact discontinuity point
    if (ip1 .eq. 1) then
      tx = xysh(1, ip2, ish2) - xysh(1, ip2 - 1, ish2)
      ty = xysh(2, ip2, ish2) - xysh(2, ip2 - 1, ish2)
    else
      tx = xysh(1, ip2 - 1, ish2) - xysh(1, ip2, ish2)
      ty = xysh(2, ip2 - 1, ish2) - xysh(2, ip2, ish2)
    end if
    dum = sqrt(tx**2 + ty**2)
    tx = tx/dum; ty = ty/dum
    xyshu(1, ip2, ish2) = xyshu(1, ip2, ish2) + eps*tx
    xyshu(2, ip2, ish2) = xyshu(2, ip2, ish2) + eps*ty
    xyshd(1, ip2, ish2) = xyshd(1, ip2, ish2) + eps*tx
    xyshd(2, ip2, ish2) = xyshd(2, ip2, ish2) + eps*ty

    xyshd(1, ip1, ish1) = xyshd(1, ip2, ish2)
    xyshd(2, ip1, ish1) = xyshd(2, ip2, ish2)

    xyshd(1, ip3, ish3) = xyshu(1, ip2, ish2)
    xyshd(2, ip3, ish3) = xyshu(2, ip2, ish2)
  end subroutine te_displace

! EP: ports co_pnt_dspl.f90's 'EP' branch verbatim (lines 804-814) --
! simply overlaps up/downstream with the un-displaced point.
  subroutine ep_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    class(end_point_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry

    integer(i4) :: i, ish, ip

    ish = this%ish(1); i = this%leg(1) - 1
    ip = 1 + i*(nshockpoints(ish) - 1)

    xyshd(1, ip, ish) = xysh(1, ip, ish)
    xyshd(2, ip, ish) = xysh(2, ip, ish)
    xyshu(1, ip, ish) = xysh(1, ip, ish)
    xyshu(2, ip, ish) = xysh(2, ip, ish)
  end subroutine ep_displace

! SP: ports co_pnt_dspl.f90's 'SP' branch verbatim (lines 835-845) --
! byte-identical to ep_displace's body in this chain (the original
! duplicates this same 4-line overlap between its EP and SP branches;
! kept duplicated here too, matching the source rather than merging).
  subroutine sonic_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    class(start_point_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry

    integer(i4) :: i, ish, ip

    ish = this%ish(1); i = this%leg(1) - 1
    ip = 1 + i*(nshockpoints(ish) - 1)

    xyshd(1, ip, ish) = xysh(1, ip, ish)
    xyshd(2, ip, ish) = xysh(2, ip, ish)
    xyshu(1, ip, ish) = xysh(1, ip, ish)
    xyshu(2, ip, ish) = xysh(2, ip, ish)
  end subroutine sonic_displace

! C (periodic=.false.) / PC (periodic=.true.): ports co_pnt_dspl.f90's
! 'C' branch verbatim (lines 816-833); 'PC' is a confirmed no-op in this
! chain (empty branch in the original), so periodic returns immediately.
  subroutine conn_displace(this, xysh, xyshu, xyshd, vshnor, zroeshu,&
  &nshockpoints, dx_carry, dy_carry)
    class(connection_t), intent(inout) :: this
    real(wp), intent(in) :: xysh(:, :, :)
    real(wp), intent(inout) :: xyshu(:, :, :), xyshd(:, :, :)
    real(wp), intent(in) :: vshnor(:, :, :), zroeshu(:, :, :)
    integer(i4), intent(in) :: nshockpoints(:)
    real(wp), intent(inout) :: dx_carry, dy_carry

    integer(i4) :: i, ish1, ish2, ip1, ip2

    if (this%periodic) return

    ish1 = this%ish(1); i = this%leg(1) - 1; ip1 = 1 + i*(nshockpoints(ish1) - 1)
    ish2 = this%ish(2); i = this%leg(2) - 1; ip2 = 1 + i*(nshockpoints(ish2) - 1)

    xyshd(1, ip2, ish2) = xyshd(1, ip1, ish1)
    xyshd(2, ip2, ish2) = xyshd(2, ip1, ish1)
    xyshu(1, ip2, ish2) = xyshu(1, ip1, ish1)
    xyshu(2, ip2, ish2) = xyshu(2, ip1, ish1)
  end subroutine conn_displace

! Shared default for remesh_boundary/relocate/correct_normal (still
! their increment-2 minimal (this)-only signature -- increments 6-8
! give each its own real signature the same way displace just did,
! replacing this default with a dedicated sp_..._noop where needed).
! Literally does nothing, matching today's empty if/elseif arms (e.g.
! fx_msh_sps.f90's TP/QP/EP/C/SP branches, fx_dps_loc.f90's TP/QP/RRX/
! EP/SP/C/TE branches, ...). This one CAN be shared as-is across all
! types: it's bound directly on special_point_t itself (not overriding
! a deferred binding in an extension), so its passed-object dummy
! legitimately is the abstract base type.
  subroutine sp_noop(this)
    class(special_point_t), intent(inout) :: this
  end subroutine sp_noop

  function new_triple_point_t() result(sp)
    type(triple_point_t) :: sp
    sp%code = 'TP'
    sp%nshe = 4
  end function new_triple_point_t

  function new_quad_point_t() result(sp)
    type(quad_point_t) :: sp
    sp%code = 'QP'
    sp%nshe = 5
  end function new_quad_point_t

  function new_trailing_edge_t() result(sp)
    type(trailing_edge_t) :: sp
    sp%code = 'TE'
    sp%nshe = 3
  end function new_trailing_edge_t

  function new_regular_reflection_t(code) result(sp)
    character(len=*), intent(in) :: code
    type(regular_reflection_t) :: sp
    select case (trim(code))
    case ('RRX')
      sp%curved = .false.
    case ('RR')
      sp%curved = .true.
    case default
      error stop 'new_regular_reflection_t: unrecognized code '//trim(code)
    end select
    sp%code = code
    sp%nshe = 2
  end function new_regular_reflection_t

  function new_wall_float_t(code) result(sp)
    character(len=*), intent(in) :: code
    type(wall_float_t) :: sp
    select case (trim(code))
    case ('IPX')
      sp%axis = 'X'; sp%restore_only = .true.
    case ('IPY')
      sp%axis = 'Y'; sp%restore_only = .true.
    case ('OPX')
      sp%axis = 'X'; sp%restore_only = .false.
    case ('OPY')
      sp%axis = 'Y'; sp%restore_only = .false.
    case ('WPNRX')
      sp%axis = 'X'; sp%restore_only = .false.
    case ('WPNRY')
      sp%axis = 'Y'; sp%restore_only = .false.
    case default
      error stop 'new_wall_float_t: unrecognized code '//trim(code)
    end select
    sp%code = code
    sp%nshe = 1
  end function new_wall_float_t

  function new_floating_wall_point_t() result(sp)
    type(floating_wall_point_t) :: sp
    sp%code = 'FWP'
    sp%nshe = 1
  end function new_floating_wall_point_t

  function new_end_point_t() result(sp)
    type(end_point_t) :: sp
    sp%code = 'EP'
    sp%nshe = 1
  end function new_end_point_t

  function new_start_point_t() result(sp)
    type(start_point_t) :: sp
    sp%code = 'SP'
    sp%nshe = 1
  end function new_start_point_t

  function new_connection_t(code) result(sp)
    character(len=*), intent(in) :: code
    type(connection_t) :: sp
    select case (trim(code))
    case ('C')
      sp%periodic = .false.
    case ('PC')
      sp%periodic = .true.
    case default
      error stop 'new_connection_t: unrecognized code '//trim(code)
    end select
    sp%code = code
    sp%nshe = 2
  end function new_connection_t

end module mod_special_point
