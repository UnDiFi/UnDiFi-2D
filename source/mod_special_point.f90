module mod_special_point
! Phase 3.2 increment 2 (merged 3.2+3.3+3.4, ROADMAP.md #14): type
! hierarchy skeleton for the special-point object model that will
! eventually replace all six typespecpoints-keyed dispatch chains
! (fx_state_dps, co_pnt_dspl, fx_msh_sps, fx_dps_loc, co_norm, and the
! re_sdw_info/wrt_sdw_info serialization pair). See
! /home/lk/.claude/plans/proud-jingling-donut.md's "Merged 3.2" section
! for the full design (type-hierarchy rationale, per-code nshe table,
! which of the five per-chain behaviors each type overrides vs. leaves
! at the shared no-op default, and the six known latent bugs that must
! be preserved bit-for-bit when each chain is actually ported).
!
! This increment defines ONLY the type hierarchy, the 9 concrete types'
! constructors (which set nshe and any type-specific discriminating
! fields -- the one piece of "real" behavior this increment delivers),
! and stub bodies for every behavior method. Nothing calls any of these
! types yet; increments 3-8 wire up the registry and port each dispatch
! chain's real logic into these stubs one at a time.
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
! special_point_t built fresh via the registry (increment 3), used, and
! discarded, matches the existing call pattern with zero pointer-bridge
! complexity.

  use mod_kinds, only: wp, i4
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
    procedure :: displace => sp_noop ! co_pnt_dspl.f90 -- default no-op (increment 5)
    procedure :: relocate => sp_noop ! fx_dps_loc.f90  -- default no-op (increment 7)
    procedure :: correct_normal => sp_noop ! co_norm.f90     -- default no-op (increment 8)
  end type special_point_t

  abstract interface
    subroutine sp_solve_if(this)
      import :: special_point_t
      class(special_point_t), intent(inout) :: this
    end subroutine sp_solve_if
  end interface

! TP -- triple point (internal), nshe=4. Newton-solved via co_utp today.
  type, extends(special_point_t) :: triple_point_t
  contains
    procedure :: solve_state => tp_solve_state
  end type triple_point_t

! QP -- quadruple point (internal), nshe=5. Newton-solved via co_uqp today.
  type, extends(special_point_t) :: quad_point_t
  contains
    procedure :: solve_state => qp_solve_state
  end type quad_point_t

! TE -- trailing edge point, nshe=3. Also Newton-solved via co_uqp today
! (fx_state_dps.f90 reuses the quadruple-point solver for TE, with two
! of the five states defaulted to placeholder values).
  type, extends(special_point_t) :: trailing_edge_t
  contains
    procedure :: solve_state => te_solve_state
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
  end type end_point_t

! SP -- start/sonic point (characteristic coalescence), nshe=1.
  type, extends(special_point_t) :: start_point_t
  contains
    procedure :: solve_state => sonic_solve_state
  end type start_point_t

! C (periodic=.false.) / PC (periodic=.true.) -- connection between two
! shocks, nshe=2. iclr is only meaningful when periodic.
  type, extends(special_point_t) :: connection_t
    logical :: periodic = .false.
    integer(i4) :: iclr(2) = 0
  contains
    procedure :: solve_state => conn_solve_state
  end type connection_t

contains

! Shared message builder behind every concrete type's solve_state stub
! below -- this increment ports the type hierarchy only; increment 4
! replaces each type's binding with real fx_state_dps.f90 logic
! (co_utp/co_uqp/co_urr for the three that Newton-solve, direct algebra
! for the rest). A single shared subroutine can't implement a deferred
! binding for every type directly: gfortran requires the passed-object
! dummy of a deferred binding's implementation to be of the exact
! extending type, not the abstract base -- hence one small per-type
! wrapper below, all funneling into this one message.
  subroutine sp_report_not_ported(code)
    character(len=*), intent(in) :: code
    error stop 'special_point_t%solve_state: '//trim(code)//&
    &' not yet ported (Phase 3.2 increment 4)'
  end subroutine sp_report_not_ported

  subroutine tp_solve_state(this)
    class(triple_point_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine tp_solve_state

  subroutine qp_solve_state(this)
    class(quad_point_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine qp_solve_state

  subroutine te_solve_state(this)
    class(trailing_edge_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine te_solve_state

  subroutine rr_solve_state(this)
    class(regular_reflection_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine rr_solve_state

  subroutine wf_solve_state(this)
    class(wall_float_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine wf_solve_state

  subroutine fwp_solve_state(this)
    class(floating_wall_point_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine fwp_solve_state

  subroutine ep_solve_state(this)
    class(end_point_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine ep_solve_state

  subroutine sonic_solve_state(this)
    class(start_point_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine sonic_solve_state

  subroutine conn_solve_state(this)
    class(connection_t), intent(inout) :: this
    call sp_report_not_ported(this%code)
  end subroutine conn_solve_state

! Shared default for the four behaviors most codes leave untouched in
! their corresponding legacy chain -- literally does nothing, matching
! today's empty if/elseif arms (e.g. fx_msh_sps.f90's TP/QP/EP/C/SP
! branches, fx_dps_loc.f90's TP/QP/RRX/EP/SP/C/TE branches, ...). This
! one CAN be shared as-is across all types: it's bound directly on
! special_point_t itself (not overriding a deferred binding in an
! extension), so its passed-object dummy legitimately is the abstract
! base type.
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
