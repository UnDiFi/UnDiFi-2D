module mod_special_point_registry
! Phase 3.2 increment 3 (merged 3.2+3.3+3.4, ROADMAP.md #14): the
! factory (make_special_point) that replaces every dispatch chain's
! typespecpoints if/elseif with a single type_code -> concrete-type
! lookup, plus sp_read/sp_write, the shared serialization step that
! replaces re_sdw_info.f90/wrt_sdw_info.f90's own per-code nshe/ispclr
! branching (those two files independently hardcoded the same 16-way
! nshe table with no cross-check between them -- make_special_point's
! constructors are now the single source of truth for nshe).
!
! sp_read/sp_write take the caller's existing shinspps(:,:,isppnts)/
! ispclr(:,isppnts) slices directly (not a new intermediate
! representation) -- they read/write exactly the same npshmax-agnostic
! per-point legacy arrays the dispatch chains already pass around, so
! re_sdw_info.f90/wrt_sdw_info.f90 (and, in later increments, the other
! four chains) need no data-model change beyond calling into this
! module instead of repeating their own if/elseif.

  use mod_kinds, only: wp, i4
  use mod_special_point, only: special_point_t, triple_point_t,&
  &quad_point_t, trailing_edge_t, regular_reflection_t, wall_float_t,&
  &floating_wall_point_t, end_point_t, start_point_t, connection_t,&
  &new_triple_point_t, new_quad_point_t, new_trailing_edge_t,&
  &new_regular_reflection_t, new_wall_float_t,&
  &new_floating_wall_point_t, new_end_point_t, new_start_point_t,&
  &new_connection_t
  implicit none(type, external)
  private
  public :: make_special_point, sp_read, sp_write, sp_unpack

contains

! Allocates sp to the concrete type matching `code`, with nshe and any
! type-specific discriminating fields (curved/axis/restore_only/
! periodic) already set by that type's constructor. Preserves the
! legacy dispatch chains' fatal-stop-on-unrecognized-code behavior
! (re_sdw_info.f90/wrt_sdw_info.f90/fx_msh_sps.f90/fx_dps_loc.f90/
! co_norm.f90 all `stop` on an unhandled typespecpoints value) -- the
! exact message text is not shared verbatim across all former call
! sites (each legacy file wrote a slightly different string to its own
! log unit before stopping), which is an acceptable simplification
! since no valid production run ever reaches this branch and message
! text isn't part of the regression checksum.
  subroutine make_special_point(code, sp)
    character(len=*), intent(in) :: code
    class(special_point_t), allocatable, intent(out) :: sp

    select case (trim(code))
    case ('TP')
      sp = new_triple_point_t()
    case ('QP')
      sp = new_quad_point_t()
    case ('TE')
      sp = new_trailing_edge_t()
    case ('RRX', 'RR')
      sp = new_regular_reflection_t(code)
    case ('WPNRX', 'WPNRY', 'IPX', 'IPY', 'OPX', 'OPY')
      sp = new_wall_float_t(code)
    case ('FWP')
      sp = new_floating_wall_point_t()
    case ('EP')
      sp = new_end_point_t()
    case ('SP')
      sp = new_start_point_t()
    case ('C', 'PC')
      sp = new_connection_t(code)
    case default
      error stop 'make_special_point: condition not implemented: '//trim(code)
    end select
  end subroutine make_special_point

! True for the 3 concrete-type/variant combinations whose sh00.dat/
! sh99.dat record carries a per-leg boundary color (ispclr) alongside
! the (ish,leg) pair -- RR (not RRX), FWP (always), PC (not C). Matches
! re_sdw_info.f90/wrt_sdw_info.f90's read/write statements exactly.
  logical function sp_has_iclr(sp) result(has_iclr)
    class(special_point_t), intent(in) :: sp
    has_iclr = .false.
    select type (sp)
    type is (regular_reflection_t)
      has_iclr = sp%curved
    type is (floating_wall_point_t)
      has_iclr = .true.
    type is (connection_t)
      has_iclr = sp%periodic
    end select
  end function sp_has_iclr

! Reads sp%nshe (ish,leg[,iclr]) records from `unit`, echoing each to
! `log_unit` -- the exact per-point body of re_sdw_info.f90's dispatch
! chain, now driven by sp%nshe/sp_has_iclr instead of a hardcoded
! per-branch nshe and a repeated ispclr read in 3 of the 16 branches.
  subroutine sp_read(unit, log_unit, sp, shinspps, ispclr)
    integer(i4), intent(in) :: unit, log_unit
    class(special_point_t), intent(in) :: sp
    integer(i4), intent(out) :: shinspps(:, :) ! (2, >=nshe) slice: shinspps(:,:,isppnts)
    integer(i4), intent(inout) :: ispclr(:) ! (>=nshe) slice: ispclr(:,isppnts)
    integer(i4) :: k
    logical :: has_iclr

    has_iclr = sp_has_iclr(sp)
    do k = 1, sp%nshe
      if (has_iclr) then
        read (unit, *) shinspps(1, k), shinspps(2, k), ispclr(k)
        write (log_unit, *) shinspps(1, k), shinspps(2, k), ispclr(k)
      else
        read (unit, *) shinspps(1, k), shinspps(2, k)
        write (log_unit, *) shinspps(1, k), shinspps(2, k)
      end if
    end do
  end subroutine sp_read

! Mirror of sp_read for wrt_sdw_info.f90's write side.
  subroutine sp_write(unit, log_unit, sp, shinspps, ispclr)
    integer(i4), intent(in) :: unit, log_unit
    class(special_point_t), intent(in) :: sp
    integer(i4), intent(in) :: shinspps(:, :) ! (2, >=nshe) slice: shinspps(:,:,isppnts)
    integer(i4), intent(in) :: ispclr(:) ! (>=nshe) slice: ispclr(:,isppnts)
    integer(i4) :: k
    logical :: has_iclr

    has_iclr = sp_has_iclr(sp)
    do k = 1, sp%nshe
      if (has_iclr) then
        write (unit, *) shinspps(1, k), shinspps(2, k), ispclr(k)
        write (log_unit, *) shinspps(1, k), shinspps(2, k), ispclr(k)
      else
        write (unit, *) shinspps(1, k), shinspps(2, k)
        write (log_unit, *) shinspps(1, k), shinspps(2, k)
      end if
    end do
  end subroutine sp_write

! Phase 3.2 increment 4: populates sp%ish/sp%leg (and, for the 3
! iclr-carrying variants, sp%iclr, when a well-shaped ispclr slice is
! available) from the shinspps(:,:,isppnts)/ispclr(:,isppnts) slices
! already read once by re_sdw_info.f90 at startup -- the behavior
! chains (fx_state_dps.f90 and, in later increments, co_pnt_dspl/
! fx_msh_sps/fx_dps_loc/co_norm) call this right after
! make_special_point, every dispatch, since each special_point_t is
! built fresh per point per call (see mod_special_point.f90's header
! for why that's the right shape here, unlike mod_discontinuity.f90's
! persistent aliased storage).
!
! ispclr is OPTIONAL because fx_state_dps.f90 itself declares its own
! ispclr dummy argument rank-1 (ispclr(*)) rather than the rank-2
! (5,*) shape every other chain uses -- a real, confirmed pre-existing
! mismatch (bug #1 in the merged-3.2 plan) reachable only via 'RR',
! which no regression fixture exercises. Rather than paper over that
! bug by reshaping a rank-1 array into something sp_unpack can slice
! per-leg, fx_state_dps.f90 omits ispclr here entirely and instead
! passes the raw ispclr(isppnts) scalar straight through to
! solve_state's ispclr_flat argument, reproducing the legacy flat
! single-subscript read bit-for-bit. Callers with a properly rank-2
! ispclr (increments 5-8) pass the real per-point slice instead.
  subroutine sp_unpack(sp, shinspps, ispclr)
    class(special_point_t), intent(inout) :: sp
    integer(i4), intent(in) :: shinspps(:, :) ! (2, >=nshe) slice: shinspps(:,:,isppnts)
    integer(i4), intent(in), optional :: ispclr(:) ! (>=nshe) slice: ispclr(:,isppnts)

    allocate (sp%ish(sp%nshe), sp%leg(sp%nshe))
    sp%ish(1:sp%nshe) = shinspps(1, 1:sp%nshe)
    sp%leg(1:sp%nshe) = shinspps(2, 1:sp%nshe)

    if (present(ispclr)) then
      select type (sp)
      type is (regular_reflection_t)
        if (sp%curved) sp%iclr(1:sp%nshe) = ispclr(1:sp%nshe)
      type is (floating_wall_point_t)
        sp%iclr = ispclr(1)
      type is (connection_t)
        if (sp%periodic) sp%iclr(1:sp%nshe) = ispclr(1:sp%nshe)
      end select
    end if
  end subroutine sp_unpack

end module mod_special_point_registry
