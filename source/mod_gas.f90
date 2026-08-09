module mod_gas
! Issue #21 (E1, ROADMAP.md Part III): the gas model this project has
! never had -- UnDiFi-2D is single-gas today, one scalar GA read from
! input.dat by re_inp_data.f90 into COMMON/PARAMT/, referenced directly
! at roughly 200 sites (fx_state_dps.f90, co_utp.f90, co_uqp.f90,
! co_urr.f90, co_shock.f90, ...). None of those ~200 sites are touched
! by this increment -- see mod_gas_registry.f90's own header for what
! is and isn't wired up yet.
!
! gas_t is deliberately minimal: gamma is everything co_shock_2gas.f90/
! co_interface.f90 (issue #23) need today. R (specific gas constant) is
! carried for a future caloric/thermal model (issue text: "a caloric
! model" is explicitly a *later* addition) -- present now so
! gas_registry_t's array shape doesn't need to change again once that
! lands, not because anything reads it yet.
  use mod_kinds, only: wp
  implicit none(type, external)
  private
  public :: gas_t

  type :: gas_t
    character(len=:), allocatable :: name
    real(wp) :: gamma = 1.4_wp
    real(wp) :: r_gas = 287.0_wp ! J/(kg K), air's value; unused until a caloric model exists
  end type gas_t

end module mod_gas

module mod_gas_registry
! A small append-only collection of gas_t, indexed by integer -- the
! same integer discontinuity_t%gas_up/%gas_dn already carry (Phase 3.1
! scaffolding, mod_discontinuity.f90, default 1/1 = "the one gas" for
! every case that never calls add_gas). Index 1 is reserved for the
! default/single gas so existing single-gas behavior is exactly index-1
! lookups everywhere, once (if) this registry is ever wired into a real
! call site -- not done in this increment. This module has no consumer
! yet, same status mod_discontinuity.f90's own `kind`/`gas_up`/`gas_dn`
! fields were left in after Phase 3.1: real storage for Phase 3.2-style
! wiring to start binding against, not itself that wiring.
  use mod_kinds, only: wp, i4
  use mod_gas, only: gas_t
  implicit none(type, external)
  private
  public :: gas_registry_t

  type :: gas_registry_t
    type(gas_t), allocatable :: gases(:)
  contains
    procedure :: add => gas_registry_add
    procedure :: gamma => gas_registry_gamma
    procedure :: count => gas_registry_count
  end type gas_registry_t

contains

! Appends a new gas, returns its 1-based index (the value to store in
! discontinuity_t%gas_up/%gas_dn). First call on an unallocated registry
! creates it; every call is a realloc-and-copy (this registry is built
! once at case setup, not in any hot loop -- simplicity over an
! amortized-growth scheme here is deliberate).
  function gas_registry_add(this, name, gamma, r_gas) result(idx)
    class(gas_registry_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(wp), intent(in) :: gamma
    real(wp), intent(in), optional :: r_gas
    integer(i4) :: idx

    type(gas_t), allocatable :: tmp(:)
    integer(i4) :: n

    if (.not. allocated(this%gases)) then
      allocate (this%gases(0))
    end if

    n = size(this%gases)
    allocate (tmp(n + 1))
    tmp(1:n) = this%gases
    tmp(n + 1)%name = name
    tmp(n + 1)%gamma = gamma
    if (present(r_gas)) tmp(n + 1)%r_gas = r_gas

    call move_alloc(tmp, this%gases)
    idx = n + 1
  end function gas_registry_add

  real(wp) function gas_registry_gamma(this, idx) result(gam)
    class(gas_registry_t), intent(in) :: this
    integer(i4), intent(in) :: idx
    gam = this%gases(idx)%gamma
  end function gas_registry_gamma

  integer(i4) function gas_registry_count(this) result(n)
    class(gas_registry_t), intent(in) :: this
    n = 0
    if (allocated(this%gases)) n = size(this%gases)
  end function gas_registry_count

end module mod_gas_registry
