module test_mod_gas
! Issue #21 (E1): mod_gas/gas_registry -- the foundational type_gas_t +
! gas_registry_t this issue's other tasks (discontinuity_t%gas_up/
! %gas_dn wiring, the ~200-site GA/GM1 replacement) depend on. Neither
! of those consumers exist yet (see mod_gas.f90's own header); this
! suite is purely about the registry's own, small contract.
  use mod_kinds, only: wp, i4
  use mod_gas_registry, only: gas_registry_t
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none(type, external)
  private
  public :: collect_mod_gas

contains

  subroutine collect_mod_gas(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
    &new_unittest("empty_registry_has_zero_count", test_empty),&
    &new_unittest("add_returns_1based_sequential_index", test_add_index),&
    &new_unittest("gamma_lookup_matches_what_was_added", test_gamma_lookup)&
    &]
  end subroutine collect_mod_gas

  subroutine test_empty(error)
    type(error_type), allocatable, intent(out) :: error
    type(gas_registry_t) :: reg

    call check(error, reg%count(), 0_i4)
    if (allocated(error)) return
  end subroutine test_empty

  subroutine test_add_index(error)
    type(error_type), allocatable, intent(out) :: error
    type(gas_registry_t) :: reg
    integer(i4) :: idx1, idx2, idx3

    idx1 = reg%add("air", 1.4_wp)
    idx2 = reg%add("helium", 1.667_wp)
    idx3 = reg%add("CO2", 1.289_wp)

    call check(error, idx1, 1_i4)
    if (allocated(error)) return
    call check(error, idx2, 2_i4)
    if (allocated(error)) return
    call check(error, idx3, 3_i4)
    if (allocated(error)) return
    call check(error, reg%count(), 3_i4)
    if (allocated(error)) return
  end subroutine test_add_index

  subroutine test_gamma_lookup(error)
    type(error_type), allocatable, intent(out) :: error
    type(gas_registry_t) :: reg
    integer(i4) :: idx_air, idx_he

    idx_air = reg%add("air", 1.4_wp)
    idx_he = reg%add("helium", 1.667_wp, r_gas=2077.0_wp)

    call check(error, reg%gamma(idx_air), 1.4_wp, thr=1.0e-14_wp)
    if (allocated(error)) return
    call check(error, reg%gamma(idx_he), 1.667_wp, thr=1.0e-14_wp)
    if (allocated(error)) return
    call check(error, reg%gases(idx_he)%r_gas, 2077.0_wp, thr=1.0e-14_wp)
    if (allocated(error)) return
    call check(error, reg%gases(idx_air)%name == "air")
    if (allocated(error)) return
  end subroutine test_gamma_lookup

end module test_mod_gas
