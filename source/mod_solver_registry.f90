module mod_solver_registry
! Phase 3.7 (ROADMAP.md #14, issue #15): the flow_solver_t factory,
! separated from mod_solver_iface.f90 itself for the same reason
! mod_special_point_registry.f90 is separate from mod_special_point.f90
! -- the factory needs to know about every concrete type (eulfs_t/
! su2_t/neo_t), while the abstract base module must NOT know about any
! of them (a concrete type's own module already has to `use
! mod_solver_iface` to extend flow_solver_t; if flow_solver_t's own
! module also `use`d the concrete modules, that would be a circular
! module dependency).

  use mod_solver_iface, only: flow_solver_t
  implicit none(type, external)
  private

  public :: make_flow_solver

contains

  ! eulfs/su2 mirror the exact CLI-derived logicals main.f90 already
  ! computes (neo = .not.eulfs .and. .not.su2, unchanged) -- this
  ! factory doesn't touch CLI parsing or the "Use eulfs?"/"Use su2?"
  ! reporting, both of which keep reading the original logicals
  ! directly.
  function make_flow_solver(eulfs, su2) result(solver)
    use mod_eulfs_solver, only: eulfs_t
    use mod_su2_solver, only: su2_t
    use mod_neo_solver, only: neo_t
    logical, intent(in) :: eulfs, su2
    class(flow_solver_t), allocatable :: solver

    if (eulfs) then
      allocate (eulfs_t :: solver)
    elseif (su2) then
      allocate (su2_t :: solver)
    else
      allocate (neo_t :: solver)
    end if
  end function make_flow_solver

end module mod_solver_registry
