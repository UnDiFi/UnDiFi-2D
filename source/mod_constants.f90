module mod_constants
! Phase 2.1 (ROADMAP.md #13): replaces the PARAMETER block of paramt.h.

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private

  public :: zero, half, one, two, pi
  public :: ndim, ndof, nshmax, npshmax, neshmax, nspmax, naddholesmax, nprdbndmax

  real(wp), parameter :: zero = 0.00d0
  real(wp), parameter :: half = 0.5d0
  real(wp), parameter :: one = 1.00d0
  real(wp), parameter :: two = 2.00d0
  real(wp), parameter :: pi = 3.141593d0

  integer(i4), parameter :: ndim = 2
  integer(i4), parameter :: ndof = 4
  integer(i4), parameter :: nshmax = 10
  integer(i4), parameter :: npshmax = 500
  integer(i4), parameter :: neshmax = npshmax - 1
  integer(i4), parameter :: nspmax = 12
  integer(i4), parameter :: naddholesmax = 10
  integer(i4), parameter :: nprdbndmax = 1

end module mod_constants
