module mod_kinds

  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none(type, external)
  private
  public :: wp, i4

  integer, parameter :: wp = real64
  integer, parameter :: i4 = int32

end module mod_kinds
