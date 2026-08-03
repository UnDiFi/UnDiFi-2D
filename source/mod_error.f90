module mod_error

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: fatal

contains

  subroutine fatal(msg, code)
    use mod_kinds, only: wp, i4
    character(len=*), intent(in) :: msg
    integer(i4), intent(in) :: code

    write (0, '(a)') trim(msg)
    error stop code
  end subroutine fatal

end module mod_error
