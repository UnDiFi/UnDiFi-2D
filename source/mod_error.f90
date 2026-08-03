module mod_error

  implicit none(type, external)
  private
  public :: fatal

contains

  subroutine fatal(msg, code)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: code

    write (0, '(a)') trim(msg)
    error stop code
  end subroutine fatal

end module mod_error
