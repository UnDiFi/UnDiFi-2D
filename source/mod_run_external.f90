module mod_run_external

  use mod_kinds, only: i4
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  implicit none(type, external)
  private
  public :: run_external

  interface
    ! The C library's system() (ISO C, so this is portable across
    ! compilers -- gfortran's SYSTEM() GNU intrinsic and ifx disagree on
    ! how to declare that extension, see ROADMAP.md Phase 5a.2 for the
    ! eventual library-coupled replacement of this call entirely).
    function c_system(command) bind(c, name="system") result(stat)
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: command(*)
      integer(c_int) :: stat
    end function c_system
  end interface

contains

  ! Runs cmd via the C library system() call. name is used only in the
  ! error message. By default a nonzero exit status is fatal; pass
  ! fatal_on_error=.false. for best-effort commands (periodic backups,
  ! cleanup) where the caller intends to proceed regardless.
  function run_external(cmd, name, fatal_on_error) result(status)
    character(len=*), intent(in) :: cmd
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: fatal_on_error
    integer(i4) :: status

    logical :: check

    check = .true.
    if (present(fatal_on_error)) check = fatal_on_error

    status = int(c_system(trim(cmd)//c_null_char), i4)
    call flush (6)

    if (check .and. status /= 0) then
      write (6, *) trim(name)//' has returned an error code ifail = ', status
      write (6, *) 'command: '//trim(cmd)
      error stop status
    end if
  end function run_external

end module mod_run_external
