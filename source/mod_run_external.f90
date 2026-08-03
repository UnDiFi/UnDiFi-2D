module mod_run_external

  use mod_kinds, only: i4
  implicit none(type, external)
  private
  public :: run_external

contains

  ! Runs cmd via the C library system() call (still the underlying
  ! transport for now -- see ROADMAP.md Phase 5a.2 for the eventual
  ! library-coupled replacement). name is used only in the error
  ! message. By default a nonzero exit status is fatal; pass
  ! fatal_on_error=.false. for best-effort commands (periodic backups,
  ! cleanup) where the caller intends to proceed regardless.
  function run_external(cmd, name, fatal_on_error) result(status)
    character(len=*), intent(in) :: cmd
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: fatal_on_error
    integer(i4) :: status

    integer(i4) :: system

    logical :: check

    check = .true.
    if (present(fatal_on_error)) check = fatal_on_error

    status = system(cmd)
    call flush (6)

    if (check .and. status /= 0) then
      write (6, *) trim(name)//' has returned an error code ifail = ', status
      write (6, *) 'command: '//trim(cmd)
      error stop status
    end if
  end function run_external

end module mod_run_external
