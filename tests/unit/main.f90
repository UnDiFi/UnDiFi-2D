program undifi2d_unit_tests
! test-drive's own driver skeleton (see its README), extended with one
! new_testsuite(...) entry per collect_* module added under tests/unit/.
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type,&
  &select_suite, run_selected, get_argument
  use test_newton_solve, only: collect_newton_solve
  use test_co_shock, only: collect_co_shock
  use test_co_norm, only: collect_co_norm
  use test_co_urr, only: collect_co_urr
  use test_co_uqp, only: collect_co_uqp
  use test_co_utp, only: collect_co_utp
  use test_rd_dps, only: collect_rd_dps
  use test_rd_dps_eq, only: collect_rd_dps_eq
  use test_interp, only: collect_interp
  use test_two_shock_theory, only: collect_two_shock_theory
  use test_three_shock_theory, only: collect_three_shock_theory
  implicit none

  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)

  stat = 0
  testsuites = [&
  &new_testsuite("newton_solve", collect_newton_solve),&
  &new_testsuite("co_shock", collect_co_shock),&
  &new_testsuite("co_norm", collect_co_norm),&
  &new_testsuite("co_urr", collect_co_urr),&
  &new_testsuite("co_uqp", collect_co_uqp),&
  &new_testsuite("co_utp", collect_co_utp),&
  &new_testsuite("rd_dps", collect_rd_dps),&
  &new_testsuite("rd_dps_eq", collect_rd_dps_eq),&
  &new_testsuite("interp", collect_interp),&
  &new_testsuite("two_shock_theory", collect_two_shock_theory),&
  &new_testsuite("three_shock_theory", collect_three_shock_theory)&
  &]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
      else
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end if
    else
      write (error_unit, '(a)') "Available testsuites"
      do is = 1, size(testsuites)
        write (error_unit, '(a)') "- "//testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, a)') stat, " test(s) failed"
    error stop 1
  end if
end program undifi2d_unit_tests
