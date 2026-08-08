# Phase 7.5 (ROADMAP.md #20): unit tests for the physics kernels, using
# test-drive (https://github.com/fortran-lang/test-drive). Fetched rather
# than vendored -- it's a small, actively-maintained single-purpose
# testing library with no transitive dependencies of its own, unlike the
# heavier third-party codes (EulFS, NEO, triangle) this project vendors
# as prebuilt archives in lib/.
include(FetchContent)
FetchContent_Declare(
  test-drive
  GIT_REPOSITORY https://github.com/fortran-lang/test-drive
  GIT_TAG v0.6.1
)
# Build only test-drive's own library, not its own test suite.
set(TESTDRIVE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(test-drive)

# Each unit-test source file compiles the real production kernel it
# exercises straight from source/ -- no reimplementation, no mocking.
# Grows one entry per increment as more kernels get coverage (co_shock,
# co_norm, special_point_t%residual, ...); paramt.h's COMMON/PARAMT/ (GA
# etc.) is shared by name with whatever the test sets it to, same as any
# other translation unit in the real driver.
add_executable(undifi2d_unit_tests
  main.f90
  test_newton_solve.f90
  test_co_shock.f90
  ${CMAKE_SOURCE_DIR}/source/mod_kinds.f90
  ${CMAKE_SOURCE_DIR}/source/mod_constants.f90
  ${CMAKE_SOURCE_DIR}/source/mod_search.f90
  ${CMAKE_SOURCE_DIR}/source/mod_newton_solve.f90
  ${CMAKE_SOURCE_DIR}/source/fnd_phps.f90
  ${CMAKE_SOURCE_DIR}/source/co_shock.f90
)
target_include_directories(undifi2d_unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(undifi2d_unit_tests PRIVATE test-drive-lib)

# One ctest entry per test-drive testsuite, selected by name via
# test-drive's own get_argument(1, suite_name) convention (see main.f90)
# -- a failing suite is then individually reportable/re-runnable, not
# just "the test binary failed".
set(undifi2d_unit_test_suites
  newton_solve
  co_shock
)
foreach(_suite ${undifi2d_unit_test_suites})
  add_test(NAME "unit.${_suite}" COMMAND undifi2d_unit_tests "${_suite}")
endforeach()
