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
  test_co_norm.f90
  test_co_urr.f90
  test_co_uqp.f90
  test_co_utp.f90
  test_rd_dps.f90
  test_rd_dps_eq.f90
  test_interp.f90
  ${CMAKE_SOURCE_DIR}/source/mod_kinds.f90
  ${CMAKE_SOURCE_DIR}/source/mod_constants.f90
  ${CMAKE_SOURCE_DIR}/source/mod_search.f90
  ${CMAKE_SOURCE_DIR}/source/mod_newton_solve.f90
  ${CMAKE_SOURCE_DIR}/source/fnd_phps.f90
  ${CMAKE_SOURCE_DIR}/source/co_shock.f90
  # co_norm.f90's own transitive closure -- mod_special_point.f90 declares
  # `external` calls to co_utp/co_uqp/co_urr/co_intr_pnt(co_pnt_dspl.f90)/
  # findbedg(fx_msh_sps.f90)/finder(interp.f90) inside deferred type-bound
  # procedures that the co_norm test itself never exercises (nspecpoints=0
  # in every test case here), but the linker still needs every symbol
  # resolved since Fortran doesn't dead-code-eliminate unused type-bound
  # procedures. Pulled in whole rather than stubbed, per this suite's own
  # "real production kernel, no mocking" convention.
  ${CMAKE_SOURCE_DIR}/source/mod_log.f90
  ${CMAKE_SOURCE_DIR}/source/mod_special_point.f90
  ${CMAKE_SOURCE_DIR}/source/mod_special_point_registry.f90
  ${CMAKE_SOURCE_DIR}/source/co_utp.f90
  ${CMAKE_SOURCE_DIR}/source/co_uqp.f90
  ${CMAKE_SOURCE_DIR}/source/co_urr.f90
  ${CMAKE_SOURCE_DIR}/source/co_pnt_dspl.f90
  ${CMAKE_SOURCE_DIR}/source/fx_msh_sps.f90
  ${CMAKE_SOURCE_DIR}/source/interp.f90
  ${CMAKE_SOURCE_DIR}/source/co_norm.f90
  ${CMAKE_SOURCE_DIR}/source/rd_dps.f90
  ${CMAKE_SOURCE_DIR}/source/rd_dps_eq.f90
)
target_include_directories(undifi2d_unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_libraries(undifi2d_unit_tests PRIVATE test-drive-lib)
# interp.f90's finder/finder_search call area/icycl, which live in the
# vendored libmylib archive (same one source/CMakeLists.txt links into the
# real UnDiFi-2D executable) rather than anywhere under source/.
undifi_link_vendor_libs(undifi2d_unit_tests mylib)

# co_norm.f90 unconditionally opens log/co_norm.log (relative to ctest's
# working directory, which defaults to this binary dir) -- create it once
# at configure time rather than depending on a real driver run having made
# it first.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/log)

# One ctest entry per test-drive testsuite, selected by name via
# test-drive's own get_argument(1, suite_name) convention (see main.f90)
# -- a failing suite is then individually reportable/re-runnable, not
# just "the test binary failed".
set(undifi2d_unit_test_suites
  newton_solve
  co_shock
  co_norm
  co_urr
  co_uqp
  co_utp
  rd_dps
  rd_dps_eq
  interp
)
foreach(_suite ${undifi2d_unit_test_suites})
  add_test(NAME "unit.${_suite}" COMMAND undifi2d_unit_tests "${_suite}"
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
