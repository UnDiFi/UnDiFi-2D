cmake_minimum_required(VERSION 3.20)

project(UnDiFi2D
  DESCRIPTION "UnDiFi-2D: an Unstructured Discontinuity Fitting code for 2D grids"
  LANGUAGES Fortran C)

# All Fortran here is legacy fixed-form F77 (source/, most of source_utils/)
# with a handful of free-form .f90 utilities. -std=legacy relaxes gfortran's
# strict-mode checks (Hollerith constants, old-style DO loops, etc.) that
# this code relies on; it is a superset of plain -O, never a restriction.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-std=legacy;-fno-range-check>")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g;-fcheck=all;-ffpe-trap=invalid,zero,overflow;-Wall>")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g;-check;all;-fpe0;-traceback>")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-g;-Mbounds;-Ktrap=fp>")
  endif()
else()
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-O2>")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Phase 4 (ROADMAP.md #16): OpenMP shared-memory parallelism, opt-in so
# every existing preset/workflow keeps building single-threaded until the
# CI/benchmark infra for the OMP variant lands. When ON, the !$omp
# directives already landed in co_norm.f90/co_state_dps.f90 (4.2) become
# real parallel regions instead of inert comments.
option(UNDIFI_ENABLE_OPENMP "Build with OpenMP shared-memory parallelism" OFF)
if(UNDIFI_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# Prebuilt vendored support libraries (libport, libmylib, libfxdr,
# libsparse-blas, libtirpc, libgeometry, libtriangulation, ...).
# These are checked in as static archives in lib/linux_intel_x86_64/ and
# rebuilt from lib/ via the existing lib/Makefile + lib/compile.sh; that
# vendoring workflow is out of scope for this CMake build (see ROADMAP.md
# Phase 0/issue #11 -- it targets source/ and source_utils/, i.e. the
# project's own code, not the third-party dependencies it vendors).
set(UNDIFI_VENDOR_LIB_DIR ${CMAKE_SOURCE_DIR}/lib/linux_intel_x86_64
    CACHE PATH "Directory containing the prebuilt vendored .a libraries")

# The vendored archives above were built with gfortran. A different Fortran
# compiler's driver (ifx, nvfortran) won't automatically pull in
# libgfortran's runtime, leaving I/O/runtime symbols (_gfortran_st_write,
# _gfortran_stop_string, ...) undefined at link time. Link libgfortran
# explicitly in that case -- a pragmatic bridge for now, not a substitute
# for eventually rebuilding lib/ itself per-compiler.
function(undifi_link_vendor_libs target)
  target_link_directories(${target} PRIVATE ${UNDIFI_VENDOR_LIB_DIR})
  target_link_libraries(${target} PRIVATE ${ARGN})
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    find_program(UNDIFI_GFORTRAN gfortran)
    if(UNDIFI_GFORTRAN)
      execute_process(
        COMMAND ${UNDIFI_GFORTRAN} -print-file-name=libgfortran.so
        OUTPUT_VARIABLE _libgfortran_path
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      get_filename_component(_libgfortran_dir "${_libgfortran_path}" DIRECTORY)
      if(_libgfortran_dir)
        target_link_directories(${target} PRIVATE ${_libgfortran_dir})
        target_link_libraries(${target} PRIVATE gfortran)
      endif()
    endif()
  endif()
endfunction()

enable_testing()

add_subdirectory(source)
add_subdirectory(source_utils)

# Phase 7.5 (ROADMAP.md #20): unit tests for the physics kernels
# (test-drive, fetched over the network at configure time -- opt-out for
# builds without network access, e.g. an offline vendored-lib rebuild).
option(UNDIFI_ENABLE_TESTS "Build the Fortran unit test suite (test-drive)" ON)
if(UNDIFI_ENABLE_TESTS)
  add_subdirectory(tests/unit)
endif()
