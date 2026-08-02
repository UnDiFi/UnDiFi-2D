# Module dependency order (r_files_mod, boundary_edge_mod before main)
# is auto-detected by CMake's Fortran support from the .mod producer/
# consumer relationships, matching the order the legacy Makefile lists.
add_executable(triangle2grd
  r_files_mod.f90
  boundary_edge_mod.f90
  main.f90
)
