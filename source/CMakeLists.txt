# UnDiFi-2D driver. Source list mirrors the OBJS list in the legacy
# Makefile (the authoritative list -- SRCS there has drifted, see the
# co_pnt_dipl.f typo fixed alongside this CMake build).
set(UNDIFI2D_SOURCES
  mod_kinds.f90
  mod_error.f90
  fx_msh_sps.f90
  co_utp.f90
  co_uqp.f90
  co_urr.f90
  chg_bnd_ptr.f90
  main.f90
  interp.f90
  rand.f90
  readmesh.f90
  readpmap.f90
  rtri.f90
  fnd_phps.f90
  rdstrshp.f90
  fixendpts.f90
  setbndrynodeptr.f90
  re_sdw_info.f90
  re_inp_data.f90
  wrt_sdw_info.f90
  co_norm.f90
  co_pnt_dspl.f90
  co_state_dps.f90
  fx_state_dps.f90
  mv_dps.f90
  fx_dps_loc.f90
  ch_sh_topology.f90
  fltr_dls.f90
  rd_dps.f90
  rd_dps_eq.f90
  co_shock.f90
  solzne.f90
  pr_sh_state.f90
  fx_sh_state.f90
  wtri.f90
  wtri0.f90
  wsh_mean.f90
  calc_vel.f90
  mv_grid.f90
)

add_executable(UnDiFi-2D ${UNDIFI2D_SOURCES})

# paramt.h / shock.com are resolved via the Fortran INCLUDE statement
# (not the C preprocessor), which gfortran/ifx/nvfortran all search via -I.
target_include_directories(UnDiFi-2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

undifi_link_vendor_libs(UnDiFi-2D port mylib fxdr sparse-blas tirpc)
