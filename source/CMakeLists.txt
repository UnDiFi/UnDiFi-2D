# UnDiFi-2D driver. Source list mirrors the OBJS list in the legacy
# Makefile (the authoritative list -- SRCS there has drifted, see the
# co_pnt_dipl.f typo fixed alongside this CMake build).
set(UNDIFI2D_SOURCES
  fx_msh_sps.f
  co_utp.f
  co_uqp.f
  co_urr.f
  chg_bnd_ptr.f
  main.f
  interp.f
  rand.f
  readmesh.f
  readpmap.f
  rtri.f
  fnd_phps.f
  rdstrshp.f
  fixendpts.f
  setbndrynodeptr.f
  re_sdw_info.f
  re_inp_data.f
  wrt_sdw_info.f
  co_norm.f
  co_pnt_dspl.f
  co_state_dps.f
  fx_state_dps.f
  mv_dps.f
  fx_dps_loc.f
  ch_sh_topology.f
  fltr_dls.f
  rd_dps.f
  rd_dps_eq.f
  co_shock.f
  solzne.f
  pr_sh_state.f
  fx_sh_state.f
  wtri.f
  wtri0.f
  wsh_mean.f
  calc_vel.f
  mv_grid.f
)

add_executable(UnDiFi-2D ${UNDIFI2D_SOURCES})

# paramt.h / shock.com are resolved via the Fortran INCLUDE statement
# (not the C preprocessor), which gfortran/ifx/nvfortran all search via -I.
target_include_directories(UnDiFi-2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

undifi_link_vendor_libs(UnDiFi-2D port mylib fxdr sparse-blas tirpc)
