# Mesh/solver format converters and pre/post-processing tools used by the
# UnDiFi-2D outer loop (see the system() calls in source/main.f) and by
# users inspecting results.
#
# fsplot is excluded by default: it links against a proprietary Tecplot
# 360 SDK plus X11/Motif/OpenGL (see fsplot/plot/Makefile), none of which
# are available in a normal dev or CI environment. Enable it explicitly
# with -DUNDIFI_BUILD_FSPLOT=ON on a machine that has the SDK installed
# and TEC_360_2008 set.
option(UNDIFI_BUILD_FSPLOT "Build fsplot (requires proprietary Tecplot SDK + X11/Motif)" OFF)

add_subdirectory(triangle-v1.6)
add_subdirectory(triangle2dat)
add_subdirectory(dat2triangle)
add_subdirectory(su22triangle)
add_subdirectory(triangle2su2)
add_subdirectory(NEO2triangle)
add_subdirectory(Na00x2vvvv)
add_subdirectory(Triangle2grd)
add_subdirectory(Grid_0)
add_subdirectory(Na_creation)
add_subdirectory(dat2paraview)

if(UNDIFI_BUILD_FSPLOT)
  add_subdirectory(fsplot/plot)
endif()
