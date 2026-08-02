# Shewchuk's Triangle mesh generator, built as a standalone executable
# (TRILIBRARY undefined). showme (the X11 mesh viewer) is not built here;
# it is a debugging aid, not part of the UnDiFi-2D outer-loop pipeline.
add_executable(triangle triangle.c)
target_compile_definitions(triangle PRIVATE LINUX NO_TIMER)
target_link_libraries(triangle PRIVATE m)
