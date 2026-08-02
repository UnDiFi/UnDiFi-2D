add_executable(su22triangle
  main.f
  seterr.f
  rtri.f
  nodcod.f
  wtri.f
  rsu2.f
)
# -std=legacy -fno-range-check applied globally at top level.
undifi_link_vendor_libs(su22triangle port mylib)
