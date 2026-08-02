# -I./common8/ from the legacy Makefile is a stale flag; no such directory
# exists in this tree, so it is dropped here.
add_executable(dat2triangle
  main.f
  readat.f
  nodcod.f
  seterr.f
  solzne.f
  wtri.f
  movep.f
  renumber.f
  wkspce.f
)
undifi_link_vendor_libs(dat2triangle mylib port fxdr tirpc)
