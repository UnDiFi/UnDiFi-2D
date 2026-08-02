add_executable(triangle2su2
  main.f
  seterr.f
  rtri.f
  wsu2.f
)
undifi_link_vendor_libs(triangle2su2 port mylib)
