add_executable(triangle2dat
  main.f
  seterr.f
  periodic.f
  xa112d.f
  rtri.f
)
undifi_link_vendor_libs(triangle2dat port fxdr mylib tirpc)
