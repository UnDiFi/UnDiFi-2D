! Fix and correct the boundary structure of the wedge

subroutine fx_bndry_wedge

  use mod_kinds, only: wp, i4
  implicit none(type, external)

!     open log file
  open (8, file='log/fx_bndry_wedge.log')

  close (8)
  return
end subroutine fx_bndry_wedge
