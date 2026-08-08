module boundary_edge_mod

! Phase 5a (ROADMAP.md #17): in-process copy of
! source_utils/Triangle2grd/boundary_edge_mod.f90, unchanged except the
! hardcoded unit 46 -> newunit (see triangle2grd.f90 for why). See
! na2vvvv.f90 for the rationale of this in-process port.

! ========================================
! Determinate the number of boundary edge
! ========================================

  implicit none

  contains

  subroutine boundary_edge(nb_ed,nb_bed,f_edge)
  integer,intent(out)::nb_bed                         ! Number of boundary edge
  integer,intent(in)::nb_ed                           ! Number of edge
  character(len=*),intent(in)::f_edge                 ! Edge file
  integer::i, temp, mark, u                           ! mark = boundary marker

  nb_bed = 0

  open(newunit=u, file=f_edge, status='old')
  read(u,*)
  do i=1,nb_ed
     read(u,*) temp, temp, temp, mark
     if ((mark/=0)) then   !.and. (mark/=10)
        nb_bed = nb_bed+1
     end if
  end do
  close(u)

  end subroutine boundary_edge

end module boundary_edge_mod
