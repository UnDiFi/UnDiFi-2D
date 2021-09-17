module boundary_edge_mod

! ========================================
! Determinate the number of boundary edge
! ========================================

  implicit none
  
  contains

  subroutine boundary_edge(nb_ed,nb_bed,f_edge)
  integer,intent(out)::nb_bed                         ! Number of boundary edge
  integer,intent(in)::nb_ed                           ! Number of edge
  character(len=*),intent(in)::f_edge                 ! Edge file
  integer::i, temp, mark                              ! mark = boundary marker
  
  nb_bed = 0

  open(46, file=f_edge, status='old')
  read(46,*)
  do i=1,nb_ed
     read(46,*) temp, temp, temp, mark
     if ((mark/=0)) then   !.and. (mark/=10) 
        nb_bed = nb_bed+1
     end if
  end do
  close(46)

  end subroutine boundary_edge

end module boundary_edge_mod
