module reading_mod

  implicit none
  
contains
  
  subroutine reading(xy_node,bnd_f,bmk_n,bmk_f,fname,k)
  
    real*8,dimension(:,:),intent(out)::xy_node
    integer,dimension(:,:),intent(out)::bnd_f

    integer,dimension(:),intent(out)::bmk_f, bmk_n
    character(len=*),intent(in)::fname
    integer,intent(in)::k

    integer::nb_node, nb_bnd_f, nb_ele, temp
    integer::i, a, b
    real*8::x1,x2,y1,y2

    xy_node(:,:) = 0
    bmk_n(:) = 0
    bnd_f(:,:) = 0
    bmk_f(:) = 0

    ! ========
    ! Reading
    ! ========
    open(1, file=fname(1:k)//".grd", status='old')
    read(1,*) temp, nb_ele, nb_node, nb_bnd_f
    read(1,*) 
    do i = 1,nb_ele
      read(1,*)
    end do
    read(1,*)
    do i = 1, nb_node
      read(1,*) xy_node(1,i), xy_node(2,i)
    end do
    read(1,*)
    do i = 1, nb_bnd_f
      read(1,*) bnd_f(1,i), bnd_f(2,i)
    end do
    close(1)

    open(2,file="bords.data",status="old")
    read(2,*)
    read(2,*) y1, x2, y2, x1
    close(2)
    
    ! ================
    ! Boundary marker
    ! ================
    do i = 1, nb_node
       if ( (xy_node(1,i)==x1) .or. (xy_node(1,i)==x2) .or. (xy_node(2,i)==y1) .or. (xy_node(2,i)==y2)) then
          bmk_n(i) = 2
       end if
    end do
    
    do i = 1, nb_bnd_f
       a = bnd_f(1,i)+1
       b = bnd_f(2,i)+1
       
       if ( (xy_node(1,a)==x1) .and. (xy_node(1,b)==x1) ) bmk_f(i) = 4
       if ( (xy_node(1,a)==x2) .and. (xy_node(1,b)==x2) ) bmk_f(i) = 2
       if ( (xy_node(2,a)==y1) .and. (xy_node(2,b)==y1) ) bmk_f(i) = 1
       if ( (xy_node(2,a)==y2) .and. (xy_node(2,b)==y2) ) bmk_f(i) = 3
       
!      if ( (xy_node(1,a)==0) .and. (xy_node(1,b)==0) ) bmk_f(i) = 4
!      if ( (xy_node(1,a)==5) .and. (xy_node(1,b)==5) ) bmk_f(i) = 2
!      if ( (xy_node(2,a)==0) .and. (xy_node(2,b)==0) ) bmk_f(i) = 1
!      if ( (xy_node(2,a)==1) .and. (xy_node(2,b)==1) ) bmk_f(i) = 3
   
    end do

  end subroutine reading

end module reading_mod
