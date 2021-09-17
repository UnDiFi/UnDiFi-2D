module writing_mod
  
  implicit none
  
contains

  subroutine writing(xy_node,bnd_f,bmk_n,bmk_f,nb_node,nb_bnd_f,Xs,Z1D,Z2D,Z3D,Z4D,Z1U,Z2U,Z3U,Z4U)

    real*8,dimension(:,:),intent(in)::xy_node
    integer,dimension(:,:),intent(in)::bnd_f
    integer,dimension(:),intent(in)::bmk_f, bmk_n

    real*8,intent(in)::Xs

    real*8,intent(in)::Z1D, Z2D, Z3D, Z4D         
    real*8,intent(in)::Z1U, Z2U, Z3U, Z4U 

    integer::nb_node, nb_bnd_f, nb_ele
    integer::i

    ! ===================
    ! Writing .node file
    ! ===================
    open(1, file="na00.node", status="unknown")
    write(1,*) nb_node, 2, 4, 1
    do i = 1, nb_node
       if (xy_node(1,i) <= Xs) then
          write(1,*) i,xy_node(1,i),xy_node(2,i),Z1D,Z2D,Z3D,Z4D,bmk_n(i) 
       end if
       if (xy_node(1,i) > Xs) then
          write(1,*) i,xy_node(1,i),xy_node(2,i),Z1U,Z2U,Z3U,Z4U,bmk_n(i)
       end if
    end do
    close(1)

    ! ===================
    ! Writing .poly file
    ! ===================
    open(2, file="na00.poly", status="unknown")
    write(2,*) 0, 2, 0, 0
    write(2,*) nb_bnd_f, 1
    do i = 1 ,nb_bnd_f
       write(2,*) i,bnd_f(1,i)+1, bnd_f(2,i)+1, bmk_f(i)
    end do
    write(2,*) 0
    close(2)

  end subroutine writing

end module writing_mod
