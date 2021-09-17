!> It creates the neogrid0.grd file

program main
implicit none

real*8, dimension(:,:), allocatable :: xy_node
integer :: NN, NN_S, i, temp

open(1, file="na00.1.node", status="old")
read(1,*) NN
close(1)

open(1, file="na00001.node", status="old")
read(1,*) NN_S
close(1)

allocate(xy_node(2,NN_S))

open(3, file="na00.1.node", status="old")
read(3,*)
do i = 1, NN
  read(3,*) temp, xy_node(1,i), xy_node(2,i)
end do
close(3)

open(7, file="na00001.node", status="old")
read(7,*)
do i = 1, NN
  read(7,*)
end do
do i = NN+1, NN_S
  read(7,*) temp, xy_node(1,i), xy_node(2,i)
end do
close(7)

open(46, file="./NEO_data/input/neogrid0.grd", status="new")
do i = 1, NN_S
  write(46,'(2F20.16,I8)') xy_node(1,i), xy_node(2,i), i-1
end do
close (46)

end program main
