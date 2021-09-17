program main

! write file vvvv_input.dat
!
! ===============================================
!  To compile : gfortran -o na2vvvv main_n2v.f90
! ===============================================

  implicit none

  real*8,dimension(:),allocatable::z1,z2,z3,z4,x,y
  real*8,dimension(:),allocatable::rho,p,u,v,H,k
  integer::lenK,temp,nb_node,i

  character(len=128)::fname

  print*," Enter the name of the file "
  read*,fname

  lenK = len(trim(fname))
  write(6,*)'Opening file ',fname(1:lenK)

  open(1,file=fname(1:lenK)//".node",status="old")
  read(1,*) nb_node
  close(1)

  allocate(x(nb_node),y(nb_node),z1(nb_node),z2(nb_node),z3(nb_node),z4(nb_node))
  allocate(rho(nb_node),p(nb_node),u(nb_node),v(nb_node),H(nb_node),k(nb_node))

  ! Reading
  open(2,file=fname(1:lenK)//".node",status="old")
  read(2,*)
  do i=1,nb_node
     read(2,*)temp,x(i),y(i),z1(i),z2(i),z3(i),z4(i)
  end do
  close(2)

  ! Computation
  do i=1,nb_node
     rho(i) = z1(i)*z1(i)
     H(i)   = z2(i)/z1(i)
     u(i)   = z3(i)/z1(i)
     v(i)   = z4(i)/z1(i)
     k(i)   = 0.50d0*(z4(i)*z4(i) + z3(i)*z3(i))/rho(i)
     p(i) = 0.40d0/1.40d0 * rho(i) * ( H(i) - k(i) )
  end do

  ! Writing
  !fname = "../../source_utils/NEO_source/output/vvvv_input.dat"
  fname = "./NEO_data/output/vvvv_input.dat"
  lenK = len(trim(fname))
  write(6,*)'Writing file ',fname(1:lenK)
  open(3,file=fname(1:lenK),status="unknown")
  write(3,"(1A37)") "TITLE      =  Unstructured grid data "
  write(3,"(1A55)") "VARIABLES  =  x  y  rho  u  v p H  Ma  s T time  sensor"
  write(3,"(1A59)") "ZONE    N  =  0    E  =  0    F = FEPOINT    ET = TRIANGLE "
  write(3,*)
  do i=1,nb_node
!    write(3,"(1F12.8,11F13.8)") x(i), y(i), rho(i), u(i), v(i), p(i), H(i), 0., 0., 0., 0., 0.
     write(3,"(12F20.16)") x(i), y(i), rho(i), u(i), v(i), p(i), H(i), 0., 0., 0., 0., 0.
!     write(3,*) x(i), y(i), rho(i), u(i), v(i), p(i), H(i), 0., 0., 0., 0., 0.
  end do
  write(3,*)
  close(3)

  deallocate(x,y,z1,z2,z3,z4,rho,p,u,v,H,k)

end program main
