program create_na00

  use reading_mod
  use writing_mod

! =====================================
!> Create na00.node and na00.poly files
!> for a steady planar shock
! =====================================

! gfortran -o exe_na physical_values.f90 reading_mod.f90 writing_mod.f90 main.f90

  implicit none

  real*8:: Xs
  real*8::Z1D, Z2D, Z3D, Z4D           ! Downstream Roe's variables
  real*8::Z1U, Z2U, Z3U, Z4U           ! Upstream Roe's variables

  real*8,dimension(:,:),allocatable::xy_node
  integer,dimension(:,:),allocatable::bnd_f
  integer,dimension(:),allocatable::bmk_f, bmk_n
  character(len=50)::fname

  integer::nb_node, nb_bnd_f, nb_ele, temp, k

  print*,"Insert the name of the file .grid to convert"
  read*,fname

  print*,'Select the shock position : Xs'
  read*,Xs

  k=len(trim(fname))

  ! ================
  ! Physical values
  ! ================
  call physical_values(Z1D,Z2D,Z3D,Z4D,Z1U,Z2U,Z3U,Z4U)

  ! =====================
  ! Reading general data
  ! =====================
  open(1, file=fname(1:k)//".grd", status='old')
  read(1,*) temp, nb_ele, nb_node, nb_bnd_f
  close(1)

  ! ===========
  ! Allocation
  ! ===========
  allocate(xy_node(2,nb_node))
  allocate(bnd_f(2,nb_bnd_f))
  allocate(bmk_f(nb_bnd_f))
  allocate(bmk_n(nb_node))

  ! ===========================================
  ! Reading & determination of boundary marker
  ! ===========================================
  call reading(xy_node,bnd_f,bmk_n,bmk_f,fname,k)

  ! =================================
  ! Writing of .node and .poly files
  ! =================================
  call writing(xy_node,bnd_f,bmk_n,bmk_f,nb_node,nb_bnd_f,Xs,Z1D,Z2D,Z3D,Z4D,Z1U,Z2U,Z3U,Z4U)

  deallocate(xy_node,bnd_f,bmk_f,bmk_n)

end program create_na00
