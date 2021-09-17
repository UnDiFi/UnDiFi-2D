program trianglegrd

  ! =============================================================
  ! Conversion of .triangle format to .grd format (Mario's code)
  ! This program needs .edge, .neigh and .poly files
  ! Allow to create .grd file
  ! Inspiré de Triangle2dat.f
  ! =============================================================

  ! ************************************************************************************
  ! To compile : gfortran -o triangle2grd boundary_edge_mod.f90 r_files_mod.f90 main.f90
  ! To run     : ./exe
  ! ************************************************************************************

  use boundary_edge_mod
  use r_files_mod

  implicit none

  integer:: nb_bed, nb_ele, nb_node, nb_att, nb_ed, K, temp, i     ! # bound_edge / # of element / # node / # attribut / # edge (whole grid) / temporary(x3)
  real*8,dimension(:,:),allocatable:: xy_node                      ! Coodonnate of nodes
  real*8,dimension(:,:),allocatable:: att_node                     ! Attribute of nodes
  integer,dimension(:,:),allocatable:: b_edge, bmk_f               ! Boundary edges / Boundary makers of faces
  integer,dimension(:,:),allocatable:: ele, neigh                  ! Triangle's nodes / neigh's node
  integer,dimension(:),allocatable::bmk_edge
  integer,dimension(10)::type
  integer,parameter:: dim=2                                        ! Dimension of the space
  integer,dimension(10):: lenfname                                 ! Array of number of character
  character(len=100),dimension(10):: fname                         ! File


  !print*,'========================================='
  !print*,'Conversion of triangle files in grid file'
  !print*,'========================================='

  ! =============
  ! General data
  ! =============
  fname(10)= "file.1"                             ! ***** Lecture du fichier automatique  *******
  print*, 'Enter fname'
  read*, fname(10)

  K = len(trim(fname(10)))                        ! Number of character of the entering file / "trim" allows to remove spaces
  FNAME(1)(1:K+5) = FNAME(10)(1:K)//".node"
  FNAME(2)(1:K+4) = FNAME(10)(1:K)//".ele"
  FNAME(3)(1:K+6) = FNAME(10)(1:K)//".neigh"      ! Not necessary, allows to verify boundary faces
  FNAME(4)(1:K+5) = FNAME(10)(1:K)//".edge"
  lenfname(1) = K+5
  lenfname(2) = K+4
  lenfname(3) = K+6
  lenfname(4) = K+5

  ! -----------------
  ! Reading of .node
  ! -----------------
  open(1,file=fname(1)(1:K+5),status='old')
  read(1,*) nb_node, temp, nb_att, temp
  if(nb_att.NE.4)then
     print*,'Number of attributes should be 4 in ',fname(1)(1:K+5),'while it is ',nb_att
     STOP
  endif
  if(temp.NE.1)then
     print*,'Number of bndry markers should be 1 in ',fname(1)(1:K+5),'while it is ',temp
     STOP
  endif
  close(1)

  ! ----------------
  ! Reading of .ele
  ! ----------------
  open(2,file=FNAME(2)(1:K+4),status="old")
  rewind(2)
  read(2,*) nb_ele
  close(2)

  ! -----------------
  ! Reading of .edge
  ! -----------------
  open(3,file=FNAME(4)(1:K+5),status="old")
  rewind(3)
  read(3,*) nb_ed
  close(3)

  ! -------------------------
  ! Reading of face_type.dat
  ! -------------------------
  open(9,file="type.dat",status="old")
  read(9,*)
  do i=1,10
     read(9,*) temp, type(i)
  end do
  close(9)

  ! =============================
  ! Determinate de boundary edge
  ! =============================

  call boundary_edge(nb_ed,nb_bed,fname(4)(1:K+5))

  ! ===========
  ! Allocation
  ! ===========
  allocate(xy_node(2,nb_node))
  allocate(att_node(nb_att,nb_node))
  allocate(b_edge(2,nb_bed))
  allocate(ele(3,nb_ele))
  allocate(neigh(3,nb_ele))
  allocate(bmk_f(3,nb_bed))
  allocate(bmk_edge(nb_bed))

  xy_node(:,:) = 0
  att_node(:,:) = 0
  b_edge(:,:) = 0
  ele(:,:) = 0
  neigh(:,:) = 0
  bmk_f(:,:) = 0
  bmk_edge(:) = 0

  ! ==========================
  ! Detailed reading of files
  ! ==========================
  call r_files(xy_node,att_node,ele,neigh,b_edge,bmk_f,nb_bed,fname,lenfname,type,bmk_edge)

  ! =======================================
  ! Writing the mesh and connectivity file
  ! =======================================
  !fname(9) = './../source_utils/NEO_source/input/neogrid.grd'
  fname(9) = './NEO_data/input/neogrid.grd'
  K = len(trim(fname(9)))                        ! Number of character of the entering file / "trim" allows to remove spaces
  write(6,*)'Triangle2grd is writing file: ',fname(9)(1:K)
  open(4,file=fname(9)(1:K),status="unknown")
  rewind(4)

  write(4,*) dim, nb_ele, nb_node, nb_bed
  write(4,*) ' '

  do i=1, nb_ele
     write(4,*) ele(1,i)-1, ele(2,i)-1, ele(3,i)-1
  end do
  write(4,*) ' '

  do i=1, nb_node
     write(4,*) xy_node(1,i), xy_node(2,i), i-1
  end do
  write(4,*) ' '

  do i=1, nb_bed
     write(4,*) b_edge(1,i)-1, b_edge(2,i)-1, bmk_f(1,i), bmk_f(2,i), bmk_f(3,i)
  end do
  close(4)

  deallocate(xy_node,att_node,b_edge,ele,neigh,bmk_f)

end program trianglegrd
