program NEO2triangle

  ! ==============================
  ! Create .node and .poly files 
  ! .poly file as the same
  ! ==============================

  ! ----------------------------------------------------
  ! To compile :  gfortran -o NEO2triangle main_d2t.f90
  ! ----------------------------------------------------
  !aldo
  !
  ! ribattezzato l'eseguibile per evitare confusione con l'altro convertitore
  !
  !
  !aldo
  ! ----------------------------------------------------
  ! To compile :  gfortran -o dat2triangle main_d2t.f90
  ! ----------------------------------------------------

  implicit none

  real*8,dimension(:,:),allocatable:: xy_node, Z                    ! coordinates of node / Roe's variable (Z1,Z2,Z3,Z4) 
  integer,dimension(:,:),allocatable:: bnd_face                     ! Boundary face
  real*8,dimension(:),allocatable::rho, u, v, H                     ! Density / Velocity (u,v) / Total enthalpy
  integer,dimension(:),allocatable::bmk ,bndf_mk                    ! Boundary marker of each node, boundary marker of faces

  integer:: NN, NBF, dim, nb_edge                                   ! # of node / # of boundary face / dimension / # of edges                a
  integer::i, temp, kspace, mark, s1, s2,AllocateStatus
  real*8::tempr
  character(len=50)::filename
  character(len=80)::EXECMD

  read(5,'(A)') filename
  kspace = len(trim(filename))

  ! ========================
  ! Reading of general data
  ! ========================

  write(6,*)"NEO2triangle: reading dims from './NEO_data/input/neogrid.grd'"
  open(1, file = './NEO_data/input/neogrid.grd', status = 'old')
  read(1,*) dim, temp, NN, NBF
  close(1)

  ! ===========
  ! Allocation
  ! ===========

  allocate(xy_node(2,NN),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  allocate(Z(4,NN),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  allocate(bmk(NN),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  allocate(rho(NN), u(NN), v(NN), H(NN),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  allocate(bnd_face(2,NBF),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  allocate(bndf_mk(NBF),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  xy_node(:,:) = 0
  Z(:,:) = 0
  bmk(:) = 0
  rho(:) = 0
  u(:) = 0
  v(:) = 0
  H(:) = 0
  bnd_face(:,:) = 0
  bndf_mk(:) = 0

  ! =================
  ! Thorough reading
  ! =================

  ! -------------------
  ! Reading of outfile
  ! -------------------

  write(6,*)"NEO2triangle: opening and reading from './NEO_data/output/vvvv.dat'"
  open(2, file = './NEO_data/output/vvvv.dat', status = 'old')
  read(2,*)
  read(2,*)
  read(2,*)
  read(2,*)
  do i = 1, NN
     read(2,*) xy_node(1,i), xy_node(2,i), rho(i), u(i), v(i), tempr, H(i)
  end do
  close(2)

  ! -----------------------------------------------
  ! Determination of boundary marker for each node
  ! -----------------------------------------------

  write(6,*)"NEO2triangle: reading file ", filename(1:kspace)//".edge"
  open(46, file = filename(1:kspace)//".edge", status = "old")
  read(46,*) nb_edge
  do i=1, nb_edge
     read(46,*) temp, s1, s2, mark
     if ((mark/=0)) then      
        bmk(s1) = 2
        bmk(s2) = 2
     end if
  end do
  close(46)


  ! --------------------------
  ! Connectivity are the same
  ! --------------------------
  
  !open(7, file = filename(1:kspace)//".poly", status = 'unknown')
  !read(7,*)
  !read(7,*)
  !do i = 1, NBF
  !   read(7,*) temp, bnd_face(1,i), bnd_face(2,i), bndf_mk(i)
  !end do
  !close(7)

  ! ==========================
  ! Convert in Roe's variable
  ! ==========================

  do i=1,NN
     Z(1,i) = sqrt(rho(i))
     Z(2,i) = sqrt(rho(i))*H(i)
     Z(3,i) = sqrt(rho(i))*u(i)
     Z(4,i) = sqrt(rho(i))*v(i)
  end do

  ! ========================
  ! Writing of the new file
  ! ========================

  EXECMD = "mv -v "//filename(1:kspace)//".node "//filename(1:kspace)//".node.BAK"
  WRITE(6,*) EXECMD
  CALL SYSTEM(EXECMD) 

  ! -----------
  ! .node file
  ! -----------

  write(6,*)"NEO2triangle: writing file ", filename(1:kspace)//".node"
  open(3, file = filename(1:kspace)//".node", status = 'unknown')    
  write(3,*) NN, dim, 4, 1
  do i=1,NN
!   write(3,'(I6,1X,F20.18,2X,F20.18,2X,F20.18,2X,F20.18,2X,F20.18,2X,F20.18,1X,I4)') i, &
!         xy_node(1,i), xy_node(2,i), Z(1,i), Z(2,i), Z(3,i), Z(4,i), bmk(i) 
    write(3,*) i, xy_node(1,i), xy_node(2,i), Z(1,i), Z(2,i), Z(3,i), Z(4,i), bmk(i)
  end do
  close(3)  
!'(I6,6F20.15,I1)'
  ! -----------
  ! .poly file
  ! -----------

  EXECMD = "cp "//filename(1:kspace)//".poly "//filename(1:kspace)//".poly.BAK"
  WRITE(6,*) EXECMD
  CALL SYSTEM(EXECMD) 

  !open(4, file = filename(1:kspace)//".poly", status = 'unknown')
  !write(4,*) 0, dim, 0, 0
  !write(4,*) NBF, 1
  !do i = 1, NBF
   !  write(4,*) i, bnd_face(1,i), bnd_face(2,i), bndf_mk(i)
  !end do
  !write(4,*) 0
  !close(4)

  deallocate(xy_node, Z, bmk, rho, u, v, H, bnd_face, bndf_mk)

end program NEO2triangle
