subroutine neo2triangle(filename)

! Phase 5a (ROADMAP.md #17): in-process port of
! source_utils/NEO2triangle/main_d2t.f90 -- see na2vvvv.f90 for the
! rationale. filename is passed directly instead of read from stdin;
! hardcoded low unit numbers -> newunit. The two CALL SYSTEM backup
! renames (.node -> .node.BAK, .poly -> .poly.BAK) are left as-is: they
! are cheap single-file mv/cp operations, not full executable spawns,
! and out of this port's scope. The original standalone executable
! (source_utils/NEO2triangle/) is left untouched and buildable.

  ! ==============================
  ! Create .node and .poly files
  ! .poly file as the same
  ! ==============================

  implicit none

  character(len=*), intent(in) :: filename

  real*8,dimension(:,:),allocatable:: xy_node, Z                    ! coordinates of node / Roe's variable (Z1,Z2,Z3,Z4)
  integer,dimension(:,:),allocatable:: bnd_face                     ! Boundary face
  real*8,dimension(:),allocatable::rho, u, v, H                     ! Density / Velocity (u,v) / Total enthalpy
  integer,dimension(:),allocatable::bmk ,bndf_mk                    ! Boundary marker of each node, boundary marker of faces

  integer:: NN, NBF, dim, nb_edge                                   ! # of node / # of boundary face / dimension / # of edges                a
  integer::i, temp, kspace, mark, s1, s2,AllocateStatus
  real*8::tempr
  character(len=80)::EXECMD
  integer :: log_unit, u1, u2, u46, u3

  open(newunit=log_unit, file='log/NEO2triangle.log', status='unknown')

  kspace = len(trim(filename))

  ! ========================
  ! Reading of general data
  ! ========================

  write(log_unit,*)"NEO2triangle: reading dims from './NEO_data/input/neogrid.grd'"
  open(newunit=u1, file = './NEO_data/input/neogrid.grd', status = 'old')
  read(u1,*) dim, temp, NN, NBF
  close(u1)

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

  write(log_unit,*)"NEO2triangle: opening and reading from './NEO_data/output/vvvv.dat'"
  open(newunit=u2, file = './NEO_data/output/vvvv.dat', status = 'old')
  read(u2,*)
  read(u2,*)
  read(u2,*)
  read(u2,*)
  do i = 1, NN
     read(u2,*) xy_node(1,i), xy_node(2,i), rho(i), u(i), v(i), tempr, H(i)
  end do
  close(u2)

  ! -----------------------------------------------
  ! Determination of boundary marker for each node
  ! -----------------------------------------------

  write(log_unit,*)"NEO2triangle: reading file ", filename(1:kspace)//".edge"
  open(newunit=u46, file = filename(1:kspace)//".edge", status = "old")
  read(u46,*) nb_edge
  do i=1, nb_edge
     read(u46,*) temp, s1, s2, mark
     if ((mark/=0)) then
        bmk(s1) = 2
        bmk(s2) = 2
     end if
  end do
  close(u46)

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
  write(log_unit,*) EXECMD
  CALL SYSTEM(EXECMD)

  ! -----------
  ! .node file
  ! -----------

  write(log_unit,*)"NEO2triangle: writing file ", filename(1:kspace)//".node"
  open(newunit=u3, file = filename(1:kspace)//".node", status = 'unknown')
  write(u3,*) NN, dim, 4, 1
  do i=1,NN
    write(u3,*) i, xy_node(1,i), xy_node(2,i), Z(1,i), Z(2,i), Z(3,i), Z(4,i), bmk(i)
  end do
  close(u3)

  ! -----------
  ! .poly file
  ! -----------

  EXECMD = "cp "//filename(1:kspace)//".poly "//filename(1:kspace)//".poly.BAK"
  write(log_unit,*) EXECMD
  CALL SYSTEM(EXECMD)

  deallocate(xy_node, Z, bmk, rho, u, v, H, bnd_face, bndf_mk)

  close(log_unit)

end subroutine neo2triangle
