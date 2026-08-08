subroutine na2vvvv(fname)

! Phase 5a (ROADMAP.md #17): in-process port of
! source_utils/Na00x2vvvv/main_n2v.f90 -- was a standalone executable
! spawned via system() every outer iteration (mod_run_external.f90);
! the benchmark in issue #17 showed the process-spawn/file-round-trip
! overhead dominates wall time far more than this converter's own
! trivial per-node computation, so this in-process version drops the
! fork/exec, keeping the exact same numerics and the same
! log/na2vvvv.log output. fname is passed directly instead of read
! from stdin. The original standalone executable
! (source_utils/Na00x2vvvv/) is left untouched and buildable, per the
! issue's own "keep the executable path working throughout" guidance.

  implicit none

  character(len=*), intent(in) :: fname

  real*8, dimension(:), allocatable::z1, z2, z3, z4, x, y
  real*8, dimension(:), allocatable::rho, p, u, v, H, k
  integer::lenK, temp, nb_node, i, log_unit, u1, u2, u3

  character(len=128)::outfname

  open (newunit=log_unit, file='log/na2vvvv.log', status='unknown')

  lenK = len(trim(fname))
  write (log_unit, *) 'Opening file ', fname(1:lenK)

  open (newunit=u1, file=fname(1:lenK)//".node", status="old")
  read (u1, *) nb_node
  close (u1)

  allocate (x(nb_node), y(nb_node), z1(nb_node), z2(nb_node), z3(nb_node), z4(nb_node))
  allocate (rho(nb_node), p(nb_node), u(nb_node), v(nb_node), H(nb_node), k(nb_node))

! Reading
  open (newunit=u2, file=fname(1:lenK)//".node", status="old")
  read (u2, *)
  do i = 1, nb_node
    read (u2, *) temp, x(i), y(i), z1(i), z2(i), z3(i), z4(i)
  end do
  close (u2)

! Computation
  do i = 1, nb_node
    rho(i) = z1(i)*z1(i)
    H(i) = z2(i)/z1(i)
    u(i) = z3(i)/z1(i)
    v(i) = z4(i)/z1(i)
    k(i) = 0.50d0*(z4(i)*z4(i) + z3(i)*z3(i))/rho(i)
    p(i) = 0.40d0/1.40d0*rho(i)*(H(i) - k(i))
  end do

! Writing
  outfname = "./NEO_data/output/vvvv_input.dat"
  lenK = len(trim(outfname))
  write (log_unit, *) 'Writing file ', outfname(1:lenK)
  open (newunit=u3, file=outfname(1:lenK), status="unknown")
  write (u3, "(1A37)") "TITLE      =  Unstructured grid data "
  write (u3, "(1A55)") "VARIABLES  =  x  y  rho  u  v p H  Ma  s T time  sensor"
  write (u3, "(1A59)") "ZONE    N  =  0    E  =  0    F = FEPOINT    ET = TRIANGLE "
  write (u3, *)
  do i = 1, nb_node
    write (u3, "(12F20.16)") x(i), y(i), rho(i), u(i), v(i), p(i), H(i), 0., 0., 0., 0., 0.
  end do
  write (u3, *)
  close (u3)

  deallocate (x, y, z1, z2, z3, z4, rho, p, u, v, H, k)

  close (log_unit)

end subroutine na2vvvv
