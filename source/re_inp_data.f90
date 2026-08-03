! Read the file input.dat containing the information concerning
! mesh generation in the proximity of the shocks shock integration
! additional hole points

subroutine re_inp_data

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_config, only: config_t, read_config
  implicit none(type, external)
  include 'paramt.h'

!     .. local scalars ..
  type(config_t) :: cfg
  integer(i4) i

!     open log file
  open (8, file='log/re_inp_data.log')

! open file input.dat containing the information concerning
!  1) mesh generation in the proximity of the shocks
!  2) shock integration
!  3) additional hole points
!  4) periodic boundarìes

  open (12, file='input.dat', status='old')
  write (8, *) ' open file input.dat'
  cfg = read_config(12)

  eps = cfg%eps
  sndmin = cfg%sndmin
  dxcell = cfg%dxcell
  shrelax = cfg%shrelax
  ibak = cfg%ibak
  ga = cfg%ga
  gm1 = cfg%gm1
  naddholes = cfg%naddholes
  do i = 1, naddholes
    caddhole(1, i) = cfg%caddhole(1, i)
    caddhole(2, i) = cfg%caddhole(2, i)
  end do
  nprdbnd = cfg%nprdbnd
  do i = 1, nprdbnd
    prdbndclr(1, i) = cfg%prdbndclr(1, i)
    prdbndclr(2, i) = cfg%prdbndclr(2, i)
    prdbndclr(3, i) = cfg%prdbndclr(3, i)
  end do
  flt_dspeed = cfg%flt_dspeed
  imtf = cfg%imtf

  close (12)
  close (8)
  return
end subroutine re_inp_data
