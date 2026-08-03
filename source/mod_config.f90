module mod_config
! Phase 2.1 (ROADMAP.md #13): replaces the COMMON/PARAMT/ runtime block of
! paramt.h. read_config() reproduces re_inp_data's input.dat parse, but
! returns a config_t instead of writing into module/COMMON-global state.

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, naddholesmax, nprdbndmax, zero, one

  implicit none(type, external)
  private

  public :: config_t, read_config

  type :: config_t
    real(wp)    :: eps = zero
    real(wp)    :: sndmin = zero
    real(wp)    :: dxcell = zero
    real(wp)    :: shrelax = zero
    real(wp)    :: ga = zero
    real(wp)    :: gm1 = zero
    real(wp)    :: flt_dspeed = zero
    integer(i4) :: ibak = 0
    integer(i4) :: naddholes = 0
    integer(i4) :: nprdbnd = 0
    integer(i4) :: imtf = 0
    real(wp)    :: caddhole(ndim, naddholesmax) = zero
    integer(i4) :: prdbndclr(3, nprdbndmax) = 0
  end type config_t

contains

  function read_config(unit) result(cfg)
    integer(i4), intent(in) :: unit
    type(config_t) :: cfg

    integer(i4) :: i

    read (unit, *) cfg%eps        !< distance between the two shock faces
    read (unit, *) cfg%sndmin     !< maximum nondimensional distance of the phantom nodes
    read (unit, *) cfg%dxcell     !< lenght of shock edges
    read (unit, *) cfg%shrelax    !< relax coefficent of shock point integration
    read (unit, *) cfg%ibak       !< number of steps between the writing of a solution file and the next
    read (unit, *) cfg%ga         !< specific heat ratio
    cfg%gm1 = cfg%ga - one

    read (unit, *) cfg%naddholes  !< number of addition hole points
    do i = 1, cfg%naddholes
      read (unit, *) cfg%caddhole(1, i), cfg%caddhole(2, i) !< coordinates of addition hole points
    end do
    read (unit, *) cfg%nprdbnd    !< number of periodic boundaries (allowed value 0=no periodic boundary 1 = one periodic boundary
    do i = 1, cfg%nprdbnd
      read (unit, *) cfg%prdbndclr(1, i), cfg%prdbndclr(2, i), cfg%prdbndclr(3, i) !< corresponding color pairs (indices 1 and 2)
      !< of each periodic boundaries and equal coordinate
      !< of periodic boundary (index 3 1=x 2=y)
    end do
    read (unit, *) cfg%flt_dspeed !< filter on discontinuity speeds [0,1.0] (0=disactive)
    read (unit, *) cfg%imtf       !< iteration of mesh topology freezing (this value must be equal to one backup iteration)
  end function read_config

end module mod_config
