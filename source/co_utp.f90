! Compute an unsteady triple point (TP)

module co_utp_ctx_m
! Phase 3.2 increment 1 (ROADMAP.md #14): closure payload for futp's
! residual evaluation, replacing the long positional argument list
! (a, b, r14, dxr14, dyr14, r23, unsh1, flag1) co_utp used to pass on
! every call into futp, now that both go through mod_newton_solve's
! generic residual_if(i, y, ctx) interface.
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: utp_ctx_t

  type :: utp_ctx_t
    real(wp) :: a(11, 20) = 0.0_wp
    real(wp) :: b(11) = 0.0_wp
    real(wp) :: r14 = 0.0_wp, dxr14 = 0.0_wp, dyr14 = 0.0_wp
    real(wp) :: r23 = 0.0_wp, unsh1 = 0.0_wp
    logical :: flag1 = .false.
  end type utp_ctx_t

end module co_utp_ctx_m

subroutine co_utp(y, r14, dxr14, dyr14, r23, unsh1, yn1, flag1, ifail)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_newton_solve, only: newton_solve, residual_if
  use co_utp_ctx_m, only: utp_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: futp
  include 'paramt.h'

  integer(i4) i, j, nn
  real(wp) gam, delta
  real(wp) r14, dxr14, dyr14, r23, unsh1
  real(wp) y, yn1
  logical flag1, ifail
  dimension y(20), yn1(20)
  type(utp_ctx_t) :: ctx

!     assign constants
  gam = ga
  delta = (gam - 1.0)/2.

  nn = 20

! build matrix a and vector b of the additional conditions
! which varies depending on the case

! calculation of oblique shock being the upstream state and deviation known

  do i = 1, nn - (2*4 + 1)
    do j = 1, nn
      ctx%a(i, j) = 0.d+0
    end do
    ctx%b(i) = 0.d+0
  end do

! state 1 known
  ctx%a(1, 1) = 1.0d+0
  ctx%a(2, 2) = 1.0d+0
  ctx%a(3, 3) = 1.0d+0
  ctx%a(4, 4) = 1.0d+0
  ctx%b(1) = y(1)
  ctx%b(2) = y(2)
  ctx%b(3) = y(3)
  ctx%b(4) = y(4)

! state 2 known
  ctx%a(5, 5) = 1.0d+0
  ctx%a(6, 6) = 1.0d+0
  ctx%a(7, 7) = 1.0d+0
  ctx%a(8, 8) = 1.0d+0
  ctx%b(5) = y(5)
  ctx%b(6) = y(6)
  ctx%b(7) = y(7)
  ctx%b(8) = y(8)

! value sx12 known
  ctx%a(9, 17) = 1.0d+0
  ctx%b(9) = y(17)

! p3=p4 and theta3=theta4 (rows 10,11) are disabled upstream in futp too
! (kept as zero rows) -- see futp's own comment on the i=10/i=11 branches

  ctx%r14 = r14
  ctx%dxr14 = dxr14
  ctx%dyr14 = dyr14
  ctx%r23 = r23
  ctx%unsh1 = unsh1
  ctx%flag1 = flag1

!     open log file
  open (8, file='log/co_utp.log')

! calculate the downstream state and shock velocity with the
! newton-raphson method
  call newton_solve(nn, y, futp, ctx, 0.5_wp, 1.0e-7_wp, 0.001_wp, yn1,&
  &ifail=ifail, log_unit=8_i4, log_iter=.true., detect_divergence=.true.)

  write (8, *) 'initial and final state'
  do i = 1, nn
    write (8, *) i, y(i), yn1(i)
  end do
  close (8)

end subroutine co_utp

real(wp) function futp(i, y, ctx) result(r)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use co_utp_ctx_m, only: utp_ctx_t
  implicit none(type, external)
  include 'paramt.h'

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  integer(i4) ii, j
  real(wp) a(11, 20), b(11), r14, dxr14, dyr14, r23, unsh1
  real(wp) wn, wt, wsh
  real(wp) ro1, ro2, p1, p2, u1, u2, gam, delta
  real(wp) v1, v2, un1, un2, ut1, ut2
  real(wp) taux12, tauy12, nx12, ny12, sx12
  real(wp) taux23, tauy23, nx23, ny23, sx23
  real(wp) taux14, tauy14, nx14, ny14, sx14
  real(wp) sx23corr
  logical flag1

  select type (ctx)
  type is (utp_ctx_t)
    a = ctx%a
    b = ctx%b
    r14 = ctx%r14
    dxr14 = ctx%dxr14
    dyr14 = ctx%dyr14
    r23 = ctx%r23
    unsh1 = ctx%unsh1
    flag1 = ctx%flag1
  end select

!     assign constants and variables
  gam = ga
  delta = (gam - 1.0)/2.

  if (i .ge. 1 .and. i .le. 4) then
    ro1 = y(5)
    p1 = y(6)
    u1 = y(7)
    v1 = y(8)
    ro2 = y(9)
    p2 = y(10)
    u2 = y(11)
    v2 = y(12)
    sx12 = y(17)
    sx23 = y(18)
    wsh = y(20)

! calculation of the tangent and normal versor to shock 12 (incident shock)
    taux12 = cos(sx12)
    tauy12 = sin(sx12)
    nx12 = -tauy12
    ny12 = taux12

! calculation of the tangent and normal versor to shock 23
    taux23 = cos(sx23)
    tauy23 = sin(sx23)
    nx23 = -tauy23
    ny23 = taux23

! calculation of the normal and tangential velocity components
    un1 = u1*nx23 + v1*ny23
    un2 = u2*nx23 + v2*ny23
    ut1 = u1*taux23 + v1*tauy23
    ut2 = u2*taux23 + v2*tauy23

! calculation of the normal velocity component to shock 23 in the triple point
    wn = wsh*(taux12*nx23 + tauy12*ny23)

    wt = wsh*(taux12*taux23 + tauy12*tauy23) +&
    &unsh1*(nx12*taux23 + ny12*tauy23)

    sx23corr = 0.

    taux23 = cos(sx23 - sx23corr)
    tauy23 = sin(sx23 - sx23corr)
    nx23 = -tauy23
    ny23 = taux23

! calculation of the normal and tangential velocity components
    un1 = u1*nx23 + v1*ny23
    un2 = u2*nx23 + v2*ny23
    ut1 = u1*taux23 + v1*tauy23
    ut2 = u2*taux23 + v2*tauy23

! calculation of the normal velocity component of the triple point
    wn = wsh*(taux12*nx23 + tauy12*ny23)

    r = 0.d0
    if (i .eq. 1) then
      r = ro1*un1 - ro2*un2 - wn*(ro1 - ro2)
    elseif (i .eq. 2) then
      r = p1 + ro1*(un1 - wn)**2 - p2 - ro2*(un2 - wn)**2

    elseif (i .eq. 3) then
      r = ut1 - ut2
    elseif (i .eq. 4) then
      r = gam/(gam - 1.0)*p1/ro1 + 0.5*(un1 - wn)**2&
      &- gam/(gam - 1.0)*p2/ro2 - 0.5*(un2 - wn)**2

    end if
  elseif (i .ge. 5 .and. i .le. 8 + 1) then
    ro1 = y(1)
    p1 = y(2)
    u1 = y(3)
    v1 = y(4)
    ro2 = y(13)
    p2 = y(14)
    u2 = y(15)
    v2 = y(16)
    sx12 = y(17)
    sx14 = y(19)
    wsh = y(20)

! calculation of the tangent and normal versor to shock 12 (incident shock)
    taux12 = cos(sx12)
    tauy12 = sin(sx12)
    nx12 = -tauy12
    ny12 = taux12

! calculation of the tangent and normal versor to shock 14
    taux14 = cos(sx14)
    tauy14 = sin(sx14)
    nx14 = -tauy14
    ny14 = taux14

! calculation of the normal and tangential velocity components
    un1 = u1*nx14 + v1*ny14
    un2 = u2*nx14 + v2*ny14
    ut1 = u1*taux14 + v1*tauy14
    ut2 = u2*taux14 + v2*tauy14

! calculation of the normal velocity component of the triple point
    wn = wsh*(taux12*nx14 + tauy12*ny14)
    wt = wsh*(taux12*taux14 + tauy12*tauy14) +&
    &unsh1*(nx12*taux14 + ny12*tauy14)

    r = 0.d0
    if (i .eq. 5) then
      r = ro1*un1 - ro2*un2 - wn*(ro1 - ro2)
    elseif (i .eq. 6) then
      r = p1 + ro1*(un1 - wn)**2 - p2 - ro2*(un2 - wn)**2

    elseif (i .eq. 7) then
      r = ut1 - ut2
    elseif (i .eq. 8) then
      r = gam/(gam - 1.0)*p1/ro1 + 0.5*(un1 - wn)**2&
      &- gam/(gam - 1.0)*p2/ro2 - 0.5*(un2 - wn)**2
    elseif (i .eq. 9) then

      ro2 = y(13)
      p2 = y(14)
      u2 = y(15)
      v2 = y(16)
      sx14 = y(19)

! calculation of the normal and tangential versor
      taux14 = cos(sx14)
      tauy14 = sin(sx14)
      nx14 = -tauy14
      ny14 = taux14

      un2 = u2*dxr14 + v2*dyr14

! calculation of the normal and tangential velocity components
      r = sqrt(gam*p2/ro2) + delta*un2 - r14

    end if
  elseif (i .eq. 10) then
    ro1 = y(9)
    p1 = y(10)
    u1 = y(11)
    v1 = y(12)
    ro2 = y(13)
    p2 = y(14)
    u2 = y(15)
    v2 = y(16)
    r = p1 - p2
  elseif (i .eq. 11) then
    ro1 = y(9)
    p1 = y(10)
    u1 = y(11)
    v1 = y(12)
    ro2 = y(13)
    p2 = y(14)
    u2 = y(15)
    v2 = y(16)
    r = (u1*v2 - v1*u2)
  elseif (i .ge. 12) then
    r = 0.d0
    ii = i - (2*4 + 1 + 1 + 1)
    do j = 1, 20
      r = r + a(ii, j)*y(j)
    end do
    r = r - b(ii)
  end if

end function futp
