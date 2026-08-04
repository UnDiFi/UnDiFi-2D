! Compute a shock point

module co_shock_ctx_m
! Phase 3.2 increment 1 (ROADMAP.md #14): closure payload for f's
! residual evaluation, replacing the COMMON/shck/ block co_shock used
! to communicate rom/pm/um/gam/delta/R2 into f, now that both go
! through mod_newton_solve's generic residual_if(i, y, ctx) interface.
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: shock_ctx_t

  type :: shock_ctx_t
    real(wp) :: rom = 0.0_wp, pm = 0.0_wp, um = 0.0_wp
    real(wp) :: gam = 0.0_wp, delta = 0.0_wp, r2 = 0.0_wp
  end type shock_ctx_t

end module co_shock_ctx_m

module co_dc_ctx_m
! Phase 3.2 increment 1 (ROADMAP.md #14): closure payload for fdc's
! residual evaluation, replacing the COMMON/dc/ block co_dc used to
! communicate gam/delta/R1/R2/S1/S2 into fdc.
  use mod_kinds, only: wp, i4
  implicit none(type, external)
  private
  public :: dc_ctx_t

  type :: dc_ctx_t
    real(wp) :: gam = 0.0_wp, delta = 0.0_wp
    real(wp) :: r1 = 0.0_wp, r2 = 0.0_wp, s1 = 0.0_wp, s2 = 0.0_wp
  end type dc_ctx_t

end module co_dc_ctx_m

subroutine co_shock(x1, x2, wshk, R14)

!     x1(1) and x2(1) upstream and downstream density
!     x1(2) and x2(2) upstream and downstream pressure
!     x1(3) and x2(3) upstream and downstream normal velocity
!     wshk shock velocity

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_newton_solve, only: newton_solve, residual_if
  use co_shock_ctx_m, only: shock_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: f
  include 'paramt.h'

  real(wp) x1, x2, wshk, R14
  dimension x1(4), x2(4)
  real(wp) rov, rom, pv, pm, uv, um, gam, delta, w
  real(wp) y0(4), yn1(4)
  type(shock_ctx_t) :: ctx

! NOTE: w is used uninitialized below (yn1(4) = w), reproducing a
! pre-existing latent bug in the pre-refactor source (its "initialization
! of downstream state and shock velocity" block is entirely commented
! out, leaving w never assigned). Not this increment's job to fix -- see
! the merged-3.2 plan's "preserve and document, do not fix" policy for
! untouched latent bugs; giving it an explicit deterministic seed here
! would itself be an unreviewed numerics change.

! constants assign
  gam = GA
  delta = (gam - 1.0)/2.

! downstream state 1
  rov = x1(1)
  pv = x1(2)
  uv = x1(3)

! upstream state 2
  rom = x2(1)
  pm = x2(2)
  um = x2(3)

! calculation of invariants
  ctx%r2 = R14
  ctx%rom = rom
  ctx%pm = pm
  ctx%um = um
  ctx%gam = gam
  ctx%delta = delta

! compute the downstream state and shock velocity with Newton-Raphson method
! initialize the vector of unknowns
  y0(1) = rov
  y0(2) = pv
  y0(3) = uv
  y0(4) = w

  call newton_solve(4_i4, y0, f, ctx, 0.2_wp, 1.0e-7_wp, 0.001_wp, yn1)

  wshk = yn1(4)
  x1(1) = yn1(1)
  x1(2) = yn1(2)
  x1(3) = yn1(3)

  return
end subroutine co_shock

! ************************************
real(wp) function f(i, y, ctx) result(r)
  use mod_kinds, only: wp, i4
  use co_shock_ctx_m, only: shock_ctx_t
  implicit none(type, external)

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  real(wp) rov, pv, uv, w, rom, pm, um, gam, delta, R2

  select type (ctx)
  type is (shock_ctx_t)
    rom = ctx%rom
    pm = ctx%pm
    um = ctx%um
    gam = ctx%gam
    delta = ctx%delta
    R2 = ctx%r2
  end select

  rov = y(1)
  pv = y(2)
  uv = y(3)
  w = y(4)

  r = 0
  if (i .eq. 1) then
    r = rov*(uv - w) - rom*(um - w)
  elseif (i .eq. 2) then
    r = pv + rov*(uv - w)**2 - pm - rom*(um - w)**2
  elseif (i .eq. 3) then
    r = gam/(gam - 1.0)*pv/rov + 0.5*(uv - w)**2&
    &- gam/(gam - 1.0)*pm/rom - 0.5*(um - w)**2

  elseif (i .eq. 4) then
    r = sqrt(gam*pv/rov) + delta*uv - R2
  end if
  return
end function f

subroutine invmat(a, b, r)
  use mod_kinds, only: wp, i4

  implicit real*8(a - h, o - z)
!     subroutine for matrix inversion
!     a     : matrix r*r to be inverted (real)
!     b     : matrix r*r inverted
!     r     : dimension (max 20) (integer)
  integer(i4) r, j, i, k, l
  real(wp) a(4, 4), b(4, 4), c(4, 4)
  logical found_pivot

  do i = 1, r
    do j = 1, r
      b(i, j) = 0.d+00
      c(i, j) = a(i, j)
    end do
  end do
  do i = 1, r
    b(i, i) = 1.d+00
  end do
  do j = 1, r
    found_pivot = .false.
    do i = j, r
      if (a(i, j) .ne. 0.) then
        found_pivot = .true.
        exit
      end if
    end do

    if (.not. found_pivot) then
      do i = 1, r
        if (a(j, i) .ne. 0.) then
          write (*, *) 'singular matrix'
          return
        end if
      end do
      cycle
    end if

    do k = 1, r
      s = a(j, k)
      a(j, k) = a(i, k)
      a(i, k) = s
      s = b(j, k)
      b(j, k) = b(i, k)
      b(i, k) = s
    end do
    t = 1/a(j, j)
    do k = 1, r
      a(j, k) = t*a(j, k)
      b(j, k) = t*b(j, k)
    end do
    do l = 1, r
      if (l .eq. j) cycle
      t = -a(l, j)
      do k = 1, r
        a(l, k) = a(l, k) + t*a(j, k)
        b(l, k) = b(l, k) + t*b(j, k)
      end do
    end do
  end do
  do i = 1, r
    do j = 1, r
      a(i, j) = c(i, j)
    end do
  end do
  return
end subroutine invmat

subroutine co_dc(x1, x2, wdc)

!     x1(1) and x2(1) upstream and downstream density
!     x1(2) and x2(2) upstream and downstream pressure
!     x1(3) and x2(3) upstream and downstream normal velocity
!     wdc  contact discontinuity velocity

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, nprdbndmax
  use mod_newton_solve, only: newton_solve, residual_if
  use co_dc_ctx_m, only: dc_ctx_t
  implicit none(type, external)
  procedure(residual_if) :: fdc
  include 'paramt.h'

  real(wp) x1, x2, wdc
  dimension x1(4), x2(4)

  real(wp) ro1, ro2, p1, p2, u1, u2, gam, delta
  real(wp) y0(7), yn1(7)
  type(dc_ctx_t) :: ctx

  ro1 = x1(1)
  p1 = x1(2)
  u1 = x1(3)
  ro2 = x2(1)
  p2 = x2(2)
  u2 = x2(3)

! assign constants
  gam = GA
  delta = (gam - 1.0)/2.

! compute invariants
  ctx%gam = gam
  ctx%delta = delta
  ctx%r1 = sqrt(gam*p1/ro1) + delta*u1
  ctx%r2 = sqrt(gam*p2/ro2) - delta*u2
  ctx%s1 = p1/ro1**gam
  ctx%s2 = p2/ro2**gam

! compute the downstream state and shock velocity with Newton-Raphson method

! initialize the vector of unknowns
  y0(1) = ro1
  y0(2) = p1
  y0(3) = u1
  y0(4) = ro2
  y0(5) = p2
  y0(6) = u2
  y0(7) = wdc

  call newton_solve(7_i4, y0, fdc, ctx, 1.0_wp, 1.0e-10_wp, 0.01_wp, yn1)

  wdc = yn1(7)
  x1(1) = yn1(1)
  x1(2) = yn1(2)
  x1(3) = yn1(3)
  x2(1) = yn1(4)
  x2(2) = yn1(5)
  x2(3) = yn1(6)

  return
end subroutine co_dc

real(wp) function fdc(i, y, ctx) result(r)
  use mod_kinds, only: wp, i4
  use co_dc_ctx_m, only: dc_ctx_t
  implicit none(type, external)

  integer(i4), intent(in) :: i
  real(wp), intent(in) :: y(:)
  class(*), intent(in) :: ctx

  real(wp) ro1, ro2, p1, p2, u1, u2, w, gam, delta
  real(wp) R1, R2, S1, S2

  select type (ctx)
  type is (dc_ctx_t)
    gam = ctx%gam
    delta = ctx%delta
    R1 = ctx%r1
    R2 = ctx%r2
    S1 = ctx%s1
    S2 = ctx%s2
  end select

  ro1 = y(1)
  p1 = y(2)
  u1 = y(3)
  ro2 = y(4)
  p2 = y(5)
  u2 = y(6)
  w = y(7)

  r = 0.d+0
  if (i .eq. 1) then
    r = sqrt(gam*p1/ro1) + delta*u1 - R1
  elseif (i .eq. 2) then
    r = p1/ro1**gam - S1
  elseif (i .eq. 3) then
    r = sqrt(gam*p2/ro2) - delta*u2 - R2
  elseif (i .eq. 4) then
    r = p2/ro2**gam - S2
  elseif (i .eq. 5) then
    r = p1 - p2
  elseif (i .eq. 6) then
    r = u1 - u2
  elseif (i .eq. 7) then
    r = w - u1
  end if

  return
end function fdc
