! Fix and correct the nodal position in all special points

subroutine fx_dps_loc(&
&xysh,&
&xyshu,&
&xyshd,&
&zroeshu,&
&zroeshd,&
&zroeshuold,&
&zroeshdold,&
&vshnor,&
&wsh,&
&nshocks,&
&nshockpoints,&
&nshockedges,&
&typesh,&
&iter,&
&nspecpoints,&
&typespecpoints,&
&shinspps,&
&ispclr,&
&ia,&
&ja,&
&iclr,&
&nclr,&
&zroe,&          !vale
&corg,&
&shtopolchanged)!vale

  use mod_kinds, only: wp, i4
  use mod_constants, only: ndim, ndof, npshmax, nshmax
  use mod_special_point, only: special_point_t
  use mod_special_point_registry, only: make_special_point, sp_unpack
  implicit none(type, external)

!     .. scalar arguments ..
  integer(i4) iter, nshocks, nshockpoints(nshmax), nshockedges(nshmax),&
  &nspecpoints, shinspps(2, 5, *), ispclr(*)
  character*1 typesh(*)
  character*5 typespecpoints(*)

!     .. array arguments ..
  real(wp)&
  &xysh(ndim, npshmax, *),&
  &xyshu(ndim, npshmax, *),&
  &xyshd(ndim, npshmax, *),&
  &zroeshu(ndof, npshmax, *),&
  &zroeshd(ndof, npshmax, *),&
  &zroeshuold(ndof, npshmax, *),&
  &zroeshdold(ndof, npshmax, *),&
  &vshnor(ndim, npshmax, *),&
  &wsh(ndim, npshmax, *)
  integer(i4) nclr
  integer(i4) ia(*), ja(*), iclr(nclr)
  real(wp) corg(ndim, *), zroe(ndof, *)

! vale
  real(wp) xywedge(2), dum1, dum2
  logical shtopolchanged
! vale

!     .. array arguments ..
!     character*(*) fname

!     .. local scalars ..
  integer(i4) isppnts

  class(special_point_t), allocatable :: sp

!     double precision dx,dy,kine,ws,hh,wws,taux,tauy,dum,wrr,wss
!     double precision wqpx,wqpy,w,cs,dum1,f1,f3,dxr14,dyr14
!     double precision help,x1(ndof),x2(ndof),taux12,tauy12,nx,ny
!     double precision r2(npshmax,nshmax),xtpi(24),xtp(24),r23,r14
!     double precision r14old
!     double precision dumm2,dumx,dumy
!     integer i,im,iv,k,totnshockpoints
!     integer isppnts
!     integer ip,ish,ip1,ip2,ip3,ip4,ip5,ish1,ish2,ish3,ish4,ish5
!     logical flag1,ifail

!     open log file
  open (8, file='log/fx_dps_loc.log')

! Phase 3.2 increment 7: each isppnts' relocate now dispatches through
! the special_point_t registry instead of a 16-way if/elseif.
! make_special_point error-stops on any code outside the 16 handled
! here, preserving the original's final 'condition not defined' stop.
  do isppnts = 1, nspecpoints
    call make_special_point(typespecpoints(isppnts), sp)
    call sp_unpack(sp, shinspps(:, :, isppnts))
    call sp%relocate(xysh(:, :, 1:nshmax), nshockpoints(1:nshmax),&
    &isppnts, ispclr, ia, ja, iclr, nclr, corg)
  end do

!    Save old upstream state and assign recomputed shock upstream state

!       do ish=1,nshocks
!        do ip=1, nshockpoints(ish)
!         do iv=1,ndof
!          zroeshdold(iv,ip,ish)=zroeshd(iv,ip,ish)
!          zroeshuold(iv,ip,ish)=zroeshu(iv,ip,ish)
!         enddo
!        enddo
!       enddo

!      do ish=1,nshocks
!       if(typesh(ish).eq.'s'.and.iter.le.1001.and.mod(iter,50).eq.0)then
!
!        j=nshockpoints(ish)
!        xysh(1,j+1,ish)=2*xysh(1,j,ish)-xysh(1,j-1,ish)
!        xysh(2,j+1,ish)=2*xysh(2,j,ish)-xysh(2,j-1,ish)
!        zroeshu(1,j+1,ish)=zroeshu(1,j,ish)
!        zroeshu(2,j+1,ish)=zroeshu(2,j,ish)
!        zroeshu(3,j+1,ish)=zroeshu(3,j,ish)
!        zroeshu(4,j+1,ish)=zroeshu(4,j,ish)
!
!        zroeshd(1,j+1,ish)=zroeshd(1,j,ish)
!        zroeshd(2,j+1,ish)=zroeshd(2,j,ish)
!        zroeshd(3,j+1,ish)=zroeshd(3,j,ish)
!        zroeshd(4,j+1,ish)=zroeshd(4,j,ish)

!        zroeshu(1,j,ish)=zroeshu(1,j-1,ish)
!        zroeshu(2,j,ish)=zroeshu(2,j-1,ish)
!        zroeshu(3,j,ish)=zroeshu(3,j-1,ish)
!        zroeshu(4,j,ish)=zroeshu(4,j-1,ish)
!
!        zroeshd(1,j,ish)=zroeshd(1,j-1,ish)
!        zroeshd(2,j,ish)=zroeshd(2,j-1,ish)
!        zroeshd(3,j,ish)=zroeshd(3,j-1,ish)
!        zroeshd(4,j,ish)=zroeshd(4,j-1,ish)
!
!        nshockpoints(ish)=j+1
!       endif
!      enddo

  return
end subroutine fx_dps_loc
