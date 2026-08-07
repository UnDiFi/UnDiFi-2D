! Updates values in the phantom nodes of the background mesh(0)
! interpolating values in the shocked mesh (1)

!    .           work,
subroutine interp(&
&ibndfac,&        !not used
&nbfac,&          !not used
&icelnod,&
&nvt,&
&nelem,&
&xy,&
&zroe,&
&xysh,&
&xyshu,&
&xyshd,&
&nphpoin,&        !not used
&xybkg,&
&zbkg,&
&nodcod,&
&npoin,&
&nshocks,&
&nshockpoints,&
&ia,&
&ja,&
&iclr,&
&nclr)

  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, nprdbndmax, npshmax, nshmax
  use mod_search, only: bin_grid_t, build_bin_grid, query_point
  implicit none(type, external)
  external finder, finder_search
  include 'paramt.h'

!     .. scalar arguments ..
  integer(i4) nelem, nvt, nbfac
  integer(i4) nshocks, nshockpoints(nshmax), nphpoin

!     .. array arguments ..
!    &                 work(ndim,npshmax,*),
  real(wp) xysh(ndim, npshmax, *),&
  &xyshu(ndim, npshmax, *),&
  &xyshd(ndim, npshmax, *),&
  &xy(ndim, *),&
  &zroe(ndof, *),&
  &xybkg(ndim, *),&
  &zbkg(ndof, *)

  integer(i4) ibndfac(3, *),&
  &icelnod(nvt, nelem),&
  &nodcod(*),&
  &npoin(0:*)

  integer(i4) nclr
  integer(i4) ia(*), ja(*), iclr(nclr)

!     .. local scalars ..
  integer(i4) ipoin, ielem, i, ii, k, n, ifail, ish, clr, bbgn, bend, j, kp1, ibc
  real(wp) x0, y0, x1, y1, x2, y2, dum, dum1, dum2

!     Phase 4.1 (ROADMAP.md #16): spatial acceleration structure --
!     narrows finder's candidate list instead of the brute-force
!     do ielem=1,nelem scan (see finder_search below). Rebuilt once per
!     call since the mesh changes every outer iteration.
  type(bin_grid_t) :: grid
  integer(i4), allocatable :: cand(:)

!     open log file
  open (8, file='log/interp.log')

  write (8, *) 'enter in interp'

!     make upstream node coordinates coincide with downstream node coordinates
!     these are the coordinates of the shocked mesh (1)
!     Note: the nof of shock points is that on the shocked mesh
!     not the one on the background mesh since this one might have
!     been updated in the shock redistribution routine called by shockmov
  do ish = 1, nshocks
    do ipoin = 1, nshockpoints(ish)
      xyshu(1, ipoin, ish) = xysh(1, ipoin, ish)
      xyshu(2, ipoin, ish) = xysh(2, ipoin, ish)
      xyshd(1, ipoin, ish) = xysh(1, ipoin, ish)
      xyshd(2, ipoin, ish) = xysh(2, ipoin, ish)
    end do
  end do

!     interpolate in the background mesh nodes, using the connectivity of the shocked mesh;
!     the interpolation is necessary only for the ghost nodes
  call build_bin_grid(grid, icelnod, nvt, xy, nelem)
  do ipoin = 1, npoin(0)
!       if((nodcod(ipoin).eq.-1) or.(nodcod(ipoin).eq.-2)) then
    if ((nodcod(ipoin) .eq. -1)) then
      write (8, *) 'trying to locate ', ipoin,&
      &'(', xybkg(1, ipoin), ',', xybkg(2, ipoin), ')'
      ifail = 0
!            if(ipoin.eq.1449)ifail=99
      call query_point(grid, xybkg(1, ipoin), xybkg(2, ipoin), cand)
      call finder_search(icelnod, xy, ndim, zroe, ndof, xybkg(1, ipoin),&
      &zbkg(1, ipoin), cand, size(cand), ielem, ifail)
      write (8, *) 'found in cell ', ielem, ifail
      if (ifail .ne. 0) then
        write (8, *) 'search failed for vertex ', ipoin
        write (8, *) (xybkg(i, ipoin), i=1, ndim)
        write (8, *) 'cell no is ', ielem
        stop
      end if
    else
!           call dcopy(ndof,zroe(1,ipoin),1,zbkg(1,ipoin),1)
      do i = 1, ndof
        zbkg(i, ipoin) = zroe(i, ipoin)
      end do
    end if
  end do

!     interpolate phantom node on the boundary
  do ipoin = 1, npoin(0)
    if ((nodcod(ipoin) .eq. -2)) then
      write (8, *) 'trying to locate bounday point', ipoin,&
      &'(', xybkg(1, ipoin), ',', xybkg(2, ipoin), ')'

      x0 = xybkg(1, ipoin)
      y0 = xybkg(2, ipoin)
!        write(*,*)'x0,y0:',x0,y0

      do i = 1, nbfac

        k = ibndfac(1, i)
        kp1 = ibndfac(2, i)
        ibc = ibndfac(3, i)
        if (ibc .lt. 10 .and. ibc .gt. 0) then

          x1 = xy(1, k)
          y1 = xy(2, k)
          x2 = xy(1, kp1)
          y2 = xy(2, kp1)

          dum = ((x2 - x1)**2 + (y2 - y1)**2)
          dum1 = ((x0 - x1)**2 + (y0 - y1)**2)
          dum2 = ((x0 - x2)**2 + (y0 - y2)**2)

          if ((dum1 + dum2) .le. dum) then
            write (8, *) 'search succesfully for vertex ', ipoin
            write (8, *) (xybkg(ii, ipoin), ii=1, ndim)

            dum = sqrt(dum1) + sqrt(dum2)
            dum1 = sqrt(dum1)/dum
            dum2 = sqrt(dum2)/dum

            write (8, *) 'dum:', dum, 'dum1:', dum1, 'dum2:', dum2
            do j = 1, ndof
              zbkg(j, ipoin) = dum2*zroe(j, k) + dum1*zroe(j, kp1)
              write (8, *) 'k:', k, 'kp1:', kp1
              write (8, *) 'z', j, ':', zbkg(j, ipoin), zroe(j, k), zroe(j, kp1)

            end do

          end if
        end if

      end do

    end if

  end do

!    goto 65
!    do ish=1,nshocks
!     do ipoin = 1, nshockpoints(ish)
!        k = ipoin ! shock point
!        n = k + nshockpoints(ish) ! duplicated shock point
!        xysh(1,n,ish) = work(1,ipoin,ish)
!        xysh(2,n,ish) = work(2,ipoin,ish)
!     enddo
!     enddo
65 continue

  write (8, *) 'exit interp'

  close (8)

  return
end subroutine interp

subroutine finder(icelnod, nelem, coor, ndim, zroe, ndof, xyin, zout,&
&ielem, info)

!     input:
!            xyin nodal coords (belonging to the background grid)
!                 to be located inside the current grid
!            icelnod cell to node pointer
!            coor nodal coordinates
!            ndim space dimension
!            ndof nof degrees of freedom
!     output:
!            ielem is the cell node xyin falls inside
!            info = 0 node found !=0 search failed
!            zout(*) is filled with the interpolated value

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  integer(i4) ielem, ndim, ndof, info, nelem, info1
  integer(i4) icelnod(3, *)
  real(wp) coor(ndim, *), zroe(ndof, *), aa
  real(wp) xyin(ndim), zout(ndof)
  real(wp) x0, y0, xp(4), yp(4), a(3), t, s, help
  integer(i4) idxs(3), iv, ipoin, ivar, ilog
  real(wp) eps
  parameter(eps=1.e-08, ilog=1)
  real(wp) area
  integer(i4) icycl
  external area, icycl

  info1 = info
  info = 0

!     x0,y0 are the coordinates of the point to be located
  x0 = xyin(1)
  y0 = xyin(2)

  do 1000 ielem = 1, nelem
    do 50 iv = 1, 3
      ipoin = icelnod(iv, ielem) ! global node number
      xp(iv) = coor(1, ipoin)
      yp(iv) = coor(2, ipoin)
50    continue
      xp(4) = x0 ! node to be located
      yp(4) = y0 ! node to be located

      idxs(1) = 1
      idxs(2) = 2
      idxs(3) = 3
      help = 1.d0/area(xp, yp, 3, idxs)
      idxs(3) = 4 ! node to be located

!        compute area coordinates (in the x-y plane)
      do iv = 1, 3
        idxs(1) = icycl(1 + iv, 3)
        idxs(2) = icycl(2 + iv, 3)
        a(iv) = area(xp, yp, 3, idxs)*help
!           if(abs(a(iv)).le.eps)a(iv)=0.d0
      end do

      s = min(a(1), a(2), a(3))
      t = max(a(1), a(2), a(3))
      aa = abs(a(1)) + abs(a(2)) + abs(a(3))
      if (info1 .eq. 99) then
!          if(t.le.2.0.and.s.ge.-1.0)then
!          if(aa.le.2.5.and.aa.gt.0.5)then
        write (*, *) 'x0,y0:', x0, y0
        write (*, *) 'ielem', ielem
        write (*, *) t, s, aa
        write (*, *) xp(1), yp(1)
        write (*, *) xp(2), yp(2)
        write (*, *) xp(3), yp(3)

        write (*, *)
!          endif
      end if
!        write(6,*)'area coords ',(a(j),j=1,3),' cell ',ielem,help
!        write(6,*)'y/z coords ',(xp(j),j=1,4),(yp(j),j=1,4)
      if ((s .ge. 0.d0 .and. s .le. 1.d0) .and.&
      &(t .ge. 0.d0 .and. t .le. 1.d0)) then
        if (ilog .eq. 0) write (6, fmt=*) (icelnod(iv, ielem), iv=1, 3)
        do ivar = 1, ndof
          zout(ivar) = 0.d0
        end do
        do iv = 1, 3
          ipoin = icelnod(iv, ielem)
          help = a(iv)
          do ivar = 1, ndof
            zout(ivar) = zout(ivar) + help*zroe(ivar, ipoin)
          end do
        end do
        if (ilog .eq. 0) then
          do ivar = 1, ndim
            write (6, fmt=400) xyin(ivar),&
            &(coor(ivar, icelnod(iv, ielem)), iv=1, 3), (a(iv), iv=1, 3), s, t
          end do
          do ivar = 1, ndof
            write (6, fmt=300) zout(ivar),&
            &(zroe(ivar, icelnod(iv, ielem)), iv=1, 3), (a(iv), iv=1, 3)
          end do
        end if
        info = 0
        return
      end if
1000  continue
      info = 1
      write (6, *) 'search failed for vertex coords ', x0, y0
      write (6, fmt=1100) (a(iv), iv=1, 3), a(1) + a(2) + a(3), s, t
      return
1100  format('a(i),s,s,t ', 6(e12.4, 1x))
300   format(7(f10.5, 1x))
400   format(9(f10.5, 1x))
      end subroutine finder

! Phase 4.1 (ROADMAP.md #16): like finder above, but scans only a given
! candidate element list (from mod_search's query_point) instead of
! 1..nelem. Used by interp's per-phantom-node loop, which calls this
! O(N_phantom) times per outer iteration -- finder itself is left
! untouched since mod_special_point.f90's sonic_interpolate_state also
! calls it directly, at most nspecpoints times per iteration, far too
! rarely for the O(nelem) scan there to matter. The debug-only
! info1==99/ilog==0 branches in finder are dropped here, not ported --
! interp.f90's call site always passes ifail=0 into finder too (the
! only line that could set it to 99 is commented out), so those
! branches are equally unreachable from this call path either way.
      subroutine finder_search(icelnod, coor, ndim, zroe, ndof, xyin, zout,&
      &candidates, ncand, ielem, info)

        use mod_kinds, only: wp, i4
        implicit none(type, external)
        integer(i4) ielem, ndim, ndof, info, ncand
        integer(i4) icelnod(3, *), candidates(ncand)
        real(wp) coor(ndim, *), zroe(ndof, *)
        real(wp) xyin(ndim), zout(ndof)
        real(wp) x0, y0, xp(4), yp(4), a(3), t, s, help
        integer(i4) idxs(3), iv, ipoin, ivar, kc
        real(wp) area
        integer(i4) icycl
        external area, icycl

        x0 = xyin(1)
        y0 = xyin(2)

        do kc = 1, ncand
          ielem = candidates(kc)
          do iv = 1, 3
            ipoin = icelnod(iv, ielem)
            xp(iv) = coor(1, ipoin)
            yp(iv) = coor(2, ipoin)
          end do
          xp(4) = x0
          yp(4) = y0

          idxs(1) = 1
          idxs(2) = 2
          idxs(3) = 3
          help = 1.d0/area(xp, yp, 3, idxs)
          idxs(3) = 4

          do iv = 1, 3
            idxs(1) = icycl(1 + iv, 3)
            idxs(2) = icycl(2 + iv, 3)
            a(iv) = area(xp, yp, 3, idxs)*help
          end do

          s = min(a(1), a(2), a(3))
          t = max(a(1), a(2), a(3))

          if ((s .ge. 0.d0 .and. s .le. 1.d0) .and.&
          &(t .ge. 0.d0 .and. t .le. 1.d0)) then
            do ivar = 1, ndof
              zout(ivar) = 0.d0
            end do
            do iv = 1, 3
              ipoin = icelnod(iv, ielem)
              help = a(iv)
              do ivar = 1, ndof
                zout(ivar) = zout(ivar) + help*zroe(ivar, ipoin)
              end do
            end do
            info = 0
            return
          end if
        end do

        info = 1
        write (6, *) 'search failed for vertex coords ', x0, y0
        write (6, fmt=1200) (a(iv), iv=1, 3), a(1) + a(2) + a(3), s, t
        return
1200    format('a(i),s,s,t ', 6(e12.4, 1x))
      end subroutine finder_search

! *************************************************************
! Updates values in the phantom nodes of the background mesh(0)
! interpolating values in the shocked mesh (1)
! *************************************************************

      subroutine interp_sp(&
      &icelnod,&
      &nvt,&
      &nelem,&
      &xy,&
      &zroe,&
      &xysh,&
      &zroeshu,& !upstream
      &zroesh,&  !downstream
      &vshnor,&
      &npoin,&
      &nshocks,&
      &nshockpoints,&
      &typesh,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps)

        use mod_kinds, only: wp, i4
        use mod_constants, only: ndim, ndof, npshmax, nshmax
        use mod_special_point, only: special_point_t
        use mod_special_point_registry, only: make_special_point, sp_unpack
        implicit none(type, external)

!     .. scalar arguments ..
        integer(i4) nelem, nvt, nbfac
        integer(i4) nshocks, nshockpoints(nshmax), nphpoin
        integer(i4) nspecpoints, shinspps(2, 5, *)

!     .. array arguments ..
        real(wp) xysh(ndim, npshmax, *),&
        &xy(ndim, *),&
        &zroesh(ndof, npshmax, *),&
        &zroeshu(ndof, npshmax, *),&
        &zroe(ndim, *)

        integer(i4) icelnod(nvt, nelem),&
        &npoin(0:*)
        real(wp) vshnor(ndim, npshmax, *)

!     .. character array arguments
        character*1 typesh(*)
        character*5 typespecpoints(*)

!     .. local scalars ..
        integer(i4) isppnts

        class(special_point_t), allocatable :: sp

!     open log file
        open (8, file='log/interp_sp.log')

        write (8, *) 'enter in interp_sp'

!     make upstream node coordinates coincide with downstream node coordinates
!     these are the coordinates of the shocked mesh (1)
!     Note: the nof of shock points is that on the shocked mesh
!     not the one on the background mesh since this one might have
!     been updated in the shock redistribution routine called by shockmov
! Phase 3.2 increment 9: a seventh, previously-undiscovered
! typespecpoints dispatch site -- this loop's single 'SP'-only branch
! (no elseif, no fatal stop for anything else) now goes through the
! same special_point_t registry as the other six, via
! interpolate_state (default no-op, overridden only by start_point_t).
        do isppnts = 1, nspecpoints
          call make_special_point(typespecpoints(isppnts), sp)
          call sp_unpack(sp, shinspps(:, :, isppnts))
          call sp%interpolate_state(icelnod, nelem, xy, zroe,&
          &xysh(:, :, 1:nshmax), zroesh(:, :, 1:nshmax),&
          &zroeshu(:, :, 1:nshmax), nshockpoints(1:nshmax))
        end do

        close (8)
        return
      end subroutine interp_sp
