! Find the cells crossed by the shock and the phantom points

subroutine fnd_phps(nface,&
&ibndfac,&       ! not used
&nbfac,&
&icelnod,&
&nvt,&
&nelem,&
&xy,&
&xysh,&
&nodcod,&
&npoin,&
&inodptr,&
&nbpoin,&
&nshocks,&
&nshockpoints,&  ! not used
&nshockedges,&
&nphpoin,&
&pmap)

  use mod_kinds, only: wp, i4
  implicit none(type, external)
  include 'paramt.h'

  integer(i4) nface, nelem, npoin, nvt, nbfac, nbpoin
  integer(i4) nshocks, nshockedges(nshmax), nshockpoints(nshmax)

!     .. array arguments ..
  real(wp) xy(ndim, npoin), xysh(ndim, npshmax, nshmax)
  integer(i4) ibndfac(3, nbfac + neshmax), inodptr(nbpoin, 3),&
  &icelnod(nvt, nelem),&
  &nodcod(npoin),&
  &iedgptr(3, nface),&
  &pmap(npoin)

!     character*(*) fname
  integer(i4) inode(2)

!     .. local scalars ..
  real(wp) xc1, xc2, xc3, yc1, yc2, yc3, xs1, xs2, ys1, ys2, rdshp
  real(wp) d1, d2, d3
  integer(i4) i, k, n1, n2, n3, ielem, ielemsh, ishel1, ishel2, ii, nphpoin,&
  &ifail, nbphp, ish, iface, last, ipos, ipoin

! open log file
  open (8, file='log/fnd_phps.log')

! set the minimum distance between shock and nodal point. The distance is normalized
! using the shock element lenght
!
!      sndmin=0.20
!      sndmin=0.40

! set to 0 the number of phantom nodes and the color of the  phantom nodes
! is set to 0
  nphpoin = 0

!      find mesh cells crossed by the shock
  do ish = 1, nshocks
    do ielemsh = 1, nshockedges(ish)
      do ielem = 1, nelem

        n1 = icelnod(1, ielem)
        n2 = icelnod(2, ielem)
        n3 = icelnod(3, ielem)
        xc1 = xy(1, n1)
        xc2 = xy(1, n2)
        xc3 = xy(1, n3)
        yc1 = xy(2, n1)
        yc2 = xy(2, n2)
        yc3 = xy(2, n3)

        xs1 = xysh(1, ielemsh, ish)
        ys1 = xysh(2, ielemsh, ish)
        xs2 = xysh(1, ielemsh + 1, ish)
        ys2 = xysh(2, ielemsh + 1, ish)

! if both the functions ishel1 and ishel2 return 0, then shock segment denoted
! by the two shock points (xs1, ys1, xs2, ys2) crosses the cell with
! the vertices xc1, yc1, xc2, yc2 ,xc3, yc3
        i = ishel1(xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1, xs2, ys2)
        if (i .eq. 0) then
          ii = ishel2(xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1, xs2, ys2)
          if (ii .eq. 0) then

! for each triangle crossed by the shock element compute the distance of
! vertices from the shock straight line.
            d1 = rdshp(xc1, yc1, xs1, ys1, xs2, ys2, sndmin)
            d2 = rdshp(xc2, yc2, xs1, ys1, xs2, ys2, sndmin)
            d3 = rdshp(xc3, yc3, xs1, ys1, xs2, ys2, sndmin)

! if distance is too small the node become a phantom node
            if (d1 .ge. 0 .and. d1 .lt. sndmin .and. nodcod(n1) .eq. 0) nodcod(n1) = -1
            if (d1 .ge. 0 .and. d1 .lt. sndmin .and. nodcod(n1) .gt. 0) nodcod(n1) = -2
            if (d2 .ge. 0 .and. d2 .lt. sndmin .and. nodcod(n2) .eq. 0) nodcod(n2) = -1
            if (d2 .ge. 0 .and. d2 .lt. sndmin .and. nodcod(n2) .gt. 0) nodcod(n2) = -2
            if (d3 .ge. 0 .and. d3 .lt. sndmin .and. nodcod(n3) .eq. 0) nodcod(n3) = -1
            if (d3 .ge. 0 .and. d3 .lt. sndmin .and. nodcod(n3) .gt. 0) nodcod(n3) = -2

          end if
        end if
      end do
    end do
  end do

! check on the periodic boundaries
! if one node on a pariodic boundary is phantom
! also the corresponding point on the other periodic boundary
! must be set phantom
  do k = 1, npoin
    if (nodcod(k) .eq. -2 .and. pmap(k) .ne. 0) then
      nodcod(pmap(k)) = -2
    end if
  end do

  nphpoin = 0
  nbphp = 0
  do k = 1, npoin
    if ((nodcod(k) .eq. -1) .or.&
    &(nodcod(k) .eq. -2)) nphpoin = nphpoin + 1
    if (nodcod(k) .eq. -2) nbphp = nbphp + 1
    if (nodcod(k) .eq. -1)&
    &write (8, *) 'node ', k, ' has become a phantom'
    if (nodcod(k) .eq. -2)&
    &write (8, *) 'node ', k, ' on the bndry has become a phantom'
  end do

  write (8, *) 'number of phantom nodes (incl. those on the bndry)',&
  &nphpoin
  if (nbphp .gt. 0) then
    write (8, *) 'uh! oh! there are ', nbphp, ' phantom nodes&
    &         on the boundary'
  end if

  close (8)

!     This part of the routine updates the boundary data
!     in order to keep into account for eventual phantom nodes
!     on the boundary
!
!     NDIM  is the space dimension =2
!     NVT = NDIM+1 is the number of vertices
!
!     NBFAC  are the boundary faces (shock segments excluded)
!     NBPOIN are the boundary points (shock points excluded)

! open log file
  open (8, file='log/ChangeBndryPtr.log')

  write (8, *) 'subr changebndryptr; nbfac was = ', nbfac
  write (8, *) 'subr changebndryptr; nbpoin was = ', nbpoin
  write (8, *) 'subr changebndryptr; npoin was = ', npoin

  do 1 ipoin = 1, npoin
1   continue
    write (8, *) 'Subr ChangeBndryPtr; NBFAC is now = ', NBFAC
! part before commented

!     call x04eaf('general',' ',3,nbfac,ibndptr,3,
!    +            'bndry pointer in changebndryptr',ifail)
!     pause
!     stop

    close (8)

    return
    end subroutine fnd_phps

! the function ishel1 returns 0 if the cell is crossed by the straight line
! passing for the two shock points. Otherwise the function ishel1 returns 1
!
! In particular the function evaluates the sign of results obtained
! replacing the vertex coordinates in the equation of the straight line
! denoted by the two shock points. If all results have the same sign then
! the straight line does not cross the triangle.

    integer(i4) function ishel1(xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1,&
    &xs2, ys2)

      use mod_kinds, only: wp, i4
      real(wp) xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1, xs2, ys2
      real(wp) eval, x, y

      eval(x, y) = (xs1 - x)*(ys1 - ys2) - (ys1 - y)*(xs1 - xs2)

      ishel1 = sign(1.0d+0, eval(xc1, yc1)) +&
      &sign(1.0d+0, eval(xc2, yc2)) +&
      &sign(1.0d+0, eval(xc3, yc3))
      ishel1 = abs(ishel1)/3
      return

    end function ishel1

! the function ishel2 is applied only to the cells where the function
! ishel1 return 0.
! the function ishel2 returns 0 if ?almost? one cell segment has the
! intersection point with the shock straight line enclosed between the
! two shock points

    integer(i4) function ishel2(xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1,&
    &xs2, ys2)

      use mod_kinds, only: wp, i4
      real(wp) xc1, yc1, xc2, yc2, xc3, yc3, xs1, ys1, xs2, ys2
      integer(i4) nn
      parameter(nn=2)
      real(wp) a(nn, nn), b(nn), x(nn)
      real(wp) xi, yi, rlsh2, rl2

      ishel2 = 0

! write the equation of the straight line passing for the two shock points
      a(1, 1) = (ys2 - ys1)
      a(1, 2) = (xs1 - xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1) = xs2*(ys2 - ys1) + ys2*(xs1 - xs2)

! write the equation of the straight line perpendicular to previous shock line
! and passing for the vertex node
      a(2, 1) = (yc2 - yc1)
      a(2, 2) = (xc1 - xc2)
!     b(2)=xc1*(yc2-yc1)+yc1*(xc1-xc2)
      b(2) = xc2*(yc2 - yc1) + yc2*(xc1 - xc2)

! solve the linear system and find the intersection point coordinates
      call solg(nn, nn, a, b, x)
      xi = x(1)
      yi = x(2)

! compute the distance between the shock points and the sum of distances
! of the intersection point from the both shock points

      rlsh2 = (xs1 - xs2)**2 + (ys1 - ys2)**2
      rl2 = (xs1 - xi)**2 + (ys1 - yi)**2 +&
      &(xs2 - xi)**2 + (ys2 - yi)**2

! if the distance between the shock points is equal to the sum of distances
! of the intersection point from the both shock points the intersection
! point is enclosed between the two shock point
!     if(rlsh2.ge.rl2) return
      if ((rlsh2 - rl2) .ge. -1e-5) return

! repeat the same algorithm for other triangle vertices

      a(1, 1) = (ys2 - ys1)
      a(1, 2) = (xs1 - xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1) = xs2*(ys2 - ys1) + ys2*(xs1 - xs2)

      a(2, 1) = (yc3 - yc1)
      a(2, 2) = (xc1 - xc3)
!     b(2)=xc1*(yc3-yc1)+yc1*(xc1-xc3)
      b(2) = xc3*(yc3 - yc1) + yc3*(xc1 - xc3)

      call solg(nn, nn, a, b, x)

      xi = x(1)
      yi = x(2)
      rlsh2 = (xs1 - xs2)**2 + (ys1 - ys2)**2
      rl2 = (xs1 - xi)**2 + (ys1 - yi)**2 +&
      &(xs2 - xi)**2 + (ys2 - yi)**2

!     if(rlsh2.ge.rl2) return
      if ((rlsh2 - rl2) .ge. -1e-5) return

      a(1, 1) = (ys2 - ys1)
      a(1, 2) = (xs1 - xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1) = xs2*(ys2 - ys1) + ys2*(xs1 - xs2)

      a(2, 1) = (yc3 - yc2)
      a(2, 2) = (xc2 - xc3)
!     b(2)=xc2*(yc3-yc2)+yc2*(xc2-xc3)
      b(2) = xc3*(yc3 - yc2) + yc3*(xc2 - xc3)

      call solg(nn, nn, a, b, x)

      xi = x(1)
      yi = x(2)
      rlsh2 = (xs1 - xs2)**2 + (ys1 - ys2)**2
      rl2 = (xs1 - xi)**2 + (ys1 - yi)**2 +&
      &(xs2 - xi)**2 + (ys2 - yi)**2

!     if(rlsh2.ge.rl2) return
      if ((rlsh2 - rl2) .ge. -1e-5) return

      ishel2 = 1
      return
    end function ishel2

    real(wp) function rdshp(xc, yc, xs1, ys1, xs2, ys2, Smin)

      use mod_kinds, only: wp, i4
      integer(i4) nn
      parameter(nn=2)
      real(wp) a(nn, nn), b(nn), x(nn)
      real(wp) xi, yi, rlsh2, rl2, Smin
      real(wp) xc, yc, xs1, ys1, xs2, ys2
      real(wp) rlsh3

      rdshp = -1.0

      a(1, 1) = (ys2 - ys1)
      a(1, 2) = (xs1 - xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1) = xs2*(ys2 - ys1) + ys2*(xs1 - xs2)
      a(2, 1) = a(1, 2)
      a(2, 2) = -a(1, 1)
      b(2) = a(2, 1)*xc + a(2, 2)*yc

      call solg(nn, nn, a, b, x)

      xi = x(1)
      yi = x(2)
      rlsh2 = (xs1 - xs2)**2 + (ys1 - ys2)**2
!     rlsh2=0.005**2
!ren ccccccccccccccccccccccccccccccccccc
! added line
      rlsh3 = ((1.0d0 + Smin)**2 + Smin**2)*rlsh2
!ren cccccccccccccccccccccccccccccccccc
      rl2 = (xs1 - xi)**2 + (ys1 - yi)**2 +&
      &(xs2 - xi)**2 + (ys2 - yi)**2
!ren ccccccccccccccccccccccccccccccccc
! fixed line
!     if(rlsh2.gt.rl2)return
      if (rlsh3 .lt. rl2) return
!ren ccccccccccccccccccccccccccccccccc
      rdshp = (xi - xc)**2 + (yi - yc)**2
      rdshp = rdshp/rlsh2
      rdshp = sqrt(rdshp)
      return
    end function rdshp

    subroutine solg(n, nmax, a, b, x)
!     subroutine : gauss method for the solution of a
!                  linear algebraic system
      use mod_kinds, only: wp, i4
      implicit none(type, external)
      real(wp) a, b, x
      integer(i4) n, nmax

      real(wp) summ, pik, app
      integer(i4) i, j, k, imax, j1, i1

      dimension a(nmax, nmax), b(nmax), x(nmax)
!     Triangularization of matrix A with partial pivot
      do k = 1, n - 1
        imax = k
        do i1 = k + 1, n
          if (abs(a(i1, k)) .gt. abs(a(imax, k))) imax = i1
        end do
        if (imax .ne. k) then
          do j1 = k, n
            app = a(imax, j1)
            a(imax, j1) = a(k, j1)
            a(k, j1) = app
          end do
          app = b(k)
          b(k) = b(imax)
          b(imax) = app
        end if
        do i = k + 1, n
          pik = a(i, k)/a(k, k)
          a(i, k) = 0.0d+00
          b(i) = b(i) - pik*b(k)
          do j = k + 1, n
            a(i, j) = a(i, j) - pik*a(k, j)
          end do
        end do
      end do
!     Calculate results with backward substitution
      x(n) = b(n)/a(n, n)
      do i = n - 1, 1, -1
        summ = 0
        do j = i + 1, n
          summ = summ + a(i, j)*x(j)
        end do
        x(i) = (b(i) - summ)/a(i, i)
      end do
      return
    end subroutine solg
