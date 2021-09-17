! Find the cells crossed by the shock and the phantom points

      subroutine fnd_phps (nface,
     +                     ibndfac,       ! not used
     +                     nbfac,
     +                     icelnod,
     +                     nvt,
     +                     nelem,
     +                     xy,
     +                     xysh,
     +                     nodcod,
     +                     npoin,
     .                     inodptr,
     .                     nbpoin,
     +                     nshocks,
     +                     nshockpoints,  ! not used
     +                     nshockedges,
     +                     nphpoin,
     +                     pmap)

      implicit none
      include 'paramt.h'

      integer nface,nelem,npoin,nvt,nbfac,nbpoin
      integer nshocks,nshockedges(nshmax),nshockpoints(nshmax)

!     .. array arguments ..
      double precision xy(ndim,npoin),xysh(ndim,npshmax,nshmax)
      integer ibndfac(3,nbfac+neshmax),inodptr(nbpoin,3),
     &icelnod(nvt,nelem),
     &nodcod(npoin),
     &iedgptr(3,nface),
     &pmap(npoin)

!     character*(*) fname
      integer inode(2)

!     .. local scalars ..
      double precision xc1,xc2,xc3,yc1,yc2,yc3, xs1,xs2,ys1,ys2,rdshp
      double precision d1,d2,d3
      integer i,k,n1,n2,n3,ielem,ielemsh,ishel1,ishel2,ii,nphpoin,
     &ifail,nbphp,ish,iface,last,ipos,ipoin

! open log file
      open(8,file='log/fnd_phps.log')

! set the minimum distance between shock and nodal point. The distance is normalized
! using the shock element lenght
!
!      sndmin=0.20
!      sndmin=0.40

! set to 0 the number of phantom nodes and the color of the  phantom nodes
! is set to 0
       nphpoin=0
caldo
caldo  The following cycle is no more necessary because
caldo  before entering in this routine NODCOD is reset
caldo  to the value which had for the background grid
caldo  that is WITHOUT phantom nodes
caldo
!      do k = 1, npoin
!       if(nodcod(k).eq.-1)nodcod(k)=0
!       if(nodcod(k).eq.-2)nodcod(k)=2
!      end do

!      find mesh cells crossed by the shock
       do ish=1, nshocks
       do ielemsh=1, nshockedges(ish)
       do ielem=1, nelem

        n1=icelnod(1,ielem)
        n2=icelnod(2,ielem)
        n3=icelnod(3,ielem)
        xc1=xy(1,n1)
        xc2=xy(1,n2)
        xc3=xy(1,n3)
        yc1=xy(2,n1)
        yc2=xy(2,n2)
        yc3=xy(2,n3)

        xs1=xysh(1,ielemsh,ish)
        ys1=xysh(2,ielemsh,ish)
        xs2=xysh(1,ielemsh+1,ish)
        ys2=xysh(2,ielemsh+1,ish)

! if both the functions ishel1 and ishel2 return 0, then shock segment denoted
! by the two shock points (xs1, ys1, xs2, ys2) crosses the cell with
! the vertices xc1, yc1, xc2, yc2 ,xc3, yc3
        i=ishel1(xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1, xs2, ys2)
        if(i.eq.0) then
          ii=ishel2(xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1, xs2, ys2)
          if(ii.eq.0) then

! for each triangle crossed by the shock element compute the distance of
! vertices from the shock straight line.
            d1=rdshp(xc1, yc1,xs1, ys1, xs2, ys2,sndmin)
            d2=rdshp(xc2, yc2,xs1, ys1, xs2, ys2,sndmin)
            d3=rdshp(xc3, yc3,xs1, ys1, xs2, ys2,sndmin)

! if distance is too small the node become a phantom node
        if(d1.ge.0.and.d1.lt.sndmin.and.nodcod(n1).eq.0) nodcod(n1)=-1
        if(d1.ge.0.and.d1.lt.sndmin.and.nodcod(n1).gt.0) nodcod(n1)=-2
caldo
        if(d2.ge.0.and.d2.lt.sndmin.and.nodcod(n2).eq.0) nodcod(n2)=-1
        if(d2.ge.0.and.d2.lt.sndmin.and.nodcod(n2).gt.0) nodcod(n2)=-2
caldo
        if(d3.ge.0.and.d3.lt.sndmin.and.nodcod(n3).eq.0) nodcod(n3)=-1
        if(d3.ge.0.and.d3.lt.sndmin.and.nodcod(n3).gt.0) nodcod(n3)=-2
caldo

          endif
         endif
       end do
       end do
       end do

! check on the periodic boundaries
! if one node on a pariodic boundary is phantom
! also the corresponding point on the other periodic boundary
! must be set phantom
        do k = 1, npoin
         if(nodcod(k).eq.-2.and.pmap(k).ne.0)then
          nodcod(pmap(k))=-2
         endif
        enddo

       nphpoin=0
       nbphp=0
       do k = 1, npoin
        if((nodcod(k).eq.-1).or.
     &     (nodcod(k).eq.-2))nphpoin=nphpoin+1
        if(nodcod(k).eq.-2)nbphp=nbphp+1
        if(nodcod(k).eq.-1)
     &   write(8,*)'node ',k,' has become a phantom'
        if(nodcod(k).eq.-2)
     &   write(8,*)'node ',k,' on the bndry has become a phantom'
       end do

       write(8,*)'number of phantom nodes (incl. those on the bndry)',
     &nphpoin
        if(nbphp .gt. 0)then
           write(8,*)'uh! oh! there are ',nbphp,' phantom nodes
     &on the boundary'
       endif

      close(8)

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
      open(8,file='log/ChangeBndryPtr.log')

      write(8,*)'subr changebndryptr; nbfac was = ',nbfac
      write(8,*)'subr changebndryptr; nbpoin was = ',nbpoin
      write(8,*)'subr changebndryptr; npoin was = ',npoin

      do 1 ipoin = 1, npoin

!     if the node was de-activated (NODCOD=-2)
!     it is searched in the INODPTR vector

! previously commented
!          if( nodcod(ipoin) .eq. -2 )then
!              call binsrc(ipoin,inodptr(1,1),nbpoin,ipos,last)
!              if(ipos.eq.0)then
!                  write(8,*)'entry not found for ',ipoin
!                  write(*,*)'entry not found for ',ipoin
!                  stop
!               endif
!               write(8,*)'removing node ',ipoin,
!     &                   ' belongs to edges ',(inodptr(ipos,k),k=2,3)
!               do 6 k = 2,3
!                  iface = inodptr(ipos,k)
!                  if( ibndfac(1,iface) .ne. ipoin )then
!                      inode(k-1) = ibndfac(1,iface)
!                  else
!                      inode(k-1) = ibndfac(2,iface)
!                  endif
!    6          continue
!    the face is modified
!            iface = inodptr(ipos,2)
!            ibndfac(1,iface) = inode(1)
!            ibndfac(2,iface) = inode(2)
!            write(8,*)'face ',iface,' has been updated with ',
!     &(inode(k),k=1,2)
!aldo
!            iface = inodptr(ipos,3)
!            write(8,*)'face ',iface,' has been removed'
!            ibndfac(1,iface) = ibndfac(1,iface)
!            ibndfac(2,iface) = ibndfac(2,iface)
!            ibndfac(3,iface) =-ibndfac(3,iface)
!         endif
    1 continue
      write(8,*)'Subr ChangeBndryPtr; NBFAC is now = ',NBFAC
! part before commented

!     call x04eaf('general',' ',3,nbfac,ibndptr,3,
!    +            'bndry pointer in changebndryptr',ifail)
!     pause
!     stop

      close(8)

      return
      end

! the function ishel1 returns 0 if the cell is crossed by the straight line
! passing for the two shock points. Otherwise the function ishel1 returns 1
!
! In particular the function evaluates the sign of results obtained
! replacing the vertex coordinates in the equation of the straight line
! denoted by the two shock points. If all results have the same sign then
! the straight line does not cross the triangle.

      integer function ishel1(xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1,
     +                       xs2, ys2)

      double precision xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1,xs2, ys2
      double precision eval,x,y

      eval(x,y)=(xs1-x)*(ys1-ys2)-(ys1-y)*(xs1-xs2)

      ishel1=sign(1.0d+0,eval(xc1, yc1))+
     +       sign(1.0d+0,eval(xc2, yc2))+
     +       sign(1.0d+0,eval(xc3, yc3))
      ishel1=abs(ishel1)/3
      return

      end

! the function ishel2 is applied only to the cells where the function
! ishel1 return 0.
! the function ishel2 returns 0 if ?almost? one cell segment has the
! intersection point with the shock straight line enclosed between the
! two shock points

      integer function ishel2(xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1,
     +                       xs2, ys2)

      double precision xc1, yc1, xc2, yc2 ,xc3, yc3, xs1, ys1,xs2, ys2
      integer nn
      parameter (nn=2)
      double precision a(nn,nn), b(nn),x(nn)
      double precision xi, yi , rlsh2,rl2

      ishel2=0

! write the equation of the straight line passing for the two shock points
      a(1,1)=(ys2-ys1)
      a(1,2)=(xs1-xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1)=xs2*(ys2-ys1)+ys2*(xs1-xs2)

! write the equation of the straight line perpendicular to previous shock line
! and passing for the vertex node
      a(2,1)=(yc2-yc1)
      a(2,2)=(xc1-xc2)
!     b(2)=xc1*(yc2-yc1)+yc1*(xc1-xc2)
      b(2)=xc2*(yc2-yc1)+yc2*(xc1-xc2)

! solve the linear system and find the intersection point coordinates
      call solg(nn,nn,a,b,x)
      xi=x(1)
      yi=x(2)


! compute the distance between the shock points and the sum of distances
! of the intersection point from the both shock points

      rlsh2=(xs1-xs2)**2+(ys1-ys2)**2
      rl2  =(xs1-xi )**2+(ys1-yi )**2+
     +      (xs2-xi )**2+(ys2-yi )**2

! if the distance between the shock points is equal to the sum of distances
! of the intersection point from the both shock points the intersection
! point is enclosed between the two shock point
!     if(rlsh2.ge.rl2) return
      if((rlsh2-rl2).ge.-1e-5) return

! repeat the same algorithm for other triangle vertices

      a(1,1)=(ys2-ys1)
      a(1,2)=(xs1-xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1)=xs2*(ys2-ys1)+ys2*(xs1-xs2)

      a(2,1)=(yc3-yc1)
      a(2,2)=(xc1-xc3)
!     b(2)=xc1*(yc3-yc1)+yc1*(xc1-xc3)
      b(2)=xc3*(yc3-yc1)+yc3*(xc1-xc3)

      call solg(nn,nn,a,b,x)

      xi=x(1)
      yi=x(2)
      rlsh2=(xs1-xs2)**2+(ys1-ys2)**2
      rl2  =(xs1-xi )**2+(ys1-yi )**2+
     +      (xs2-xi )**2+(ys2-yi )**2

!     if(rlsh2.ge.rl2) return
      if((rlsh2-rl2).ge.-1e-5) return

      a(1,1)=(ys2-ys1)
      a(1,2)=(xs1-xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1)=xs2*(ys2-ys1)+ys2*(xs1-xs2)

      a(2,1)=(yc3-yc2)
      a(2,2)=(xc2-xc3)
!     b(2)=xc2*(yc3-yc2)+yc2*(xc2-xc3)
      b(2)=xc3*(yc3-yc2)+yc3*(xc2-xc3)

      call solg(nn,nn,a,b,x)

      xi=x(1)
      yi=x(2)
      rlsh2=(xs1-xs2)**2+(ys1-ys2)**2
      rl2  =(xs1-xi )**2+(ys1-yi )**2+
     +      (xs2-xi )**2+(ys2-yi )**2

!     if(rlsh2.ge.rl2) return
      if((rlsh2-rl2).ge.-1e-5) return

      ishel2=1
      return
      end


      double precision function rdshp(xc,yc,xs1,ys1,xs2,ys2,Smin)

      integer nn
      parameter (nn=2)
      double precision a(nn,nn), b(nn),x(nn)
      double precision xi, yi, rlsh2,rl2,Smin
      double precision xc, yc, xs1, ys1, xs2, ys2
      double precision rlsh3

      rdshp=-1.0

      a(1,1)=(ys2-ys1)
      a(1,2)=(xs1-xs2)
!     b(1)=xs1*(ys2-ys1)+ys1*(xs1-xs2)
      b(1)=xs2*(ys2-ys1)+ys2*(xs1-xs2)
      a(2,1)=a(1,2)
      a(2,2)=-a(1,1)
      b(2)= a(2,1)*xc+a(2,2)*yc

      call solg(nn,nn,a,b,x)

      xi=x(1)
      yi=x(2)
      rlsh2=(xs1-xs2)**2+(ys1-ys2)**2
!     rlsh2=0.005**2
!ren ccccccccccccccccccccccccccccccccccc
! added line
      rlsh3=((1.0d0+Smin)**2+Smin**2)*rlsh2
!ren cccccccccccccccccccccccccccccccccc
      rl2  =(xs1-xi )**2+(ys1-yi )**2+
     +      (xs2-xi )**2+(ys2-yi )**2
!ren ccccccccccccccccccccccccccccccccc
! fixed line
!     if(rlsh2.gt.rl2)return
      if(rlsh3.lt.rl2)return
!ren ccccccccccccccccccccccccccccccccc
      rdshp=(xi-xc)**2+(yi-yc)**2
      rdshp=rdshp/rlsh2
      rdshp=sqrt(rdshp)
      return
      end


      subroutine solg(n,nmax,a,b,x)
!     subroutine : gauss method for the solution of a
!                  linear algebraic system
      implicit none
      double precision a,b,x
      integer n, nmax

      double precision summ,pik,app
      integer i,j,k,imax,j1,i1

      dimension a(nmax,nmax),b(nmax),x(nmax)
!     Triangularization of matrix A with partial pivot
      do 10 k=1,n-1
      imax=k
      do 3 i1=k+1,n
      if (abs(a(i1,k)).gt.abs(a(imax,k))) imax=i1
3     continue
      if (imax.eq.k) goto 7
      do 5 j1=k,n
      app=a(imax,j1)
      a(imax,j1)=a(k,j1)
5     a(k,j1)=app
      app=b(k)
      b(k)=b(imax)
      b(imax)=app
7     do 10 i=k+1,n
      pik=a(i,k)/a(k,k)
      a(i,k)=0.0d+00
      b(i)=b(i)-pik*b(k)
      do 10 j=k+1,n
10    a(i,j)=a(i,j)-pik*a(k,j)
!     Calculate results with backward substitution
      x(n)=b(n)/a(n,n)
      do 30 i=n-1,1,-1
      summ=0
      do 20 j=i+1,n
20    summ=summ+a(i,j)*x(j)
30    x(i)=(b(i)-summ)/a(i,i)
      return
      end
