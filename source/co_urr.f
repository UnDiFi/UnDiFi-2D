! Compute an unsteady regular reflection point

      subroutine co_urr(y,tauwx,tauwy,yn1)

      implicit none
      include 'paramt.h'

      integer i,j,k,nn
      double precision gam,delta,a,b,bb
      double precision tauwx,tauwy
      double precision y,futp1,g,yn,yn1,g1,dum,dum1,dum2,dyn1
      double precision theta,dyn,sx14,taux14,tauy14,nx14,ny14
      double precision store13,store14,store15,store16
      logical flag1
      dimension y(15),yn(15),yn1(15),g(15,15),g1(15,15),a(10,15),b(10)
      dimension bb(15),dyn(15)

!     assign constants
      gam=ga
      delta=(gam-1.0)/2.

      nn=15

! build the matrix a and the vecotr b of the additional conditions
! (which varies according to the case)

!     compute oblique shock given the upsteam state and deviation angle
      do i=1,nn-(4+1)
        do j=1,nn
          a(i,j)=0.d+0
        enddo
        b(i)=0.d+0
      enddo

!     state 1 is known
      a(1,1)=1.0d+0
      a(2,2)=1.0d+0
      a(3,3)=1.0d+0
      a(4,4)=1.0d+0
      b(1)=y(1)
      b(2)=y(2)
      b(3)=y(3)
      b(4)=y(4)

!     state 2 is known
      a(5,5)=1.0d+0
      a(6,6)=1.0d+0
      a(7,7)=1.0d+0
      a(8,8)=1.0d+0
      b(5)=y(5)
      b(6)=y(6)
      b(7)=y(7)
      b(8)=y(8)

!     value of sx12 is known
      a(9,13)=1.0d+0
      b(9)=y(13)

!     velocity of the reflection point
      a(10,15)=1.0d+0
      b(10)=y(15)

! compute downstream state and shock velocity with the Newton-Raphson method

!     initialize vector of unknowns
      do i=1,nn
        yn1(i) = y(i)
      enddo

10    do i=1,nn
        yn(i)=yn1(i)
        bb(i)=futp1(i,yn1,a,b,tauwx,tauwy)
!       write(*,*)i,bb(i)
      enddo

!     compute jacobian
      do i=1,nn
        do j=1,nn
          do k=1,nn
            yn1(k)=yn(k)
          enddo
          dyn1=abs(yn1(j))*.01
          if(dyn1.lt.1.0d-7)dyn1=1.0d-7
          yn1(j)=yn(j)+dyn1
!         yn1(j)=yn(j)
          dum2=futp1(i,yn1,a,b,tauwx,tauwy)
          yn1(j)=yn(j)-dyn1
!         yn1(j)=yn(j)
          dum1=futp1(i,yn1,a,b,tauwx,tauwy)
          g(i,j)=(dum2-dum1)/(2d0*dyn1)
!         g(i,j)=(dum2-dum1)/(1d0*dyn1)
        enddo
      enddo

!     write(*,100)((g(i,j),j=1,nn),i=1,nn)
100   format(15(1x,f10.5))
      call solg(nn,15,g,bb,dyn)
      do i=1,nn
        yn1(i)=yn(i)-dyn(i)
!       write(*,*)i,yn(i),yn1(i),dyn(i)
      enddo

!     compute and check residual
      dum=0.
      do i=1,nn
        dum=dum+abs(yn1(i)-yn(i))
      end do
!     write(*,*)dum
!     pause
!     continue

      if(dum.gt.1e-06)goto 10

      return
      end

      double precision function futp1(i,y,a,b,tauwx,tauwy)

      implicit none
      include 'paramt.h'

      integer i,ii,j
      double precision y,a,b,tauwx,tauwy,wrr,wn
      dimension y(15),a(10,15),b(10)
      double precision ro2,ro3,p2,p3,u2,u3,gam,delta
      double precision v2,v3,un2,un3,ut2,ut3,e1,e2
      double precision nwx,nwy
      double precision taux23,tauy23,nx23,ny23,sx23
      double precision taux14,tauy14,nx14,ny14,sx14
      double precision sx23corr
      logical flag1

!     assign constants and variables
      gam=ga
      delta=(gam-1.0)/2.

      ro2=y(5)
      p2 =y(6)
      u2 =y(7)
      v2 =y(8)
      ro3=y(9)
      p3 =y(10)
      u3 =y(11)
      v3 =y(12)
      sx23=y(14)
      wrr=y(15)

!     compute normal and tangent versor
      taux23=cos(sx23)
      tauy23=sin(sx23)
      nx23=-tauy23
      ny23=taux23
!     write(*,*)'nx23,ny23',nx23,ny23

!     compute the normal and tangential velocity components
      un2=u2*nx23+v2*ny23
      un3=u3*nx23+v3*ny23
      ut2=u2*taux23+v2*tauy23
      ut3=u3*taux23+v3*tauy23

!     compute the normal velocity component of the reflection point
      wn=wrr*(tauwx*nx23+tauwy*ny23)

      futp1=0.d0
      if(i.eq.1)then
        futp1=ro2*un2-ro3*un3-wn*(ro2-ro3)
      elseif(i.eq.2)then
        futp1=p2+ro2*(un2-wn)**2-p3-ro3*(un3-wn)**2
      elseif(i.eq.3)then
        futp1=ut2-ut3
      elseif(i.eq.4)then
       futp1=gam/(gam-1.0)*p2/ro2+0.5*(un2-wn)**2
     &      -gam/(gam-1.0)*p3/ro3-0.5*(un3-wn)**2

      elseif(i.eq.5)then
        nwx=-tauwy
        nwy= tauwx
        futp1=u3*nwx+v3*nwy

      elseif(i.ge.6)then
        futp1=0.d0
        ii=i-(1*4+1)
        do j=1,15
          futp1=futp1+a(ii,j)*y(j)
        enddo
        futp1=futp1-b(ii)
      endif

      return
      end
