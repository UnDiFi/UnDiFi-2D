! Compute an unsteady quadruple point (QP)

      subroutine co_uqp(y,wqpx,wqpy,yn1)

      implicit none
      include 'paramt.h'

      integer i,j,k,nn
      double precision gam,delta,a,b,bb
      double precision wqpx,wqpy
      double precision y,futp2,g,yn,yn1,g1,dum,dum1,dum2,dyn1
      double precision theta,dyn,sx14,taux14,tauy14,nx14,ny14
      double precision store13,store14,store15,store16
      logical flag1
      dimension y(24),yn(24),yn1(24),g(24,24),g1(24,24),a(14,24),b(14)
      dimension bb(24),dyn(24)

      open(8,file='log/co_uqp.log')

!     assign costans
      gam=ga
      delta=(gam-1.0)/2.

      nn=24

! build matrix a and vector b of the additional conditions
! which varies depending on the case

! calculation of incident shock 1, incident shock 2 being the upstream
! state and deviation known
       do i=1,(3*4+2)
         do j=1,nn
           a(i,j)=0.d+0
         enddo
         b(i)=0.d+0
       enddo

! state 1 known
       a(1,1)=1.0d+0
       a(2,2)=1.0d+0
       a(3,3)=1.0d+0
       a(4,4)=1.0d+0
       b(1)=y(1)
       b(2)=y(2)
       b(3)=y(3)
       b(4)=y(4)

! state 2 known
       a(5,5)=1.0d+0
       a(6,6)=1.0d+0
       a(7,7)=1.0d+0
       a(8,8)=1.0d+0
       b(5)=y(5)
       b(6)=y(6)
       b(7)=y(7)
       b(8)=y(8)

! state 3 known
       a( 9, 9)=1.0d+0
       a(10,10)=1.0d+0
       a(11,11)=1.0d+0
       a(12,12)=1.0d+0
       b( 9)=y(9)
       b(10)=y(10)
       b(11)=y(11)
       b(12)=y(12)

! value sx12 known
       a(13,21)=1.0d+0
       b(13)=y(21)

! value sx13 known
       a(14,23)=1.0d+0
       b(14)=y(23)

! calcuate the downstream state and shock velocity
! with the newton-raphson method

! initialization of the vector of unknowns
      do i=1,nn
       yn1(i) = y(i)
      enddo

10    do i=1,nn
        yn(i)=yn1(i)
        bb(i)=futp2(i,yn1,a,b,wqpx,wqpy)
      enddo

! jacobian calculation
      do i=1,nn
        do j=1,nn
          do k=1,nn
           yn1(k)=yn(k)
          enddo
          dyn1=abs(yn1(j))*.001
          if(dyn1.lt.1.0d-7)dyn1=1.0d-7
          yn1(j)=yn(j)+dyn1
          dum2=futp2(i,yn1,a,b,wqpx,wqpy)
          yn1(j)=yn(j)-dyn1
          dum1=futp2(i,yn1,a,b,wqpx,wqpy)
          g(i,j)=(dum2-dum1)/(2*dyn1)
        enddo
      enddo

!     write(*,100)((g(i,j),j=1,24),i=1,24)
100   format(24(1x,f10.5))
      call solg(nn,nn,g,bb,dyn)
      do i=1,nn
       yn1(i)=yn(i)-0.5*dyn(i)
!      write(*,*)i,yn(i),yn1(i),dyn(i)
      enddo

! calculation and check of residual
      dum=0.
      do i=1,nn
       dum=dum+abs(yn1(i)-yn(i))
      end do
      write(8,*)'conv--->',dum
!     if(dum.gt.10.)then
!        write(*,*)'change r23 --> r14'
!        flag1=.false.
!      endif

      if(dum.gt.1e-11)goto 10

      write(8,*)
      do i=1,nn
        write(8,*)i,y(i),yn1(i)
      enddo

      close(8)

      return
      end

      double precision function futp2(i,y,a,b,wqpx,wqpy)

      implicit none
      include 'paramt.h'

      integer i,ii,j
      double precision y,a,b,wqpx,wqpy,wn,wt, wsh,unsh1
      dimension y(24),a(14,24),b(14)
      double precision ro1,ro2,p1,p2,u1,u2,gam,delta
      double precision v1,v2,un1,un2,ut1,ut2,e1,e2
      double precision taux12,tauy12,nx12,ny12,sx12
      double precision taux23,tauy23,nx23,ny23,sx23
      double precision taux14,tauy14,nx14,ny14,sx14
      double precision sx23corr
      logical flag1

! assign constants and variables
      gam=ga
      delta=(gam-1.0)/2.

      if(i.ge.1.and.i.le.4)then
       ro1=y(5)
       p1 =y(6)
       u1 =y(7)
       v1 =y(8)
       ro2=y(13)
       p2 =y(14)
       u2 =y(15)
       v2 =y(16)
       sx12=y(22)

! calculation of the tangent and normal versor to shock 24 (reflected shock 1)
       taux12=cos(sx12)
       tauy12=sin(sx12)
       nx12=-tauy12
       ny12=taux12

! calculation of the normal and tangential velocity components
       un1=u1*nx12+v1*ny12
       un2=u2*nx12+v2*ny12
       ut1=u1*taux12+v1*tauy12
       ut2=u2*taux12+v2*tauy12

! calculation of the normal velocity component to shock 24 in the quadruple point
       wn=wqpx*nx12+wqpy*ny12

       futp2=0.d0
       if(i.eq.1)then
         futp2=ro1*un1-ro2*un2-wn*(ro1-ro2)
       elseif(i.eq.2)then
         futp2=p1+ro1*(un1-wn)**2-p2-ro2*(un2-wn)**2
       elseif(i.eq.3)then
         futp2=ut1-ut2
       elseif(i.eq.4)then
        futp2=gam/(gam-1.0)*p1/ro1+0.5*(un1-wn)**2
     &       -gam/(gam-1.0)*p2/ro2-0.5*(un2-wn)**2

       endif
      elseif(i.ge.5.and.i.le.8)then
       ro1=y(9)
       p1 =y(10)
       u1 =y(11)
       v1 =y(12)
       ro2=y(17)
       p2 =y(18)
       u2 =y(19)
       v2 =y(20)
       sx12=y(24)

! calculation of the tangent and normal versor to shock 35 (reflected shock 2)
       taux12=cos(sx12)
       tauy12=sin(sx12)
       nx12=-tauy12
       ny12=taux12

! calculation of the normal velocity component to shock 35 in the quadruple point
       wn=wqpx*nx12+wqpy*ny12

! calculation of the normal and tangential velocity components
       un1=u1*nx12+v1*ny12
       un2=u2*nx12+v2*ny12
       ut1=u1*taux12+v1*tauy12
       ut2=u2*taux12+v2*tauy12

       futp2=0.d0
       if(i.eq.5)then
         futp2=ro1*un1-ro2*un2-wn*(ro1-ro2)
       elseif(i.eq.6)then
         futp2=p1+ro1*(un1-wn)**2-p2-ro2*(un2-wn)**2

       elseif(i.eq.7)then
         futp2=ut1-ut2
       elseif(i.eq.8)then
        futp2=gam/(gam-1.0)*p1/ro1+0.5*(un1-wn)**2
     &       -gam/(gam-1.0)*p2/ro2-0.5*(un2-wn)**2
       endif
       elseif(i.eq.9)then
       ro1=y(13)
       p1 =y(14)
       u1 =y(15)
       v1 =y(16)
       ro2=y(17)
       p2 =y(18)
       u2 =y(19)
       v2 =y(20)
       futp2=p1-p2
      elseif(i.eq.10)then
       ro1=y(13)
       p1 =y(14)
       u1 =y(15)
       v1 =y(16)
       ro2=y(17)
       p2 =y(18)
       u2 =y(19)
       v2 =y(20)
       futp2=(u1*u2+v1*v2)**2-(u1*u1+v1*v1)*(u2*u2+v2*v2)
      elseif(i.ge.11)then
        futp2=0.d0
        ii=i-(2*4+1+1)
        do j=1,24
          futp2=futp2+a(ii,j)*y(j)
        enddo
        futp2=futp2-b(ii)
      endif

      return
      end
