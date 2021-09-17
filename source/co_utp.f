! Compute an unsteady triple point (TP)

      subroutine co_utp(y,r14,dxr14,dyr14,r23,unsh1,yn1,flag1,ifail)

      implicit none
      include 'paramt.h'

      integer i,j,k,nn,icont
      double precision gam,delta,a,b,bb
      double precision r14,dxr14,dyr14,r23,unsh1
      double precision y,futp,g,yn,yn1,g1,dum,dum1,dum2,dyn1,dumold
      double precision theta,dyn,sx14,taux14,tauy14,nx14,ny14
      double precision store13,store14,store15,store16
      logical flag1,ifail
      dimension y(20),yn(20),yn1(20),g(20,20),g1(20,20),a(11,20),b(11)
      dimension bb(20),dyn(20)

!     open log file
      open(8,file='log/co_utp.log')

!     assign constants
      gam=ga
      delta=(gam-1.0)/2.

      nn=20

!     memorize state 4
      store13=y(13)
      store14=y(14)
      store15=y(15)
      store16=y(16)

! build matrix a and vector b of the additional conditions
! which varies depending on the case

! calculation of oblique shock being the upstream state and deviation known

       do i=1,nn-(2*4+1)
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

! value sx12 known
       a(9,17)=1.0d+0
       b(9)=y(17)

! p3=p4
!      a(10,10)=1.0d+0
!      a(10,14)=-1.0d+0
!      b(10)=0.d+0

! theta3=theta4 --> u3v4-u4v3=0
!      a(11,16)=y(11)
!      a(11,12)=-y(15)
!      b(11)=0.0d+0

!     do i=1,nn
!      write(*,*)i,futp(i,y,a,b,r14)
!     enddo

! calcuate the downstream state and shock velocity
! with the newton-raphson method

! initialization of the vector of unknowns
      do i=1,nn
       yn1(i) = y(i)
      enddo

! initialization of the iteration counter
      icont=0

10    do i=1,nn
        yn(i)=yn1(i)
        bb(i)=futp(i,yn1,a,b,r14,dxr14,dyr14,r23,unsh1,flag1)
!       write(*,*)i,bb(i)
      enddo

! jacobian calculation
      do i=1,nn
!       bb(i)=futp(i,yn,a,b,r14,r23,unsh1,flag1)
        do j=1,nn
          do k=1,nn
           yn1(k)=yn(k)
          enddo
!         dum1=futp(i,yn1,a,b,r14,r23,unsh1,flag1)
          dyn1=abs(yn1(j))*.001
          if(dyn1.lt.1.0d-7)dyn1=1.0d-7
          yn1(j)=yn(j)+dyn1
          dum2=futp(i,yn1,a,b,r14,dxr14,dyr14,r23,unsh1,flag1)
          yn1(j)=yn(j)-dyn1
          dum1=futp(i,yn1,a,b,r14,dxr14,dyr14,r23,unsh1,flag1)
          g(i,j)=(dum2-dum1)/(2*dyn1)
!         g(i,j)=(dum2-dum1)/(dyn1)
        enddo
      enddo

!     write(*,100)((g(i,j),j=1,20),i=1,20)
100   format(20(1x,f10.5))
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
      icont=icont+1
      write(8,*)'conv--->',dum, icont
!     if(dum.gt.10.)then
!        write(*,*)'change r23 --> r14'
!        flag1=.false.
!      endif

! update coefficients matrix a with calculated values
!     a(11,16)=yn1(11)
!     a(11,12)=-yn1(15)

! recalculate r14 using the new vector n14
       sx14=yn1(19)
       taux14=cos(sx14)
       tauy14=sin(sx14)
       nx14=-tauy14
       ny14=taux14

!      r14=sqrt(ga*y(14)/y(13))+
!    &        gm1*0.5d0*(y(15)*nx14+y(16)*ny14)
!      r14=sqrt(ga*store14/store13)+
!    &        gm1*0.5d0*(store15*nx14+store16*ny14)

!     do i=1,nn
!      write(*,*)i,futp(i,yn1,a,b)
!     enddo

      if(icont.eq.1)then
        dumold=dum
        goto 10
      endif
      if(dum.gt.dumold.and.icont.gt.1)then
        ifail=.true.
        return
      endif
      if(dum.gt.1e-07)then
         dumold=dum
         goto 10
      endif

      write(8,*)'initial and final state'
      do i=1,nn
       write(8,*)i,y(i),yn(i)
      enddo
      close(8)
      return
      end

      double precision function futp(i,y,a,b,r14,dxr14,dyr14,r23,unsh1,
     &                               flag1)

      implicit none
      include 'paramt.h'

      integer i,ii,j
      double precision y,a,b,r14,r23,wn,wt, wsh,unsh1
      dimension y(20),a(11,20),b(11)
      double precision ro1,ro2,p1,p2,u1,u2,gam,delta
      double precision v1,v2,un1,un2,ut1,ut2,e1,e2
      double precision taux12,tauy12,nx12,ny12,sx12
      double precision taux23,tauy23,nx23,ny23,sx23
      double precision taux14,tauy14,nx14,ny14,sx14
      double precision sx23corr,dxr14,dyr14
      logical flag1

!     assign constants and variables
      gam=ga
      delta=(gam-1.0)/2.

      if(i.ge.1.and.i.le.4)then
       ro1=y(5)
       p1 =y(6)
       u1 =y(7)
       v1 =y(8)
       ro2=y(9)
       p2 =y(10)
       u2 =y(11)
       v2 =y(12)
       sx12=y(17)
       sx23=y(18)
       wsh=y(20)

! calculation of the tangent and normal versor to shock 12 (incident shock)
       taux12=cos(sx12)
       tauy12=sin(sx12)
       nx12=-tauy12
       ny12=taux12

! calculation of the tangent and normal versor to shock 23
       taux23=cos(sx23)
       tauy23=sin(sx23)
       nx23=-tauy23
       ny23=taux23

! calculation of the normal and tangential velocity components
       un1=u1*nx23+v1*ny23
       un2=u2*nx23+v2*ny23
       ut1=u1*taux23+v1*tauy23
       ut2=u2*taux23+v2*tauy23

! calculation of the normal velocity component to shock 23 in the triple point
       wn=wsh*(taux12*nx23+tauy12*ny23)
!    &  +unsh1*(nx12  *nx23+ny12  *ny23)

       wt=wsh*(taux12*taux23+tauy12*tauy23)+
     &  unsh1*(nx12  *taux23+ny12  *tauy23)

!       sx23corr=atan(wt/un1)
        sx23corr=0.

! modification
!      u1=u1-wt*taux23
!      v1=v1-wt*tauy23
!      u2=u2-wt*taux23
!      v2=v2-wt*tauy23

!      un1=u1*nx23+v1*ny23
!      un2=u2*nx23+v2*ny23
!      ut1=u1*taux23+v1*tauy23
!      ut2=u2*taux23+v2*tauy23

       taux23=cos(sx23-sx23corr)
       tauy23=sin(sx23-sx23corr)
       nx23=-tauy23
       ny23=taux23

! calculation of the normal and tangential velocity components
       un1=u1*nx23+v1*ny23
       un2=u2*nx23+v2*ny23
       ut1=u1*taux23+v1*tauy23
       ut2=u2*taux23+v2*tauy23

! calculation of the normal velocity component of the triple point
       wn=wsh*(taux12*nx23+tauy12*ny23)
!    &  +unsh1*(nx12  *nx23+ny12  *ny23)

       futp=0.d0
       if(i.eq.1)then
         futp=ro1*un1-ro2*un2-wn*(ro1-ro2)
       elseif(i.eq.2)then
!        futp=p1+ro1*(un1)**2-p2-ro2*(un2)**2-wn*(ro1*un1-ro2*un2)
         futp=p1+ro1*(un1-wn)**2-p2-ro2*(un2-wn)**2

       elseif(i.eq.3)then
!        futp=ro1*un1*ut1-ro2*un2*ut2-wn*(ro1*ut1-ro2*ut2)
!        futp=ut1-ut2-wn*(ut1-ut2)
         futp=ut1-ut2
       elseif(i.eq.4)then
!        e1=ro1*(1.0/(gam-1.0)*p1/ro1+(u1**2+v1**2)*0.5)
!        e2=ro2*(1.0/(gam-1.0)*p2/ro2+(u2**2+v2**2)*0.5)
!        futp=gam/(gam-1.0)*p1*un1+0.5*ro1*un1*(u1**2+v1**2)
!    +    -gam/(gam-1.0)*p2*un2-0.5*ro2*un2*(u2**2+v2**2)
!    +    -wn*(e1-e2)
        futp=gam/(gam-1.0)*p1/ro1+0.5*(un1-wn)**2
     &      -gam/(gam-1.0)*p2/ro2-0.5*(un2-wn)**2

       endif
      elseif(i.ge.5.and.i.le.8+1)then
       ro1=y(1)
       p1 =y(2)
       u1 =y(3)
       v1 =y(4)
       ro2=y(13)
       p2 =y(14)
       u2 =y(15)
       v2 =y(16)
       sx12=y(17)
       sx14=y(19)
       wsh=y(20)

! calculation of the tangent and normal versor to shock 12 (incident shock)
       taux12=cos(sx12)
       tauy12=sin(sx12)
       nx12=-tauy12
       ny12=taux12

! calculation of the tangent and normal versor to shock 14
       taux14=cos(sx14)
       tauy14=sin(sx14)
       nx14=-tauy14
       ny14=taux14

! calculation of the normal and tangential velocity components
       un1=u1*nx14+v1*ny14
       un2=u2*nx14+v2*ny14
       ut1=u1*taux14+v1*tauy14
       ut2=u2*taux14+v2*tauy14

! calculation of the normal velocity component of the triple point
       wn=wsh*(taux12*nx14+tauy12*ny14)
!    &  +unsh1*(nx12  *nx14+ny12  *ny14)
      wt=wsh*(taux12*taux14+tauy12*tauy14)+
     &  unsh1*(nx12  *taux14+ny12  *tauy14)

! modification
!      u1=u1+wt*taux14
!      v1=v1+wt*tauy14
!      u2=u2+wt*taux14
!      v2=v2+wt*tauy14

!      un1=u1*nx14+v1*ny14
!      un2=u2*nx14+v2*ny14
!      ut1=u1*taux14+v1*tauy14
!      ut2=u2*taux14+v2*tauy14

       futp=0.d0
       if(i.eq.5)then
         futp=ro1*un1-ro2*un2-wn*(ro1-ro2)
       elseif(i.eq.6)then
!        futp=p1+ro1*(un1)**2-p2-ro2*(un2)**2-wn*(ro1*un1-ro2*un2)
         futp=p1+ro1*(un1-wn)**2-p2-ro2*(un2-wn)**2

       elseif(i.eq.7)then
!        futp=ro1*un1*ut1-ro2*un2*ut2-wn*(ro1*ut1-ro2*ut2)
!        futp=ut1-ut2-wn*(ut1-ut2)
         futp=ut1-ut2
       elseif(i.eq.8)then
!        e1=ro1*(1.d0/(gam-1.0)*p1/ro1+(u1**2+v1**2)*0.5)
!        e2=ro2*(1.d0/(gam-1.0)*p2/ro2+(u2**2+v2**2)*0.5)
!        futp=gam/(gam-1.0)*p1*un1+0.5*ro1*un1*(u1**2+v1**2)
!    +    -gam/(gam-1.0)*p2*un2-0.5*ro2*un2*(u2**2+v2**2)
!    +    -wn*(e1-e2)
        futp=gam/(gam-1.0)*p1/ro1+0.5*(un1-wn)**2
     &      -gam/(gam-1.0)*p2/ro2-0.5*(un2-wn)**2
       elseif(i.eq.9)then

!        if (flag1)then

!         ro2=y(9)
!         p2 =y(10)
!         u2 =y(11)
!         v2 =y(12)
!         sx23=y(18)

! calculation of the normal and tangential versor
!         taux23=cos(sx23)
!         tauy23=sin(sx23)
!         nx23=-tauy23
!         ny23=taux23

! calculation of the normal and tangential velocity components
!         un2=u2*nx23+v2*ny23

!         futp=sqrt(gam*p2/ro2)+delta*un2-r23

!        else

          ro2=y(13)
          p2 =y(14)
          u2 =y(15)
          v2 =y(16)
          sx14=y(19)

! calculation of the normal and tangential versor
          taux14=cos(sx14)
          tauy14=sin(sx14)
          nx14=-tauy14
          ny14=taux14

caldo
!         nx14= 0.9
!         ny14=-sqrt(1.00-nx14**2)
          un2=u2*dxr14+v2*dyr14
caldo

! calculation of the normal and tangential velocity components
!         un2=u2*nx14+v2*ny14

          futp=sqrt(gam*p2/ro2)+delta*un2-r14
!        endif

       endif
      elseif(i.eq.10)then
       ro1=y(9)
       p1 =y(10)
       u1 =y(11)
       v1 =y(12)
       ro2=y(13)
       p2 =y(14)
       u2 =y(15)
       v2 =y(16)
       futp=p1-p2
      elseif(i.eq.11)then
       ro1=y(9)
       p1 =y(10)
       u1 =y(11)
       v1 =y(12)
       ro2=y(13)
       p2 =y(14)
       u2 =y(15)
       v2 =y(16)
!      futp=(u1*u2+v1*v2)**2-(u1*u1+v1*v1)*(u2*u2+v2*v2)
       futp=(u1*v2-v1*u2)
      elseif(i.ge.12)then
        futp=0.d0
        ii=i-(2*4+1+1+1)
        do j=1,20
          futp=futp+a(ii,j)*y(j)
        enddo
        futp=futp-b(ii)
      endif

      return
      end
