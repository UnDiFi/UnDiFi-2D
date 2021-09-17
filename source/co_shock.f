! Compute a shock point

      subroutine  co_shock(x1,x2,wshk,R14)

!     x1(1) and x2(1) upstream and downstream density
!     x1(2) and x2(2) upstream and downstream pressure
!     x1(3) and x2(3) upstream and downstream normal velocity
!     wshk shock velocity

      implicit none
      include 'paramt.h'

      double precision x1,x2,wshk,R14
      dimension x1(4),x2(4)
      integer I,J,K,NN
      double precision rov,rom,pv,pm,uv,um,w,gam,delta,R1,R2,a,dyn
      double precision f,g,yn,yn1,g1,dum,dum1,dum2,dyn1,help,Mm,b
      dimension yn(4),yn1(4),g(4,4),g1(4,4),dyn(4),b(4)
      common /shck/ rom,pm,um,gam,delta,R2

! constants assign
      gam=GA
      delta=(gam-1.0)/2.

! determination of the upstream and downstream state
!     if(x1(2).gt.x2(2))then

! downstream state 1
          rov=x1(1)
          pv =x1(2)
          uv =x1(3)

! upstream state 2
          rom=x2(1)
          pm =x2(2)
          um =x2(3)
!     else
!         rov=x2(1)
!         pv =x2(2)
!         uv =x2(3)

!         rom=x1(1)
!         pm =x1(2)
!         um =x1(3)
!     endif

! calculation of invariants
!     R2 = sqrt(gam*pv/rov)+delta*uv
      R2=R14
!     R1 = sqrt(gam*pm/rom)-delta*um

! initialization of downstream state and shock velocity
! with isentropic flow hypotesis
!     uv = (R2-R1)/2./delta
!     a  = (R2+R1)/2.
!     w  = (um+sqrt(gam*pm/rom)+uv+a)/2.
!     dum=pm/rom**(gam)
!     rov=(dum*gam/a**2)**(1./(1.-gam))
!     pv = a**2*rov/gam

!     help = dsign(1.d0,w)

! initialization of downstream state and shock velocity
! with steady shock hypotesis
! N.B w must be an arbitrarily small value but <> 0
! otherwise the jacobian calculation fails
!     Mm=abs(um)/sqrt(gam*pm/rom)
!     pv=pm*(1.d0+2.d0*gam/(gam+1.d0)*(Mm**2-1.d0))
!     rov=rom*(gam+1.d0)*Mm**2/(2.d0+(gam-1.d0)*Mm**2)
!     uv=um*rom/rov
!     w=-0.001
!     if (R14.eq.0.0d0)then
!        R14=sqrt(gam*pv/rov)+delta*uv
!        R2=R14
!     endif

! compute the downstream state and shock velocity with Newton-Raphson method
! initialize the vector of unknowns
      yn1(1) = rov
      yn1(2) = pv
      yn1(3) = uv
      yn1(4) = w

10    do i=1,4
        yn(i)=yn1(i)
      enddo

! compute jacobian
      do i=1,4
        b(i)=f(i,yn)
        do j=1,4
          do k=1,4
           yn1(k)=yn(k)
          enddo
!         dum1=f(i,yn1)
          dyn1=abs(yn1(j))*.001
          if(dyn1.le.1.0d-7)dyn1=1.0d-7

          yn1(j)=yn(j)-dyn1
          dum1=f(i,yn1)

          yn1(j)=yn(j)+dyn1
          dum2=f(i,yn1)
          g(i,j)=(dum2-dum1)/(2.0d0*dyn1)
!         g(i,j)=(dum2-dum1)/(dyn1)

        enddo
      enddo

! invert jacobian matrix
!     call invmat(g,g1,4)

! compute solution at step n+1
!     do i=1,4
!       dyn(i)=0.
!       do j=1,4
!         dyn(i)=dyn(i)+g1(i,j)*f(j,yn)
!       end do
!       yn1(i)=yn(i)-dyn(i)
!     enddo

      nn=4
      call solg(nn,nn,g,b,dyn)
      do i=1,4
       yn1(i)=yn(i)-0.2d0*dyn(i)
!      yn1(i)=yn(i)-dyn(i)
      enddo

! compute and check residual
      dum=0.
      do i=1,4
       dum=dum+abs(yn1(i)-yn(i))
      end do
      if(dum.gt.1d-07)goto 10

      wshk=yn1(4)
!     wshk=help*ABS(yn1(4))
!     if(x1(2).gt.x2(2))then
          x1(1)=yn1(1)
          x1(2)=yn1(2)
          x1(3)=yn1(3)
!     else
!         x2(1)=yn1(1)
!         x2(2)=yn1(2)
!         x2(3)=yn1(3)
!     endif

      return
      end

! ************************************
      double precision function f(i,y)
      integer i
      double precision y
      dimension y(4)
      double precision rov,rom,pv,pm,uv,um,w,gam,delta,R2
      common /shck/ rom,pm,um,gam,delta,R2

      rov=y(1)
      pv =y(2)
      uv =y(3)
      w  =y(4)

      f=0
      if(i.eq.1)then
        f=rov*(uv-w)-rom*(um-w)
!       f=1.0d0-rom*(um-w)/(rov*(uv-w))
      elseif(i.eq.2)then
        f=pv+rov*(uv-w)**2-pm-rom*(um-w)**2
!       f=1.0+rov*(uv-w)**2/pv-pm/pv-rom*(um-w)**2/pv
      elseif(i.eq.3)then
        f=gam/(gam-1.0)*pv/rov+0.5*(uv-w)**2
     &   -gam/(gam-1.0)*pm/rom-0.5*(um-w)**2
!       f=gam/(gam-1.0)*pv*rom+0.5*(uv-w)**2*rov*rom
!    &   -gam/(gam-1.0)*pm*rov-0.5*(um-w)**2*rov*rom

      elseif(i.eq.4)then
        f=sqrt(gam*pv/rov)+delta*uv-R2
!       f=sqrt(gam*pv)+(delta*uv-R2)*sqrt(rov)
      endif
      return
      end

      subroutine invmat(a,b,r)

      implicit real*8 (a-h,o-z)
!     subroutine for matrix inversion
!     a     : matrix r*r to be inverted (real)
!     b     : matrix r*r inverted
!     r     : dimension (max 20) (integer)
      integer r,j,i,k,l
      real*8 a(4,4),b(4,4),c(4,4)

      do 2 i=1,r
        do 3 j=1,r
          b(i,j)=0.d+00
          c(i,j)=a(i,j)
3       continue
2     continue
      do 1 i=1,r
        b(i,i)=1.d+00
1     continue
      do 10 j=1,r
        do 20 i=j,r
          if(a(i,j).ne.0.) goto 210
20      continue
        do 21 i=1,r
          if(a(j,i).ne.0.) goto 211
21      continue
        goto 10
211     write(*,*)'singular matrix'
        return
210     do 30 k=1,r
         s=a(j,k)
         a(j,k)=a(i,k)
         a(i,k)=s
         s=b(j,k)
         b(j,k)=b(i,k)
         b(i,k)=s
30      continue
        t=1/a(j,j)
        do 40 k=1,r
          a(j,k)=t*a(j,k)
          b(j,k)=t*b(j,k)
40      continue
        do 50 l=1,r
          if(l.eq.j) goto 50
          t=-a(l,j)
          do 60 k=1,r
            a(l,k)=a(l,k)+t*a(j,k)
            b(l,k)=b(l,k)+t*b(j,k)
60        continue
50      continue
10    continue
      do 110 i=1,r
        do 120 j=1,r
          a(i,j)=c(i,j)
120     continue
110   continue
100   return
      end

      subroutine  co_dc(x1,x2,wdc)

!     x1(1) and x2(1) upstream and downstream density
!     x1(2) and x2(2) upstream and downstream pressure
!     x1(3) and x2(3) upstream and downstream normal velocity
!     wdc  contact discontinuity velocity

      implicit none
      include 'paramt.h'

      double precision x1,x2,wdc
      dimension x1(4),x2(4)

      integer I,J,K,NN
      double precision ro1,ro2,p1,p2,u1,u2,w,gam,delta,R1,R2,S1,S2,dyn
      double precision fdc,g,yn,yn1,g1,dum,dum1,dum2,dyn1,help,Mm,b
      dimension yn(7),yn1(7),g(7,7),g1(7,7),dyn(7),b(7)

      common /dc/ gam,delta,R1,R2,S1,S2

      ro1=x1(1)
      p1 =x1(2)
      u1 =x1(3)
      ro2=x2(1)
      p2 =x2(2)
      u2 =x2(3)
      wdc=wdc

! assign constants
      gam=GA
      delta=(gam-1.0)/2.

! compute invariants
      R1 = sqrt(gam*p1/ro1)+delta*u1
      R2 = sqrt(gam*p2/ro2)-delta*u2
      S1 = p1/ro1**gam
      S2 = p2/ro2**gam

! compute the downstream state and shock velocity with Newton-Raphson method

! initialize the vector of unknowns
      yn1(1) = ro1
      yn1(2) = p1
      yn1(3) = u1
      yn1(4) = ro2
      yn1(5) = p2
      yn1(6) = u2
      yn1(7) = wdc

10    do i=1,7
        yn(i)=yn1(i)
      enddo

! compute jacobian
      do i=1,7
        b(i)=fdc(i,yn)
        do j=1,7
          do k=1,7
           yn1(k)=yn(k)
          enddo
!         dum1=fdc(i,yn1)
          dyn1=abs(yn1(j))*.01
          if(dyn1.le.1.0d-7)dyn1=1.0d-7

          yn1(j)=yn(j)-dyn1
          dum1=fdc(i,yn1)

          yn1(j)=yn(j)+dyn1
          dum2=fdc(i,yn1)
          g(i,j)=(dum2-dum1)/(2.0d0*dyn1)
!         g(i,j)=(dum2-dum1)/(dyn1)

        enddo
      enddo

      nn=7
      call solg(nn,nn,g,b,dyn)
!     write(*,*)'******'
      do i=1,7
!      yn1(i)=yn(i)-0.5d0*dyn(i)
       yn1(i)=yn(i)-dyn(i)
!      write(*,*)i,yn1(i)
      enddo

! compute and check residual
      dum=0.
      do i=1,7
       dum=dum+abs(yn1(i)-yn(i))
      end do
!     write(*,*)dum
      if(dum.gt.1d-10)goto 10

      wdc=yn1(7)
      x1(1)=yn1(1)
      x1(2)=yn1(2)
      x1(3)=yn1(3)
      x2(1)=yn1(4)
      x2(2)=yn1(5)
      x2(3)=yn1(6)

      return
      end

      double precision function fdc(i,y)
      integer i
      double precision y
      dimension y(7)
      double precision ro1,ro2,p1,p2,u1,u2,w,gam,delta
      double precision R1,R2,S1,S2
      common /dc/ gam,delta,R1,R2,S1,S2

      ro1=y(1)
      p1 =y(2)
      u1 =y(3)
      ro2=y(4)
      p2 =y(5)
      u2 =y(6)
      w  =y(7)

      fdc=0.d+0
      if(i.eq.1)then
        fdc=sqrt(gam*p1/ro1)+delta*u1-R1
      elseif(i.eq.2)then
        fdc=p1/ro1**gam-S1
      elseif(i.eq.3)then
        fdc=sqrt(gam*p2/ro2)-delta*u2-R2
      elseif(i.eq.4)then
        fdc=p2/ro2**gam-S2
      elseif(i.eq.5)then
        fdc=p1-p2
      elseif(i.eq.6)then
        fdc=u1-u2
      elseif(i.eq.7)then
        fdc=w-u1
      endif

      return
      end
