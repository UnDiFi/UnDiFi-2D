      subroutine calc_vel(npoin,varray,dt,xy,wsh,i,ispredictor,
     +                    nowtime,testcase)

!     This subroutine moves the grid-points according to the motion of the piston
!     x=a*t^3*(1-x/l); where "a" depends on the acceleration of the piston 
!     itself and "l" is the total grid length along the x-coordinate

      implicit none
      include 'paramt.h'
      
      integer npoin,i
      double precision xy(ndim,npoin),varray(ndim,npoin),l0
      double precision wsh(2,npshmax*nshmax)
      double precision xnew(ndim,npoin),xold(ndim,npoin)
      integer ipoin,ilist,ish
      double precision dt,t,l,vp,x0, nowtime
      character c5*5, xyfile*11, velfile*12,ispredictor
      character(len=20) :: testcase

!     distinguish unsteady testcases      
      if (testcase == 'ShockExpansion') then
        vp=0.734611368103595d+0
        l0=5.0d0
      else
        vp=0.0d0
        l0=2.0d0
      end if

!     predictor step
      if (ispredictor .eq. 'y') then
        t=2.0d+0*dt*(i-1)+dt
        open(46,file="timesteps.dat",status="old")
        write(46,*) dt
        close(46)
      else 
!     corrector step
        t=dt*i
        open(46,file="timesteps.dat",status="old")
        write(46,*) dt
        close(46)
      end if

!     actual simulation time        
      nowtime = t
         
!     piston 	
      x0=xy(1,1)

!     distance of the piston from the right wall
      l=l0-x0
      xold=xy
      xnew=xy

      write(c5,fmt="(i5.5)")i
      xyfile="co"//c5//".dat"
      velfile="vel"//c5//".dat"
      open(7,file='log/calc_vel.log')

      ! creation of velocity file for neo
      open(10, file='./NEO_data/input/vel.dat', status='unknown')  
      write(10,'(1a35)')'        u                         v'
      if(mod(i,ibak).eq.0)then
        open(8,file='vel/'//xyfile)
        open(9,file='vel/'//velfile)
        write(8,*)'%        ','x                         ','y'
        write(9,*)'%        ','u                         ','v'
      endif 
      write(7,*)'dt=',dt
      write(7,*)'t=',t
      write(7,*)'l=',l
      write(7,*)'npoin=',npoin
      write(7,*)'xold=',x0

!     implementation of the grid velocity
      ilist=npoin+2*nshmax*npshmax
      write(7,*)'ilist=',ilist

      do ipoin=1,npoin
!      undeformed initial condition
       if(xnew(1,ipoin).eq.x0)then
!         piston law 
          xnew(1,ipoin)=vp*(1-1/exp(t))
!         no motion in y
          xnew(2,ipoin)=xnew(2,ipoin)
       else
        xnew(1,ipoin)=vp*(1-1/exp(t))+
     +               (xold(1,ipoin)-x0)*(l0-vp*(1-1/exp(t)))/l
        xnew(2,ipoin)=xnew(2,ipoin)
        varray(1,ipoin)=varray(1,ipoin)*(1-(xold(1,ipoin)-x0)/l)
       endif
        varray(1,ipoin)=(xnew(1,ipoin)-xold(1,ipoin))/dt
        varray(2,ipoin)=0d+0  
        write(10,*) varray(1,ipoin),varray(2,ipoin)
       if(mod(i,ibak).eq.0)then 
         write(8,*) xnew(1,ipoin),xnew(2,ipoin)
         write(9,*) varray(1,ipoin),varray(2,ipoin)
       endif
      end do

      do ipoin=npoin+1,npoin+nshmax*npshmax
        ish=ipoin-npoin
        varray(1,ipoin)=wsh(1,ish)
        varray(2,ipoin)=wsh(2,ish)
        write(10,*) varray(1,ipoin),varray(2,ipoin)
       if(mod(i,ibak).eq.0)then
         write(9,*) varray(1,ipoin),varray(2,ipoin)
       endif
      end do

      do ipoin=npoin+nshmax*npshmax+1,npoin+2*nshmax*npshmax
        ish=ipoin-nshmax*npshmax-npoin
        varray(1,ipoin)=wsh(1,ish)
        varray(2,ipoin)=wsh(2,ish)
        write(10,*) varray(1,ipoin),varray(2,ipoin)
       if(mod(i,ibak).eq.0)then
         write(9,*) varray(1,ipoin),varray(2,ipoin)
       endif
      end do

      close(10)
      close(7)

      if(mod(i,ibak).eq.0)then
        close(8)
        close(9)
      endif

      return
      end
