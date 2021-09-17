      subroutine mv_grid(npoin,varray,dt,xy,wsh,i,testcase)

!     Purpose:
!        this subroutine moves the grid-points according to the motion of the piston
!        x=a*t^3*(1-x/l); where "a" depends on the acceleration of the pistion itself and
!        "l" is the total grid length along the x-coordinate

      implicit none
      include 'paramt.h'
      
      integer npoin,i
      double precision xy(ndim,npoin),varray(ndim,npoin),l0
      double precision wsh(2,nshmax*npshmax)
      double precision xold(ndim,npoin)
      integer ipoin,ilist,status
      double precision dt,t,l,vp,x0
      character c5*5, xyfile*11, velfile*12
      character(len=20) :: testcase

      if (testcase == 'ShockExpansion') then
        vp=0.734611368103595d+0 ! piston's speed
        l0=5.0d0                ! domain's length
      else
        vp=0.0d0
        l0=2.0d0
      end if

!     this subroutine is only called in the corrector step
      t=dt*i 
!     piston (constant value)
      x0=xy(1,1)
!     distance of the piston from the right wall
      l=l0-x0
      xold=xy
      
      write(c5,fmt="(i5.5)")i
      xyfile="xy"//c5//".dat"
      open(7,file='log/mv_grid.log')
      if(mod(i,ibak).eq.0)then
        open(8,file='vel/'//xyfile)
        write(8,*)'%        ','x                         ','y'
      endif 
      write(7,*)'dt=',dt
      write(7,*)'t=',t
      write(7,*)'l=',l
      write(7,*)'npoin=',npoin
! check file     
!       open(unit=6,file='log/check_mv_grid.log',status='new',
!     & action='write',iostat=status)
!       openif: if (status == 0) then
!      open was ok. write values.	       
!       write(6,*) x0,l,xold
!       else openif 
!       write(*,1040) status
!1040   format (' ','error opening file: iostat = ', i6)
!       end if openif

!     law of motion of the grid

      ilist=npoin+2*nshmax*npshmax
      write(7,*)'ilist=',ilist

      do ipoin=1,npoin
       if(xy(1,ipoin).eq.x0)then
          xy(1,ipoin)=vp*(1-1/exp(t))

!         xy(1,ipoin)=xy(1,ipoin)	!shock-vortex interaction
          xy(2,ipoin)=xy(2,ipoin)

       else
! 
        xy(1,ipoin)=xy(1,1)+(xold(1,ipoin)-x0)*(l0-xy(1,1))/l
!       xy(1,ipoin)=xy(1,ipoin)
        xy(2,ipoin)=xy(2,ipoin)

!       write(*,*) xold(1,npoin)-x0, (l0-xy(1,1))/l, xy(1,1)
!       write(*,*)(xold(1,npoin)-x0)*(l0-xy(1,1))/l, xy(1,1),xy(1,npoin)
!       write(*,*)(xold(1,ipoin)-x0)*(l0-xy(1,1))/l, xy(1,1),xy(1,ipoin)

       endif
       if(mod(i,ibak).eq.0)then 
         write(8,*) xy(1,ipoin),xy(2,ipoin)
       endif
      end do

      close(7)
      if(mod(i,ibak).eq.0)then
        close(8)
!       close(6)
      endif

      return
      end
