    subroutine electcond(T,chis,Elnorm,sigma)
!
!   Compute electrical conductivity
!
!   $Id: electcond.f90,v 1.2 2014/04/10 09:15:19 tesistim Exp $
!
    implicit none

    include 'plasma.h'
    include 'paramchem.h'

!   Variables
    integer i,j,ir,fitsigma
    double precision chis(NSP),Elnorm,T
    double precision fs1,fs2,fs4
    double precision sigma

!   Choice a fit for sigma
!   fitsigma = 1 -> Fit of Boyd AIAA article data (only for p=0.013atm)
!   fitsigma = 2 -> Raizer (only for p=1bar)
!   fitsigma = 3 -> Chapman & Cowling
!   fitsigma = 4 -> Fit Colonna & D'Angola data(only p=1bar) 
!   else -> sigma=1

    fitsigma = 5

    if (fitsigma.eq.1)then
       sigma = fs1(T)        !(Ohm/cm)
       sigma = sigma*1.0d2   !(Ohm/m)
    elseif (fitsigma.eq.2)then
       if (T.le.8.d3)then
          T = 8.d3
       endif 
       if(T.ge.1.4d4)then
          T = 1.4d4
       endif
       sigma = 8.3d3*exp(-3.6d5/T)  !(Ohm/m)
    elseif (fitsigma.eq.3)then
       sigma = fs2(chis,T)   !(Ohm/cm) 
       sigma = sigma*1.0d2   !(Ohm/m)
    elseif(fitsigma.eq.4)then              
       sigma = fs4(T)        !(Ohm/m)        
    else
!       write(6,*)'Choice a model for the electrical conductivity'
!       stop
       sigma = 1.d0
    endif
!
   return
   end   
   
!   Analytical fitting function fs1
    double precision function fs1(x)
    implicit none     
    double precision x,m(5)
    !implicit none     
    m(1) = 2.433d-15
    m(2) = -9.35d-11
    m(3) = 1.241d-6
    m(4) = -6.192d-3 
    m(5) = 8.533d0
 
    fs1 = m(1)*x**4 + m(2)*x**3 + m(3)*x**2 + m(4)*x + m(5)
    fs1 = 1.d1**fs1 
!
    return
    end

!   Analytical fitting function fs2
    double precision function fs2(chi,T)
    implicit none
    double precision T,chi(4),Q,alpha
!    implicit none

    Q=5.d-15
!    alpha=chi(3)/(1-chi(3))
    alpha=chi(3)
    fs2 = 3.34d-10 * alpha/Q/sqrt(T)

    return
    end

!    Analytical fitting function fs3
!    double precision function fs3(chi,E)
!    double precision chis(4),E

!    implicit none

!    fs3 =

!    return
!    end
    
!   Analytical fitting function fs4
    double precision function fs4(x)
    implicit none
    double precision x,x1,m1(6),m2(6)

!   INTERVALLO TEMPERATURE (1000 - 5000 K)
    m1(1) =  1.232d-015;
    m1(2) = -1.962d-011;
    m1(3) =  1.202d-007;
    m1(4) =  -0.0003545d0;
    m1(5) =     0.5128d0;
    m1(6) =        -313.d0;

!   INTERVALLO TEMPERATURE (5000 - 13000 K)
    m2(1) =   7.34d-019;
    m2(2) = -3.793d-014;
    m2(3) =  7.856d-010;
    m2(4) = -8.212d-006;
    m2(5) =     0.04402d0;
    m2(6) =      -90.63d0;

    x1 = x
    if(x1.gt.1.3d4)then
        x1=1.3d4;
    endif  
     
    if(x1.lt.1.d3)then
        fs4=0.0d0;
    elseif((x1.ge.1.d3).and.(x1.le.5.d3))then
        fs4=m1(1)*x1**5 + m1(2)*x1**4 + m1(3)*x1**3 + m1(4)*x1**2 + m1(5)*x1 + m1(6)
        fs4 = exp(fs4)
    elseif(x.gt.5.d3)then
        fs4=m2(1)*x1**5 + m2(2)*x1**4 + m2(3)*x1**3 + m2(4)*x1**2 + m2(5)*x1 + m2(6)
        fs4 = exp(fs4)    
    endif
!
    return
    end


