!**********************************************************************************
!  PURPOSE:  To compute the adimensional heat source due to joule effect.
!
!  $Id: ohmsource.f90,v 1.4 2014/04/10 09:14:21 tesistim Exp $
!
!*****************************************************************************************

       double precision function ohmheat(zroe,pstar,esqr,NDIM,NOFVAR)

       implicit none

       ! Variables   

       include 'paramt.h' 
       include 'plasma.h'
       include 'chem.h'
       include 'paramchem.h'
       include 'constants.h'
       include 'commonv.inc'
       include 'commonchem.inc'
       include 'streamplasma.com'
       include 'electric.com'

      integer i,j,isp,ir,NDIM,NOFVAR
      double precision zroe(MAXNOFVAR),pstar,press,sqrtr,rho_R,rho_Rstar
      double precision T,rho,rhostar(NSP),rhos(NSP),mols(NSP),chis(NSP)
      double precision Tstar,esqr,sigma,ohmheat2
      double precision elnorm,elnormsqrt,nmol
      double precision ddot
      logical verbose
!      parameter(verbose=.TRUE.)
     parameter(verbose=.FALSE.)
    
!      Z_rho (sqrt(rho))  
       sqrtr = ZERO
       do isp = 1,NSP
           sqrtr = sqrtr + zroe(isp) 
       enddo   
  
!     Chemical species densities (kg/m^3)  
      do i = 1,NSP
          rhostar(i) = sqrtr * zroe(i)
          rhos(i) = rhostar(i) * RREFP   !(kg/m^3)
      enddo

!      Pressure
       press = pstar * PREFP                    !(Pa)
       
!      \rho* R* = \sum(rhoi* Ri*) 
       rho_Rstar = ZERO 
       do i = 1,NSP
           rho_Rstar = rho_Rstar + rhostar(i) * RSTARP(i)
       enddo 
!
!      Temperature
       Tstar = pstar/rho_Rstar 
       T = Tstar * TREFP                           !(K)
!
!      Molar density (mol/cm^3)     
       nmol=ZERO
       do i = 1 , NSP
           mols(i) = rhos(i)/Ms(i)     !(kmol/m^3)
           nmol = nmol + mols(i)
       enddo
!      kmol/m^3 -> mol/cm^3
!       mols = mols * 1.0d-3 
!
!      Molar fractions         
       do i = 1,NSP
          chis(i) = mols(i)/nmol
          if(verbose)write(6,*)'chi(',i,')=,',chis(i) 
       enddo 
       
!     Scalar product of electric field
      if(verbose)write(6,*)'Elect field^2 ADIM=',esqr

      elnormsqrt = esqr*PHIREF*PHIREF/(LREFP*LREFP)  !Dimensional electric field squared norm (V^2/m^2)              
      if(verbose)write(6,*)'Elect field^2 DIM=',elnormsqrt

!     Electric field norm
      elnorm = sqrt(elnormsqrt)   ! (V/m)                    
      if(verbose)write(6,*)'Elect field DIM=',elnorm

!     Electrical conductivity 
      if(verbose)then
       write(6,*)'Temperature=',T
      endif
      call electcond(T,chis,elnorm,sigma)  !(1/(Ohm m))      
      if(verbose)then
       write(6,*)'sigma DIM (OHM/M)=',sigma
       write(6,*)'Temperature=',T
      endif
  
!pepe elnorm = PHIREF/1.d-1
!pepe elnormsqrt = elnorm*elnorm
!
!     Ohmic source term
      Ohmheat = sigma*elnormsqrt           !(J/(m^3 s)
      ohmheat2 = sigma*esqr
     
      if(verbose)write(6,*)'OHMHEAT DIM=',Ohmheat
!
!     Non-dimensional ohmic source term
      Ohmheat = Ohmheat * LREFP / (RREFP*UREFP**3) 
      ohmheat2 = ohmheat2 * PHIREF**2/(LREFP*RREFP*UREFP**3)

      if(verbose)then
         write(6,*)'OHMHEAT ADIM=',Ohmheat
         write(6,*)'OHMHEAT2 ADIM=',ohmheat2
         pause    
      endif
  
      return
      end
