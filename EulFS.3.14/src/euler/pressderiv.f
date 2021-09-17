      subroutine pressderiv(Z,chis,kappa)

!*********************************************************************
!
!     PURPOSE:  To compute the pressure derivatives chis and k
!      
!     Input: Z,Temp 
!     Ooutput: chis, k
!
!*********************************************************************

      INCLUDE'paramt.h'
      INCLUDE'paramchem.h'  
      INCLUDE'constants.h'
      INCLUDE'chem.h'
      INCLUDE'plasma.h'
      INCLUDE'streamplasma.com'
      INCLUDE'commonv.inc'
 
      double precision Z(MAXNOFEQN),chis(NSP),kappa
      double precision en(NSP),Cv(NSP),Temp,pres
      double precision ZrR,ZrCv,ZrHf,Zrho,gams(nsp),KINETIC     
      integer isp,id      
           
      ZrR = ZERO
      ZrCv = ZERO
      ZrHf = ZERO
      Zrho = ZERO 
      do isp = 1,NSP
!         Cv(isp) = GOGM1 * RSTARP(isp) !ideal mixture of caloric perfect gases (monoatomic or diatomic) 
         gams(isp) = GAM  !ideal mixture of caloric perfect gases (monoatomic or diatomic)
         Cv(isp) = ONE/(gams(isp)-1)*RSTARP(isp)                 
         ZrR = ZrR + Z(isp)*RSTARP(isp) 
         ZrCv = ZrCv + Z(isp)*Cv(isp)
         ZHf = ZHf + Z(isp)*HF0(isp)/HREFP
         Zrho = Zrho + Z(isp) 
      enddo    
      
      KINETIC = HALF * (Z(IX)*Z(IX) + Z(IY)*Z(IY))
      if(NDIM.eq.3)then
         KINETIC = KINETIC + HALF*Z(IZ)*Z(IZ)
      endif
 
!     Temperature
      Temp = (Z(IE) - KINETIC/Zrho - ZHf)/(ZrR + ZrCv) !ideal mixture of caloric perfect gases
!      write(6,*)'TEMP=',Temp
       
!     Pressure  
      pres = Zrho*ZrR*Temp
!      write(6,*)'pres=',pres
        
!     Pressure derivatives      
      kappa = ZrR/ZrCv   

      do isp = 1,NSP
         en(isp) = Cv(isp)*Temp + HF0(isp)/HREFP    !caloric perfect gas 
         chis(isp) = RSTARP(isp)*Temp-kappa*en(isp)      
      enddo                                                            

      return  
      end
