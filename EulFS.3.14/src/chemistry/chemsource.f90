!**********************************************************************************
!  PURPOSE:  To compute the adimensional chemical source term for an Argon mixture.
!
!  CHEMICAL SPECIES
!  Ar: Argon ground level
!  Ar*: Argon metastable (4s level)
!  Ar+: Ionized Argon (first ionizzation)
!  e-: Electron
!
!  CHEMICAL REACTIONS
!  Electronic processes (kf e kb by Colonna)
!  p1) Ar + e <-> Ar+ + e + e
!  p2) Ar + e <-> Ar* + e
!  p3) Ar* + e <-> Ar+ + e + e
!
!  Atom-atom processes (kf e kb by Bacri & Vlcek)
!  p1) Ar + Ar <-> Ar + e + Ar+
!  p2) Ar + Ar <-> Ar + Ar*
!  p3) Ar + Ar* <-> Ar + e + Ar+
! 
!   VARIABLES
!
!   Ns,Nreac: species & chemical reactions number 
!   Nrmax: max chemical reactions number     
!   Nc: number of electronic processes (Colonna fit)
!   Nb: number of atom-atom processes (Bacri & Vlcek fit)
!   scr,scp: stoichiometric coefficient (reagents,products)
!   kf,kb: forward and backward rates
!   Keq: chemical equilibrium constant (particles concentration) 
!   rhos(i): density species i
!   rhoc(i): molar concentration species i
!   rhoN(i): particle concentration species i
!   Ss(i): source term species i
! 
!   INPUT
!    
!   zroe: adimensional parameter vector 
!   pstar: adimensional pressure   
!   
!   OUTPUT 
!  
!   Sstar: adimensional source term (Ss* = Ss *(L/rho/u))     
!
!*****************************************************************************************

       subroutine chemsource(zroe,pstar,Sstar)

       implicit none

       ! Variables   

       include 'paramt.h' 
       include 'plasma.h'
       include 'chem.h'
       include 'paramchem.h'
       include 'commonv.inc'
       include 'commonchem.inc'
       include 'streamplasma.com'
       include 'tauchem.com'

      integer i,j,isp,ir
      double precision zroe(MAXNOFVAR),pstar,press,sqrtr,rho_R,rho_Rstar
      double precision T,rho,rhostar(NSP),rhos(NSP),rhoc(NSP),rhoN(NSP)
      double precision k1,k2,kf(NRMAX),kb(NRMAX),Tstar
      double precision Ss(NSP),S,Sstar(NSP),Shelp !,tau(NSP)
      double precision SsSU2(NSP)
      integer fitrate,Stanford
     
!      Z_rho (sqrt(rho))  
       sqrtr = 0d0
       do isp = 1,NSP
           sqrtr = sqrtr + zroe(isp) 
       enddo   
  
!       write(6,*)'NREAC=',NREAC 
!       write(6,*)'CHEM=',CHEM 

!     Chemical species densities (kg/m^3)  
      do i = 1,NSP
          rhostar(i) = sqrtr * zroe(i)
          rhos(i) = rhostar(i) * RREFP   !(kg/m^3)
      enddo

!      Pressure
!      External: PREFP=RREFP*UREFP^2        
!      Internal: PREFP=P10

       press = pstar * PREFP                    !(Pa)
       
!      \rho* R* = \sum(rhoi* Ri*) 
!      External: RSTARP(i) = Tinf/Uinf^2 * R(i)
!      Internal: RSTARP(i) = rho0 * T0 /p0 * R(i)
 
       rho_Rstar = 0d0  
       do i = 1,NSP
           rho_Rstar = rho_Rstar + rhostar(i) * RSTARP(i)
       enddo 
!
!      Temperature
!      External: TREFP = Tin
!      Internal: TREFP = T0

       Tstar = pstar/rho_Rstar 
       T = Tstar * TREFP                           !(K)
!       write(6,*)'T=',T
!
!      Molar concentration (mol/cm^3)     
       do i = 1 , NSP
           rhoc(i) = rhos(i)/Ms(i)     !(kmol/m^3)
       enddo
!      kmol/m^3 -> mol/cm^3
       rhoc = rhoc * 1.0d-3 

!      Number of particles concentration (1/cm^3)     
       do i = 1 , NSP
           rhoN(i) = rhoc(i)*Na 
       enddo

!      Chemical reaction rates  (kf & kb) (cm^3/s) or (cm^6/s)
       do i=1,NREAC
           !Fit type
           if (i.le.3)then 
               fitrate=1       !Colonna
           else
               fitrate=2       !Bacri &  Vlcek
           endif        
           call rate(T,fitrate,i,k1,k2)
           kf(i) = k1
           kb(i) = k2
!           write(6,*)'kf',i,k1
!           write(6,*)'kb',i,k2
!!           pause
       enddo

!      Chemical source terms
       do isp = 1,NSP
           call massaction(rhoN,kf,kb,S,isp)     !  g/(mol*cm^3*s)
           Ss(isp) = S           
       enddo

       do isp = 1,NSP
!          g/(mol*cm^3*s) -> g/(cm^3*s)  (Na=Avogadro number 1/mol)              
           Ss(isp) = Ss(isp)/Na
!          g/(cm^3*s) -> kg/(m^3*s)
           Ss(isp) = Ss(isp) * 1.0d3    
!          Non-dimensional chemical source
!           Sstar(isp) = Ss(isp) * LREFP / (RREFP*UREFP) 
!           write(6,*)'S*=',isp,Sstar(isp)
        enddo      
!
!####################################################################   
!
!     Stanford SU2 rates (Hofferd, Lien and Itikawa)      
      !If Stanford = 1 compute source terms with SU2 rate

!       Stanford = 0  
      
!        if (Stanford.eq.1)then               
!           call massactionSU2(T,rhos,SsSU2)      
!           do isp = 1,NSP
!              Sstar(isp)=SsSU2(isp) * LREFP / (RREFP*UREFP)
!           enddo 
!        endif
!
!###################################################################
!
!     External UREFP = Uinf
!     Internal UREFP = sqrt(RMIX*T0)

!   check \sum_{s=1}^{N_s} S_s = 0
      Shelp = 0.d0
      do i = 1,NSP
         Shelp = Shelp + Ss(i)         
      enddo  
!      Shelp = Shelp/RREFP
!
!      write(6,*)'SUM Ss* =',Shelp
!           
      if(abs(Shelp).gt.1d-10)then
!          write(6,*)'THE SUM OF SOURCE TERMS IS DIFFERENT FROM ZERO'
!          write(6,*)'SUM Ss* =',Shelp
          if(abs(Ss(3)/Ms(3)-Ss(4)/Ms(4))*1.0d3.gt.1d-10)then
!             write(6,*)'S(+)-(-)=',1.0d3*(Ss(3)/Ms(3)-Ss(4)/Ms(4))
!             Ss(4)=Ss(3)*Ms(4)/Ms(3)
!             pause   
          endif
!          Ss(1) = -Ss(2) - Ss(3) - Ss(4)
!         stop                    
      endif  

      DO isp = 1,NSP
           Sstar(isp) = Ss(isp) * LREFP / (RREFP*UREFP) 
      ENDDO

!     Dimensional characteristic chemical time       
      CALL tauchem4Ar(rhoN,kf,kb,tau,Nreac)
      DO isp = 1,NSP
         tau(isp) = tau(isp)*UREFP/LREFP
      ENDDO

       return
       end

      SUBROUTINE tauchem4Ar(rhoN,kf,kb,tau,NREAC)

      IMPLICIT NONE

      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'

      DOUBLE PRECISION rhoN(NSP),kf(NREAC),kb(NREAC),tau(NSP),dSdr(NSP)
      INTEGER NREAC,I

      IF(NREAC.eq.6)THEN   
        dSdr(1)=-rhoN(1)*(kf(4)+kf(5))-rhoN(4)*(kf(1)+kf(2))
        dSdr(2)=-rhoN(1)*(kf(5)+kf(6))-rhoN(4)*(kb(2)+kf(3))
        dSdr(3)=-rhoN(4)*(rhoN(1)*(kb(4)+kb(6)) + rhoN(4)*(kb(1)+kb(3)))
        dSdr(4)=rhoN(1)*kf(1) + rhoN(2)*kf(3)
        dSdr(4)=dSdr(4) - 2.d0*rhoN(3)*rhoN(4)*(kb(1)+kb(3))
        dSdr(4)=dSdr(4) - rhon(1)*rhoN(3)*(kb(4)+kb(6))
        if(abs(dSdr(4)-dSdr(3)).gt.1.d-12)then
          dSdr(4)=dSdr(3)
        endif
      ELSEIF(NREAC.eq.3)THEN
        dSdr(1)=-rhoN(4)*(kf(1)+kf(2))
        dSdr(2)=-rhoN(4)*(kb(2)+kf(3))
        dSdr(3)=-rhoN(4)*(rhoN(4)*(kb(1)+kb(3)))
        dSdr(4)=rhoN(1)*kf(1) + rhoN(2)*kf(3)
        dSdr(4)=dSdr(4) - rhon(1)*rhoN(3)*(kb(4)+kb(6))
        if(abs(dSdr(4)-dSdr(3)).gt.1.d-12)then
           dSdr(4)=dSdr(3)
        endif
       ENDIF 


      DO I=1,NSP
         tau(I) = -1.d0/dSdr(I)
      ENDDO

      RETURN
      END   
       



