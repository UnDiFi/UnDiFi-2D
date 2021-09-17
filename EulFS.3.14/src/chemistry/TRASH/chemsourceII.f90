!****************************************************************************
!  PURPOSE:  To compute the chemical source term for an Argon mixture.
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
!   VARIABLES
!
!   Ns,Nr: species & chemical reactions number 
!   scr,scp: stoichiometric coefficient (reagents,products)
!   kf,kb: forward and backward rates
!   Keq: chemical equilibrium constant (particles concentration) 
!   rhos(i): density species i
!   rhoc(i): molar concentration species i
!   rhoN(i): particle concentration species i
!   Ss(i): source term species i
!
!   INPUT
!   T: dimensional temperature
!   press: dimensional pressure
!   rhos(i): chemical species density
!
!   OUTPUT
!   Ss(i): dimensional source term     
!    
!****************************************************************************

       subroutine chemsourceII(T,press,rhos,Ss)

       implicit none

       ! Variables   

       include 'paramt.h' 
       include 'plasma.h'
       include 'chem.h'
       include 'paramchem.h'
       include 'commonv.inc'
       include 'streamplasma.com'

      integer i,j,isp,ir
      double precision press,T
      double precision rhos(NSP),rhoc(NSP),rhoN(NSP)
      double precision k1,k2,kf(NR),kb(NR)
      double precision Ss(NSP),S,Shelp
      integer fitrate

!       write(6,*)'SEI IN CHEMSOURCEII'
!       write(6,*)'scr',scr         
!       write(6,*)'scp',scp
!       write(6,*)'mc',mc         
!       write(6,*)'mbv',mbv         
       
!     Molar masses (g/mol)
      Ms(1) = MAr              !Ar
      Ms(2) = MAr              !Ar*
      Ms(3) = MAr - Me         !Ar+
      Ms(4) = Me               !e-

!       write(6,*)'molar masses',Ms
         
!      Molar concentration (mol/cm^3)     
       do i = 1 , NSP
           rhoc(i) = rhos(i)/Ms(i)     !(kmol/m^3)
       enddo

!      mol/m^3 -> mol/cm^3
       rhoc = rhoc * 1.0d-3 
       
!       write(6,*)'molar densities',rhoc

!      Number of particles concentration (1/cm^3)     
       do i = 1 , NSP
           rhoN(i) = rhoc(i)*Na 
       enddo

!       write(6,*)'particles densities',rhoN

!      Chemical reaction rates  (kf & kb) (cm^3/s) or (cm^6/s)
       do i=1,NR
           !Fit type
           if (i.le.3)then 
               fitrate=1       !Colonna
           else
               fitrate=2       !Bacri &  Vlcek
           endif        
           call rate(T,fitrate,i,k1,k2)
           kf(i) = k1
           kb(i) = k2
!       write(*,*)'kf reaction',i,'=',k1
!       write(*,*)'kb reaction',i,'=',k2
       enddo

!      Chemical source 
       do i=1,NSP
           call massaction(rhoN,kf,kb,S,i)     !  g/(mol*cm^3*s)
           Ss(i)=S
       enddo

!       write(6,*),'chemical source',Ss        
  
!      g/(mol*cm^3*s) -> g/(cm^3*s)  (Na=Avogadro number 1/mol)      
       Ss = Ss/Na

!      g/(cm^3*s) -> kg/(m^3*s)
       Ss = Ss * 1.0d3
        
!       write(6,*),'chemical source',Ss   
     
!     check \sum_{s=1}^{N_s} S_s = 0
      Shelp = 0.d0
      do i = 1,NSP
         Shelp = Shelp + Ss(i)
      enddo  

      if(abs(Shelp).gt.1d-5)then
!          write(6,*)'THE SUM OF SOURCE TERMS IS DIFFERENT FROM ZERO'
!          write(6,*)'SUM Ss* =',Shelp
!          pause
      endif  

       return
       end

