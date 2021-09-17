    subroutine rate(T,fitrate,ir,kf,kb)

!   Compute rate coefficients

    implicit none

    include 'plasma.h'
    include 'chem.h'
    include 'paramchem.h'
    include 'commonv.inc'
    include 'commonchem.inc'

!   Variables
    integer i,j,ir,fitrate
    double precision T,Te,kf,kb,Keq
    double precision mf(Nfit),mb(Nfit),meq(Nfit)
    double precision Ep,bp,Qp,Keqr
    double precision f1,f2,f3,fbv
    double precision Tfilter,kbthresh
    integer fitORdetbal,reacoff(NREAC),zeldovich

!   Filter for T
    Tfilter = 100.0d0
    if (T.lt.Tfilter)then 
      T = Tfilter
    endif

!   Threshold for kb2 computed with detailed balance
    kbthresh = 10.d0

!   Flag to select the type of calculation for the backward rate
!   fitORdetbal=1: fit Colonna
!   fitORdetbal=2: detailed balance 
    fitORdetbal = 1


!   Flag to compute electronic excitation and deexcitation rates, with Zel'dovich formula
!   zeldovich=0: fit Colonna for el_Mol processes
!   zeldovich=1: Zeldovich & Raizer formula for el_Mol excitation (G -> 4P) (reaction 2)
    zeldovich = 0
    if (zeldovich.eq.1)then
       fitORdetbal = 2
    endif

!   Flag to turn off chemical ractions i (reacoff=0)
    reacoff(1) = 1      !G -> I   el_Mol 
    reacoff(2) = 1      !G -> 4P  el_Mol
    reacoff(3) = 1      !4P -> I  el_Mol
    reacoff(4) = 1      !G -> I   Mol_Mol
    reacoff(5) = 1      !G -> 4P  Mol_Mol
    reacoff(6) = 1      !4P -> I  Mol_Mol
    
!   fit coefficients    
    do i = 1 , Nfit
        mf(i) = mc(i,ir)
        mb(i) = mc(i,ir+Ncol)
        meq(i) = mc(i,2*Ncol+1)
    enddo       

!   Equilibrium constant 
    Keq = f1(T,meq)             !log(Keq) 
    Keq = exp(Keq)              !Keq (1/cm^3)     

!   Partition function    
    Qp=1+6*exp(-hArm/kBol/T)
          
!   Electronic temperature (K)
    Te=T    

!   Forward e Backward Rate    
    if (fitrate.eq.1)then  
    ! Electron_Atom impact (Colonna fit)
        if (ir.eq.1)then  
            kf = f1(Te,mf)       !log(kf1) 
            kf = exp(kf)        !kf1 (cm^3/s)    
            if (fitORdetbal.eq.1)then            
               kb = f3(Te,mb)       !kb1fit (cm^6/s)
            elseif(fitORdetbal.eq.2)then 
               Keqr = Keq*Qp       
               kb = kf/Keqr       !detailed balance
            endif
            if (reacoff(ir).eq.0)then 
               kf = 0d0
               kb = 0d0 
            endif 
        elseif (ir.eq.2)then   
            if (zeldovich.eq.0)then       !Colonna
               kf = f1(Te,mf)       !log(kf2) 
               kf = exp(kf)        !kf2 (cm^3/s)
            elseif (zeldovich.eq.1)then    !Zeldovich & Raizer              
               kf = 4.348758d-12*sqrt(Te)*(11.55d0+2.0d0*kBol*Te)*exp(-11.55d0/kBol/T)
            endif      
            if (fitORdetbal.eq.1)then                                 
               kb = f1(Te,mb)       !log(kb2) 
               kb = exp(kb)        !kb2fit (cm^3/s)                                
            elseif(fitORdetbal.eq.2)then
               Keqr = 6*exp(-hArm/kBol/T)       
               kb = kf/Keqr        !detailed balance                 
!               if(kb.ge.kbthresh)then
!                  kb = kbthresh              
!               endif 
            endif
            if (reacoff(ir).eq.0)then
               kf = 0d0
               kb = 0d0
            endif         
        elseif (ir.eq.3)then    
            kf = f1(Te,mf)       !log(kf3)
            kf = exp(kf)        !kf3 (cm^3/s)
            if (fitORdetbal.eq.1)then         
               kb = f2(Te,mb)       !log(kb3)
               kb = exp(kb)        !kb3 (cm^6/s)
            elseif(fitORdetbal.eq.2)then
               Keqr = Keq*Qp/(6*exp(-hArm/kBol/T)) 
               kb = kf/Keqr       !detailed balance
            endif
            if (reacoff(ir).eq.0)then
               kf = 0d0
               kb = 0d0
            endif
        endif
    elseif(fitrate.eq.2)then
    ! Atom_Atom impact (Bacri & Vlcek)
        if (ir.eq.4)then
             bp = mbv(1,1)      !(cm^3/(s*eV)/K^0.5)   
             Ep = mbv(2,1)      !eV
             kf = fbv(T,bp,Ep)  !kf4 (cm^3/s)
             Keqr = Keq*Qp
             kb = kf/Keqr        !kb4 (cm^6/s)                                    
             if (reacoff(ir).eq.0)then
               kf = 0d0
               kb = 0d0
            endif
        elseif (ir.eq.5)then
             bp = mbv(1,2)      !(cm^3/(s*eV)/K^0.5)   
             Ep = mbv(2,2)      !eV
             kf = fbv(T,bp,Ep)  !kf4 (cm^3/s)
             Keqr = 6*exp(-hArm/kBol/T) 
             kb = kf/Keqr        !kb4 (cm^6/s)                     
             if (reacoff(ir).eq.0)then
                kf = 0d0
                kb = 0d0
             endif
        elseif (ir.eq.6)then
             bp = mbv(1,3)      !(cm^3/(s*eV)/K^0.5)   
             Ep = mbv(2,3)      !eV
             kf = fbv(T,bp,Ep)  !kf4 (cm^3/s)
             Keqr = Keq*Qp/(6*exp(-hArm/kBol/T)) 
             kb = kf/Keqr        !kb4 (cm^6/s)
            if (reacoff(ir).eq.0)then
               kf = 0d0
               kb = 0d0
            endif
        endif
    endif
!
   return
   end   
   
!   Analytical fitting function f1
    double precision function f1(T,m)
    double precision T,m(5)
!    implicit none     

    f1 = m(1) + m(2)*(1000/T)**m(3) + m(4)*exp(-T/m(5))
!
    return
    end

!   Analytical fitting function f2
    double precision function f2(T,m)
    double precision T,m(5)
    integer im
!    implicit none

    f2 = m(1) + m(2) * log(T) + m(3) * (log(T))**2
!
    return
    end

!   Analytical fitting function f3
    double precision function f3(T,m)
    double precision T,m(5)
!    implicit none

    f3 = m(1)*(T/1000)**m(2) + m(3)*exp((-T/m(4))**m(5)) !+ m(6)*exp(-((T-m(7))/m(8))**2)
!
    return
    end

!   Analytical fitting function fBV 
    double precision function fbv(T,bp,Ep)
    double precision T,bp,Ep,k
    parameter(k=8.6173324d-5)
!    implicit none

    fbv = bp*sqrt(T)*(Ep+2*k*T)*exp(-Ep/(k*T))

    return
    end

