    subroutine massactionSU2(T,rhos,Ss)

!######################################################################################
!
!   Compute source term for each chemical species using the model implemented in SU2  
!
!   Chemical rates are given by:
!   Hoffert, Lien - "Quasi-one dimensional nonequilibium gas dynamics of 
!   partially ionized 2-temperature Argon", The Physiscs of Fluids, 10, 1967
!
!######################################################################################    

    implicit none
    include 'plasma.h'
    include 'paramchem.h'

    integer isp
    double precision T,rhos(NSP),Ss(NSP),Chis(NSP),Ms(NSP),RR,qr,qp
    double precision kf(NSP),kb(NSP),Keq,Cs(NSP)
    double precision Csm,Cse,eta,theta,Cks,zeta,phi
    double precision Tkr,Tkeq
    parameter(Csm=10.12d0,Cse=22.59d4)
!    parameter(Csm=10.12d0,Cse=22.59)
    parameter(eta=1.5d0,theta=1.353d6,Cks=2.9d22,zeta=1.5,phi=1.831d6)

!   Molecular weight
    Ms(1) = MAr!*1.0d-3
    Ms(2) = MAr!*1.0d-3
    Ms(3) = (MAr - Me)!*1.0d-3
    Ms(4) = Me!*1.0d-3

!   Rate constants
    Cs(1) = Csm
    Cs(2) = Csm
    Cs(3) = Csm
    Cs(4) = Cse

!   Adimensional Temperature
    Tkr = theta/T
    Tkeq = phi/T

!   Equilibrium constant (m^-3)
    Keq = Cks*T**zeta*exp(-TKeq)    
 
    do isp = 1,NSP    
!      Forward rates (m^3/s)
       kf(isp) = Cs(isp)*T**eta*(Tkr + 2.d0)*exp(-Tkr) 
!      Backward rates(m^6/s) 
       kb(isp) = kf(isp)/Keq                     

!      Molar densities
       Chis(isp) = rhos(isp)/Ms(isp)        

   enddo

!   Rate of reaction
    RR = 0.d0
    do isp = 1,NSP
       qp = kf(isp)*(Chis(1)+Chis(2))*Chis(isp)
       qr = kb(isp)*(Chis(3)*Chis(4)*Chis(isp))
       RR = RR + (qr-qp)
    enddo

!   Source terms
    Ss(1) = Ms(1)*RR
    Ss(2) = 0.d0
    Ss(3) = -Ms(3)*RR     
    Ss(4) = -Ms(4)*RR

    return
    end
