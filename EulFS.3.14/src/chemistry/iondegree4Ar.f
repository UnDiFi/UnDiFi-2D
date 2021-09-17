       subroutine iondegree4Ar(P,T,chi)
! 
!      Implementation of Saha equation to compute the degree of ionization 
!      Mixture of Ar,Ar+ and e- at equilibrium  
!
!      $Id: iondegree4Ar.f,v 1.1 2013/01/25 09:48:16 abonfi Exp $
!
       implicit none
       
       include 'plasma.h'
       include 'paramchem.h'
       
       double precision P,T,chi
       double precision ntot,nion,beta,const
       parameter (const=2.4d15)
C
       ntot = P/(kBolSI*T)
       ntot = ntot/1.0d6
C
       beta = const * T**(3/2) * exp(-hArI/(kBol*T))/ntot
       chi = (-beta + sqrt(beta**2 + 4*beta))/2
C 
       return
       end    
