       subroutine boltzeq4Ar(P,T,chi,cm)
! 
!      Compute the equilibrium concentrations of Argon internal level species (Ar0, Ar*)
!      
!      $Id: boltzeq4Ar.f,v 1.1 2013/01/25 09:47:36 abonfi Exp $      
!
       implicit none
       
       include 'plasma.h'
       include 'paramchem.h'
       
       double precision P,T,chi,Sm
       double precision ntot,nion,ng,n0,n1,cm
C
       ntot = P/(kBolSI*T)
       ntot = ntot/1.0d6
C
       Sm=6.0d0*exp(-hArm/kBol/T)
C
       nion = chi/(1.0d0+chi)*ntot
       ng = ntot - 2*nion
       n1 = Sm/(1.0d0+Sm)*ng
       n0 = ng - n1
C     
       cm = n1/ntot        
C      
       return
       end    
