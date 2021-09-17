C
      SUBROUTINE PLASFLOW_INF(P,T,CHI,RHOS,RHO,RMIX,RS,HFTOT,CM)
C
C     p,T,chi,cm ---> rhos,rho,Rmix,Hf0 
C
C     Purpose:
C     1) computing global and single species density for a mixture of Ar,Ar*,Ar+,e-
C     2) computing constant of the mixture & global formation enthalpy
C     
C     Input
C     p: pressure
C     T: temperature
C     chi: ionization degree
C     cm: metastable concentration
C
C     Output
C     Rhos(NSP): chemical species denities
C     Rho: density
C     Rmix: mixture gas constant 
C     Rs(NSP): single species gas constant
C     HFTOT: mixture formation enthalpy
C                
C      
      IMPLICIT NONE     
C
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'paramchem.h'
C
      INTEGER ISP
      DOUBLE PRECISION P,T,CHI,RHOS(NSP),RHO,RMIX,HF0(NSP),HFTOT,CM
      DOUBLE PRECISION NION,NTOT,N0,NS(NSP),Ms(NSP),RS(NSP)
      DOUBLE PRECISION eV2KJoMOL
      PARAMETER (eV2KJoMOL=96.485d0)

C     Molar masses (g/mol)
      Ms(1) = MAr              !Ar
      Ms(2) = MAr              !Ar*
      Ms(3) = MAr - Me         !Ar+
      Ms(4) = Me               !e-

C     Heat of formation (eV)
      HF0(1) = ZERO            !Ar
      HF0(2) = hArm            !Ar*
      HF0(3) = hArI            !Ar+            
      HF0(4) = ZERO            !e-

C     eV -> J/kg
      DO ISP = 1,NSP
          HF0(ISP) = HF0(ISP) * eV2KJoMOL / Ms(ISP)*1.0d6
      ENDDO

C     Rgas (J * K^-1 * kg^-1)
      DO ISP = 1,NSP
          RS(ISP) = Rgp / Ms(ISP) * 1.0d3
      ENDDO
C
C     Molar concentration (mol/m^3)
      NTOT = P / (Rgp * T)
      NION = CHI / (CHI + ONE) * NTOT
!      NION = CHI * NTOT        
C
      NS(2) = CM*NTOT
!      NS(2) = 0.d0
      NS(1) = NTOT - TWO * NION - NS(2)
      NS(3) = NION
      NS(4) = NION      
!      write(6,*)'HF0=',HF0
C
C     Densitiy of chemical species (kg/m^3)
      DO ISP = 1,NSP
         RHOS(ISP) = NS(ISP) * MS(ISP) * 1.0d-3
      ENDDO
C
C     Total density (kg/m^3)
      RHO = 0.d0
      DO ISP=1,NSP
        RHO = RHO + RHOS(ISP)
      ENDDO
C
C     Total formation energy (J/kg)
      HFTOT = 0.d0
      DO ISP =  1,NSP
        HFTOT = HFTOT + RHOS(ISP) * HF0(ISP)
      ENDDO
      HFTOT = HFTOT / RHO

!      write(6,*)'HFTOT=',HFTOT
!      pause
C
C     Constant of the mixture   (J/kg K)
      RMIX = ZERO
      DO ISP = 1,NSP
        RMIX = RMIX + RHOS(ISP) * Rgp / MS(ISP) * 1.0d3
C        RMIX = RMIX + RHOS(ISP) * RS(ISP) 
      ENDDO
      RMIX = RMIX / RHO 
!      write(6,*)'RMIX=',RMIX
C
C
      RETURN
      END                
