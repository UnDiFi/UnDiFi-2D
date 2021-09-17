      SUBROUTINE readArgon()

      IMPLICIT NONE
      INCLUDE 'plasma.h'
      INCLUDE 'chem.h'
      INCLUDE 'paramchem.h'
      INCLUDE 'commonv.inc'

      INTEGER ISP,IR,IFIT 
      DOUBLE PRECISION eV2KJoMOL  
      PARAMETER (eV2KJoMOL=96.485d0)
      CHARACTER(len=255) :: homedir,fname 

      CALL get_environment_variable("FSPL_DIR", homedir)

!     Stoichiometric coefficients (reagents)
      fname = TRIM(homedir)//'/src/chemistry/data/Ar/creagent.dat'
!     OPEN(10,FILE='creagent.dat')
      OPEN(10,FILE=fname)
        DO ISP = 1,NSP
            READ(10,*)(scr(ISP,IR),IR=1,NRMAX)
        ENDDO
      CLOSE(10)

!     Stoichiometric coefficients (products)
      fname = TRIM(homedir)//'/src/chemistry/data/Ar/cproduct.dat'
!     OPEN(10,FILE='cproduct.dat')
      OPEN(10,FILE=fname)
        DO ISP = 1,NSP
            READ(10,*)(scp(ISP,IR),IR=1,NRMAX)
        ENDDO
      CLOSE(10)

!     Rate fit coefficients (atom-electron collisions - Colonna) 
      fname = TRIM(homedir)//'/src/chemistry/data/Ar/fitdata.dat'
!     OPEN(10,FILE='fitdata.dat')
      OPEN(10,FILE=fname)
        DO IFIT = 1,NFIT
             READ(10,*)(mc(IFIT,IR),IR=1,2*Ncol+1)
        ENDDO
      CLOSE(10)

!      Rate fit coefficients (atom-atom collisions - Bacri & Vlcek) 
      fname = TRIM(homedir)//'/src/chemistry/data/Ar/fitdatabacri.dat'
!     OPEN(10,FILE='fitdatabacri.dat')
      OPEN(10,FILE=fname)
         DO IFIT = 1,2
              READ(10,*)(mbv(IFIT,IR),IR=1,Nbac)
         ENDDO
      CLOSE(10)

!     Molar masses (g/mol)
      Ms(1) = MAr              !Ar
      Ms(2) = MAr              !Ar*
      Ms(3) = MAr - Me         !Ar+
      Ms(4) = Me               !e-

!     Heat of formation (eV)
      HF0(1) = 0.d0            !Ar
      HF0(2) = hArm            !Ar*
      HF0(3) = hArI            !Ar+            
      HF0(4) = 0.d0            !e-

!     eV/part -> J/kg
      DO ISP = 1,NSP
          HF0(ISP) = HF0(ISP) * eV2KJoMOL / Ms(ISP) * 1.0d6
      ENDDO     

!     Rgas (J * K^-1 * kg^-1)
      DO ISP = 1,NSP
          RGASS(ISP) = Rgp / Ms(ISP) * 1.0d3
      ENDDO           

!      write(6,*)'HF0=',HF0
!      write(6,*)'RGASS=',RGASS
!      pause
      RETURN
      END
