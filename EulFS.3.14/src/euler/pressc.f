      DOUBLE PRECISION FUNCTION PRESSC( NDIM , ZROE )
C
C    .. This function computes PRESSURE from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'pfcgas.com'
C
      INTEGER NDIM
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION TEMP
C
      TEMP       = ZROE(3)*ZROE(3) + ZROE(4)*ZROE(4)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(5)*ZROE(5)
      TEMP = HALF * TEMP
      PRESSC = GM1OG * ( ZROE(1)*ZROE(2) - TEMP )
C
      RETURN
      END
C
       DOUBLE PRECISION FUNCTION PRESS4AR( NDIM , ZROE )
C
C    .. This function computes PRESSURE from Roe's
C       parameter vector ..
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
      INCLUDE 'pfcgas.com'
      INCLUDE 'dofs.com'
C
      INTEGER NDIM,ISP
      DOUBLE PRECISION ZROE(*)
      DOUBLE PRECISION SQRTR1,HELP1,PRESSC,ZRDR1,INVDENS,DR1(NSP)
      DOUBLE PRECISION HELP,ZRDR,TEMP,ZRCHI,ERRORP
      DOUBLE PRECISION SQRTR,DR(NSP),DE,DM(3),CHI(NSP),KAPPA
C
C
      KAPPA = GM1   
      DE = KAPPA
C
      HELP1 = KAPPA/(ONE + KAPPA)
      HELP = ONE/(ONE+DE)
C
      TEMP       = ZROE(IX)*ZROE(IX) + ZROE(IY)*ZROE(IY)
      IF( NDIM .EQ. 3 )TEMP = TEMP + ZROE(IZ)*ZROE(IZ)
      TEMP = HALF * TEMP
C
C     VERIFICA
C
      SQRTR = ZERO
      DO ISP = 1,NSP
         SQRTR = SQRTR + ZROE(ISP)
      ENDDO
C
      INVDENS = ONE/SQRTR
      INVDENS = INVDENS*INVDENS 
C
      ZRDR = ZERO
      ZRCHI = ZERO
      DO ISP = 1,NSP         
!         write(*,*) DR(ISP)
         DR(ISP)=CHI(ISP)+KAPPA*TEMP*INVDENS  
         ZRDR = ZRDR + ZROE(ISP)*DR(ISP)
!         write(*,*) DR(ISP)
         ZRCHI = ZRCHI + ZROE(ISP)*CHI(ISP)
      ENDDO 
      ZRDR = ZRDR * SQRTR  
      ZRCHI = ZRCHI * SQRTR
C      
      PRESS4Ar = HELP *(ZRDR + DE*(SQRTR*ZROE(IE) - TWO*TEMP))  
!      write(*,*)'press4Ar',PRESS4Ar
C   
!rpepe manca il termine dovuto a ZRDR  PRESSC = HELP1 * (ZROE(IE)*SQRTR - TEMP)
!      write(*,*)'pressc Corretta',PRESSC      
C      
      PRESSC = ZRCHI/(ONE+KAPPA)+ HELP1 * (ZROE(IE)*SQRTR - TEMP)
!      write(*,*)'pressc',PRESSC      

!      ERRORP = (PRESS4Ar-PRESSC)/PRESS4Ar
!      write(*,*)'errorP ',ERRORP 
!      pause

      RETURN
      END
C 
