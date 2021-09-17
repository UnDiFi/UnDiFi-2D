
C///////////////////////////////////////////////////////////////////////

C       Identification Keywords:

C       @(#)tfv1.f      1.0.0.1   95/08/29        12:00:00
C       /alpha1/users/mubadmin/scalar/sccs/s.tsc.f
C

C///////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////
C
       DOUBLE PRECISION FUNCTION TFT1(TD,TTD,DU)
*  ==============================================================
*
*
*      FILENAME: tfv1.f
*
*
*
*      DESCRIPTION:
*      -----------
*
*      turbulence model function
*      (Spalart & Allmaras model: version)
*
*      FILE HISTORY:
*      ------------
*
*      LANGUAGE : Fortran 77
*      LAST CHANGE : April 97
*
*  ==============================================================

c2345678
      IMPLICIT NONE
      INCLUDE 'turb.com'
      INTEGER IDUM
      DOUBLE PRECISION TD,TTD
      DOUBLE PRECISION DU,TGT
      DOUBLE PRECISION DUM1
      DUM1=TD**2+TGT(DU)**2*TTD**2
      TFT1=TCT1*TGT(DU)*EXP(-TCT2*((TST/DU)**2)*DUM1)
      RETURN
      END
