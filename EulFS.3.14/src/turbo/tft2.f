
C///////////////////////////////////////////////////////////////////////

C       Identification Keywords:

C       @(#)tfv2.f      1.0.0.1   95/09/21        12:00:00
C       /alpha1/users/mubadmin/scalar/sccs/s.tft2.f
C

C///////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////
C
       REAL*8 FUNCTION TFT2(TCHI)
*  ==============================================================
*
*
*      FILENAME: tft2.f
*
*
*
*      DESCRIPTION:
*      -----------
*
*      turbulence model function
*      (Spalart & Allmaras model: 2D mgship version)
*
*      FILE HISTORY:
*      ------------
*
*      LANGUAGE : Fortran 77
*      LAST CHANGE : April 97
*
*  ==============================================================

c2345678
      INCLUDE 'turb.com'
      DOUBLE PRECISION TCHI
      TFT2=TCT3*EXP(-TCT4*TCHI**2)
      RETURN
      END
