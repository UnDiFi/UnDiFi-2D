
C///////////////////////////////////////////////////////////////////////

C

C///////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////
C
       DOUBLE PRECISION FUNCTION TFV1(TCHI)
       IMPLICIT NONE
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
*      LAST CHANGE : October 97
*
*  ==============================================================

c2345678
      INCLUDE 'turb.com'
      DOUBLE PRECISION TCHI

      TFV1=TCHI**3/(TCHI**3+TCV1**3)
      RETURN
      END
