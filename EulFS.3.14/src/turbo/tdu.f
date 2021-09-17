
C///////////////////////////////////////////////////////////////////////

C       Identification Keywords:

C       @(#)tfv1.f      1.0.0.1   95/09/21        12:00:00
C       /alpha1/users/mubadmin/scalar/sccs/s.tdu.f
C

C///////////////////////////////////////////////////////////////////////

C//////////////////////////////////////////////////////////////////////
C
       DOUBLE PRECISION FUNCTION TDU(UX,UY,UZ)
*  ==============================================================
*
*
*      FILENAME: tdu.f
*
*
*
*      DESCRIPTION:
*      -----------
*
*      function that compute the module difference between the volcity
*      vector in the trip point (if trip point is on wall velocity 
*      in trip point is zero) end the velocity vector in a flowfield
*      point 
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
      DOUBLE PRECISION UX,UY,UZ,UXT,UYT,UZT
C
C     (UXT,UYT,UZT) is the velocity at trip point 
C
      UXT=0.
      UYT=0.
      UZT=0. 
      TDU=SQRT( (UXT-UX)**2+(UYT-UY)**2+(UZT-UZ)**2 )
      RETURN
      END
