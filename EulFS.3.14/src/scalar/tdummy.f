      SUBROUTINE TDUMMY(iarg1,arg2,arg3,iarg4,iarg5,iarg6,iarg7,
     &                  arg8,arg9,arg10,arg11,larg12,sarg13,sarg14,
     &                  arg15,larg16,VISCL,VISCT)
C
C     $Id: tdummy.f,v 1.4 2013/01/26 11:14:27 abonfi Exp $
C
C     A dummy turbulence model that returns the
C     constant diffusion coefficient for scalar
C     advection-diffusion problems
C
C     This routine is called by the higher level routine
C     noname()
C

      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'visco.com'
C
C     .. Scalar Arguments ..
      INTEGER iarg1,iarg4,iarg5,iarg6,iarg7
      DOUBLE PRECISION VISCL,VISCT,arg2,arg3,arg8,arg9
      DOUBLE PRECISION arg10,arg11,arg15
      LOGICAL larg12,larg16
      EXTERNAL sarg13,sarg14
C     ..
C     ..
      VISCL = REINV
      VISCT = ZERO

      RETURN

      END
