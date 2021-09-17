      SUBROUTINE LAMINAR(iarg1,arg2, arg3,iarg4, iarg5, iarg6,iarg7,
     &                    arg8,arg9,arg10,arg11,larg12,sarg13,sarg14,
     &                   arg15,COMPRESSIBLE,VISCL,VISCT)
c
c     $Id: laminar.f,v 1.4 2013/01/26 11:56:44 abonfi Exp $
c
c     this is a dummy turbulence model:
c     returns laminar viscosity
c
c
      implicit none
c
      include 'paramt.h'
      include 'constants.h'
      include 'three.com'
c
c     input:
c     ------
      INTEGER iarg1,iarg4,iarg5,iarg6,iarg7
      DOUBLE PRECISION arg2,arg3,arg8,arg9
      DOUBLE PRECISION arg10,arg11,arg15
      LOGICAL larg12,COMPRESSIBLE
      EXTERNAL sarg13,sarg14
c
c     output:
c     ------
      double precision viscl,visct
c
c     local:
c     -----
c
      double precision dummy
      double precision sutherlaw
c
      if(compressible)then
      viscl = sutherlaw(dummy,abar,asqr)
      else
      viscl = ONE
      endif
      visct = ZERO
c
      return
      end
