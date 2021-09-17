!> \brief \b TFW
!> \par Purpose
!>
!> (Spalart & Allmaras model: version: all)
!>
!> Computes the function:
!>
!> \f[
!> f_{w} = g \, \left(\frac{1+c_{w_3}^6}{g^6+c_{w_3}^6}\right)^{1/6}
!> \f]
!> where:
!> \f[
!> g = r + c_{w_2} \left(r^6-r\right)
!> \f]
!> and
!> \f[
!> r = \frac{\tilde{\nu}}{\tilde{S}\,\kappa^2\,d^2}
!> \f]
!>
!>
!> @param[in] R is \f$ \frac{\tilde{\nu}}{\tilde{S}\,\kappa^2\,d^2} \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.2 $
!> \date $Date: 2015/05/11 07:21:49 $
       DOUBLE PRECISION FUNCTION TFW(R)
c2345678
      INCLUDE 'constants.h'
      INCLUDE 'turb.com'
      DOUBLE PRECISION R
      DOUBLE PRECISION GG,HELP
!     GG=R+TCW2*(R**6.-R)
      GG=R*(ONE+TCW2*(R*R*R*R*R-ONE))
      HELP = TCW3*TCW3*TCW3*TCW3*TCW3*TCW3
      TFW=GG*((ONE+HELP)/(GG**6+HELP))**(ONE/6.d0)
      RETURN
      END
