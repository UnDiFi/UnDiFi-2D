!> \brief \b TFV2
!> \par Purpose
!>
!> compute the function \f$ f_{v_2} = \frac{1}{\left(1+\chi/c_{v_2}\right)^3} \f$
!> (Spalart & Allmaras model: version: SA-fv3)
!>
!>
!> @param[in] TCHI is \f$ \tilde{\nu}/\nu \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2020/02/05 15:03:28 $
      DOUBLE PRECISION FUNCTION TFV2(TCHI)
c2345678
      IMPLICIT NONE
      INCLUDE 'turb.com'
      INCLUDE 'constants.h'
      DOUBLE PRECISION TCHI

C     new version from the Ashford's thesis; this is also known as the SA-fv3 variant of the model
C
      TFV2=ONE/(ONE+TCHI/TCV2)**3
      RETURN
      END
!> \brief \b TFV2SA
!> \par Purpose
!>
!> compute the function \f$ f_{v_2} = 1-\frac{\chi}{\left(1+\chi\,f_{v_1}\right)} \f$
!> (Spalart & Allmaras model: version: SA)
!>
!>
!> @param[in] TCHI is \f$ \tilde{\nu}/\nu \f$
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2020/02/05 15:03:28 $
      DOUBLE PRECISION FUNCTION TFV2SA(TCHI)
c2345678
      IMPLICIT NONE
      INCLUDE 'turb.com'
      INCLUDE 'constants.h'
      DOUBLE PRECISION TCHI
      DOUBLE PRECISION TFV1
C
      TFV2SA=ONE-TCHI/(ONE+TCHI*TFV1(TCHI))
      RETURN
      END
