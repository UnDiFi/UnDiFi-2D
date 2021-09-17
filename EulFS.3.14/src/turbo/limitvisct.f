      SUBROUTINE limitvisct(zroe,nofvar,nturb,npoin)
      IMPLICIT NONE
      INCLUDE 'constants.h'
      INTEGER nofvar,nturb,npoin
      DOUBLE PRECISION ZROE(nofvar,*)
      INTEGER ipoin
      DOUBLE PRECISION S
      DO ipoin = 1,NPOIN
         S = ZROE(nofvar,ipoin) 
         ZROE(nofvar,ipoin) = MAX(S,ZERO)
      ENDDO
      RETURN
      END
