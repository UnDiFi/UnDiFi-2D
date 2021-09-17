C
C     $Id: bodyf.com,v 1.1 2013/01/25 08:13:34 abonfi Exp $
C
      DOUBLE PRECISION GRAV(3)
      COMMON /BDYCOM/ GRAV
