      DOUBLE PRECISION	BETA,BETASQR,nu_p,nu_m,X
      COMMON / blkEulerII / BETA,BETASQR,nu_p,nu_m,X
C
C     $Id: blkEulerII.com,v 1.1 2013/01/25 08:12:46 abonfi Exp $
C
C     Variables used in the Hyperbolic-Elliptic Splitting
C     for compressible flows
