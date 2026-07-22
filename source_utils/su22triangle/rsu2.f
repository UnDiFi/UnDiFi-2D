      SUBROUTINE RSU2(FNAME,LENFNAM,ZROE,NDOF,NPOIN,GAM)
C
C     Reads a SU2 restart CSV (as written by SU2's own solver, in the
C     same "<FNAME>_restart.csv" layout WSU2/triangle2su2 wrote as
C     RESTART_SOL input) and converts the conservative variables back
C     to the Roe parameter vector, per
C     codes/UnDiFi-2D/doc/userguide/040_solver.md:
C       Z = ( sqrt(rho), sqrt(rho)*H, sqrt(rho)*u, sqrt(rho)*v )
C
C     Row order is trusted to match the mesh's own point order (SU2
C     preserves it), the same convention RTRI already relies on for
C     the Triangle .node file -- the leading PointID/x/y columns are
C     read and discarded, not used to re-index.
C
      IMPLICIT NONE
C
      INTEGER LENFNAM,NDOF,NPOIN
      CHARACTER*(*) FNAME
      DOUBLE PRECISION ZROE(NDOF,*),GAM
C
      INTEGER IPOIN,IDUM
      DOUBLE PRECISION RHO,RHOU,RHOV,RHOE,U,V,P,H,GM1,XDUM,YDUM
      CHARACTER FWORK*255,LINE*512
C
      GM1 = GAM - 1.D0
C
      FWORK(1:LENFNAM+12) = FNAME(1:LENFNAM)//'_restart.csv'
      OPEN(UNIT=23,FILE=FWORK(1:LENFNAM+12),STATUS='OLD')
      READ(23,FMT='(A)')LINE
C
      DO 10 IPOIN = 1,NPOIN
         READ(23,*)IDUM,XDUM,YDUM,RHO,RHOU,RHOV,RHOE
         U = RHOU/RHO
         V = RHOV/RHO
         P = GM1*(RHOE - 0.5D0*RHO*(U*U+V*V))
         H = (RHOE + P)/RHO
         ZROE(1,IPOIN) = SQRT(RHO)
         ZROE(2,IPOIN) = SQRT(RHO)*H
         ZROE(3,IPOIN) = SQRT(RHO)*U
         ZROE(4,IPOIN) = SQRT(RHO)*V
   10 CONTINUE
      CLOSE(23)
C
      RETURN
      END
