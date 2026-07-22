      SUBROUTINE RSU2(FNAME,LENFNAM,ZROE,NDOF,NPOIN,ICELNOD,NELEM,GAM)
C
C     Reads a SU2 restart CSV (as written by SU2's own solver, in the
C     same "<FNAME>_restart.csv" layout WSU2/triangle2su2 wrote as
C     RESTART_SOL input) and converts the conservative variables back
C     to the Roe parameter vector, per
C     codes/UnDiFi-2D/doc/userguide/040_solver.md:
C       Z = ( sqrt(rho), sqrt(rho)*H, sqrt(rho)*u, sqrt(rho)*v )
C
C     WSU2 wrote a restart row only for points ICELNOD actually
C     references, skipping Triangle's unreferenced phantom-node
C     padding (see wsu2.f) -- rebuild the same ICELNOD-driven
C     compaction here so CSV row order (SU2 preserves it) lines back
C     up with the original point it came from. Points IMAP leaves
C     unmarked never went to SU2, so their ZROE is left as RTRI
C     already loaded it.
C
      IMPLICIT NONE
C
      INTEGER LENFNAM,NDOF,NPOIN,NELEM
      CHARACTER*(*) FNAME
      DOUBLE PRECISION ZROE(NDOF,*),GAM
      INTEGER ICELNOD(3,*)
C
      INTEGER IPOIN,IDUM,IELEM
      INTEGER IMAP(NPOIN)
      DOUBLE PRECISION RHO,RHOU,RHOV,RHOE,U,V,P,H,GM1,XDUM,YDUM
      CHARACTER FWORK*255,LINE*512
C
      GM1 = GAM - 1.D0
C
      DO 5 IPOIN = 1,NPOIN
         IMAP(IPOIN) = 0
    5 CONTINUE
      DO 6 IELEM = 1,NELEM
         IMAP(ICELNOD(1,IELEM)) = 1
         IMAP(ICELNOD(2,IELEM)) = 1
         IMAP(ICELNOD(3,IELEM)) = 1
    6 CONTINUE
C
      FWORK(1:LENFNAM+12) = FNAME(1:LENFNAM)//'_restart.csv'
      OPEN(UNIT=23,FILE=FWORK(1:LENFNAM+12),STATUS='OLD')
      READ(23,FMT='(A)')LINE
C
      DO 10 IPOIN = 1,NPOIN
         IF(IMAP(IPOIN).EQ.0)GOTO 10
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
