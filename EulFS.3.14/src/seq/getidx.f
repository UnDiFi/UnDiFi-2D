      SUBROUTINE GetIdx(NGHOST,NOFVAR,GHOSTS,GHOSTS2)
c
      IMPLICIT NONE
c
c     nghost is the number of ghost meshpoints
c     ghosts stores the nodal indexes (1-based indexing)
c            of the ghost meshpoints
c     ghosts2 stores the adresses (0-based indexing)
c              of the unknowns of the ghost meshpoints
c

C     .. Scalar Arguments ..
      INTEGER NGHOST,NOFVAR
C     ..
C     .. Array Arguments ..
      INTEGER GHOSTS(NGHOST),GHOSTS2(NOFVAR,NGHOST)
C     ..
C     .. Local Scalars ..
      INTEGER I,IPOIN,J
C     ..
      DO 1 I = 1,NGHOST
          IPOIN = GHOSTS(I)
          DO 1 J = 1,NOFVAR
              GHOSTS2(J,I) = (IPOIN-1)*NOFVAR + J - 1
    1 CONTINUE
      RETURN

      END
