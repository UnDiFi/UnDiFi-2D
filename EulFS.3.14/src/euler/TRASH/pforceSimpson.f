      SUBROUTINE PFORCE(ICLR,IVERT,VCN,NDIM,VCZ,NOFVAR,NOFVERT,PRESSURE)
C
C     $Id: pforce.f,v 1.5 2008/12/03 11:16:18 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'paramt.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
C
      INTEGER ICLR,IVERT,NDIM,NOFVERT,NOFVAR
C
C     This routine computes the pressure force acting on a
C     boundary face as the aritmetic mean of the face pressure
C     values, ie. assuming linear variation of pressure;
C     while this is o.k. for the INcompressible flow eqns., for
C     COmpressible flows one should compute
C     the EXACT integral of pressure which is a quadratic function
C
C     IVERT is the local nodenumber of the vertex opposite
C           the boundary face
C
C
      DOUBLE PRECISION VCZ(NOFVAR,NOFVERT),VCN(NDIM,NOFVERT)
C
      DOUBLE PRECISION P
      INTEGER I,IV
C
      DOUBLE PRECISION PRESSURE
      INTEGER ICYCL
      EXTERNAL ICYCL,PRESSURE
C
C  .. Loop over the nodes of the face
C
      P = ZERO
      DO 1 I = 1,NOFVERT - 1
          IV = ICYCL(IVERT+I,NOFVERT)
          P = P + PRESSURE(NDIM,VCZ(1,IV))
    1 CONTINUE
C
      PRESF(1,ICLR) = PRESF(1,ICLR) - P*VCN(1,IVERT)/NDIM
      PRESF(2,ICLR) = PRESF(2,ICLR) - P*VCN(2,IVERT)/NDIM
      IF(NDIM.EQ.3)PRESF(3,ICLR) = PRESF(3,ICLR) - P*VCN(3,IVERT)/NDIM
C
      RETURN

      END
