C
C     $Id: periodic.com,v 1.3 2013/04/30 07:12:17 abonfi Exp $
C
      DOUBLE PRECISION PITCH,COSALPHA,SINALPHA
      DOUBLE PRECISION CYY,CZZ,CYZ,CZY
      DOUBLE PRECISION QMAT(MAX_NOFVAR_SQR)
      LOGICAL PERIODIC_MESH,ANNULAR
      LOGICAL PFLAG(MAXNOFVERT)
      INTEGER NBLADES

      COMMON/PCOML/PERIODIC_MESH,ANNULAR,PFLAG
      COMMON/PCOMR/PITCH,COSALPHA,SINALPHA,CYY,CZZ,CYZ,CZY,
     &QMAT
      COMMON/PCOMI/NBLADES
