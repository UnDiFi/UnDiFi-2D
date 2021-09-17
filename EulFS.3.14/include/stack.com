      INTEGER NVA
      PARAMETER(NVA=140000000)
      DOUBLE PRECISION DSTAK
      COMMON /CSTAK/DSTAK(NVA)
C
C     The stack DSTAK(1:NVA) is the unique workarray where
C     all global variables are stored.
C     To increase the size of the problem increase NVA here ONLY.
C     Memory allocation is handled by the port library.
C
