C
C    $Id: datatype.h,v 1.1 2013/01/25 08:20:21 abonfi Exp $
C
C    1 - LOGICAL
C    2 - INTEGER
C    3 - REAL
C    4 - DOUBLE PRECISION
C    5 - COMPLEX
C
      INTEGER KIND_LOGICAL,KIND_INTEGER,KIND_REAL4,KIND_REAL8,KIND_CMPLX
      PARAMETER (KIND_LOGICAL=1,KIND_INTEGER=2,KIND_REAL4=3,
     &           KIND_REAL8  =4,KIND_CMPLX  =5)
C
C     refer to the port library for the different datatypes
C
