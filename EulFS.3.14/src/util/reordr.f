      SUBROUTINE REORDR(N,DATA,INCX,INDEX,WORK)

C
C     This routine is used to reorder the array DATA
C     according to the permutation array INDEX
C     i.e. A(I) := A(INDEX(I)) 
C     It should follow a previous call to SORTRX or QSORT
C
C     On entry
C     WORK must contain a copy of DATA
C     WORK(i) = DATA(j) i=1,N j = 1,INCX*N,INCX
C     DATA is allowed to have non unit stride
C
C

C     .. Scalar Arguments ..
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      INTEGER DATA(*),INDEX(N),WORK(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,XADDR
C     ..
      XADDR = 1
C
      IF  ( INCX .LT. 0 )  THEN
          XADDR = (-N+1)*INCX + 1
      ENDIF
C
      DO 1 I = 1,N

          DATA(XADDR) = WORK(INDEX(I))
          XADDR       = XADDR + INCX
    1 CONTINUE

      RETURN

      END
