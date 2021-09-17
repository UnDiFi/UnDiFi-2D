      SUBROUTINE IGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  IGTHR -- INTEGER GATHER                                 ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A INTEGER VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A INTEGER VECTOR  X  IN COMPRESSED FORM (X,INDX).
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C
      INTEGER             NZ, INDX (*)
C
      INTEGER             Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I) = Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
