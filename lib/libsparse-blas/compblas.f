      SUBROUTINE   CAXPYI   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  CAXPYI -- INDEXED COMPLEX ELEMENTARY VECTOR OPERATION   ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CAXPYI ADDS A COMPLEX SCALAR MULTIPLE OF 
C             A COMPLEX SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       COMPLEX     SCALAR MULTIPLIER OF  X.
C         X       COMPLEX     ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
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
      INTEGER             NZ, INDX (*)
C         
      COMPLEX             Y (*), X (*), A
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
      IF  ( A .EQ. ( 0.0E0, 0.0E0 ) )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      COMPLEX FUNCTION   CDOTCI   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  CDOTCI -- COMPLEX CONJUGATED INDEXED DOT PRODUCT        ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CDOTCI COMPUTES THE CONJUGATED VECTOR INNER PRODUCT OF 
C             A COMPLEX SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       COMPLEX     ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         CDOTCI   COMPLEX    COMPLEX FUNCTION VALUE EQUAL TO THE 
C                             CONJUGATED VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  CDOTCI IS SET TO ZERO.
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
      INTEGER             NZ, INDX (*)
C
      COMPLEX             X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      CDOTCI = ( 0.0E0, 0.0E0 )
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          CDOTCI = CDOTCI  +  CONJG ( X(I) ) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      COMPLEX FUNCTION CDOTUI ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  CDOTUI -- COMPLEX UNCONJUGATED INDEXED DOT PRODUCT      ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CDOTUI COMPUTES THE UNCONJUGATED VECTOR INNER PRODUCT OF 
C             A COMPLEX SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       COMPLEX     ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         CDOTUI   COMPLEX    COMPLEX FUNCTION VALUE EQUAL TO THE 
C                             UNCONJUGATED VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  CDOTCI IS SET TO ZERO.
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
      INTEGER             NZ, INDX (*)
C
      COMPLEX             X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      CDOTUI = ( 0.0E0, 0.0E0 )
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          CDOTUI = CDOTUI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
          RETURN
      END
      SUBROUTINE CGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  CGTHR -- COMPLEX GATHER                                 ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A COMPLEX VECTOR  X  IN COMPRESSED FORM (X,INDX).
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
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       COMPLEX     ARRAY CONTAINING THE VALUES GATHERED INTO
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
      COMPLEX             Y (*), X (*)
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
      SUBROUTINE CGTHRZ ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  CGTHRZ -- COMPLEX GATHER AND ZERO                       ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A COMPLEX VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     UPDATED ...
C
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       COMPLEX     ARRAY CONTAINING THE VALUES GATHERED INTO
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
      INTEGER             NZ, INDX (*)
C
      COMPLEX             Y (*), X (*)
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
          X(I)       = Y(INDX(I))
          Y(INDX(I)) = ( 0.0E0, 0.0E0 )
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE CSCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  CSCTR -- COMPLEX SCATTER                                ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         CSCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A COMPLEX VECTOR  Y  
C             IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE MODIFIED.  THE VALUES IN  INDX  MUST BE DISTINCT TO
C         ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE SCATTERED FROM 
C                             COMPRESSED FORM.
C         X       COMPLEX     ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       COMPLEX     ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
C                             HAVE BEEN SET TO THE CORRESPONDING 
C                             ENTRIES OF  X.  ONLY THE ELEMENTS  
C                             CORRESPONDING TO THE INDICES IN  INDX
C                             HAVE BEEN MODIFIED.
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
      INTEGER             NZ, INDX (*)
C
      COMPLEX             X (*), Y (*)
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
          Y(INDX(I)) = X(I)
   10 CONTINUE
C
      RETURN
      END
