      SUBROUTINE   SAXPYI   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SAXPYI -- INDEXED REAL ELEMENTARY VECTOR OPERATION      ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SAXPYI ADDS A REAL SCALAR MULTIPLE OF 
C             A REAL SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A REAL VECTOR  Y  IN FULL STORAGE FORM.
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
C         A       REAL        SCALAR MULTIPLIER OF  X.
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
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
      REAL                Y (*), X (*), A
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
      IF  ( A .EQ. 0.0E0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      REAL FUNCTION   SDOTI   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SDOTI -- REAL INDEXED DOT PRODUCT                       ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SDOTI COMPUTES THE VECTOR INNER PRODUCT OF 
C             A REAL SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM.
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
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         SDOTI   REAL        REAL FUNCTION VALUE EQUAL TO THE 
C                             VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  SDOTI IS SET TO ZERO.
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
      REAL                X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      SDOTI = 0.0E0
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          SDOTI = SDOTI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SGTHR ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  SGTHR -- REAL GATHER                                    ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SGTHR GATHERS THE SPECIFIED ELEMENTS FROM 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A REAL VECTOR  X  IN COMPRESSED FORM (X,INDX).
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
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     OUTPUT ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES GATHERED INTO
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
      REAL                Y (*), X (*)
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
      SUBROUTINE SGTHRZ ( NZ, Y, X, INDX )
C
C     ==================================================================
C     ==================================================================
C     ====  SGTHRZ -- REAL GATHER AND ZERO                          ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A REAL VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A REAL VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
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
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES GATHERED INTO
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
      REAL                Y (*), X (*)
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
          Y(INDX(I)) = 0.0E0
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE   SROTI   ( NZ, X, INDX, Y, C, S )
C
C     ==================================================================
C     ==================================================================
C     ====  SROTI  --  APPLY INDEXED REAL GIVENS ROTATION           ====
C     ==================================================================
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C     PURPOSE
C     -------
C
C         SROTI APPLIES A GIVENS ROTATION TO 
C             A SPARSE VECTOR   X  STORED IN COMPRESSED FORM  (X,INDX) 
C         AND 
C             ANOTHER VECTOR  Y  IN FULL STORAGE FORM.
C
C         SROTI DOES NOT HANDLE FILL-IN IN  X  AND THEREFORE, IT IS
C         ASSUMED THAT ALL NONZERO COMPONENTS OF  Y  ARE LISTED IN 
C         INDX.  ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN
C         INDX  ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST
C         BE DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
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
C         INDX    INTEGER     ARRAY CONTAINING THE INDICIES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C         C,S     REAL        THE TWO SCALARS DEFINING THE GIVENS
C                             ROTATION.
C
C     UPDATED ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE 
C                             SPARSE VECTOR IN COMPRESSED FORM.
C         Y       REAL        ARRAY WHICH CONTAINS THE VECTOR  Y
C                             IN FULL STORAGE FORM.  ONLY THE
C                             ELEMENTS WHOSE INDICIES ARE LISTED IN
C                             INDX  HAVE BEEN REFERENCED OR MODIFIED.
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
      REAL                X (*), Y (*), C, S
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
      REAL                TEMP
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( ( C .EQ. 1.0E0 ) .AND. ( S .EQ. 0.0E0 ) )  RETURN
C
      DO 10 I = 1, NZ
          TEMP        = - S * X (I)  +  C * Y (INDX(I))
          X (I)       =   C * X (I)  +  S * Y (INDX(I))
          Y (INDX(I)) = TEMP
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SSCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  SSCTR -- REAL SCATTER                                   ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         SSCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A REAL VECTOR  Y  
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
C         X       REAL        ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       REAL        ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
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
      REAL                X (*), Y (*)
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
