      SUBROUTINE ISCTR ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  ISCTR -- INTEGER SCATTER                                   ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         ISCTR SCATTERS THE COMPONENTS OF 
C             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX) 
C         INTO 
C             SPECIFIED COMPONENTS OF A INTEGER VECTOR  Y  
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
C         X       INTEGER        ARRAY CONTAINING THE VALUES TO BE 
C                             SCATTERED FROM COMPRESSED FORM INTO FULL 
C                             STORAGE FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
C                             TO BE SCATTERED FROM COMPRESSED FORM.  
C                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX 
C                             ARE DISTINCT.
C
C     OUTPUT ...
C
C         Y       INTEGER        ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
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
      INTEGER             X (*), Y (*)
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
