      IF(MATRIX_ASSEMBLY.EQV..FALSE.)GOTO 109
      IF(ICHECK.NE.0)THEN
      GOTO 69
C
C        print the element stiffness matrix
C
      DO 65 J = 1,NOFVERT
          DO 65 I = 1,NOFVERT
              WRITE (6,FMT=*) 'IELEM = ',IELEM,'I,J = ',I,J
              CALL R8Mat_Print('General',' ',NORDER,NORDER,
     +        STIFC(1,1,I,J),NDOF,' Element stiffness matrix ',INFO)
   65 CONTINUE
      PAUSE

   69 CONTINUE
C
C     ... Check that \sum_{j} C_{i,j} = 0
C
      DO 61 I = 1,NOFVERT
          CALL DINIT(MAXNOFVAR*MAXNOFVAR,ZERO,WKSP2,1)
          DO 63 J = 1,NOFVERT
              CALL MATSUM(WKSP2,MAXNOFVAR,STIFC(1,1,I,J),NDOF,NORDER,
     +                    NORDER)
   63     CONTINUE
          LFLAG = NULLMAT(WKSP2,NORDER,NORDER,MAXNOFVAR,1.D-14)
          IF (.NOT.LFLAG) THEN
              WRITE (6,FMT=*) IELEM,I
              CALL R8Mat_Print('General',' ',NORDER,NORDER,WKSP2,
     1             MAXNOFVAR,
     2             ' Sum_j C_{i,j} equals the null matrix ? (LDA)',
     3             INFO)
              PAUSE

          ENDIF

   61 CONTINUE
C
C     ... CELRES[NDOF+1:2*NDOF] stores the residual in 
C         the "V" variable computed as - \sum_{i,j} C_{i,j} V_j
C
      DO 59 J = 1,NOFVERT
          DO 59 I = 1,NOFVERT
              CALL DGEMV('N',NORDER,NORDER,MONE,STIFC(1,1,I,J),NDOF,
     +                   W(1,J),1,ONE,CELRES(NDOF+1),1)
   59 CONTINUE
C
C     check that the two residuals are the same
C
      LFLAG = .FALSE.
      DO 57 IVAR = 1,NORDER
          TEMP = CELRES(IVAR) - CELRES(NDOF+IVAR)
          IF (DABS(TEMP).GT.1.D-14) LFLAG = .TRUE.
   57 CONTINUE
      IF (LFLAG .EQV. .TRUE.) THEN
          WRITE (6,FMT=*) 'System scheme, elem ',IELEM
          WRITE (6,FMT=*) 'var #     implicit        "true"'
          DO 55 IVAR = 1,NORDER
              WRITE (6,FMT=*) IVAR,CELRES(NDOF+IVAR),CELRES(IVAR)
   55     CONTINUE
      ENDIF
      ENDIF
  109 CONTINUE
