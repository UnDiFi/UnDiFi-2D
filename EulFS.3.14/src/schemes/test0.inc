      IF (ICHECK.NE.0) THEN
C
C       .. Sum up the distribution matrices ..
C
         CALL MATSUM(WKSP2,MAXNOFVAR,BETA(IADD),NORDER,NORDER,NORDER)
C
      ENDIF
