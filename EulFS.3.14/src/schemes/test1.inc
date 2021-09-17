      IF (ICHECK.NE.0) THEN
          DO IVERT = 1,NOFVERT 
             IADD = (IVERT-1)*ORDSQR+1
             CALL R8Mat_Print('G',' ',NORDER,NORDER,BETA(IADD),NORDER,
     +                    'Beta matrix',INFO)
          ENDDO
C
C       .. Check that the BETA's sum up to zero ..
C
          LFLAG = UNITMAT(WKSP2,NORDER,NORDER,MAXNOFVAR,1.D-14)
          IF (.NOT.LFLAG) THEN
              WRITE (6,FMT=*) 'Element # ',IELEM
              CALL R8Mat_Print('G',' ',NORDER,NORDER,WKSP2,MAXNOFVAR,
     +                    'Sum of the distribution matrices',INFO)
              PAUSE

          ENDIF
C
      ENDIF
