          IF( PERIODIC_MESH .AND. ANNULAR )THEN
             DO 90 JV = 1, NOFVERT
                IF(ICELNOD(JV,NELEM+IELEM).GT.(NPOIN+NGHOST))THEN
                    PFLAG(JV)=.TRUE.
                ELSE
                    PFLAG(JV)=.FALSE.
                ENDIF
   90        CONTINUE
          ENDIF
C
