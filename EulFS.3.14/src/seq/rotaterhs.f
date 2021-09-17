      IF( PERIODIC_MESH .AND. ANNULAR )THEN
          DO 19 IVERT = 1, NOFVERT
C     check if a periodic node, if not jump
C     to the following vertex
             IF( PFLAG(IVERT) )THEN
C     pre-multiply by Qt
                CYY = NODRES(IY,IVERT)
                CZZ = NODRES(IZ,IVERT)
                NODRES(IY,IVERT) =   CYY*COSALPHA+CZZ*SINALPHA
                NODRES(IZ,IVERT) = - CYY*SINALPHA+CZZ*COSALPHA
             ENDIF
   19     CONTINUE
      ENDIF
