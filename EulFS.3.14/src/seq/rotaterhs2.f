      IF( PERIODIC_MESH .AND. ANNULAR )THEN
          DO 19 IV = 1, NOFVERT-1
C     check if a periodic node, if not skip the
C     following loop
             IF( PFLAG(IV) )THEN
C     pre-multiply by Qt
                CYY = NODRES(IY,IV)
                CZZ = NODRES(IZ,IV)
                NODRES(IY,IV) =  CYY*COSALPHA+CZZ*SINALPHA
                NODRES(IZ,IV) = -CYY*SINALPHA+CZZ*COSALPHA
             ENDIF
   19     CONTINUE
      ENDIF
