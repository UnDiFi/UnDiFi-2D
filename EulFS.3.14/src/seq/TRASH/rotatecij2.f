      IF( PERIODIC_MESH .AND. ANNULAR )THEN
          DO 9 IV = 1, NOFVERT
C     check if a periodic node, if not skip the
C     following loop
             IF( PFLAG(IV) )THEN
C     pre-multiply C_il by Qt
             DO 8 JV = 1, NOFVERT
                 CYY = STIFEL(IY,IY,IV,JV)
                 CYZ = STIFEL(IY,IZ,IV,JV)
                 CZY = STIFEL(IZ,IY,IV,JV)
                 CZZ = STIFEL(IZ,IZ,IV,JV)
                 STIFEL(IY,IY,IV,JV) = CYY*COSALPHA-CZY*SINALPHA
                 STIFEL(IY,IZ,IV,JV) = CYZ*COSALPHA-CZZ*SINALPHA
                 STIFEL(IZ,IY,IV,JV) = CYY*SINALPHA+CZY*COSALPHA
                 STIFEL(IZ,IZ,IV,JV) = CYZ*SINALPHA+CZZ*COSALPHA
    8        CONTINUE
C     post-multiply C_ii by Q
             CYY = STIFEL(IY,IY,IV,IV)
             CYZ = STIFEL(IY,IZ,IV,IV)
             CZY = STIFEL(IZ,IY,IV,IV)
             CZZ = STIFEL(IZ,IZ,IV,IV)
             STIFEL(IY,IY,IV,IV) = CYY*COSALPHA-CYZ*SINALPHA
             STIFEL(IY,IZ,IV,IV) = CYY*SINALPHA+CYZ*COSALPHA
             STIFEL(IZ,IY,IV,IV) = CYZ*COSALPHA-CZZ*SINALPHA
             STIFEL(IZ,IZ,IV,IV) = CZY*SINALPHA+CZZ*COSALPHA
!            write(6,*)ielem
!            write(6,*)(icn(k)+1,k=1,nofvert)
!            write(6,*)(icelnod(k,nelem+ielem),k=1,nofvert)
             ENDIF
caldo        pause
    9     CONTINUE
      ENDIF
