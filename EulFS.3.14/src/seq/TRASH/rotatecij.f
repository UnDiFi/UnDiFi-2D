      IF( PERIODIC_MESH .AND. ANNULAR )THEN
          DO 9 IVERT = 1, NOFVERT
C     check if a periodic node, if not skip the
C     following loop
             IF( PFLAG(IVERT) )THEN
C     pre-multiply C_il by Qt
             DO 8 JVERT = 1, NOFVERT
                 CYY = STIFEL(IY,IY,IVERT,JVERT)
                 CYZ = STIFEL(IY,IZ,IVERT,JVERT)
                 CZY = STIFEL(IZ,IY,IVERT,JVERT)
                 CZZ = STIFEL(IZ,IZ,IVERT,JVERT)
                 STIFEL(IY,IY,IVERT,JVERT) = CYY*COSALPHA-CZY*SINALPHA
                 STIFEL(IY,IZ,IVERT,JVERT) = CYZ*COSALPHA-CZZ*SINALPHA
                 STIFEL(IZ,IY,IVERT,JVERT) = CYY*SINALPHA+CZY*COSALPHA
                 STIFEL(IZ,IZ,IVERT,JVERT) = CYZ*SINALPHA+CZZ*COSALPHA
    8        CONTINUE
C     post-multiply C_ii by Q
             CYY = STIFEL(IY,IY,IVERT,IVERT)
             CYZ = STIFEL(IY,IZ,IVERT,IVERT)
             CZY = STIFEL(IZ,IY,IVERT,IVERT)
             CZZ = STIFEL(IZ,IZ,IVERT,IVERT)
             STIFEL(IY,IY,IVERT,IVERT) = CYY*COSALPHA-CYZ*SINALPHA
             STIFEL(IY,IZ,IVERT,IVERT) = CYY*SINALPHA+CYZ*COSALPHA
             STIFEL(IZ,IY,IVERT,IVERT) = CYZ*COSALPHA-CZZ*SINALPHA
             STIFEL(IZ,IZ,IVERT,IVERT) = CZY*SINALPHA+CZZ*COSALPHA
!            write(6,*)ielem
!            write(6,*)(icn(k)+1,k=1,nofvert)
!            write(6,*)(icelnod(k,nelem+ielem),k=1,nofvert)
caldo        pause
             ENDIF
    9     CONTINUE
      ENDIF
