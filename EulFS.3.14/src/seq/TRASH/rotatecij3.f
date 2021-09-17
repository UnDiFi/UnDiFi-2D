      IF( PERIODIC_MESH .AND. ANNULAR )THEN
          DO 39 IV = 1, NOFVERT-1
C     check if a periodic node, if not skip the
C     following loop
             IF( PFLAG(IV) )THEN
C     pre-multiply C_il by Qt
             DO 38 JV = 1, 1
                 CYY = STIFEL(IY,IY,IV,JV)
                 CYZ = STIFEL(IY,IZ,IV,JV)
                 CZY = STIFEL(IZ,IY,IV,JV)
                 CZZ = STIFEL(IZ,IZ,IV,JV)
                 STIFEL(IY,IY,IV,JV) = CYY*COSALPHA-CZY*SINALPHA
                 STIFEL(IY,IZ,IV,JV) = CYZ*COSALPHA-CZZ*SINALPHA
                 STIFEL(IZ,IY,IV,JV) = CYY*SINALPHA+CZY*COSALPHA
                 STIFEL(IZ,IZ,IV,JV) = CYZ*SINALPHA+CZZ*COSALPHA
   38        CONTINUE
C     post-multiply C_ii by Q
             CYY = STIFEL(IY,IY,IV,IV)
             CYZ = STIFEL(IY,IZ,IV,IV)
             CZY = STIFEL(IZ,IY,IV,IV)
             CZZ = STIFEL(IZ,IZ,IV,IV)
             STIFEL(IY,IY,IV,IV) = CYY*COSALPHA-CYZ*SINALPHA
             STIFEL(IY,IZ,IV,IV) = CYY*SINALPHA+CYZ*COSALPHA
             STIFEL(IZ,IY,IV,IV) = CYZ*COSALPHA-CZZ*SINALPHA
             STIFEL(IZ,IZ,IV,IV) = CZY*SINALPHA+CZZ*COSALPHA
!            write(6,*)ielem,npoin+nghost
!            write(6,*)(icn(i)+1,i=1,nofvert)
!            write(6,*)(icelnod(i,ielem),i=1,nofvert)
!            write(6,*)(icelnod(i,nelem+ielem),i=1,nofvert)
!            write(6,*)(pflag(i),i=1,nofvert)
!            pause
             ENDIF
   39     CONTINUE
      ENDIF
