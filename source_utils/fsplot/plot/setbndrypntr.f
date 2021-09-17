      SUBROUTINE SetBndryPntr(ICELNOD,ICELCEL,NOFVERT,NELEM,IBNDFAC,
     &NBFAC,IBNDPTR,INODPTR,NBPOIN)
C
      IMPLICIT NONE
C
C     $Id:$ 
C
C
C
C     .. Parameters ..
C     ..
C     .. Local Scalars ..
      INTEGER NOFVERT,NELEM,NBPOIN,NBFAC
      INTEGER IFAIL,IELEM,IVERT,I,J,IBFAC,IPOIN,NL,NR,LAST,IPOS,IEDGE
C     ..
C     .. Local Arrays ..
C
      INTEGER ICELNOD(NOFVERT,NELEM),INODPTR(NBPOIN,5)
      INTEGER ICELCEL(NOFVERT,NELEM),IBNDFAC(3,NBFAC),IBNDPTR(6,NBFAC)
C     ..
C     .. External Subroutines ..
C     ..
C
C
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Equivalences ..
C     ..
C     .. External Functions ..
C
      INTEGER JCYCL,ICYCL
C
C
C     this routine fills the entries of IBNDPTR
C
C        IEL          IE         IER
C   +-----------+-----------+-----------+----
C  NNL          NL          NR          NNR 
C
C
C        IBNDPTR(1,IE) = NL
C        IBNDPTR(2,IE) = NR
C        IBNDPTR(3,IE) = NNL
C        IBNDPTR(4,IE) = NNR
C        IBNDPTR(5,IE) = IEL
C        IBNDPTR(6,IE) = IER
C
C
      DO IBFAC = 1,NBFAC
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
         DO J = 1,2 ! the nodes that make up the edge
            IBNDPTR(J,IBFAC) = ICELNOD(JCYCL(IVERT+J),IELEM)
         ENDDO
C
         DO J = 1,2 ! pick up the nodes once at a time
            NL = IBNDPTR(J,IBFAC)
            NR = IBNDPTR(ICYCL(J+1,2),IBFAC)
C           find node
            CALL BINSRC(NL,INODPTR(1,1),NBPOIN,IPOS,LAST)
            IF(IPOS.EQ.0)THEN
                WRITE(6,*)'Subr. SetBndryPntr: Entry NOT found for node'
     &, IPOIN,' face = ',IBFAC
                CALL X04EAF('General',' ',3,NBFAC,IBNDFAC,3,
     +            'Face Bndry pointer',IFAIL)
                CALL EXIT(LAST)
            ENDIF
            DO I = 2,3 ! loop over the edges on both side of node NL
               IEDGE = INODPTR(IPOS,I) ! the edge on one side of NL
               IPOIN = INODPTR(IPOS,I+2) ! the node on one side of NL
               IF( IEDGE .EQ. IBFAC )then
                    IF( IPOIN .NE. NR )then
                        write(6,*)'test #1 failed ',IPOIN,NR
                        call exit(IPOIN-NR) 
                    endif
               ELSE
                    IBNDPTR(J+2,IBFAC) = IPOIN
                    IBNDPTR(J+4,IBFAC) = IEDGE
               ENDIF
               ENDDO ! end loop over the edges on both side of node NL
         ENDDO ! end loop on nodes that make up the bndry edge
      ENDDO ! end loop over bndry faces
      CALL X04EAF('General',' ',6,NBFAC,IBNDPTR,6,
     +            'Bndry edge pointer',IFAIL)
      stop
      RETURN
      END
