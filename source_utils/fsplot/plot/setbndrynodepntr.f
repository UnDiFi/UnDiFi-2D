      SUBROUTINE SetBndryNodePtr(LBNDFAC,LNODCOD,LNODPTR,LCELNOD,NBFAC,
     &NPOIN,NBPOIN)
C
      IMPLICIT NONE
C
C     $Id:$ 
C
C
C
C     NDIM  is the space dimension =2
C     NVT = NDIM+1 is the number of vertices
C
C
C
C     .. Parameters ..
C     ..
C     .. Local Scalars ..
      INTEGER LBNDFAC,LNODCOD,LCELNOD
      INTEGER LNODPTR,LWKSP
      INTEGER NBFAC,NPOIN,IPOIN,NBPOIN
      INTEGER IFAIL,IER
C     ..
C     .. Local Arrays ..
C
C     ..
C     .. External Subroutines ..
      EXTERNAL DINIT,IINIT,ISTKIN,ISTKRL
C     ..
C
C
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DSTAK(1)
      INTEGER ISTAK(1)
      COMMON/CSTAK/DSTAK
C     ..
C     .. Equivalences ..
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C     ..
C     .. External Functions ..
      INTEGER  ISTKGT
      EXTERNAL ISTKGT
C
C     count the nof boundary vertices
C
      NBPOIN = 0
      DO 1 IPOIN = 0, NPOIN-1
          IF( ISTAK(LNODCOD+IPOIN) .GT. 0 )NBPOIN = NBPOIN+1
    1 CONTINUE
      write(6,*)'Found ',NBPOIN,' boundary points in the nodcod array'
      LNODPTR = ISTKGT(5*NBPOIN,2)
      LWKSP   = ISTKGT(2*NBPOIN,2)
      CALL IINIT(5*NBPOIN,0,ISTAK(LNODPTR),1)
      CALL IINIT(2*NBPOIN,0,ISTAK(LWKSP),1)
      CALL MYROUTINE(ISTAK(LNODCOD),ISTAK(LBNDFAC),NBFAC,ISTAK(LNODPTR),
     &NPOIN,NBPOIN,ISTAK(LCELNOD),ISTAK(LWKSP))
      CALL ISTKRL(1)
      RETURN
      END
C
      SUBROUTINE MYROUTINE(NODCODE,IBNDPTR,NBFAC,INODPTR,NPOIN,NBPOIN,
     &ICELNOD,IWORK)
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NPOIN,NBPOIN
      INTEGER NODCODE(*),ICELNOD(3,*)
      INTEGER IBNDPTR(3,NBFAC),INODPTR(NBPOIN,5),IWORK(2*NBPOIN)
      INTEGER IPOIN,IPOS,LAST,IFAIL,J,K,IFACE,JPOIN,IELEM,IVERT,IER
      INTEGER iarray(2),INODE(2)
      INTEGER  ICYCL,JCYCL
      EXTERNAL ICYCL,JCYCL
C
C     INODPTR is a nodal pointer for boundary nodes
C     INODPTR(*,1) addresses the global nodenumber
C     INODPTR(*,2) addresses one of the two edges it belongs to
C     INODPTR(*,3) addresses the other edge it belongs to
C     INODPTR(*,4) is the other vertex number of meshpoint
C     INODPTR(*,5) is the other
C
C     fill the first entry of the pointer with the node number
C
      LAST = 0
      DO 1 IPOIN = 1, NPOIN
          IF( NODCODE(IPOIN) .GT. 0 )THEN
              LAST = LAST+1
              IWORK(LAST) = IPOIN
          ENDIF 
    1 CONTINUE
      IF(LAST.EQ.NBPOIN)THEN
         write(6,*)'Found ',NBPOIN,' boundary points in array NODCODE()'
      ELSE
         STOP 'LAST .NE. NBPOIN'
      ENDIF 
C
C     IWORK(1:NBPOIN) stores the NBPOIN node numbers
C
      CALL ISORTRX(NBPOIN,IWORK,IWORK(NBPOIN+1))
      DO 2 IPOIN = 1, NBPOIN
         INODPTR(IPOIN,1) = IWORK(IWORK(NBPOIN+IPOIN))
    2 CONTINUE
!     write(6,*)(inodptr(ipoin,1),ipoin=1,nbpoin)
!     pause
C
!     CALL X04EAF('General',' ',NBPOIN,3,INODPTR,NBPOIN,
!    +            'Nodal Bndry pointer',IFAIL)
!     CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!    +            'Face Bndry pointer',IFAIL)
C
      IFAIL = 0
      DO 4 IFACE = 1, NBFAC
         IELEM = IBNDPTR(1,IFACE)
         IVERT = IBNDPTR(2,IFACE)
         DO J = 1,2 ! loop over the nodes of the bndry face
            INODE(J) = ICELNOD(JCYCL(IVERT+J),IELEM)
         ENDDO
         DO 5 J = 1,2 ! loop over the nodes of the bndry face
            IPOIN = INODE(J)
            JPOIN = INODE(ICYCL(J+1,2))
            CALL BINSRC(IPOIN,INODPTR(1,1),NBPOIN,IPOS,LAST)
            IF(IPOS.EQ.0)THEN
            write(6,*)'Subr. CheckBndryPntr: Entry NOT found for node',
     &      IPOIN,' face = ',iface
            CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
     +            'Face Bndry pointer',IFAIL)
            STOP
            ENDIF
            DO 6 K = 2,3
               IF(INODPTR(IPOS,K).EQ.0)THEN
                  INODPTR(IPOS,K)=IFACE
                  INODPTR(IPOS,K+2)=JPOIN
                  GOTO 5
               ENDIF 
    6       CONTINUE
C
C    il nodo sembra appartenere a piu di 2 facce di frontiera
C    l'array e' gia' pieno
C
         IFAIL = IPOIN
         write(6,*)'Il nodo sembra appartenere a piu di 2 facce di front
     &iera'
         write(6,*)'Face no. ',IFACE,(IBNDPTR(K,IFACE),K=1,3)
         write(6,*)'Node no. ',IPOIN,(INODPTR(IPOS,K),K=1,3)
         GOTO 7
    5 CONTINUE ! end loop over the nodes of the bndry face (J)

    4 CONTINUE ! loop over bndry faces
    7 CONTINUE
C
C     Make a simple check: all nodes in INODPTR(*,4:5) should be boundary nodes
C
      DO J = 1, NBPOIN
         DO K = 4,5
            IPOIN = INODPTR(J,K)
            CALL BINSRC(IPOIN,INODPTR(1,1),NBPOIN,IPOS,LAST)
            IF(IPOS.EQ.0)THEN
            write(6,*)'Subr. CheckBndryPntr: CHECK failed for node',
     &      IPOIN
            CALL X04EAF('General',' ',NBPOIN,5,INODPTR,NBPOIN,
     +            'Nodal Bndry pointer',IFAIL)
            CALL EXIT(1)
            ENDIF 
         ENDDO
      ENDDO
C
      IF(IFAIL.NE.0)THEN 
          IER = IFAIL
          CALL X04EAF('General',' ',NBPOIN,5,INODPTR,NBPOIN,
     +            'Nodal Bndry pointer',IFAIL)
          CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
     +            'Face Bndry pointer',IFAIL)
          WRITE(6,FMT=*)'Unrecoverable error in CheckBndryPntr'
          CALL EXIT(IER) 
      ENDIF
C
      CALL X04EAF('General',' ',NBPOIN,5,INODPTR,NBPOIN,
     +            'Nodal Bndry pointer',IFAIL)
!     CALL X04EAF('General',' ',3,NBFAC,IBNDPTR,3,
!    +            'Face Bndry pointer',IFAIL)
      return
      do iface = 1, nbfac
         ielem = ibndptr(1,iface)
         ivert = ibndptr(2,iface)
         do j = 1,2
            iarray(j) = icelnod(jcycl(ivert+j),ielem)
         enddo
         write(6,*)iface,ielem,(iarray(j),j=1,2),icelnod(ivert,ielem)
      enddo
!     stop 
      RETURN
      END

