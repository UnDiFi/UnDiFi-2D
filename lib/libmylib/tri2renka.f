      SUBROUTINE TRI2RENKA(NPOIN,NELEM,NOFVERT,ICELNOD,ICELCEL,JA,IA,
     +                     NODCOD,LIST,LPTR,LEND,LNEW,COOR,NDIM,IERR)
C
      IMPLICIT NONE
C
C     This routine returns the linKed list of the elements each node
C     is surrounded by
C     IA(K) points to the position of JA() where the list of node K starts
C
C
C     .. Scalar Arguments ..
      INTEGER IERR,LNEW,NELEM,NOFVERT,NPOIN,NDIM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COOR(NDIM,*)
      INTEGER IA(*),ICELNOD(NOFVERT,*),JA(*),ICELCEL(NOFVERT,*),
     &NODCOD(*),LIST(*),LPTR(*),LEND(NPOIN)
!    &,IAO(*),JAO(*)
C     ..
C     .. Local Scalars ..
      INTEGER IELEM,IC,J,JBGN,JEND,K,KNOD,KSAV,KSAVN,NEXT,
     &IPOIN,NITEMS,NN,NEIGHB,N1,N2,IFAIL
      INTEGER IWRK(10)
C     ..
C     .. External Subroutines ..
      LOGICAL LEFT
      INTEGER ICYCL
      EXTERNAL IINIT,ICYCL,LEFT
C     ..
C
C     Count the number of elements surrounding each node
C
      LNEW = 0
      DO 3 IPOIN = 1,NPOIN
          JBGN = IA(IPOIN)
          JEND = IA(IPOIN+1)-1
          NITEMS = JEND-JBGN+1 ! nof triangles surrounding gridpoint IPOIN
          NN = NITEMS  ! NN is the number of gridpoints neighbouring gridpoint IPOIN
          IF(NODCOD(IPOIN).NE.0)NN=NN+1 ! add 1 for boundary nodes
!         IF(NODCOD(IPOIN).NE.0)GOTO 3
          KSAV = 0
          DO 10 J = JBGN,JEND ! loop over neighbouring triangles
              IELEM = JA(J)
C         K is J mapped into [1,NITEMS] 
              K = J-JBGN+1
              K = JBGN-1+ICYCL(K+1,NITEMS) ! K:=J+1
              NEXT = JA(K) ! this is the next triangle in the list (cyclically)
!             WRITE(6,*)'IP = ',IPOIN,' IE = ',IELEM,' NEXT = ',NEXT
              DO 8 K = 1, NOFVERT
                 IF( ICELCEL(K,IELEM) .EQ. NEXT )THEN
                     KSAV = KSAV + 1
                     LNEW = LNEW+1
                     LIST(LNEW) = ICELNOD(K,IELEM)
                     GOTO 10
                 ENDIF
    8         CONTINUE
C
C     If we arrive here: either IPOIN is a boundary node or there is an error
C
C
              IF(NODCOD(IPOIN).NE.0)THEN ! boundary nodes require special treatment
!                IF(NITEMS.EQ.1)THEN ! bndry node belongs to only one triangle
                    DO K = 1, NOFVERT
                       NEIGHB = ICELNOD(K,IELEM)
                       IF( NEIGHB .NE. IPOIN )THEN
                           LNEW = LNEW+1
                           LIST(LNEW) = NEIGHB
                       ENDIF
                    ENDDO 
!                ELSE
!                   STOP
!                ENDIF
!                WRITE(6,*)IPOIN,' is a bndry node with ',NITEMS,' neigh
!    &bouring triangles'
!                WRITE(6,*)'found ',KSAV,' vertex neighbours out of ',NN
!                WRITE(6,*)IPOIN,' IE = ',IELEM,' NEXT = ',NEXT
!                WRITE(6,*)'ICN = ',(ICELNOD(K,IELEM),K=1,NOFVERT)
!                WRITE(6,*)'ICC = ',(ICELCEL(K,IELEM),K=1,NOFVERT)
!                WRITE(6,*)'JA = ',(JA(K),K=JBGN,JEND)
!                WRITE(6,*)'Should not arrive here'
              ELSE ! not a bndry node: then there is an error
                 WRITE(6,*)IPOIN,' IE = ',IELEM,' NEXT = ',NEXT
                 WRITE(6,*)'ICN = ',(ICELNOD(K,IELEM),K=1,NOFVERT)
                 WRITE(6,*)'ICC = ',(ICELCEL(K,IELEM),K=1,NOFVERT)
                 WRITE(6,*)'JA = ',(JA(K),K=JBGN,JEND)
                 WRITE(6,*)'Should not arrive here'
                 STOP
              ENDIF
   10     CONTINUE ! end loop over triangles neighbouring gridpoint IPOIN
          LEND(IPOIN) = LNEW
          IFAIL = 0 
          IF(IPOIN.EQ.1)THEN
             JBGN = 1
          ELSE
             JBGN = LEND(IPOIN-1)
          ENDIF 
          JEND = LEND(IPOIN)
          DO J = JBGN,JEND
             n1 = LIST(J)
             NITEMS = JEND-JBGN+1 ! nof gridpoints surrounding gridpoint IPOIN
             K = J-JBGN+1
             K = JBGN-1+ICYCL(K+1,NITEMS) ! K:=J+1
             n2 = LIST(K)
             IF( .NOT. LEFT(COOR(1,N1),COOR(2,N1),COOR(1,N1),COOR(2,N1),
     &                      COOR(1,N1),COOR(2,N1)) )IFAIL = IFAIL+1
          ENDDO
          IF(IFAIL.NE.0)THEN
             WRITE(6,*)'LIST = ',(LIST(K),K=JBGN,JEND)
          ELSE
             WRITE(6,*)'IPOIN = ',IPOIN,NITEMS,(LIST(K),K=JBGN,JEND)
          ENDIF
C
C Here we need to set up LPTR
C
!         WRITE(6,*)'Success!'
!         WRITE(6,*)IPOIN
!         WRITE(6,*)'JA = ',(JA(K),K=JBGN,JEND)
!         WRITE(6,*)NITEMS,IC,(IWRK(K),K=1,IC)
!         LNEW = LNEW + NN
    3 CONTINUE ! loop over gridpoints
      LNEW = LNEW+1
      WRITE(6,*)LNEW,6*NPOIN-12
      RETURN
C
      IERR = 0
      RETURN

      END
