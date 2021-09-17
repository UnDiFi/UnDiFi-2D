      SUBROUTINE SETBNDFLXPTR(BFLX,IFACNOD,IBPTR,NOFVAR,NDIM,NOFVERT,
     &IBNDPTR,NBFAC,ICELNOD,RV,IREAD,JOB)
c
c     $Id: setbndflxptr.f,v 1.3 2020/03/28 09:46:12 abonfi Exp $
c
      IMPLICIT NONE
c
      INCLUDE 'bnd.h'       
      INCLUDE 'bnd.com'       
      INCLUDE 'io.com'       
c
      INTEGER NOFVAR,NDIM,NOFVERT,NBFAC,IREAD,JOB
      INTEGER IFACNOD(NDIM,*),IBPTR(*),ICELNOD(NOFVERT,*),
     &IBNDPTR(3,NBFAC)
      DOUBLE PRECISION BFLX(NOFVAR,*),RV(*)
c
      INTEGER IFACE,I,J,IFAIL,IE,IV,IB,IC
c
      INTEGER IARR(3,2) 
c
      INTEGER ICYCL
c
c     there are NBFLX(2) faces on the current processor where a
c     prescribed flux is enforced
c
c     these faces are defined by their nodes, stored in FACNOD(NDIM,NBFLX(2))
c     this vector is no longer needed upon return
c
c     IBPTR(NBFLX(2)) addresses the corresponding entry in IBNDPTR 
c
      GOTO(1,2)JOB
      WRITE(6,*)'Wrong JOB ',JOB,' in subroutine SETBNDFLXPTR'
      CALL EXIT(128)
    1 CONTINUE
      WRITE(NOUT,*)'Reading ',NBFLX(2),' boundary fluxes'
      DO I = 1,NBFLX(2)
         READ(IREAD,*)(IFACNOD(J,I),J=1,(NOFVERT-1)),IC
         READ(IREAD,*)(BFLX(J,I),J=1,NOFVAR)
      ENDDO
      CLOSE(IREAD)
      CALL I4Mat_Print('General',' ',NOFVERT-1,NBFLX(2),IFACNOD,
     +            NOFVERT-1,'Boundary vertices ',IFAIL)
      CALL R8Mat_Print('General',' ',NOFVAR,NBFLX(2),BFLX,
     +            NOFVAR,'Boundary fluxes ',IFAIL)
      RETURN
    2 CONTINUE
      CALL I4Mat_Print('General',' ',NOFVERT-1,NBFLX(2),IFACNOD,
     +            NOFVERT-1,'Boundary vertices (2) ',IFAIL)
      DO 5 I = 1,NBFLX(2)
         DO J = 1,(NOFVERT-1)
            IARR(J,1) = IFACNOD(J,I)
         ENDDO
         CALL QSORT2(NOFVERT-1,IARR(1,1),IFAIL)
C
C        loop over all boundary faces to find the one that matches
C
         DO IFACE = 1, NBFAC
            IE = IBNDPTR(1,IFACE)
            IV = IBNDPTR(2,IFACE)
            IB = IBNDPTR(3,IFACE)
            DO J = 1,NOFVERT-1
               IARR(J,2) = ICELNOD(ICYCL(IV+J,NOFVERT),IE)
            ENDDO
            CALL QSORT2(NOFVERT-1,IARR(1,2),IFAIL)
            IFAIL = 0
            DO J = 1,NOFVERT-1
               IF(IARR(J,1).NE.IARR(J,2))IFAIL = IFAIL+1
            ENDDO
            IF(IFAIL.EQ.0)THEN
               IBPTR(I) = IFACE 
               GOTO 5
            ENDIF
         ENDDO ! end loop over faces
         WRITE(NOUT,*)(IARR(J,1),J=1,NOFVERT-1)
         WRITE(NOUT,*)(IFACNOD(J,I),J=1,NOFVERT-1)
         WRITE(NOUT,*)'item is ',I
         STOP 'should not arrive here'
    5 CONTINUE ! loop over boundary fluxes
C
C     here we use IFACNOD as a 1-dimensional work array
C     it is 3 times the required length
C
cnag  CALL M01DBF(IBPTR,1,NBFLX(2),'Ascending',IFACNOD,IFAIL)
      CALL QSORTI(IFACNOD,NBFLX(2),IBPTR)
      CALL RNKIDX(IFACNOD,1,NBFLX(2),IFAIL)
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
cnag  CALL M01EBF(IBPTR,1,NBFLX(2),IFACNOD,IFAIL) ! re-arrange the integer pointer
      CALL I4RANK(IBPTR,1,NBFLX(2),IFACNOD,IFAIL) ! re-arrange the integer pointer
      IF(IFAIL.NE.0)CALL EXIT(IFAIL)
      DO I = 1,NOFVAR ! sort fluxes
         CALL DCOPY(NBFLX(2),BFLX(I,1),NOFVAR,RV,1)
cnag     CALL M01EAF(RV,1,NBFLX(2),IFACNOD,IFAIL) ! re-arrange the i-th flux component
         CALL R8RANK(RV,1,NBFLX(2),IFACNOD,IFAIL) ! re-arrange the i-th flux component
         IF(IFAIL.NE.0)CALL EXIT(IFAIL)
         CALL DCOPY(NBFLX(2),RV,1,BFLX(I,1),NOFVAR)
      ENDDO
      WRITE(6,*)(IBPTR(I),I=1,NBFLX(2))
      CALL R8Mat_Print('General',' ',NOFVAR,NBFLX(2),BFLX,
     +            NOFVAR,'Sorted Boundary fluxes ',IFAIL)
      RETURN
      END
