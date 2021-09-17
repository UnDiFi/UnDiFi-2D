      SUBROUTINE BURKARDT(NDOF,NDIM,NPOIN,NPNOD,NOFVAR,NOFVERT,
     &NELEM,ICELNOD,ICELCEL,XYBKGR,XY,ZBKGR,ZOUT,ORD)
C
C     Purpose: interpolate from the 2D grid 1 onto the 3D grid 2
C     using the triangulation library
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER NDOF
      INTEGER NDIM,NELEM(*),NOFVAR,NPOIN(*),NOFVERT,NPNOD(*)
      INTEGER ICELNOD(NOFVERT,*),ICELCEL(NOFVERT,*),ORD(*)
      LOGICAL LFLAG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION XYBKGR(NDIM,*),XY(NDIM,*),
     &                  ZBKGR(NDOF,*),ZOUT(NDOF,*)
C     ..
C     .. Arrays in Common ..
C     ..
C     .. Local Scalars ..
      INTEGER EDGE,IELEM,N,IPOIN,I,IFAIL,NITEMS,INFO,NOTFOUND
      INTEGER IVAR,IFREQ,J,NAIVEFAIL,CLEVERFAIL
      DOUBLE PRECISION P(2),WGHT(3)
C     ..
C     .. Local Arrays ..
      INTEGER SEED(3)
      COMMON/RANDOM/SEED
C     ..
C     .. External Functions ..
      INTEGER ICYCL
      DOUBLE PRECISION AREA,DNRM2,DASUM
      EXTERNAL ICYCL,AREA,DNRM2,DASUM
C     ..
C     .. Intrinsic Functions ..
C     INTRINSIC SQRT
C     ..
C     .. External Functions ..
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DINIT,IINIT,ISTKIN
C     ..
C     .. Common blocks ..
C     ..
C     .. Equivalences ..
C     ..
      WRITE(6,*)'Background mesh info'
      WRITE(6,*)'has ',NPOIN(2)+NPNOD(2),' gridpoints'
      WRITE(6,*)'has ',NOFVAR,' degrees of freedom'
      DO I = 1,NOFVAR
         WRITE(6,*)I,DNRM2(NPOIN(2)+NPNOD(2),ZBKGR(I,1),NOFVAR),
     &DASUM(NPOIN(2)+NPNOD(2),ZBKGR(I,1),NOFVAR)/
     &REAL(NPOIN(2)+NPNOD(2))
      ENDDO
      WRITE(6,*)'New mesh info'
      WRITE(6,*)'has ',NPOIN(3)+NPNOD(3),' gridpoints'
      WRITE(6,*)'has ',NOFVAR,' degrees of freedom'
      DO I = 1,NOFVAR
         WRITE(6,*)I,DNRM2(NPOIN(3)+NPNOD(3),ZOUT(I,1),NOFVAR),
     &DASUM(NPOIN(3),ZOUT(I,1),NOFVAR)/REAL(NPOIN(3)+NPNOD(3))
      ENDDO
C     ..
C
      LFLAG = .TRUE.
      SEED(1) = 33
      DO IELEM = 1,NELEM(2)
         DO J = 1,NOFVERT
            N = ICELCEL(J,IELEM) 
            IF( (N .LT. 1) .OR. (N .GT. NELEM(2)) )THEN
                N = -(3*IELEM + J-1)
                ICELCEL(J,IELEM) = N
            ENDIF
         ENDDO
      ENDDO
C
      CALL triangulation_order3_check ( NPOIN(2)+NPNOD(2), NELEM(2),
     &  ICELNOD, IFAIL )
      IF(IFAIL.EQ.0) WRITE(6,*)'triangulation_order3_check passed'
      CALL triangulation_order3_neighbor_triangles (
     &  NELEM(2), ICELNOD, ICELCEL )
      WRITE(6,*)'triangulation_order3_neighbor_triangles passed'
C
      CALL X04CAF('General',' ',NDIM,NPOIN(3)+NPNOD(3),XY,
     +            NDIM,'Nodal coordinates of the probes',IFAIL)
C
C
      WRITE(6,*)" Now interpolating; log file is trpack.log"
      OPEN(12,FILE="trpack.log",FORM="formatted",STATUS="unknown")
      CLEVERFAIL = 0
      NAIVEFAIL = 0
      NITEMS = NPOIN(3)+NPNOD(3)
      IFREQ = MAX0( NITEMS/100 , 1 )
      WRITE(6,*)'Interpolating values in ',NITEMS,' gridpoints'
      DO 93 IPOIN = 1,NITEMS
         IF((IPOIN/IFREQ)*IFREQ .EQ. IPOIN)WRITE(*,111)IPOIN,
     &      IPOIN/REAL(NITEMS)*100,CLEVERFAIL,
     &      CLEVERFAIL/REAL(IPOIN)*100,NAIVEFAIL
         P(1) = XY(1,IPOIN) 
         P(2) = XY(2,IPOIN) 
         CALL triangulation_order3_search ( NPOIN(2)+NPNOD(2) , XYBKGR ,
     &        NOFVERT, NELEM(2), ICELNOD, ICELCEL,
     &        p, seed, IELEM , edge )
C        IELEM = -1 
         IF( IELEM .EQ. -1 .OR. EDGE .LT. 0 )THEN
C
C           the clever search has failed; increment IFAIL
C
            CLEVERFAIL = CLEVERFAIL+1
            CALL triangulation_order3_search_naive ( 
     3           NPOIN(2)+NPNOD(2) , XYBKGR ,  NOFVERT, NELEM(2), 
     3           ICELNOD, p,  IELEM )
                 IF(IELEM.GT.NELEM(2).OR.IELEM.LT.1)THEN
C
C           also the naive search has failed
C
                   NAIVEFAIL = NAIVEFAIL + 1 
                   LFLAG = .FALSE.
!                  WRITE(6,*)'Naive search **FAILED** '
                   WRITE(12,*)'Using inverse distance at (x,y) ',
     &          (P(i),i=1,2)
                   call inversed(ord,XYBKGR,NDIM,npoin(2),ZBKGR,NDOF,
     &                           p,ZOUT(1,IPOIN),info)
                 ELSE
                   LFLAG = .TRUE.
!                  WRITE(6,*)'Naive search has ret ',ipoin,ielem
                 ENDIF
         ELSE
!                  WRITE(6,*)'Clever search has ret ',ipoin,ielem
                   LFLAG = .TRUE.
         ENDIF
C
         IF(LFLAG)THEN ! the inverse distance has already interpolated
C
                call interp(icelnod,xybkgr,NDIM,ZBKGR,NDOF,P,
     &                      ZOUT(1,IPOIN),ielem,info,wght)
                IF( info .NE. 0 )THEN
!                   WRITE(6,*)'Error locating ',(P(i),i=1,2),ielem,
!    +                        (wght(i),i=1,3)
                    WRITE(12,*)'Error locating ',(P(i),i=1,2),ielem,
     +                        (wght(i),i=1,3)
                    CALL EXIT(1)
                ENDIF
         ENDIF
   93 CONTINUE ! loop over nodes
C
      WRITE(6,*)'Clever search failed in ',CLEVERFAIL,' nodes' ,
     &(CLEVERFAIL/REAL(NITEMS))*100.
      WRITE(6,*)NAIVEFAIL,' nodes' ,
     &(NAIVEFAIL/REAL(NITEMS))*100.,' where extrapolation occurs'
C
      DO I = 1,NOFVAR
         WRITE(6,*)I,DNRM2(NPOIN(3)+NPNOD(3),ZOUT(I,1),NOFVAR),
     &DASUM(NPOIN(3)+NPNOD(3),ZOUT(I,1),NOFVAR)/
     3REAL(NPOIN(3)+NPNOD(3))
      ENDDO
 111  FORMAT(I7,2X,F6.2,'%',1X,I7,2X,F6.2,'%',1X,I6)
 1    FORMAT(A)
      RETURN
      END
      SUBROUTINE RPROBE(IUNIT,XY,NDIM,NITEMS)
      IMPLICIT NONE
      INTEGER IUNIT,NDIM,NITEMS
      DOUBLE PRECISION XY(NDIM,NITEMS)
      INTEGER I,J
      DO I = 1,NITEMS
         READ(IUNIT,*)(XY(J,I),J=1,NDIM)
      ENDDO
      END
      SUBROUTINE WPROBE(IUNIT,XY,NDIM,NITEMS,Z,NDOF)
      IMPLICIT NONE
      INTEGER IUNIT,NDIM,NITEMS,NDOF
      DOUBLE PRECISION XY(NDIM,NITEMS),Z(NDOF,NITEMS)
      INTEGER I,J
      DO I = 1,NITEMS
         WRITE(IUNIT,*)I,(XY(J,I),J=1,NDIM),
     &       (Z(J,I),J=1,NDOF)
      ENDDO
      END
