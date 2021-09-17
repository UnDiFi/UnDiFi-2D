      SUBROUTINE PROBEIN(IELEM,WEIGHTS,NITEMS,NOFVERT,IUNIT)
C
C     $Id: probe.f,v 1.1 2011/03/30 09:16:20 abonfi Exp $
C
      IMPLICIT NONE
      INTEGER NOFVERT,NITEMS,IUNIT
      INTEGER IELEM(NITEMS)
      DOUBLE PRECISION WEIGHTS(NOFVERT,NITEMS)
      INTEGER I,J
      DO I = 1, NITEMS
         READ(IUNIT,*) IELEM(I),(WEIGHTS(J,I),J=1,NOFVERT)
!        WRITE(6,*) IELEM(I),(WEIGHTS(J,I),J=1,NOFVERT)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PROBEOUT(IELEM,WEIGHTS,NITEMS,
     &ZROE,NDOF,ICELNOD,NOFVERT,TIME,IUNIT)
C
C     $Id: probe.f,v 1.1 2011/03/30 09:16:20 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INTEGER NOFVERT,NITEMS,NDOF,IUNIT
      DOUBLE PRECISION TIME
      INTEGER IELEM(NITEMS),ICELNOD(NOFVERT,*)
      DOUBLE PRECISION WEIGHTS(NOFVERT,NITEMS),ZROE(NDOF,*)
C
      DOUBLE PRECISION ZPROBE(MAXNOFVAR),WF
      INTEGER I,J,K,L,IPOIN
      DO I = 1, NITEMS ! loop over probes
         DO J = 1,NDOF ! initialize interpolated value
            ZPROBE(J) = ZERO 
         ENDDO
         K = IELEM(I) ! pick up the element the i-th probe belongs to
         DO L = 1,NOFVERT  ! loop over vertices of the element
            IPOIN = ICELNOD(L,K) ! pick up the gridpoint
            WF = WEIGHTS(L,I) ! pick up the weight
            DO J = 1,NDOF
               ZPROBE(J) = ZPROBE(J) + ZROE(J,IPOIN) * WF
            ENDDO
         ENDDO ! loop over vertices
         WRITE(IUNIT,*) I,TIME,(ZPROBE(J),J=1,NDOF)
      ENDDO
      RETURN
      END
