      SUBROUTINE RTRI(XY,NPOIN,ZROE,NDOF,ICELNOD,ICELCEL,IBNDPTR,NELEM,
     &NBFAC,FNAME,LENFNAM)

      IMPLICIT NONE

      INTEGER NPOIN,NELEM,NBFAC,NDOF
      CHARACTER*(*) FNAME(*)
      DOUBLE PRECISION XY(2,*),ZROE(NDOF,NPOIN)
      INTEGER ICELNOD(3,*),ICELCEL(3,*),IBNDPTR(3,*),LENFNAM(*)
      INTEGER I,J,K,IC,ID,IELEM,IVERT,IBC,IFACE,N1,N2,IN1,IN2
      INTEGER IBFAC,IPOIN,NFACE
      INTEGER  LENSTR,JCYCL
      EXTERNAL LENSTR,JCYCL
      INTEGER ICLR(0:10)
C
      DO K = 0, 10
         ICLR(K) = 0
      ENDDO
C
      k = LENFNAM(1)
      OPEN(17,FILE=FNAME(1)(1:K))
      READ(17,*)NPOIN
      write(6,*)'Opening and reading ',npoin,' meshpoints from ',
     &fname(1)(1:K)
      DO 1 IPOIN = 1, NPOIN
         READ(17,*)I,XY(1,IPOIN),XY(2,IPOIN),
     &   (ZROE(K,IPOIN),K=1,NDOF),ID
    1 CONTINUE
      CLOSE(17)
      write(6,*)'Done ! '

      k = LENFNAM(2)
      OPEN(17,FILE=FNAME(2)(1:K))
      READ(17,*)NELEM
      write(6,*)'Opening and reading ',nelem,' triangles from ',
     &fname(2)(1:LENFNAM(2))
      DO 3 IELEM = 1, NELEM
         READ(17,*)I,(ICELNOD(J,IELEM),J=1,3)
    3 CONTINUE
      CLOSE(17)
      write(6,*)'Done ! '

!     return

      k = LENFNAM(3)
      write(6,*)'Opening and reading ',nelem,' neighbours from ',
     &fname(3)(1:K)
      OPEN(17,FILE=FNAME(3)(1:K))
      READ(17,*)
      DO 5 IELEM = 1, NELEM
         READ(17,*)I,(ICELCEL(J,IELEM),J=1,3)
    5 CONTINUE
      CLOSE(17)
      write(6,*)'Done ! '

      IBFAC = 0
      DO 7 IELEM = 1, NELEM
         DO 7 I = 1,3
         IF( ICELCEL(I,IELEM) .EQ. -1 )THEN
            IBFAC = IBFAC+1
            ICELCEL(I,IELEM) = 0
            IN1 = ICELNOD(JCYCL(I+1),IELEM)
            IN2 = ICELNOD(JCYCL(I+2),IELEM)
            write(17,*)IN1,IN2
            IBNDPTR(1,IBFAC) = IELEM
            IBNDPTR(2,IBFAC) = I
! still undefined
            IBNDPTR(3,IBFAC) = -1
         ENDIF
    7 CONTINUE
      write(6,*)'There appear to be ',IBFAC,' bndry faces '
      NBFAC = IBFAC

      k = LENFNAM(4)
      OPEN(17,FILE=fname(4)(1:K))
      READ(17,*)NFACE,IC
      write(6,*)'Opening and reading ',nface,' faces from ',
     &fname(4)(1:K)
      IF( IC.NE.1 )THEN
          write(6,*)'There should be 1 bndry marker in ',fname(4)(1:K),
     &    'while there appear to be ',IC,' run triangle with -e'
          CALL EXIT(1)
      ENDIF
! scan the edges
      IBFAC = 0
      DO 9 IFACE = 1, NFACE
         READ(17,*)I,N1,N2,IBC
! if a bndry face, look for the parent element in ibndptr
         IF(IBC .EQ. 0 )THEN
              ICLR(IBC) = ICLR(IBC)+1
         ELSE
              IBFAC = IBFAC + 1
              DO 11 I = 1,NBFAC
                 IELEM = IBNDPTR(1,I)
                 IVERT = IBNDPTR(2,I)
                 IN1 = ICELNOD(JCYCL(IVERT+1),IELEM)
                 IN2 = ICELNOD(JCYCL(IVERT+2),IELEM)
                 IF( (IN1 .EQ. N1 .AND. IN2 .EQ. N2 ) .OR.
     &               (IN1 .EQ. N2 .AND. IN2 .EQ. N1 ) )THEN
                     IF( IBNDPTR(3,I) .NE. -1 )then
                 WRITE(6,*)'We have a problem !! ',I,IBNDPTR(3,I)
                 CALL EXIT(1)
                     ENDIF
                     IBNDPTR(3,I) = IBC
                     ICLR(IBC) = ICLR(IBC)+1
                     GOTO 9
                 ENDIF
   11         CONTINUE
              WRITE(6,*)'Cannot match boundary face ',iface,n1,n2,
     &'within subroutine rtri'
!             CALL EXIT(1)
         ENDIF
    9 CONTINUE
      CLOSE(17)
       
      WRITE(6,*)'Boundary faces are ',IBFAC,NBFAC
      DO K = 0, 10
         IF(ICLR(K).NE.0)WRITE(6,*)'Bndry edges coloured ',K,
     &' are ',ICLR(K)
      ENDDO
      IF(IBFAC.NE.NBFAC)THEN
         WRITE(6,*) 'Non matching boundary faces'
         CALL EXIT(1)
      ENDIF
      RETURN
      END
