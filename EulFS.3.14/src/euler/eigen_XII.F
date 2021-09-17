      SUBROUTINE Eigen_XII(Matrix,LDA,DVDZ,DUDV,NDIM,LDB)
C
C     $Id: eigen_XII.F,v 1.4 2020/03/28 09:51:15 abonfi Exp $
C
C     Conserved variables 
C
C     define FSPL_USE_DEBUG for debugging purposes, i.e. when using
C     -check
C
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h'
C
      INCLUDE 'dofs.com'
      INCLUDE 'three.com'
      INCLUDE 'four.com'
      INCLUDE 'flags.com'
C
C
      INTEGER LDA,NDIM,LDB
C
      DOUBLE PRECISION Matrix(LDA,LDA,*),DVDZ(LDB,LDB),
     +DUDV(LDB,LDB)
#ifdef FSPL_USE_DEBUG
      DOUBLE PRECISION WORK(NMAX*NMAX)
      LOGICAL NULLMAT
#endif
C
      DOUBLE PRECISION DENSABARINV,DENSINVABARINV
      INTEGER IDIM,ifail,I,J
      DOUBLE PRECISION TMPA(NSP)
C
      DOUBLE PRECISION DDOT 
crpepe      ,PIR,PIE,PIM
      EXTERNAL         DDOT 
crpepe      ,PIR,PIE,PIM 
C
C
C     Conserved variables for a plasma ..
C
      IF(ICHECK.NE.0)THEN 
C
C     N.B. GRAD_CHAR  non serve, tranne quando ICHECK <> 0
C          cosi` la matrice, tranne che si utilizzi la
C          decomposizione numerica. 
C
C     here GRAD_CHAR are the conserved variables
C
      CALL PARM2CONS4Ar(ZAVG,DVDZ,LDB,NDIM)
      CALL DGEMM('NoTranspose','NoTranspose',LDB,NDIM,LDB,
     +           ONE,DVDZ,LDB,GRAD_PARM,NMAX,ZERO,GRAD_CHAR,LDW)
CREM  CALL R8Mat_Print('G',' ',LDB,NDIM,GRAD_CHAR,LDW,
CREM +      'CONSV variables within eigenXI ',IFAIL)
CREM  pause
C
C
C     Jacobian matrix in conserved variables
C  
      DO 15 I = 1 , NSP
        TMPA(I) = DR(I) - UAVG(IE)    
   15 CONTINUE
C
      IDIM = 1
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            Matrix(I,J,IDIM) = - ALPHA(I) * UAVG(IE+IDIM) 
            IF (I.EQ.J)THEN
                Matrix(I,J,IDIM) = Matrix(I,J,IDIM) + UAVG(IE+IDIM)
            ENDIF
        ENDDO
      ENDDO
C
      DO I = 1 , NSP
        Matrix(I,IE,IDIM) = ZERO
        Matrix(I,IX,IDIM) = ALPHA(I)
        Matrix(I,IY,IDIM) = ZERO
        Matrix(I,IZ,IDIM) = ZERO
      ENDDO
C
      DO I = 1 , NSP
        Matrix(IE,I,IDIM) = UAVG(IE+IDIM) * TMPA(I)
      ENDDO
      Matrix(IE,IE,IDIM) = UAVG(IX) * (ONE + DE)
      Matrix(IE,IX,IDIM) = UAVG(IE) + DM(1)*UAVG(IX)
      Matrix(IE,IY,IDIM) = DM(2) * UAVG(IX)
      Matrix(IE,IZ,IDIM) = DM(3) * UAVG(IX)
C
      DO I = 1 , NSP
        Matrix(IX,I,IDIM) = -UAVG(IX)*UAVG(IX) + DR(I)
      ENDDO
      Matrix(IX,IE,IDIM) = DE
      Matrix(IX,IX,IDIM) = DM(1) + TWO*UAVG(IX)
      Matrix(IX,IY,IDIM) = DM(2)
      Matrix(IX,IZ,IDIM) = DM(3)
C
      DO I = 1 , NSP
        Matrix(IY,I,IDIM) = -UAVG(IX)*UAVG(IY)
      ENDDO
      Matrix(IY,IE,IDIM) = ZERO
      Matrix(IY,IX,IDIM) = UAVG(IY)
      Matrix(IY,IY,IDIM) = UAVG(IX)
      Matrix(IY,IZ,IDIM) = ZERO
C
      DO I = 1 , NSP
        Matrix(IZ,I,IDIM) = -UAVG(IX)*UAVG(IZ)
      ENDDO
      Matrix(IZ,IE,IDIM) = ZERO
      Matrix(IZ,IX,IDIM) = UAVG(IZ)
      Matrix(IZ,IY,IDIM) = ZERO
      Matrix(IZ,IZ,IDIM) = UAVG(IX)
C
      IDIM = 2
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            Matrix(I,J,IDIM) = - ALPHA(I) * UAVG(IE+IDIM)
            IF (I.EQ.J)THEN
                Matrix(I,J,IDIM) = Matrix(I,J,IDIM) + UAVG(IE+IDIM)
            ENDIF
        ENDDO
      ENDDO
C
      DO I = 1 , NSP
        Matrix(I,IE,IDIM) = ZERO
        Matrix(I,IX,IDIM) = ZERO
        Matrix(I,IY,IDIM) = ALPHA(I)
        Matrix(I,IZ,IDIM) = ZERO
      ENDDO
C
      DO I = 1 , NSP
        Matrix(IE,I,IDIM) = UAVG(IE+IDIM) * TMPA(I)
      enddo
      Matrix(IE,IE,IDIM) = UAVG(IY) * (ONE + DE)
      Matrix(IE,IX,IDIM) = UAVG(IY) * DM(1)
      Matrix(IE,IY,IDIM) = UAVG(IE) + UAVG(IY) * DM(2)
      Matrix(IE,IZ,IDIM) = UAVG(IY) * DM(3)
C
      DO I = 1 , NSP
        Matrix(IX,I,IDIM) = -UAVG(IX)*UAVG(IY)
      ENDDO
      Matrix(IX,IE,IDIM) = ZERO
      Matrix(IX,IX,IDIM) = UAVG(IY)
      Matrix(IX,IY,IDIM) = UAVG(IX)
      Matrix(IX,IZ,IDIM) = ZERO
C
      DO I = 1 , NSP
        Matrix(IY,I,IDIM) = -UAVG(IY)*UAVG(IY) + DR(I)
      ENDDO
      Matrix(IY,IE,IDIM) = DE
      Matrix(IY,IX,IDIM) = DM(1)
      Matrix(IY,IY,IDIM) = DM(2) + TWO*UAVG(IY)
      Matrix(IY,IZ,IDIM) = DM(3)
C
      DO I = 1 , NSP
        Matrix(IZ,I,IDIM) = -UAVG(IY)*UAVG(IZ)
      ENDDO
      Matrix(IZ,IE,IDIM) = ZERO
      Matrix(IZ,IX,IDIM) = ZERO
      Matrix(IZ,IY,IDIM) = UAVG(IZ)
      Matrix(IZ,IZ,IDIM) = UAVG(IY)
C
C
      IDIM = 3
C
C
      DO I = 1 , NSP
        DO J = 1 , NSP
            Matrix(I,J,IDIM) = - ALPHA(I) * UAVG(IE+IDIM)
            IF (I.EQ.J)THEN
                Matrix(I,J,IDIM) = Matrix(I,J,IDIM) + UAVG(IE+IDIM)
            ENDIF
        ENDDO
      ENDDO
C      
      DO I = 1 , NSP
        Matrix(I,IE,IDIM) = ZERO
        Matrix(I,IX,IDIM) = ZERO
        Matrix(I,IY,IDIM) = ZERO
        Matrix(I,IZ,IDIM) = ALPHA(I)
      ENDDO
C
      DO I = 1 , NSP
        Matrix(IE,I,IDIM) = UAVG(IZ) * TMPA(I) 
      ENDDO
      Matrix(IE,IE,IDIM) = UAVG(IZ) * (ONE+DE)
      Matrix(IE,IX,IDIM) = DM(1) * UAVG(IZ)
      Matrix(IE,IY,IDIM) = DM(2) * UAVG(IZ)
      Matrix(IE,IZ,IDIM) = UAVG(IE) + DM(3) * UAVG(IZ)
C
      DO I = 1 , NSP
        Matrix(IX,I,IDIM) = -UAVG(IX)*UAVG(IZ)
      ENDDO
      Matrix(IX,IE,IDIM) = ZERO
      Matrix(IX,IX,IDIM) = UAVG(IZ)
      Matrix(IX,IY,IDIM) = ZERO
      Matrix(IX,IZ,IDIM) = UAVG(IX)
C
      DO I = 1 , NSP
        Matrix(IY,I,IDIM) = -UAVG(IY)*UAVG(IZ)
      ENDDO
      Matrix(IY,IE,IDIM) = ZERO
      Matrix(IY,IX,IDIM) = ZERO
      Matrix(IY,IY,IDIM) = UAVG(IZ)
      Matrix(IY,IZ,IDIM) = UAVG(IY)
C
      DO I = 1 , NSP
        Matrix(IZ,I,IDIM) = -UAVG(IZ)*UAVG(IZ)+DR(I)
      ENDDO
      Matrix(IZ,IE,IDIM) = DE
      Matrix(IZ,IX,IDIM) = DM(1)
      Matrix(IZ,IY,IDIM) = DM(2)
      Matrix(IZ,IZ,IDIM) = DM(3) + TWO*UAVG(IZ)
C
C
C
      ENDIF
C
C
      RETURN
      END
