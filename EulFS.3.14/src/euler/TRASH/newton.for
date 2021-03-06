      SUBROUTINE NEW(NDIM,NOFVAR,NOFVERT,JAC,LDJ,VCN,VCZ,ZEPS,VOLUME)

      IMPLICIT NONE

C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=5)
      DOUBLE PRECISION ROOT_MACHINE_EPS
      PARAMETER (ROOT_MACHINE_EPS=1.d-07)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION VOLUME
      INTEGER LDJ,NDIM,NOFVAR,NOFVERT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION JAC(LDJ,LDJ),VCN(NDIM,NOFVERT),
     +                 VCZ(NOFVAR,NOFVERT),ZEPS(NOFVAR,NOFVERT)
      double precision dndq,area 
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPS
      INTEGER I,IELEM,IFAIL,IVERT,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DKNEG(NMAX*NMAX),DKPOS(NMAX*NMAX),
     +                 KMAT(NMAX*NMAX),KMAT2(NMAX*NMAX),KNEG(NMAX*NMAX),
     +                 KNEG2(NMAX*NMAX),KPOS(NMAX*NMAX),
     +                 KPOS2(NMAX*NMAX),LNEG(NMAX),LPOS(NMAX),
     +                 VLEFT(NMAX*NMAX),VRIGHT(NMAX*NMAX),WR(NMAX),
     +                 APOS(NMAX*NMAX),ANEG(NMAX*NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DSCAL,LINEARIZE,MATSPLITVIII,X04CAF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,SIGN
C     ..
      CALL DINIT(NMAX*NMAX,0.d0,APOS,1) 
      CALL DINIT(NMAX*NMAX,0.d0,ANEG,1) 
      IELEM = -1
      DO 1 IVERT = 1,NOFVERT
          CALL MATSPLITVIII(IELEM,NDIM,NOFVAR,VCN(1,IVERT),JAC,LDJ,KMAT,
     +                      KPOS,KNEG,VLEFT,VRIGHT,NMAX,WR,LPOS,LNEG,
     +                      .TRUE.)
          DO 4 J = 1,NOFVERT
              area = vcn(1,j)**2 + vcn(2,j)**2
              if( ndim .eq. 3 )area = area + vcn(3,j)**2
              area = sqrt(area) 
              DO 3 I = 1,NOFVAR
                  CALL DCOPY(NOFVAR*NOFVERT,VCZ,1,ZEPS,1)
                  EPS = ROOT_MACHINE_EPS*MAX(ABS(ZEPS(I,J)),1.d0)*
     +                  SIGN(1.d0,ZEPS(I,J))
                  ZEPS(I,J) = ZEPS(I,J)+EPS
                  CALL LINEARIZE(IELEM,.FALSE.,VCN,NDIM,NOFVERT,ZEPS,
     +                           NOFVAR,VOLUME)
                  CALL MATSPLITVIII(IELEM,NDIM,NOFVAR,VCN(1,IVERT),JAC,
     +                              LDJ,KMAT2,KPOS2,KNEG2,VLEFT,VRIGHT,
     +                              NMAX,WR,LPOS,LNEG,.TRUE.)
                  CALL DCOPY(NOFVAR*NOFVAR,KPOS,1,DKPOS,1)
                  CALL DCOPY(NOFVAR*NOFVAR,KNEG,1,DKNEG,1)
                  CALL DAXPY(NOFVAR*NOFVAR,-1.d0,KPOS2,1,DKPOS,1)
                  CALL DAXPY(NOFVAR*NOFVAR,-1.d0,KNEG2,1,DKNEG,1)
                  CALL DSCAL(NOFVAR*NOFVAR,1.d0/EPS,DKPOS,1)
                  CALL DSCAL(NOFVAR*NOFVAR,1.d0/EPS,DKNEG,1)
                  write(6,*)i,j
C                 CALL X04CAF('General',' ',NOFVAR,NOFVERT,VCZ,NOFVAR,
C    +                        'U',IFAIL)
C                 CALL X04CAF('General',' ',NOFVAR,NOFVERT,ZEPS,NOFVAR,
C    +                        'dU',IFAIL)
                  CALL X04CAF('General',' ',NOFVAR,NOFVAR,KPOS,NOFVAR,
     +                        'K+',IFAIL)
                  CALL X04CAF('General',' ',NOFVAR,NOFVAR,KPOS2,NOFVAR,
     +                        'K+(U+eps)',IFAIL)
                  CALL X04CAF('General',' ',NOFVAR,NOFVAR,DKPOS,NOFVAR,
     +                        'd(K+)',IFAIL)
                  if( i .eq. 1 )then
                      dndq = 0.d0
                  else
                      dndq = vcn(i-1,j)/area
                  endif
                  CALL JACVIII(dndq,VCN(1,J),APOS,ANEG,NOFVAR,NDIM)
                  CALL X04CAF('General',' ',NOFVAR,NOFVAR,APOS,NOFVAR,
     +                        'd(K+)/dQ',IFAIL)
                  pause

    3         CONTINUE
    4     CONTINUE
    1 CONTINUE
      RETURN

      END
