      PROGRAM SU2TRIANGLE
C
C     converts SU2's restart-state output back into a triangle .node
C     file, the SU2 analogue of dat2triangle. Re-reads the same
C     Triangle .node/.ele/.neigh/.edge/.poly files triangle2su2 read
C     (to recover connectivity/boundary structure exactly, the same
C     way dat2triangle relies on EulFS's own file001/file002.dat
C     rather than re-deriving it), then overlays the state SU2
C     computed by reading "<fname>_restart.csv" (see rsu2.f).
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NELEM,NPOIN,ns,nfr,nt,K,NHOLE,mdim,nbc,i,NBPOIN
C
      INTEGER NVA
      PARAMETER (NVA=39800000)
C
      DOUBLE PRECISION DSTAK(NVA)
      INTEGER ISTAK(1)
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      COMMON /CSTAK/DSTAK
C
      INTEGER LCORG,LCELCEL,LCELNOD,LBNDPTR,LOCZ,LNODCOD,NITEMS
C
      INTEGER   NOFVERT,NDIM,NDOF
      PARAMETER(NDIM=2,NOFVERT=(NDIM+1),NDOF=4)
C
      INTEGER ILEN(10),ICHOLE(1),KSU2
C
      INTEGER  ISTKGT,LENSTR
      EXTERNAL ISTKGT,LENSTR
C
      DOUBLE PRECISION EPS,GAM
      COMMON /EPSCOM/EPS
      CHARACTER FNAME(10)*255,SU2NAME*255
C
      equivalence(nt,nelem)
      equivalence(ns,npoin)
      equivalence(nfr,nbfac)
C
C     .. Executable Statements ..
C
      GAM = 1.4D0
C
      WRITE(6,*)
      WRITE(6,*)
     &'This code converts SU2 restart files into triangle files '
      WRITE(6,*)
      FNAME(10) = "file.1"
      WRITE(6,*)'Enter fname'
      READ(5,*)FNAME(10)
      SU2NAME = "file"
      WRITE(6,*)'Enter SU2 input basename'
      READ(5,*)SU2NAME
      KSU2 = LENSTR(SU2NAME)
      K = LENSTR(FNAME(10))
      FNAME(1)(1:K+5) = FNAME(10)(1:K)//".node"
      FNAME(2)(1:K+4) = FNAME(10)(1:K)//".ele"
      FNAME(3)(1:K+6) = FNAME(10)(1:K)//".neigh"
      FNAME(4)(1:K+5) = FNAME(10)(1:K)//".edge"
      FNAME(5)(1:K+5) = FNAME(10)(1:K)//".poly"
      ILEN(1) = K+5
      ILEN(2) = K+4
      ILEN(3) = K+6
      ILEN(4) = K+5
      ILEN(5) = K+5
C
      open(1,file=FNAME(1)(1:K+5),status="old")
      rewind(1)
      read(1,*) ns,i,i,i
      close(1)
C
      open(1,file=FNAME(2)(1:K+4),status="old")
      rewind(1)
      read(1,*) nt
      close(1)
C
      open(1,file=FNAME(4)(1:K+5),status="old")
      rewind(1)
      read(1,*) nfr
      close(1)
C
      open(1,file=FNAME(5)(1:K+5),status="old")
      rewind(1)
      read(1,*)nitems,mdim,i,nbc
      do i = 1,nitems
         read(1,*)
      enddo
      read(1,*)nitems,nbc
      do i = 1,nitems
         read(1,*)
      enddo
      read(1,*)nhole
      close(1)
C
      write(6,50)ns,nt,nfr,nhole
   50 format(5X,'NPOIN = ',I8,5X,'NELEM ',I8,5X,'NFACE = ',I8,' NHOLE
     &= ',I2)
C
      CALL ISTKIN(NVA,4)
C
      LCORG = ISTKGT(NPOIN*NDIM,4)
      CALL DINIT(NPOIN*NDIM,0.d0,DSTAK(LCORG),1)
C
      LOCZ = ISTKGT(NDOF*NPOIN,4)
      CALL DINIT(NPOIN*NDOF,0.d0,DSTAK(LOCZ),1)
C
      LCELNOD = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELNOD),1)
C
      LCELCEL = ISTKGT(NELEM*NOFVERT,2)
      CALL IINIT(NELEM*NOFVERT,0,ISTAK(LCELCEL),1)
C
      LBNDPTR = ISTKGT(3*NBFAC,2)
      CALL IINIT(3*NBFAC,0,ISTAK(LBNDPTR),1)
C
      LNODCOD = ISTKGT(NPOIN,2)
      CALL IINIT(NPOIN,0,ISTAK(LNODCOD),1)
C
      WRITE(6,*)'Reading mesh..........'
      CALL RTRI(DSTAK(LCORG),NPOIN,DSTAK(LOCZ),NDOF,ISTAK(LCELNOD),
     &ISTAK(LCELCEL),ISTAK(LBNDPTR),NELEM,NBFAC,FNAME,ILEN)
      WRITE(6,*)'Done'
C
      WRITE(6,*)'Reading SU2 restart file..........'
      CALL RSU2(SU2NAME,KSU2,DSTAK(LOCZ),NDOF,NPOIN,ISTAK(LCELNOD),
     &NELEM,GAM)
      WRITE(6,*)'Done'
C
      NBPOIN = 0
      CALL NODCOD(ISTAK(LNODCOD),NPOIN,NBPOIN,ISTAK(LCELNOD),NOFVERT,
     &NELEM,ISTAK(LBNDPTR),NBFAC)
C
C     back up the pre-solve .node file before WTRI overwrites it in
C     place, same convention as dat2triangle/main.f
C
      CALL SYSTEM('mv '//FNAME(1)(1:K+5)//' '//FNAME(1)(1:K+5)//
     &'.BAK')
C
      CALL WTRI(ISTAK(LBNDPTR),NBFAC,ISTAK(LCELNOD),NOFVERT,
     &DSTAK(LCORG),NDIM,DSTAK(LOCZ),NDOF,ISTAK(LNODCOD),NPOIN,
     &FNAME(10),ICHOLE,0)
C
      CALL EXIT(0)
      END
