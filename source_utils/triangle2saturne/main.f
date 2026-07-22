      PROGRAM TRIANGLESATURNE
C
C     converts a Triangle mesh into an I-DEAS Universal file (.unv)
C     that codes/code_saturne's preprocessor can read, extruding the
C     2D triangulation into a single-cell-thick 3D wedge mesh (see
C     wunv.f). The code_saturne analogue of triangle2dat/triangle2su2;
C     see that program's header for the required
C     "triangle -n -e -p filename.poly" invocation.
C
C     No periodic-surface support (same scope note as triangle2su2).
C     No solution-state conversion here either, unlike triangle2su2:
C     code_saturne's restart/checkpoint format isn't a plain text
C     format like SU2's CSV, so seeding/reading back Roe-vector state
C     is left as a follow-up (see issue #4) -- this program is
C     mesh-only for now.
C
      IMPLICIT NONE
C
      INTEGER NBFAC,NELEM,NPOIN,ns,nfr,nt,K,NHOLE,mdim,nbc,i
C
      INTEGER NVA
      PARAMETER (NVA=39800000)
C
      DOUBLE PRECISION DSTAK(NVA)
      INTEGER ISTAK(1)
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      COMMON /CSTAK/DSTAK
C
      INTEGER LCORG,LCELCEL,LCELNOD,LBNDPTR,LOCZ,NITEMS
C
      INTEGER   NOFVERT,NDIM,NDOF
      PARAMETER(NDIM=2,NOFVERT=(NDIM+1),NDOF=4)
C
      INTEGER ILEN(10)
C
      INTEGER  ISTKGT,LENSTR
      EXTERNAL ISTKGT,LENSTR
C
      DOUBLE PRECISION EPS
      COMMON /EPSCOM/EPS
      CHARACTER FNAME(10)*255
C
      equivalence(nt,nelem)
      equivalence(ns,npoin)
      equivalence(nfr,nbfac)
C
C     .. Executable Statements ..
C
      WRITE(6,*)
      WRITE(6,*)'This code converts triangle files into a '//
     &'code_saturne (.unv) mesh '
      WRITE(6,*)
      FNAME(10) = "file.1"
      WRITE(6,*)'Enter fname'
      READ(5,*)FNAME(10)
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
      WRITE(6,*)'Trying to open ',FNAME(1)(1:K+5)
      open(1,file=FNAME(1)(1:K+5),status="old")
      rewind(1)
      read(1,*) ns,i,i,i
      close(1)
C
      WRITE(6,*)'Trying to open ',FNAME(2)(1:K+4)
      open(1,file=FNAME(2)(1:K+4),status="old")
      rewind(1)
      read(1,*) nt
      close(1)
C
      WRITE(6,*)'Trying to open ',FNAME(4)(1:K+5)
      open(1,file=FNAME(4)(1:K+5),status="old")
      rewind(1)
      read(1,*) nfr
      close(1)
C
      WRITE(6,*)'Trying to open ',FNAME(5)(1:K+5)
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
      WRITE(6,*)'Reading mesh..........'
      CALL RTRI(DSTAK(LCORG),NPOIN,DSTAK(LOCZ),NDOF,ISTAK(LCELNOD),
     &ISTAK(LCELCEL),ISTAK(LBNDPTR),NELEM,NBFAC,FNAME,ILEN)
      WRITE(6,*)'Done'
C
      CALL WUNV(DSTAK(LCORG),NDIM,ISTAK(LCELNOD),ISTAK(LBNDPTR),
     &NPOIN,NELEM,NBFAC,FNAME(10),K)
C
      CALL EXIT(0)
      END
