      LENW    = NBFAC*(NOFVERT-1)
c     max facce di un dato colore X (nofvert-1)
      LNLIST  = ISTKGT((NOFVERT-1)*

      SUBROUTINE SUBXAA(ICELNOD,ICELCEL,ICELFAC,NOFVERT,NELEM,CORG,NDIM,
     +NPOIN,LENBC,NFACE,NBFAC,NBINT)
      CHKBND(IBNDFAC,NBFAC,ICELFAC,ICELNOD,NOFVERT,NELEM,
     +FACENORM,COORD,NDIM,NFACE,NWFAC,NBODY4,NBODY6)
C
      IMPLICIT NONE
      INCLUDE 'bnd.h'
C
      INCLUDE 'constants'
      INCLUDE 'bnd'
      INCLUDE 'io.com'
C
C     This routines 
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NFACE,NOFVERT,NPOIN
C     ..
C     .. Array Arguments ..
      INTEGER IBNDFAC(3,*),ICELNOD(NOFVERT,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUMMY
      INTEGER IELEM,IFACE,IFAIL,IFREQ,I,J,K,L,N,NERR,IOPT
      CHARACTER*72 ERRMSG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL,I1MACH
      EXTERNAL DDOT,DNRM2,ICYCL,I1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL CROSS_PROD,DAXPY,DSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,INT,ISIGN,MAX0,SIGN
C     ..
C     .. Data statements ..
C     ..
C
      INTEGER NBFAC,NOFVERT,NELEM,NDIM,NFACE
      INTEGER NN,NERR,IOPT
      PARAMETER(NN=NBTYPE+1)
      CHARACTER ERRMSG*72
C
      INTEGER IBC,J,JCOLOR,IFAIL,i,IFRST
      INTEGER IELEM,IBFAC,IFREQ,IVERT,IFACE
      DOUBLE PRECISION TEMP,TEMP1
C
      INTEGER LENBC(0:NBTYPE)
C
      EXTERNAL DAXPY
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
C
      INTRINSIC MAX0
C
C
      DATA (LENBC(J),J=0,NBTYPE),IFAIL/ NN*0,0 /
      DATA ERRMSG(1:7)/'BNDCHK '/
C
C
C
      DO 14 IBC = 0 , NCOLOR
         LENBC(IBC) = 0
   14 CONTINUE
C
C
      DO 16 IBFAC = 1 , NBFAC
C
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
C
C     NON il colore !!!
C
         IBC   = IBNDFAC(3,IBFAC)
C
C     .. pick up the nodes on the boundary face
C
         DO 14 J = 1,(NOFVERT-1)
             JVERT = ICYCL(IVERT+J,NOFVERT) 
             IPOIN = ICELNOD(JVERT,IELEM) 
             LENBC(IBC) = LENBC(IBC) + 1
             NLIST(LENBC(IBC),IBC) = IPOIN
   14    CONTINUE
   16 CONTINUE
C
C  .. remove duplicated entries from the list of bnd nodes
C
      LENW = 0
      DO 12 IBC = 0 , NCOLOR
          IF(LENBC(IBC).EQ.0)GOTO 12
          N = LENBC(IBC)
          CALL SORTSP(N,NLIST(1,IBC),LENBC(IBC))
          LENW = LENW + LENBC(IBC)
          WRITE(NOUT,355)LENBC(IBC),IBC
C
C  qui dovrei spostare le liste
C
   12 CONTINUE
C
C
      RETURN
C
C     I/O FORMATS
C
  100 FORMAT(I3.3)
  110 FORMAT(I2.2)
  111 FORMAT('.',$)
  114 FORMAT(10X,'CHECKING BOUNDARY FACES ',$)
  115 FORMAT(I6,/)
  350 FORMAT(10X,A20,' BOUNDARY FACES : ',I6,2X,'(',D9.4,')')
  355 FORMAT(10X,'THERE ARE ',I5,' BOUNDARY POINTS COLOURED ',I2)
  454 FORMAT(57X,9("-"),/57X,D9.4)
  500 FORMAT('Inconsistent data : IBNDFAC(',I1,',',I6,') = ',I6)
C
      END

      LVTXNOR = ISTKGT(LENW*NDIM,4) 
      LVTXPTR = ISTKGT(LENW,2) 
      LDEGREE = ISTKGT(LENW,2) 

      SUBROUTINE GETNOR(ICELNOD,ICELCEL,ICELFAC,NOFVERT,NELEM,
     +              NDIM,NPOIN,FACNOR,NFACE,VTXNOR,NNORM,NBFAC,NBINT)
C
      IMPLICIT NONE
C
      INCLUDE 'constants'
      INCLUDE 'io.com'
C
C     This routines finds all faces of a 3D mesh
C
C     .. Scalar Arguments ..
      INTEGER NBFAC,NDIM,NELEM,NFACE,NOFVERT,NPOIN,NBINT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION VTXNOR(NDIM,NNORM),FACNOR(NDIM,NFACE)
      INTEGER VTXLST(NNORM),ICELFAC(NOFVERT,NELEM),
     +        ICELNOD(NOFVERT,NELEM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DUMMY,S,W
      INTEGER IELEM,IFACE,IFAIL,IFREQ,II,IV,I,J,K,L,LL,N,JV,JVERT
     +        NEIGHB,NERR,IOPT
      CHARACTER*72 ERRMSG
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WKSP(3)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      INTEGER ICYCL,I1MACH
      EXTERNAL DDOT,DNRM2,ICYCL,I1MACH
C     ..
C     .. External Subroutines ..
      EXTERNAL CROSS_PROD,DAXPY,DSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,IABS,INT,ISIGN,MAX0,SIGN
C     ..
C     .. Data statements ..
C     ..
C
      IFACE = 0
      IFREQ = MAX0(1,INT(NFACE)/20)
      WRITE (NOUT,FMT=113)
C
C Loop over the cells to compute node based normals
C
      DO 16 IBFAC = 1 , NBFAC
C
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
         IBC   = IBNDFAC(3,IBFAC)
C
C     .. pick up nodes on the boundary face
C
         DO 14 J = 1,(NOFVERT-1)
             JVERT = ICYCL(IVERT+J,NOFVERT) 
             IPOIN = ICELNOD(JVERT,IELEM) 
             IFACE = ICELCEL(JVERT,IELEM)
             ALPHA = DSIGN(ONE,IFACE)
             IFACE = IABS(IFACE)
             CALL BINSRC(IPOIN,NLIST(1,IBC),NLEN(IBC),IPOS,LAST)
             IF( IPOS .LT. 0 )THEN 
             ENDIF
             DEGREE(IPOS) = DEGREE(IPOS) + 1
             CALL DAXPY(NDIM,ALPHA,FACNOR(1,IFACE),1,VTXNOR(1,IPOS),1)
   14    CONTINUE
   16 CONTINUE
      DO 26 I = 1,NNORM
         ALPHA = ONE/REAL(DEGREE(I))
         CALL DSCAL(NDIM,ALPHA,VTXNOR(1,I),1)
   26 CONTINUE
C
      RETURN
C
C     I/O FORMATS
C
  111 FORMAT ('.',$)
  113 FORMAT (10X,'COMPUTING NORMALS ',$)
  115 FORMAT (I6,/)
  207 FORMAT (5X,10 ('*'),' WARNING ',10 ('*'),/,5X,I7,
     +       ' BOUNDARY FACES(EDGES) WERE EXPECTED',/,5X,I7,
     +       ' INTER-PROC FACES(EDGES) HAVE BEEN FOUND',/,5X,I7,
     +       ' BOUNDARY FACES(EDGES) HAVE BEEN FOUND',/)
  208 FORMAT ('FF -- CHECK ON BOUNDARY FACES(EDGES) FAILED')
  205 FORMAT (5X,10 ('*'),' ERROR ',10 ('*'),/,5X,I7,
     +       ' FACES(EDGES) WERE EXPECTED',/,5X,I7,
     +       ' FACES(EDGES) HAVE BEEN FOUND',/)
  206 FORMAT ('FF -- TOO MUCH ROOM TO STORE THE NORMALS') 
  210 FORMAT (3X,'WARNING !! ZERO NORMAL FOR FACE ',I5,' NODES ',
     +       3 (3X,I6))
  405 FORMAT (/,25X,10 ('*'),' ERROR ',10 ('*'),/,5X,
     +       'THE ESTIMATED NUMBER OF FACES IS INSUFFICIENT',/,5X,
     +       'CHECK THE CONNECTIVITY OF THE MESH')
  406 FORMAT (' FF -- NOT ENOUGH ROOM TO STORE THE NORMALS')
  999 FORMAT (5X,'ERROR : Scaled normals in cell ',I6,' sum up to ',
     +       E9.4)
      END
