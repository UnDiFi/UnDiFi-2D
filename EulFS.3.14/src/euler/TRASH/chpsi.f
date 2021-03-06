      SUBROUTINE CHPSI(ICELNOD,ICELFAC,VOL,NELEM,ZROE,NPOIN,VCZ,VFACNOR,
     &                 NDIM,NOFVERT,NOFVAR,ICHPSI)
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'three'
      INCLUDE 'constants'
C
      INTEGER NDIM,NOFVERT,NOFVAR,NELEM,NPOIN
C
C     On entry:
C     --------
C
C     NDIM    dimension of the space (2 or 3)
C     NOFVERT number of vertices per element (=NDIM+1, since
C             only triangles or tetrahedra are allowed)
C     NOFVAR  number of variables (degrees of freedom)
C             in each meshpoint
C     NELEM   no. of processor owned elements (triangles/tetrahedra);
C             global number of elements in the uni-processor case
C
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      INTEGER ICHPSI(NPOIN)
C
C     ICELNOD(1:NOFVERT,1:NELEM)
C            Cell to Node pointer : ICELNOD(i,ielem) gives the
C            global node number of the i-th vertex of the ielem-th cell
C
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM), ZROE(NOFVAR,*),
     +VCZ(NOFVAR,NOFVERT)
C
C     FACNOR(1:NDIM,1:NFACE)  cartesian components
C                             of the NFACE faces/edges
      INTEGER IVAR,JVERT,IELEM
      INTEGER ICN(4)
      DOUBLE PRECISION VCN(12),VOLUME,SUM
C
C     Compute CHPSI
C
      CALL IINIT(NPOIN,0,ICHPSI,1)
C
      DO 1999 IELEM = 1,NELEM

          CALL CELPTR(IELEM, ICELNOD, ICELFAC, VOL, ZROE, VFACNOR, NDIM,
     +    NOFVERT, NOFVAR, ICN, VCZ, VCN, VOLUME)

c     Compute Mach number:

      DO 10 IVAR = 1 , NOFVAR
       SUM = ZERO
            DO 12 JVERT = 1 , NOFVERT
               SUM = SUM + VCZ( IVAR , JVERT )
   12       CONTINUE
      ZAVG(IVAR) = SUM / NOFVERT
   10 CONTINUE

      IF(NDIM.EQ.3)
     &UAVG(5) = ZAVG(5)/ZAVG(1) ! z componenet of the velocity vector
      UAVG(4) = ZAVG(4)/ZAVG(1) ! y componenet of the velocity vector
      UAVG(3) = ZAVG(3)/ZAVG(1) ! x componenet of the velocity vector
      UAVG(2) = ZAVG(2)/ZAVG(1) ! Total Enthalpy
      UAVG(1) = ZAVG(1)*ZAVG(1) ! Density
C
      KINETIC = UAVG(3)*UAVG(3) + UAVG(4)*UAVG(4)
      IF( NDIM .EQ. 3 )KINETIC = KINETIC + UAVG(5)*UAVG(5)
      KINETIC = HALF * KINETIC
      ASQR = GM1 * ( UAVG(2) - KINETIC )
      IF( ASQR .LT. ZERO )THEN
         WRITE(6,*)'Negative averaged pressure in element ',IELEM
         STOP
      ENDIF
      ABAR      = SQRT(ASQR)
c     MACHSQR   = DMAX1( TWO * KINETIC / ASQR, EPSQR_STAG) ! Barth & Deconinck
      MACHSQR   = TWO * KINETIC / ASQR
      MACH      = SQRT( MACHSQR )

      if (MACH.ge.ONE) then

            DO JVERT = 1 , NOFVERT
               ICHPSI(ICN(JVERT)+1) = 1
            END DO

      end if

 1999 CONTINUE
C     do jvert=1,npoin
C        write(16,*)ichpsi(jvert)
C     enddo 
C     stop
C
      RETURN
      END
