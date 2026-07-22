      SUBROUTINE WSU2(XY,NDIM,ZROE,NDOF,ICELNOD,IBNDPTR,NPOIN,NELEM,
     &NBFAC,GAM,FNAME,LENFNAM)
C
C     Writes a SU2-native ASCII mesh (<FNAME>.su2) and a restart-state
C     file (<FNAME>_restart.csv) from Triangle connectivity plus the
C     Roe-parameter-vector state already carried on the mesh, mirroring
C     what triangle2dat's RTRI/main.f do for the EulFS .dat format.
C
C     .su2 mesh format (NDIME/NELEM/NPOIN/NMARK sections, VTK element
C     type codes: 5=triangle, 3=line, zero-based numbering) checked
C     against the vendored codes/SU2/QuickStart/mesh_NACA0012_inv.su2
C     example.
C
C     Roe parameter vector convention, per
C     codes/UnDiFi-2D/doc/userguide/040_solver.md:
C       Z = ( sqrt(rho), sqrt(rho)*H, sqrt(rho)*u, sqrt(rho)*v )
C     converted here to the conservative variables (rho,rho*u,rho*v,
C     rho*E) that SU2's restart/solution files expect.
C
      IMPLICIT NONE
C
      INTEGER NDIM,NDOF,NPOIN,NELEM,NBFAC,LENFNAM
      CHARACTER*(*) FNAME
      DOUBLE PRECISION XY(NDIM,*),ZROE(NDOF,*),GAM
      INTEGER ICELNOD(3,*),IBNDPTR(3,*)
C
      INTEGER MAXCLR
      PARAMETER(MAXCLR=20)
C
      INTEGER I,J,K,IPOIN,IELEM,IVERT,IN1,IN2,NCOLR,ICOLR(MAXCLR)
      INTEGER NPOUT
      INTEGER IMAP(NPOIN)
      DOUBLE PRECISION RHO,H,U,V,Q2,P,RHOE,GM1
      CHARACTER FWORK*255,MSHNAME*255
C
      INTEGER JCYCL
      EXTERNAL JCYCL
C
      GM1 = GAM - 1.D0
C
C     Triangle's phantom-node preallocation leaves large numbers of
C     unreferenced, coincident points (marker -99) in the .node file
C     that ICELNOD never touches; EulFS/NEO tolerate the padding, but
C     SU2 requires every declared NPOIN entry to be referenced by some
C     element (CPhysicalGeometry::DistributeColoring's "Mismatch
C     between NPOIN and number of points listed in mesh file" check).
C     Compact the numbering down to only the points ICELNOD actually
C     uses, and remap element/boundary connectivity through IMAP.
C
      DO 2 IPOIN = 1,NPOIN
         IMAP(IPOIN) = 0
    2 CONTINUE
      DO 4 IELEM = 1,NELEM
         IMAP(ICELNOD(1,IELEM)) = 1
         IMAP(ICELNOD(2,IELEM)) = 1
         IMAP(ICELNOD(3,IELEM)) = 1
    4 CONTINUE
      NPOUT = 0
      DO 6 IPOIN = 1,NPOIN
         IF(IMAP(IPOIN).NE.0)THEN
            NPOUT = NPOUT + 1
            IMAP(IPOIN) = NPOUT
         ENDIF
    6 CONTINUE
C
C     ---------------------------------------------------------------
C     mesh file
C     ---------------------------------------------------------------
C
      FWORK(1:LENFNAM+4) = FNAME(1:LENFNAM)//'.su2'
      MSHNAME(1:LENFNAM+4) = FWORK(1:LENFNAM+4)
      OPEN(UNIT=21,FILE=FWORK(1:LENFNAM+4))
C
      WRITE(21,FMT='(A,I0)')'NDIME= ',NDIM
C
      WRITE(21,FMT='(A,I0)')'NELEM= ',NELEM
      DO 10 IELEM = 1,NELEM
         WRITE(21,FMT='(I0,3(1X,I0),1X,I0)')5,
     &   IMAP(ICELNOD(1,IELEM))-1,IMAP(ICELNOD(2,IELEM))-1,
     &   IMAP(ICELNOD(3,IELEM))-1,IELEM-1
   10 CONTINUE
C
      WRITE(21,FMT='(A,I0)')'NPOIN= ',NPOUT
      DO 20 IPOIN = 1,NPOIN
         IF(IMAP(IPOIN).EQ.0)GOTO 20
         WRITE(21,FMT='(1X,ES24.16E3,1X,ES24.16E3,1X,I0)')
     &   XY(1,IPOIN),XY(2,IPOIN),IMAP(IPOIN)-1
   20 CONTINUE
C
C     group boundary faces by their Triangle/.poly colour (IBNDPTR(3,*))
C     into SU2 markers named "bc<colour>"; colour <= 0 is not expected
C     for a genuine boundary edge (see RTRI) and is skipped with a
C     warning rather than silently mis-tagged.
C
      DO 30 I = 1,MAXCLR
         ICOLR(I) = 0
   30 CONTINUE
      DO 40 I = 1,NBFAC
         J = IBNDPTR(3,I)
         IF(J.GE.1.AND.J.LE.MAXCLR)THEN
            ICOLR(J) = ICOLR(J) + 1
         ELSE
            WRITE(6,*)'WSU2: WARNING boundary face ',I,
     &      ' has an unexpected colour ',J,' -- skipped'
         ENDIF
   40 CONTINUE
C
      NCOLR = 0
      DO 50 I = 1,MAXCLR
         IF(ICOLR(I).GT.0)NCOLR = NCOLR + 1
   50 CONTINUE
C
      WRITE(21,FMT='(A,I0)')'NMARK= ',NCOLR
      DO 70 I = 1,MAXCLR
         IF(ICOLR(I).LE.0)GOTO 70
         WRITE(21,FMT='(A,I0)')'MARKER_TAG= bc',I
         WRITE(21,FMT='(A,I0)')'MARKER_ELEMS= ',ICOLR(I)
         K = 0
         DO 60 J = 1,NBFAC
            IF(IBNDPTR(3,J).NE.I)GOTO 60
            IELEM = IBNDPTR(1,J)
            IVERT = IBNDPTR(2,J)
            IN1 = IMAP(ICELNOD(JCYCL(IVERT+1),IELEM))
            IN2 = IMAP(ICELNOD(JCYCL(IVERT+2),IELEM))
            WRITE(21,FMT='(I0,2(1X,I0),1X,I0)')3,IN1-1,IN2-1,K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
      CLOSE(21)
C
C     ---------------------------------------------------------------
C     restart-state file (conservative variables, SU2 CSV convention)
C     ---------------------------------------------------------------
C
      FWORK(1:LENFNAM+12) = FNAME(1:LENFNAM)//'_restart.csv'
      OPEN(UNIT=22,FILE=FWORK(1:LENFNAM+12))
      WRITE(22,FMT='(A)')
     &'"PointID","x","y","Density","Momentum_x","Momentum_y",'//
     &'"Energy"'
      DO 80 IPOIN = 1,NPOIN
         IF(IMAP(IPOIN).EQ.0)GOTO 80
         RHO  = ZROE(1,IPOIN)**2
         H    = ZROE(2,IPOIN)/ZROE(1,IPOIN)
         U    = ZROE(3,IPOIN)/ZROE(1,IPOIN)
         V    = ZROE(4,IPOIN)/ZROE(1,IPOIN)
         Q2   = U*U + V*V
         P    = RHO*GM1/GAM*(H - 0.5D0*Q2)
         RHOE = P/GM1 + 0.5D0*RHO*Q2
         WRITE(22,FMT='(I0,6('','',ES24.16E3))')IMAP(IPOIN)-1,
     &   XY(1,IPOIN),XY(2,IPOIN),RHO,RHO*U,RHO*V,RHOE
   80 CONTINUE
      CLOSE(22)
C
      WRITE(6,*)'WSU2: wrote ',MSHNAME(1:LENFNAM+4),' and a matching ',
     &'restart file with ',NCOLR,' boundary marker(s)'
C
      RETURN
      END
