      SUBROUTINE WUNV(XY,NDIM,ICELNOD,IBNDPTR,NPOIN,NELEM,NBFAC,
     &FNAME,LENFNAM)
C
C     Extrudes the 2D Triangle mesh into a single-cell-thick 3D
C     unstructured mesh (linear wedge/prism cells) and writes it as
C     an I-DEAS Universal file (<FNAME>.unv), the format
C     codes/code_saturne's preprocessor reads natively -- see
C     preprocessor/pre-post/ecs_pre_ideas.c: dataset 2411 = nodes,
C     2412 = elements (fe_des_id 112 = linear wedge, 94 = linear
C     shell quad, 91 = linear shell tri), 2430 = groups (entity
C     type 8 = element, tag = element label). Record layouts (field
C     counts, 4-per-line wrapping for 2430) were read directly out of
C     that parser, not from external documentation -- untested against
C     an actual code_saturne preprocessor run (not built in this
C     environment); verify before trusting it end to end.
C
C     code_saturne is a genuinely 3D FV code (no OpenFOAM-style
C     "empty" 2D direction), so the extruded front/back cap faces are
C     written as their own "front"/"back" boundary groups for the
C     caller to mark SYMMETRY in the case setup; the original 2D
C     boundary edges become quad side faces, one boundary group per
C     Triangle/.poly colour (bc<n>), matching triangle2su2's
C     convention.
C
      IMPLICIT NONE
      INTEGER NDIM,NPOIN,NELEM,NBFAC,LENFNAM
      CHARACTER*(*) FNAME
      DOUBLE PRECISION XY(NDIM,*)
      INTEGER ICELNOD(3,*),IBNDPTR(3,*)
C
      INTEGER MAXCLR
      PARAMETER(MAXCLR=20)
C
      INTEGER I,J,K,IPOIN,IELEM,IVERT,IN1,IN2,NGRP,LBL,ICOLR(MAXCLR)
      DOUBLE PRECISION XMIN,XMAX,YMIN,YMAX,DIAG,DZ
      CHARACTER FWORK*255
C
      INTEGER JCYCL
      EXTERNAL JCYCL
C
      FWORK(1:LENFNAM+4) = FNAME(1:LENFNAM)//'.unv'
      OPEN(UNIT=21,FILE=FWORK(1:LENFNAM+4))
C
C     ---- extrusion thickness: 1% of the mesh's bounding-box diagonal ----
C
      XMIN = XY(1,1)
      XMAX = XY(1,1)
      YMIN = XY(2,1)
      YMAX = XY(2,1)
      DO 5 IPOIN = 2,NPOIN
         XMIN = MIN(XMIN,XY(1,IPOIN))
         XMAX = MAX(XMAX,XY(1,IPOIN))
         YMIN = MIN(YMIN,XY(2,IPOIN))
         YMAX = MAX(YMAX,XY(2,IPOIN))
    5 CONTINUE
      DIAG = SQRT((XMAX-XMIN)**2 + (YMAX-YMIN)**2)
      DZ = 0.01D0*DIAG
C
C     ---- dataset 2411: nodes (front layer 1..NPOIN @ z=0, back layer
C          NPOIN+1..2*NPOIN @ z=DZ) ----
C
      WRITE(21,FMT='(A)')'    -1'
      WRITE(21,FMT='(A)')'  2411'
      DO 10 IPOIN = 1,NPOIN
         WRITE(21,FMT='(4I10)')IPOIN,1,1,1
         WRITE(21,FMT='(3(1X,ES24.16E3))')XY(1,IPOIN),XY(2,IPOIN),0.D0
   10 CONTINUE
      DO 20 IPOIN = 1,NPOIN
         WRITE(21,FMT='(4I10)')NPOIN+IPOIN,1,1,1
         WRITE(21,FMT='(3(1X,ES24.16E3))')XY(1,IPOIN),XY(2,IPOIN),DZ
   20 CONTINUE
      WRITE(21,FMT='(A)')'    -1'
C
C     ---- dataset 2412: elements ----
C     labels 1..NELEM            = wedge cells (extruded triangles)
C           NELEM+1..NELEM+NBFAC = side quads (extruded bndry edges)
C           next NELEM           = front cap triangles (z=0)
C           next NELEM           = back  cap triangles (z=DZ)
C
      WRITE(21,FMT='(A)')'    -1'
      WRITE(21,FMT='(A)')'  2412'
C
      DO 30 IELEM = 1,NELEM
         WRITE(21,FMT='(6I10)')IELEM,112,1,1,1,6
         WRITE(21,FMT='(6I10)')ICELNOD(1,IELEM),ICELNOD(2,IELEM),
     &   ICELNOD(3,IELEM),NPOIN+ICELNOD(1,IELEM),
     &   NPOIN+ICELNOD(2,IELEM),NPOIN+ICELNOD(3,IELEM)
   30 CONTINUE
C
      DO 40 I = 1,NBFAC
         LBL = NELEM + I
         J = IBNDPTR(3,I)
         IELEM = IBNDPTR(1,I)
         IVERT = IBNDPTR(2,I)
         IN1 = ICELNOD(JCYCL(IVERT+1),IELEM)
         IN2 = ICELNOD(JCYCL(IVERT+2),IELEM)
         WRITE(21,FMT='(6I10)')LBL,94,1,1,J,4
         WRITE(21,FMT='(4I10)')IN1,IN2,NPOIN+IN2,NPOIN+IN1
   40 CONTINUE
C
      DO 50 IELEM = 1,NELEM
         LBL = NELEM + NBFAC + IELEM
         WRITE(21,FMT='(6I10)')LBL,91,1,1,MAXCLR+1,3
         WRITE(21,FMT='(3I10)')ICELNOD(1,IELEM),ICELNOD(2,IELEM),
     &   ICELNOD(3,IELEM)
   50 CONTINUE
      DO 60 IELEM = 1,NELEM
         LBL = 2*NELEM + NBFAC + IELEM
         WRITE(21,FMT='(6I10)')LBL,91,1,1,MAXCLR+2,3
         WRITE(21,FMT='(3I10)')NPOIN+ICELNOD(1,IELEM),
     &   NPOIN+ICELNOD(2,IELEM),NPOIN+ICELNOD(3,IELEM)
   60 CONTINUE
C
      WRITE(21,FMT='(A)')'    -1'
C
C     ---- dataset 2430: groups (one per boundary colour + front/back) ----
C     obsolete-but-simple layout: each group entity is 2 ints (type,
C     tag), 4 entities per 8I10 line; type=8 marks an element label
C
      DO 70 I = 1,MAXCLR
         ICOLR(I) = 0
   70 CONTINUE
      DO 80 I = 1,NBFAC
         J = IBNDPTR(3,I)
         IF(J.GE.1.AND.J.LE.MAXCLR)ICOLR(J) = ICOLR(J) + 1
   80 CONTINUE
C
      WRITE(21,FMT='(A)')'    -1'
      WRITE(21,FMT='(A)')'  2430'
C
      NGRP = 0
      DO 100 I = 1,MAXCLR
         IF(ICOLR(I).LE.0)GOTO 100
         NGRP = NGRP + 1
         WRITE(21,FMT='(8I10)')NGRP,0,0,0,0,0,0,ICOLR(I)
         WRITE(21,FMT='(A,I0)')'bc',I
         K = 0
         DO 90 J = 1,NBFAC
            IF(IBNDPTR(3,J).NE.I)GOTO 90
            LBL = NELEM + J
            WRITE(21,FMT='(2I10)',ADVANCE='NO')8,LBL
            K = K + 1
            IF(MOD(K,4).EQ.0)WRITE(21,FMT='()')
   90    CONTINUE
         IF(MOD(K,4).NE.0)WRITE(21,FMT='()')
  100 CONTINUE
C
      NGRP = NGRP + 1
      WRITE(21,FMT='(8I10)')NGRP,0,0,0,0,0,0,NELEM
      WRITE(21,FMT='(A)')'front'
      K = 0
      DO 110 IELEM = 1,NELEM
         LBL = NELEM + NBFAC + IELEM
         WRITE(21,FMT='(2I10)',ADVANCE='NO')8,LBL
         K = K + 1
         IF(MOD(K,4).EQ.0)WRITE(21,FMT='()')
  110 CONTINUE
      IF(MOD(K,4).NE.0)WRITE(21,FMT='()')
C
      NGRP = NGRP + 1
      WRITE(21,FMT='(8I10)')NGRP,0,0,0,0,0,0,NELEM
      WRITE(21,FMT='(A)')'back'
      K = 0
      DO 120 IELEM = 1,NELEM
         LBL = 2*NELEM + NBFAC + IELEM
         WRITE(21,FMT='(2I10)',ADVANCE='NO')8,LBL
         K = K + 1
         IF(MOD(K,4).EQ.0)WRITE(21,FMT='()')
  120 CONTINUE
      IF(MOD(K,4).NE.0)WRITE(21,FMT='()')
C
      WRITE(21,FMT='(A)')'    -1'
      CLOSE(21)
C
      WRITE(6,*)'WUNV: wrote ',FWORK(1:LENFNAM+4),' -- ',NELEM,
     &' wedge cells, ',NGRP,' boundary groups, dz = ',DZ
C
      RETURN
      END
