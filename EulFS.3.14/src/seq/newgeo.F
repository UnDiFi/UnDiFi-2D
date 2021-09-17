!> \par Purpose
!>
!> here in NEWGEO() we re-compute some grid metrics, more precisely:
!>
!> a) measures of the triangle/tetrahedra (LVOL) are stored in VOL(NELEM+1:2*NELEM)
!>
!> b) measures of the median dual control volumes (LMEDIAN) are stored in MEDIAN(NP+1:2*NP)
!>
!> c) face normals (LFACNOR) are overwritten on the current ones 
!>
!> d) cell to face pointers (LCELFAC) are overwritten on the current ones, but should not change
!>        anyway, since it is assumed that the connectivity is un-changed 
!>
!> e) grid velocities LXYZDOT
!>
!>   this has to be done BEFORE calculating the residual i.e. before calling RHSFunction
!>
!>   in NEWGEO all metrics are updated at some time which is intermediate btw n and n+1 
!>   the areas/volumes of the median dual CVs are however those at n+1 since these
!>   are only needed to update the solution, not to compute the residual
!>
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELCEL Cell to cell pointer: \c ICELCEL(i,j) gives the element number that shares the face opposite the i-th vertex of the j-th element if 0 or > \c NELEM, that face is a boundary face
!> @param[in] PMAP Map for periodic nodes
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in,out] FACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in,out] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in,out] VMEDIAN area/volume of the median dual cells
!> @param[in,out] CORG Cartesian coordinates of the meshpoints; CORG(1,2*NP+1) stores the grid at time level 0, i.e. the fixed grid; new grid coordinates are put into CORG(1,NP+1) and the grid velocity is computed accordingly
!> @param[out] XYZDOT Cartesian components of the nodal grid velocities
!> @param[in] IBNDPTR Boundary informations: \c IBNDFAC(1,j) gives the element the j-th face belongs to; \c IBNDFAC(2,j) gives the local vertex number of element \c IBNDFAC(1,j) opposite the boundary face; \c IBNDFAC(3,j) gives the color of the boundary face
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NELEM nof boundary faces
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in] NGHOST nof ghost nodes in the mesh
!> @param[in] NPNOD nof periodic nodes in the mesh
!> @param[in] NFACE nof faces in the mesh
!> @param[in] NBFAC nof boundary faces
!> @param[in] NBINT nof inter-processor faces
!> @param[in] NHOLE nof holes in the mesh (should be 0 in 3D)
!> @param[in] TIME physical time to be reached, i.e. time \f$(n+1)\Delta t\f$.
!>
!>    Upon entry:
!>
!> \verbatim
!>
!>    LCORG       |---- n   ----| |---- XXX ----| |---- 0 ----|
!>    LMEDIAN     |---- n   ----| |---- XXX ----| |---- n ----|
!>    LVOL        |---- n   ----| |---- XXX ----| |----n-1----|
!>    LFACNOR     |---- n   ----|
!>
!>    XXX means that the entry is irrelevant
!> \endverbatim
!>
!>    Upon return:
!>
!> \verbatim
!>
!>    LCORG       |---- n+a ----| |---- n+1 ----| |---- 0 ----|
!>    LMEDIAN     |---- n   ----| |---- n+1 ----| |---- n ----|
!>    LVOL        |---- n+a ----| |---- n+1 ----| |----n-1----|
!>    LFACNOR     |---- n+a ----|
!> \endverbatim
!>
!> where \c n stands for the time level and \c a is the \c ALFALE parameter which is set through \c -ale_grid_weight
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.18 $
!> \date $Date: 2020/04/25 10:18:01 $
!>
      SUBROUTINE NEWGEO(MotionSolver,PetscMatC,NodalBcs,ICELNOD,ICELCEL,
     2PMAP,ICELFAC,FACNOR,VOL,VMEDIAN,CORG,XYZDOT,IBNDPTR,NDIM,NOFVERT,
     3NELEM,NPOIN,NGHOST,NPNOD,NFACE,NBFAC,NBINT,NHOLE,TIME)
C
C     $Id: newgeo.F,v 1.18 2020/04/25 10:18:01 abonfi Exp $
C
C
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscpc.h"
      use petscpc
      use petscksp
C
      IMPLICIT NONE
C
      INCLUDE "paramt.h"
      INCLUDE "constants.h"
      INCLUDE "time.h"
      INCLUDE "time.com"
      INCLUDE "conv.com"
      INCLUDE "io.com"
C
      INTEGER NDIM,NOFVERT,NELEM,NPOIN,NGHOST,NPNOD,NFACE,NBFAC,NBINT,
     2NHOLE
      DOUBLE PRECISION TIME
      DOUBLE PRECISION CORG(NDIM,*)
      INTEGER ICELNOD(NOFVERT,NELEM),ICELCEL(NOFVERT,NELEM),PMAP(*)
C
      Mat              PetscMatC
      IS, dimension(0:*) :: NodalBcs
      KSP MotionSolver
C
      DOUBLE PRECISION FACNOR(NDIM,NFACE) ,VOL(NELEM), VMEDIAN(*), 
     &XYZDOT(*)
      INTEGER I,INDX,IFAIL
      INTEGER ICELFAC(NOFVERT,NELEM),IBNDPTR(3,NBFAC)
C
      INTEGER NP,NITEMS,iw
      DOUBLE PRECISION TNOW,maxdxyzdt(3),avgdxyzdt(3)
      LOGICAL DEBUG
      PARAMETER(DEBUG=.FALSE.)
      CHARACTER*80 FNAME,GNAME,BNAME
C
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      DOUBLE PRECISION DASUM
      EXTERNAL DASUM
C
      DATA FNAME(1:14)/"grid000000.dat"/
      DATA GNAME(1:14)/"step000000.dat"/
      DATA BNAME(1:18)/"gridvel_000000.dat"/
C
      IFAIL = IBAK
      NP = NPOIN+NGHOST+NPNOD
C
      IF(DEBUG)THEN
         iw = 97
         call print4check(corg,vol,vmedian,facnor,ndim,np,nelem,
     &time,'After entering NEWgeo',iw)
      ENDIF 
!
!     We read the grid velocity at all gridpoints whenever ALE_READ_GRIDVEL = .TRUE.
!
      IF(ALE_READ_GRIDVEL)THEN
         WRITE(BNAME(9:14),FMT="(I6.6)")ITIM
         CALL SOLZNE(BNAME(1:18),XYZDOT,NDIM,NP,"read")
         IALE = 0 ! move the gridpoints using the grid velocities that have just been read
         CALL DCOPY(NDIM*NP,CORG(1,1),1,CORG(1,NP+1),1) ! put grid at t=n into t=n+1
         CALL DAXPY(NDIM*NP,DELT,XYZDOT,1,CORG(1,NP+1),1) ! update the grid at t=n+1
      ENDIF
C
      IF(ALE_MOVE_GRID)THEN
C
C     Here we move the grid using ad-hoc subroutines, see grid_motion.f
C
C
         IF(IALE.NE.0)
     &   CALL MOVEXYZ(CORG(1,1+2*NP),CORG(1,NP+1),NDIM,NP,TIME) ! move the grid
C                     ^              ^
C                     |              |
C                     |              |
C
C     this is the grid * this is the grid
C     at t=0 i.e. the  * at time t=TIME i.e. n=n+1
C     "fixed" grid     *
C
C
C     compute the grid velocity as (x^{n+1}-x^n)/dt
C
         CALL DCOPY(NDIM*NP,CORG(1,NP+1),1,XYZDOT,1) ! XYZDOT := x^{n+1}
         CALL DAXPY(NDIM*NP,MONE,CORG(1,1),1,XYZDOT,1) ! XYZDOT := x^{n+1} - x^n
         CALL DSCAL(NDIM*NP,ONE/DELT,XYZDOT,1) ! XYZDOT := (x^{n+1}-x^n)/dt
C
         DO I = 1, NDIM
            avgdxyzdt(I) = DASUM(NP,XYZDOT(I),NDIM)/NP
            INDX = IDAMAX(NP,XYZDOT(I),NDIM)
            INDX = (INDX-1)*NDIM+I 
            maxdxyzdt(I) = ABS(XYZDOT(INDX))
         ENDDO
C
         WRITE(NOUT,*)
         WRITE(NOUT,*)' avg |b_x,b_y,b_z| = ',(avgdxyzdt(I),I=1,NDIM)
         WRITE(NOUT,*)' max |b_x,b_y,b_z| = ',(maxdxyzdt(I),I=1,NDIM)
         WRITE(NOUT,*)
         WRITE(IOALE,*)TIME,(avgdxyzdt(I),I=1,NDIM),
     &                      (maxdxyzdt(I),I=1,NDIM)
!
C
      ELSEIF(ALE_LAPLACE_GRID)THEN
C
C     Here we move the grid by solving Laplace's equation
C
         WRITE(GNAME(5:10),FMT="(I6.6)")ITIM
!
!        Move the grid: this is only needed to compute the grid velocity
!        the way the grid is move is hardwired into the code; see grid_motion.f
!
         IF(IALE.NE.0)
     &   CALL MOVEXYZ(CORG(1,1+2*NP),CORG(1,NP+1),NDIM,NP,TIME) ! move the grid
C                     ^              ^
C                     |              |
C                     |              |
C
C     this is the grid   this is the grid
C     at t=0 i.e. the    at time t=TIME i.e. n=n+1
C     "fixed" grid
C
C     compute the grid velocity as (x^{n+1}-x^n)/dt:
C     this is only needed in order to have the grid velocity on all boundaries
C
         CALL DCOPY(NDIM*NP,CORG(1,NP+1),1,XYZDOT,1) ! XYZDOT := x^{n+1}
         CALL DAXPY(NDIM*NP,MONE,CORG(1,1),1,XYZDOT,1) ! XYZDOT := x^{n+1} - x^n
         CALL DSCAL(NDIM*NP,ONE/DELT,XYZDOT,1) ! XYZDOT := (x^{n+1}-x^n)/dt
C
C     solve Laplace's Eqn to get the grid velocity at all gridpoints  
C
         CALL GridDeformation(MotionSolver,PetscMatC,NodalBcs,NDIM,
     2                        NOFVERT,NPOIN,NGHOST,NPNOD,NELEM,NBFAC,
     3                        NFACE)
!
!     compute the average grid velocity once all velocities have been re-computed inside GridDeformation
!
         DO I = 1, NDIM
            avgdxyzdt(I) = DASUM(NP,XYZDOT(I),NDIM)/NP
            INDX = IDAMAX(NP,XYZDOT(I),NDIM)
            INDX = (INDX-1)*NDIM+I 
            maxdxyzdt(I) = ABS(XYZDOT(INDX))
         ENDDO
         WRITE(NOUT,*)
         WRITE(NOUT,*)' avg |b_x,b_y,b_z| = ',(avgdxyzdt(I),I=1,NDIM)
         WRITE(NOUT,*)' max |b_x,b_y,b_z| = ',(maxdxyzdt(I),I=1,NDIM)
         WRITE(NOUT,*)
         WRITE(IOALE,*)TIME,(avgdxyzdt(I),I=1,NDIM),
     &                      (maxdxyzdt(I),I=1,NDIM)
!
!        we must now move the nodes again
!
         NITEMS = NDIM*NP
!
!        copy the grid at t=n into the location where the grid at t=n+1 will be stored
!
         CALL DCOPY(NITEMS,CORG(1,1),1,CORG(1,NP+1),1)
C                                           ^
C                                           |
C                                           |
C                                           |
C                        now     this is the grid
C                        at time t=TIME-DELT i.e. t=n*DELT
!
!        update using the grid velocity
!
         CALL DAXPY(NITEMS,DELT,XYZDOT,1,CORG(1,NP+1),1)
!
      ENDIF 
C
C    At this stage we have computed the mesh coordinates at time level n+1 i.e. t+dt
C    now CORG(1,NP+1) stores the grid at time level n+1 
C
C
      IF(DEBUG)THEN
         call print4check(corg,vol,vmedian,facnor,ndim,np,nelem,
     &time,'After updating the mesh coordinates ',iw)
      ENDIF
C
C
C
C
C
      IF((ITIM/IFAIL)*IFAIL.EQ.ITIM)THEN
            WRITE(GNAME(5:10),FMT="(I6.6)")ITIM
            WRITE(6,*)'Writing mesh at time ',TIME,' in file = ',
     &GNAME(1:12)
            CALL GRIDPTS(NPOIN,NGHOST,NPNOD,NELEM,NBFAC,NBINT,NHOLE,
     2             NOFVERT,NDIM,CORG(1,NP+1),ICELNOD,IBNDPTR, ! icn might use the periodic copy?
     3             GNAME(1:14),"noname.dat","w") ! write the NEW grid
C
          WRITE(FNAME(5:10),FMT="(I6.6)")ITIM
          CALL TECPLOT(FNAME(1:14),CORG(1,NP+1),NDIM,NP,ICELNOD,NOFVERT,
     &                 NELEM,XYZDOT,NDIM)
      ENDIF
CC
CC    LCORG       |---- n   ----| |---- n+1 ----| |---- 0 ----|
CC                                       ^
CC                                       |
CC                                       +-------------- we have updated this one ONLY
CC
CC    LMEDIAN     |---- n   ----| |---- XXX ----| |---- n ----|
CC    LVOL        |---- n   ----| |---- XXX ----| |----n-1----|
CC    LFACNOR     |---- n   ----|
CC
C
C
! compute the face normals using the updated grid coordinates
! i.e. those at time  "n+1"; this is actually required only to compute the volume
! at time n+1
!
      CALL FF(ICELNOD,ICELCEL,ICELFAC,NOFVERT,NELEM,CORG(1,NP+1),NDIM,
     +        NP,FACNOR,NFACE,NBFAC,NBINT)
caldo
!         CALL I4Mat_Print('General',' ',NOFVERT,NELEM,ICELFAC,
!    +    NOFVERT,'New ICELFAC array ',IFAIL)
!         CALL R8Mat_Print('General',' ',NDIM,NFACE,FACNOR,
!    +            NDIM,'New face normals ',IFAIL)
C
C
C     The volume of EACH triangular/tetrahedral element is computed and stored in VOL(NELEM+1:2*NELEM)
C     this is at time level n+1
C
      CALL CMPVOL(ICELNOD,ICELFAC,NOFVERT,NELEM,
     +            FACNOR,NFACE,CORG(1,NP+1),NDIM,NP,
     +            VOL(NELEM+1),TIME)
C
C     this is the area/vol of the median dual cell at time level n+1
C
      CALL MEDIAN(VMEDIAN(NP+1),NPOIN,NGHOST,NPNOD,NP,VOL(NELEM+1),
     &   ICELNOD,PMAP,NOFVERT,NELEM)
CC
CC    LCORG       |---- n   ----| |---- n+1 ----| |---- 0 ----|
CC    LMEDIAN     |---- n   ----| |---- n+1 ----| |---- n ----|
CC    LVOL        |---- n   ----| |---- n+1 ----| |----n-1----|
CC    LFACNOR     |---- n+1 ----|        ^
CC                       ^               |
CC                       +---------------+-------------- we have updated these three (LMEDIAN,LVOL,LFACNOR)
CC
      IF(DEBUG)THEN
         call print4check(corg,vol,vmedian,facnor,ndim,np,nelem,
     &time,'Before computing the mesh coordinates at tnow ',iw)
      ENDIF
C
C     compute the coordinates of the grid to be used in the residual calculation
C     at some intermediate time btw. n and n+1:
C     x^{n+a} = a*x^n + (1-a)*x^{n+1}
C
      TNOW = ALFALE*(TIME-DELT)+(ONE-ALFALE)*TIME
      CALL DSCAL(NDIM*NP,ALFALE,CORG(1,1),1)
      CALL DAXPY(NDIM*NP,(ONE-ALFALE),CORG(1,NP+1),1,CORG(1,1),1)
CC
CC    LCORG       |---- n+a ----| |---- n+1 ----| |---- 0 ----|
CC                       ^                
CC                       |                
CC                       +---------------+-------------- we have updated this one
CC    LMEDIAN     |---- n   ----| |---- n+1 ----| |---- n ----|
CC    LVOL        |---- n   ----| |---- n+1 ----| |----n-1----|
CC    LFACNOR     |---- n+1 ----|
CC
!
!     compute the face normals at this intermediate grid location (n+a);
!     the face normals are overwritten
!
      CALL FF(ICELNOD,ICELCEL,ICELFAC,NOFVERT,NELEM,CORG(1,1),NDIM,
     +        NP,FACNOR,NFACE,NBFAC,NBINT)
caldo
!         CALL I4Mat_Print('General',' ',NOFVERT,NELEM,ICELFAC,
!    +    NOFVERT,'New ICELFAC array ',IFAIL)
!         CALL R8Mat_Print('General',' ',NDIM,NFACE,FACNOR,
!    +            NDIM,'New face normals ',IFAIL)
C
C     The volume of EACH triangular/tetrahedral element is computed and stored in VOL(1:NELEM)
C     notice that we do NOT compute the area/volume of the median dual cells at this intermediate time
C
      CALL CMPVOL(ICELNOD,ICELFAC,NOFVERT,NELEM,
     +            FACNOR,NFACE,CORG(1,1),NDIM,NP,
     +            VOL(1),TNOW)
CC
CC    LCORG       |---- n+a ----| |---- n+1 ----| |---- 0 ----|
CC    LMEDIAN     |---- n   ----| |---- n+1 ----| |---- n ----|
CC    LVOL        |---- n+a ----| |---- n+1 ----| |---- n ----| |----n-1----|
CC                       ^                
CC                       |                
CC                       +---------------+-------------- we have updated this one
CC    LFACNOR     |---- n+a ----|
CC                       ^                
CC                       |                
CC                       +---------------+-------------- we have updated this one
CC
      IF(DEBUG)THEN
         call print4check(corg,vol,vmedian,facnor,ndim,np,nelem,
     &time,'Before leaving NEWGEO ',iw)
      ENDIF
C
C     VOL(1:NELEM) stores the areas/volumes at time n+a
C     VOL(NELEM+1:2*NELEM) stores the areas/volumes at time n+1
C     VOL(2*NELEM+3*NELEM) stores the areas/volumes at time n
C     VOL(3*NELEM+4*NELEM) stores the areas/volumes at time n-1
C
      IF( GCLCHK .AND. ((ITIM/IFAIL)*IFAIL.EQ.ITIM) )THEN
      WRITE(6,*)'Checking GCL'
      CALL CHKGCL(ICELNOD,ICELFAC,FACNOR,XYZDOT,VOL,CORG,
     2            NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NDIM,
     3            itim)
      ENDIF 
C
      RETURN
      END
