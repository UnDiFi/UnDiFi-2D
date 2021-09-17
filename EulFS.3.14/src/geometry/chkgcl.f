!> \par Purpose
!>
!> Checks the GCL on a cell basis when the following runtime option is set:
!> \verbatim
!> -ale_check_gcl
!> \endverbatim
!> and writes some statistics to file \c volXXXXX.log where \c XXXXX is the iteration number.
!>
!> In particular, we dump to file:
!> 1. the current iteration \c ITER
!> 2. \f$T_e|^{n+1}-|T_e|^n\f$
!> 3. \f$ \left(\nabla \cdot \mathrm{b} \right) |T_e| \Delta t \f$
!> 4. The ratio btw. items 2 and 3 which should be 1 as long as Eq.(1) is satisfied to machine accuracy
!> 5. \f$ |T_e|^{n+\alpha}\f$
!> 6. \f$ |T_e|^{n+1}\f$
!> 7. \f$|T_e|^n\f$
!> 8. \f$|T_e|^{n-1}\f$
!>
!> To be more precise, we check that:
!> \f[
!> \frac{\mathrm{d}|T_e|}{\mathrm{d}t} = \frac{|T_e|^{n+1}-|T_e|^n}{\Delta t} = 
!> \nabla \cdot \mathbf{b}\, |T_e|^{n+\alpha} = \frac{1}{d} \sum_{j=1}^{d+1} \mathbf{b}_j \cdot \mathbf{n}_j^{n+\alpha} \quad\quad \mbox{(1)}
!> \f]
!> where:
!> \f[
!> \mathbf{b} = \frac{\mathbf{x}^{n+1}-\mathbf{x}^{n}}{\Delta t}
!> \f]
!> and the cell normals \f$\mathbf{n}_j^{n+\alpha}\f$ are computed using the geometry at some intermediate time:
!> \f[
!> \mathbf{x}^{n+\alpha} = \alpha \mathbf{x}^{n} + \left(1-\alpha\right) \mathbf{x}^{n+1}
!> \f]
!> where \f$\alpha = {\tt ALFALE}\f$ is a constant which is set using the runtime option
!> \verbatim
!> -ale_grid_weight ALFALE
!> \endverbatim
!> It may be readily verified that only the choice \f$\alpha = 0.5\f$ makes Eq.(1) an identity.
!>
!>
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j} gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] VFACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] XYZDOT Cartesian components of the nodal grid velocities
!> @param[in] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] ZROE nodal values of the dependent variable
!> @param[in] NELEM nof boundary faces
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in] NGHOST nof ghost nodes in the mesh
!> @param[in] NPNOD nof periodic nodes in the mesh
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NOFVAR nof dofs
!> @param[in] ITER physical time iteration counter
!> \author $Author: abonfi $
!> \version $Revision: 1.4 $
!> \date $Date: 2014/04/15 10:08:13 $
      SUBROUTINE CHKGCL(ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE,
     2                   NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     3                   ITER)
C
      IMPLICIT NONE
C
C     $Id: chkgcl.f,v 1.4 2014/04/15 10:08:13 abonfi Exp $
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'flags.com'
      INCLUDE 'nloc.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
C
      INTEGER NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NELEM,iter
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM), ZROE(NOFVAR,*),
     2                 XYZDOT(NDIM,*)
C
C
      integer I,J
 
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C
      INTEGER IELEM,NP
      INTEGER IFAIL
      DOUBLE PRECISION EPS,DVDT,DIVB
C
C     ..
C     .. Local Arrays ..
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION VCN(3*MAXNOFVERT),VOLUME(MAXTIMLEVS+1)
      DOUBLE PRECISION VCZ(MAXNOFVAR*MAXNOFVERT*MAXTIMLEVS),
     &VCB(3*MAXNOFVERT)
      DOUBLE PRECISION DIV,volnew,volold,tnew,tnow,told,taux
      CHARACTER*31 stringa
C
C
C     Some initializations ....
C
      stringa = "nodal coordinates of cell 00000"
      stringa = "vol00000.log"
C
      NP = NPOIN + NGHOST + NPNOD
C
      write(stringa(4:8),FMT="(I5.5)")iter
C
      tnew =  itim*delt
      tnow =  tnew-delt
      told =  tnow-delt
      taux = alfale*told+(ONE-ALFALE)*tnew
C
      OPEN(125,FILE=stringa(1:12))
      write(125,FMT=55)tnew,tnow,taux,tnew,tnow,told
C
      DO 2000 IELEM = 1,NELEM
C
          write(stringa(27:31),FMT="(I5.5)")ielem
C
          CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC,  VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
          volnew = VOLUME(2)
          volold = VOLUME(3)
          DVDT = VOLNEW-VOLOLD
          DIVB = DIV(NDIM,NOFVERT,VCN,VCB)*DELT
C
!         do i= 1,ndim
!            write(125,*)i,(vcb(i+(j-1)*ndim),j=1,nofvert)
!         enddo
          write(125,FMT=50)ielem,dvdt,divb,dvdt/divb,(volume(i),i=1,4) 
C
C
 2000 CONTINUE ! end loop over elements
      WRITE(6,*)
      CLOSE(125)
      RETURN
   55 FORMAT('Cell no.; V(',F11.4,')-V(',F11.4,'); div(b)*dt;  dvdt/div(
     &b)'4(1X,('V(',F11.4,')')))
   50 FORMAT(I6,7(1X,E12.5))
      END
