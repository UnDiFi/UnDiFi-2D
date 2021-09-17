!> \par Purpose
!>
!> Compute the integral of the conserved quantities over the entire computational domain \f$\Omega\f$:
!> the routine has been written assuming that parameter vector
!> is the dependent variable and therefore works for compressible flows only
!>
!> The integral is written to \c UNIT = \c IHST4
!>
!> In order to be consistent with the use of parameter vector as the set of dependent variables,
!> the mean value of the conservative variables should be computed as follows:
!> \f{eqnarray*}{
!> <U>\,|\Omega| &=& \int_{\Omega} U \, \mathrm{d}V \\\
!> &=& \frac{1}{2} \, \sum_{e=1}^{NE} \int_{T_e} \left(\frac{\partial U}{\partial Z}\right) \, Z \, \mathrm{d}V
!> = \frac{1}{2} \sum_{e=1}^{NE} \, \sum_{j \in e} \sum_{k \in e} \gamma_{jk}^e \left(\frac{\partial U}{\partial Z}\right)_{Z=Z_j} Z_k
!> \f}
!> In deriving the equation above we have used the fact that the following relation holds:
!> \f[
!> U = \frac{1}{2} \left(\frac{\partial U}{\partial Z}\right) \, Z
!> \f]
!> where both matrix \f$\partial U/\partial Z\f$ and vector \f$Z\f$ are linear functions of the independent variable \f$\mathbf{x}\f$.
!> the coefficients \f$\gamma_{jk}^e\f$ that appear in the above equation are:
!> \f[
!> \gamma^e_{jk} = |T_e| \times \left\{ \begin{array}{cc} \frac{1}{20} & j \ne k \\ \frac{1}{10} & j = k \end{array} \right.
!> \f]
!>
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] ICELFAC Cell to face pointer: \c ICELFAC(i,j) gives the global face number of the face opposite the i-th vertex of the j-th element
!> @param[in] VFACNOR Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] XYZDOT Cartesian components of the nodal grid velocities
!> @param[in] VOL area/volume of the simplicial elements (triangles,tetrahedra)
!> @param[in] ZROE is parameter vector
!> @param[in] NELEM nof boundary faces
!> @param[in] NPOIN nof interior nodes in the mesh
!> @param[in] NGHOST nof ghost nodes in the mesh
!> @param[in] NPNOD nof periodic nodes in the mesh
!> @param[in] NDIM dimension of the space
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NOFVAR is the nof dofs
!> @param[in] NTURB nof turbulent variables
!> @param[in] TIME physical time to be reached, i.e. time \f$(n+1)\Delta t\f$.
!>
!> \author $Author: abonfi $
!> \version $Revision: 1.3 $
!> \date $Date: 2013/09/17 10:01:18 $
!> \bug Only works for compressible (perfect gas) flows
!> \warning The integrals are area/volume multiplied; in a moving grid case (when the size of \f$\Omega\f$ changes with time) it might be more appropriate to divide by the current area/volume
!
      SUBROUTINE CHKINT(ICELNOD,ICELFAC,VFACNOR,XYZDOT,VOL,ZROE,
     3                  NELEM,NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,
     4                  NTURB,TIME)
C
      IMPLICIT NONE
C
C     $Id: chkint.f,v 1.3 2013/09/17 10:01:18 abonfi Exp $
C
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'conv.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'nloc.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
C
      INTEGER NPOIN,NGHOST,NPNOD,NDIM,NOFVERT,NOFVAR,NELEM,NTURB
      INTEGER ICELNOD(NOFVERT,NELEM),ICELFAC(NOFVERT,NELEM)
      DOUBLE PRECISION TIME
      DOUBLE PRECISION VFACNOR(NDIM,*),VOL(NELEM), ZROE(NOFVAR,*),
     2                 XYZDOT(NDIM,*)
C
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
C
      INTEGER IVAR,I,J,IADD,JADD,IELEM,NP
      INTEGER IFAIL,N4
C
C     ICN stores the vertices of the current element (0-based indexing)
C
C     ..
C     .. Local Arrays ..
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION VCN(3*MAXNOFVERT),VOLUME(MAXTIMLEVS+1),
     &WKSP(MAXNOFVAR,2)
      DOUBLE PRECISION VCZ(MAXNOFVAR*MAXNOFVERT),VCB(3*MAXNOFVERT)
      DOUBLE PRECISION DUDZ(MAX_NOFVAR_SQR*MAXNOFVERT)
C
      DOUBLE PRECISION WEIGHT,DIAG,OFFD
C
C     Some initializations ....
C
      NP = NPOIN + NGHOST + NPNOD
      OFFD = ONE/(REAL(NDIM+ONE)*REAL(NDIM+TWO))
      DIAG = TWO*OFFD
C
C     set work array equal to zero
C
      CALL DINIT(NOFVAR,ZERO,WKSP(1,1),1)
      CALL DINIT(NOFVAR,ZERO,WKSP(1,2),1)
C
      DO 2000 IELEM = 1,NELEM
C
         CALL DINIT(NOFVAR,ZERO,WKSP(1,2),1)
C
         CALL CELPTR(IELEM, NELEM, ICELNOD, ICELFAC, VOL, ZROE,
     +                VFACNOR, XYZDOT, NDIM, NOFVERT, NOFVAR, NP, ICN,
     3                VCZ, VCN, VCB, VOLUME)
C
         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +               VCZ,NOFVAR,VOLUME(1))
C
         DO I = 1,NOFVERT
            IADD = (I-1)*NOFVAR+1
            JADD = (I-1)*(NOFVAR*NOFVAR)+1
            CALL PARM2CONS(VCZ(IADD),DUDZ(JADD),NOFVAR,NDIM) ! dUdZ
         ENDDO ! end loop over vertices
C
C     compute the integral
C
         DO I = 1,NOFVERT
            IADD = (I-1)*(NOFVAR*NOFVAR)+1 ! (dU/dZ)_i
            DO J = 1,NOFVERT
               JADD = (I-1)*NOFVAR+1 ! Z_j
               IF(J.EQ.I)THEN
                  WEIGHT = VOLUME(1)*DIAG
               ELSE
                  WEIGHT = VOLUME(1)*OFFD
               ENDIF
               CALL DGEMV('No',NOFVAR,NOFVAR,WEIGHT,DUDZ(IADD),NOFVAR,
     &                    VCZ(JADD),1,ONE,WKSP(1,2),1) ! 
            ENDDO ! end loop over vertices
         ENDDO ! end loop over vertices
C
C     add cell integral to the global one
C
         CALL DAXPY(NOFVAR,HALF,WKSP(1,2),1,WKSP(1,1),1)
C
 2000 CONTINUE ! end loop over elements
      WRITE(IHST4,*)ITIM,ITER,TIME,(WKSP(J,1),J=1,NOFVAR)
      RETURN
      END
