!> \brief \b CHKFLX
!
!> \par Purpose
!>
!> \verbatim
!>
!> This routine checks the mass flux through the boundaries of the computational domain
!> \endverbatim
!>
!> @param[in] IBNDFAC Boundary informations: \c IBNDFAC(1,j) gives the element the j-th face belongs to; \c IBNDFAC(2,j) gives the local vertex number of element \c IBNDFAC(1,j) opposite the boundary face; \c IBNDFAC(3,j) gives the color of the boundary face
!> @param[in] NBFAC nof boundary faces
!> @param[in] ICELFAC Cell to face pointer: \c iface = \c ICELFAC(i,j) gives the global face number of the face opposite the i-th vertex of the j-th element; if \c iface >0 then the normal points inside the cell; otherwise it is the outward normal
!> @param[in] ICELNOD Cell to node pointer: \c ICELNOD(i,j) gives the global node number of the i-th vertex of the j-th element
!> @param[in] NOFVERT nof boundary faces
!> @param[in] NELEM nof boundary faces
!> @param[in] FACENORM Cartesian components of the normals to a face, multiplied by the face area
!> @param[in] ZROE Nodal values of the dependent variable;  compressible equations: \f$ Z = \sqrt{\rho} \left( 1, H, \mathbf{u} \right) \f$; incompressible equations: \f$ Z = \left( p, \mathbf{u} \right) \f$
!> @param[in] XYZDOT the cartesian components of the grid velocity
!> @param[in] NOFVAR nof dofs in ZROE
!> @param[in] NDIM dimension of the space
!> @param[in] NFACE nof faces in the mesh
!> @param[in] COMPRESSIBLE is .TRUE. when solving compressible flows
!> @param[in] LWRITE is .TRUE. when the boundary fluxes are to be written to a file
!> @param[in] ITER is current iteration number
C
C
!> \author $Author: abonfi $
!> \version $Revision: 1.16 $
!> \date $Date: 2021/03/28 16:27:21 $
C
      SUBROUTINE CHKFLX(IBNDFAC,NBFAC,ICELFAC,ICELNOD,NOFVERT,NELEM,
     &FACENORM,ZROE,XYZDOT,NOFVAR,NDIM,NFACE,COMPRESSIBLE,LWRITE,ITER)
C
C     $Id: chkflx.f,v 1.16 2021/03/28 16:27:21 abonfi Exp $
C
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      INCLUDE 'bnd.h'
      INCLUDE 'bnd.com'
      INCLUDE 'periodic.com'
      INCLUDE 'stream.com'
      INCLUDE 'time.com'
      INCLUDE 'io.com'
C
      INTEGER NBFAC,NOFVERT,NELEM,NDIM,NFACE,NOFVAR,ITER
      INTEGER NN,IOPT,IDOF
      PARAMETER(NN=NBTYPE+1)
      CHARACTER ERRMSG*72
      LOGICAL COMPRESSIBLE,VERBOSE,LWRITE
C
      INTEGER IBNDFAC(3,*),ICELFAC(NOFVERT,*),ICELNOD(NOFVERT,*)
      DOUBLE PRECISION FACENORM(NDIM,*),ZROE(NOFVAR,*),XYZDOT(NDIM,*)
C
      INTEGER I,IC,J,IBC,JVERT,NITEMS
      INTEGER IELEM,IBFAC,IFREQ,IVERT,IFACE,IPOIN
      DOUBLE PRECISION TEMP
C
      CHARACTER*12 FNAME
C
      DOUBLE PRECISION WKSP(MAXNOFVAR,0:NCOLOR),
     &FLXN(MAXNOFVAR*MAXNOFVERT),WORK(MAXNOFVAR),
     &VCZ(MAXNOFVAR*MAXNOFVERT),VCB(3*MAXNOFVERT)
      INTEGER IWKSP(0:NCOLOR),IV(VMAX)
C
      PARAMETER(VERBOSE=.FALSE.)
!     PARAMETER(VERBOSE=.TRUE.)
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION DDOT
      INTEGER ICYCL
      EXTERNAL         DDOT,ICYCL,FLUX2,FLUX4
C
      INTRINSIC MAX0
C
      DATA ERRMSG(1:7)/'BNDCHK '/
C
      IF(VERBOSE)
     &OPEN(100,FILE="bflux.log",FORM="FORMATTED",STATUS="unknown")
C
      IOPT = 1
      TEMP = RREF*UREF*LREF*LREF
C
      IF( PERIODIC_MESH .AND. ANNULAR )TEMP=TEMP*REAL(NBLADES)
C
      WRITE(NOUT,2000)MY_PE
2000  FORMAT(//' CHECKING BOUNDARY FLUXES ON PE #',I4,/' ',19('=')/)
C
      CALL DINIT(NOFVAR,ZERO,WORK,1)
      CALL DINIT(NOFVAR*NOFVERT,ZERO,FLXN,1)
      CALL DINIT(NOFVAR*NOFVERT,ZERO,VCZ,1)
      DO 14 J = 0 , NCOLOR
          IWKSP(J) = 0
          DO IDOF = 1,NOFVAR
             WKSP(IDOF,J) = ZERO
          ENDDO
   14 CONTINUE
C
C     each boundary face has a color, which is stored in IBNDFAC(3,*);
C     the option "-color" assigns a boundary type (see bnd.h) to each color
C
      IF(LWRITE)THEN
         FNAME = 'bflux000.dat'
         WRITE(FNAME(6:8),FMT="(I3.3)")MY_PE
         OPEN(UNIT=120,FILE=FNAME)
         NITEMS = 0
         DO I = 1, NBFLX(1) ! there are NBFLX(1) patches where the flux has to be computed
            IC = IBFLX(I)
            NITEMS = NITEMS + MCOLOR(IC)
         ENDDO
         WRITE(120,*)NITEMS
      ENDIF
C
      IFREQ = MAX0( 1 , NBFAC / 20 )
      IF(.NOT.LALE)THEN
          CALL DINIT(NOFVERT*NDIM,ZERO,VCB,1)
      ENDIF
      DO 16 IBFAC = 1 , NBFAC ! main loop over all boundary faces
         IBC = IBNDFAC(3,IBFAC)
         IWKSP(IBC) = IWKSP(IBC) + 1
C
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
C
         IFACE = ICELFAC(IVERT,IELEM)
         IF(IFACE.LT.0)STOP 'negative face in chkflx; should NOT occur'
         DO 12 JVERT = 1,(NOFVERT-1)
            IPOIN = ICELNOD(ICYCL(IVERT+JVERT,NOFVERT),IELEM+NELEM)
            IV(JVERT) = IPOIN
            CALL DCOPY(NOFVAR,ZROE(1,IPOIN),1,VCZ((JVERT-1)*NOFVAR+1),1)
            IF(LALE)THEN
               CALL DCOPY(NDIM,XYZDOT(1,IPOIN),1,
     &                   VCB((JVERT-1)*NDIM+1),1)
            ENDIF
   12    CONTINUE ! end loop over the vertices of cell IELEM
         IF(COMPRESSIBLE)THEN
            IF(NDIM.EQ.2)THEN
               CALL SIMPSON(NDIM,NOFVERT,NOFVAR,VCZ,VCB,
     &                   FACENORM(1,IFACE),FLXN,FLUX4)
            ELSEIF(NDIM.EQ.3)THEN
               CALL QUADRATURE(NDIM,NOFVERT,NOFVAR,VCZ,VCB,
     &                   FACENORM(1,IFACE),FLXN,FLUX4)
            ELSE
                    STOP 'Error with NDIM in CHKFLX'
            ENDIF
         ELSE ! incompressible
            IF(NDIM.EQ.2)THEN
               CALL SIMPSON(NDIM,NOFVERT,NOFVAR,VCZ,VCB,
     &                   FACENORM(1,IFACE),FLXN,FLUX2)
            ELSEIF(NDIM.EQ.3)THEN
               CALL QUADRATURE(NDIM,NOFVERT,NOFVAR,VCZ,VCB,
     &                   FACENORM(1,IFACE),FLXN,FLUX2)
            ELSE
                    STOP 'Error with NDIM in CHKFLX'
            ENDIF
         ENDIF
         IF(LWRITE)THEN
caldo               write(6,*)'nbflx(1) = ', nbflx(1) 
caldo               write(6,*)'ibflx = ', ibflx
            DO I = 1, NBFLX(1) ! there are NBFLX(1) patches where the flux has to be computed
               IC = IBFLX(I)
               IF(IC.EQ.IBC)THEN 
                  WRITE(120,*)(IV(JVERT),JVERT=1,(NOFVERT-1)),IBC
                  WRITE(120,*)(FLXN(J),J=1,NOFVAR)
               ENDIF
            ENDDO
         ENDIF
         CALL DAXPY(NOFVAR,ONE,FLXN((NOFVERT-1)*NOFVAR+1),1,
     &                                        WKSP(1,IBC),1)
C
         IF(VERBOSE)THEN
            WRITE(100,FMT=340)IBFAC,IBC,(ICELNOD(ICYCL(IVERT+JVERT,NOFVE
     &RT),IELEM+NELEM),JVERT=1,(NOFVERT-1))
            WRITE(100,FMT=345)(FLXN((NOFVERT-1)*NOFVAR+IDOF),IDOF=1,NOFV
     &AR)
         ENDIF
C
   16 CONTINUE ! main loop over all boundary faces
C
      DO 18 J = 0 , NCOLOR
         IF(IWKSP(J).NE.0)THEN ! does smthg iff there is at least a bndry face for color J
            CALL DAXPY(NOFVAR,ONE,WKSP(1,J),1,WORK,1)
            WRITE(NOUT,FMT=330)J,(WKSP(IDOF,J),IDOF=1,NOFVAR)!*TEMP
            WRITE(IFUNIT(J),FMT=*)ITER,(WKSP(IDOF,J),IDOF=1,NOFVAR)!*TEMP
         ENDIF
   18 CONTINUE
C
      WRITE(NOUT,FMT=335)(WORK(IDOF),IDOF=1,NOFVAR)
C
      IF(VERBOSE)CLOSE(100)
      IF(LWRITE)CLOSE(120)
C
      RETURN
C
C     I/O FORMATS
C
  330 FORMAT(10X,'INVISCID FLUX THROUGH BOUNDARY ',I2,6(1X,E12.6))
  335 FORMAT(10X,'INVISCID FLUX THROUGH',/,28X,'ALL BOUNDARIES ',
     &6(1X,E12.6))
  340 FORMAT(/,10X,'FLUXES THROUGH BOUNDARY ',I6,' IBC = ',I3,
     &' NODES ARE: ',3(1X,I6))
  345 FORMAT(10X,6(1X,E12.6))
C
      END
!> \brief \b FLUX4
!
!> \par Purpose
!>
!>    compute the Eulerian fluxes through the face of normal FACN from parameter vector
!> \f[ 
!> F = \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} ,
!> \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) H,
!> \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) \mathbf{u} + p \mathbf{n}
!> \right)
!> \f] 
!
!> @param[in] NDIM dimension of the space
!> @param[in] FACN the cartesian components of the scaled inward normal to the face
!> @param[in] ZROE Roe's parameter vector
!> @param[in] XYZDOT the cartesian components of the grid velocity
!> @param[out] FLXN the eulerian (inviscid) flux through the face
!>
C
!>    \author $Author: abonfi $
!>    \version $Revision: 1.16 $
!>    \date $Date: 2021/03/28 16:27:21 $
C
      SUBROUTINE FLUX4(NDIM,FACN,ZROE,XYZDOT,FLXN)
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INTEGER NDIM
C
      DOUBLE PRECISION FACN(NDIM),ZROE(*),XYZDOT(*),FLXN(*)
C
      DOUBLE PRECISION DOTM,PRES
      DOUBLE PRECISION PRESSC
C
C
      DOTM = FACN(1)*(ZROE(3)-ZROE(1)*XYZDOT(1))+
     &       FACN(2)*(ZROE(4)-ZROE(1)*XYZDOT(2))
      IF(NDIM.EQ.3)DOTM = DOTM + FACN(3)*(ZROE(5)-ZROE(1)*XYZDOT(3))
C
      PRES = PRESSC( NDIM , ZROE )
C
      FLXN(1) = ZROE(1)*DOTM
      FLXN(2) = ZROE(2)*DOTM
      FLXN(3) = ZROE(3)*DOTM+PRES*FACN(1)
      FLXN(4) = ZROE(4)*DOTM+PRES*FACN(2)
      IF(NDIM.EQ.3)FLXN(5) = ZROE(5)*DOTM+PRES*FACN(3)
C
      RETURN
      END 
!> \brief \b FLUX2
!
!> \par Purpose
!>
!>    compute the incompressible Eulerian fluxes through the face of normal FACN from primitive variables
!> \f[ 
!> F = \left( \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} ,
!> \left( \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) \mathbf{u} + p \mathbf{n}
!> \right)
!> \f] 
!
!> @param[in] NDIM dimension of the space
!> @param[in] FACN the cartesian components of the scaled inward normal to the face
!> @param[in] ZROE p,u,v,w
!> @param[in] XYZDOT the cartesian components of the grid velocity
!> @param[out] FLXN the eulerian (inviscid) flux through the face
!>
C
!>    \author $Author: abonfi $
!>    \version $Revision: 1.16 $
!>    \date $Date: 2021/03/28 16:27:21 $
C
      SUBROUTINE FLUX2(NDIM,FACN,ZROE,XYZDOT,FLXN)
      IMPLICIT NONE
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
C
      INTEGER NDIM
C
      DOUBLE PRECISION FACN(NDIM),ZROE(*),XYZDOT(*),FLXN(*)
C
      DOUBLE PRECISION DOTM,PRES
C
      DOTM = FACN(1)*(ZROE(2)-XYZDOT(1))+
     &       FACN(2)*(ZROE(3)-XYZDOT(2))
      IF(NDIM.EQ.3)DOTM = DOTM + FACN(3)*(ZROE(4)-XYZDOT(3))
C
      PRES = ZROE(1)
C
      FLXN(1) = DOTM
      FLXN(2) = ZROE(2)*DOTM+PRES*FACN(1)
      FLXN(3) = ZROE(3)*DOTM+PRES*FACN(2)
      IF(NDIM.EQ.3)FLXN(4) = ZROE(4)*DOTM+PRES*FACN(3)
C
      RETURN
      END 
