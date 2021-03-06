!> \par Purpose
!>
!> subroutine GETDF2CORRDU computes matrix: 
!> \f[
!> \frac{ \partial F}{\partial U}
!> = \left( \begin{array}{cccc}
!>  0 & n_x & n_y & n_z  \\
!>  0 & u_n + u n_x & u n_y & u n_z  \\
!>  0 & v n_x & u_n + v n_y & v n_z  \\
!>  0 & w n_x & w n_y & u_n + w n_z
!> \end{array} \right)
!> \f]
!> where
!> \f$ F \f$ is given in subroutine INWLLI
!
!> @param[in] UVW the x,y,z components of the velocity
!> @param[in] VB the NDIM cartesian component of nodal grid velocity
!> @param[in] VN the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR nof dofs
!> @param[out] DFCORRDU the jacobian matrix
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2013/08/21 10:49:44 $
!> \warning NOFVAR should be set equal to NDIM+1 in the calling routine
!> \warning This routine has never been tested with a non-zero grid velocity
C     
      SUBROUTINE GETDF2CORRDU(UVW,VB,VN,NDIM,NOFVAR,DFCORRDU)
C
C     $Id: wall.f,v 1.6 2013/08/21 10:49:44 abonfi Exp $
C
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'constants.h'
C     ..
C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFCORRDU(NOFVAR,NOFVAR),UVW(NDIM),VN(NDIM),
     &VB(NDIM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION UDOTN
C     ..
C     We account for the grid velocity
C
      UDOTN = (UVW(1)-VB(1))*VN(1) + (UVW(2)-VB(2))*VN(2)
      IF (NDIM.EQ.3) UDOTN = UDOTN + (UVW(3)-VB(3))*VN(3)

      DFCORRDU(1,1) = ZERO
      DFCORRDU(1,2) = VN(1)
      DFCORRDU(1,3) = VN(2)

      DFCORRDU(2,1) = ZERO
      DFCORRDU(2,2) = UDOTN + UVW(1)*VN(1)
      DFCORRDU(2,3) = UVW(1)*VN(2)

      DFCORRDU(3,1) = ZERO
      DFCORRDU(3,2) = UVW(2)*VN(1)
      DFCORRDU(3,3) = UDOTN + UVW(2)*VN(2)

      IF (NDIM.EQ.3) THEN
          DFCORRDU(1,4) = VN(3)
          DFCORRDU(2,4) = UVW(1)*VN(3)
          DFCORRDU(3,4) = UVW(2)*VN(3)
          DFCORRDU(4,1) = ZERO
          DFCORRDU(4,2) = UVW(3)*VN(1)
          DFCORRDU(4,3) = UVW(3)*VN(2)
          DFCORRDU(4,4) = UDOTN + UVW(3)*VN(3)
      ENDIF

      RETURN

      END
!> \par Purpose
!>
!> subroutine GETDF4CORRDU computes matrix 
!> \f[
!> \frac{\partial F_{n}}{\partial Z} =
!>   \left( \begin{array}{ccccc}
!>  \sqrt{\rho} u_n - 2z_1b_n & 0 & z_1 n_x & z_1 n_y & z_1 n_z  \\
!>  -b_n\,z_2 & \sqrt{\rho}\left(u_n-b_n\right) & z_2 n_x & z_2 n_y & z_2 n_z  \\
!>  -b_n\,z_3 & 0 & \sqrt{\rho}\left(u_n-b_n\right)+ z_3 n_x & z_3 n_y & z_3 n_z  \\
!>  -b_n\,z_4 & 0 & z_4 n_x & \sqrt{\rho}\left(u_n-b_n\right) + z_4 n_y & z_4 n_z  \\
!>  -b_n\,z_5 & 0 & z_5 n_x & z_5 n_y & \sqrt{\rho}\left(u_n-b_n\right) + z_5 n_z
!> \end{array} \right)
!> \f]
!> where:
!> \f[
!> \sqrt{\rho}\left(u_n-b_n\right) =
!> \sqrt{\rho}\left(\mathbf{u}-\mathbf{b}\right)\cdot\mathbf{n} = z_3\,n_x+z_4\,n_y+z_5\,n_z-z_1\,\left(\mathbf{b}\cdot\mathbf{n}\right)
!> \f]
!> and
!> \f$ F \f$ is given in subroutine INWLL
!
!> @param[in] ZROE Roe's parameter vector in the gridpoint
!> @param[in] B the NDIM cartesian component of nodal grid velocity
!> @param[in] VN the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVAR nof dofs
!> @param[out] DFCORRDU the jacobian matrix
!
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2013/08/21 10:49:44 $
!> \warning NOFVAR should be set equal to NDIM+2 in the calling routine
C     
      SUBROUTINE GETDF4CORRDU(ZROE,B,VN,NDIM,NOFVAR,DFCORRDU)

      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'constants.h'

C     .. Scalar Arguments ..
      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(NDIM),VN(NDIM),ZROE(NOFVAR),
     &DFCORRDU(NOFVAR,NOFVAR)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION UDOTN,BDOTN
C     ..
      UDOTN = ZROE(3)*VN(1) + ZROE(4)*VN(2) ! UDOTN is \sqrt{\rho} (u\,n_x+v\,n_y+w\,n_z)
      BDOTN = B(1)*VN(1) + B(2)*VN(2) ! BDOTN is b_in_i i.e.\ b \cdot \mathbf{n}
      IF (NDIM.EQ.3) THEN
          UDOTN = UDOTN + ZROE(5)*VN(3)
          BDOTN = BDOTN + B(3)*VN(3)
      ENDIF
      UDOTN = UDOTN - BDOTN*ZROE(1) ! now UDOTN is \sqrt{\rho} [(u-b_x)n_x + etc

      DFCORRDU(1,1) = UDOTN-BDOTN*ZROE(1)
      DFCORRDU(1,2) = ZERO
      DFCORRDU(1,3) = VN(1)*ZROE(1)
      DFCORRDU(1,4) = VN(2)*ZROE(1)

      DFCORRDU(2,1) = -BDOTN*ZROE(2)
      DFCORRDU(2,2) = UDOTN
      DFCORRDU(2,3) = ZROE(2)*VN(1)
      DFCORRDU(2,4) = ZROE(2)*VN(2)

      DFCORRDU(3,1) = -BDOTN*ZROE(3)
      DFCORRDU(3,2) = ZERO
      DFCORRDU(3,3) = UDOTN + ZROE(3)*VN(1)
      DFCORRDU(3,4) = ZROE(3)*VN(2)

      DFCORRDU(4,1) = -BDOTN*ZROE(4)
      DFCORRDU(4,2) = ZERO
      DFCORRDU(4,3) = ZROE(4)*VN(1)
      DFCORRDU(4,4) = UDOTN + ZROE(4)*VN(2)

      IF (NDIM.EQ.3) THEN

          DFCORRDU(1,5) = VN(3)*ZROE(1)
          DFCORRDU(2,5) = ZROE(2)*VN(3)
          DFCORRDU(3,5) = ZROE(3)*VN(3)
          DFCORRDU(4,5) = ZROE(4)*VN(3)

          DFCORRDU(5,1) = -BDOTN*ZROE(5)
          DFCORRDU(5,2) = ZERO
          DFCORRDU(5,3) = ZROE(5)*VN(1)
          DFCORRDU(5,4) = ZROE(5)*VN(2)
          DFCORRDU(5,5) = UDOTN + ZROE(5)*VN(3)

      ENDIF

      RETURN

      END
C
!> \details
!> \par Purpose
!>
!>  INVISCID WALL boundary condition (INCompressible flows)
!>
!> add a correction flux having the form:
!> \f[ 
!> F = -\left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} ,
!> \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) \mathbf{u}
!> \right)
!> \f]
!
!> @param[in] NDIM the dimension of the space
!> @param[in] VNOR stores the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] B stores the NDIM cartesian component of nodal grid velocity
!> @param[in] Z is the vector of primitive variables (p,u,v,w)
!> @param[out] F the inviscid flux for an inviscid wall (eventually moving)
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2013/08/21 10:49:44 $
!
C     
      SUBROUTINE INVWLLI( NDIM , VNOR , B , Z , F )
C
C
      IMPLICIT NONE
C
      INTEGER NDIM
      DOUBLE PRECISION VNOR(NDIM),Z(*),B(NDIM),F(*)
      DOUBLE PRECISION VDOTN
C
C    | U_n         |    | 0   |
C    | U_n U + p n |    | p n |
C
C     We account for the possibility that the grid is moving....
C
      VDOTN      = VNOR(1)*(Z(2)-B(1)) + VNOR(2)*(Z(3)-B(2))
      IF(NDIM.EQ.3)VDOTN = VDOTN + VNOR(3)*(Z(4)-B(3))
C
      F(1) =-VDOTN
      F(2) =-VDOTN * Z(2)
      F(3) =-VDOTN * Z(3)
      IF(NDIM.EQ.3) F(4) =-VDOTN * Z(4)
C
      RETURN
      END
C
!> \details
!> \par Purpose
!>
!> INVISCID WALL boundary condition (Compressible)
!>
!> We add a correction in the form of a boundary flux:
!> \f[ 
!> F = -\left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} ,
!> \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) H,
!> \left( \rho \left( \mathbf{u}-\mathbf{b} \right) \cdot \mathbf{n} \right) \mathbf{u}
!> \right)
!> \f]
!> so that only the pressure flux is (or should be) left afterwords
!>
!> @param[in] NDIM the dimension of the space
!> @param[in] VNOR stores the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] B stores the NDIM cartesian component of nodal grid velocity \f$ \mathbf{b} \f$
!> @param[in] Z is Roe's parameter vector 
!> @param[out] F the inviscid flux for an inviscid wall (eventually moving)
!
!> \author $Author: abonfi $
!> \version $Revision: 1.6 $
!> \date $Date: 2013/08/21 10:49:44 $
C     
      SUBROUTINE INVWLL( NDIM , VNOR , B, Z , F )
C
C
      IMPLICIT NONE
C
      INTEGER NDIM
      DOUBLE PRECISION VNOR(NDIM),Z(*),B(NDIM),F(*)
      DOUBLE PRECISION VDOTN
C
C    | rho U_n         |    | 0   |
C    | rho U_n H       |  - | 0   |
C    | rho U_n U + p n |    | p n |
C
C
C     We account for the possibility that the grid is moving....
C
      VDOTN      = VNOR(1)*(Z(3)-Z(1)*B(1)) + VNOR(2)*(Z(4)-Z(1)*B(2))
      IF(NDIM.EQ.3)VDOTN = VDOTN + VNOR(3)*(Z(5)-Z(1)*B(3))
C
      F(1) =-VDOTN * Z(1)
      F(2) =-VDOTN * Z(2)
      F(3) =-VDOTN * Z(3)
      F(4) =-VDOTN * Z(4)
      IF(NDIM.EQ.3)F(5) =-VDOTN * Z(5)
C
      RETURN
      END
C
      SUBROUTINE GETDF4CORRDU4Ar(ZROE,VN,NDIM,NOFVAR,DFCORRDU)

      IMPLICIT NONE
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h' 
      INCLUDE 'dofs.com' 

C     .. Scalar Arguments ..

      INTEGER NDIM,NOFVAR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DFCORRDU(NOFVAR,NOFVAR),VN(NDIM),ZROE(NOFVAR)
      INTEGER I,J
!      INTEGER IE,IX,IY,IZ
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION UDOTN
C     ..
C
      UDOTN = ZROE(IX)*VN(1) + ZROE(IY)*VN(2)
      IF (NDIM.EQ.3) UDOTN = UDOTN + ZROE(IZ)*VN(3)
      DO I = 1 , NSP
        DO J = 1 , NSP
            DFCORRDU(I,J) = ZERO
            IF (I.EQ.J) THEN  
                DFCORRDU(I,J) = UDOTN            
            ENDIF
        ENDDO
      ENDDO  

      DO I = 1 , NSP    
        DFCORRDU(I,IE) = ZERO
        DFCORRDU(I,IX) = VN(1)*ZROE(I)
        DFCORRDU(I,IY) = VN(2)*ZROE(I)
      ENDDO

      DO J = 1 , NSP
        DFCORRDU(IE,J) = ZERO
        DFCORRDU(IX,J) = ZERO
        DFCORRDU(IY,J) = ZERO
      ENDDO

      DFCORRDU(IE,IE) = UDOTN
      DFCORRDU(IE,IX) = ZROE(IE)*VN(1)
      DFCORRDU(IE,IY) = ZROE(IE)*VN(2)

      DFCORRDU(IX,IE) = ZERO
      DFCORRDU(IX,IX) = UDOTN + ZROE(IX)*VN(1)
      DFCORRDU(IX,IY) = ZROE(IX)*VN(2)
      
      DFCORRDU(IY,IE) = ZERO
      DFCORRDU(IY,IX) = ZROE(IY)*VN(1)
      DFCORRDU(IY,IY) = UDOTN + ZROE(IY)*VN(2)

      IF (NDIM.EQ.3) THEN
          DO I=1,NSP
            DFCORRDU(I,IZ) = VN(3)*ZROE(I)
            DFCORRDU(IZ,I) = ZERO
          ENDDO
          DFCORRDU(IE,IZ) = ZROE(IE)*VN(3)
          DFCORRDU(IX,IZ) = ZROE(IX)*VN(3)
          DFCORRDU(IY,IZ) = ZROE(IY)*VN(3)
          
          DFCORRDU(IZ,IE) = ZERO
          DFCORRDU(IZ,IX) = ZROE(IZ)*VN(1)
          DFCORRDU(IZ,IY) = ZROE(IZ)*VN(2)
          DFCORRDU(IZ,IZ) = UDOTN + ZROE(IZ)*VN(3)
      ENDIF

      RETURN

      END
      SUBROUTINE INVWLL4Ar( NDIM , VNOR , Z , F )
C
C    .. INVISCID WALL boundary condition ..
C
      IMPLICIT NONE
C
      INCLUDE 'constants.h'
      INCLUDE 'plasma.h' 
      INCLUDE 'dofs.com' 
C
      INTEGER NDIM,ISP
      DOUBLE PRECISION VNOR(*),Z(*),F(*)
      DOUBLE PRECISION VDOTN
C
C
C    | rho U_n         |    | 0   |
C    | rho U_n H       |  - | 0   |
C    | rho U_n U + p n |    | p n |
C
      VDOTN      = VNOR(1)*Z(IX) + VNOR(2)*Z(IY)
      IF(NDIM.EQ.3)VDOTN = VDOTN + VNOR(3)*Z(IZ)
C
      DO ISP = 1 , NSP
        F(ISP) =-VDOTN * Z(ISP)
      ENDDO
      F(IE) =-VDOTN * Z(IE)
      F(IX) =-VDOTN * Z(IX)
      F(IY) =-VDOTN * Z(IY)
      IF(NDIM.EQ.3)F(IZ) =-VDOTN * Z(IZ)
C
      RETURN
      END
C
