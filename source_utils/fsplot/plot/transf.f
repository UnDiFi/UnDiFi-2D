      SUBROUTINE parm_to_cons(ZROE)
C
      IMPLICIT NONE
C
      INCLUDE'constants'
      include'dim_flags'
      include'int_flags'
      include'mesh_i4'
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION ZROE(NOFVAR,1)
C
C     .. Local Scalars ..
C
      INTEGER INODE 
C
C     .. External Functions ..
C
      DOUBLE PRECISION	DDOT
      EXTERNAL	DDOT 
C
C     .. Executable Statements ..
C
      do inode = 1 , NPOIN
         ZROE(2,inode) = GINV * ZROE(1,inode) * ZROE(2,inode) + HALF *
     .  GM1OG * DDOT(DIM,ZROE(3,INODE),1,ZROE(3,INODE),1)
C
         IF(DIM.EQ.3)ZROE(5,inode) = ZROE(1,inode) * ZROE(5,inode)
         ZROE(4,inode) = ZROE(1,inode) * ZROE(4,inode)
         ZROE(3,inode) = ZROE(1,inode) * ZROE(3,inode)
         ZROE(1,inode) = ZROE(1,inode) * ZROE(1,inode)
      enddo
c
      return
      end

      SUBROUTINE cons_to_parm(ZROE)
C
      IMPLICIT NONE
C
c
c This routine transforms the conservative vector into the
c paramter vector overwriting the array Z.
c
      INCLUDE'constants'
      include'dim_flags'
      include'int_flags'
      include'mesh_i4'
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION	ZROE
      DIMENSION ZROE(NOFVAR,1)
C
C     .. Local Scalars ..
C
      INTEGER*4 INODE
C
C     .. External Functions ..
C
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT 
C
C     .. Executable Statements ..
C
      do inode = 1 , NPOIN
c
         if(ZROE(1,inode) .le. zero )then
           write(*,100)inode,ZROE(1,inode)!,(point(j,inode),j=1,3)
  100 format(1X,'Neg. density in node ',I5,4(2X,F7.3))
           CALL EXIT(1)
         endif
c
         ZROE(1,inode) = sqrt(ZROE(1,inode))
         ZROE(3,inode) = ZROE(3,inode) / ZROE(1,inode)
         ZROE(4,inode) = ZROE(4,inode) / ZROE(1,inode)
         IF(DIM.EQ.3)ZROE(5,inode) = ZROE(5,inode) / ZROE(1,inode)
         ZROE(2,inode) = ( ZROE(2,inode) - HALF * GM1OG *
     1 ( DDOT(DIM,ZROE(3,INODE),1,ZROE(3,INODE),1) ) ) * GAM /
     2   ZROE(1,inode)
      enddo
c
      return
      end
C
C ------------------------------+------------------------------
C
      SUBROUTINE Conserved2Primitive(ZROE)
C
      IMPLICIT NONE
C
      INCLUDE'constants'
      INCLUDE'mesh_i4'
      INCLUDE'dim_flags'
      INCLUDE'int_flags'
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION	ZROE(NOFVAR,1)
C
C     .. Local Scalars ..
C
      INTEGER*4 IPOIN
C
C     .. External Functions ..
C
      DOUBLE PRECISION	DDOT
      EXTERNAL	DDOT
C
C     .. Executable Statements ..
C
C       Conserved variables are : velocities and static Pressure
C       Primitive variables are : velocities and Total Pressure
C
      DO 1 IPOIN = 1 , NPOIN
	ZROE(2,IPOIN) = ZROE(2,IPOIN) + HALF * 
     &  DDOT(DIM,ZROE(3,IPOIN),1,ZROE(3,IPOIN),1)
    1 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE Primitive2Conserved(ZROE)
C
      IMPLICIT NONE
C
      INCLUDE'constants'
      INCLUDE'mesh_i4'
      INCLUDE'dim_flags'
      INCLUDE'int_flags'
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION	ZROE(NOFVAR,1)
C
C     .. Local Scalars ..
C
      INTEGER	IPOIN
C
C     .. External Functions ..
C
      DOUBLE PRECISION	DDOT
      EXTERNAL	DDOT
C
C     .. Executable Statements ..
C
C       Primitive variables are : velocities and Total Pressure
C       Conserved variables are : velocities and static Pressure
C
      DO 1 IPOIN = 1 , NPOIN
	ZROE(2,IPOIN) = ZROE(2,IPOIN) - HALF * 
     &  DDOT(DIM,ZROE(3,IPOIN),1,ZROE(3,IPOIN),1)
    1 CONTINUE
C
      RETURN
      END
