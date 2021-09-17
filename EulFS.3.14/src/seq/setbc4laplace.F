      SUBROUTINE SETBC4LAPLACE(A,X,RHS,U,NodalBcs)
C
C     This routine SETs Boundary Conditions for Laplace's equation
C
C     $Id: setbc4laplace.F,v 1.5 2020/04/23 09:46:19 abonfi Exp $
C
C
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
      use petscvec
      use petscmat
      use petscis
C
      IMPLICIT NONE
C
C
      PetscBool IFLAG
C
      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      include 'bnd.h'
      include 'iset.h'
C
      INCLUDE 'bnd.com'
      INCLUDE 'dofs.com'
      INCLUDE 'flags.com'
      INCLUDE 'io.com'
      INCLUDE 'pfcgas.com'
      INCLUDE 'verbose.com'
      INCLUDE 'time.com'
C
      INTEGER NP,NDOF
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      INTEGER IFAIL
      INTEGER I,NI,ICLR,IPOIN,NOFVAR,IOFF
C
C
C
      DOUBLE PRECISION TIME(-1:1)
      DOUBLE PRECISION T,DV,HELP,s
      LOGICAL COMPRESSIBLE,INCOMPRESSIBLE,VISCOUS,LFLAG,PLASMA
     
C
      Mat A
      Vec X
      Vec RHS
      Vec U ! vector of conserved variables
      IS, dimension(0:*) :: NodalBcs
      integer     is_array(1)
      PetscOffset i_is,iaddr,jaddr
      PetscScalar x_array(1)
      PetscOffset i_x
      PetscScalar u_array(1)
      PetscOffset i_u
C
      DOUBLE PRECISION SETVOLTAGE
C
C     Executable Statements ..
C
      COMPRESSIBLE = (ABS(KAN).EQ.4)
      INCOMPRESSIBLE = (.NOT.COMPRESSIBLE)
      VISCOUS = (KAN.EQ.2).OR.(KAN.EQ.4)
      PLASMA = (KAN.EQ.-3)
C
      IF( .NOT. PLASMA )THEN
          WRITE(6,*)" Laplace's equation can only be set for plasma"
          CALL EXIT(1)
      ENDIF
      WRITE (NOUT,FMT=2000)
C
      TIME( 1) = ITIM*DELT
      TIME( 0) = TIME(1) - DELT
      TIME(-1) = TIME(0) - DELT
C
C     set voltage at times n-1,n,n+1 
C
      DO ICLR = 0, NCOLOR
         DO I = -1,1
            T = TIME(I)
            VOLTAGE(ICLR,I) = SETVOLTAGE(ICLR,T)
         ENDDO
!        write(6,*)iclr,(VOLTAGE(ICLR,I),i=-1,1)
      ENDDO
C
      call VecGetBlockSize(X,NOFVAR,ifail) 
#ifdef DEBUG
      call VecNorm(X,NORM_2,s,ifail) 
      write(6,*)'X norm before VecGetArray ',s
#endif
      call VecGetArray(X,x_array,i_x,ifail) ! vector of the unknowns
      call VecGetArray(U,u_array,i_u,ifail) ! vector of conserved variables
      DO 10 ICLR = 0, NCOLOR ! loop over all boundary patches
         IF( MCOLOR(ICLR) .EQ. 0 )GOTO 10 ! the boundary is EMPTY: do nothing
         IF( ICOLOR(ICLR,3) .EQ. BC_TYPE_SUPS_INLET )THEN ! an electrode if its boundary type is 3 
             call ISGetSize(NodalBcs(ICLR),NI,IFAIL)
             DV = VOLTAGE(ICLR,1) ! voltage for bndry patch ICLR at time level n+1
             WRITE (NOUT,FMT=2100)DV,ICLR,NI 
             call ISGetIndices(NodalBcs(ICLR),is_array,i_is,IFAIL)
C be careful: these index sets adress gridpoints, NOT unknowns and are 1-based indexed
             DO I = 1, NI
                IPOIN = is_array(i_is+I) ! these IS are 1-based indexed
                IOFF = IPOIN*NOFVAR
                iaddr = i_x + IOFF ! (IPOIN-1)*NOFVAR + NOFVAR ! address the location of the electric potential
                jaddr = i_u + IOFF ! (IPOIN-1)*NOFVAR + NOFVAR ! address the location of the electric potential
                x_array(iaddr) = DV-u_array(jaddr)
#ifdef DEBUG
                write(6,*)i,x_array(iaddr)
#endif
             ENDDO
             call ISRestoreIndices(NodalBcs(ICLR),is_array,i_is,IFAIL)
         ENDIF
   10 CONTINUE 
      call VecRestoreArray(X,x_array,i_x,ifail)
      call VecGetArray(U,u_array,i_u,ifail)
#ifdef DEBUG
      call VecNorm(X,NORM_2,s,ifail) 
      write(6,*)'X norm after VecRestoreArray ',s
#endif
C
C     Here we take advantage of the functionality of MatZeroRowsIS
C     to modify the rhs
C
      CALL MatZeroRowsIS(A,NodalBcs(Dirichlet4Poisson),ONE,X,RHS,IFAIL)
C
#ifdef DEBUG
      call VecNorm(X,NORM_2,s,ifail) 
      write(6,*)'X norm after MatZeroRowsIS ',s
      call VecNorm(RHS,NORM_2,s,ifail) 
      write(6,*)'RHS norm after MatZeroRowsIS ',s
C
      call ISGetSize(NodalBcs(Dirichlet4Poisson),NI,IFAIL)
      call ISGetIndices(NodalBcs(Dirichlet4Poisson),is_array,i_is,IFAIL)
      call VecGetArray(X,x_array,i_x,ifail)
      call VecGetArray(RHS,u_array,i_u,ifail)
             DO I = 1, NI
                iaddr = is_array(i_is+I) + 1 ! +1 because IS are 0-based indexed
                help = x_array(iaddr)
                x_array(iaddr) = DV
                write(6,*)i,is_array(i_is+I),x_array(i_x+iaddr),
     &          u_array(i_u+iaddr)
             ENDDO
      call ISRestoreIndices(NodalBcs(Dirichlet4Poisson),is_array,i_is,
     &                      IFAIL)
      call VecrestoreArray(X,x_array,i_x,ifail)
      call VecrestoreArray(RHS,u_array,i_u,ifail)
#endif
C
 2000 FORMAT (/,/," SETTING-BOUNDARY CONDITIONS 4 LAPLACE's EQN ",/,1X,
     &28('='),/)
 2100 FORMAT (5X,' SETTING VOLTAGE = ',F12.6,' ON BOUNDARY COLOURED ',
     &I2,' WITH ',I6,' POINTS')
C
C
C
      RETURN
      END
C
C
C
      DOUBLE PRECISION FUNCTION SETVOLTAGE(IBC,T)
      IMPLICIT NONE
      INCLUDE 'constants.h' 
      INCLUDE 'bnd.h' 
      INCLUDE 'bnd.com' 
      INCLUDE 'time.com'
      INCLUDE 'electric.com' 
      INTEGER IBC,IVOLT,TESTCASE
      DOUBLE PRECISION T
      DOUBLE PRECISION V,OMEGA,IC,VDIS,ALPHA,VG,BETA,BETAT,VG0      
C
C     TESTCASE identify the testcase
C     TESTCASE = 1 -> duct
C     TESTCASE = 2 -> electrodes
C
      TESTCASE = 1       
C 
C     IVOLT identify the time evolution law for the generator potential 
C     IVOLT = 1   ->   VG = CONST
C     IVOLT = 2   ->   VG = alpha*t if VG<1
C     IVOLT = 3   ->   VG = 1 + 0.5*SIN(OMEGA*T)
C       
      IVOLT = 1     
C    
C     The voltage should be hardwired here, can use CFLUX to get the current through
C     the electrodes
C
      IF(TESTCASE.EQ.1)THEN       
         IC=(ABS(CFLUX(3))+ABS(CFLUX(7)))/2
      ELSEIF(TESTCASE.EQ.2)THEN
         IC=(ABS(CFLUX(6))+ABS(CFLUX(7)))/2
      ELSE
         WRITE(6,*)'ERROR: TESTCASE must be 1 or 2'
         STOP
      ENDIF
C         
C     Generator potential setting
C
      VG0=1.0d0
      IF(IVOLT.EQ.1)THEN 
         VG = VG0
      ELSEIF(IVOLT.EQ.2)THEN
         BETA = 0.05d0
         BETAT = BETA*ITSTEP*DELT
         ALPHA = 1/BETAT     
         IF(T.lt.BETAT)THEN
            VG=ALPHA*T*VG0
         ELSEIF(T.ge.BETAT)THEN
            VG=VG0
         ENDIF
      ELSEIF(IVOLT.EQ.3)THEN
         OMEGA = TWO*PI/(ITSTEP*DELT)
         VG=VG0*(ONE + HALF*SIN(OMEGA*T))
       ELSE
         VG = VG0
      ENDIF  
C
C     Ohm law for the external circuit
C
      VDIS = VG - IC*RC/PHIREF      
C
C     Boundary condition setting 
C     (5,6,7,8) -> Immersed electrodes
C     (3,7) -> Duc
C     
      IF(TESTCASE.EQ.1)THEN
           IF(IBC .EQ. 7)THEN
               V = ONE
            ELSEIF(IBC .EQ. 3)THEN
               V = ONE + VDIS
            ELSE
               V = ZERO
            ENDIF
       ELSEIF(TESTCASE.EQ.2)THEN
            IF( (IBC .EQ. 5).OR.(IBC .EQ. 6) )THEN
               V = ONE
            ELSEIF( (IBC .EQ. 7).OR.(IBC .EQ. 8) )THEN
               V = ONE + VDIS
            ELSE
               V = ZERO
            ENDIF
       ENDIF


      SETVOLTAGE = V

      RETURN
      END

