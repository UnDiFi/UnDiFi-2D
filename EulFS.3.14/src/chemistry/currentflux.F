      SUBROUTINE CURRENTFLUX(ICELNOD,ICELFAC,NOFVERT,NELEM,IBNDFAC,
     &NBFAC,FACNOR,XYZDOT,ZROE,NDIM,NOFVAR,NP,VOL,TIME)
C
      IMPLICIT NONE
C
C     $Id: currentflux.F,v 1.4 2013/09/17 10:06:31 abonfi Exp $
C
C     Here we compute the 
C
      INTEGER NOFVERT,NELEM,NBFAC,NDIM,NOFVAR,NP
      INTEGER ICELNOD(*),ICELFAC(*),IBNDFAC(3,NBFAC)
      DOUBLE PRECISION FACNOR(*),VOL(*),XYZDOT(*),ZROE(*)
      DOUBLE PRECISION TIME
C

      INCLUDE 'paramt.h'
      INCLUDE 'constants.h'
      include 'bnd.h'
      include 'plasma.h'
      include 'bnd.com'
      include 'streamplasma.com'
      include 'electric.com'
      INCLUDE 'io.com'
      INCLUDE 'time.com'
      INCLUDE 'three.com'
C
#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
      integer isum,source,dest,total,NPES,tag
#endif
C
      INTEGER ICN(MAXNOFVERT)
      DOUBLE PRECISION VCZ(MAXNOFVAR*MAXNOFVERT*MAXTIMLEVS)
      DOUBLE PRECISION VCB(3*MAXNOFVERT) 
      DOUBLE PRECISION VCN(3*MAXNOFVERT) 
      DOUBLE PRECISION WKSP(0:NCOLOR) 
C
      INTEGER IELEM,IBFAC
      DOUBLE PRECISION VOLUME(MAXTIMLEVS+1),ESQR,HELP,SIGMA,PRESS
      INTEGER IBC,IVERT
C
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
C
      DOUBLE PRECISION PRESS4Ar,DDOT,GETSIGMA
C
      CALL DINIT(NCOLOR+1,ZERO,CFLUX,1)
C
      DO IBFAC = 1,NBFAC
         IELEM = IBNDFAC(1,IBFAC)
         IVERT = IBNDFAC(2,IBFAC)
         IBC   = IBNDFAC(3,IBFAC)

         CALL CELPTR(IELEM,NELEM,ICELNOD,ICELFAC,VOL,ZROE,FACNOR,
     2               XYZDOT,NDIM,NOFVERT,NOFVAR,NP,ICN,VCZ,VCN,VCB,
     3               VOLUME)

         CALL LINEARIZE(IELEM,LALE,VCN,VCB,NDIM,NOFVERT,
     +VCZ,NOFVAR,VOLUME(1))
C
         PRESS = PRESS4Ar(NDIM,ZAVG)
         ESQR = DDOT(NDIM,GRAD_PARM(NOFVAR,1),NOFVAR,
     &                    GRAD_PARM(NOFVAR,1),NOFVAR)
         SIGMA = GETSIGMA(ZAVG,PRESS,ESQR)
         HELP = DDOT(NDIM,GRAD_PARM(NOFVAR,1),NOFVAR,
     &          VCN((IVERT-1)*NDIM+1),1)
         HELP = HELP * PHIREF/LREFP

         CFLUX(IBC) = CFLUX(IBC) - HELP
         
      ENDDO
C
#ifdef MPI
      CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, NPES, ifail)
      TAG = 50
      DEST = 0
      if (MY_PE .EQ. 0) then
          do source = 1, NPES-1
              call MPI_RECV(WKSP,NCOLOR+1,MPI_DOUBLE_PRECISION,source,
     +               tag, PETSC_COMM_WORLD, istatus, ifail)
              DO IBC = 0,NCOLOR
                 CFLUX(IBC) = CFLUX(IBC) + WKSP(IBC)
              ENDDO
          enddo
      else
          call MPI_SEND(CFLUX, NCOLOR+1, MPI_DOUBLE_PRECISION, dest,
     +          tag, PETSC_COMM_WORLD, ifail)
      endif
#endif
C
      IF(MY_PE.EQ.0)THEN
         write(IWEFLX,FMT=100)TIME,(CFLUX(IBC),IBC=1,10)
      ENDIF
  100 FORMAT(F12.6,10(1X,F12.6))
      RETURN
      END
