      SUBROUTINE RHSBC6( NOFVAR, RHS, DT ) 
C
C     This routine sets the periodic boundary conditions in the
C     "EXPLICIT" residual ..
C
CCCCC#include "include/finclude/petsc.h"
#include "include/finclude/vec.h"
#include "include/finclude/is.h"
      use petscvec
      use petscis
C
      IMPLICIT NONE
C
      Vec RHS,DT
      VecScatter PeriodicScatter
C
      INCLUDE 'iset.com'
C
      INTEGER NOFVAR
C
C
      INTEGER IFAIL,na,nb

      CALL VecScatterCreate(RHS,PeriodicSetA,RHS,PeriodicSetB,
     +                      PeriodicScatter,IFAIL) 
C
C     nodal residuals of nodes in PeriodicSetB are updated
C     with nodal residuals from nodes in PeriodicSetA.
C
      CALL VecScatterBegin(RHS,RHS,ADD_VALUES,SCATTER_FORWARD,
     +                     PeriodicScatter,IFAIL) 
      CALL VecScatterEnd(RHS,RHS,ADD_VALUES,SCATTER_FORWARD,
     +                   PeriodicScatter,IFAIL) 
C
C     nodal residuals of nodes in PeriodicSetB are copied
C     onto nodal residuals in nodes in PeriodicSetA.
C
      CALL VecScatterBegin(RHS,RHS,INSERT_VALUES,SCATTER_REVERSE,
     +                     PeriodicScatter,IFAIL) 
      CALL VecScatterEnd(RHS,RHS,INSERT_VALUES,SCATTER_REVERSE,
     +                   PeriodicScatter,IFAIL) 
      CALL VecScatterDestroy(PeriodicScatter,IFAIL)
C
C     done with residuals, do with timestep 
C
      CALL VecScatterCreate(DT,PeriodicNodesA,DT,PeriodicNodesB,
     +                      PeriodicScatter,IFAIL) 
C
C     nodal timestep of nodes in PeriodicNodesB are updated
C     with nodal residuals from nodes in PeriodicNodesA.
C
      CALL VecScatterBegin(DT,DT,ADD_VALUES,SCATTER_FORWARD,
     +                     PeriodicScatter,IFAIL) 
      CALL VecScatterEnd(DT,DT,ADD_VALUES,SCATTER_FORWARD,
     +                   PeriodicScatter,IFAIL) 
C
C     nodal timestep of nodes in PeriodicNodesB are copied
C     onto nodal residuals in nodes in PeriodicNodesA.
C
      CALL VecScatterBegin(DT,DT,INSERT_VALUES,SCATTER_REVERSE,
     +                     PeriodicScatter,IFAIL) 
      CALL VecScatterEnd(DT,DT,INSERT_VALUES,SCATTER_REVERSE,
     +                   PeriodicScatter,IFAIL) 
      CALL VecScatterDestroy(PeriodicScatter,IFAIL)
C
      RETURN
C
C     write(6,*)ipoin,ni,(y(i),i=1,ni),(ix(i),i=1,ni)
C     write(6,*)'VecSetValues has returned ',ifail
C     write(6,*)'VecAssemblyBegin has returned ',ifail
C     write(6,*)'VecAssemblyEnd has returned ',ifail
C     write(6,*)'Leaving SETBC1'
C
C     write(6,*)'VecScatterCreate has returned ',ifail
C     write(6,*)'VecScatterDestroy has returned ',ifail
C     write(6,*)'VecScatterEnd has returned ',ifail
C     write(6,*)'VecScatterBegin has returned ',ifail
C     write(6,*)'VecScatterEnd has returned ',ifail
C     write(6,*)'VecScatterBegin has returned ',ifail
C     CALL ISGetSize(PeriodicSetA,na,ifail)
C     CALL ISGetSize(PeriodicSetB,nb,ifail)
C     write(6,*)na,nb
C     CALL VecGetSize(RHS,nb,ifail)
C     write(6,*)nb
      END 
