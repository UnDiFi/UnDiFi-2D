!> \par Purpose
!> 
!> Define some parameters to be used for statically dimensioned arrays
!> \verbatim 
!>    NMAX is the maximum number of variables currently allowed
!>    VMAX is the maximum number of vertices  currently allowed
!>    MAXNOFEQN is the maximum number of mean flow equations
!>    MAXNOFVAR same as NMAX
!>    MAXNOFVERT is the max number of vertices = d+1 (= VMAX)
!> \endverbatim 
!>    \author $Author: abonfi $
!>    \version $Revision: 1.5 $
!>    \date $Date: 2013/08/21 07:56:05 $
!> \warning MAXNOFVAR = MAXNOFEQN + no of turbulence eqn
!
!     $Id: paramt.h,v 1.5 2013/08/21 07:56:05 abonfi Exp $
!
      INTEGER NMAX,VMAX,MAXNOFEQN,MAXNOFVAR,MAXNOFVERT,MAX_NOFVERT_SQR
      INTEGER MAX_NOFVAR_SQR,MAXTIMLEVS
!
      PARAMETER(NMAX=8,VMAX=4,MAXNOFVAR=NMAX,MAXNOFEQN=NMAX)
      PARAMETER(MAXNOFVERT=VMAX,MAX_NOFVAR_SQR=MAXNOFVAR*MAXNOFVAR)
      PARAMETER(MAX_NOFVERT_SQR=MAXNOFVERT*MAXNOFVERT,MAXTIMLEVS=3)
!
!       
