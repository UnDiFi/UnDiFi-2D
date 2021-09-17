      subroutine circum_av(a,ja,ia,degree,rcline,zroe,zcl,wksp,
     &                     nofvar,ncl,lda)
c
c     $Id: qb.F,v 1.8 2020/04/23 09:37:34 abonfi Exp $
c
#include "petsc/finclude/petscsys.h"
      use petscsys
C
      implicit none
C
      include 'paramt.h'
      include 'bnd.h'
      include 'constants.h'
      include 'periodic.com'
      include 'stream.com'
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
      integer nofvar,ncl,lda
      integer ja(lda,*),ia(ncl+1),degree(ncl)
      double precision a(lda,*),rcline(ncl),zcl(nofvar,ncl),
     &ZROE(nofvar,*),wksp(*)
c
c     rcline radial coord for the current c-line
c     zcl    parameter vector for the current c-line
c     ia(3,j) indices of the three meshpoints surrounding a node
c     on a c-line
c     a(3,j)  area coordinates of c-line vertex j
c     a(4,j)  theta coordinate of c-line vertex j
c
      integer j,j1,j2,k,il,jbgn,jend,ipoin,ivar,IFAIL
      double precision denom,ux,uy,uz,ur,ut,sint,cost,wght,theta,dth,
     &help
c
!     write(6,*)'Hello world! This is ia on # ',MY_PE,
!    &(ia(il),il=1,ncl+1)
!     write(6,*)'Hello world! This is deg on # ',MY_PE,
!    &(degree(il),il=1,ncl)
c
      do 1 il = 1, ncl ! loop over clines
         do 7 ivar = 1, nofvar
              zcl(ivar,il) = ZERO
    7    continue
         jbgn = ia(il)     ! first vertex for the current cline
         jend = ia(il+1)-1 ! last  vertex for the current cline
c        write(6,*)il,jbgn,jend
         if( jend .LT. jbgn )goto 1
c
c        the integration formula to average over a c-line is:
c
c        int = 1/2*[z(1)+z(2)]*(th(2)-th(1))+
c              1/2*[z(2)+z(3)]*(th(3)-th(2))+
c              ..........
c              1/2*[z(n-1)+z(n)]*(th(n)-th(n-1))
c
c        which can be re-written as: 
c
c        int = 1/2*[z(1)*(th(2)-th(1))+z(2)*(th(3)-th(1))+..
c                   z(i)*(th(i+1)-th(i-1))+..+z(n)*(th(n)-th(n-1)]
c
c
         do 3 j = jbgn,jend ! loop over vertices of the current clines
         if    (j.EQ.jbgn)then
            j2 = j+1
            j1 = j
         elseif(j.EQ.jend)then
            j2 = jend
            j1 = jend-1
         else
            j2 = j+1
            j1 = j-1
         endif
c
c        we get the absolute value since the c-line might be traversed
c        for decreasing values of theta, but at the end we average by
c        the positive pitch
c
         dth = HALF*ABS(a(4,j2)-a(4,j1))
             do 4 ivar = 1, nofvar
                wksp(ivar) = ZERO
    4        continue
             theta = a(4,j)
             do 5 k = 1,3 ! loop over the vertices of the bndry face
                 ipoin = ja(k,j)
c     write( 6,*)ipoin,(zroe(ivar,ipoin),ivar=1,nofvar)
                 wght  =  a(k,j)
                 do 5 ivar = 1,nofvar
                    wksp(ivar) = wksp(ivar) + ZROE(ivar,ipoin)*wght
    5        continue
c
c     perform rotation 
c
!            wksp(1) = poutlet
             ux = wksp(2)
             uy = wksp(3)
             uz = wksp(4)
             cost = cos(theta)
             sint = sin(theta)
             ur = uy*cost+uz*sint
             ut =-uy*sint+uz*cost
             zcl(1,il) = zcl(1,il) + wksp(1) * dth
             zcl(2,il) = zcl(2,il) + ux * dth
             zcl(3,il) = zcl(3,il) + ur * dth
             zcl(4,il) = zcl(4,il) + ut * dth
    3    continue ! end loop over vertices of the current c-line
    1 continue ! end loop over c-line
#ifdef MPI
      CALL MPI_Allreduce(zcl,wksp,ncl*NOFVAR,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,PETSC_COMM_WORLD,IFAIL)
      CALL DCOPY(ncl*NOFVAR,wksp,1,zcl,1)
#endif
c
c
      call dscal(ncl*nofvar,ONE/pitch,zcl,1)
!     write(6,*)'p on all c-lines is ',(zcl(1,il),il=1,ncl)
c
      return
      end
c
      subroutine isre(rcline,zcl,nofvar,ncl)
c
c     purpose: compute pressure distribution over the clines
c              based on simple radial equilibrium
c
      implicit none
      include 'paramt.h'
      include 'bnd.h'
      include 'constants.h'
      include 'stream.com'
      INTEGER MY_PE
      COMMON/MPICOM/MY_PE
      integer nofvar,ncl
      double precision rcline(ncl),zcl(nofvar,ncl)
      double precision ut
      integer il,ilm1
c
c     rcline radial coord for the current c-line
c     zcl    parameter vector for the current c-line
c
      do 1 il = 1, ncl ! loop over clines
         ut = zcl(4,il)
         zcl(2,il) = ut*ut/ rcline(il)
    1 continue ! end loop over clines
!        write(6,*)'ut**2 is ',(zcl(2,il),il=1,ncl)
         zcl(1,1) = poutlet
      do 3 il = 2, ncl ! loop over clines
         ilm1 = il-1
         zcl(1,il) = zcl(1,ilm1) + HALF*(zcl(2,il)+zcl(2,ilm1))*
     &   (rcline(il)-rcline(ilm1)) 
!        zcl(1,il) = poutlet
    3    continue
!        write(6,*)(rcline(il),il=1,ncl)
!        write(16+MY_PE,*)(zcl(1,il),il=1,ncl)
!        write(6,*)'press is ',(zcl(1,il),il=1,ncl)
c
      return
      end
