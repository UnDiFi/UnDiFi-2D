      subroutine movep(xy,ndim,nbkgrd,npoin,npnod,zroe,wksp,ndof)
      implicit none
      integer ndim,nbkgrd,npoin,npnod,ndof
      double precision xy(ndim,*),zroe(ndof,*),wksp(*) 
      integer nitems
      integer i,j
C
C     +--------------NPOIN--------------+---NPNOD-----+
C     |--------NBKGRD----------|
C
      do i = 1,npoin+npnod
         write(36,*)i,(xy(j,i),j=1,ndim)
      enddo
      nitems = npnod*ndim
      call dcopy(nitems,xy(1,npoin+1),1,wksp,1) ! store periodic nodes into temporary location
      nitems = (npoin-nbkgrd)*ndim
      call dcopy(nitems,xy(1,nbkgrd+1),-1,xy(1,nbkgrd+npnod+1),-1)
      nitems = npnod*ndim
      call dcopy(nitems,wksp,1,xy(1,nbkgrd+1),1)
      do i = 1,npoin+npnod
         write(37,*)i,(xy(j,i),j=1,ndim)
      enddo
C
      nitems = npnod*ndof
      call dcopy(nitems,zroe(1,npoin+1),1,wksp,1) ! store periodic nodes into temporary location
      nitems = (npoin-nbkgrd)*ndof
      call dcopy(nitems,zroe(1,nbkgrd+1),-1,zroe(1,nbkgrd+npnod+1),-1)
      nitems = npnod*ndof
      call dcopy(nitems,wksp,1,zroe(1,nbkgrd+1),1)
      return
      end
