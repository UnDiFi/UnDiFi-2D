! Shock point redistribution procedure
   
      subroutine rdstrshpnt(xysh,zold,shpnt,dx1,
     &                     rxysh,znew,rshpnt,nmax,ish)

      implicit none
      include 'paramt.h'

      integer shpnt,rshpnt,nmax
      double precision dx1,dx,s,sr,alpha,beta,ds,dsj
      dimension s(1000),sr(1000)
      double precision zold(ndof,*),znew(ndof,*)
      double precision rxysh(ndim,*),xysh(ndim,*)
      double precision kk
      integer i,j,ivar,ish,nn,mm

      dx=dx1

!     compute curvilinear abscissa of the original vector
      s(1)=0.
      do i=2,shpnt
        ds=(xysh(1,i)-xysh(1,i-1))**2+(xysh(2,i)-xysh(2,i-1))**2
        ds=sqrt(ds)
        s(i)=s(i-1)+ds
!       write(*,*)i,s(i)
      enddo

!     compute curvilinear abscissa of the redistributed vector
!     compute the number of redistributed points
      rshpnt=s(shpnt)/dx+1
!ren cccccccccccccccccccccc
!     the shock points of a shock must be always odd
      rshpnt=rshpnt/2
      rshpnt=rshpnt*2+1
      if(rshpnt.lt.3)rshpnt=3
!ren cccccccccccccccccccccc

      if(2*rshpnt.gt.npshmax)then
         write(6,*)' too many shock points! increase npshmax in paramt
     & up to ',2*rshpnt,' at least'
         call exit(1)
      endif
      if(rshpnt.gt.nmax)then
         write(6,*)' subr. rdstrshpnt must be called with nmax >= ',
     &rshpnt
         call exit(1)
      endif

!     compute the redistribution step
      dx=s(shpnt)/(rshpnt-1)
      sr(1)=0.
      do i=2,rshpnt
!       sr(i)=sr(i-1)+dx+dx/4.d0*(-1)**i
        sr(i)=sr(i-1)+dx
!       write(*,*)i, sr(i)
      enddo

!     if(ish.eq.3)then
!       rshpnt=rshpnt+1
!       sr(rshpnt)=sr(rshpnt-1)
!       sr(rshpnt-1)=sr(rshpnt-2)+dx*.70
!     endif

!     if(ish.eq.3)then
!!       mm=10
!!       nn=4
!!       kk=1.30242240774664
!        mm= 8
!        nn=2
!        kk=1.99196419704968
!        rshpnt=rshpnt+(mm-nn)
!        do i=1,mm
!         sr(rshpnt-mm+i)=sr(rshpnt-mm+i-1)+dx/(kk)**(i-1)
!         write(*,*)sr(rshpnt-mm+i)
!      enddo
!     endif

!     if(ish.eq.4)then
!       rshpnt=rshpnt+1
!       sr(rshpnt)=sr(rshpnt-1)
!       sr(rshpnt-1)=sr(rshpnt-2)+dx*.85
!     endif

!      if(ish.eq.4)then
!!       mm=10
!!       nn=4
!!       kk=1.30242240774664

!        mm= 8
!        nn=2
!        kk=1.99196419704968

!        rshpnt=rshpnt+(mm-nn)
!        do i=1,mm
!         sr(rshpnt-mm+i)=sr(rshpnt-mm+i-1)+dx/(kk)**(i-1)
!         write(*,*)sr(rshpnt-mm+i)
!        enddo

!      endif

!     if(ish.eq.2)then
!       rshpnt=rshpnt+1
!       sr(rshpnt)=sr(rshpnt-1)
!       sr(rshpnt-1)=sr(rshpnt-2)+dx*.70
!     endif

!     if(ish.eq.2)then
!!       mm=10
!!       nn=4
!!       kk=1.30242240774664

!        mm= 8
!        nn=2
!        kk=1.99196419704968

!        rshpnt=rshpnt+(mm-nn)
!        do i=1,mm
!         sr(rshpnt-mm+i)=sr(rshpnt-mm+i-1)+dx/(kk)**(i-1)
!         write(*,*)sr(rshpnt-mm+i)
!        enddo
!
!      endif

!     interpolation of shock points in the new distribution
!     assignment of the first and last point
      rxysh(1,1)=xysh(1,1)
      rxysh(2,1)=xysh(2,1)

      rxysh(1,rshpnt)=xysh(1,shpnt)
      rxysh(2,rshpnt)=xysh(2,shpnt)

      do ivar = 1, ndof
         znew(ivar,1)       =zold(ivar,1)
         znew(ivar,rshpnt)  =zold(ivar,shpnt)
         znew(ivar,1+rshpnt)=zold(ivar,shpnt+1)
         znew(ivar,2*rshpnt)=zold(ivar,2*shpnt)
      enddo

!     compute internal points
      do i=2,rshpnt-1
!       search in the original distribution the values to use in the interpolation
        do j=2,shpnt
          if(sr(i).ge.s(j-1).and.sr(i).lt.s(j))then
            ds=s(j)-s(j-1)
            dsj=sr(i)-s(j-1)
            alpha=dsj/ds
            beta=1.0d0-alpha
            rxysh(1,i)=beta*xysh(1,j-1)+alpha*xysh(1,j)
            rxysh(2,i)=beta*xysh(2,j-1)+alpha*xysh(2,j)
            do ivar = 1, ndof
               znew(ivar,i)=beta*zold(ivar,j-1)+alpha*zold(ivar,j)
               znew(ivar,i+rshpnt)=beta*zold(ivar,j+shpnt-1)+
     &                            alpha*zold(ivar,j+shpnt)
            enddo
          endif
        enddo
      enddo

      return
      end
