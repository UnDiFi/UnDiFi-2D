      subroutine set_bc_laplace(np,xyz,ndim,zroe,ndof,nodecode)
!
!     $Id: set_bc_laplace.f,v 1.2 2013/06/04 16:08:13 abonfi Exp $
!
!
      implicit none
      include 'constants.h'
      integer np,ndim,ndof
      double precision xyz(ndim,*),zroe(ndof,*)
      integer nodecode(*)
      integer ipoin,ivar
      write(6,*)'Setting dirichlet bcs for laplace' 
      do ipoin = 1,np
         if( nodecode(ipoin) .NE. 0 )then
            zroe(ndof,ipoin) = xyz(1,ipoin)**2-xyz(2,ipoin)**2
         else
            zroe(ndof,ipoin) = ZERO
         endif
      enddo
      return 
      end 
