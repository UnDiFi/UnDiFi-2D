!> \par Purpose
!>
!> Compute the ALE flux balance over a cell by calculating the integral of the ALE flux:
!>
!> \f[
!> F_n^{ALE} = - \left(\mathbf{b}\cdot\mathbf{n}\right)
!> \left(\begin{array}{c}
!> \rho \\\ \rho E \\\ \rho \mathbf{u} \end{array} \right)
!> = - \frac{1}{2} \left(\mathbf{b}\cdot\mathbf{n}\right)
!> \left( \frac{\partial U}{\partial Z} \right) Z
!> \f]
!> over the \c d+1 sides/faces of the simplex
!>
!> the calculation of the ALE flux balance is needed when using the ALE approach that is alternative to
!> Deconinck's approach
!>
!> Assuming that \f$Z\f$ and \f$\mathbf{b}\f$ vary linearly in space and since \f$\frac{\partial U}{\partial Z}\f$
!> ia also linear in \f$Z\f$, the ALE flux is a cubic function of the space coordinates and can therefore
!> be integrated exactly using Simpson's rule
!>
!> Performing the integration in 2D, we get:
!> \f[
!> I = I_{12} + I_{23} + I_{31} = \frac{1}{24} \, \sum_{i=1}^{3} \left(\frac{\partial U}{\partial Z}\right)_i \left[ \sum_{j=1}^{3} c_{ij} Z_j \right]
!> \f]
!>
!> where the entries \f$c_{ij}\f$ of the symmetric matrix \f$C\f$ are:
!> 
!> \f[
!> C =
!> \left(
!> \begin{array}{ccc}
!> 3\left(b_1^{n_3} + b_1^{n_2} \right) + b_2^{n_3} + b_3^{n_2} & b_1^{n_3} + b_2^{n_3} & b_3^{n_2} + b_1^{n_2}  \\\
!> c_{12} & 3\left(b_2^{n_3}+b_2^{n_1}\right) + b_1^{n_3} + b_3^{n_1} & b_2^{n_1} + b_3^{n_1} \\\
!> c_{13} & c_{23} & 3\left(b_3^{n_1}+b_3^{n_2}\right) + b_2^{n_1} + b_1^{n_2}
!> \end{array} \right) \quad = \quad
!> \left(
!> \begin{array}{ccc}
!> -2\,b_1^{n_1} + c_{12} + c_{13} & b_1^{n_3} + b_2^{n_3} & b_3^{n_2} + b_1^{n_2}  \\\
!> c_{12} & -2\,b_2^{n_2} + c_{21} + c_{23} & b_2^{n_1} + b_3^{n_1} \\\
!> c_{13} & c_{23} & -2\,b_3^{n_3} + c_{13} + c_{23}
!> \end{array} \right)
!> \f]
!> The following notation has been used:
!> \f[
!> b_i^{n_j} = \mathbf{b}_i\cdot\mathbf{n}_j
!> \f]
!> where the subscripts address the vertices of a cell; we have also used the fact that:
!> \f[
!> \sum_j \mathbf{n}_j = 0 \quad \rightarrow \quad \sum_j b_i^{n_j} = 0
!> \f]
!>
!> @param[in] NDIM the dimension of the space
!> @param[in] NOFVERT the nof vertices of the  simplicial element
!> @param[in] NDOF the nof dofs
!> @param[in] VCN the NDIM cartesian component of the (NOFVERT) inward face normals, scaled by its measure
!> @param[in] VCB the NDIM cartesian component of the grid velocity
!> @param[in] VCZ the NDOF dofs of the dependent variable within the NOFVERT vertices of the cell
!> @param[out] FLUXALE is the ALE flux through the face
!> @param[in] TEST set to .TRUE. when using the \c -check runtime option
C
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2020/03/28 09:51:15 $
!> \bug This routine does NOT work yet in 3D
C
      subroutine aleflux(ndim,nofvert,ndof,vcn,vcb,vcz,fluxale,test)
C
C
C     $Id: flxale.f,v 1.5 2020/03/28 09:51:15 abonfi Exp $
C
C
      implicit none
      include 'paramt.h' 
      include 'constants.h' 
      include 'three.com' 
      integer ndim,nofvert,ndof
      double precision b(3,3),c(3,3)
      double precision vcn(ndim,nofvert),vcb(ndim,nofvert)
      double precision vcz(ndof,nofvert),fluxale(ndof)
c
      double precision dudz(max_nofvar_sqr*maxnofvert)
      double precision work(maxnofvar),flux2(maxnofvar),flux1(maxnofvar)
      double precision f(maxnofvar,maxnofvert),divb,help
      double precision zmid(maxnofvar),bmid(maxnofvar)
      double precision db,dz
      double precision ddot,dasum
      integer i,j,k,ivar,ordsqr,iaddr,ifail,ng,n,irule
!
!  http://en.wikipedia.org/wiki/Gaussian_quadrature#Computation_of_Gaussian_quadrature_rules
!
      parameter(ng=3,irule=1)
!     parameter(ng=2)
!     parameter(ng=1)
      double precision gp(ng),weight(ng),zg(MAXNOFVAR,ng),bg(3,ng)
      integer jcycl
      logical test,verbose
      parameter(verbose=.FALSE.)
      character*12 label
c
      if(ng.EQ.2)then
         gp(1) = 1.d0/sqrt(3.d0)
         gp(2) = -gp(2)
         weight(1) = ONE
         weight(2) = ONE
      elseif(ng.EQ.3)then 
         gp(1) = 0.d0
         gp(2) = sqrt(3.d0/5.d0)
         gp(3) = -gp(2)
         weight(1) = 8.d0/9.d0
         weight(2) = 5.d0/9.d0
         weight(3) = 5.d0/9.d0
      endif
c
      ordsqr = ndof*ndof
      IF(ndim.NE.2)call exit(ndim)
c
!     build the b_{ij} = \mathbf{b}_i \cdot \mathbf{n}_j
c
      do j = 1,nofvert
         do i = 1,nofvert
            b(i,j) = ddot(ndim,vcb(1,i),1,vcn(1,j),1)
         enddo
      enddo
      IF(VERBOSE)
     &CALL R8Mat_Print('G',' ',Nofvert,Nofvert,b,nofvert,
     +      'B matrix ',IFAIL)
      do i = 1, nofvert
         iaddr = (i-1)*ordsqr +1
         call parm2cons(vcz(1,i),dudz(iaddr),ndof,ndim)
      enddo
      IF(VERBOSE)
     &CALL R8Mat_Print('G',' ',Ndof,Ndof*nofvert,dudz,ndof,
     +      'dUdZ matrix ',IFAIL)
c
c     build the C matrix
c
      c(1,2) = b(1,3)+b(2,3) 
      c(1,3) = b(3,2)+b(1,2) 
!     c(1,1) = 3.d0*(b(1,3)+b(1,2))
      c(1,1) = c(1,2)+c(1,3)+TWO*(b(1,3)+b(1,2))
c
      c(2,1) = c(1,2)
      c(2,3) = b(2,1)+b(3,1) 
!     c(2,2) = 3.d0*(b(2,3)+b(2,1))
      c(2,2) = c(2,1)+c(2,3)+TWO*(b(2,3)+b(2,1))
c
      c(3,1) = c(1,3)
      c(3,2) = c(2,3)
!     c(3,3) = 3.d0*(b(3,1)+b(3,2))
      c(3,3) = c(3,1)+c(3,2)+TWO*(b(3,1)+b(3,2))
c
c     call DSCAL(9,1.d0/24.d0,c,1)
      divb = (c(1,1) + c(2,2) + c(3,3))/3.d0
c
      IF(VERBOSE)THEN 
         CALL R8Mat_Print('G',' ',Nofvert,Nofvert,c,nofvert,
     +      'C matrix ',IFAIL)
         write(6,*)'div(b) = ',divb
      ENDIF
c
      call dinit(ndof,ZERO,fluxale,1)
      do i = 1, nofvert
         call dinit(ndof,ZERO,work,1)
         do j = 1, nofvert
            call daxpy(ndof,c(i,j),vcz(1,j),1,work,1)
         enddo ! j
         iaddr = (i-1)*ordsqr +1
         call dgemv('No',ndof,ndof,ONE/24.d0,dudz(iaddr),ndof,
     &              work,1,ONE,fluxale,1)
      enddo ! i
!     help = dasum(ndof,fluxale,1)
!
!     test should be TRUE when using -check as runtime option
!
      if(.NOT.test)return
   10 continue
c
c     here I test with either simpson, midpoint or trapezoidal, depending on the value of irule
c
      call dinit(ndof,ZERO,flux1,1) ! used for Gaussian quadrature
      call dinit(ndof,ZERO,flux2,1) ! used for either Simpson, midpoint or trapezoidal quadrature
c
      do i = 1,nofvert ! loop over the edges of the triangle
         call dinit(maxnofvar*maxnofvert,ZERO,f(1,1),1)
         j = jcycl(i+1)
         k = jcycl(i+2)
         do ivar = 1,ndof
            zmid(ivar) = HALF*(vcz(ivar,k)+vcz(ivar,j)) ! Z at midpoint
         enddo
         do ivar = 1,ndim
            bmid(ivar) = HALF*(vcb(ivar,k)+vcb(ivar,j)) ! b at midpoint
         enddo
         call intale(vcz(1,j),ndof,vcb(1,j),vcn(1,i),ndim,f(1,1))
         call intale(vcz(1,k),ndof,vcb(1,k),vcn(1,i),ndim,f(1,2))
         call intale(zmid    ,ndof,bmid    ,vcn(1,i),ndim,f(1,3))
!     Simpson's rule
         if(irule.EQ.1)then
!        write(6,*)"Using Simpson's rule "
            label = "Simpson's  "
            call daxpy(ndof,ONE/6.d0,f(1,1),1,flux2,1)
            call daxpy(ndof,ONE/6.d0,f(1,2),1,flux2,1)
            call daxpy(ndof,4.d0/6.d0,f(1,3),1,flux2,1)
         elseif(irule.EQ.2)then
!     trapezoidal rule
!           write(6,*)"Using Trapeziodal rule "
            label = "Trapeziodal"
            call daxpy(ndof,ONE/2.d0,f(1,1),1,flux2,1)
            call daxpy(ndof,ONE/2.d0,f(1,2),1,flux2,1)
         elseif(irule.EQ.3)then
!     midpoint rule
!     write(6,*)"Using Midpoint rule "
            label = "Midpoint   "
            call daxpy(ndof,ONE,f(1,3),1,flux2,1)
         endif
C
C        Use Gaussian quadrature and put it into flux1
C
c
c        Z and b at Gauss points
c
         do n = 1,ng 
            do ivar = 1,ndof
               dz = vcz(ivar,k)-vcz(ivar,j)
               zg(ivar,n) = HALF*dz*gp(n)+zmid(ivar)
!        write(6,*)'Gauss po ',n,' zg = ',zg(ivar,n),zmid(ivar)
            enddo
            do ivar = 1,ndim
               db = vcb(ivar,k)-vcb(ivar,j)
               bg(ivar,n) = HALF*db*gp(n)+bmid(ivar)
!        write(6,*)'Gauss po ',n,' bg = ',bg(ivar,n),bmid(ivar)
            enddo
         enddo
         do n = 1,ng
            call intale(zg(1,n),ndof,bg(1,n),vcn(1,i),ndim,f(1,1))
            call dscal(ndof,HALF*weight(n),f(1,1),1)
!        do ivar = 1,ndof
!        write(6,*)'Gauss po ',n,' flx(',ivar,') = ',f(1,1),f(1,3)
!        enddo
            call daxpy(ndof,ONE,f(1,1),1,flux1,1)
         enddo ! end loop over gauss points
      enddo ! loop over vertices (i)
      write(6,*)'ALE flux =  analytical,',label,'  Gauss quadrature'
      do i = 1,ndof
         write(6,*)'flux 2; i = ',i,fluxale(i) ,flux2(i),flux1(i)
!    &fluxale(i)/flux2(i)
c
c        copy the integral (computed using either Simpson or Midpoint or Trapezoidal) into the ALE flux
c
         fluxale(i) = flux2(i)
      enddo
!     pause
      return
      end
!> \par Purpose
!>
!> Compute the ALE flux at a point in space in the direction \f$\mathbf{n}\f$
!>
!> The ALE flux is computed as follows:
!>
!> \f[
!> F_n^{ALE} = - \left(\mathbf{b}\cdot\mathbf{n}\right)
!> \left(\begin{array}{c}
!> \rho \\\ \rho E \\\ \rho \mathbf{u} \end{array} \right)
!> = - \frac{1}{2} \left(\mathbf{b}\cdot\mathbf{n}\right)
!> \left( \frac{\partial U}{\partial Z} \right) Z
!> \f]
!>
!>
!>
!> @param[in] Z the NDOF dofs of the dependent variable
!> @param[in] NDOF the nof dofs
!> @param[in] B the NDIM cartesian components of the grid velocity
!> @param[in] VN the NDIM cartesian component of the inward face normal, scaled by its measure
!> @param[in] NDIM the dimension of the space
!> @param[out] FLUX is the ALE flux \f$F_n^{ALE}\f$
C
!> \author $Author: abonfi $
!> \version $Revision: 1.5 $
!> \date $Date: 2020/03/28 09:51:15 $
C
      subroutine intale(z,ndof,b,vn,ndim,flux)
      implicit none 
      include 'paramt.h' 
      include 'constants.h' 
      integer ndof,ndim
      double precision b(ndim),vn(ndim),z(ndof),flux(ndof)
      double precision dudz(max_nofvar_sqr*maxnofvert)
      double precision alpha
      double precision ddot
      call parm2cons(z,dudz,ndof,ndim)
      alpha = ddot(ndim,vn,1,b,1) * HALF
         call dgemv('No',ndof,ndof,alpha,dudz,ndof,
     &              z,1,ZERO,flux,1)
      return
      end
