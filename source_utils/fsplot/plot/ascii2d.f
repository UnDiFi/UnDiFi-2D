      subroutine ascii2d(coor2,ndim,q2,nnu,nnodes,q2fr,nnufr,nnodesfr,
     +iconnec,nelem,nelemtypes,nufac,logfac,nfac)
c
      implicit none
c 
c NOTE In this example, there are only two variables (named q 
c      and qfr hereafter). q contains the set of variable in 
c      the whole field, it has nnu unknown for each grid point.
c      The variable qfr contains other variable on specific 
c      surface S of the mesh (inflow, outflow, wall,etc~...). 
c      It contains nnufr variables to be stored at the grid 
c      points of the surface S.
c       The unit number(s) below (unit1 and unit2) are arbitrary
c      and to be selected by the participant.
c
c  **************************************************************
c
      integer ndim,nnodes,nelem,nelemtypes,nnu,nfac,nnufr,nnodesfr
      real*8 coor2(ndim,nnodes)
      real*8 q2(nnu,nnodes)
      real q2fr(6,nnufr,nfac)
      integer nnperelement(6),iconnec(*),nufac(6,*)
      integer itypel,itype,nblock,npatch,i,l,nu,iunit,nb,nn,iel,iptr
      integer unit1,ifac
      integer logfac(*)
c
c parameters
c ----------
c
c ndimax               : space dimension (1,2 or 3)
c idimax,jdimax,kdimax : maximum of indices in the i,j,k direction,
c                        structured meshes
c idifrmax,jdifrmax,
c kdifrmax             : maximum of indices in the ifr,jfr,kfr direction,
c                      : structured mesh for S
c nnodesmax            : maximum number of nodes, unstructured meshes
c nnodesfrmax          : maximum number of nodes on S, unstructured
c                        meshes
c nelemax              : maximum number of element, ustructured meshes
c nnumax,
c nnufrmax             : maximum number of variables
c nfacmax              : maximum number of elements on S
c npatchmax            : maximum number of patches (structured meshes)
c 
c arrays
c ------
c coor1,coor2          : x,y,z coordinates 
c q1,q2,q1fr,q2fr      : array of variables
c nnperlement          : number of nodes per element
c ietyp                : define for a given element its type
c iconnec              : connectivity table
c iconnecp(nelemax)    : connectivity table pointer
c iconnecpfac(nfacmax) : connectivity table pointer for the frontier
c nufac(6,nfacmax)     : connectivity table for the frontier
c logfac(nfacmax)      : logic type of the frontier elements
c curv                 : values of variables in file unit2
c npoint               : number of value on the curves curv(*,*)
c ibegfr,...kendfr     : extreme indices of surfaces where surface
c                        values are asked for.
c
c note : the dimension 6*nelemax for iconnec has to be related to the
c size of nnperlement.
c
c   ***************************************************************
c
      unit1=22
      open (unit1,file='soldata',form='formatted')
c
      if(ndim.eq.2)then
         itypel = 3
      else
         itypel = 2
      endif
c
c     Number of blocks. Note that blocks should be either structured
c     or unstructured, but one database file can contain both
c     structured and unstructured blocks
c 
      nblock = 1
c
C     write (unit1) nblock
c
      do 103 nb=1,nblock
c
c     Read type of block, and number of unknowns stored for this block 
c
c     itype = 0: Unstructured mesh
c             1: Structured mesh
c     nnu   = number of unknowns stored
c     nnufr = number of unknown stored on the frontier
c
      itype =0
      npatch =0
c
c
c     ******************************************************************
c     Structured mesh
c     ******************************************************************
c
c
c     ******************************************************************
c     Unstructured mesh
c     ******************************************************************
c
c
c     nelemtypes = number of element types
c     nelem      = number of elements
c     nnodes     = number of nodes
c     nfac       = number of facets on the frontier
c     nnodesfr   = number of nodes on the frontier where data is stored
c                  they correspond to type of frontiers that are given
c                  by the next paragraph, problem per problem
c  remark : only on frontier is assumed to be of interest
c           same remarks as for structured
c
         write (unit1,FMT=*) 'points ',nnodes
         do 1 i=1,nnodes
            write (unit1,FMT=*) (sngl(coor2(l,i)),l=1,ndim)
    1 continue
         write (unit1,FMT=*) 'triangles ',nelem
c
c        write the connectivity table
c

            do 203 iel=1,nelem
               iptr = (iel-1)*3
               write (unit1,*) (iconnec(iptr+i)-1,i=1,3)
  203       continue
         write (unit1,FMT=*) 'scalars rho '
         do 3 i=1,nnodes
            write (unit1,FMT=*) sngl(q2(1,i))
    3 continue
            if(nnu.EQ.1)goto 103
         write (unit1,FMT=*) 'scalars rho*E '
         do 5 i=1,nnodes
            write (unit1,FMT=*) sngl(q2(2,i))
    5 continue
         write (unit1,FMT=*) 'scalars rho*u '
         do 7 i=1,nnodes
            write (unit1,FMT=*) sngl(q2(3,i))
    7 continue
         write (unit1,FMT=*) 'scalars rho*v '
         do 9 i=1,nnodes
            write (unit1,FMT=*) sngl(q2(4,i))
    9 continue
c
      close(unit1)
c
  103 continue
      close(unit1)
c
      return
      end
