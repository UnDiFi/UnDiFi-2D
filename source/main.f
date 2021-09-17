      program undifi_2d

      implicit none

! ********************************************************************************************************************************
!  This program reads a triangle grid that already contains more than a shock, moves the shocks and re-grids using triangle
! ********************************************************************************************************************************

!>    ndim                        space dimension
!>    ndof                        number of degree of freedom
!>    nvt = ndim+1                number of vertices (=3) of a cell
!>    npoin                       number of nodes in the mesh
!>    nbfac                       number of boundary faces in the mesh
!>    nelem                       number of elements in the mesh
!>    nhole                       number of bodies(holes) in the mesh
!>    nholes                      number of holes used to remove cells in between the shock layers
!>    nedge                       number of edges in the mesh

!>    nshocks                     number of schocks and contact discontinuities
!>    nshockpoints(nshmax)        number of shock points for each shock
!>    nshocksegs(nshmax)          number of shock segments fore each shock(nshocksegs=nshockpoints-1)
!>    nphampoints                 number of phantom points. (points of the original mesh that are too near to the shock segments)

!>    typeshocks()                discontinuity type 's'=shock  'd'=contact discontinuity
!>    nspecpoints                 number of special points
!>    typespecpoints(nspmax)      type of special points
!>    shinspps(5,2,nspmax)        shocks in the special points
!>    ispclr(nspmax)              color of special point in other words, the special point moves along a specified colored boundary

!>    nshmax                      max number of shocks
!>    npshmax                     max number of shock points for each shock.
!>    neshmax                     max number of shock element for each shocks

!>    xysh(ndim,npshmax,nshmax)   coordinates of shock and discontinuity points
!>    xysh(ndim,npshmax,nshmax)   coordinates of shock and discontinuity points
!>    xyshu(ndim,npshmax,nshmax)  coordinates of shock and discontinuity points at upstream
!>    xyshd(ndim,npshmax,nshmax)  coordinates of shock and discontinuity points at downstream
!>    roeshu(ndof,npshmax,nshmax) upstream variable status of shocks and discontinuities
!>    roeshd(ndof,npshmax,nshmax) downstream variable status of shocks and discontinuties
!>    norsh(ndim,npshmax,nshmax)  normal unit vectors shocks and discontinuities
!>    wsh(npshmax,nshmax)         shock/discontinuity speed
!>    nodcodsh(npshmax,nshmax )   code characterizing the shock and discontinuity points
!                                 nodcodsh= -99 means that (i, ish) is not a shock/discontinuity point
!                                 nodcodsh= 10 means that (i, ish) is a shock/discontinuity point

!    .. parameters ..
      include 'paramt.h'
      integer*4 nin,nout
      parameter (nin=5,nout=6)
      integer nva
      parameter (nva=9990000)

!     .. array definitions
      double precision xysh  ( ndim,  npshmax, nshmax),
     +                 xyshu ( ndim,  npshmax, nshmax),
     +                 xyshd ( ndim,  npshmax, nshmax),
     +                 zroeshuold( ndof,  npshmax, nshmax),
     +                 zroeshdold( ndof,  npshmax, nshmax),
     +                 norsh ( ndim,  npshmax, nshmax),
     +                 wsh   ( ndim,  npshmax, nshmax)

!     arrays for unsteady predictor-corrector time accurate integration
      double precision xyshnew       ( ndim,  npshmax, nshmax),
     +                 norshnew      ( ndim,  npshmax, nshmax),
     +                 wshnew        ( ndim,  npshmax, nshmax),
     +                 wshmean       ( ndim,  npshmax, nshmax),
     +                 zroeshuoldnew ( ndof,  npshmax, nshmax),
     +                 zroeshdoldnew ( ndof,  npshmax, nshmax),
     +                 varray(ndim,30000) ! modify with the logic dstak/istak

      integer*4        nodcodsh(npshmax,nshmax),
     +                 nshocksegs(nshmax),
     +                 nshockpoints(nshmax),
     +                 shinspps(2,5,nspmax),
     +                 ispclr(5,nspmax)

      logical shtopolchanged
      logical neo, eulfs

      character        typespecpoints*5,
     +                 typeshocks*1

      dimension        typespecpoints(nspmax),
     +                 typeshocks(nshmax)

!     .. scalar definition
      INTEGER*4       nShocks,
     +                nPhamPoints,
     +                nSpecPoints

!     .. arrays in common ..
      double precision dstak(nva)
      character execmd*255,
     +          fname*255,
     +          fname2*255,
     +          fnameback*255,
     +          backdir*255,
     +          bindir*255,
     +          gastype*4,
     +          hostype*10,
     +          color1*2,
     +          color2*2,
     +          VELFILE*18,
     +          MODE,
     +          ISPREDICTOR

!     .. local scalars ..
      integer i,
     +        nshockpointsold(nshmax),ish,
     +        nholes,totshockpoints,ii,
     +        nvt,ifail,nsteps,nitems,nbegin

!     .. local arrays ..
      integer istak(1),lout(0:2)

!     pointers in 0 refer to the background mesh
      integer lbndfac(0:2),lcelcel(0:2),lcelnod(0:2),lcorg(0:2),
     & lnodcod(0:2),lzroe(0:2),ledgptr(0:2),lnodptr(0:2),
     & lshnor,lxyshold,lxyshnew,lwork,lpmap(0:2)
      integer nbfac(0:2),nelem(0:2),nhole(0:2),nedge(0:2),npoin(0:2),
     &nbpoin(0:2),nbfac_sh,npnod(0:2)
      integer lia(0:2),lja(0:2),liclr(0:2),nclr(0:2)
      logical fndbnds

!     .. external functions ..
      integer  initxdr,istkgt,istkst,system
!     external initxdr,istkgt,istkst,system
      external initxdr,istkgt,istkst
      double precision rand
      external rand

!     .. external subroutines ..
      external dinit,iinit,istkin,istkrl

!     .. common blocks ..
      common /cstak/dstak

!     .. equivalences ..
      equivalence (dstak(1),istak(1))

!     Time steps for predictor-corrector
      double precision dtpr, dtco, nowtime

!     Read command line arguments
      integer           :: no, n_args
      character(len=20) :: testcase
      character(len=20) :: args(5)
      logical           :: steady, unsteady

      n_args = command_argument_count();
      if (n_args /= 5) then
        write(*,*) 'Usage: ../../bin/UnDiFi-2D_x86_64
     +              0 501 false true "TestCaseName"'
        call abort()
      end if
      do i = 1, n_args
        call get_command_argument(i, args(i))
        args(i) = trim(adjustl(args(i)))
      end do
      read(args(1),*) nbegin
      read(args(2),*) nsteps
      read(args(3),*) eulfs
      read(args(4),*) steady
      read(args(5),*) testcase

      write(*,*) 'nbegin: ',    nbegin
      write(*,*) 'nsteps: ',    nsteps
      write(*,*) 'Use eulfs? ', eulfs
      write(*,*) 'Is steady? ', steady
      write(*,*) 'testcase: ',  testcase

!     flag to select the shock-capturing solver, eulfs or neo
      NEO = (.not. EULFS)

!     flag to select the type of simulation: steady or unsteady
      UNSTEADY = (.not. STEADY)

! --------- set character variables

      bindir  = "../../bin/"
      hostype = "x86_64"
!     hostype = 'i386'
      write(*,*)'1 ',hostype

!     gastype = "g140"
!     gastype = 'g167'

!     necessary for EulFS ALE
      velfile = "gridvel_000001.dat"
      mode    = 'w'

! Vale
!     flag used in fx_dps_loc.f for the shock reflection on a wedge
      ShTopolChanged=.false.
! Vale

!caldo
!      eulfs = .true.
!      !eulfs = .false.
!      neo = (.not. eulfs)
!caldo

!     NDIM = 2
!     NDOF = 4 ! will be reset within rtri

      call istkin(nva,4)

! ---------- allocate space

!     write (nout,fmt=4000)

4000  format (/,/,' memory allocation   ',/,' ',19 ('='),/)

      ifail = system("echo 'Running on' `uname -a` > triangle.log")
      ifail = system("date >> triangle.log")

      if(eulfs)then
         execmd = "rm -fv convergenza.dat"
         ifail = system(execmd)
         if(ifail.ne.0)call exit(ifail)

!        copy file .petsrc in home
!        for UNSTEADY EulFS simulations this file
!        will be overwritten with .petsrc_predictor and
!        .petsrc_corrector in their respective steps
!        They only differ in the dt value
         if (.not.UNSTEADY) then
         execmd = "cp -fv .petscrc .petscrc"
         endif
!        ifail = system(execmd)
         if(ifail.ne.0)call exit(ifail)
      endif

!     the initial grid is stored in a file called na00.1
      fname = "na00.1"
      fnameback = "na99"

! **********************************************************************
!  Read the original mesh, allocate space for the mesh variables and
!  allocate additional space for mesh points and segments
! **********************************************************************

      write(*,1001,advance='no')'readmesh               -->  '
      fndbnds=.true.
!     fndbnds=.false.
      call readmesh(
     +     lbndfac(0),
     +     lcelcel(0),
     +     lcelnod(0),
     +     lcorg(0),
     +     ledgptr(0),
     +     lnodcod(0),
     +     lnodptr(0),
     +     lzroe(0),
     +     nbfac(0),
     +     npoin(0),
     +     nelem(0),
     +     nhole(0),
     +     nbpoin(0),
     +     nvt,
     +     nedge(0),
     +     fname,
     +     lia(0),
     +     lja(0),
     +     liclr(0),
     +     nclr(0),
     +     fndbnds)
      write(*,1002)' ok'
1001  format(a)
1002  format(a)

! **********************************************************************
!  Read pmap
! **********************************************************************

      write(*,1001,advance='no')'readpmap               -->  '
      call readpmap(npoin(0),
     +     npnod(0),
     +     lpmap(0))
      write(*,1002)' ok'

! **********************************************************************
!  Make backups of some of the arrays of the background mesh:
!  the nodal flag (nodcode), the boundary data structure (bndfac),
!  the boundary node pointer (nodptr)
! **********************************************************************

      write(*,1001,advance='no')'copy mesh(0) in mesh(2)-->  '
      nitems = nbfac(0)+2*nshmax*neshmax
      lbndfac(2) = istkgt(3*nitems   ,2)
      lnodptr(2) = istkgt(3*nbpoin(0),2)
      lnodcod(2) = istkgt(npoin(0),2)
      call icopy(3*nbfac(0), istak(lbndfac(0)),1,istak(lbndfac(2)),1)
      call icopy(3*nbpoin(0),istak(lnodptr(0)),1,istak(lnodptr(2)),1)
      call icopy(npoin(0),istak(lnodcod(0)),1,istak(lnodcod(2)),1)
      nbfac(2)  = nbfac(0)
      nbpoin(2) = nbpoin(0)
      write(*,1002)' ok'

      lout(0) = istkst(1) ! number of arrays allocated in the stack
!     write(6,*)lout(0),' arrays have been allocated on the background mesh'

!     call x04eaf('general',' ',3,nbfac,istak(lbndfac(2)),3,
!    +            'bndry pointer(2) in main',ifail)
!     call x04eaf('general',' ',3,nbfac,istak(lbndfac(0)),3,
!    +            'bndry pointer(0) in main',ifail)
!
!     allocate an array to store the normal to the shock
!     only one normal is stored for each pair of shock points
!
!     lshnor = istkgt(ndim*nshmax*npshmax,4)
!
!     add the shock points and find the cells crossed by the shock and
!     the phantom points

! **********************************************************************
!  Read file input.dat containing information about mesh generation,
!  shock/discontinuity integration and additional hole point
! **********************************************************************

      write(*,1001,advance='no')'re_inp_data            -->  '
      call re_inp_data
      write(*,1002)' ok'

! **********************************************************************
!  Read the timesteps.dat file containing information about the dt
!  to be applied in the predictor/corrector steps (unsteady case)
!  Note: the timesteps.dat file should be already present in the
!        working directory
! **********************************************************************

      if (UNSTEADY) then

        write(*,1001,advance='no')'re_dt_data             -->  '
        open(unit=12,file='timesteps.dat',status='old',action='read')
        read(12,*) dtpr
        read(12,*) dtco
        close(12)
        write(*,1002)' ok'
        write(*,'(1x,f7.4, 1x, f7.4)') dtpr, dtco

      end if ! UNSTEADY

! **********************************************************************
!  Read file sh00.dat containing information about shock/discontinuity
! **********************************************************************

      write(*,1001,advance='no')'re_sdw_info            -->  '
      call re_sdw_info(
     +   xysh,
     +   dstak(lzroe(0)+npoin(0)*ndof),                        !upstream state
     +   dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof),    !downstream state
     +   zroeshuold,
     +   zroeshdold,
     +   istak(lnodcod(0)+npoin(0)),
     +   dstak(lcorg(0)),   !vale
     +   istak(lbndfac(0)), !vale
     +   nbfac(0),          !vale
     +   npoin(0),          !vale
     +   nshocks,
     +   nshockpoints,
     +   nshocksegs,
     +   typeshocks,
     +   nspecpoints,
     +   typespecpoints,
     +   shinspps,
     +   ispclr)
      write(*,1002)' ok'

! **********************************************************************
!  Shock points (equally-spaced) redistribution strategy
!  Note: points redistribution may not be mandatory but beneficial
! **********************************************************************

      if (STEADY) then

         write(*,1001,advance='no')'rd_sps_eq              -->  '
         call rd_dps_eq(
     +        xysh,
     +        dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +        dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +        nshocks,
     +        nshockpoints,
     +        nshocksegs)
         write(*,1002)' ok'

      call dcopy(nshmax*npshmax*ndof,dstak(lzroe(0)+npoin(0)*ndof),1,
     &           zroeshuold,1)
      call dcopy(nshmax*npshmax*ndof,
     &           dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof),1,
     &           zroeshdold,1)

      end if ! STEADY

!     call fx_sh_state(
!    +     dstak(lzroe(0)+npoin(0)*ndof), ! upstream state
!    +     nshocks,
!    +     nshockpoints,
!    +     nshocksegs)

!     call pr_sh_state(
!    +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
!    +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
!    +     nshocks,
!    +     nshockpoints,
!    +     nshocksegs)

! **********************************************************************
!  Compute the normal vectors
! **********************************************************************

      if (UNSTEADY) then

        write(*,1001,advance='no')'co_norm                -->  '
        call co_norm(
     +       xysh,
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       norsh,
     +       nshocks,
     +       nshockpoints,
     +       typeshocks,
     +       nspecpoints,
     +       typespecpoints,
     +       shinspps,
     +       ispclr,
     +       istak(lia(0)),
     +       istak(lja(0)),
     +       istak(liclr(0)),
     +       nclr(0),
     +       dstak(lcorg(0)))
        write(*,1002)' ok'

!       fix the normal orientation which otherwise
!       creates problems due to steady flow upstream
        if (testcase == 'ShockExpansion') then
          do no = 1, nshockpoints(1)
            norsh(1,no,1)=-norsh(1,no,1)
            norsh(2,no,1)=-norsh(2,no,1)
          enddo
        end if

! **********************************************************************
!  Compute the velocity of the shock (wsh)
! **********************************************************************

        write(*,1001,advance='no')'co_state_dps           -->  '
        call co_state_dps(
     +       xysh,
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       zroeshuold,
     +       zroeshdold,
     +       norsh,
     +       wsh,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks,
     +       i)
        write(*,1002)' ok'

! **********************************************************************
!  Fix the states in the discontinuity points
! **********************************************************************

        write(*,1001,advance='no')'fx_state_dps           -->  '
        call fx_state_dps(
     +       xysh,
     +       dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +       dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       zroeshuold,
     +       zroeshdold,
     +       norsh,
     +       wsh,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks,
     +       i,
     +       nspecpoints,
     +       typespecpoints,
     +       shinspps,
     +       ispclr,
     +       istak(lia(0)),
     +       istak(lja(0)),
     +       istak(liclr(0)),
     +       nclr(0),
     +       dstak(lcorg(0)))
        write(*,1002)' ok'

      end if ! UNSTEADY testcases

! **********************************************************************
!  Start the time loop
! **********************************************************************

      backdir = 'stepxyzk'

!     write(6,*)
!     write(6,*)' enter the  first  step '
!     write(6,*)
!     read(5,*)nbegin
!     write(6,*)
!     write(6,*)' enter the no steps '
!     write(6,*)
!     read(5,*)nsteps
      write(6,*)
      write(6,*)' starting the time loop'
      write(6,*)
      do 1000 i=1+nbegin,nsteps+nbegin
      write(6,*)'***********************************'
      write(6,*)' time level is ',i
      write(6,*)'                                   '

! **********************************************************************
!  Find the cells crossed by the shock and the phantom points and
!  update data structure on the boundary to take into account of the
!  presence of phantom points
! **********************************************************************

      write(*,1001,advance='no')'fnd_phps               -->  '
      call fnd_phps(
     +     nedge(0),
     +     istak(lbndfac(0)),
     +     nbfac(0),
     +     istak(lcelnod(0)),
     +     nvt,
     +     nelem(0),
     +     dstak(lcorg(0)),
     +     xysh,
     +     istak(lnodcod(0)),
     +     npoin(0),
     +     istak(lnodptr(0)),
     +     nbpoin(0),
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     nphampoints,
     +     istak(lpmap(0)))
      write(*,1002)' ok'

! **********************************************************************
!  Compute the normal unit vector to shocks and discontinuities
! **********************************************************************

      write(*,1001,advance='no')'co_norm                -->  '
      call co_norm(
     +     xysh,
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     norsh,
     +     nshocks,
     +     nshockpoints,
     +     typeshocks,
     +     nspecpoints,
     +     typespecpoints,
     +     shinspps,
     +     ispclr,
     +     istak(lia(0)),
     +     istak(lja(0)),
     +     istak(liclr(0)),
     +     nclr(0),
     +     dstak(lcorg(0)))
      write(*,1002)' ok'

!     fix the normal orientation which otherwise
!     creates problems due to steady flow upstream
      if (testcase == 'ShockExpansion') then
        do no = 1, nshockpoints(1)
          norsh(1,no,1)=-norsh(1,no,1)
          norsh(2,no,1)=-norsh(2,no,1)
        enddo
      end if

! **********************************************************************
!  Interp() updates values in the special points
! **********************************************************************

      if (STEADY) then

!       goto 2340
        write(*,1001,advance='no')'interp_sp              -->  '
        call interp_sp(
     +       istak(lcelnod(0)),
     +       nvt,
     +       nelem(0),
     +       dstak(lcorg(0)),
     +       dstak(lzroe(0)),
     +       xysh,
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       norsh,
     +       npoin(0),
     +       nshocks,
!    +       nshockpointsold)
     +       nshockpoints,
     +       typeshocks,
     +       nspecpoints,
     +       typespecpoints,
     +       shinspps)
        write(*,1002)' ok'

      end if ! STEADY

! **********************************************************************
!  Compute downstream and upstream mesh coordinate of each shock and
!  discontinuity point
! **********************************************************************

      write(*,1001,advance='no')'co_pnt_dspl            -->  '
      call co_pnt_dspl(
     +     xysh,
     +     dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +     dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     istak(lnodcod(0)+npoin(0)),
     +     norsh,
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     typeshocks,
     +     nspecpoints,
     +     typespecpoints,
     +     shinspps,
     +     ispclr)
      write(*,1002)' ok'

! **********************************************************************
!  Fix the mesh around the special points
! **********************************************************************

      write(*,1001,advance='no')'fx_msh_sps             -->  '
      call fx_msh_sps(
     +     istak(lbndfac(0)),
     +     istak(lnodcod(0)),
     +     nbfac(0),
     +     nbfac_sh,
     +     nvt,
     +     nelem(0),
     +     dstak(lcorg(0)),
     +     xysh,
     +     dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +     dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +     npoin(0),
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     nspecpoints,
     +     typespecpoints,
     +     shinspps,
     +     ispclr)
      write(*,1002)' ok'

! **********************************************************************
!  Copy shock(0) in shock(1)
!  Note: it calls mylib library dcopy(n,x,incx,y,incy)
!  The first vector is copied (in the second)
!  The unity increment defines the standard "direction" of the copy
!
!  The following vector fields are being copied:
!
!  xysh -------> xyshnew
!  zroeshdold -> zroeshdoldnew
!  zroeshuold -> zroeshuoldnew
!  norsh ------> norshnew
!  wsh --------> wshnew
!
!  This copy is required for the implementation of the
!  predictor/corrector algorithm: t+dt/2 --> old state
!                                 t+dt   --> new state
! **********************************************************************

      if (UNSTEADY) then

!       Copy here, so that if steady simulation is performed,
!       out of the corrector step, updated arrays are used,
!       and we can use "new" ones, at the cost of additional
!       memory allocation.

        write(*,1001,advance='no')'copy sh(0) in sh(1)    -->  '
        call dcopy(ndim*npshmax*nshmax,xysh,      1,xyshnew,      1)
        call dcopy(ndim*npshmax*nshmax,zroeshdold,1,zroeshdoldnew,1)
        call dcopy(ndim*npshmax*nshmax,zroeshuold,1,zroeshuoldnew,1)
        call dcopy(ndim*npshmax*nshmax,norsh,     1,norshnew,     1)
        call dcopy(ndim*npshmax*nshmax,wsh,       1,wshnew,       1)
        write(*,1002)' ok'

! **********************************************************************
!  Compute the grid velocity
!  Note: also the shock velocity vector is required
! **********************************************************************

        write(*,1001,advance='no')'calc_vel               -->  '
        call calc_vel(
     +       npoin(0),
     +       varray,
     +       dtpr,
     +       dstak(lcorg(0)),
     +       wsh,
     +       i,
     +       'y',
     +       nowtime,
     +       testcase)
        write(*,1002)' ok'

! **********************************************************************
!  It gives to EulFS info about grid velocity needed for the ALE
! **********************************************************************

        if (EULFS) then
        write(*,1001, advance='no')'solzne                -->   '
        call solzne(
     +       velfile,
     +       varray,
     +       ndim,
     +       npoin(0)+2*npshmax*nshmax,
     +       mode)
        write(*,1002)' ok'
        end if ! EULFS

      end if ! END UNSTEADY

! **********************************************************************
!  Write new poly file including all mesh points except phantom points
! **********************************************************************

      write(fname(3:7),fmt="(i5.5)")i
      write(*,1001,advance='no')'wtri                   -->  '
      call wtri(
     +     istak(lbndfac(0)),
     +     nbfac(0),
     +     nbfac_sh,
     +     istak(lcelnod(0)),
     +     nvt,
     +     dstak(lcorg(0)),
     +     xysh,
     +     dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +     dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +     dstak(lzroe(0)),
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     istak(lnodcod(0)),
     +     istak(lnodcod(0)+npoin(0)),
     +     npoin(0),
     +     fname(1:7),
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     nphampoints)
      write(*,1002)' ok'

! **********************************************************************
!  If use NEO and it is the 1st iteration, creates the neogrid0.grd
!  grid file in /NEO_data/input/
! **********************************************************************

      if (NEO) then ! NEO solver
        if ( i == 1+nbegin .and. testcase=="ShockExpansion" ) then
          write(*,1001,advance='no')'neogrid0               -->  '
          execmd = bindir(1:10) // 'neogrid0'
          ifail  = system(execmd)
          call flush(6)
          if (ifail.ne.0) then
            write(6,*)'neogrid0 has returned an error code ifail = ',
     +      ifail
            call exit(1)
          endif
          write(*,1002)' ok'
        end if
      end if

! **********************************************************************
!  Generate the new mesh
! **********************************************************************

!        write(6,*)
!        write(6,*)' meshing with triangle; input file is ',fname(1:7)
!        write(6,*)
         write(*,1001,advance='no')'triangle               -->  '

         execmd = bindir(1:10) // 'triangle_'//hostype(1:6)// ' -nep '
     +                  //fname(1:7)//' > log/triangle.log'
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
            write(6,*)'triangle has returned an error code ifail = ',
     +      ifail
            call exit(ifail)
         endif

         write(*,1002)' ok'

! **********************************************************************
!  Freeze mesh topology
! **********************************************************************

       if (STEADY) then

         if (imtf.ne.0.and.i.gt.imtf)then
           write(*,1001,advance='no')'mesh topology freezing -->  '

           fname2='stepXXXXX/naXXXXX.1'
           write(fname2(5:9),fmt="(i5.5)")imtf
           write(fname2(13:17),fmt="(i5.5)")imtf
           execmd ='cp '//fname2(1:19)//'.ele '//fname(1:7)//'.1.ele'
!          write(*,*)execmd
           ifail = system(execmd)
           call flush(6)
           if (ifail.ne.0) then
             write(6,*)'cp has returned an error code ifail = ',
     +       ifail
             call exit(1)
           endif
           execmd ='cp '//fname2(1:19)//'.neigh '//
     +                    fname(1:7)//'.1.neigh'
!          write(*,*)execmd
           ifail = system(execmd)
           call flush(6)
           if (ifail.ne.0) then
             write(6,*)'cp has returned an error code ifail = ',
     +       ifail
             call exit(1)
           endif
           execmd ='cp '//fname2(1:19)//'.edge '//fname(1:7)//'.1.edge'
!          write(*,*)execmd
           ifail = system(execmd)
           call flush(6)
           if (ifail.ne.0) then
             write(6,*)'cp has returned an error code ifail = ',
     +       ifail
             call exit(1)
           endif

           write(*,1002)' ok'

         endif

       end if ! STEADY

! ***********************************
      if (EULFS) then ! EulFS SOLVER
! ***********************************

! **********************************************************************
!  Convert the triangle files into a fmt readable by the code using:
!  echo na0x.1 | triangle2dat
! **********************************************************************

         write(*,1001,advance='no')'triangle2dat           -->  '
!        execmd = "echo " // fname(1:7)
!     +   // ".1 |" // bindir(1:10) // "triangle2dat_" // hostype(1:6)
!     +   // " > log/triangle2dat.log"
         if (nprdbnd .eq. 0) then                            ! for the cases without periodic BCs
           execmd = "printf '" // fname(1:7)
     +     // ".1\nn'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"
         elseif(nprdbnd.eq.1.and.prdbndclr(3,1).eq.1)then    ! for the cases with only one periodic boundary
           write(color1,fmt="(i2.2)")prdbndclr(1,1)          ! with points having the same x
           write(color2,fmt="(i2.2)")prdbndclr(2,1)
           execmd = "printf '" // fname(1:7)
     +     // ".1\ny\n"
     +     // color1 // "\n"
     +     // color2 // "\nx'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"
         elseif (nprdbnd.eq.1 .and. prdbndclr(3,1).eq.2) then ! for the cases with only one periodic boundary
           write(color1,fmt="(i2.2)")prdbndclr(1,1)           ! with points having the same y
           write(color2,fmt="(i2.2)")prdbndclr(2,1)
           execmd = "printf '" // fname(1:7)
     +     // ".1\ny\n"
     +     // color1 // "\n"
     +     // color2 // "\ny'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"

         else ! for cases with more thatn one periodic boundary
           write(*,*)' case not implemented!'
         endif

!        write(*,*)execmd

!        execmd = "printf '" // fname(1:7)
!    +   // ".1\ny\n2\n4\nx'|" ! for case cascade
!    +   // ".1\ny\n1\n3\nx'|" ! for case nacapar2
!    +   // ".1\nn'|"          ! for cases w/o periodic BCs'
!    +   // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
!    +   // " > log/triangle2dat.log"
!         write(*,*)execmd

         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)execmd
           write(6,*)'triangle2dat has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

         write(*,1002)' ok'

! ****************************
!  Run one step of EulFS code
! ****************************

         if (UNSTEADY) then
!          It runs the predictor step of the EulFS code (we need to use dt/2)
           execmd = "cp -f .petscrc_predictor .petscrc"
           ifail = system(execmd)
           if(ifail.ne.0)call exit(1)
         end if

         write(*,1001,advance='no')'eulfs                  -->  '
!        execmd = bindir(1:10) // "eulfs11.13_"
!     +                        //gastype(1:4)//"_"//hostype(1:6)
!     +  // " -itmax 1 > log/eulfs.log"
         execmd = bindir(1:10) // "EulFS_"//hostype(1:6)
     +   // " -itmax 1 > log/eulfs.log"

         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
            write(6,*)'eulfs has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

         if (unsteady) then
           execmd = "cp step000001.dat file001.dat"
           ifail = system(execmd)
           if(ifail.ne.0)call exit(1)
         endif

!        execmd = "cp file003.dat file010.dat"
!        ifail = system(execmd)
!        if(ifail.ne.0)call exit(1)

         write(*,1002)' ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0x.1 | dat2triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.bak
! **********************************************************************

         write(*,1001,advance='no')'dat2triangle           -->  '
!        execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
!    +   // "dat2triangle_" // hostype(1:6)
!    +   // ">log/dat2triangle.log"
         execmd = "printf '"//fname(1:7)//".1' | "// bindir(1:10)
     +   // "dat2triangle-NEW-" // hostype(1:6)
     +   // ">log/dat2triangle.log"
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)'dat2triangle has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

         write(*,1002)' ok'

! ***********************************
      elseif (NEO) then ! NEO SOLVER
! ***********************************

!      neogrid0 works only for the 1st iteration of ShockVortex
!      but in all other cases we need this conversion
       if( STEADY .or. i /= 1+nbegin .or. testcase=="ShockVortex") then

         write(*,1001,advance='no')'na00xTovvvv            -->   '
         execmd = "echo " // fname(1:7)
     +          // ".1 |" // bindir(1:10) // "na2vvvv"
     +          // " > log/na2vvvv.log"
         ifail = system(execmd)
         write(*,1002)'ok'

       end if

! ****************************************************************
!       convert the triangle files into a format readable by the
!       code using: echo na0X.1 | triangle2dat
! ****************************************************************

         write(*,1001,advance='no')'triangle2grd           -->   '
         execmd = "echo " // fname(1:7)
     +   // ".1 |" // bindir(1:10) // "triangle2grd"
     +   // " > log/triangle2grd.log"
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)'triangle2grd has returned an error code ifail = ',
     &ifail
           call exit(ifail)
         endif
         write(*,1002)'ok'

!     ******************
      if (UNSTEADY) then
!     ******************

!       if it is the 1st iteration
!       **************************
        if (i == 1+nbegin) then

          write(*,1001,advance='no')'NEO 1st iteration      -->  '

          execmd = bindir(1:10) // "CRD_euler"
     +                          // "> log/neo.log"
          ifail = system(execmd)
          call flush(6)
          if (ifail.ne.0) then
            write(6,*)'NEO (1st iteration)
     +      has returned an error code ifail = ',ifail
          call exit(1)
          endif

          execmd = "cp ./NEO_data/output/vvvv.dat "//
     +                "./NEO_data/output/vvvv0.dat "
          ifail  = system(execmd)

          execmd = "mv ./NEO_data/output/vvvv.dat "//
     +               "./NEO_data/output/vvvv_input.dat "
          ifail  = system(execmd)

!         Here the following happens (for unsteady cases):
!         - the 1st iteration uses NEO_data/textinput/inputfile-exp.txt
!         - After the 1st NEO call, inputfile-exp.txt is moved in BAK
!         - Then, the inputfile-exp.txt in the testcase folder is
!           copied in NEO_data/textinput/inputfile-exp.txt
!         This is done because the two inputfile-exp.txt files differ
!         for the "Initial state" value. In the first case, it is 14
!         which means that the NEO function initial_solution()
!         writes the initial solution for the centered expansion test,
!         while after the 1st iteration should be 0, since we don't need
!         to initialize it again but instead just read_solution() which
!         happens if the variable intial_solution = 0.

          execmd = "mv ./NEO_data/textinput/inputfile-exp.txt "//
     +                "./NEO_data/textinput/inputfile-exp.txt.BAK "
          ifail  = system(execmd)

          execmd = "cp inputfile-exp.txt "//"./NEO_data/textinput/"
          ifail  = system(execmd)

          write(*,1002)' ok'

        end if ! 1ST ITERATION

      end if ! UNSTEADY


!     for all the other iterations
!     ****************************
         write(*,1001,advance='no')'NEO                    -->   '
         execmd = bindir(1:10) // "CRD_euler"
     +                         // "> log/neo.log"
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
            write(6,*)'neo has returned an error code ifail = ',
     &ifail
            call exit(ifail)
         endif
         write(*,1002)'ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0X.1 | NEO2triangle ex (dat2triangle)
!  Note: the file na0X.1.node will be overwritten with the values
!        updated by the code and a copy with "old" values is copied in
!        na0X.1.node.BAK
! **********************************************************************

         write(*,1001,advance='no')'NEO2triangle           -->   '
         execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
     +   // "NEO2triangle" // ">log/NEO2triangle.log"
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)'neo2triangle has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif
         write(*,1002)'ok'

      ELSE

         write(*,*)'should be running either EULFS ',eulfs,' or NEO ',
     &              neo
         call exit(10)

      ENDIF ! IF-THEN-ELSE ON THE CFD CODE

! **********************************************************************
!  Here the corrector step starts
! **********************************************************************

!     ******************
      if (UNSTEADY) then
!     ******************

! **********************************************************************
!  Read the mesh generated by triangle, with nodal values updated by
!  the code, allocate space for the mesh variables and allocate
!  additional space for mesh points and segments.
!  Phantom nodes will be read as well, with wrong nodal coordinates
!  and values, but we do not care since phantom nodes are not addressed
!  in the connectivity
! **********************************************************************

        write(*,1001,advance='no')'readmesh               -->  '
        fname(1:9) = fname(1:7)//".1"
        fndbnds = .false.
        call readmesh(
     +       lbndfac(1),
     +       lcelcel(1),
     +       lcelnod(1),
     +       lcorg(1),
     +       ledgptr(1),
     +       lnodcod(1),
     +       lnodptr(1),
     +       lzroe(1),
     +       nbfac(1),
     +       npoin(1),
     +       nelem(1),
     +       nhole(1),
     +       nbpoin(1),
     +       nvt,
     +       nedge(1),
     +       fname,
     +       lia(1),
     +       lja(1),
     +       liclr(1),
     +       nclr(1),
     +       fndbnds)

        lout(1) = istkst(1) ! number of arrays allocated in the stack
        write(*,1002)' ok'

        ! TODO: check whether FX_USTATE should be added here ...

! **********************************************************************
!  Update the nodal values on the backgroud grid (0) using values of
!  the shocked grid (1); the shocked grid contains "wrong" values in
!  the phantom nodes but these will be changed at a later stage in
!  interp() we need to perform this copy here, since the shockmov()
!  routine works on nodal values of grid (0)
! **********************************************************************

        totshockpoints=2*nshmax*npshmax

        write(*,1001,advance='no')'zroe(1)->zroe(0)       -->  '

        if ( npoin(1) .eq. (npoin(0)+totshockpoints)) then

          call dcopy(ndof*npoin(1),dstak(lzroe(1)),1,dstak(lzroe(0)),1)
          write(*,1002)' ok'

        else

!         the nof gridpoints in grid(1) must equal the number of
!         gridpoints on the background mesh + 2 * nshockpoints

          write(6,*) 'there is a mismatch in the nof gridpoints'
          write(6,*) 'btw grid(0) and grid(1)'
          write(*,*) npoin(0), totshockpoints
          write(*,*) npoin(1), totshockpoints
          call exit(1)

        endif

!       work is a work array used to store nodal values in the shock points
!       work is used in interp() and shockmov()

        lwork    = istkgt(2*ndof*nshmax*npshmax,4) ! work array
        lxyshold = istkgt(2*ndim*nshmax*npshmax,4) ! work array
        lxyshnew = istkgt(2*ndim*nshmax*npshmax,4) ! work array

! **********************************************************************
!  Updates nodal values in all the shock points of grid (0) using R-H
!  relations and compute the shock speed
!  Note: xysh coordinates are used only to write tecplot file but not
!        elsewhere
! **********************************************************************

        write(*,1001,advance='no')'co_state_dps           -->  '
        call co_state_dps(
     +       xyshnew,
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       zroeshuoldnew,
     +       zroeshdoldnew,
     +       norshnew,
     +       wshnew,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks,
     +       i)
        write(*,1002)' ok'

! **********************************************************************

        write(*,1001,advance='no')'fx_state_dps           -->  '
        call fx_state_dps(
     +       xyshnew,                                           ! not used
     +       dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream   coord.
     +       dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream   state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       zroeshuoldnew,
     +       zroeshdoldnew,
     +       norshnew,
     +       wshnew,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks,
     +       i,
     +       nspecpoints,
     +       typespecpoints,
     +       shinspps,
     +       ispclr,
     +       istak(lia(0)),
     +       istak(lja(0)),
     +       istak(liclr(0)),
     +       nclr(0),
     +       dstak(lcorg(0)))
        write(*,1002)' ok'

! **********************************************************************

        write(*,1001,advance='no')'zroesh(0)->zroesh(1)   -->  '
        call dcopy(ndof*totshockpoints,
     +       dstak(lzroe(0)+npoin(0)*ndof),1,
     +       dstak(lzroe(1)+npoin(0)*ndof),1)
        write(*,1002)' ok'

! **********************************************************************

        write(*,1001,advance='no')'calc_vel               -->  '
        call calc_vel(
     +       npoin(0),
     +       varray,
     +       dtco,
     +       dstak(lcorg(0)),
     +       wsh, !WSHnew?
     +       i,
     +       'n',
     +       nowtime,
     +       testcase)
        write(*,1002)' ok'

! **********************************************************************
!  It gives to eulfs information about grid velocity (corrector step)
! **********************************************************************

        if (EULFS) then
        write(*,1001,advance='no')'solzne                 -->   '
        call solzne(
     +       velfile,
     +       varray,
     +       ndim,
     +       npoin(0)+2*npshmax*nshmax,
     +       mode)
        write(*,1002)' ok'
        end if

! **********************************************************************
!  It generates the new mesh
! **********************************************************************

        write(*,1001,advance='no')'triangle               -->  '
        execmd = bindir(1:10) // 'triangle_'//hostype(1:6)// ' -nep '
     +                        //fname(1:7)//' > log/triangle.log'
        ifail  = system(execmd)
        call flush(6)
        if (ifail.ne.0) then
          write(6,*)'triangle has returned an error code ifail = ',
     +    ifail
          call exit(1)
        endif
        write(*,1002)' ok'

! ***********************************
      if (EULFS) then ! EulFS SOLVER
! ***********************************

! **********************************************************************
!  Convert the triangle files into a fmt readable by the code using:
!  echo na0x.1 | triangle2dat
! **********************************************************************

         write(*,1001,advance='no')'triangle2dat           -->  '
!        execmd = "echo " // fname(1:7)
!     +   // ".1 |" // bindir(1:10) // "triangle2dat_" // hostype(1:6)
!     +   // " > log/triangle2dat.log"
         if (nprdbnd .eq. 0) then                            ! for the cases without periodic BCs
           execmd = "printf '" // fname(1:7)
     +     // ".1\nn'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"
         elseif(nprdbnd.eq.1.and.prdbndclr(3,1).eq.1)then    ! for the cases with only one periodic boundary
           write(color1,fmt="(i2.2)")prdbndclr(1,1)          ! with points having the same x
           write(color2,fmt="(i2.2)")prdbndclr(2,1)
           execmd = "printf '" // fname(1:7)
     +     // ".1\ny\n"
     +     // color1 // "\n"
     +     // color2 // "\nx'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"
         elseif (nprdbnd.eq.1 .and. prdbndclr(3,1).eq.2) then ! for the cases with only one periodic boundary
           write(color1,fmt="(i2.2)")prdbndclr(1,1)           ! with points having the same y
           write(color2,fmt="(i2.2)")prdbndclr(2,1)
           execmd = "printf '" // fname(1:7)
     +     // ".1\ny\n"
     +     // color1 // "\n"
     +     // color2 // "\ny'|"
     +     // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
     +     // " > log/triangle2dat.log"

         else ! for cases with more thatn one periodic boundary
           write(*,*)' case not implemented!'
         endif

!        write(*,*)execmd

!        execmd = "printf '" // fname(1:7)
!    +   // ".1\ny\n2\n4\nx'|" ! for case cascade
!    +   // ".1\ny\n1\n3\nx'|" ! for case nacapar2
!    +   // ".1\nn'|"          ! for cases w/o periodic BCs'
!    +   // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
!    +   // " > log/triangle2dat.log"
!         write(*,*)execmd

         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)execmd
           write(6,*)'triangle2dat has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

         write(*,1002)' ok'

! ****************************
!  Run one step of EulFS code
! ****************************

!        if (UNSTEADY) then
!          It runs the corrector step of the EulFS code (now we use the full dt)
           execmd = "cp -f .petscrc_corrector .petscrc"
           ifail = system(execmd)
           if(ifail.ne.0)call exit(1)
!        end if

         write(*,1001,advance='no')'eulfs                  -->  '
!        execmd = bindir(1:10) // "eulfs11.13_"
!     +                        //gastype(1:4)//"_"//hostype(1:6)
!     +  // " -itmax 1 > log/eulfs.log"
         execmd = bindir(1:10) // "EulFS_"//hostype(1:6)
     +   // " -itmax 1 > log/eulfs.log"

         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
            write(6,*)'eulfs has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

!        if (UNSTEADY) then
           execmd = "cp step000001.dat file001.dat"
           ifail = system(execmd)
           if(ifail.ne.0)call exit(1)
!        endif

!        execmd = "cp file003.dat file010.dat"
!        ifail = system(execmd)
!        if(ifail.ne.0)call exit(1)

         write(*,1002)' ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0x.1 | dat2triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.bak
! **********************************************************************

         write(*,1001,advance='no')'dat2triangle           -->  '
!        execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
!    +   // "dat2triangle_" // hostype(1:6)
!    +   // ">log/dat2triangle.log"
         execmd = "printf '"//fname(1:7)//".1' | "// bindir(1:10)
     +   // "dat2triangle-NEW-" // hostype(1:6)
     +   // ">log/dat2triangle.log"
         ifail = system(execmd)
         call flush(6)
         if(ifail.ne.0)then
           write(6,*)'dat2triangle has returned an error code ifail = ',
     &ifail
            call exit(1)
         endif

         write(*,1002)' ok'

      elseif (NEO) then ! NEO SOLVER for UNSTEADY predictor step

! **********************************************************************
!   Update of vvvv.dat (input) for NEO
! **********************************************************************

        write(*,1001,advance='no')'na2vvvv                -->  '
        execmd = "echo " // fname(1:7)
     +    // ".1 |" // bindir(1:10) // "na2vvvv"
     +    // " > log/na2vvvv.log"
        ifail  = system(execmd)
        write(*,1002)' ok'

! **********************************************************************
!  convert the triangle files into a fmt readable by the code
!  using: echo na0x.1 | triangle2dat
! **********************************************************************

        write(*,1001,advance='no')'triangle2grd           -->  '
        execmd = "echo " // fname(1:7)
     +    // ".1 |" // bindir(1:10) // "triangle2grd"
     +    // " > log/triangle2grd.log"
        ifail = system(execmd)
        call flush(6)
        if (ifail.ne.0) then
          write(6,*)'triangle2grd has returned an error code ifail = ',
     &    ifail
          call exit(1)
        endif
        write(*,1002)' ok'

! **********************************************************************
!  One step with NEO
! **********************************************************************

        write(*,1001,advance='no')'NEO                    -->  '
        execmd = bindir(1:10) // "CRD_euler"
     +                        // " > log/neo.log"
        ifail  = system(execmd)
        call flush(6)
        if (ifail.ne.0) then
          write(6,*)'NEO has returned an error code ifail = ',
     +    ifail
          call exit(1)
        endif
        write(*,1002)' ok'

! **********************************************************************
!  Convert dat file to triangle file
! **********************************************************************

        write(*,1001,advance='no')'NEO2triangle           -->  '
        execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
     +                  // "NEO2triangle" // ">log/NEO2triangle.log"
        ifail = system(execmd)
        call flush(6)
        if (ifail.ne.0) then
          write(6,*)'NEO2triangle has returned an error code ifail = ',
     +    ifail
          call exit(1)
        endif
        write(*,1002)' ok'

      end if ! end SOLVER (EULFS/NEO) for UNSTEADY (corrector step)

! **********************************************************************

        write(*,1001,advance='no')'mv_grid                -->  '
        call mv_grid(
     +       npoin(0),
     +       varray,
     +       dtco,
     +       dstak(lcorg(0)),
     +       wshnew,
     +       i,
     +       testcase) ! as in calc_vel, added arg to switch case
        write(*,1002)' ok'

! **********************************************************************

!       write(*,1001,advance='no')'copy sh(0) in sh(1)    -->  '
!       call dcopy(ndim*npshmax*nshmax,xysh,1,xyshnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshdold,1,zroeshdoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshuold,1,zroeshuoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,norsh,1,norshnew,1)
        call dcopy(ndim*npshmax*nshmax,wshnew,1,wsh,1)
!       write(*,1002)' ok'

      end if ! UNSTEADY TODO: check the variables in this part ...
!     *****************

! **********************************************************************
!  Read the mesh generated by triangle, with nodal values updated by
!  the code, allocate space for the mesh variables and allocate
!  additional space for mesh points and segments.
!  Phantom nodes will be read as well, with wrong nodal coordinates
!  and values, but we do not care since phantom nodes are not addressed
!  in the connectivity
! **********************************************************************

      write(*,1001,advance='no')'readmesh               -->  '
      fname(1:9) = fname(1:7)//".1"
      fndbnds=.false.
      call readmesh(
     +     lbndfac(1),
     +     lcelcel(1),
     +     lcelnod(1),
     +     lcorg(1),
     +     ledgptr(1),
     +     lnodcod(1),
     +     lnodptr(1),
     +     lzroe(1),
     +     nbfac(1),
     +     npoin(1),
     +     nelem(1),
     +     nhole(1),
     +     nbpoin(1),
     +     nvt,
     +     nedge(1),
     +     fname,
     +     lia(1),
     +     lja(1),
     +     liclr(1),
     +     nclr(1),
     +     fndbnds)

      lout(1) = istkst(1) ! number of arrays allocated in the stack
      write(*,1002)' ok'

! **********************************************************************
!  TODO: comment this procedure
! **********************************************************************

!     call pr_sh_state(
!    +   dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
!!   +   dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
!    +   nshocks,
!    +   nshockpoints,
!    +   nshocksegs)

! **********************************************************************
!  Update the nodal values on the backgroud grid (0) using values of
!  the shocked grid (1); the shocked grid contains "wrong" values in
!  the phantom nodes but these will be changed at a later stage in
!  interp() we need to perform this copy here, since the shockmov()
!  routine works on nodal values of grid (0)
! **********************************************************************

      totshockpoints=2*nshmax*npshmax

      write(*,1001,advance='no')'zroe(1)->zroe(0)       -->  '
      if (npoin(1) .eq. (npoin(0)+totshockpoints)) then
        call dcopy(ndof*npoin(1),dstak(lzroe(1)),1,
     +                           dstak(lzroe(0)),1)

      write(*,1002)' ok'

      else

!       the nof gridpoints in grid(1) must equal the number of
!       gridpoints on the background mesh + 2 x nshockpoints

        write(6,*)'there is a mismatch in the nof gridpoints'
        write(6,*)'btw grid(0) and grid(1)'
        write(*,*)npoin(0),totshockpoints
        write(*,*)npoin(1),totshockpoints
        call exit(1)
      endif

!     work is a work array used to store nodal values in the shock points
!     work is used in interp() and shockmov()

      lwork    = istkgt(2*ndof*nshmax*npshmax,4)
      lxyshold = istkgt(2*ndim*nshmax*npshmax,4) ! work array
      lxyshnew = istkgt(2*ndim*nshmax*npshmax,4) ! work array

!     SHOCKmov updates nodal values in the shock points of grid (0)
!     computes R-H relations, moves the shock
!     Nodal values in the shockpoints of the shocked mesh (1) are
!     updated too, so that the correct downstream values are available
!     for interpolation in interp()
!     REM: nShockPoints might be changed by SHOCKmov()

!        do ish=1,nshocks
!          nshockpointsold(ish)=nshockpoints(ish)
!        enddo

!     REcreate shock edges near the triple point

!     write(6,*)' calling chktpnt in main '
!     call chktpnt2(istak(lbndfac(0)),nbfac(0),nbfac_sh,
!    &     istak(lcelnod(0)),nvt,
!    &     nelem(0),
!    &     dstak(lcorg(0)),
!    &     dstak(lcorg(0)+npoin(0)*ndim),
!    &     dstak(lcorg(1)+npoin(0)*ndim),
!    &     dstak(lzroe(0)+npoin(0)*ndof),ndof,
!    &     ndim,
!    &     npoin(0),
!    &     nshocks,nshockpoints,nshocksegs)

! **********************************************************************
!  Update nodal values in all the shock points of grid (0) using R-H
!  relations and compute the shock speed
! **********************************************************************

      write(*,1001,advance='no')'co_state_dps           -->  '
      call co_state_dps(
     +     xysh,
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     zroeshuold, !ZROESHuOLDnew?
     +     zroeshdold, !ZROESHdOLDnew?
     +     norsh,      !NORSHnew?
     +     wsh,        !WSHnew?
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     typeshocks,
     +     i)
      write(*,1002)' ok'

! **********************************************************************
!  Fix and correct the nodal values and shock speed in all special
!  point using the correct s-s interaction relation
! **********************************************************************

      write(*,1001,advance='no')'fx_state_dps           -->  '
      call fx_state_dps(
     +     xysh, !XYSHnew?
     +     dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +     dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     zroeshuold, !ZROESHuOLDnew?
     +     zroeshdold, !ZROESHdOLDnew?
     +     norsh,      !NORSHnew?
     +     wsh,        !WSHnew?
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     typeshocks,
     +     i,
     +     nspecpoints,
     +     typespecpoints,
     +     shinspps,
     +     ispclr,
     +     istak(lia(0)),
     +     istak(lja(0)),
     +     istak(liclr(0)),
     +     nclr(0),
     +     dstak(lcorg(0)))
      write(*,1002)' ok'

! **********************************************************************
!  Update the nodal values of shocks on the grid (1)
!  Note: we need to perform this copy here, since the fx_state_sps and
!        co_state_dps routines work on nodal values of grid (0)
! **********************************************************************

      write(*,1001,advance='no')'zroesh(0)->zroesh(1)   -->  '
      call dcopy(ndof*totshockpoints,
     +     dstak(lzroe(0)+npoin(0)*ndof),1,
     +     dstak(lzroe(1)+npoin(0)*ndof),1)
      write(*,1002)' ok'

! **********************************************************************
!  Calculate the mean shock velocity
! **********************************************************************

      if (UNSTEADY) then
!     ******************

        write(*,1001,advance='no')'wsh_mean               -->  '
        call wsh_mean(
     +       wsh,
     +       wshnew,
     +       wshmean)
        write(*,1002)' ok'

      end if ! UNSTEADY

! **********************************************************************
!  Move the shocks
! **********************************************************************

      write(*,1001,advance='no')'mv_dps                 -->  '
      call mv_dps(
     +     xysh,
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     wsh,
     +     i,
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     typeshocks)
      write(*,1002)' ok'

! **********************************************************************
!  Fix and correct the nodal position in all special point
! **********************************************************************

      if (STEADY) then
!     ****************

        write(*,1001,advance='no')'fx_dps_loc             -->  '
        call fx_dps_loc(
     +       xysh,
     +       dstak(lcorg(0)+npoin(0)*ndim),                     ! upstream coord.
     +       dstak(lcorg(0)+npoin(0)*ndim+nshmax*npshmax*ndim), ! downstream coord.
     +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       zroeshuold,
     +       zroeshdold,
     +       norsh,
     +       wsh,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks,
     +       i,
     +       nspecpoints,
     +       typespecpoints,
     +       shinspps,
     +       ispclr,
     +       istak(lia(0)),
     +       istak(lja(0)),
     +       istak(liclr(0)),
     +       nclr(0),
     +       dstak(lzroe(0)), ! vale
     +       dstak(lcorg(0)),
     +       shtopolchanged)  ! vale
       write(*,1002)' ok'

! **********************************************************************
!  Filters the shocks
! **********************************************************************

        write(*,1001,advance='no')'fltr_dls               -->  '
        call fltr_dls(
     +       xysh,
     +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +       wsh,
     +       i,
     +       nshocks,
     +       nshockpoints,
     +       nshocksegs,
     +       typeshocks)
        write(*,1002)' ok'

      end if ! STEADY

! **********************************************************************
!  Interp() updates values in the phantom nodes of the background
!  mesh (0) using values in the shocked mesh (1)
!  Note: npoin is passed as array
! **********************************************************************

!      goto 2340
       write(*,1001,advance='no')'interp                 -->  '
       call interp(
     +      istak(lbndfac(1)),
     +      nbfac(1),
     +      istak(lcelnod(1)),
     +      nvt,
     +      nelem(1),
     +      dstak(lcorg(1)),
     +      dstak(lzroe(1)),
     +      xysh,
     +      dstak(lcorg(1)+npoin(0)*ndim),                     !upstream coord.
     +      dstak(lcorg(1)+npoin(0)*ndim+nshmax*npshmax*ndim), !downstream coord.
     +      nphampoints,
     +      dstak(lcorg(0)),
     +      dstak(lzroe(0)),
     +      istak(lnodcod(0)),
     +      npoin(0),
     +      nshocks,
!    +      nshockpointsold)
     +      nshockpoints,
     +      istak(lia(1)),
     +      istak(lja(1)),
     +      istak(liclr(0)),
     +      nclr(0))
       write(*,1002)' ok'

2340   continue

! **********************************************************************
!  Redistribute the shock nodes
!  Note: this redistribution is different from the previous one
! **********************************************************************

!     if(i.gt.1000)  goto 3450
      write(*,1001,advance='no')'rd_dps                 -->  '
      call rd_dps(
     +     xysh,
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs)
      write(*,1002)' ok'

3450  continue

! **********************************************************************
!  Redistribute the shock nodes
!  Note: shock node redistribution may be not necessary but beneficial
! **********************************************************************

!     if (mod(i,17).ne.0) goto 3451
!       write(*,1001)'rd_sps_eq              -->   '
!       call rd_dps_eq(
!    +       xysh,
!    +       dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
!    +       dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
!    +       nshocks,
!    +       nshockpoints,
!    +       nshocksegs)

!       write(*,1002)'rd_sps_eq              --> ok'

3451  continue

! **********************************************************************
!  Write the triangle node file of the background with the updated
!  value of zroe
! **********************************************************************

      write(*,1001,advance='no')'wtri0                  -->  '
      call wtri0(
     +     dstak(lcorg(0)),
!    +     ndim,
     +     dstak(lzroe(0)),
!    +     ndof,
     +     istak(lnodcod(0)),
     +     npoin(0),
     +     fnameback(1:4))
      write(*,1002)' ok'

! **********************************************************************
!  Write file sh99.dat containing information about shock/discontinuity
! **********************************************************************

      write(*,1001,advance='no')'wrt_sdw_info           -->  '
      call wrt_sdw_info(
     +     xysh,
     +     dstak(lzroe(0)+npoin(0)*ndof),                     ! upstream state
     +     dstak(lzroe(0)+npoin(0)*ndof+nshmax*npshmax*ndof), ! downstream state
     +     istak(lnodcod(0)+npoin(0)),
     +     nshocks,
     +     nshockpoints,
     +     nshocksegs,
     +     typeshocks,
     +     nspecpoints,
     +     typespecpoints,
     +     shinspps,
     +     ispclr)
      write(*,1002)' ok'

! **********************************************************************
!  Copy new shock on the old one
! **********************************************************************

      if (UNSTEADY) then

        write(*,1001,advance='no')'copy sh(1) in sh(0)    -->  '
        call dcopy(ndim*npshmax*nshmax,zroeshdoldnew,1, zroeshdold,1)
        call dcopy(ndim*npshmax*nshmax,zroeshuoldnew,1, zroeshuold,1)
        call dcopy(ndim*npshmax*nshmax,norshnew,1, norsh,1)
        call dcopy(ndim*npshmax*nshmax,wshnew,1, wsh,1)
        write(*,1002)' ok'

      end if ! UNSTEADY

      call istkrl(3)

!     call solzne("file004.dat",dstak(lzroe(0)),ndof,npoin(0),"w")

!     call x04eaf('general',' ',3,nbfac,istak(lbndfac(2)),3,
!    +            'bndry pointer(2) in main',ifail)
!     call x04eaf('general',' ',3,nbfac,istak(lbndfac(0)),3,
!    +            'bndry pointer(0) in main',ifail)
!     pause

! **********************************************************************
!  Restore the original arrays of the background grid
! **********************************************************************

      nbfac(0)  = nbfac(2)
      nbpoin(0) = nbpoin(2)
      nitems    = nbfac(2)
      call icopy(3*nitems,   istak(lbndfac(2)),1,istak(lbndfac(0)),1)
      call icopy(3*nbpoin(2),istak(lnodptr(2)),1,istak(lnodptr(0)),1)
      call icopy(npoin(0),istak(lnodcod(2)),1,istak(lnodcod(0)),1)

! **********************************************************************
!  Release all pointers allocated for the shocked mesh (1)
! **********************************************************************

      call istkrl(lout(1)-lout(0))

! **********************************************************************
!  Create a directory to backup files
! **********************************************************************

      write(backdir(5:9),fmt="(i5.5)")i
      if(mod(i-1,ibak).eq.0)then
          execmd = "mkdir -v "//backdir(1:9)
          ifail = system(execmd)
!     execmd = "mv shocknor.dat shock.log file00[1-3].dat file010.dat fs
!    &pl.out shocks.dat "//fname(1:7)//".* "//backdir(1:9)
          if(eulfs)then
             execmd = "mv -v shocknor.dat file00[1-4].dat file010.dat sh
     &99.dat  "//fname(1:7)//".* "//fnameback(1:4)//".node
     &                 "//backdir(1:9)
             ifail = system(execmd)
          elseif(neo)then
             execmd = "mv -v shocknor.dat sh99.dat
     &                  "//fname(1:7)//".* "//fnameback(1:4)//".node "//
     &                backdir(1:9)
             ifail = system(execmd)
             execmd =
     &       "cp -vp ./NEO_data/input/neogrid.grd
     &       ./NEO_data/input/vel.dat "//backdir(1:9)

!    &       "mv -v ../../source_utils/NEO_source/input/neogrid.grd
!    &       ../../source_utils/NEO_source/input/vel.dat "//backdir(1:9)

! Note: make attention with cp/mv because of vel.dat file

          ifail = system(execmd)
!
!         if (i == 1+nbegin) then
!            execmd = "mv ../../../NEO_source/output/vvvv0.dat "
!    &             //backdir(1:9)
!            ifail = system(execmd)
!         end if
!
             execmd =
     &       "cp -v ./NEO_data/output/vvvv.dat "// backdir(1:9)
             ifail = system(execmd)

             execmd =
     &       "mv -v ./NEO_data/output/vvvv_input.dat "//backdir(1:9)
             ifail = system(execmd)

! copy the solution fie in Tec folder
!         if (i == 1+nbegin) then
!            execmd = "cp "//backdir(1:9)//"/vvvv0.dat "
!    &            //"tec/vvvv00000.dat"
!            ifail = system(execmd)
!         end if

!         execmd = "cp ../../source_utils/NEO_source/output/vvvv.dat "//"Tec/vvvv"
!    &         //backdir(5:9)//".dat"
!         ifail = system(execmd)

          endif
!         ifail = system(execmd)
      else ! backing up or not ...
!     execmd = "rm shocknor.dat shock.log file00[1-3].dat file010.dat fs
!    &pl.out shocks.dat "//fname(1:7)//".*"
          if(EULFS)then
      execmd = "rm shocknor.dat file00[1-4].dat file010.dat
     &                  "//fname(1:7)//".* "//fnameback(1:4)//".node "//
     &  "sh99.dat "
          elseif(NEO)then
      execmd = "rm shocknor.dat "//fname(1:7)//".* "//fnameback(1:4)//
     &".node "//"sh99.dat "
          endif
        ifail = system(execmd)
      endif

      if (EULFS) then
        execmd = "cut -c34- convhst.l2 >> convergenza.dat"
        ifail = system(execmd)
      elseif (NEO) then
!     .. can we do something similar with NEO?
      endif

 1000 continue

! 100 format (12x,a)
! 120 format (12x,i6)

      write(*,*)'End program'

      stop
      end
