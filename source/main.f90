program undifi_2d

  use mod_error, only: fatal
  use mod_run_external, only: run_external
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, neshmax, nprdbndmax, npshmax, nshmax, nspmax
  use mod_mesh, only: mesh_t
  implicit none(type, external)

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
  integer(i4) nin, nout
  parameter(nin=5, nout=6)

!     .. array definitions
  real(wp) xysh(ndim, npshmax, nshmax),&
  &xyshu(ndim, npshmax, nshmax),&
  &xyshd(ndim, npshmax, nshmax),&
  &zroeshuold(ndof, npshmax, nshmax),&
  &zroeshdold(ndof, npshmax, nshmax),&
  &norsh(ndim, npshmax, nshmax),&
  &wsh(ndim, npshmax, nshmax)

!     arrays for unsteady predictor-corrector time accurate integration
  real(wp) xyshnew(ndim, npshmax, nshmax),&
  &norshnew(ndim, npshmax, nshmax),&
  &wshnew(ndim, npshmax, nshmax),&
  &wshmean(ndim, npshmax, nshmax),&
  &zroeshuoldnew(ndof, npshmax, nshmax),&
  &zroeshdoldnew(ndof, npshmax, nshmax),&
  &varray(ndim, 30000)

  integer(i4) nodcodsh(npshmax, nshmax),&
  &nshocksegs(nshmax),&
  &nshockpoints(nshmax),&
  &shinspps(2, 5, nspmax),&
  &ispclr(5, nspmax)

  logical shtopolchanged
  logical neo, eulfs, su2

  character typespecpoints*5,&
  &typeshocks*1

  dimension typespecpoints(nspmax),&
  &typeshocks(nshmax)

!     .. scalar definition
  integer(i4) nShocks,&
  &nPhamPoints,&
  &nSpecPoints

  character execmd*255,&
  &fname*255,&
  &fname2*255,&
  &fnameback*255,&
  &backdir*255,&
  &bindir*255,&
  &gastype*4,&
  &hostype*10,&
  &color1*2,&
  &color2*2,&
  &VELFILE*18,&
  &MODE,&
  &ISPREDICTOR

!     .. local scalars ..
  integer(i4) i,&
  &nshockpointsold(nshmax), ish,&
  &nholes, totshockpoints, ii,&
  &nvt, ifail, nsteps, nbegin

!     background(0)/fitting(1)/backup(2) meshes -- see mod_mesh (issue #13)
  type(mesh_t) :: bkg, fit, bak
  integer(i4) nbfac_sh
  logical fndbnds

!     .. external functions ..
  integer(i4) initxdr
  external initxdr
  real(wp) rand
  external rand

!     .. external subroutines ..
  external calc_vel, co_norm, co_pnt_dspl, co_state_dps, dcopy, fltr_dls,&
  &fnd_phps, fx_dps_loc, fx_msh_sps, fx_state_dps, interp,&
  &interp_sp, mv_dps, mv_grid, rd_dps, rd_dps_eq, re_inp_data,&
  &re_sdw_info, readmesh, readpmap, solzne, wrt_sdw_info, wsh_mean,&
  &wtri, wtri0

!     Time steps for predictor-corrector
  real(wp) dtpr, dtco, nowtime

!     Read command line arguments
  integer(i4)           :: no, n_args
  character(len=20) :: testcase
  character(len=20) :: args(6)
  character(len=20) :: solvername
  logical           :: steady, unsteady

  n_args = command_argument_count(); 
  if (n_args /= 5 .and. n_args /= 6) then
    write (*, *) 'Usage: ../../bin/UnDiFi-2D_x86_64&
    &                                 0 501 false true "TestCaseName" [su2]'
    call fatal('wrong number of command-line arguments', 1)
  end if
  do i = 1, n_args
    call get_command_argument(i, args(i))
    args(i) = trim(adjustl(args(i)))
  end do
  read (args(1), *) nbegin
  read (args(2), *) nsteps
  read (args(3), *) eulfs
  read (args(4), *) steady
  read (args(5), *) testcase

!     optional 6th arg: a solver name that overrides args(3) when it
!     is a solver EULFS/NEO's boolean can't express (currently only
!     "su2"); omitted (n_args==5) is fully backward compatible with
!     every existing eulfs/neo caller (scripts/run_steady.sh, etc.)
  solvername = ""
  su2 = .false.
  if (n_args == 6) then
    solvername = args(6)
    su2 = (trim(solvername) == "su2")
  end if
  if (su2) then
    eulfs = .false.
  end if

  write (*, *) 'nbegin: ', nbegin
  write (*, *) 'nsteps: ', nsteps
  write (*, *) 'Use eulfs? ', eulfs
  write (*, *) 'Use su2? ', su2
  write (*, *) 'Is steady? ', steady
  write (*, *) 'testcase: ', testcase

!     flag to select the shock-capturing solver, eulfs, neo, or su2
  NEO = (.not. EULFS) .and. (.not. SU2)

!     flag to select the type of simulation: steady or unsteady
  UNSTEADY = (.not. STEADY)

! --------- set character variables

  bindir = "../../bin/"
  hostype = "x86_64"
!     hostype = 'i386'
  write (*, *) '1 ', hostype

!     gastype = "g140"
!     gastype = 'g167'

!     necessary for EulFS ALE
  velfile = "gridvel_000001.dat"
  mode = 'w'

! Vale
!     flag used in fx_dps_loc.f for the shock reflection on a wedge
  ShTopolChanged = .false.
! Vale

!     NDIM = 2
!     NDOF = 4 ! will be reset within rtri

! ---------- allocate space

!     write (nout,fmt=4000)

4000 format(/, /, ' memory allocation   ', /, ' ', 19('='),/)

  ifail = run_external("echo 'Running on' `uname -a` > triangle.log", 'echo', fatal_on_error=.false.)
  ifail = run_external("date >> triangle.log", 'date', fatal_on_error=.false.)

  if (eulfs) then
    execmd = "rm -fv convergenza.dat"
    ifail = run_external(execmd, 'rm')

!        copy file .petsrc in home
!        for UNSTEADY EulFS simulations this file
!        will be overwritten with .petsrc_predictor and
!        .petsrc_corrector in their respective steps
!        They only differ in the dt value
    if (.not. UNSTEADY) then
      execmd = "cp -fv .petscrc .petscrc"
    end if
!        ifail = system(execmd)
    if (ifail .ne. 0) call fatal('system command failed', ifail)
  end if

!     the initial grid is stored in a file called na00.1
  fname = "na00.1"
  fnameback = "na99"

! **********************************************************************
!  Read the original mesh, allocate space for the mesh variables and
!  allocate additional space for mesh points and segments
! **********************************************************************

  write (*, 1001, advance='no') 'readmesh               -->  '
  fndbnds = .true.
!     fndbnds=.false.
  call readmesh(bkg, fname, fndbnds)
  nvt = bkg%nvt
  write (*, 1002) ' ok'
1001 format(a)
1002 format(a)

! **********************************************************************
!  Read pmap
! **********************************************************************

  write (*, 1001, advance='no') 'readpmap               -->  '
  call readpmap(bkg)
  write (*, 1002) ' ok'

! **********************************************************************
!  Make backups of some of the arrays of the background mesh:
!  the nodal flag (nodcode), the boundary data structure (bndfac),
!  the boundary node pointer (nodptr)
! **********************************************************************

  write (*, 1001, advance='no') 'copy mesh(0) in mesh(2)-->  '
  bak%nbfac = bkg%nbfac
  bak%nbpoin = bkg%nbpoin
!     only the real (unpadded) entries are backed up -- the shock-edge/
!     shock-point tail that fx_msh_sps/fnd_phps write into bndfac/nodcod's
!     padded region is deliberately left alone by the restore below
  bak%bndfac = bkg%bndfac(:, 1:bkg%nbfac)
  bak%nodptr = bkg%nodptr
  bak%nodcod = bkg%nodcod(1:bkg%npoin)
  write (*, 1002) ' ok'

!     call x04eaf('general',' ',3,nbfac,istak(lbndfac(2)),3,
!    +            'bndry pointer(2) in main',ifail)
!     call x04eaf('general',' ',3,nbfac,bkg%bndfac,3,
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

  write (*, 1001, advance='no') 're_inp_data            -->  '
  call re_inp_data
  write (*, 1002) ' ok'

! **********************************************************************
!  Read the timesteps.dat file containing information about the dt
!  to be applied in the predictor/corrector steps (unsteady case)
!  Note: the timesteps.dat file should be already present in the
!        working directory
! **********************************************************************

  if (UNSTEADY) then

    write (*, 1001, advance='no') 're_dt_data             -->  '
    open (unit=12, file='timesteps.dat', status='old', action='read')
    read (12, *) dtpr
    read (12, *) dtco
    close (12)
    write (*, 1002) ' ok'
    write (*, '(1x,f7.4, 1x, f7.4)') dtpr, dtco

  end if ! UNSTEADY

! **********************************************************************
!  Read file sh00.dat containing information about shock/discontinuity
! **********************************************************************

  write (*, 1001, advance='no') 're_sdw_info            -->  '
  call re_sdw_info(&
  &xysh,&
  &bkg%zroe(1, bkg%npoin + 1),&                        !upstream state
  &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),&    !downstream state
  &zroeshuold,&
  &zroeshdold,&
  &bkg%nodcod(bkg%npoin + 1),&
  &bkg%xy,&   !vale
  &bkg%bndfac,& !vale
  &bkg%nbfac,&          !vale
  &bkg%npoin,&          !vale
  &nshocks,&
  &nshockpoints,&
  &nshocksegs,&
  &typeshocks,&
  &nspecpoints,&
  &typespecpoints,&
  &shinspps,&
  &ispclr)
  write (*, 1002) ' ok'

! **********************************************************************
!  Shock points (equally-spaced) redistribution strategy
!  Note: points redistribution may not be mandatory but beneficial
! **********************************************************************

  if (STEADY) then

    write (*, 1001, advance='no') 'rd_sps_eq              -->  '
    call rd_dps_eq(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &nshocks,&
    &nshockpoints,&
    &nshocksegs)
    write (*, 1002) ' ok'

    call dcopy(nshmax*npshmax*ndof, bkg%zroe(1, bkg%npoin + 1), 1,&
    &zroeshuold, 1)
    call dcopy(nshmax*npshmax*ndof,&
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax), 1,&
    &zroeshdold, 1)

  end if ! STEADY

!     call fx_sh_state(
!    +     dstak(lzroe(0)+bkg%npoin*ndof), ! upstream state
!    +     nshocks,
!    +     nshockpoints,
!    +     nshocksegs)

!     call pr_sh_state(
!    +     dstak(lzroe(0)+bkg%npoin*ndof),                     ! upstream state
!    +     dstak(lzroe(0)+bkg%npoin*ndof+nshmax*npshmax*ndof), ! downstream state
!    +     nshocks,
!    +     nshockpoints,
!    +     nshocksegs)

! **********************************************************************
!  Compute the normal vectors
! **********************************************************************

  if (UNSTEADY) then

    write (*, 1001, advance='no') 'co_norm                -->  '
    call co_norm(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &norsh,&
    &nshocks,&
    &nshockpoints,&
    &typeshocks,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr,&
    &bkg%ia,&
    &bkg%ja,&
    &bkg%iclr,&
    &bkg%nclr,&
    &bkg%xy)
    write (*, 1002) ' ok'

!       fix the normal orientation which otherwise
!       creates problems due to steady flow upstream
    if (testcase == 'ShockExpansion') then
      do no = 1, nshockpoints(1)
        norsh(1, no, 1) = -norsh(1, no, 1)
        norsh(2, no, 1) = -norsh(2, no, 1)
      end do
    end if

! **********************************************************************
!  Compute the velocity of the shock (wsh)
! **********************************************************************

    write (*, 1001, advance='no') 'co_state_dps           -->  '
    call co_state_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &zroeshuold,&
    &zroeshdold,&
    &norsh,&
    &wsh,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &i)
    write (*, 1002) ' ok'

! **********************************************************************
!  Fix the states in the discontinuity points
! **********************************************************************

    write (*, 1001, advance='no') 'fx_state_dps           -->  '
    call fx_state_dps(&
    &xysh,&
    &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
    &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &zroeshuold,&
    &zroeshdold,&
    &norsh,&
    &wsh,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &i,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr,&
    &bkg%ia,&
    &bkg%ja,&
    &bkg%iclr,&
    &bkg%nclr,&
    &bkg%xy)
    write (*, 1002) ' ok'

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
  write (6, *)
  write (6, *) ' starting the time loop'
  write (6, *)
  do 1000 i = 1 + nbegin, nsteps + nbegin
    write (6, *) '***********************************'
    write (6, *) ' time level is ', i
    write (6, *) '                                   '

! **********************************************************************
!  Find the cells crossed by the shock and the phantom points and
!  update data structure on the boundary to take into account of the
!  presence of phantom points
! **********************************************************************

    write (*, 1001, advance='no') 'fnd_phps               -->  '
    call fnd_phps(&
    &bkg%nedge,&
    &bkg%bndfac,&
    &bkg%nbfac,&
    &bkg%celnod,&
    &nvt,&
    &bkg%nelem,&
    &bkg%xy,&
    &xysh,&
    &bkg%nodcod,&
    &bkg%npoin,&
    &bkg%nodptr,&
    &bkg%nbpoin,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &nphampoints,&
    &bkg%pmap)
    write (*, 1002) ' ok'

! **********************************************************************
!  Compute the normal unit vector to shocks and discontinuities
! **********************************************************************

    write (*, 1001, advance='no') 'co_norm                -->  '
    call co_norm(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &norsh,&
    &nshocks,&
    &nshockpoints,&
    &typeshocks,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr,&
    &bkg%ia,&
    &bkg%ja,&
    &bkg%iclr,&
    &bkg%nclr,&
    &bkg%xy)
    write (*, 1002) ' ok'

!     fix the normal orientation which otherwise
!     creates problems due to steady flow upstream
    if (testcase == 'ShockExpansion') then
      do no = 1, nshockpoints(1)
        norsh(1, no, 1) = -norsh(1, no, 1)
        norsh(2, no, 1) = -norsh(2, no, 1)
      end do
    end if

! **********************************************************************
!  Interp() updates values in the special points
! **********************************************************************

    if (STEADY) then

!       goto 2340
      write (*, 1001, advance='no') 'interp_sp              -->  '
!    +       nshockpointsold)
      call interp_sp(&
      &bkg%celnod,&
      &nvt,&
      &bkg%nelem,&
      &bkg%xy,&
      &bkg%zroe,&
      &xysh,&
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &norsh,&
      &bkg%npoin,&
      &nshocks,&
      &nshockpoints,&
      &typeshocks,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps)
      write (*, 1002) ' ok'

    end if ! STEADY

! **********************************************************************
!  Compute downstream and upstream mesh coordinate of each shock and
!  discontinuity point
! **********************************************************************

    write (*, 1001, advance='no') 'co_pnt_dspl            -->  '
    call co_pnt_dspl(&
    &xysh,&
    &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
    &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%nodcod(bkg%npoin + 1),&
    &norsh,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr)
    write (*, 1002) ' ok'

! **********************************************************************
!  Fix the mesh around the special points
! **********************************************************************

    write (*, 1001, advance='no') 'fx_msh_sps             -->  '
    call fx_msh_sps(&
    &bkg%bndfac,&
    &bkg%nodcod,&
    &bkg%nbfac,&
    &nbfac_sh,&
    &nvt,&
    &bkg%nelem,&
    &bkg%xy,&
    &xysh,&
    &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
    &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
    &bkg%npoin,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr)
    write (*, 1002) ' ok'

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

      write (*, 1001, advance='no') 'copy sh(0) in sh(1)    -->  '
      call dcopy(ndim*npshmax*nshmax, xysh, 1, xyshnew, 1)
      call dcopy(ndim*npshmax*nshmax, zroeshdold, 1, zroeshdoldnew, 1)
      call dcopy(ndim*npshmax*nshmax, zroeshuold, 1, zroeshuoldnew, 1)
      call dcopy(ndim*npshmax*nshmax, norsh, 1, norshnew, 1)
      call dcopy(ndim*npshmax*nshmax, wsh, 1, wshnew, 1)
      write (*, 1002) ' ok'

! **********************************************************************
!  Compute the grid velocity
!  Note: also the shock velocity vector is required
! **********************************************************************

      write (*, 1001, advance='no') 'calc_vel               -->  '
      call calc_vel(&
      &bkg%npoin,&
      &varray,&
      &dtpr,&
      &bkg%xy,&
      &wsh,&
      &i,&
      &'y',&
      &nowtime,&
      &testcase)
      write (*, 1002) ' ok'

! **********************************************************************
!  It gives to EulFS info about grid velocity needed for the ALE
! **********************************************************************

      if (EULFS) then
        write (*, 1001, advance='no') 'solzne                -->   '
        call solzne(&
        &velfile,&
        &varray,&
        &ndim,&
        &bkg%npoin + 2*npshmax*nshmax,&
        &mode)
        write (*, 1002) ' ok'
      end if ! EULFS

    end if ! END UNSTEADY

! **********************************************************************
!  Write new poly file including all mesh points except phantom points
! **********************************************************************

    write (fname(3:7), fmt="(i5.5)") i
    write (*, 1001, advance='no') 'wtri                   -->  '
    call wtri(&
    &bkg%bndfac,&
    &bkg%nbfac,&
    &nbfac_sh,&
    &bkg%celnod,&
    &nvt,&
    &bkg%xy,&
    &xysh,&
    &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
    &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
    &bkg%zroe,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &bkg%nodcod,&
    &bkg%nodcod(bkg%npoin + 1),&
    &bkg%npoin,&
    &fname(1:7),&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &nphampoints)
    write (*, 1002) ' ok'

! **********************************************************************
!  If use NEO and it is the 1st iteration, creates the neogrid0.grd
!  grid file in /NEO_data/input/
! **********************************************************************

    if (NEO) then ! NEO solver
      if (i == 1 + nbegin .and. testcase == "ShockExpansion") then
        write (*, 1001, advance='no') 'neogrid0               -->  '
        execmd = bindir(1:10)//'neogrid0'
        ifail = run_external(execmd, 'neogrid0')
        write (*, 1002) ' ok'
      end if
    end if

! **********************************************************************
!  Generate the new mesh
! **********************************************************************

!        write(6,*)
!        write(6,*)' meshing with triangle; input file is ',fname(1:7)
!        write(6,*)
    write (*, 1001, advance='no') 'triangle               -->  '

    execmd = bindir(1:10)//'triangle_'//hostype(1:6)//' -nep '&
    &//fname(1:7)//' > log/triangle.log'
    ifail = run_external(execmd, 'triangle')

    write (*, 1002) ' ok'

! **********************************************************************
!  Freeze mesh topology
! **********************************************************************

    if (STEADY) then

      if (imtf .ne. 0 .and. i .gt. imtf) then
        write (*, 1001, advance='no') 'mesh topology freezing -->  '

        fname2 = 'stepXXXXX/naXXXXX.1'
        write (fname2(5:9), fmt="(i5.5)") imtf
        write (fname2(13:17), fmt="(i5.5)") imtf
        execmd = 'cp '//fname2(1:19)//'.ele '//fname(1:7)//'.1.ele'
        ifail = run_external(execmd, 'cp')
        execmd = 'cp '//fname2(1:19)//'.neigh '//&
        &fname(1:7)//'.1.neigh'
        ifail = run_external(execmd, 'cp')
        execmd = 'cp '//fname2(1:19)//'.edge '//fname(1:7)//'.1.edge'
        ifail = run_external(execmd, 'cp')

        write (*, 1002) ' ok'

      end if

    end if ! STEADY

! ***********************************
    if (EULFS) then ! EulFS SOLVER
! ***********************************

! **********************************************************************
!  Convert the triangle files into a fmt readable by the code using:
!  echo na0x.1 | triangle2dat
! **********************************************************************

      write (*, 1001, advance='no') 'triangle2dat           -->  '
!        execmd = "echo " // fname(1:7)
!     +   // ".1 |" // bindir(1:10) // "triangle2dat_" // hostype(1:6)
!     +   // " > log/triangle2dat.log"
      if (nprdbnd .eq. 0) then                            ! for the cases without periodic BCs
        execmd = "printf '"//fname(1:7)&
        &//".1\nn'|"&
        &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
        &//" > log/triangle2dat.log"
      elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 1) then    ! for the cases with only one periodic boundary
        write (color1, fmt="(i2.2)") prdbndclr(1, 1)          ! with points having the same x
        write (color2, fmt="(i2.2)") prdbndclr(2, 1)
        execmd = "printf '"//fname(1:7)&
        &//".1\ny\n"&
        &//color1//"\n"&
        &//color2//"\nx'|"&
        &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
        &//" > log/triangle2dat.log"
      elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 2) then ! for the cases with only one periodic boundary
        write (color1, fmt="(i2.2)") prdbndclr(1, 1)           ! with points having the same y
        write (color2, fmt="(i2.2)") prdbndclr(2, 1)
        execmd = "printf '"//fname(1:7)&
        &//".1\ny\n"&
        &//color1//"\n"&
        &//color2//"\ny'|"&
        &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
        &//" > log/triangle2dat.log"

      else ! for cases with more thatn one periodic boundary
        write (*, *) ' case not implemented!'
      end if

!        write(*,*)execmd

!        execmd = "printf '" // fname(1:7)
!    +   // ".1\ny\n2\n4\nx'|" ! for case cascade
!    +   // ".1\ny\n1\n3\nx'|" ! for case nacapar2
!    +   // ".1\nn'|"          ! for cases w/o periodic BCs'
!    +   // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
!    +   // " > log/triangle2dat.log"
!         write(*,*)execmd

      ifail = run_external(execmd, 'triangle2dat')

      write (*, 1002) ' ok'

! ****************************
!  Run one step of EulFS code
! ****************************

      if (UNSTEADY) then
!          It runs the predictor step of the EulFS code (we need to use dt/2)
        execmd = "cp -f .petscrc_predictor .petscrc"
        ifail = run_external(execmd, 'cp')
      end if

      write (*, 1001, advance='no') 'eulfs                  -->  '
!        execmd = bindir(1:10) // "eulfs11.13_"
!     +                        //gastype(1:4)//"_"//hostype(1:6)
!     +  // " -itmax 1 > log/eulfs.log"
      execmd = bindir(1:10)//"EulFS_"//hostype(1:6)&
      &//" -itmax 1 > log/eulfs.log"

      ifail = run_external(execmd, 'eulfs')

      if (unsteady) then
        execmd = "cp step000001.dat file001.dat"
        ifail = run_external(execmd, 'cp')
      end if

!        execmd = "cp file003.dat file010.dat"
!        ifail = system(execmd)
!        if(ifail.ne.0)call exit(1)

      write (*, 1002) ' ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0x.1 | dat2triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.bak
! **********************************************************************

      write (*, 1001, advance='no') 'dat2triangle           -->  '
!        execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
!    +   // "dat2triangle_" // hostype(1:6)
!    +   // ">log/dat2triangle.log"
      execmd = "printf '"//fname(1:7)//".1' | "//bindir(1:10)&
      &//"dat2triangle-NEW-"//hostype(1:6)&
      &//">log/dat2triangle.log"
      ifail = run_external(execmd, 'dat2triangle')

      write (*, 1002) ' ok'

! ***********************************
    elseif (SU2) then ! SU2 SOLVER
! ***********************************

! **********************************************************************
!  Convert the triangle files into SU2's native mesh + restart state:
!  echo na0x.1 / su2case | triangle2su2
!  "su2case" is a fixed basename (not fname) so su2case.cfg's
!  MESH_FILENAME/SOLUTION_FILENAME/RESTART_FILENAME never have to
!  change across outer iterations even though the Triangle basename
!  (fname) does. No periodic-BC support yet (see source_utils/
!  triangle2su2/main.f) -- fine for now, CircularCylinder has none.
! **********************************************************************

      write (*, 1001, advance='no') 'triangle2su2           -->  '
      execmd = "printf '"//fname(1:7)&
      &//".1\nsu2case'|"&
      &//bindir(1:10)//"triangle2su2-"//hostype(1:6)&
      &//" > log/triangle2su2.log"
      ifail = run_external(execmd, 'triangle2su2')

      write (*, 1002) ' ok'

! **************************
!  Run one step of SU2 code
! **************************
!  su2case.cfg sets ITER=1 with RESTART_SOL=YES: one implicit step
!  per outer UNDIFI iteration, the SU2 analogue of EulFS's -itmax 1.

      write (*, 1001, advance='no') 'su2                    -->  '
      execmd = bindir(1:10)//"SU2_CFD"&
      &//" su2case.cfg > log/su2.log"

      ifail = run_external(execmd, 'su2')

      write (*, 1002) ' ok'

! **********************************************************************
!  Convert su2case's restart state back into triangle fmt:
!  echo na0x.1 / su2case | su22triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.BAK
! **********************************************************************

      write (*, 1001, advance='no') 'su22triangle           -->  '
      execmd = "printf '"//fname(1:7)&
      &//".1\nsu2case'|"&
      &//bindir(1:10)//"su22triangle-"//hostype(1:6)&
      &//" > log/su22triangle.log"
      ifail = run_external(execmd, 'su22triangle')

      write (*, 1002) ' ok'

! ***********************************
    elseif (NEO) then ! NEO SOLVER
! ***********************************

!      neogrid0 works only for the 1st iteration of ShockVortex
!      but in all other cases we need this conversion
      if (STEADY .or. i /= 1 + nbegin .or. testcase == "ShockVortex") then

        write (*, 1001, advance='no') 'na00xTovvvv            -->   '
        execmd = "echo "//fname(1:7)&
        &//".1 |"//bindir(1:10)//"na2vvvv"&
        &//" > log/na2vvvv.log"
        ifail = run_external(execmd, 'na2vvvv', fatal_on_error=.false.)
        write (*, 1002) 'ok'

      end if

! ****************************************************************
!       convert the triangle files into a format readable by the
!       code using: echo na0X.1 | triangle2dat
! ****************************************************************

      write (*, 1001, advance='no') 'triangle2grd           -->   '
      execmd = "echo "//fname(1:7)&
      &//".1 |"//bindir(1:10)//"triangle2grd"&
      &//" > log/triangle2grd.log"
      ifail = run_external(execmd, 'triangle2grd')
      write (*, 1002) 'ok'

!     ******************
      if (UNSTEADY) then
!     ******************

!       if it is the 1st iteration
!       **************************
        if (i == 1 + nbegin) then

          write (*, 1001, advance='no') 'NEO 1st iteration      -->  '

          execmd = bindir(1:10)//"CRD_euler"&
          &//"> log/neo.log"
          ifail = run_external(execmd, 'NEO (1st iteration)')

          execmd = "cp ./NEO_data/output/vvvv.dat "//&
          &"./NEO_data/output/vvvv0.dat "
          ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

          execmd = "mv ./NEO_data/output/vvvv.dat "//&
          &"./NEO_data/output/vvvv_input.dat "
          ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

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

          execmd = "mv ./NEO_data/textinput/inputfile-exp.txt "//&
          &"./NEO_data/textinput/inputfile-exp.txt.BAK "
          ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

          execmd = "cp inputfile-exp.txt "//"./NEO_data/textinput/"
          ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

          write (*, 1002) ' ok'

        end if ! 1ST ITERATION

      end if ! UNSTEADY

!     for all the other iterations
!     ****************************
      write (*, 1001, advance='no') 'NEO                    -->   '
      execmd = bindir(1:10)//"CRD_euler"&
      &//"> log/neo.log"
      ifail = run_external(execmd, 'neo')
      write (*, 1002) 'ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0X.1 | NEO2triangle ex (dat2triangle)
!  Note: the file na0X.1.node will be overwritten with the values
!        updated by the code and a copy with "old" values is copied in
!        na0X.1.node.BAK
! **********************************************************************

      write (*, 1001, advance='no') 'NEO2triangle           -->   '
      execmd = "echo "//fname(1:7)//".1 | "//bindir(1:10)&
      &//"NEO2triangle"//">log/NEO2triangle.log"
      ifail = run_external(execmd, 'neo2triangle')
      write (*, 1002) 'ok'

    else

      write (*, *) 'should be running either EULFS ', eulfs, ' or NEO ',&
      &neo
      error stop 10

    end if ! IF-THEN-ELSE ON THE CFD CODE

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

      write (*, 1001, advance='no') 'readmesh               -->  '
      fname(1:9) = fname(1:7)//".1"
      fndbnds = .false.
      call readmesh(fit, fname, fndbnds)
      write (*, 1002) ' ok'

      ! TODO: check whether FX_USTATE should be added here ...

! **********************************************************************
!  Update the nodal values on the backgroud grid (0) using values of
!  the shocked grid (1); the shocked grid contains "wrong" values in
!  the phantom nodes but these will be changed at a later stage in
!  interp() we need to perform this copy here, since the shockmov()
!  routine works on nodal values of grid (0)
! **********************************************************************

      totshockpoints = 2*nshmax*npshmax

      write (*, 1001, advance='no') 'zroe(1)->zroe(0)       -->  '

      if (fit%npoin .eq. (bkg%npoin + totshockpoints)) then

        call dcopy(ndof*fit%npoin, fit%zroe, 1, bkg%zroe, 1)
        write (*, 1002) ' ok'

      else

!         the nof gridpoints in grid(1) must equal the number of
!         gridpoints on the background mesh + 2 * nshockpoints

        write (6, *) 'there is a mismatch in the nof gridpoints'
        write (6, *) 'btw grid(0) and grid(1)'
        write (*, *) bkg%npoin, totshockpoints
        write (*, *) fit%npoin, totshockpoints
        error stop 1

      end if

! **********************************************************************
!  Updates nodal values in all the shock points of grid (0) using R-H
!  relations and compute the shock speed
!  Note: xysh coordinates are used only to write tecplot file but not
!        elsewhere
! **********************************************************************

      write (*, 1001, advance='no') 'co_state_dps           -->  '
      call co_state_dps(&
      &xyshnew,&
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuoldnew,&
      &zroeshdoldnew,&
      &norshnew,&
      &wshnew,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &i)
      write (*, 1002) ' ok'

! **********************************************************************

      write (*, 1001, advance='no') 'fx_state_dps           -->  '
      call fx_state_dps(&
      &xyshnew,&                                           ! not used
      &bkg%xy(1, bkg%npoin + 1),&                     ! upstream   coord.
      &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream   state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuoldnew,&
      &zroeshdoldnew,&
      &norshnew,&
      &wshnew,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &i,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps,&
      &ispclr,&
      &bkg%ia,&
      &bkg%ja,&
      &bkg%iclr,&
      &bkg%nclr,&
      &bkg%xy)
      write (*, 1002) ' ok'

! **********************************************************************

      write (*, 1001, advance='no') 'zroesh(0)->zroesh(1)   -->  '
      call dcopy(ndof*totshockpoints,&
      &bkg%zroe(1, bkg%npoin + 1), 1,&
      &fit%zroe(1, bkg%npoin + 1), 1)
      write (*, 1002) ' ok'

! **********************************************************************

      write (*, 1001, advance='no') 'calc_vel               -->  '
      call calc_vel(&
      &bkg%npoin,&
      &varray,&
      &dtco,&
      &bkg%xy,&
      &wsh,& !WSHnew?
      &i,&
      &'n',&
      &nowtime,&
      &testcase)
      write (*, 1002) ' ok'

! **********************************************************************
!  It gives to eulfs information about grid velocity (corrector step)
! **********************************************************************

      if (EULFS) then
        write (*, 1001, advance='no') 'solzne                 -->   '
        call solzne(&
        &velfile,&
        &varray,&
        &ndim,&
        &bkg%npoin + 2*npshmax*nshmax,&
        &mode)
        write (*, 1002) ' ok'
      end if

! **********************************************************************
!  It generates the new mesh
! **********************************************************************

      write (*, 1001, advance='no') 'triangle               -->  '
      execmd = bindir(1:10)//'triangle_'//hostype(1:6)//' -nep '&
      &//fname(1:7)//' > log/triangle.log'
      ifail = run_external(execmd, 'triangle')
      write (*, 1002) ' ok'

! ***********************************
      if (EULFS) then ! EulFS SOLVER
! ***********************************

! **********************************************************************
!  Convert the triangle files into a fmt readable by the code using:
!  echo na0x.1 | triangle2dat
! **********************************************************************

        write (*, 1001, advance='no') 'triangle2dat           -->  '
!        execmd = "echo " // fname(1:7)
!     +   // ".1 |" // bindir(1:10) // "triangle2dat_" // hostype(1:6)
!     +   // " > log/triangle2dat.log"
        if (nprdbnd .eq. 0) then                            ! for the cases without periodic BCs
          execmd = "printf '"//fname(1:7)&
          &//".1\nn'|"&
          &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
          &//" > log/triangle2dat.log"
        elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 1) then    ! for the cases with only one periodic boundary
          write (color1, fmt="(i2.2)") prdbndclr(1, 1)          ! with points having the same x
          write (color2, fmt="(i2.2)") prdbndclr(2, 1)
          execmd = "printf '"//fname(1:7)&
          &//".1\ny\n"&
          &//color1//"\n"&
          &//color2//"\nx'|"&
          &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
          &//" > log/triangle2dat.log"
        elseif (nprdbnd .eq. 1 .and. prdbndclr(3, 1) .eq. 2) then ! for the cases with only one periodic boundary
          write (color1, fmt="(i2.2)") prdbndclr(1, 1)           ! with points having the same y
          write (color2, fmt="(i2.2)") prdbndclr(2, 1)
          execmd = "printf '"//fname(1:7)&
          &//".1\ny\n"&
          &//color1//"\n"&
          &//color2//"\ny'|"&
          &//bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)&
          &//" > log/triangle2dat.log"

        else ! for cases with more thatn one periodic boundary
          write (*, *) ' case not implemented!'
        end if

!        write(*,*)execmd

!        execmd = "printf '" // fname(1:7)
!    +   // ".1\ny\n2\n4\nx'|" ! for case cascade
!    +   // ".1\ny\n1\n3\nx'|" ! for case nacapar2
!    +   // ".1\nn'|"          ! for cases w/o periodic BCs'
!    +   // bindir(1:10)//"triangle2dat-NEW-"//hostype(1:6)
!    +   // " > log/triangle2dat.log"
!         write(*,*)execmd

        ifail = run_external(execmd, 'triangle2dat')

        write (*, 1002) ' ok'

! ****************************
!  Run one step of EulFS code
! ****************************

!        if (UNSTEADY) then
!          It runs the corrector step of the EulFS code (now we use the full dt)
        execmd = "cp -f .petscrc_corrector .petscrc"
        ifail = run_external(execmd, 'cp')
!        end if

        write (*, 1001, advance='no') 'eulfs                  -->  '
!        execmd = bindir(1:10) // "eulfs11.13_"
!     +                        //gastype(1:4)//"_"//hostype(1:6)
!     +  // " -itmax 1 > log/eulfs.log"
        execmd = bindir(1:10)//"EulFS_"//hostype(1:6)&
        &//" -itmax 1 > log/eulfs.log"

        ifail = run_external(execmd, 'eulfs')

!        if (UNSTEADY) then
        execmd = "cp step000001.dat file001.dat"
        ifail = run_external(execmd, 'cp')
!        endif

!        execmd = "cp file003.dat file010.dat"
!        ifail = system(execmd)
!        if(ifail.ne.0)call exit(1)

        write (*, 1002) ' ok'

! **********************************************************************
!  Convert the code files into triangle fmt using:
!  echo na0x.1 | dat2triangle
!  The file na0x.1.node will be overwritten with the values updated by
!  the code and a copy with "old" values is copied in na0x.1.node.bak
! **********************************************************************

        write (*, 1001, advance='no') 'dat2triangle           -->  '
!        execmd = "echo "//fname(1:7)//".1 | "// bindir(1:10)
!    +   // "dat2triangle_" // hostype(1:6)
!    +   // ">log/dat2triangle.log"
        execmd = "printf '"//fname(1:7)//".1' | "//bindir(1:10)&
        &//"dat2triangle-NEW-"//hostype(1:6)&
        &//">log/dat2triangle.log"
        ifail = run_external(execmd, 'dat2triangle')

        write (*, 1002) ' ok'

      elseif (NEO) then ! NEO SOLVER for UNSTEADY predictor step

! **********************************************************************
!   Update of vvvv.dat (input) for NEO
! **********************************************************************

        write (*, 1001, advance='no') 'na2vvvv                -->  '
        execmd = "echo "//fname(1:7)&
        &//".1 |"//bindir(1:10)//"na2vvvv"&
        &//" > log/na2vvvv.log"
        ifail = run_external(execmd, 'na2vvvv', fatal_on_error=.false.)
        write (*, 1002) ' ok'

! **********************************************************************
!  convert the triangle files into a fmt readable by the code
!  using: echo na0x.1 | triangle2dat
! **********************************************************************

        write (*, 1001, advance='no') 'triangle2grd           -->  '
        execmd = "echo "//fname(1:7)&
        &//".1 |"//bindir(1:10)//"triangle2grd"&
        &//" > log/triangle2grd.log"
        ifail = run_external(execmd, 'triangle2grd')
        write (*, 1002) ' ok'

! **********************************************************************
!  One step with NEO
! **********************************************************************

        write (*, 1001, advance='no') 'NEO                    -->  '
        execmd = bindir(1:10)//"CRD_euler"&
        &//" > log/neo.log"
        ifail = run_external(execmd, 'NEO')
        write (*, 1002) ' ok'

! **********************************************************************
!  Convert dat file to triangle file
! **********************************************************************

        write (*, 1001, advance='no') 'NEO2triangle           -->  '
        execmd = "echo "//fname(1:7)//".1 | "//bindir(1:10)&
        &//"NEO2triangle"//">log/NEO2triangle.log"
        ifail = run_external(execmd, 'NEO2triangle')
        write (*, 1002) ' ok'

      end if ! end SOLVER (EULFS/NEO) for UNSTEADY (corrector step)

! **********************************************************************

      write (*, 1001, advance='no') 'mv_grid                -->  '
      call mv_grid(&
      &bkg%npoin,&
      &varray,&
      &dtco,&
      &bkg%xy,&
      &wshnew,&
      &i,&
      &testcase) ! as in calc_vel, added arg to switch case
      write (*, 1002) ' ok'

! **********************************************************************

!       write(*,1001,advance='no')'copy sh(0) in sh(1)    -->  '
!       call dcopy(ndim*npshmax*nshmax,xysh,1,xyshnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshdold,1,zroeshdoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshuold,1,zroeshuoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,norsh,1,norshnew,1)
      call dcopy(ndim*npshmax*nshmax, wshnew, 1, wsh, 1)
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

    write (*, 1001, advance='no') 'readmesh               -->  '
    fname(1:9) = fname(1:7)//".1"
    fndbnds = .false.
    call readmesh(fit, fname, fndbnds)
    write (*, 1002) ' ok'

! **********************************************************************
!  TODO: comment this procedure
! **********************************************************************

!     call pr_sh_state(
!    +   dstak(lzroe(0)+bkg%npoin*ndof),                     ! upstream state
!!   +   dstak(lzroe(0)+bkg%npoin*ndof+nshmax*npshmax*ndof), ! downstream state
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

    totshockpoints = 2*nshmax*npshmax

    write (*, 1001, advance='no') 'zroe(1)->zroe(0)       -->  '
    if (fit%npoin .eq. (bkg%npoin + totshockpoints)) then
      call dcopy(ndof*fit%npoin, fit%zroe, 1,&
      &bkg%zroe, 1)

      write (*, 1002) ' ok'

    else

!       the nof gridpoints in grid(1) must equal the number of
!       gridpoints on the background mesh + 2 x nshockpoints

      write (6, *) 'there is a mismatch in the nof gridpoints'
      write (6, *) 'btw grid(0) and grid(1)'
      write (*, *) bkg%npoin, totshockpoints
      write (*, *) fit%npoin, totshockpoints
      error stop 1
    end if

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
!     call chktpnt2(bkg%bndfac,bkg%nbfac,nbfac_sh,
!    &     bkg%celnod,nvt,
!    &     bkg%nelem,
!    &     bkg%xy,
!    &     dstak(lcorg(0)+bkg%npoin*ndim),
!    &     dstak(lcorg(1)+bkg%npoin*ndim),
!    &     dstak(lzroe(0)+bkg%npoin*ndof),ndof,
!    &     ndim,
!    &     bkg%npoin,
!    &     nshocks,nshockpoints,nshocksegs)

! **********************************************************************
!  Update nodal values in all the shock points of grid (0) using R-H
!  relations and compute the shock speed
! **********************************************************************

    write (*, 1001, advance='no') 'co_state_dps           -->  '
    call co_state_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &zroeshuold,& !ZROESHuOLDnew?
    &zroeshdold,& !ZROESHdOLDnew?
    &norsh,&      !NORSHnew?
    &wsh,&        !WSHnew?
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &i)
    write (*, 1002) ' ok'

! **********************************************************************
!  Fix and correct the nodal values and shock speed in all special
!  point using the correct s-s interaction relation
! **********************************************************************

    write (*, 1001, advance='no') 'fx_state_dps           -->  '
    call fx_state_dps(&
    &xysh,& !XYSHnew?
    &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
    &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &zroeshuold,& !ZROESHuOLDnew?
    &zroeshdold,& !ZROESHdOLDnew?
    &norsh,&      !NORSHnew?
    &wsh,&        !WSHnew?
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &i,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr,&
    &bkg%ia,&
    &bkg%ja,&
    &bkg%iclr,&
    &bkg%nclr,&
    &bkg%xy)
    write (*, 1002) ' ok'

! **********************************************************************
!  Update the nodal values of shocks on the grid (1)
!  Note: we need to perform this copy here, since the fx_state_sps and
!        co_state_dps routines work on nodal values of grid (0)
! **********************************************************************

    write (*, 1001, advance='no') 'zroesh(0)->zroesh(1)   -->  '
    call dcopy(ndof*totshockpoints,&
    &bkg%zroe(1, bkg%npoin + 1), 1,&
    &fit%zroe(1, bkg%npoin + 1), 1)
    write (*, 1002) ' ok'

! **********************************************************************
!  Calculate the mean shock velocity
! **********************************************************************

    if (UNSTEADY) then
!     ******************

      write (*, 1001, advance='no') 'wsh_mean               -->  '
      call wsh_mean(&
      &wsh,&
      &wshnew,&
      &wshmean)
      write (*, 1002) ' ok'

    end if ! UNSTEADY

! **********************************************************************
!  Move the shocks
! **********************************************************************

    write (*, 1001, advance='no') 'mv_dps                 -->  '
    call mv_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &wsh,&
    &i,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks)
    write (*, 1002) ' ok'

! **********************************************************************
!  Fix and correct the nodal position in all special point
! **********************************************************************

    if (STEADY) then
!     ****************

      write (*, 1001, advance='no') 'fx_dps_loc             -->  '
      call fx_dps_loc(&
      &xysh,&
      &bkg%xy(1, bkg%npoin + 1),&                     ! upstream coord.
      &bkg%xy(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream coord.
      &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &zroeshuold,&
      &zroeshdold,&
      &norsh,&
      &wsh,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks,&
      &i,&
      &nspecpoints,&
      &typespecpoints,&
      &shinspps,&
      &ispclr,&
      &bkg%ia,&
      &bkg%ja,&
      &bkg%iclr,&
      &bkg%nclr,&
      &bkg%zroe,& ! vale
      &bkg%xy,&
      &shtopolchanged)  ! vale
      write (*, 1002) ' ok'

! **********************************************************************
!  Filters the shocks
! **********************************************************************

      write (*, 1001, advance='no') 'fltr_dls               -->  '
      call fltr_dls(&
      &xysh,&
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &wsh,&
      &i,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks)
      write (*, 1002) ' ok'

    end if ! STEADY

! **********************************************************************
!  Interp() updates values in the phantom nodes of the background
!  mesh (0) using values in the shocked mesh (1)
!  Note: npoin is passed as array
! **********************************************************************

!      goto 2340
    write (*, 1001, advance='no') 'interp                 -->  '
!    +      nshockpointsold)
    call interp(&
    &fit%bndfac,&
    &fit%nbfac,&
    &fit%celnod,&
    &nvt,&
    &fit%nelem,&
    &fit%xy,&
    &fit%zroe,&
    &xysh,&
    &fit%xy(1, bkg%npoin + 1),&                     !upstream coord.
    &fit%xy(1, bkg%npoin + 1 + nshmax*npshmax),& !downstream coord.
    &nphampoints,&
    &bkg%xy,&
    &bkg%zroe,&
    &bkg%nodcod,&
    &bkg%npoin,&
    &nshocks,&
    &nshockpoints,&
    &fit%ia,&
    &fit%ja,&
    &bkg%iclr,&
    &bkg%nclr)
    write (*, 1002) ' ok'

2340 continue

! **********************************************************************
!  Redistribute the shock nodes
!  Note: this redistribution is different from the previous one
! **********************************************************************

!     if(i.gt.1000)  goto 3450
    write (*, 1001, advance='no') 'rd_dps                 -->  '
    call rd_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &nshocks,&
    &nshockpoints,&
    &nshocksegs)
    write (*, 1002) ' ok'

3450 continue

! **********************************************************************
!  Redistribute the shock nodes
!  Note: shock node redistribution may be not necessary but beneficial
! **********************************************************************

!     if (mod(i,17).ne.0) goto 3451
!       write(*,1001)'rd_sps_eq              -->   '
!       call rd_dps_eq(
!    +       xysh,
!    +       dstak(lzroe(0)+bkg%npoin*ndof),                     ! upstream state
!    +       dstak(lzroe(0)+bkg%npoin*ndof+nshmax*npshmax*ndof), ! downstream state
!    +       nshocks,
!    +       nshockpoints,
!    +       nshocksegs)

!       write(*,1002)'rd_sps_eq              --> ok'

3451 continue

! **********************************************************************
!  Write the triangle node file of the background with the updated
!  value of zroe
! **********************************************************************

    write (*, 1001, advance='no') 'wtri0                  -->  '
!    +     ndim,
!    +     ndof,
    call wtri0(&
    &bkg%xy,&
    &bkg%zroe,&
    &bkg%nodcod,&
    &bkg%npoin,&
    &fnameback(1:4))
    write (*, 1002) ' ok'

! **********************************************************************
!  Write file sh99.dat containing information about shock/discontinuity
! **********************************************************************

    write (*, 1001, advance='no') 'wrt_sdw_info           -->  '
    call wrt_sdw_info(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &bkg%nodcod(bkg%npoin + 1),&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks,&
    &nspecpoints,&
    &typespecpoints,&
    &shinspps,&
    &ispclr)
    write (*, 1002) ' ok'

! **********************************************************************
!  Copy new shock on the old one
! **********************************************************************

    if (UNSTEADY) then

      write (*, 1001, advance='no') 'copy sh(1) in sh(0)    -->  '
      call dcopy(ndim*npshmax*nshmax, zroeshdoldnew, 1, zroeshdold, 1)
      call dcopy(ndim*npshmax*nshmax, zroeshuoldnew, 1, zroeshuold, 1)
      call dcopy(ndim*npshmax*nshmax, norshnew, 1, norsh, 1)
      call dcopy(ndim*npshmax*nshmax, wshnew, 1, wsh, 1)
      write (*, 1002) ' ok'

    end if ! UNSTEADY

!     call solzne("file004.dat",bkg%zroe,ndof,bkg%npoin,"w")

! **********************************************************************
!  Restore the original arrays of the background grid
! **********************************************************************

    bkg%nbfac = bak%nbfac
    bkg%nbpoin = bak%nbpoin
    bkg%bndfac(:, 1:bak%nbfac) = bak%bndfac
    bkg%nodptr = bak%nodptr
    bkg%nodcod(1:bkg%npoin) = bak%nodcod

! **********************************************************************
!  Release all arrays allocated for the shocked mesh (1)
! **********************************************************************

    call fit%free()

! **********************************************************************
!  Create a directory to backup files
! **********************************************************************

    write (backdir(5:9), fmt="(i5.5)") i
    if (mod(i - 1, ibak) .eq. 0) then
      execmd = "mkdir -v "//backdir(1:9)
      ifail = run_external(execmd, 'mkdir', fatal_on_error=.false.)
!     execmd = "mv shocknor.dat shock.log file00[1-3].dat file010.dat fs
!    &pl.out shocks.dat "//fname(1:7)//".* "//backdir(1:9)
      if (eulfs) then
        execmd = "mv -v shocknor.dat file00[1-4].dat file010.dat sh&
        &99.dat  "//fname(1:7)//".* "//fnameback(1:4)//".node&
        &                               "//backdir(1:9)
        ifail = run_external(execmd, 'mv', fatal_on_error=.false.)
      elseif (neo) then
        execmd = "mv -v shocknor.dat sh99.dat&
        &                                        "//fname(1:7)//".* "//fnameback(1:4)//".node "//&
        &backdir(1:9)
        ifail = run_external(execmd, 'mv', fatal_on_error=.false.)
        execmd =&
        &"cp -vp ./NEO_data/input/neogrid.grd&
        &                              ./NEO_data/input/vel.dat "//backdir(1:9)

!    &       "mv -v ../../source_utils/NEO_source/input/neogrid.grd
!    &       ../../source_utils/NEO_source/input/vel.dat "//backdir(1:9)

! Note: make attention with cp/mv because of vel.dat file

        ifail = run_external(execmd, 'cp', fatal_on_error=.false.)
!
!         if (i == 1+nbegin) then
!            execmd = "mv ../../../NEO_source/output/vvvv0.dat "
!    &             //backdir(1:9)
!            ifail = system(execmd)
!         end if
!
        execmd =&
        &"cp -v ./NEO_data/output/vvvv.dat "//backdir(1:9)
        ifail = run_external(execmd, 'cp', fatal_on_error=.false.)

        execmd =&
        &"mv -v ./NEO_data/output/vvvv_input.dat "//backdir(1:9)
        ifail = run_external(execmd, 'mv', fatal_on_error=.false.)

! copy the solution fie in Tec folder
!         if (i == 1+nbegin) then
!            execmd = "cp "//backdir(1:9)//"/vvvv0.dat "
!    &            //"tec/vvvv00000.dat"
!            ifail = system(execmd)
!         end if

!         execmd = "cp ../../source_utils/NEO_source/output/vvvv.dat "//"Tec/vvvv"
!    &         //backdir(5:9)//".dat"
!         ifail = system(execmd)

      end if
!         ifail = system(execmd)
    else ! backing up or not ...
!     execmd = "rm shocknor.dat shock.log file00[1-3].dat file010.dat fs
!    &pl.out shocks.dat "//fname(1:7)//".*"
      if (EULFS) then
        execmd = "rm shocknor.dat file00[1-4].dat file010.dat&
        &                               "//fname(1:7)//".* "//fnameback(1:4)//".node "//&
        &"sh99.dat "
      elseif (NEO) then
        execmd = "rm shocknor.dat "//fname(1:7)//".* "//fnameback(1:4)//&
        &".node "//"sh99.dat "
      end if
      ifail = run_external(execmd, 'rm', fatal_on_error=.false.)
    end if

    if (EULFS) then
      execmd = "cut -c34- convhst.l2 >> convergenza.dat"
      ifail = run_external(execmd, 'cut', fatal_on_error=.false.)
    elseif (NEO) then
!     .. can we do something similar with NEO?
    end if

1000 continue

! 100 format (12x,a)
! 120 format (12x,i6)

    write (*, *) 'End program'

    stop
    end program undifi_2d
