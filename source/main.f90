program undifi_2d

  use mod_error, only: fatal
  use mod_run_external, only: run_external
  use mod_kinds, only: wp, i4
  use mod_constants, only: naddholesmax, ndim, ndof, neshmax, nprdbndmax, npshmax, nshmax, nspmax
  use mod_mesh, only: mesh_t
  use mod_shock_system, only: xysh, zroeshuold, zroeshdold, norsh, wsh,&
  &xyshnew, norshnew, wshnew, wshmean, zroeshuoldnew, zroeshdoldnew,&
  &shock_system_init, shock_system_refresh
  use mod_shock_advance, only: advance
  use mod_solver_iface, only: flow_solver_t, solver_run_ctx_t
  use mod_solver_registry, only: make_flow_solver
  use mod_timer, only: timer_tic, timer_toc
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
!     xysh, zroeshuold, zroeshdold, norsh, wsh, nodcodsh and their *new
!     predictor/corrector shadow counterparts now live in mod_shock_system
!     (Phase 3.1, ROADMAP.md #14) -- use-associated above.
  real(wp) varray(ndim, 30000)

  integer(i4) nshocksegs(nshmax),&
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
  &VELFILE*18,&
  &MODE,&
  &ISPREDICTOR

!     .. local scalars ..
  integer(i4) i,&
  &nshockpointsold(nshmax), ish,&
  &nholes, ii,&
  &nvt, ifail, nsteps, nbegin

!     background(0)/fitting(1)/backup(2) meshes -- see mod_mesh (issue #13)
!     target: bkg%zroe backs disc(:)%zu/%zd's pointer association, see
!     mod_shock_system (Phase 3.1, ROADMAP.md #14)
  type(mesh_t), target :: bkg, fit, bak
  integer(i4) nbfac_sh
  logical fndbnds

!     Phase 3.7 (ROADMAP.md #14, issue #15): one polymorphic solver
!     replacing the if(EULFS)/elseif(SU2)/elseif(NEO) dispatch that used
!     to be inlined at 7 separate call sites -- see mod_solver_iface.f90.
  class(flow_solver_t), allocatable :: solver
  type(solver_run_ctx_t) :: ctx

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
  solver = make_flow_solver(eulfs, su2)

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

!     Phase 3.7 (ROADMAP.md #14, issue #15): was `if (eulfs) then ...
!     end if` inline here -- see mod_eulfs_solver.f90's eulfs_setup for
!     the ported logic (including a preserved pre-existing bug).
  call solver%setup(unsteady)

!     the initial grid is stored in a file called na00.1
  fname = "na00.1"
  fnameback = "na99"

! **********************************************************************
!  Read the original mesh, allocate space for the mesh variables and
!  allocate additional space for mesh points and segments
! **********************************************************************

  write (*, 1001, advance='no') 'readmesh               -->  '
  call timer_tic()
  fndbnds = .true.
!     fndbnds=.false.
  call readmesh(bkg, fname, fndbnds)
  nvt = bkg%nvt
  write (*, 1002) ' ok'//timer_toc()
1001 format(a)
1002 format(a)

! **********************************************************************
!  Read pmap
! **********************************************************************

  write (*, 1001, advance='no') 'readpmap               -->  '
  call timer_tic()
  call readpmap(bkg)
  write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Make backups of some of the arrays of the background mesh:
!  the nodal flag (nodcode), the boundary data structure (bndfac),
!  the boundary node pointer (nodptr)
! **********************************************************************

  write (*, 1001, advance='no') 'copy mesh(0) in mesh(2)-->  '
  call timer_tic()
  bak%nbfac = bkg%nbfac
  bak%nbpoin = bkg%nbpoin
!     only the real (unpadded) entries are backed up -- the shock-edge/
!     shock-point tail that fx_msh_sps/fnd_phps write into bndfac/nodcod's
!     padded region is deliberately left alone by the restore below
  bak%bndfac = bkg%bndfac(:, 1:bkg%nbfac)
  bak%nodptr = bkg%nodptr
  bak%nodcod = bkg%nodcod(1:bkg%npoin)
  write (*, 1002) ' ok'//timer_toc()

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
  call timer_tic()
  call re_inp_data
  write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Read the timesteps.dat file containing information about the dt
!  to be applied in the predictor/corrector steps (unsteady case)
!  Note: the timesteps.dat file should be already present in the
!        working directory
! **********************************************************************

  if (UNSTEADY) then

    write (*, 1001, advance='no') 're_dt_data             -->  '
    call timer_tic()
    open (unit=12, file='timesteps.dat', status='old', action='read')
    read (12, *) dtpr
    read (12, *) dtco
    close (12)
    write (*, 1002) ' ok'//timer_toc()
    write (*, '(1x,f7.4, 1x, f7.4)') dtpr, dtco

  end if ! UNSTEADY

! **********************************************************************
!  Read file sh00.dat containing information about shock/discontinuity
! **********************************************************************

  write (*, 1001, advance='no') 're_sdw_info            -->  '
  call timer_tic()
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
  write (*, 1002) ' ok'//timer_toc()

!     mod_shock_system's disc(:)/disc_new(:) (Phase 3.1, ROADMAP.md #14):
!     no consumer yet, kept in sync for Phase 3.2 to build on
  call shock_system_init(nshocks, bkg%npoin, bkg%zroe)
  call shock_system_refresh(nshockpoints, nshocksegs, typeshocks)

! **********************************************************************
!  Shock points (equally-spaced) redistribution strategy
!  Note: points redistribution may not be mandatory but beneficial
! **********************************************************************

  if (STEADY) then

    write (*, 1001, advance='no') 'rd_sps_eq              -->  '
    call timer_tic()
    call rd_dps_eq(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &nshocks,&
    &nshockpoints,&
    &nshocksegs)
    write (*, 1002) ' ok'//timer_toc()

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
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

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
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Fix the states in the discontinuity points
! **********************************************************************

    write (*, 1001, advance='no') 'fx_state_dps           -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

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
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Compute the normal unit vector to shocks and discontinuities
! **********************************************************************

    write (*, 1001, advance='no') 'co_norm                -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

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
      call timer_tic()
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
      write (*, 1002) ' ok'//timer_toc()

    end if ! STEADY

! **********************************************************************
!  Compute downstream and upstream mesh coordinate of each shock and
!  discontinuity point
! **********************************************************************

    write (*, 1001, advance='no') 'co_pnt_dspl            -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Fix the mesh around the special points
! **********************************************************************

    write (*, 1001, advance='no') 'fx_msh_sps             -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

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
      call timer_tic()
!     zroeshdold/zroeshuold are (ndof, npshmax, nshmax), not (ndim, ...)
!     like xysh/norsh/wsh -- these dcopy calls used to read
!     `ndim*npshmax*nshmax` (element count) for zroeshdold/zroeshuold
!     too. Since ndof is the fastest-varying dimension,
!     `ndim*npshmax*nshmax` elements is exactly (ndim/ndof)*nshmax = 5
!     whole per-shock slices (out of the nshmax=10 max), not "half of
!     every point's dof" -- so this silently left shock slots 6-10
!     uncopied (stale) on every UNSTEADY predictor step, for any case
!     with more than 5 shocks. Fixed here (Phase 4, ROADMAP.md #16,
!     found while scoping 4.4's dcopy cleanup) alongside 4.4 itself:
!     each array's declared bound is the compile-time max (npshmax=500
!     points, nshmax=10 shocks), but real cases use far fewer of
!     either, so a flat dcopy over the whole declared extent copies
!     mostly padding. Looping per real shock and copying only its real
!     nshockpoints(ish) removes that waste.
      do ish = 1, nShocks
        call dcopy(ndim*nshockpoints(ish), xysh(1, 1, ish), 1, xyshnew(1, 1, ish), 1)
        call dcopy(ndof*nshockpoints(ish), zroeshdold(1, 1, ish), 1, zroeshdoldnew(1, 1, ish), 1)
        call dcopy(ndof*nshockpoints(ish), zroeshuold(1, 1, ish), 1, zroeshuoldnew(1, 1, ish), 1)
        call dcopy(ndim*nshockpoints(ish), norsh(1, 1, ish), 1, norshnew(1, 1, ish), 1)
        call dcopy(ndim*nshockpoints(ish), wsh(1, 1, ish), 1, wshnew(1, 1, ish), 1)
      end do
      write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Compute the grid velocity
!  Note: also the shock velocity vector is required
! **********************************************************************

      write (*, 1001, advance='no') 'calc_vel               -->  '
      call timer_tic()
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
      write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  It gives to EulFS info about grid velocity needed for the ALE
! **********************************************************************

      if (solver%supports_ale()) then
        write (*, 1001, advance='no') 'solzne                -->   '
        call timer_tic()
        call solzne(&
        &velfile,&
        &varray,&
        &ndim,&
        &bkg%npoin + 2*npshmax*nshmax,&
        &mode)
        write (*, 1002) ' ok'//timer_toc()
      end if ! supports_ale (EULFS only)

    end if ! END UNSTEADY

! **********************************************************************
!  Write new poly file including all mesh points except phantom points
! **********************************************************************

    write (fname(3:7), fmt="(i5.5)") i

!     Phase 3.7 (ROADMAP.md #14, issue #15): populate the context passed
!     to every solver%... call for the rest of this outer iteration.
!     corrector is set .false. here and only flipped .true. around the
!     UNSTEADY-corrector-step calls further down.
    ctx%fname = fname(1:7)
    ctx%bindir = bindir(1:10)
    ctx%hostype = hostype(1:6)
    ctx%testcase = testcase
    ctx%iter = i
    ctx%nbegin = nbegin
    ctx%unsteady = unsteady
    ctx%corrector = .false.

    write (*, 1001, advance='no') 'wtri                   -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  If use NEO and it is the 1st iteration, creates the neogrid0.grd
!  grid file in /NEO_data/input/
! **********************************************************************

!     Phase 3.7 (ROADMAP.md #14, issue #15): was `if (NEO) then ... end
!     if` inline here -- see mod_neo_solver.f90's neo_pre_mesh_setup.
    call solver%pre_mesh_setup(ctx)

! **********************************************************************
!  Generate the new mesh
! **********************************************************************

!        write(6,*)
!        write(6,*)' meshing with triangle; input file is ',fname(1:7)
!        write(6,*)
    write (*, 1001, advance='no') 'triangle               -->  '
    call timer_tic()

    execmd = bindir(1:10)//'triangle_'//hostype(1:6)//' -nep '&
    &//fname(1:7)//' > log/triangle.log'
    ifail = run_external(execmd, 'triangle')

    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Freeze mesh topology
! **********************************************************************

    if (STEADY) then

      if (imtf .ne. 0 .and. i .gt. imtf) then
        write (*, 1001, advance='no') 'mesh topology freezing -->  '
        call timer_tic()

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

        write (*, 1002) ' ok'//timer_toc()

      end if

    end if ! STEADY

! ***********************************
!     Phase 3.7 (ROADMAP.md #14, issue #15): was the three-way
!     if(EULFS)/elseif(SU2)/elseif(NEO)/else-error-stop dispatch inline
!     here (~260 lines) -- see mod_eulfs_solver.f90/mod_neo_solver.f90/
!     mod_su2_solver.f90 for the ported per-solver logic. The `else`
!     fallback is dropped: eulfs/su2/neo are mutually exclusive by
!     construction (see the CLI parsing above), so make_flow_solver
!     always allocates exactly one of the three concrete types and this
!     branch was already structurally unreachable.
! ***********************************
    call solver%prepare(ctx)
    call solver%run(ctx)
    call solver%harvest(ctx)


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

      ! TODO: check whether FX_USTATE should be added here ...

! **********************************************************************
!  Update the nodal values on the backgroud grid (0) using values of
!  the shocked grid (1); the shocked grid contains "wrong" values in
!  the phantom nodes but these will be changed at a later stage in
!  interp() we need to perform this copy here, since the shockmov()
!  routine works on nodal values of grid (0). Updates nodal values in
!  all the shock points of grid (0) using R-H relations and compute the
!  shock speed. Prepares grid velocity for the corrector solve.
!  (mod_shock_advance.f90, ROADMAP.md #14 3.6 -- was duplicated inline
!  with the loop-end block below, differing only in which shadow-array
!  set is targeted and in this grid-velocity prep)
! **********************************************************************

      fname(1:9) = fname(1:7)//".1"
      call advance(bkg, fit, fname, i, nshocks, nshockpoints, nshocksegs,&
      &typeshocks, nspecpoints, typespecpoints, shinspps, ispclr,&
      &new_shadow=.true., eulfs=solver%supports_ale(), varray=varray,&
      &velfile=velfile,&
      &mode=mode, testcase=testcase, dt=dtco, velflag='n')

! **********************************************************************
!  It generates the new mesh
! **********************************************************************

      write (*, 1001, advance='no') 'triangle               -->  '
      call timer_tic()
      execmd = bindir(1:10)//'triangle_'//hostype(1:6)//' -nep '&
      &//fname(1:7)//' > log/triangle.log'
      ifail = run_external(execmd, 'triangle')
      write (*, 1002) ' ok'//timer_toc()

! ***********************************
!     Phase 3.7 (ROADMAP.md #14, issue #15): was the UNSTEADY-corrector
!     if(EULFS)/elseif(NEO) dispatch inline here (no SU2 arm existed --
!     su2_t's prepare/run/harvest are no-ops when ctx%corrector is
!     .true., reproducing that gap exactly, see mod_su2_solver.f90).
! ***********************************
      ctx%corrector = .true.
      call solver%prepare(ctx)
      call solver%run(ctx)
      call solver%harvest(ctx)
      ctx%corrector = .false.


! **********************************************************************

      write (*, 1001, advance='no') 'mv_grid                -->  '
      call timer_tic()
      call mv_grid(&
      &bkg%npoin,&
      &varray,&
      &dtco,&
      &bkg%xy,&
      &wshnew,&
      &i,&
      &testcase) ! as in calc_vel, added arg to switch case
      write (*, 1002) ' ok'//timer_toc()

! **********************************************************************

!       write(*,1001,advance='no')'copy sh(0) in sh(1)    -->  '
!       call dcopy(ndim*npshmax*nshmax,xysh,1,xyshnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshdold,1,zroeshdoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,zroeshuold,1,zroeshuoldnew,1)
!       call dcopy(ndim*npshmax*nshmax,norsh,1,norshnew,1)
      do ish = 1, nShocks
        call dcopy(ndim*nshockpoints(ish), wshnew(1, 1, ish), 1, wsh(1, 1, ish), 1)
      end do
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

! **********************************************************************
!  TODO: comment this procedure
! **********************************************************************

!     call pr_sh_state(
!    +   dstak(lzroe(0)+bkg%npoin*ndof),                     ! upstream state
!!   +   dstak(lzroe(0)+bkg%npoin*ndof+nshmax*npshmax*ndof), ! downstream state
!    +   nshocks,
!    +   nshockpoints,
!    +   nshocksegs)

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
!  Read the mesh generated by triangle, update the nodal values on the
!  background grid (0) from the shocked grid (1), update nodal values in
!  all the shock points of grid (0) using R-H relations, compute the
!  shock speed, fix and correct the nodal values/shock speed in all
!  special points, and copy the updated shock state back onto grid (1).
!  (mod_shock_advance.f90, ROADMAP.md #14 3.6 -- was duplicated inline
!  with the mid-loop block above; no grid-velocity prep follows here,
!  since there is no further solve this iteration)
! **********************************************************************

    fname(1:9) = fname(1:7)//".1"
    call advance(bkg, fit, fname, i, nshocks, nshockpoints, nshocksegs,&
    &typeshocks, nspecpoints, typespecpoints, shinspps, ispclr,&
    &new_shadow=.false., eulfs=solver%supports_ale(), varray=varray,&
    &velfile=velfile,&
    &mode=mode, testcase=testcase)

! **********************************************************************
!  Calculate the mean shock velocity
! **********************************************************************

    if (UNSTEADY) then
!     ******************

      write (*, 1001, advance='no') 'wsh_mean               -->  '
      call timer_tic()
      call wsh_mean(&
      &wsh,&
      &wshnew,&
      &wshmean)
      write (*, 1002) ' ok'//timer_toc()

    end if ! UNSTEADY

! **********************************************************************
!  Move the shocks
! **********************************************************************

    write (*, 1001, advance='no') 'mv_dps                 -->  '
    call timer_tic()
    call mv_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &wsh,&
    &i,&
    &nshocks,&
    &nshockpoints,&
    &nshocksegs,&
    &typeshocks)
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Fix and correct the nodal position in all special point
! **********************************************************************

    if (STEADY) then
!     ****************

      write (*, 1001, advance='no') 'fx_dps_loc             -->  '
      call timer_tic()
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
      write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Filters the shocks
! **********************************************************************

      write (*, 1001, advance='no') 'fltr_dls               -->  '
      call timer_tic()
      call fltr_dls(&
      &xysh,&
      &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
      &wsh,&
      &i,&
      &nshocks,&
      &nshockpoints,&
      &nshocksegs,&
      &typeshocks)
      write (*, 1002) ' ok'//timer_toc()

    end if ! STEADY

! **********************************************************************
!  Interp() updates values in the phantom nodes of the background
!  mesh (0) using values in the shocked mesh (1)
!  Note: npoin is passed as array
! **********************************************************************

!      goto 2340
    write (*, 1001, advance='no') 'interp                 -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

2340 continue

! **********************************************************************
!  Redistribute the shock nodes
!  Note: this redistribution is different from the previous one
! **********************************************************************

!     if(i.gt.1000)  goto 3450
    write (*, 1001, advance='no') 'rd_dps                 -->  '
    call timer_tic()
    call rd_dps(&
    &xysh,&
    &bkg%zroe(1, bkg%npoin + 1),&                     ! upstream state
    &bkg%zroe(1, bkg%npoin + 1 + nshmax*npshmax),& ! downstream state
    &nshocks,&
    &nshockpoints,&
    &nshocksegs)
    write (*, 1002) ' ok'//timer_toc()

!     mod_shock_system's disc(:)/disc_new(:) (Phase 3.1, ROADMAP.md #14):
!     rd_dps may have changed nshockpoints/nshocksegs above
    call shock_system_refresh(nshockpoints, nshocksegs, typeshocks)

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
    call timer_tic()
!    +     ndim,
!    +     ndof,
    call wtri0(&
    &bkg%xy,&
    &bkg%zroe,&
    &bkg%nodcod,&
    &bkg%npoin,&
    &fnameback(1:4))
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Write file sh99.dat containing information about shock/discontinuity
! **********************************************************************

    write (*, 1001, advance='no') 'wrt_sdw_info           -->  '
    call timer_tic()
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
    write (*, 1002) ' ok'//timer_toc()

! **********************************************************************
!  Copy new shock on the old one
! **********************************************************************

    if (UNSTEADY) then

      write (*, 1001, advance='no') 'copy sh(1) in sh(0)    -->  '
      call timer_tic()
!     see the matching comments at the mirror-direction copy above
!     (~line 649): zroeshdold/zroeshuold need ndof, not ndim, and every
!     array here only needs its real nshockpoints(ish) columns, not the
!     full npshmax-per-shock declared bound.
      do ish = 1, nShocks
        call dcopy(ndof*nshockpoints(ish), zroeshdoldnew(1, 1, ish), 1, zroeshdold(1, 1, ish), 1)
        call dcopy(ndof*nshockpoints(ish), zroeshuoldnew(1, 1, ish), 1, zroeshuold(1, 1, ish), 1)
        call dcopy(ndim*nshockpoints(ish), norshnew(1, 1, ish), 1, norsh(1, 1, ish), 1)
        call dcopy(ndim*nshockpoints(ish), wshnew(1, 1, ish), 1, wsh(1, 1, ish), 1)
      end do
      write (*, 1002) ' ok'//timer_toc()

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
    ctx%backdir = backdir(1:9)
    ctx%fnameback = fnameback(1:4)
!     Phase 3.7 (ROADMAP.md #14, issue #15): was `if (eulfs)/elseif
!     (neo)` (backing-up branch) and `if (EULFS)/elseif (NEO)`
!     (non-backing-up branch) inline here -- see mod_eulfs_solver.f90/
!     mod_neo_solver.f90's archive() (su2_t doesn't override archive();
!     SU2 never had an arm at either branch, see mod_su2_solver.f90).
    ctx%backing_up = (mod(i - 1, ibak) .eq. 0)
    if (ctx%backing_up) then
      execmd = "mkdir -v "//backdir(1:9)
      ifail = run_external(execmd, 'mkdir', fatal_on_error=.false.)
    end if
    call solver%archive(ctx)

!     Phase 3.7 (ROADMAP.md #14, issue #15): was `if (EULFS)/elseif
!     (NEO)` inline here -- see mod_eulfs_solver.f90's
!     eulfs_log_convergence (neo_t/su2_t inherit the shared no-op
!     default, exactly matching NEO's empty arm and SU2's absent one).
    call solver%log_convergence()

1000 continue

! 100 format (12x,a)
! 120 format (12x,i6)

    write (*, *) 'End program'

    stop
    end program undifi_2d
