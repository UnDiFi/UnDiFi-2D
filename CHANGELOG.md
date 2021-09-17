History of changes among different versions of **UnDiFi-2D**
============================================================

## Version v1.4
- added `compile_all.sh` script for automatic compilation, linking and installation of all required packages
- unified steady and unsteady code algorithm
- added shock-expansion and shock-vortex interaction unsteady testcases
- added a "testcase" argument in `calc_vel` and `mv_grid`
- substituted `start_run_x86.sh` with `run.sh` to run simulation
- Updated to EulFS-3.14 version
- Updated to petsc-3.14.6 version
- Update to new version of NEO
- Moved NEO in /UnDiFi-2D (same as EulFS)
- Bug for Intel 2019 when running with EulFS
  https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/New-MPI-error-with-Intel-2019-1-unable-to-run-MPI-hello-world/td-p/1158382
  BugFix: export FI_PROVIDER=tcp or export FI_SOCKETS_IFACE=eth0

## Version v1.3
- The following end-point conditions are implemented for shocks/discontinuities:
  SP    starting shock point originated by coalescence of characteristics
  TE    trailing edge point in supersonic flow (from the point two shocks and one contact discontinuity depart)
  FWP   normal shock point on generic wall (also curved)
  RR    regular reflection point on generic wall (also curved)
  PC    periodic connection point of a shock which crosses periodic boundary
- added in mv_dps an optional filter (see file input.dat) on the discontinuity velocity
- added freezing (from a fixed specified iteration in the input.dat file) of the grid topology

## Version v1.2
- added subroutine rd_dps_eq for the redistribution of shock points at the beginning and the simulation
- increased the number of digits in the IO files to reduce truncation errors
- implemented connection between two shocks condition:
   C connection between two extremal shock points (also of the same shock, see Q1D example)
- inserted additional hole points for non-simply connected geometries (i.e. cylinder and profile)
- file "input.dat" with data relative to:
    1) grid generation around the shock
    2) shock integration
    3) additional hole points
- inserted experimental version of the flux computation along the shock edge instead of along shock point

## Version v1.1
- added solution save on background grid (file na99.node)
  subroutine wtri0 and file shock/disc. (sh99.dat) subroutine wrt_sdw_info
- iterative TP solution initialized from previous solution
- manual restart from a specific iteration

## Version v1.0
Integrates shocks and contact discontinuities and the following interactions:
- shock-shock of the same family
- shock-shock of the opposite family
- shock-wall regular reflection
- shock-wall Mach reflection
Implements the following end-point conditions for shocks

  IPX   shock point entering from the domain on the X or Y side
  IPY

  OPX   shock point exiting from the domain on the X or Y side
  OPY

  WPNRX shock point normal to X or Y wall
  WPNRY

  TP    triple point

  QP    quadruple point

  RRX   regular reflection point on X wall
        (RRY not implemented)

  EP    end point in the supersonic zone
