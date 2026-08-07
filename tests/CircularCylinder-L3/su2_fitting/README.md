# CircularCylinder-1, SU2 shock-fitting run (all 501 outer iterations)

Runs SU2 through `source/main.f`'s outer UNDIFI loop
(`UnDiFi-2D_x86_64 0 501 false true CircularCylinder su2`), as opposed
to `../su2_capturing/`'s hand-driven single long SU2 run. Each outer
iteration calls `triangle2su2 -> SU2_CFD -> su22triangle` once
(`su2case.cfg`'s `ITER=1`/`RESTART_SOL=YES`, the SU2 analogue of
EulFS's `-itmax 1`), same freestream/marker setup as
`../su2_capturing/circularcylinder.cfg`.

Run from a scratch copy of this test case (`bindir="../../bin/"` in
`main.f` is hardcoded two levels below the UnDiFi-2D root, so it can't
be run from a subdirectory of `tests/CircularCylinder/` -- only
archived here after the fact, same as `su2_capturing/`).

## Bugs fixed to get this far

- **`wsu2.f` wrote every declared Triangle point into the `.su2`
  mesh, including thousands of unreferenced, coincident phantom-node
  placeholders** (marker `-99`) that UnDiFi's fixed-size phantom-node
  preallocation leaves in every `.node` file. EulFS/NEO tolerate the
  padding; SU2 doesn't -- `CPhysicalGeometry::DistributeColoring`
  requires every declared `NPOIN` entry to be referenced by some
  element and aborts otherwise ("Mismatch between NPOIN and number of
  points listed in mesh file"). Fixed by compacting the point
  numbering in `wsu2.f` down to only the points `ICELNOD` actually
  references, remapping element/boundary connectivity through it;
  `rsu2.f` rebuilds the identical compaction when reading the restart
  CSV back so row order lines up with the point it came from.

- **The shock-front double-point seam (colour 10) needs periodic
  pairing, not "ignore it".** EulFS's own `.petscrc`
  (`-colors -1,5,5,4,5,-1,-1,-1,-1,-1,0,-1,-1,-1`, indexed from colour
  0) assigns colour 10 `BC_TYPE_PERIODIC` (bc type `0`, see
  `EulFS.3.14/include/bnd.h`) -- it pairs the coincident
  upstream/downstream point chains rather than leaving them uncoupled.
  The first working version of this coupling used SU2's
  `MARKER_INTERNAL` (a true no-op, "parse and ignore"), which let the
  two sides of the shock drift apart over many outer iterations.
  `wsu2.f` now splits colour 10 into its two connected components (a
  small union-find over the colour's own boundary edges -- the
  up-/down-side chains are disjoint point sets even though
  geometrically coincident) and writes them as separate markers
  `bc10_1`/`bc10_2`; `su2case.cfg` pairs them with `MARKER_PERIODIC`
  and a zero (coincident, not offset) transform.

- **`bc1`/`bc2`/`bc4` should be `MARKER_FAR` (SU2's Riemann/
  characteristic far-field BC), not a fixed supersonic-inlet +
  supersonic-outlet split.** Same `.petscrc` colours array: all three
  get `BC_TYPE_FAR_FIELD` (=5) uniformly. At this case's M=20 this
  made no measurable difference (bc1 is unconditionally supersonic
  inflow and bc2/bc4 unconditionally supersonic outflow, where
  far-field's characteristic decomposition degenerates to the same
  thing) but it's the structurally correct match to EulFS and is
  cheap insurance for cases where the local flow direction isn't so
  clear-cut.

- **`CFL_NUMBER` needed to come down a lot: 0.1 -> 0.002.** Diagnosed
  by dumping `sh99.dat` every outer iteration (both here and from a
  matching EulFS reference run started from the same seed mesh) and
  comparing shock stand-off distance and downstream pressure/density
  at the stagnation point. EulFS's own downstream state stays
  essentially converged (`rho~5.92`, matching the M=20 normal-shock
  value) from iteration 1; SU2's grows without bound at CFL=0.1
  (`rho`: 5.93 -> 11.27 by iteration 50), a real numerical instability
  in the per-outer-iteration cold-restarted implicit CFD sub-solve
  near a very strong (M=20) stagnation-point shock, not a boundary or
  coupling bug -- confirmed because *more* SU2 sub-iterations per
  outer step (`ITER=20`) made it fail *sooner*, not later, and
  switching to `TIME_DISCRE_FLOW=EULER_EXPLICIT` (matching EulFS's own
  `-timestepping explicit`) also made it worse. Lowering CFL
  monotonically bought more iterations (0.1: 50, 0.02: 117, 0.005:
  396) until 0.002 cleared all 501. The actual crash at the old CFL
  values was silent (gdb: "exited normally", not a signal) --
  `fx_msh_sps.f`'s own `log/fx_msh_sps.log` (not stdout) had the real
  message, `failed matching 1st shock point of the shock n. 1`, at a
  boundary-intersection special point ("OPY" type) with a wildly
  out-of-domain coordinate, i.e. a real divergence reaching the
  fitting algorithm's own consistency checks, not a crash in SU2 or
  the converters.

## Result

**All 501 outer iterations complete.** `su2fitting_run.log` is the
UnDiFi-2D driver's own stdout (trimmed; every elided iteration prints
the same step sequence ending `su2 --> ok`). Final SU2 sub-solve:
`rms[Rho]=-1.14`, `CD=2.50` -- in the same range as `su2_capturing`'s
converged `CD=2.458`. Shock point count grew only modestly over the
whole run (60 -> 68 points, well inside `NPSHMAX=500`), consistent
with a numerically well-behaved run rather than one skating on the
edge of the next failure.

`capturing_vs_fitting.png` (`plot_compare.py`) compares the converged
shock-capturing solution against this run's final (iteration 501)
state, Mach number capped at 3 to make the post-shock structure
visible (freestream is M=20) with the sonic line (M=1) dashed. The
fitting panel's shock is visibly sharp (a true jump at the fitted
front, overlaid in white from `sh99.dat`) versus capturing's smeared
transition over several cells -- the qualitative difference the two
modes are supposed to produce. The sonic-line pocket shape is
consistent between the two. The fitted shock sits measurably farther
out than the captured one (stand-off ~1.7-1.9 vs ~1.4-1.5) -- the
residue of the CFL-driven outward drift described above, now slow
enough not to diverge over 501 iterations but not eliminated; revisit
if a closer match to EulFS's own stand-off distance matters for a
given use case (a CFL ramp -- very low for the first
O(100) iterations while the state is farthest from equilibrium, then
relaxed once established -- is the likely next lever, untried here).

## Files
- `su2case.cfg` -- the working SU2 config (`ITER=1`, `RESTART_SOL=YES`, `CFL_NUMBER=0.002`, `MARKER_PERIODIC` for the shock seam, `MARKER_FAR` for the outer boundary).
- `su2fitting_run.log` -- UnDiFi-2D's own stdout, trimmed (head + final iteration), all 501 iterations completed.
- `capturing_vs_fitting.png`, `plot_compare.py` -- side-by-side Mach-number comparison against `../su2_capturing/`'s converged solution, generated from this run's final state (`na00501.1.node/.ele`) and `sh99.dat`.
