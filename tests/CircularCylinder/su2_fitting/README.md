# CircularCylinder-1, SU2 shock-fitting smoke test

First run of SU2 through `source/main.f`'s outer UNDIFI loop itself
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
  references, remapping element/boundary connectivity through it.
- **`rsu2.f`** rebuilds the identical compaction (from the same
  `ICELNOD`) when reading SU2's restart CSV back, so row order lines
  up with the point it was written from; points the compaction
  dropped (the phantom padding) simply keep whatever `ZROE` `RTRI`
  already loaded them with.
- **`su2case.cfg` needs `MARKER_INTERNAL= ( bc5, bc6, bc7, bc8, bc9,
  bc10 )`** -- fitting mode re-triangulates every iteration with the
  fitted shock front as an internal double-point seam, coloured 5-10
  per `type.dat`'s colour>4 "-1" rows. Without `MARKER_INTERNAL`, SU2
  refuses to parse a mesh marker the config doesn't define.

## Result

Ran cleanly for **42 outer iterations** at `CFL_NUMBER=0.1` (`22`
iterations at `CFL_NUMBER=0.5` before the same failure below hit
sooner), with physically reasonable output each step -- e.g. iteration
42: `rms[Rho]=-0.62`, `CD=2.68`, in the same range as
`su2_capturing`'s converged `CD=2.458` (`su2fitting_run.log`).

**Not yet resolved**: by iteration ~43 the run fails, either in
UnDiFi's own point-location search (`source/interp.f`, "search failed
for vertex coords") or in `triangle2su2`'s boundary-face matching
(`source_utils/su22triangle/rtri.f`-style mismatch), depending on the
CFL used. The reference EulFS fitting run on this exact test case
completes all 501 iterations without either failure
(`../checksum_fitting_eulfs`), so this looks like the fitted shock
front drifting further per SU2 step than UnDiFi's mesh-update
kinematics expect -- a robustness gap between SU2's single-step state
and what EulFS's specifically-tuned stepping produces, not a
converter/wiring bug. Left as a follow-up: likely needs either a
gentler SU2 stepping strategy (sub-iterations, CFL ramp) or looking at
whether the shock-velocity estimate `co_norm`/`fx_msh_sps` derive from
the CFD state is more sensitive to inter-step state jumps than
intended.

## Files
- `su2case.cfg` -- the working SU2 config (`ITER=1`, `RESTART_SOL=YES`, `MARKER_INTERNAL` for the shock seam).
- `su2fitting_run.log` -- UnDiFi-2D's own stdout through iteration 42, ending at the iteration-43 failure.
