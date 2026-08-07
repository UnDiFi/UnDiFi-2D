# CircularCylinder-1, SU2 shock-capturing smoke test

First real run of the `triangle2su2`/`su22triangle` converters (issue
#3) against a full test case, in shock-capturing mode (no fitted front
-- SU2 captures the bow shock on its own, the same role EulFS/NEO
capturing-mode runs play per `doc/userguide/060_tutorials.md`'s
"Creation of shock file" step). `main.f` is not wired up to call SU2
yet, so this was run by hand: `triangle2su2` on `na00.1.*`, `SU2_CFD`
built from `codes/SU2` (meson+ninja, MPI/CGNS/TecIO/autodiff/CoolProp/
Mutation++/MLPCpp all disabled), `su22triangle` to convert back.

Case parameters (`circularcylinder.cfg`) reproduce the documented
freestream exactly: M_1=20, gamma=1.4, derived from
`tests/CircularCylinder/na00.1.node`'s own Roe-vector state via
Z=(sqrt(rho), sqrt(rho)*H, sqrt(rho)*u, sqrt(rho)*v) --
p=T=0.0017857142857142857, a=0.05, M=20. Boundary markers per
`tests/CircularCylinder/type.dat` (1=inflow, 2/4=outflow, 3=wall).

**Converged**: 82 iterations, `rms[Rho]=-8.04`, `CD` stable at 2.458
(`su2_run.log`). Sanity check: peak density ratio across the bow shock
in the solution is 5.55, vs. the theoretical M=20 normal-shock value
of 5.93 (`(gamma+1)*M^2 / ((gamma-1)*M^2+2)`) -- expected to be a bit
lower on this ~351-node mesh since the shock is smeared over a few
cells and the peak may not land exactly on the stagnation streamline
node. A real bow shock was captured, not a numerical artifact.

**Non-obvious fix along the way**: `FLUID_MODEL` must be set to
`IDEAL_GAS` explicitly. SU2 defaults to `STANDARD_AIR`, which silently
*ignores* `GAS_CONSTANT`/`GAMMA_VALUE` and uses real air's
R=287.058 -- this desynced the restart file's nondimensional units
(R=1) from the freestream/marker state and produced a run that
"converged" to a physically inconsistent mixed-unit solution before
this was caught and fixed.

**Also non-obvious**: at M=20, `MUSCL_FLOW=YES` (2nd-order) starting
from a uniform freestream everywhere never converged (parked around
`rms[Rho]=-1.6..-1.7` with oscillating CD and non-physical-state
warnings) -- switched to first-order (`MUSCL_FLOW=NO`) with a lower
CFL (0.1, adapted up to 1.5) for robustness, standard practice for a
strong shock forming from scratch. Revisit 2nd order once this is
wired into `main.f`'s per-iteration loop, where each step starts from
an already-close solution rather than uniform freestream.

## Files
- `circularcylinder.cfg` -- the working SU2 config.
- `na00.1.su2` -- converted mesh (from `na00.1.node/.ele/.neigh/.edge/.poly`).
- `na00.1_restart_converged.csv` -- SU2's converged solution (conservative variables).
- `na00.1.node.su2_capturing` -- same solution converted back to Roe-vector `.node` format via `su22triangle`.
- `su2_run.log` -- full solver log.
