# Solver-side assessment: species/mass-fraction transport

Issue #21 (E1), task: *"Solver-side assessment (deliverable in its own right):
the shock-capturing solver must transport a species / mass-fraction field
... This assessment may redirect which solver leads the shock-bubble work,
so it should be done early and written up."* This is that write-up.

## 1. A constraint the issue text doesn't mention: UnDiFi-2D's own `ndof`

Before asking which *solver* can carry a species field, there's a more basic
question: can UnDiFi-2D's own discontinuity-fitting data structures carry
one at all? `ndof=4` is a compile-time `parameter` in `mod_constants.f90`,
and it is not a soft default — it sizes `discontinuity_t%zu`/`%zd`
(`mod_discontinuity.f90`), every `zroesh`/`zroeshu`/`zroeshd` array passed
through `co_norm.f90`/`co_state_dps.f90`/`interp.f90`/`rd_dps.f90`, and the
Roe-variable encode/decode arithmetic (`z1=sqrt(rho)`, `z2=sqrt(rho)*H`,
`z3=sqrt(rho)*u`, `z4=sqrt(rho)*v`) baked directly into each of those files.
A 5th "species mass fraction" component riding alongside `(rho,H,u,v)`
through fitting, redistribution, and background-mesh interpolation would
touch the same files the ~200-site `GA`/`GM1` replacement (this issue's own
task 3) touches, for the same "mechanical but wide" reason. **This means
the shock-bubble exit criterion needs BOTH a solver capable of carrying a
species field AND a widened `ndof` on the fitting side — the solver
assessment below only answers half the question the issue frames.**

## 2. Per-solver capability, checked against actual source/current docs

### SU2
Not vendored in this repo (`mod_su2_solver.f90` is UnDiFi's own file-based
coupling layer — `.su2` mesh export, no SU2 source tree to inspect
directly). Checked against SU2's current public documentation instead of
relying on possibly-stale training knowledge: SU2 has a **general "species
transport" capability** — passive species-fraction transport equations
addable to both incompressible and compressible flow, no chemistry/
nonequilibrium machinery required — which is exactly the two-inert-gas
shock-bubble case (a helium bubble in air is not reacting, not
thermochemically out of equilibrium). SU2 also has a separate, much
heavier **SU2-NEMO** extension (finite-rate chemistry, thermal
nonequilibrium, a link to the Mutation++ physicochemical library) built
for hypersonic/re-entry/plasma flows — not the right tool for this
problem, but worth knowing it exists so a future contributor doesn't
reach for it by default. **Verdict: SU2's plain species-transport
capability is the closest match to what shock-bubble actually needs, with
no chemistry overhead to strip out.**

### EulFS
Fully vendored (`EulFS.3.14/`), so checked directly rather than by
inference. `EulFS.3.14/src/chemistry/` is real and substantial (~40
files) — but reading the file names and one source file's own header
comment (`massactionSU2.f90`: *"Compute source term for each chemical
species using the model implemented in SU2 ... Hoffert, Lien - Quasi-one
dimensional nonequilibrium gas dynamics of partially ionized 2-temperature
Argon"*) shows this is a **specific ionized-Argon plasma chemistry
package** (`iondegree4Ar.f`, `boltzeq4Ar.f`, `electcond.f90`,
`ohmsource.f90`, `currentflux.F` — ionization degree, Boltzmann equation,
electrical conductivity, Ohmic heating source terms), not a general
inert-gas mixture model. It would be **misleading to report "EulFS has
multi-species support"** without this qualification — the actual
infrastructure is purpose-built for MHD/plasma flows, a different
physics regime than a shock-bubble interaction with two neutral gases.
`EulFS.3.14/include/paramt.h` does confirm the state-vector framework
itself is sized generously (`MAXNOFVAR=NMAX=8`, vs. the 4 mean-flow
equations a plain 2-D Euler case uses), so there IS headroom in EulFS's
own variable-count ceiling for extra transported quantities — the gap is
a general passive-scalar transport equation and its own convective/
diffusive discretization, not a framework-level variable-count limit.
**Verdict: EulFS has real multi-species-adjacent infrastructure, but for
the wrong kind of "multi-species" (ionized plasma chemistry, not inert
mixture transport) — reusing it for shock-bubble would mean either
adapting the Argon-specific chemistry package to a generic 2-species
inert case (unclear how much of `chemistry/`'s machinery is actually
chemistry-agnostic under the Argon-specific parts) or writing a new,
simpler passive-scalar equation from scratch within EulFS's existing
NOFVAR=8 headroom.**

### NEO
Fully vendored (`NEO/`), checked directly. No species, chemistry, or
passive-scalar files anywhere in `NEO/src/` (`common.h`/
`common_variables.h`, the natural place for such a global degree-of-freedom
count, has nothing named anything like it). NEO is a compact, single-
purpose Euler solver (the project's own smallest/fastest of the three,
matching prior sessions' benchmarking against it) with no apparent
species-transport groundwork at all. **Verdict: matches the issue's own
assumption exactly — NEO needs a passive scalar added from scratch, with
no existing infrastructure to build on, or the case restricted to
single-species work when NEO is the chosen solver.**

## 3. Recommendation

**SU2 should lead the shock-bubble work**, on two independent grounds
that happen to agree: (a) it has the least-adapted-needed capability of
the three (a real, general, currently-maintained passive species-
transport feature, not a repurposed plasma-chemistry package or a
from-scratch addition), and (b) issue #17 (Phase 5a, solver-as-library)
already flagged SU2 as a `5a.3` target for its own reasons (a proper
Python/C++ API instead of file round-trips) — species-transport work and
the library-coupling work could land together rather than adding a
second, separate SU2 integration effort later. **EulFS is the second
choice**, worth a real look specifically if the shock-bubble work ever
needs its existing MHD/plasma infrastructure for an ionized-flow variant
(not the near-term inert-gas case) — its chemistry package is a real,
non-trivial asset for THAT future case, just not this one. **NEO is not
a near-term candidate** for shock-bubble specifically, given zero
existing species infrastructure, though nothing here changes its
standing for the project's other (single-species) work.

## 4. What this assessment does not answer (flagged, not attempted here)

- Whether SU2's species-transport feature interoperates cleanly with a
  discontinuity-FITTING code's moving internal boundary (an interface
  curve tracked explicitly, not a capturing scheme's smeared mixture
  fraction) — SU2's species transport was designed for capturing-mode
  use cases; whether/how the two-way UnDiFi<->SU2 coupling (`mod_su2_solver.f90`)
  would carry the species field across that boundary is unexamined.
- The `ndof` widening itself (§1) — sizing, which files need touching,
  and whether it can ride along with or must precede the `GA`/`GM1`
  replacement (issue #21 task 3) is not designed here.
- Whether EulFS's `chemistry/` package has any chemistry-agnostic
  substrate worth reusing (a generic passive-scalar transport/advection
  step underneath the Argon-specific reaction-rate source terms) or
  whether adapting it would mean rewriting the transport equation itself
  too — not read closely enough to say.
