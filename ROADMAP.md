# UnDiFi-2D Development Roadmap

**Status:** draft v1.0 — 2026-08-02
**Working branch:** `modernization-f2018`
**Scope:** short-term = language/architecture/performance modernisation;
mid-term = physics extension to shock/bubble interaction and the complete
Ben-Dor shock-reflection and shock–shock interaction nomenclature.

---

## 0. Executive summary

UnDiFi-2D is a working, published (CPC 2021), physically sophisticated
unstructured shock-fitting driver. Its weakness is not the algorithm — it is
the *software substrate*: fixed-form Fortran 77, a hand-rolled `istak/dstak`
stack allocator, compile-time array bounds, process-level coupling through
`system()` + files, zero parallelism, and a 1438-line `if/elseif` chain that
encodes the entire special-point physics.

That substrate is now the binding constraint on **both** objectives:

* it caps performance (single-threaded, O(N²) kernels, ~50 fork/exec + full
  file round-trips per time step), and
* it caps extensibility — adding one new Ben-Dor configuration today means
  editing three 1000+ line files and a hardcoded 20×20 Newton system.

The roadmap therefore front-loads an **architectural refactor to Fortran 2018**
(Part II) that is deliberately designed around the *shape of the mid-term
physics* (Part III): a polymorphic `discontinuity_t` / `special_point_t` class
hierarchy makes triple points, quadruple points, von Neumann/Guderley
reflections and shock–bubble refraction points all instances of one dispatch
mechanism rather than new branches in a chain.

Sequencing rule that governs everything below: **no refactor step may change
numerical output**. Every phase in Part II is gated by bit-comparable (or
tolerance-bounded) reproduction of the existing regression checksums in
`tests/*/checksum_*`.

---

# Part I — Assessment of the current code

## I.1 Repository anatomy

| Component | Path | LOC | Language | Role |
|---|---|---:|---|---|
| **UnDiFi-2D core** | `source/` | 12 717 | Fortran 77 fixed-form | shock-fitting driver — **the target of this roadmap** |
| Converters/utilities | `source_utils/` | 30 772 | F77 / C | `triangle2dat`, `dat2triangle`, `triangle2su2`, `su22triangle`, `NEO2triangle`, `Na_creation`, plotting |
| NEO solver | `NEO/src/` | 6 292 | C | RD shock-capturing solver (vendored) |
| EulFS solver | `EulFS.3.14/src/` | 86 904 | F77 + PETSc | RD/FV shock-capturing solver (vendored) |
| Support libs | `lib/` | 225 973 | F77/C | `libport` (stack allocator), `libmylib`, `libfxdr`, `libtirpc`, SPARSKIT2, `triangle` |
| Test cases | `tests/` | — | shell + data | 12 cases, checksum-based regression |
| Docs | `doc/userguide/` | — | Markdown → LaTeX | 7-chapter user guide |

**Core module map** (`source/`, all fixed-form `.f`):

```
main.f            2089   driver: outer time loop + all solver orchestration
fx_state_dps.f    1438   special-point jump relations (IPX/OPX/WPNR/FWP/TP/QP/TE/RR/SP/C/PC)
fx_msh_sps.f      1040   mesh repair around special points
co_pnt_dspl.f      893   shock-point displacement / hole construction
co_norm.f          817   shock normal & tangent unit vectors (upwind-biased FD)
co_utp.f           435   unsteady TRIPLE point   — 20×20 Newton, FD Jacobian
co_uqp.f           262   unsteady QUADRUPLE point — Newton, FD Jacobian
co_urr.f           190   unsteady REGULAR REFLECTION point
co_shock.f         ~230  single shock point: Rankine-Hugoniot Newton solve (4×4)
ch_sh_topology.f   281   topology change: FWP -> RR on a wedge (only case implemented)
interp.f           ~430  phantom-node interpolation (background <-> fitting mesh)
fnd_phps.f         ~420  phantom point detection
rd_dps.f / rd_dps_eq.f   shock-point redistribution
+ ~25 I/O, geometry and utility files
paramt.h / shock.com     compile-time parameters + COMMON blocks
```

## I.2 Algorithm as implemented

Per outer iteration (extracted from `main.f` step trace):

```
 co_norm  -> co_state_dps -> fx_state_dps        (jump relations, t)
 fnd_phps -> co_norm -> interp_sp -> co_pnt_dspl -> fx_msh_sps
 [unsteady: calc_vel -> solzne (ALE grid velocity to EulFS)]
 wtri     -> triangle (CDT)                      (fork/exec)
 -> {triangle2dat -> EulFS -> dat2triangle}      (fork/exec x3)
 or {triangle2su2 -> SU2   -> su22triangle}      (fork/exec x3)
 or {na2vvvv -> triangle2grd -> NEO -> NEO2triangle}
 readmesh -> zroe(1)->zroe(0) -> co_state_dps -> fx_state_dps
 [unsteady: predictor/corrector repeats the whole block]
 wsh_mean -> mv_dps -> fx_dps_loc -> fltr_dls -> interp -> mv_grid
```

The predictor–corrector path duplicates the entire block a second time inside
`main.f` — a large fraction of `main.f`'s 2089 lines is literal copy-paste of
the steady path.

## I.3 Findings

### F1 — Fixed-form F77 with F90 syntax mixed in (blocker for everything)
All of `source/*.f` is fixed-form (column 6 continuation, `!` comments in the
free part, `+`/`.` continuation markers). Yet `main.f` already uses
`character(len=20) ::`, `command_argument_count()`, `advance='no'` — the code
is *de facto* Fortran 95 written in F77 clothing. Positive: `implicit none` is
present in every file except `coords.f` and `rand.f`; only 28 `GOTO`s remain in
12.7 kLOC. **The code is much cleaner than its file extension suggests — a
mechanical conversion to free-form is low-risk and high-value.**

### F2 — `istak/dstak` PORT stack allocator
`main.f` declares `double precision dstak(9 990 000)` in a `COMMON /cstak/`
with `equivalence (dstak(1), istak(1))`, then hands out integer offsets
(`lcorg(0:2)`, `lzroe(0:2)`, `lnodcod(0:2)`, …) passed to subroutines as
`dstak(lcorg(0)+npoin(0)*ndim)`. Consequences:
* no bounds checking, no type safety, aliasing that defeats every optimiser;
* memory is a hard 9.99 M-word compile-time ceiling;
* `equivalence` of `REAL*8` and `INTEGER*4` arrays is not standard-conforming
  and blocks any move to `-fcheck=all` or to a parallel model.
This must be replaced by `allocatable` derived types before any parallelism.

### F3 — Compile-time problem-size limits
`paramt.h`: `NSHMAX=10`, `NPSHMAX=500`, `NESHMAX=499`, `NSPMAX=12`,
`NADDHOLESMAX=10`, `NprdBndMAX=1`. Every shock array is
`(ndim, npshmax, nshmax)` — i.e. **10 MB of mostly-empty static arrays** are
carried and `dcopy`-ed in full each step regardless of the actual shock count.
`NSPMAX=12` in particular is a hard ceiling that a double-Mach-reflection or a
shock–bubble simulation will hit. `paramt.h` also still carries ~90 lines of
commented-out per-test-case parameter blocks (now superseded by `input.dat`,
read in `re_inp_data.f`) — dead weight that misleads new users.

### F4 — Process-level solver coupling
50 `system()` calls in `main.f`. Each outer iteration forks `triangle`, a
format converter, the flow solver, and a converter back — with the full mesh
and solution serialised to disk in between. For the CircularCylinder case at
501 iterations this is ~2000 process spawns and ~2000 full mesh round-trips.
This is also why the code cannot be MPI-parallel today: the solver is a
*separate executable*.

### F5 — O(N_elem × N_shockedge) geometric searches
`fnd_phps.f:68-70` is a bare triple loop `ish × nshockedges × nelem`; `interp.f`
does a full `do ielem = 1, nelem` element walk per phantom node. Both are
brute-force with no spatial acceleration structure. On refined grids
(`DXCELL=0.00625` for the level-4 cylinder) these dominate the fitting-side
cost. They are also *embarrassingly parallel* — the natural first OpenMP target.

### F6 — Special-point physics as a flat `if/elseif` chain
`fx_state_dps.f` dispatches on `typespecpoints(isppnts)` — a `character*5`
literal — across 1438 lines: `IPX`,`IPY`,`OPX`,`OPY`,`WPNRX`,`WPNRY`,`FWP`,
`TP`,`QP`,`TE`,`RRX`,`RRY`,`SP`,`C`,`PC`. Each branch inlines its own state
extraction (30–60 lines of `zroeshu(4,ip1,ish1)/zroeshu(1,ip1,ish1)` index
arithmetic), then calls a bespoke solver. `fx_msh_sps.f` (1040 lines) and
`co_pnt_dspl.f` (893 lines) contain *parallel* dispatch chains on the same
type codes. **Adding one Ben-Dor configuration therefore requires coordinated
edits in three places with no compiler check that they stay consistent.**

### F7 — Hardcoded dense Newton systems with finite-difference Jacobians
`co_utp.f`: 20 unknowns, `a(11,20)`/`b(11)` constraint block assembled by ~200
lines of literal `a(1,1)=1.0d+0` statements; the Jacobian is built by
central differences (`dyn1 = |y_j|*0.001`, 2 residual evaluations per entry →
800 residual calls per Newton step). `co_shock.f` does the same for the 4×4
per-shock-point system. Costs: slow, fragile near-singular Jacobians, and no
way to express a *new* configuration without writing another such file.

### F8 — Zero parallelism
No OpenMP, MPI, coarray or OpenACC anywhere in `source/`. The available
toolchain on this machine already supports all four:
gfortran 13.3, `ifx` 2025.3, `nvfortran` 26.1 (OpenACC), OpenCoarrays `caf` 2.9.2,
`mpif90` (NVHPC). **The gap is entirely in the code, not the environment.**

### F9 — Build system
`source/Makefile` is a 1980s `mkmf` template with a stale `SRCS` list
(`co_pnt_dipl.f` — a typo for `co_pnt_dspl.f` — and `dump_sdw_info.f` missing
from `OBJS`), no dependency tracking, no `-O2`/`-march`, hardcoded
`-L../lib/linux_intel_x86_64`. `compile_all.sh` `wget`s and builds PETSc 3.14.6
from scratch on every invocation. Object files and binaries are committed into
the working tree (`source/*.o`, `bin/*`), polluting `git status`.

### F10 — Testing is checksum-on-final-output only
`tests/*/checksum_*` + `checking_fitting_runs.sh` compare end-state files.
There are **no unit tests** — not for `co_shock` (R-H solve), not for `co_norm`
(unit vectors), not for the triple-point solver. Any refactor is therefore
currently validated only at the integration level, which is exactly backwards
for the work proposed here.

### F11 — Physics gaps relative to the mid-term goal
* **Single gas, single γ.** `GA` is one scalar in `input.dat`, `GM1=GA-1`
  in a COMMON block, used directly in all state conversions. No species
  transport, no per-side γ. **Shock–bubble interaction is impossible without
  this.**
* **Open polylines only.** A discontinuity is `xysh(ndim, 1:nshockpoints, ish)`
  with both ends required to be special points. A bubble interface is a
  *closed* curve — the data structure has no cyclic mode.
* **Only two discontinuity kinds:** `typeshocks = 's' | 'd'`. No slip line
  distinct from contact, no expansion fan, no refracted-shock kind.
* **One topology transition implemented:** `ch_sh_topology.f` handles only
  FWP → RR on a wedge. RR↔MR transition, MR→DMR, refraction-pattern changes are
  all absent, so a run cannot follow a configuration through a transition.

### F12 — Documentation and provenance
The user guide is good on algorithm and workflow but has broken figure
references (`Fig. undifi_algorithm_steps`, `Fig. fig:algorithm.e`) and no API
documentation. Numerous `caldo` debug markers and large commented-out blocks
remain in the physics kernels (`fx_state_dps.f`, `co_shock.f`), obscuring which
formulation is actually live.

---

# Part II — Short term: Fortran 2018 + parallelism

**Objective:** convert `source/` into a modern, modular, parallel Fortran 2018
library + thin driver, with *no change in numerical results*, and with an
object model that makes Part III cheap.

## II.0 Guiding constraints

1. **Regression-gated.** Each phase must reproduce all 12 `tests/` checksums.
   Where bitwise reproduction is impossible (reassociation under `-O2`,
   parallel reductions), a tolerance-based comparator is introduced *first*
   (Phase 0) and the tolerance is justified and recorded.
2. **One concern per phase.** Never combine a form conversion with a semantic
   change in the same commit.
3. **The public data model is designed once (Phase 3) and not re-litigated.**
4. **Parallelism is added only after allocatable memory (Phase 2)** — `dstak`
   aliasing makes any threading unsound.

## Phase 0 — Foundations (no source changes to physics)

| # | Work item | Deliverable |
|---|---|---|
| 0.1 | `.gitignore` for `*.o`, `*.mod`, `bin/`, build dirs; purge tracked build artefacts | clean `git status` |
| 0.2 | Fix `Makefile` `SRCS`/`OBJS` (`co_pnt_dipl` typo, missing `dump_sdw_info.o`) | correct build graph |
| 0.3 | **CMake build** for `source/` + `source_utils/`, with `gfortran`/`ifx`/`nvfortran` presets, `Debug` preset carrying `-fcheck=all -ffpe-trap=invalid,zero,overflow -g` | `cmake --preset debug && cmake --build` |
| 0.4 | **Tolerance-aware regression harness**: replace byte-checksum with a field comparator (L∞/L2 on `xysh`, `zroesh`, `wsh`, nodal `zroe`) + machine-readable report | `tests/compare.py`, `ctest` targets |
| 0.5 | **Baseline capture**: run all 12 cases on all 3 solvers, archive reference fields + wall-clock + `gprof`/`perf` profile | `tests/baseline/` + profile report identifying true hotspots |
| 0.6 | **CI** (GitHub Actions): build matrix (gfortran/ifx) × (Debug/Release), run the fast subset of `tests/` | green CI on `modernization-f2018` |

*Exit criterion:* CI green; baseline archived; a deliberate 1-ulp perturbation in
`co_shock.f` is detected by the harness.

**This phase is non-negotiable and must complete before any Phase-1 commit.**
Findings F9, F10 are closed here.

## Phase 1 — Mechanical modernisation (form + syntax, semantics frozen)

| # | Work item | Notes |
|---|---|---|
| 1.1 | Fixed-form → free-form (`.f` → `.f90`), automated + reviewed | `findent`/`fprettify` then per-file diff review |
| 1.2 | Uniform lower-case keywords, 2-space indent, `end subroutine <name>` | enforced by `fprettify` in CI |
| 1.3 | `REAL*8`/`INTEGER*4` → `real(wp)` / `integer(i4)` via a `mod_kinds` module (`wp = real64`) | single point of precision control; enables a future `real128` verification build |
| 1.4 | `implicit none` in `coords.f`, `rand.f`; `implicit none (type, external)` everywhere (F2018) | closes F1 |
| 1.5 | Remove the 28 `GOTO`s → `do while` / `exit` / `cycle` / `select case` | `setbndrynodeptr` (6) and `co_shock` (6) first |
| 1.6 | Delete `caldo` debug markers and dead commented blocks; the ~90 dead lines in `paramt.h` | closes half of F12 |
| 1.7 | Replace `system()` string-building with a typed `run_external(cmd, args, log) result(status)` helper; keep `system` under the hood for now | removes 50 duplicated error-handling blocks from `main.f` (~250 lines) |
| 1.8 | Replace `call exit(n)` with `error stop` and a central `fatal(msg, code)` | F2018 conformance |

*Exit criterion:* all regressions pass at Phase-0 tolerance; `main.f90` under
1200 lines; `-std=f2018 -Wall -Wextra` clean.

## Phase 2 — Memory: kill `istak/dstak`

The single highest-leverage change. Introduce `mod_mesh`:

```fortran
module mod_mesh
  use mod_kinds, only: wp, i4
  implicit none (type, external)
  private

  type, public :: mesh_t
    integer(i4) :: npoin = 0, nelem = 0, nbfac = 0, nbpoin = 0
    integer(i4) :: nedge = 0, nhole = 0, nvt = 3
    real(wp),    allocatable :: xy(:,:)        ! (ndim, npoin)
    real(wp),    allocatable :: zroe(:,:)      ! (ndof, npoin)
    integer(i4), allocatable :: celnod(:,:)    ! (nvt, nelem)
    integer(i4), allocatable :: celcel(:,:)
    integer(i4), allocatable :: bndfac(:,:)
    integer(i4), allocatable :: nodcod(:)
    integer(i4), allocatable :: nodptr(:,:), edgptr(:,:)
    integer(i4), allocatable :: ia(:), ja(:), iclr(:)
  contains
    procedure :: alloc, free, copy_from
    final     :: mesh_finalize
  end type
end module
```

| # | Work item |
|---|---|
| 2.1 | `mod_kinds`, `mod_constants` (replace `paramt.h` parameters), `mod_config` (replace the `PARAMT` COMMON; `input.dat` reader returns a `config_t`) |
| 2.2 | `mod_mesh` as above; the `(0:2)` background/fitting/backup triple becomes `type(mesh_t) :: bkg, fit, bak` |
| 2.3 | Convert leaf routines first (`co_norm`, `interp`, `fnd_phps`, `wtri`, `readmesh`) to take `type(mesh_t), intent(in/out)` |
| 2.4 | Remove `COMMON /cstak/`, `equivalence`, `istkin/istkgt/istkrl` and the `libport` dependency |
| 2.5 | Retire `shock.com` COMMON `/SHOCKCOM/` into `mod_freestream` |
| 2.6 | Assumed-shape dummies (`real(wp), intent(in) :: xy(:,:)`) + explicit interfaces everywhere; delete all `external` declarations |

*Exit criterion:* `nm` shows no `cstak_`; `-fcheck=bounds` runs clean on all 12
cases; memory ceiling gone (verified on a 4× refined cylinder). Closes F2, F3.

## Phase 3 — The object model (designed for Part III)

This is the phase that determines whether the mid-term physics is cheap or
expensive. Design target: **every Ben-Dor configuration is a type, not a
branch.**

```fortran
module mod_discontinuity
  type, public :: discontinuity_t
    character(len=:), allocatable :: kind      ! 'shock' | 'contact' | 'slip' | 'interface'
    logical     :: closed = .false.            ! <-- enables bubble interfaces (F11)
    integer(i4) :: npoints = 0, nedges = 0
    real(wp), allocatable :: xy(:,:)           ! (ndim, npoints)  -- allocated to need
    real(wp), allocatable :: zu(:,:), zd(:,:)  ! (ndof, npoints)
    real(wp), allocatable :: nor(:,:)          ! (ndim, npoints)
    real(wp), allocatable :: w(:,:)            ! (ndim, npoints)  speed
    integer(i4), allocatable :: nodcod(:)
    integer(i4) :: gas_up = 1, gas_dn = 1      ! <-- multi-species index (F11)
  contains
    procedure :: normals          ! was co_norm
    procedure :: redistribute     ! was rd_dps / rd_dps_eq
    procedure :: displace         ! was mv_dps
    procedure :: solve_jumps      ! was co_state_dps (per-point R-H)
  end type
end module

module mod_special_point
  ! Abstract base: ONE dispatch point replacing three if/elseif chains (F6)
  type, abstract, public :: special_point_t
    integer(i4) :: id
    real(wp)    :: xy(ndim), vel(ndim)
    integer(i4), allocatable :: disc(:), endp(:)   ! incident discontinuities
    integer(i4) :: colour = 0
  contains
    procedure(nunk_i),     deferred :: n_unknowns
    procedure(residual_i), deferred :: residual      ! F(y) = 0
    procedure                       :: jacobian      ! default: FD; overridable analytic
    procedure(pack_i),     deferred :: pack_state    ! disc arrays -> y
    procedure(unpack_i),   deferred :: unpack_state  ! y -> disc arrays
    procedure                       :: solve         ! shared Newton + line search
    procedure                       :: repair_mesh   ! was the fx_msh_sps branch
    procedure                       :: check_transition => no_transition  ! -> Part III
  end type
end module
```

| # | Work item |
|---|---|
| 3.1 | `mod_discontinuity` + `mod_shock_system` (a growable `type(discontinuity_t), allocatable :: disc(:)` replacing the `NSHMAX`-sized statics) |
| 3.2 | `mod_special_point` abstract base + shared `solve` (damped Newton, line search, convergence logging, failure diagnostics) |
| 3.3 | Port existing types as concrete extensions: `sp_inlet_t` (IPX/IPY), `sp_outlet_t` (OPX/OPY), `sp_wall_t` (WPNRX/WPNRY), `sp_free_wall_t` (FWP), `sp_triple_t` (TP ← `co_utp`), `sp_quadruple_t` (QP ← `co_uqp`), `sp_regular_reflection_t` (RRX/RRY ← `co_urr`), `sp_trailing_edge_t` (TE), `sp_contact_t` (C/PC), `sp_shock_point_t` (SP) |
| 3.4 | **Registry + factory**: `type_code -> allocate(sp_XXX_t :: p)`; `class(special_point_t), allocatable :: sp(:)` array. Adding a configuration = one new module + one registry line |
| 3.5 | Replace the FD Jacobian with **analytic Jacobians** where the residual is algebraic (R-H relations are polynomial — this is the single largest arithmetic saving, F7), keeping FD as the fallback the base class provides |
| 3.6 | Collapse the duplicated predictor/corrector block in `main.f` into `subroutine advance(system, mesh, dt, stage)` |
| 3.7 | `mod_solver_iface`: abstract `flow_solver_t` with `eulfs_t`, `neo_t`, `su2_t` extensions, each owning its own convert-in/run/convert-out. Removes the three-way `if(EULFS)/if(SU2)/if(NEO)` blocks (~700 lines of `main.f`) and makes issues #4–#10 (other solvers) a *plug-in*, not a fork |

*Exit criterion:* `main.f90` < 400 lines; `fx_state_dps` chain eliminated;
regressions pass. Closes F6, F7.

**Status as of 2026-08-05**: 3.1 done as designed above (`mod_discontinuity`,
`mod_shock_system`). 3.2/3.3/3.4 were merged into one step (user's explicit
choice, to avoid landing `special_point_t` with no registry to dispatch
through it) and are done, 9 increments/commits, as `mod_special_point.f90`
(abstract base + 9 concrete types — `triple_point_t`, `quad_point_t`,
`trailing_edge_t`, `regular_reflection_t`, `wall_float_t`,
`floating_wall_point_t`, `end_point_t`, `start_point_t`, `connection_t`,
covering all 16 type codes via merged codes + discriminating fields, not
the design sketch's 10-type breakdown) and `mod_special_point_registry.f90`
(`make_special_point` factory, `sp_read`/`sp_write`). Scope ended up wider
than F6's original three chains: **seven** dispatch chains were actually
unified — `fx_state_dps.f90`, `co_pnt_dspl.f90`, `fx_msh_sps.f90`,
`fx_dps_loc.f90`, `co_norm.f90`'s correction pass, the
`re_sdw_info.f90`/`wrt_sdw_info.f90` serialization pair, and a previously-
undiscovered seventh (`interp.f90`'s `interp_sp`, a single 'SP'-only branch
missed by the initial research pass and found only once real end-to-end
regression testing actually worked — see the harness-fix commits from
2026-08-05). `main.f90` needed zero edits for any of this (every dispatch
chain's external signature was unchanged) and is therefore nowhere near the
400-line target — that depends on 3.7, tracked separately as issue #15.
Six known pre-existing latent bugs in the
original dispatch chains were found and preserved bit-for-bit rather than
fixed (documented in code comments at each site); two separate, unrelated
Newton-solve infinite-loop bugs (`co_shock`'s uninitialized seed,
`co_dc`'s unreachable tolerance) were found and fixed, since — unlike the
preserved latent bugs — they made the code non-terminating rather than
merely numerically different.

**Update, 2026-08-07: 3.5 and 3.6 are also done, closing issue #14.** 3.5
(analytic Jacobians) gave all five Newton-solve call sites
(`co_shock` x2, `co_dc`, `co_urr`, `co_uqp`, `co_utp`) hand-derived analytic
Jacobians wired through a shared `mod_newton_solve.f90`, each cross-checked
at runtime against its own FD Jacobian (`verify_jac`) over tens of thousands
of real Newton iterations before being trusted; FD remains the fallback
whenever a caller doesn't supply `jac`. 3.6 collapsed `main.f90`'s two
duplicated predictor/corrector "absorb the mesh into the shock state"
sequences into `mod_shock_advance.f90`'s `subroutine advance`. All 11
roadmap-tracked regression fixtures (10 checksum-comparable, `ShockVortex`
excluded per a pre-existing harness/case mismatch) were re-verified after
every increment; zero unexplained numeric drift throughout. Issue #14 is
closed.

**Update, 2026-08-07: 3.7 (issue #15) is also done.** `main.f90`'s
three-way `if(EULFS)/elseif(SU2)/elseif(NEO)` dispatch — actually spread
across **7** call sites, not the two obvious ~250-line blocks the initial
scoping expected — is replaced by one polymorphic `class(flow_solver_t),
allocatable :: solver` (`mod_solver_iface.f90`'s abstract base +
`mod_eulfs_solver.f90`/`mod_neo_solver.f90`/`mod_su2_solver.f90`'s
concrete types, allocated via `mod_solver_registry.f90`'s
`make_flow_solver` factory — split from the base module for the same
circular-dependency reason `mod_special_point_registry.f90` is split from
`mod_special_point.f90`). `main.f90` drops from 1647 to 1215 lines — real
progress, but **still short of the <400-line exit criterion**: the
solver dispatch was the largest single contributor, not the only one,
and further reduction is out of this issue's scope. SU2's three
pre-existing gaps (no corrector-step support, no backup/cleanup
archiving) are preserved as explicit no-ops, not fixed — SU2_CFD isn't
installed in this environment and zero regression-baseline coverage
exists for any su2_* case, so a fix would ship unverified; this matches
the project's standing policy of not opportunistically fixing latent
gaps a mechanical port happens to touch. Regression-verified via the
usual 10-fixture NEO sweep **plus**, for the first time in this
project's Phase-3 work, a EulFS-backed run (`CircularCylinder`) — every
prior increment's sweeps only ever exercised NEO. Real day-to-day
environmental drift (not caused by this change) showed up on 7 of the
11 runs; confirmed via the project's established decisive test
(rebuilding the pre-3.7 commit in an isolated worktree and re-running
two of the drifted fixtures today — the unchanged old binary reproduced
the identical drifted values). Phase 3 (3.1-3.7) is now fully complete;
issue #15 closure is a separate explicit go-ahead from the user, same
as #14.

## Phase 4 — Performance and shared-memory parallelism (OpenMP)

Order set by the Phase-0.5 profile, but the expected ranking:

| # | Work item | Expected effect |
|---|---|---|
| 4.1 | **Spatial acceleration structure** (uniform bin grid or k-d tree in `mod_search`) for `fnd_phps` and `interp` — replaces the O(N_elem·N_shockedge) scans (F5) | asymptotic; dominant on refined grids |
| 4.2 | OpenMP `do concurrent` / `!$omp parallel do` over shock points in `co_norm`, `co_state_dps`, `disc%solve_jumps` — each point's R-H solve is independent | near-linear in points |
| 4.3 | OpenMP over phantom-node search and interpolation (after 4.1) | near-linear |
| 4.4 | Remove the full-array `dcopy(ndim*npshmax*nshmax, …)` copies — with allocatables these become `move_alloc` or size-`npoints` copies | removes ~10 MB/step of pure waste (F3) |
| 4.5 | Thread-safety audit: every `open(8, file='log/*.log')` in a kernel is a shared-unit hazard → per-thread or buffered logging via `mod_log` | required for 4.2/4.3 |
| 4.6 | `do concurrent (… ) local(…) shared(…)` (F2018 locality specifiers) for the kernels that are clean, so the same loop can be offloaded in Phase 6 | one loop form, three back-ends |

*Exit criterion:* ≥ 6× on 8 cores for the fitting-side cost of the level-3
cylinder case; identical results within Phase-0 tolerance; OMP and serial
builds both in CI.

**Status as of 2026-08-07: 4.1 done** (`mod_search.f90`, a uniform bin grid
-- not a k-d-tree, since mesh cell sizes are roughly uniform within one
`DXCELL`-controlled refinement level; the module's `query_point`/
`query_segment` interface can host a k-d-tree later if a non-uniform case
needs it). Both named consumers retrofitted: `fnd_phps.f90`'s shock-
segment/cell-crossing scan and `interp.f90`'s `finder`-based phantom-
point/cell-containment scan (via a new `finder_search`, since `finder`
itself is also called directly by `mod_special_point.f90` far too rarely
to need the optimization, and touching it would only add risk). Both
build a fresh grid per call rather than sharing one across `main.f90` --
mesh regenerates every outer iteration, and redundant O(nelem)
construction is trivial next to the O(N²)-ish scans it replaces.
Regression-verified (10-fixture sweep, candidate values checked
bit-identical to already-characterized day-to-day environmental drift,
not new drift). No timing/speedup measurement performed — the exit
criterion's "level-3 cylinder case" isn't checked into the repo (only
level-0, `DXCELL=0.10`, exists; level-3 would be `DXCELL=0.0125` by
inference from the documented level-4 value) and generating it plus
scaling-benchmark tooling is separate follow-up work.

**Scoping research also found two things this table's text didn't
anticipate**, left for whoever picks up 4.2:
- **`disc%solve_jumps` doesn't exist.** It's aspirational text from
  Phase 3's original design sketch; `discontinuity_t` shipped with zero
  type-bound procedures (Phase 3 built the special-point side instead).
  The real 4.2 target is the still-standalone `co_state_dps.f90` (calls
  `co_shock`/`co_dc` per point, both already reentrant via Phase 3.5's
  per-call Newton contexts).
- **Two real thread-safety hazards sit directly in 4.2's named files**,
  beyond the `open(8,...)` issue 4.5 already anticipates: `co_state_dps.f90`
  writes every loop iteration to `mod_freestream`'s module-global scratch
  (`z1m/z1v/.../z4m/z4v`) -- harmless serially (write-only, silently
  overwritten), a real data race under `!$omp parallel do`. Phase 3
  already fixed the identical pattern in `mod_special_point.f90`'s
  `tp_solve_state` by localizing the variables; `co_state_dps.f90` itself
  was never converted and needs the same fix before 4.2 lands. Separately,
  `co_norm.f90`'s sign-flip pass is a count-then-bulk-flip reduction
  (needs `omp reduction` + a barrier), not a naively parallel loop like
  its main per-point normal-computation loop (which has no such hazard).

**4.2 done**: `!$omp parallel do` over `co_norm.f90`'s per-point normal
loop (+ `reduction(+:ii)` on its count-then-flip sign pass) and
`co_state_dps.f90`'s per-point R-H-solve loop, after fixing the
`mod_freestream` scratch-variable data race both hazards above flagged.
CMake: OpenMP is opt-in (`UNDIFI_ENABLE_OPENMP`, OFF by default, new
`gfortran-openmp` preset) -- every other preset still builds the `!$omp`
lines as inert comments. Verified bit-identical output at
`OMP_NUM_THREADS=1` vs `=8` on `CircularCylinder neo/fitting/steady`. No
speedup measured yet -- expected, since these loops are a small slice of
this level-0 case's runtime; the exit criterion needs the still-missing
level-3 mesh.

**4.3 done**: `!$omp parallel do` over `interp.f90`'s two per-point
phantom-interpolation loops, `fnd_phps.f90`'s candidate-triangle geometry
scan (with a named critical section around the shared `nodcod(n1/n2/n3)`
update -- a node can be touched by candidate lists from different shock
segments), and `co_pnt_dspl.f90`'s shock-node displacement loops (no
shared state, parallelize directly).

**4.4 done**: narrowed `main.f90`'s three UNSTEADY predictor/corrector
`dcopy` blocks from the compile-time-max element count down to real
point/shock counts (~10 MB/step of pure waste removed); also fixed a
real pre-existing bug found along the way (two of the three blocks
under-copied `zroeshdold`/`zroeshuold` using the wrong per-element
count, silently leaving shock slots 6-10 stale for any UNSTEADY case
with >5 shocks -- doesn't fire on any current fixture, all have ≤5).

**4.5 partial**: `mod_log.f90` exists (`log_line`, one named
`!$omp critical` section) and every in-loop write inside the loops 4.2
and 4.3 actually parallelize is converted. The other ~24 files'
`open(8, file='log/*.log')` calls are untouched -- a full rollout is
still open work. `mod_freestream.f90`'s now-dead
`z1m/z1v/.../z4m/z4v` declarations (nothing `use`s them since the 4.2
hazard fix) are also still flagged for cleanup here.

**Not started**: 4.6 (`do concurrent` conversion), plus a CI job
building+running the OMP variant, a thread-count-sweep benchmark
script, and generating the level-3 cylinder mesh the exit criterion
needs. Issue #16 stays open until those land.

## Phase 5 — Distributed memory (MPI / coarrays)

Two distinct parallelisations, do not conflate them:

**5a — Solver-side (the real cost).** The flow solver already has an MPI build
(`EulFS.3.14/src/mpi/`, PETSc). The blocker is F4: coupling is `system()` +
files.
| 5a.1 | `mpirun` the solver from the driver, driver on rank 0 (cheapest, ships first) |
| 5a.2 | Replace file round-trips with an in-memory handoff: link EulFS/NEO as a **library** with an `init/step/finalize` API rather than an executable. This is the structural fix for F4 and removes ~2000 process spawns per run |
| 5a.3 | For SU2, use its Python/C++ API rather than the `.su2` file round-trip |

**5b — Fitting-side.** The shock front is 1-D in a 2-D domain, so the fitting
work is O(N_shockpoints) — distributing it is only worthwhile for many-shock
configurations (which is exactly what Part III produces).
| 5b.1 | Domain decomposition by *discontinuity*: each rank owns a subset of `disc(:)`; halo = the special points shared between them |
| 5b.2 | **Coarray implementation** for the special-point coupling: `real(wp) :: sp_state(:)[*]` with `sync images` at the jump-relation stage — natural fit, since a special point couples exactly the 2–5 discontinuities that meet there |
| 5b.3 | MPI alternative behind the same `mod_parallel` interface for portability (OpenCoarrays and NVHPC are both present; ifx coarrays too) |
| 5b.4 | Parallel I/O: one `wtri`/`readmesh` per rank is wrong — use rank-0 gather initially, MPI-IO later |

*Exit criterion:* a 4-shock case runs correctly on 4 images/ranks with results
matching the serial run; the library-coupled solver path (5a.2) demonstrated on
at least one solver.

## Phase 6 — GPU offload (OpenACC)

Lowest priority — the fitting kernels are small and the solver is the cost
centre. Do this only after 5a.2.

| 6.1 | `!$acc` on the Phase-4.6 `do concurrent` kernels (`nvfortran -acc`) |
| 6.2 | `!$acc data` regions around the outer loop so shock arrays stay resident |
| 6.3 | Evaluate: batched R-H Newton solves across all shock points as one GPU kernel — the natural GPU shape for this algorithm |
| 6.4 | Benchmark honestly against the OpenMP path and **document if it loses** |

*Exit criterion:* a measured, published comparison. A negative result is an
acceptable and useful deliverable here.

## Phase 7 — Compactness and flexibility (continuous)

| 7.1 | `mod_io` with a single named-unit registry; retire the bare `open(8, …)`/`open(12, …)` pattern (~30 occurrences) |
| 7.2 | Structured input: keep `input.dat` readable but add a **namelist** or TOML front-end with defaults and validation; `input.dat` becomes one supported dialect |
| 7.3 | HDF5/XDMF output alongside the Tecplot/`.dat` writers → direct ParaView, no `dat2paraview` step |
| 7.4 | FORD/Doxygen API docs generated in CI; fix the broken figure refs in `doc/userguide` (F12) |
| 7.5 | Unit tests (`test-drive` or `pFUnit`) for: `co_shock` R-H against analytic oblique-shock relations; `co_norm` against analytic curves; each `special_point_t%residual` against a manufactured solution (F10) |

---

# Part III — Mid term: shock/bubble and the full Ben-Dor nomenclature

**Objective:** UnDiFi-2D should be able to fit, and follow through transitions,
the complete catalogue of 2-D shock-wave reflection and interaction
configurations in Ben-Dor, *Shock Wave Reflection Phenomena* (2nd ed., 2007),
plus shock–bubble (shock–inhomogeneity) interaction.

**Dependency:** Part III presumes Phase 3. Each new configuration below is then
a new `special_point_t` extension + a registry entry + a verification case —
typically 200–400 lines, not a three-file surgery.

## III.1 Capability gap analysis

| Ben-Dor configuration | Abbrev. | Status today | Gap |
|---|---|---|---|
| Regular reflection | RR | ✅ `co_urr` + `RRX/RRY`, tests `RegularReflection-1/2` | port to `sp_regular_reflection_t` |
| Single Mach reflection | SMR | ✅ `co_utp` + `TP`, tests `MachReflection-1/2` | port to `sp_triple_t` |
| Direct Mach reflection | DiMR | ⚠️ implicitly covered by TP | needs explicit classification + test |
| Stationary Mach reflection | StMR | ❌ | triple point with zero trajectory angle; degenerate case of `sp_triple_t` |
| Inverse Mach reflection | InMR | ❌ | triple point moving *toward* the surface; requires the InMR→TRR transition |
| Transitioned regular reflection | TRR | ❌ | topology change MR→RR |
| Transitional Mach reflection | TMR | ❌ | triple point + reflected-shock *band* (kink); new type |
| Double Mach reflection | DMR (DMR⁺/DMR⁻) | ❌ | **two** interacting triple points + second slipstream; `NSPMAX=12` and the topology machinery both block it |
| Terminal double Mach reflection | TerDMR | ❌ | extension of DMR |
| von Neumann reflection | vNR | ❌ | weak-shock domain; 3-shock theory has no solution → needs a distinct closure |
| Vasilev reflection | VR | ❌ | von Neumann paradox regime |
| Guderley reflection | GR | ❌ | supersonic patch + expansion fan behind the triple point; needs a *fan* discontinuity kind |
| Shock–shock, same family (overtaking) | — | ⚠️ tests `SSInteractions1-2`, `2-1`, `2-2` exist | classify; verify against analytic |
| Shock–shock, opposite family (crossing) | — | ✅ `co_uqp` + `QP` | port to `sp_quadruple_t` |
| Edney Type I–VI (shock-on-bow-shock) | — | ⚠️ Type IV referenced in `paramt.h` comments and `CircularCylinder` tuning | no explicit Edney classification or test matrix |
| Shock–contact refraction (RRR/IRR, slow-fast / fast-slow) | — | ❌ | needs per-side γ and a refraction special point |
| Shock–expansion interaction | — | ❌ | needs an expansion-fan discontinuity kind |
| **Shock–bubble interaction** | SBI | ❌ | needs closed interfaces + multi-species (below) |

## III.2 Enabling work (prerequisites for the table above)

### E1 — Multi-species / variable-γ gas model
Blocks SBI and all refraction problems (F11).
* `mod_gas`: `type gas_t` with γ, R, and (later) a caloric model; a
  `gas_registry` so each side of each discontinuity names its gas.
* `discontinuity_t%gas_up / %gas_dn` (already in the Phase-3 sketch).
* Replace every direct `GA`/`GM1` COMMON reference (~200 sites) with
  `gas(i)%gamma`.
* Solver side: the shock-capturing solver must transport a species/mass-fraction
  field. **Assess per solver** — SU2 has multi-species; EulFS and NEO need
  either an added passive scalar or a case restriction. This assessment is
  itself a roadmap deliverable and may redirect which solver leads SBI work.
* R-H with different γ on the two sides (contact/interface) → generalise
  `co_shock`'s residual.

### E2 — Closed discontinuity curves
Blocks SBI (F11).
* `discontinuity_t%closed = .true.`: cyclic indexing in `normals`,
  `redistribute`, `displace` — index arithmetic modulo `npoints`, no endpoint
  special points.
* `wtri`/`fnd_phps`/`co_pnt_dspl`: a closed curve digs an *annular* hole with
  an interior sub-domain — the hole/region marking logic in `wtri.f` must
  handle a nested region, and `triangle` needs the right `-A` region attributes.
* Interface-specific redistribution: a bubble stretches and folds; point
  density control must be curvature-based, not just `DXCELL`-based.

### E3 — Discontinuity-kind extension
* `'slip'` (distinct from `'contact'`: zero pressure jump, tangential velocity
  jump — needed for DMR's second slipstream and for the triple-point slipstream
  to be treated correctly).
* `'fan'` (centred expansion): needed for GR and shock–expansion interaction.
  This is the largest single physics addition — a fan is not a curve with two
  states but a region.
* `'interface'` (material interface, γ jump): needed for SBI and refraction.

### E4 — Automatic configuration detection and transition
Generalises `ch_sh_topology.f` (currently FWP→RR on a wedge only).
* `mod_classify`: given local incident-shock strength, wedge/deflection angle
  and freestream Mach, evaluate the **detachment criterion** and the
  **von Neumann (mechanical-equilibrium) criterion**; report the configuration
  and whether the state is in the **dual-solution domain**.
* `special_point_t%check_transition` (deferred hook from Phase 3) returns a
  requested new type; a `mod_topology` transaction applies it: allocate/free
  discontinuities, re-map `shinspps`, rebuild the local mesh.
* Transitions to implement, in order: RR↔MR, MR→TRR, SMR→TMR→DMR,
  MR→InMR→TRR.
* **Hysteresis study** (RR↔MR in the dual-solution domain) is a natural
  validation *and* publication target once this exists.

### E5 — Robustness for many interacting discontinuities
* Remove `NSPMAX=12`, `NSHMAX=10` (Phase 2/3 do this) — DMR alone needs ~2
  triple points, 5 discontinuities and several special points; SBI generates
  many more.
* The shared Newton in `special_point_t%solve` needs damping, line search,
  bounds enforcement (ρ, p > 0) and a diagnostic dump on failure — today a
  failed `co_utp` silently returns `ifail`.
* Continuation/homotopy for the ill-conditioned weak-shock configurations
  (vNR/VR/GR), where the 3-shock system is near-singular by construction.

## III.3 Physics work packages

**WP-A — Consolidate what exists (unblocks everything, no new physics)**
Port RR, SMR, QP, TE, wall/inlet/outlet points to `special_point_t`; add
analytic-solution unit tests (oblique shock relations, 2-shock and 3-shock
theory); build an explicit Edney Type I–VI test matrix on the existing
`CircularCylinder` geometry.

**WP-B — Classification and transition (E4)**
`mod_classify` + `check_transition` + RR↔MR. Validation: Ivanov's wedge cases
(already partly present as commented parameter sets in `paramt.h`), and the
dual-solution/hysteresis loop.

**WP-C — Double Mach reflection (E3 slip + E5)**
Two coupled triple points, second slipstream, DMR⁺/DMR⁻/TerDMR classification.
Validation: standard pseudo-steady wedge DMR cases with published triple-point
trajectory angles.

**WP-D — Weak-shock domain: vNR, VR, GR (E3 fan + E5 continuation)**
The hardest package. GR requires the supersonic patch and an embedded
expansion fan. **This is a genuine research contribution if fitted rather than
captured** — the von Neumann paradox is precisely where capturing schemes are
ambiguous, so a fitting code that resolves it cleanly is a strong result.

**WP-E — Shock–bubble interaction (E1 + E2)**
Planar shock impinging on a cylindrical gas inhomogeneity (the canonical
Haas & Sturtevant helium/R22 configurations). Requires closed interface +
two gases + refraction special points where the shock meets the interface, and
the interface's later Richtmyer–Meshkov roll-up. Staged:
1. shock–planar-interface refraction (regular refraction, one special point);
2. shock–circular-interface, early time (two moving refraction points);
3. full SBI including transition to irregular refraction and interface roll-up
   (the point at which fitting the interface may need to hand off to capturing —
   **an honest scope boundary to establish, not to paper over**).

**WP-F — Bubble/shock beyond the canonical case**
Multiple bubbles, bubble–wall, bubble collapse. Deliberately left open pending
WP-E outcomes.

## III.4 Verification strategy for Part III

Every configuration lands with: (i) an analytic or reference-data comparison
(2-/3-shock theory, published triple-point trajectory angles, Haas–Sturtevant
interface positions); (ii) a grid-convergence study demonstrating that fitting
preserves the design order across the discontinuity; (iii) a comparison against
the same case run in capturing mode with the same solver. **(iii) is the
argument for the whole project and should be in every paper.**

---

# Part IV — Milestones and issue map

| Milestone | Content | Gate |
|---|---|---|
| **M0** Foundations | Phase 0 | CI green, baseline archived, perturbation detected |
| **M1** Modern Fortran | Phases 1–2 | free-form F2018, no `dstak`, bounds-clean, regressions pass |
| **M2** Object model | Phase 3 | `main.f90` < 400 lines, one dispatch point, regressions pass |
| **M3** Shared-memory parallel | Phase 4 | ≥6× on 8 cores, results unchanged |
| **M4** Distributed / library coupling | Phase 5 | solver-as-library demonstrated; multi-image fitting correct |
| **M5** Consolidation + classification | WP-A, WP-B | Edney matrix + RR↔MR transition + hysteresis |
| **M6** Complex reflections | WP-C, WP-D | DMR validated; vNR/VR/GR attempted with honest reporting |
| **M7** Shock–bubble | WP-E | Haas–Sturtevant configuration reproduced |
| **M8** GPU (optional) | Phase 6 | measured comparison, positive or negative |

Ordering note: M3 and M5 are independent and can proceed in parallel once M2
lands. M7 depends on M2 (object model) and E1/E2, **not** on M3/M4.

## Issue map

| Issue | Milestone | Title |
|---|---|---|
| [#11](../../issues/11) | M0 | Foundations: CMake, tolerance-aware regression harness, CI, profiling baseline |
| [#12](../../issues/12) | M1 | Phase 1 — fixed-form F77 → free-form Fortran 2018 |
| [#13](../../issues/13) | M1 | Phase 2 — replace `istak/dstak` and compile-time size limits |
| [#14](../../issues/14) | M2 | Phase 3 — `discontinuity_t` / `special_point_t` object model |
| [#15](../../issues/15) | M2 | Phase 3.7 — abstract `flow_solver_t` interface (unblocks #4–#10) |
| [#16](../../issues/16) | M3 | Phase 4 — OpenMP + spatial acceleration structures |
| [#17](../../issues/17) | M4 | Phase 5a — solver-as-library coupling |
| [#18](../../issues/18) | M4 | Phase 5b — coarray / MPI distributed fitting |
| [#19](../../issues/19) | M8 | Phase 6 — OpenACC GPU offload |
| [#20](../../issues/20) | — | Phase 7.5 — unit test suite for the physics kernels |
| [#21](../../issues/21) | E1 | Multi-species / variable-γ gas model |
| [#22](../../issues/22) | E2 | Closed discontinuity curves |
| [#23](../../issues/23) | E3 | Slip lines, expansion fans, material interfaces |
| [#24](../../issues/24) | E4 / M5 | Configuration classification and topology transition |
| [#25](../../issues/25) | M5 | WP-A — consolidate existing configurations, Edney I–VI matrix |
| [#26](../../issues/26) | M6 | WP-C — double Mach reflection |
| [#27](../../issues/27) | M6 | WP-D — von Neumann / Vasilev / Guderley reflection |
| [#28](../../issues/28) | M7 | WP-E — shock–bubble interaction |
| [#29](../../issues/29) | — | WP-F — multiple bubbles, bubble–wall, collapse |
| [#30](../../issues/30) | — | Phase 7 — I/O, input format, HDF5/XDMF, documentation |

Pre-existing solver-coupling issues #3–#10 are subsumed by #15: once
`flow_solver_t` exists, each becomes one module implementing three deferred
procedures.

## Risk register

| Risk | Impact | Mitigation |
|---|---|---|
| Refactor silently changes results | High | Phase 0 harness *before* any refactor; one concern per commit |
| Vendored solvers (EulFS/NEO) diverge from upstream | Medium | Keep them untouched in M0–M3; treat solver-as-library (5a.2) as a separate, revertible track |
| Fan/GR physics (WP-D) proves intractable in a fitting framework | Medium | Time-boxed; a documented negative result is a valid deliverable |
| SBI needs a multi-species capturing solver none of the three provide | High | E1 assessment is an early, explicit deliverable that can redirect the solver choice |
| Interface roll-up (RMI) exceeds what a fitted polyline can represent | High | Establish the scope boundary in WP-E stage 2; consider hybrid fitting/capturing hand-off |
| Coarray portability (OpenCoarrays vs ifx vs NVHPC) | Low | `mod_parallel` abstraction with an MPI back-end |

## Non-goals (explicitly out of scope for this roadmap)

* 3-D extension (a separate F90 3-D code already exists per `doc/userguide`).
* Rewriting or replacing EulFS/NEO themselves.
* Viscous/turbulent shock interaction (shock–boundary-layer) — a natural
  successor, not part of this plan.
* Real-gas / reacting flow beyond the variable-γ needed for E1.

---

## Appendix A — Target directory layout after M2

```
src/
  core/      mod_kinds, mod_constants, mod_config, mod_log, mod_io, mod_error
  mesh/      mod_mesh, mod_search (bins/kd-tree), mod_triangulate (triangle iface)
  physics/   mod_gas, mod_thermo, mod_rankine_hugoniot, mod_characteristics
  disc/      mod_discontinuity, mod_shock_system, mod_normals, mod_redistribute
  points/    mod_special_point (abstract)
             sp_inlet, sp_outlet, sp_wall, sp_free_wall,
             sp_triple, sp_quadruple, sp_regular_reflection, sp_contact,
             sp_double_mach, sp_von_neumann, sp_guderley, sp_refraction   <- Part III
  solvers/   mod_solver_iface, solver_eulfs, solver_neo, solver_su2, solver_...
  parallel/  mod_parallel (OpenMP / coarray / MPI back-ends)
  driver/    undifi.f90   (< 400 lines)
test/
  unit/      shock relations, normals, per-special-point residuals
  regression/ the 12 existing cases + the Part III verification matrix
```

## Appendix B — Files ranked by refactor priority

| Rank | File | LOC | Why |
|---:|---|---:|---|
| 1 | `main.f` | 2089 | duplication, 50 `system()` calls, all orchestration — biggest reduction |
| 2 | `fx_state_dps.f` | 1438 | the dispatch chain; blocks all Part III |
| 3 | `fx_msh_sps.f` | 1040 | parallel dispatch chain on the same type codes |
| 4 | `co_pnt_dspl.f` | 893 | third parallel dispatch chain |
| 5 | `co_norm.f` | 817 | hot, parallelisable, needs closed-curve support (E2) |
| 6 | `co_utp.f` | 435 | template for every future Ben-Dor solver — get its abstraction right |
| 7 | `interp.f` / `fnd_phps.f` | ~850 | the O(N²) hotspots (F5) |
| 8 | `paramt.h` / `shock.com` | 137 | COMMON + compile-time limits — removed in Phase 2 |
