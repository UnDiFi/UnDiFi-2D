# Solver-coupling landscape: what's in `codes/` and what it's for

Requested alongside issue #21 (E1)'s solver-side species-transport
assessment, but broader in scope: an inventory of every codebase vendored
as a sibling of `UnDiFi-2D` in `codes/`, assessed for shock-fitting
coupling suitability — not only the shock-bubble task, but possible
extensions generally. Seven of these thirteen directories already have a
thorough, existing GitHub-issue-level assessment (#2-#10); this note
summarizes those and does the same depth of work for the six that don't:
**ADPDIS3D, HTR, MUM, ramses, amrvac, UnDiFi-3D**.

## 1. Already assessed (issues #2-#10) — summary only, not re-litigated

| Solver | Issue | Mesh type | Verdict |
|---|---|---|---|
| SU2 | #3 | unstructured, native 2D | Easiest of the eight reviewed there — `.su2` format closely matches Triangle's own output |
| code_saturne | #4 | unstructured, 3D-internally (needs 2D-as-3D extrusion) | Second easiest; compressible path is secondary to its incompressible core, spike recommended first |
| OpenFOAM | #5 | unstructured, 3D-internally (needs extrusion) | Similar tier to code_saturne; more converter surface area (full case-directory tree, not one mesh file) |
| Nektar++ | #6 | unstructured, native 2D, high-order spectral/hp | Topologically matches Triangle output but needs a linear-simplex-to-high-order-expansion converter |
| COOLFluiD | #7 | unstructured | Architecturally the best conceptual fit (partitioned advance/exchange/resume is a first-class feature) but heaviest to stand up — large C++/CMake plugin framework, sparse external docs |
| CFDWARP | #8 | structured multiblock | **Blocked on #2** (no body-fitted mesh capability at all); test cases compiled into the binary |
| MFC | #9 | Cartesian only | **Blocked on #2**; otherwise strong (GPU-capable, Python-scriptable) once the representation question is settled |
| FLOW3D_MS | #10 | structured | **Blocked on #2**; hardest of all eight — no general solver entry point, per-case compiled executables, thin documentation |

Issue #2 itself is the load-bearing prerequisite for the three structured
candidates: it asks for a real design decision (embedded/immersed-
boundary-style representation on a background Cartesian/structured grid,
informed by Assonitis et al. 2022's extrapolated shock-fitting work)
before any of CFDWARP/MFC/FLOW3D_MS can be scoped concretely. Nothing in
this note changes that — it's still the right next step before touching
any structured candidate, including the two new structured-ish ones found
below (ADPDIS3D).

## 2. Not yet assessed

### UnDiFi-3D — not a coupling target, a sibling implementation

The most important finding in this note. `UnDiFi-3D` (formerly
`2Dfit-3DEulFS`) is not an external flow solver to couple with — it is
**the same shock-fitting algorithm this project implements, in 3D**,
with a startlingly direct structural correspondence:

- Same 5-step S-F cycle (remesh shock surface -> re-tetrahedralize cavity
  around it -> run the flow solver one iteration -> solve R-H via Newton
  for the new shock state -> move the shock points -> loop), matching
  UnDiFi-2D's own outer loop step-for-step (Triangle -> flow solver ->
  `co_shock`/`co_state_dps`-equivalent -> `mv_dps`-equivalent).
- Same kernel names and structure: `co_shock.f90` there implements
  **the identical three solvers this session worked on directly** --
  `co_shock` (4x4 R-H + one downstream Riemann invariant, exactly
  UnDiFi-2D's own `co_shock.f90`), `co_dc` (7x7 contact discontinuity,
  **"implemented but never called from main.f90"** in UnDiFi-3D --
  worth checking whether UnDiFi-2D's own `co_dc` is fully wired or has
  the same gap), and `co_oblkshock` (an 8x8 oblique-shock variant with
  the normal itself as an unknown, which UnDiFi-2D has no equivalent of
  at all -- a genuinely different closure UnDiFi-2D might want to look
  at if `co_norm`'s normal-then-solve split ever needs revisiting).
- Same Roe parameter vector convention (`Z = sqrt(rho)*(1,H,u,v,w)`,
  i.e. UnDiFi-2D's own `zroesh` encoding with a third velocity
  component).
- Same "drive the flow solver as an external black box via `system()`"
  architecture, and **already supports EulFS or SU2** interchangeably --
  meaning there is likely real, transferable experience/code for the
  EulFS<->SU2 swap that issue #3 (SU2 coupling for UnDiFi-2D) could
  learn from directly, not just by analogy.
- Its own `ROADMAP.md` is written in the same style as this project's
  (draft status header, per-item GitHub issue prefixes, a structured
  "what's wrong / prioritized plan / where this can go" breakdown) --
  strongly suggesting a comparable modernization effort is either
  planned or already underway there too.

**Recommendation**: before implementing anything new for E1/E2/E3 that
UnDiFi-3D's `co_shock.f90` already has a 3D analogue of (particularly
`co_oblkshock`'s different closure, and whatever explains `co_dc` being
unwired there), read that file directly for a second, independent
implementation of the same physics to cross-check against -- cheaper
than re-deriving from scratch, and a mismatch between the two would
itself be worth understanding. More broadly, `UnDiFi-3D` is the natural
long-term answer to "possible extensions" beyond the Ben-Dor 2D
catalogue this roadmap's Part III targets -- a 3D generalization already
has a working starting point in this same `codes/` tree, not a blank
page.

### HTR (Prometeo) — modern GPU multi-species reacting-flow solver

A Legion/Regent task-based-runtime solver (`src/prometeo.rg` and
friends, `.rg`/`.cc`/`.cu` per module -- CPU and CUDA implementations
side by side throughout). Directly relevant to **two** separate open
items, not just E1:

- **E1 (multi-species)**: `prometeo_mixture.rg`/`.cc`/`.cu` and
  `prometeo_chem.rg`/`.cc`/`.cu` are real, current, GPU-native mixture
  and chemistry modules -- built around species transport from the
  ground up, unlike SU2 (species transport is an addable feature) or
  EulFS (chemistry exists but is Argon-plasma-specific, per this
  session's own `species-transport-solver-assessment.md`). If shock-
  bubble ever needs GENUINE reacting-flow capability (not just an inert
  two-gas mixture), HTR is a stronger candidate than either SU2 or EulFS
  on paper -- untested here, no build/run attempted, but the module
  names and per-task CPU/CUDA structure both point to a solver
  seriously built for multi-species flows.
- **Phase 6 (issue #19, GPU offload)**: HTR's CUDA implementations
  throughout (`prometeo_cfl.cu`, `prometeo_chem.cu`, `prometeo_bc.cu`,
  `prometeo_mixture.cu`, ...) mean it is already a working example of
  exactly the OpenACC/GPU question #19 asks about, for a Legion-based
  rather than a directive-based approach -- worth a look when #19 is
  scoped, purely as prior art on what a GPU port of this class of solver
  looks like, not as a coupling target for that issue.

**Not assessed**: mesh type/format (unstructured vs. structured), how
its case setup works, or how a UnDiFi-2D-style file-round-trip coupling
would even attach to a Legion task-based runtime (this may be a much
larger integration lift than any of the file-based solvers in #2-#10's
review -- Legion's execution model doesn't obviously support "run one
outer iteration, hand control back" the way EulFS/NEO/SU2's plain
executables do). Flagging the capability, not claiming it's easy to
integrate.

### MUM — sibling architecture (Triangle+EulFS+ALE), not a coupling target either

Like `UnDiFi-3D`, `MUM` ("Moving Unstructured Meshes") is not an
external solver to couple with -- it's a Fortran driver combining
**the exact same Triangle + EulFS pairing UnDiFi-2D itself uses**, but
for a different moving-boundary problem: rigid-body motion (pitching/
translating/orbiting airfoils) via an Arbitrary-Lagrangian-Eulerian
(ALE) formulation, re-triangulating the domain around the body's new
position every physical time step. Its own README already documents
its academic-legacy-code caveats plainly ("largely unmaintained legacy
Fortran... glued together with `system()` calls").

**Relevant how**: MUM's ALE handling of a moving internal boundary
through repeated re-triangulation is structurally the same category of
problem UnDiFi-2D's own shock motion is (a boundary that moves and
forces remeshing every step), just for a rigid body instead of a shock
front. Its `re_mesh`/`mk_ini_mesh`/`mk_blck_ins`/`mk_blck_ext`/
`mk_phpnt_interp`/`finder.f` pipeline (per its own README) is worth a
comparison read against UnDiFi-2D's `wtri`/`fnd_phps`/`interp`/`finder`
(same names!) if E2's closed-curve/hole-topology work (a bubble
interface is also a moving internal boundary requiring an "annular
hole" mesh treatment) hits a design question MUM's rigid-body case may
already have solved -- MUM's bodies are also literally holes in the
mesh that move and must stay validly triangulated. Also directly
relevant to issue #7 (COOLFluiD): the paper #7 cites for MUM's own
origin (Bonfiglioli & Paciorri 2015) is explicitly named there as part
of the same COOLFluiD/UnDiFi lineage.

### ADPDIS3D — structured/curvilinear high-order solver with real chemistry (Mutation)

A high-order, low-dissipation shock-capturing solver (`docs/` cites Yee
2000/2006 -- H.C. Yee's well-known adaptive/low-dissipation high-order
scheme work -- plus scramjet and shock-vortex-interaction test cases),
MPI-parallel (`linux-x86_64-openmpi_{intel,gnu}` build variants), and
vendoring **Mutation 1.3** (`ADPDIS3D/mutation.1.3/`) -- the same
physicochemical library SU2-NEMO links against for nonequilibrium
multi-species chemistry (confirmed via this session's earlier SU2 web
search). A VKI-associated presentation (`docs/PresentationVKI/`) and a
Lani-authored slide deck (`docs/Lani-AIAA-slides.pdf` -- A. Lani leads
COOLFluiD, issue #7) both point to the same VKI/COOLFluiD research
lineage as several codes already in this tree.

**Mesh type**: not confirmed from a quick pass (no clear README stating
structured vs. unstructured), but the "converter" directory
(`bin2vtk`, `plot3d`, `resample`) and the high-order/low-dissipation
scheme literature it cites both suggest a structured/curvilinear code in
the same family as CFDWARP -- i.e. likely **blocked on #2** the same
way, not assessed further here since confirming that needs a real read
of `ADPDIS3D/docs/manual.pdf`/`ADPDIS3D-intro.tex`, not attempted in
this pass. If it does turn out to be structured, it would be a genuinely
strong post-#2 candidate specifically FOR the multi-species/chemistry
question, on the same grounds as HTR above -- real, current chemistry
infrastructure (Mutation) rather than something to build from scratch.

### ramses, amrvac — astrophysics AMR codes, not real coupling candidates

Both are well-known, actively-maintained **astrophysical** simulation
codes: RAMSES (cosmological structure formation, `amr`/`hydro`/`mhd`/
`poisson`/`pm`/`rhd` modules) and MPI-AMRVAC (adaptive-mesh-refinement
magnetohydrodynamics). Neither has any apparent connection to this
project's aerodynamic/gasdynamic shock-fitting domain, engineering
boundary conditions, or the kind of body-fitted/structured mesh coupling
#2-#10 already reviews. **Not recommended as coupling targets.** The one
place they could matter is indirect and much longer-term: both are
mature, real-world examples of **adaptive mesh refinement** applied to
shock-dominated flows (structure-formation shocks, MHD shocks) -- if
this project ever wants a serious "fitted vs. state-of-the-art AMR
capturing" comparison (a stronger, more current benchmark than the
fitted-vs-plain-capturing comparisons this session already did for
issue #25/WP-A), these are real examples of what that comparison would
be up against, not something to integrate. Flagged for completeness
since they were in scope of the ask, not because there's a concrete next
step here.

## 3. Summary recommendation

For the **shock-bubble / E1 multi-species task specifically**: this
note's own finding doesn't overturn this session's earlier
`species-transport-solver-assessment.md` recommendation (SU2 leads), but
adds two real alternatives worth knowing about if SU2's species-
transport feature turns out not to interoperate cleanly with a fitted
(not captured) moving interface: **HTR** (genuine ground-up multi-
species/chemistry, GPU-native, integration effort unknown) and,
pending a mesh-type confirmation, **ADPDIS3D** (real Mutation-based
chemistry, likely structured/blocked-on-#2 like CFDWARP).

For **possible extensions beyond the immediate roadmap**: **UnDiFi-3D**
is the standout finding of this whole survey -- not a candidate to
couple with, but a working 3D implementation of this exact algorithm,
with independently-derived (and in `co_dc`'s case, differently-scoped)
versions of kernels this session touched directly. Worth reading before,
not after, any future 3D-extension work, and worth a direct comparison
read for `co_dc`/`co_oblkshock` regardless of 3D plans, purely as a
second implementation to check UnDiFi-2D's own physics against.
