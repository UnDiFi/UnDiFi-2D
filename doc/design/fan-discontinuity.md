# Design note: the `'fan'` discontinuity kind

Issue #23 (E3), task: *"`'fan'`: write a design note first, then implement.
Time-boxed."* This is that note. No implementation in this increment.

## 1. Why a fan doesn't fit the existing model

Every discontinuity kind in this codebase today (`'S'` shock, `'D'` contact,
`'L'` slip, issue #23) shares one data model: `discontinuity_t` is a curve —
a 1-D polyline of points, each carrying exactly two states (`zu`, `zd`, one
on each side) and a jump relation connecting them (`co_shock`/`co_dc`, both
dispatched from `co_state_dps.f90`). A centred (Prandtl-Meyer) expansion fan
is not that. It is a **region**: an angular sector swept between two
characteristic lines (the *head*, where the flow first turns and becomes
supersonic-relative-to-the-corner, and the *tail*, where it finishes turning),
inside which the flow state varies **continuously** with angle, not jumps
once at a single curve. Forcing a fan into the existing "curve with a
two-state jump" model is a category error, not just a missing kind, this is
why the issue calls it "the largest single physics addition in the roadmap."

## 2. Physics recap: what's actually inside a fan

A centred expansion is a self-similar simple wave: every state property is
constant along each ray from the corner, and the whole region is fully
determined, in closed form, by three things:

- the upstream state immediately ahead of the head characteristic
  (`ρ₁, p₁, M₁`, flow direction `θ₁`),
- the total flow deflection `Δθ` (or equivalently the downstream Mach `M₂`,
  related by the **Prandtl-Meyer function** `ν(M) = √((γ+1)/(γ-1)) ·
  atan√((γ-1)/(γ+1)·(M²-1)) − atan√(M²-1)`, with `Δθ = ν(M₂) − ν(M₁)`),
- `γ`.

Given those, **every point inside the fan is a closed-form evaluation, not a
Newton solve**: at local ray angle `φ` between the head and tail angles, the
local Prandtl-Meyer angle is `ν(φ) = ν(M₁) + (φ − φ_head)`, invertible for
`M(φ)`, and `p(φ)/p₁`, `ρ(φ)/ρ₁`, `T(φ)/T₁` follow from the isentropic
relations at that `M(φ)`. This is the opposite of every other discontinuity
kind in this codebase: **the interior needs no iterative solve at all** — the
entire difficulty is at the two bounding curves (getting the head/tail
*angles* right and fitting them as evolving fronts) and at the mesh/hole
topology inside the sector, not at "solving the fan" per se.

## 3. Two discretisation strategies (the issue's own framing)

**(A) Head and tail characteristics as two coupled fronts.** Represent the
fan by two ordinary `discontinuity_t`-shaped curves (the head and tail
characteristics), each carrying **one** state per point (not two — a
characteristic isn't a jump, so `zu`/`zd` collapse to the same value on each
curve, or more precisely: `zd` on the head curve is the local fan-edge state
just inside the fan, and there is no meaningful "upstream of the tail" state
distinct from "just outside the fan" on that side). The angular sector
between them is *not* fitted point-by-point; it is filled analytically:
given the two bounding characteristic angles and the upstream state at a
ray's origin, every interior background-mesh node whose polar angle (from
the fan's origin special point) falls between head and tail gets its state
set directly from the Prandtl-Meyer closed form, the same way `interp.f90`
today sets phantom-node states by lookup rather than by solving anything.

**(B) Fan-angle parameterisation, no explicit interior curves.** Track only
the fan's origin, head angle, and tail angle as scalar state on the special
point that spawns it (e.g. a Guderley triple point); there is no
`discontinuity_t` for the fan at all. Any code that needs to know "am I
inside a fan, and what's the state there" queries the owning special point
directly by angle.

**Recommendation: (A), with a twist — model the fan as ONE `discontinuity_t`
pair (head + tail), not as a first-class region type.** Reasoning:

- The rest of the architecture (mesh regeneration via `wtri`/Triangle, hole
  detection in `fnd_phps`, phantom-node handling in `interp`, redistribution
  in `rd_dps`) is built entirely around "boundaries are curves, meshing
  reacts to curves." Two curves — head and tail — slot into that machinery
  with the least new surface area: the mesh still just needs two more
  internal boundary curves to respect, the same shape of problem it already
  solves for a shock or a slip line's own curve. Option (B) would need a
  parallel, curve-free mechanism grafted alongside the existing one for one
  discontinuity kind only, which is more new machinery, not less.
- Option (A)'s interior-state evaluation is a **pure function of local
  angle and the upstream state carried on the head curve** — no coupling
  back into the outer Newton-solve machinery (`mod_newton_solve.f90`) at
  all, unlike every other kind. This is a genuine simplification: no new
  `residual_if`/`jacobian_if` pair, no new entry in `mod_special_point.f90`'s
  `solve_state` dispatch. The only new numerics is a closed-form Prandtl-
  Meyer evaluation (elementary functions, easily unit-testable against
  tabulated values).
- The real work is bookkeeping, matching the issue's own assessment: the two
  curves must stay geometrically consistent (head and tail both anchored at
  the same origin special point, tail always at a larger or equal deflection
  than head, both advecting with the local characteristic speed
  `u ± a` rather than a shock's jump-condition-derived speed), and the mesh
  region between them needs a **region marker**, not a hole — `wtri.f`
  already distinguishes "hole" (nothing meshed inside, e.g. a solid body)
  from "meshed normally," and would need a third case here: "meshed
  normally, but every node's state gets overwritten analytically afterward,"
  which is closer to E2's already-anticipated "annular hole with an interior
  sub-domain" region-marking work (issue #22) than to anything fan-specific.

## 4. Concrete architecture sketch

- `discontinuity_t%kind`: two more values, `'fan_head'` / `'fan_tail'`
  (or reuse the existing scaffolding's word-based `kind` string finally,
  since the legacy single-char `typesh` realistically can't grow much
  further without collisions — see open question in §6).
- A `fan_head`/`fan_tail` pair share an owning special point (the corner:
  a Guderley triple point today, potentially a plain convex-corner point
  for the general shock-expansion case later). That special point stores
  `M_upstream`, `θ_upstream`, `Δθ` (or `M_downstream`) — the three numbers
  that fully determine the fan — and exposes a pure function
  `fan_state_at_angle(φ) -> (ρ, p, u, v)` implementing §2's closed form.
  This is a natural type-bound procedure on whatever special-point type
  ends up owning a fan (`sp_guderley_t` or similar, not designed here —
  out of scope for this note, which is about the discontinuity kind only).
- Redistribution (`rd_dps`/`rd_dps_eq`) and normal computation (`co_norm`)
  for `'fan_head'`/`'fan_tail'` curves: **treat them as ordinary `'S'`-like
  curves for these two purposes** — a characteristic has a well-defined
  local direction (the local Mach angle relative to the local flow), so the
  existing upwind-biased `shp_dpndnc` dependency test and edge-length-based
  redistribution need no fan-specific logic, only a fan-specific *jump*
  relation (which doesn't exist — see §2, there is no jump to solve on
  these curves, only a state to read off analytically). This is the one
  place `co_state_dps.f90`'s dispatch needs a real new branch: not a new
  Newton solve like `co_shock`/`co_dc`, but a call to
  `fan_state_at_angle` using the point's own polar angle from the fan's
  origin.

## 5. How the shock-capturing solver sees a fan

This question, worth answering explicitly since the issue asks how the
*solver* sees it: **it doesn't need to, and that's a genuine asymmetry worth
recording.** A capturing scheme (EulFS/NEO/SU2 in capturing mode) resolves a
smooth expansion as an ordinary smooth region of the flow field — no special
treatment, no internal boundary, nothing to fit. This is the opposite of a
shock, where capturing smears a true discontinuity and fitting's entire value
proposition is resolving it sharply. For a fan, fitting's job is different:
not to sharpen something capturing blurs (there is nothing sharp to begin
with — a fan is smooth by construction), but to **place the head and tail
characteristics at the geometrically and physically correct angles**, which
a capturing scheme also gets right automatically in principle but can smear
across several cells in practice, especially near the corner. **The
fitted-vs-capturing comparison for a Guderley case is therefore not "sharper
fan edges" (arguably a fan doesn't have edges to sharpen) but "correct
supersonic-patch EXTENT and correct triple-point trajectory in the
von Neumann paradox regime, where the difference between fitting and
capturing is a first-order disagreement about the solution topology, not a
smearing-vs-sharp-edge accuracy question."** This reframing matters for
scoping WP-D's own validation strategy later — worth flagging back to that
issue when GR work actually starts.

## 6. Open questions / risks (not resolved here, flagged for the implementation increment to actually decide)

1. **`typesh` single-character ceiling.** `'S'`/`'D'`/`'L'` (issue #23) are
   already 3 of a very small alphabet (`character*1`, read/written verbatim
   by `re_sdw_info.f90`/`wrt_sdw_info.f90` and every dispatch site). A fan
   needs at minimum 2 more distinguishable values (head/tail). Two
   single-character codes (e.g. `'H'`/`'T'`) fit without widening the type,
   but `'interface'` (also pending, see the companion note) will want at
   least one more on top of that, and Part III's full Ben-Dor catalogue
   (Types III-VI shear layers/jets) will likely want several more still.
   **Recommendation for whoever implements this: budget one increment to
   widen `typesh` from `character*1` to something that can hold a short
   word (or finally wire up `discontinuity_t%kind`'s already-scaffolded
   `character(len=:)` field end-to-end) BEFORE adding the fan kind, not
   after** — retrofitting a width change once 5-6 single-char codes exist
   and are baked into every committed `sh00.dat` is strictly more painful
   than doing it now while only 3 codes exist.
2. **Region marking for `triangle`.** §3/§4 assumes `wtri.f` can mark the
   fan's angular sector as "mesh normally, tag for analytic overwrite."
   Confirming `triangle`'s `-A` region-attribute mechanism (already needed
   for E2/#22's closed-curve annular holes) actually supports this
   "meshed-but-tagged" case, rather than only "hole vs. not-hole," needs a
   short spike before committing to this design — not investigated here.
3. **Fan collapse/degeneracy.** If `Δθ → 0`, head and tail coincide and the
   fan degenerates to a single characteristic (no region at all) — the
   special point owning it needs to detect this and either drop the fan
   curves entirely or keep a zero-width placeholder; not designed here.
4. **This note does not design the owning special point** (`sp_guderley_t`
   or equivalent) or the triple-point-to-fan topology transition that would
   create one — that's WP-D/#27 and #24 (E4 topology transitions)
   territory, deliberately out of scope here per the issue's own framing
   ("this design question should be settled before the Guderley work
   starts, not during it" — this note is that settling, not the Guderley
   work itself).

## 7. Recommended implementation order (for whoever picks this up next)

1. Resolve open question 1 (widen the kind representation) as its own
   small increment — blocks everything else cleanly and is valuable even
   for `'interface'` alone.
2. Implement `fan_state_at_angle` (the Prandtl-Meyer closed form) as a
   standalone, unit-tested pure function against tabulated
   Prandtl-Meyer-angle values — zero mesh/architecture risk, can happen
   before anything else here even lands.
3. Spike open question 2 (region marking) in isolation, on a trivial
   synthetic case, before touching any real mesh-generation code path.
4. Only then: the `'fan_head'`/`'fan_tail'` curve kinds, the owning
   special-point type, and the topology transition that creates one —
   genuinely WP-D/#27 work, not E3/#23's.
