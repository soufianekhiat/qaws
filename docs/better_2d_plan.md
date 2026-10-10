# better_2D: Clipper2 parity for curves, float and exact

Branch `better_2D`, based on `qaws_2_0_0`.

Goal: full feature parity with Clipper2 (v2.0.1, `C:\git\dunya\extern\clipper2`) in
two code paths, on curves as well as polygons:

- **numerical**: `qaws_curve` inputs, float arithmetic, tolerance relative to the
  scene extent;
- **exact**: integer (int64) or lattice inputs, certified predicates, exact topology.

Plus 2.5D stacks of parallel contours and batching everywhere. Scope is 2D and 2.5D only.

---

## 1. Where both libraries stand

### Clipper2 in one paragraph

A Vatti scanline sweep over int64 polygons. Orientation and collinearity
predicates are exact (128-bit), but every new vertex (crossing, horizontal clamp,
self-intersection fix) is **rounded or truncated to the integer grid**, and the
topology is then repaired with tolerances: `CleanCollinear`, `FixSelfIntersects`
(adjacent edges only), collinear joins at 0.25 / 0.35 units, area thresholds 1 and 2,
and a voting heuristic for PolyTree ownership (`Path2ContainsPath1`). Its own tests
accept ±1..5 polygons and 1..50 % area on some cases. So Clipper2 is *robust*, but
its output topology is not exact; ours can be.

### qaws today

| Area | Float | Exact |
|---|---|---|
| Boolean | `qaws_boolean_2d`: one closed curve per side, union / intersection / difference, output sampled 48 pts per arc and refit, always 1 boundary (disjoint loops chained wrongly) | `qaws_exact_boolean_2d`: one simple closed curve per side, pieces of the inputs (curve, t enclosure, reversed) linked into loops |
| Degenerate cases | none handled; hits have no cross / touch / overlap kind | tangency and overlap → `CERTIFICATION_FAILED`, collinear end-to-end touch too |
| Fill rules, XOR, multi-path, holes, tree, open paths | no | no |
| Offset | single curve, sampled 512 in parameter, refit to Catmull-Rom; no joins, no caps, closed in → open out | none |
| Area / orientation / point-in | winding by 256-sample atan2, no on-boundary result | certified winding, linear scan, fails on the boundary |
| Intersections | pairwise: 256×256 samples, absolute 0.1 candidate tolerance; batch: grid + Gauss-Newton, overlaps dropped | pairwise O(nA·nB), self O(n²), batch on the grid; segment hit is an exact rational but returned as a double enclosure |
| Coordinates | f32 or f64 | lattice 2^-20, **coordinates capped at 32 bits**, no int input API, no int64 / i128 fast path |
| Splitting | exact for Bezier, B-spline, NURBS, cubic trajectory; **Hermite interior split is wrong** (tangents not scaled); Catmull-Rom only at knots; rational Bezier, arc, polynomial, clothoid, composite, Yuksel, subdivision unsupported | no sub-curve object; pieces = parent + double enclosure |
| Closedness | Bezier, B-spline, NURBS always report *not closed*, so winding and boolean see them as empty | exact closure test |

Bugs found on the way (fixed in the first step): `qaws_boolean_2d` can read past its
256-hit buffer (the hit call returns the total count), assumes `arcs[0]` comes from A,
chains disjoint loops into one composite, ignores the winding status; the exact
no-crossing difference emits the hole with the outer orientation.

---

## 2. Parity matrix

Target: every row yes in both columns. "R" = roadmap item in section 6.

| Clipper2 feature | Float | Exact | R |
|---|---|---|---|
| Intersection, Union, Difference | partial | partial | 3, 5 |
| Xor | no | no | 3 |
| Fill rules EvenOdd, NonZero, Positive, Negative | no | no | 3 |
| Many subjects + many clips, holes, self-intersecting inputs | no | no | 3, 5 |
| Touching / overlapping / collinear edges | no | refused | 1, 5, 6 |
| Open subjects (polyline clipping) | no | no | 3 |
| Paths output and PolyTree (Parent, Level, IsHole, Area) | no | no | 3 |
| PreserveCollinear, ReverseSolution | no | no | 3 |
| Reusable prepared input (`AddReuseableData`) | prepared sets exist for hits only | same | 3, 9 |
| Z callback | no | no | 3 |
| ClipperD precision / scaling | n/a (float native) | lattice exp | 5 |
| int64 coordinate range (±2^62) | n/a | 32 bits | 5 |
| ClipperOffset: JoinType Square, Bevel, Round, Miter | no | no | 7 |
| EndType Polygon, Joined, Butt, Square, Round | no | no | 7 |
| Miter limit, arc tolerance, delta callback, offset to PolyTree | no | no | 7 |
| RectClip, RectClipLines | no | no | 8 |
| MinkowskiSum, MinkowskiDiff | no | no | 8 |
| Area, IsPositive, PointInPolygon (IsOn / Inside / Outside) | partial | partial | 2 |
| GetBounds, Translate, Scale, MakePath, Ellipse, Length, Distance | partial | no | 2 |
| StripDuplicates, TrimCollinear, SimplifyPath, RamerDouglasPeucker | no | no | 2 |
| PolyTreeToPaths, CheckPolytreeFullyContainsChildren | no | no | 3 |
| Triangulate (beta in Clipper2) | no | no | 10 |
| Flat-array C export API | our API is C already | | 3 |

Beyond Clipper2, kept or added: curved edges never flattened, exact output
topology (no snapping repairs), certified exact mode, snap-rounded integer output
that is guaranteed simple, 2.5D contour stacks, batched / threaded execution, and
derivatives.

---

## 3. Core design

### 3.1 One engine, three geometry backends

Clipper's sweep needs y-monotone straight edges; for curves we use a **planar
arrangement** instead, which also gives exact faces, exact hole ownership and
trivially all fill rules and clip types.

```
inputs ──► edges (curve pieces) ──► intersections (batch grid, executor)
       ──► vertices (identified) ──► split + merge overlaps (edge carries Δwinding per operand)
       ──► half-edge graph (angular order at vertices)
       ──► faces + component nesting ──► winding per face (BFS from the unbounded face)
       ──► classify faces (fill rule × clip type) ──► boundary loops ──► tree / paths
```

The topology core (graph, faces, windings, classification, loop walk, tree) is
written once. It asks a **backend** only these questions:

| Question | float curves | exact polylines (int64) | exact curves |
|---|---|---|---|
| all crossings / touches / overlaps of the edges | float batch + hit kinds | i128 segment tests on the grid | exact batch + overlaps + multiplicity |
| are two vertices the same point | union-find within tolerance | rational equality | algebraic equality (gcd test) |
| order of edges around a vertex | tangent angle, then curvature | exact orient of integer directions | certified sign of cross(C′a, C′b), then curvature, then higher order |
| order of vertices along an edge | parameter compare | rational compare | root compare (refine, gcd on tie) |
| which face contains a point (component nesting) | ray / winding | exact orient crossings | certified winding |

The float polygon path and the exact polygon path then share every line of
topology code, and so does the curve path. A Clipper-style float polygon path
(`ClipperD`: scale to int, run exact, scale back) is just "exact polylines" with a
scale, which gives literal ClipperD parity for free.

### 3.2 Face classification (the Clipper truth tables, exactly)

Each edge carries `(Δs, Δc)`: the change of subject and clip winding when crossing it
from right to left; merged overlapping edges sum their deltas. The unbounded face is
`(0, 0)`; faces get windings by BFS across edges.

```
filled(w) : EvenOdd w odd | NonZero w != 0 | Positive w > 0 | Negative w < 0
inside    : Union fs || fc | Intersection fs && fc | Difference fs && !fc | Xor fs != fc
```

Output boundary = edges between an inside face and an outside face, oriented inside
on the left (positive area), or reversed with `REVERSE_SOLUTION`. Loops are walked
turning most to one side at pinch vertices, so every output loop is simple (Clipper
splits touching polygons the same way). Ownership comes from face adjacency, not a
point-in-polygon vote: a hole's parent is the inside face that surrounds it.

Open subjects add edges with no winding delta; an open piece is kept when its face
satisfies Clipper's open rule (`Intersection: in clip`, `Union: not in subject and
not in clip`, others: `not in clip`).

### 3.3 Output curves: the slicing problem

Every output edge is a **real new curve**, owned by the result, of the simplest
exact kind for the piece. Internally the engine works on `{source, t0, t1}` intervals;
the new curve is built once, at the end. After an intersection the sliced section must
stay on the original curve, 1:1, so the rules are:

1. **Never refit.** The new curve is an exact restriction of the source to
   [t0, t1], built per kind:

   | Source kind | New curve for [t0, t1] |
   |---|---|
   | Bezier | Bezier, de Casteljau twice, parameter re-normalized to [0, 1] |
   | Rational Bezier | rational Bezier, homogeneous de Casteljau (to add) |
   | B-spline, NURBS | same kind, knot insertion at t0 and t1, unused spans dropped |
   | Hermite | Hermite, tangents scaled by the local span length (**fixes** today's interior split) |
   | Catmull-Rom (any alpha) | each span is a cubic Hermite: converted exactly to Bezier spans, then cut |
   | Arc | arc with new start and sweep angles |
   | Polynomial | polynomial on the narrowed domain |
   | Clothoid | clothoid with new origin, heading and κ0, same sharpness |
   | Trajectory | cubic: exact; higher degree: keep key accelerations (to fix) |
   | Composite | composite of the cut first / last segments plus the whole middle ones |
   | Reparameterized | reparameterized wrapper of the cut source |
   | Yuksel, subdivision | no closed form: exact Bezier or B-spline form of each span where the scheme allows (uniform cubic subdivision is a uniform B-spline), otherwise the reparameterized wrapper, which stays exact |

   Straight pieces of polylines come out as degree-1 curves; a run of them from one
   source is one degree-1 B-spline (the polygon of Clipper).

2. **Exact junctions.** The new curve's end points are set to the vertex position
   (the polished intersection in float, the rounded exact point in exact mode) with
   the smallest change to the end control point, so consecutive output curves meet
   bit for bit and loops are watertight. The vertex is also reported.
3. **Re-fusion.** Adjacent pieces of the same source with contiguous parameters are
   merged before the new curve is built, so a curve the operation does not touch comes
   back as an identical copy (identity round trip, tested).
4. **Overlaps between different curves.** A shared stretch maps [s0, s1] on A to
   [u0, u1] on B, possibly non-linearly (a polyline segment lying on a B-spline's
   linear run). The output keeps one representative, subject first (as Clipper's
   `SetZ`), and records which other edge it covers.
5. **Seams.** Closed curves wrap their parameter; a piece across the seam is cut on
   both sides and joined into one curve (same kind when joinable, composite otherwise).
6. **Cusps and zero derivatives** at a vertex fall back to the first non-zero
   derivative for the angular order.
7. **Parameter map.** Each output curve reports its source index and source
   interval, so hit parameters and attributes (z, user data) can be carried over.

---

## 4. Exact path specifics

- **Integer input API**: `qaws_exact_path_create_i64(points, count, closed)` and an
  integer lattice for curves (control points as int64 directly), no float on the way.
- **Fixed-width fast arithmetic**: i128 / i256 with `__int128` or `_mul128`
  intrinsics and a portable fallback. Segment × segment for int64 points: cross
  products in i128, intersection point = rational with i192 numerator / i128
  denominator, comparisons in i256 / i384. The 2048-bit `qaws_exact_int` stays the
  fallback for curves.
- **Polyline degeneracies, all exact**: collinear overlap (interval of both
  segments), collinear end-to-end touch, vertex on edge, shared vertices, zero-length
  edges, spikes.
- **Curve degeneracies**: overlap via `common_spans` + `invert_point` (exists for
  self hits) wired into curve × curve and batch hits as `HIT_OVERLAP` intervals;
  tangency via square-free part gcd(R, R′) (`qaws_exact_poly_gcd` exists), classified
  touch vs crossing by the sign of the other curve's implicit form on both sides;
  algebraic vertex identity and ordering with interval refinement plus gcd on ties.
- **Remove caps**: 256 crossings in the boolean, 32-bit coordinates (62 bits for
  degree 1), O(n²) self hits (put them on the grid).
- **Output**: exact vertices (rational, or curve + root) with double enclosures,
  and an optional **snap-rounded int64 output** (iterated snap rounding in hot
  pixels): topology-preserving and guaranteed simple, which Clipper's rounding is not.
- **Offsets in exact mode**: an offset needs unit normals (√), so the construction
  is float (or certified enclosure) snapped to the lattice, and the clean-up union is
  exact. Same contract as Clipper (approximate construction, exact boolean).

---

## 5. 2.5D: stacks of parallel contours

A 2.5D shape is a **stack**: regions (sets of closed curves, holes allowed) at
heights z0 < z1 < ... , not uniformly spaced, with the shape between two levels
interpolated. Examples:

- a cone: a circle at z = 0 and a point (degenerate contour) at the apex;
- a sphere: 64 circles at non-uniform heights;
- a height map: its contour lines, the region at level z being { h ≥ z }.

```c
typedef struct qaws_stack_level { qaws_scalar z; unsigned int first_path, path_count; } qaws_stack_level;
typedef struct qaws_stack_2d { qaws_stack_level const* levels; unsigned int level_count;
                               qaws_path_2d const* paths; qaws_fill_rule fill_rule;
                               qaws_stack_interp interp; } qaws_stack_2d;
```

**Boolean of two stacks** A op B. The output levels are the union of the input
levels: A at {0, 1} and B at {0, 0.1, 0.5, 0.7} give results at {0, 0.1, 0.5, 0.7, 1}.

1. **Section**: at every output z, each stack gives its region. At its own level it
   is the input region; between levels it is interpolated; outside its z range it is
   empty.
2. **Layer solve**: one 2D boolean per output level, run as independent jobs on the
   executor.
3. **Finalize** (post process): assemble the result stack, keep levels whose region
   is empty only where the shape starts or ends (so the extent stays right), drop
   repeated identical levels if asked, and carry z on every output path.

**Interpolation between levels.** Two modes, chosen per stack:

| Mode | How | Pros | Cons |
|---|---|---|---|
| **Correspondence blend** | match contours between consecutive levels, make them compatible (same kind, degree elevation, knot merging, aligned start point and orientation), blend control points: C(t, s) = (1 − s)·C_k(t) + s·C_k+1(t), with s = (z − z_k)/(z_k+1 − z_k) | output is an exact curve of the input kind; exact in the exact path (rational s) | needs one-to-one contours; fails when a contour splits or merges between levels (height-map saddles) |
| **Distance-field blend** (shape-based interpolation) | the section is the zero set of F(p) = (1 − s)·d_k(p) + s·d_k+1(p), with d the signed distance to each level's region (batch closest points); traced by predictor-corrector, fitted with a certified distance bound | handles splits, merges, holes appearing; a circle → point blend gives exactly the cone's circles (|p| − (1 − s)·r) | output is a fitted curve, not an input kind; float path only (√ in distances) |

The default is correspondence when the contours match one-to-one and the
distance-field blend otherwise.

**Non-linear interpolation, the easy way.** Both modes blend values level by level,
so swap the linear weight for a cubic Hermite in z with non-uniform Catmull-Rom
tangents from the neighbouring levels. A sphere given as 64 circles then follows the
curved profile between levels instead of a faceted one, at no extra structural cost.
Linear stays the default; `interp = QAWS_STACK_CUBIC` turns it on.

**Exact path**: correspondence blend at rational z only. Distance-field sections
are float.

**Degenerate contours**: a point (cone apex) or a segment is a closed contour of
zero area. It takes part in interpolation; a 2D boolean on it alone gives an empty
region, but the level is kept as the end of the shape.

**Also from the stack**: a surface mesh between consecutive levels (contour
stitching, see item 10), and the volume and area of the result from the levels.

---

## 6. Roadmap

Each item lands with tests, a showcase PNG, and docs; float and exact move
together wherever the item applies to both.

1. **Foundations**
   - Fix the boolean bugs listed in section 1.
   - `is_closed` for Bezier, B-spline, NURBS; `qaws_curve_extract(curve, t0, t1)`
     building the exact new curve of section 3.3 for every kind; Hermite interior
     split; reverse for every kind.
   - Replace the Python amalgamation with a CMake script.
   - Hit kinds `CROSSING / TOUCH / OVERLAP(interval)` in pairwise and batch float
     hits; relative candidate tolerance instead of the absolute 0.1.
2. **Paths and utilities** (float + exact)
   - `qaws_path_2d` (ordered edges, closed flag), regions = path lists.
   - Area (exact for polynomial curves by Green's theorem; rational for int64
     polygons), orientation, point-in-path with IsOn, tight bounds, length,
     translate / scale, ellipse, make path, strip duplicates, trim collinear,
     simplify (Visvalingam-style as Clipper), Ramer–Douglas–Peucker; for curves also
     re-fusion and tolerance merge via `qaws_curve_merge_chain`.
3. **Float arrangement engine**: all clip types including Xor, four fill rules,
   multi subjects / clips, self-intersecting inputs, open subjects, paths + tree
   output, PreserveCollinear, ReverseSolution, Z callback, prepared regions,
   executor. `qaws_boolean_2d` becomes a thin wrapper.
4. **Parity corpus**: C loader for Clipper2's `Polygons.txt` (195 cases),
   `Lines.txt`, `Offsets.txt`, `PolytreeHoleOwner*.txt`; their tolerances, their 750
   random cases (seed 42) re-expressed in C; the polytree regression list (#618, #942,
   #957, #973, #987, #720, #777, #831). Mathematica areas for curved cases.
5. **Exact polylines**: int64 input, i128 / i256 arithmetic, exact segment
   arrangement with all collinear cases, the same topology core, snap-rounded output,
   ClipperD-style scaled float entry. Runs the corpus with exact topology.
6. **Exact curves**: overlaps, tangencies, algebraic vertex identity and angular
   order, caps removed, grid in self hits.
7. **Offsets** (float, then exact clean-up): polygon / open path with JoinType
   Square, Bevel, Round, Miter and EndType Polygon, Joined, Butt, Square, Round,
   miter limit, arc tolerance, delta callback (per vertex, and per parameter for
   curves), clean-up union Positive / Negative, orientation kept, tree output. For
   curves: lines and arcs offset exactly in their own kind, other kinds by the
   existing offset within a tolerance; round joins are true arcs.
8. **RectClip, RectClipLines, Minkowski**: rect clip as a fast path (axis-line
   roots per curve, curve × line overlaps already exist in exact); Minkowski sum /
   diff of polygons (quads + NonZero union as Clipper), then curved paths with a
   polygon pattern, then convex curved patterns (convolution at parallel tangents).
9. **Batching and 2.5D stacks**: many independent clip / offset jobs in one call, union of N
   regions in one arrangement, prepared regions queried repeatedly, all on the
   executor with bit-identical results. Stacks: correspondence blend, distance-field
   blend, cubic option, stack booleans level by level, finalize.
10. **Triangulation and meshing**: constrained Delaunay triangulation of regions
    (exact orient / incircle in the exact path), boundary samples placed on the true
    curves to a tolerance; Delaunay refinement (minimum angle, maximum area, size
    field); isotropic remeshing (split, collapse, flip, smooth with boundary vertices
    kept on their curves); contour stitching between stack levels into a surface mesh.

Benchmarks against Clipper2 on its own random-polygon benchmark go in once item 5 is
in, as an optional C++ program outside the library (the library stays C) that links
the local checkout `C:\git\dunya\extern\clipper2` and compares results and timings.

## 6.1 Decisions

- 2D output: real new curves owned by the result, exact restrictions of the sources
  (section 3.3), with the source index and interval reported for each.
- 2.5D output: one result region per output level, z carried on each path.
- Builds: Sharpmake and CMake only, no Python anywhere.

---

## 7. Risks

- Exact tangencies and algebraic vertex identity on curves (item 6) are the
  research-heavy part; polylines (item 5) are not.
- Arrangement memory and time on large polygon sets must stay near Clipper's sweep:
  the grid batch is already the broad phase; the face pass is linear.
- Float vertex identification is the classic weak spot of arrangement codes: one
  union-find tolerance, relative to the scene, used by every predicate.
- The CMake amalgamation still calls `buildsystem/amalgamate.py`; item 1 replaces it
  with a CMake script (`cmake -P`) so the build has no Python.
- Correspondence matching between stack levels (which contour goes to which, start
  point alignment) is where the 2.5D results can surprise; the distance-field blend
  is the fallback and the matching is reported.

---

## 8. Status (2026-10-10)

| Item | State |
|---|---|
| 1. Foundations | done: `qaws_curve_extract` (exact pieces of every kind), reverse for every kind, geometric `is_closed`, Hermite interior split fixed, hit kinds CROSSING / TOUCH / OVERLAP in the curve batch, CMake amalgamation (no Python) with the single-file build fixed |
| 2. Paths and utilities | float done (`qaws_path.h`: area, bounds, point location with fill rules, winding, transforms, ellipse, polyline, Clipper2's strip / trim / simplify / RDP); exact counterparts come with the int64 types (`qaws_path64_area2` so far) |
| 3. Float engine | done: `qaws_clip.h` (every clip type and fill rule, open subjects, nesting, z callback, reverse / preserve collinear); `qaws_boolean_2d` wraps it |
| 4. Parity corpus | Polygons.txt and Lines.txt in test 85: 195 / 195 and 16 / 16, oracle clean in f64 and f32; PolytreeHoleOwner and Offsets.txt pending (with offsets) |
| 5. Exact polylines | done: `qaws_clip64.h` (rational vertices, exact predicates, 195 / 195 on the corpus with an exact oracle); benchmark against Clipper2 in `bench/`; lazy exact arithmetic pending (9-23 times slower than Clipper2 today) |
| 6-10 | next |

Decisions made on the way:

- The float engine splits polylines and composites into one source per
  segment, so every corner is a vertex and doubled-back stretches meet as
  overlaps; line x line hits are computed in double from the end points
  (one rounding), lines evaluated in double whatever the scalar type.
- Loops are cut where they pass a vertex twice (parts touching at a point),
  as Clipper2 does.
- Clipper2's stored path counts depend on its integer rounding; parity is
  checked as: areas within Clipper2's tolerances, its count between ours
  with touching pieces joined (slivers under a unit dropped) and ours, and
  every result checked against the definition by an oracle.
