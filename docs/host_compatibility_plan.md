# Host-compatible curves, surfaces and evaluation kernels

Plan for making qaws the mathematical foundation of Atelier's Curves &
Surfaces package (and a dependency of its Procedural Geometry package, of
Alyaf, and of CPU/Halide/HLSL kernels), written 2026-10-19 after a review of
HEAD. qaws stays dependency-free and usable without Atelier or Alyaf.
Priorities: zero-copy parametric evaluation, reusable kernels, safe prepared
objects.

## 0. Where things stand (review of HEAD)

**Branches.** Development is on `meshing` (2.0.0), which contains
`repulsive`, which contains `better_2D` (which `qaws_2_0_0` holds up to
1aceb17). `parameter_convention` (3.0.0: the Alwan/Suwar/Abaad argument
order, f32 AND f64 in one binary) branched from `repulsive` at 7a0c04c and is
not merged: 14 commits on it, 63 on `meshing` since. A trial merge conflicts
in 12 files (18 hunks); the real cost is the code written since in the old
order: the repulsive additions (~2.5k lines) and meshing (92 files, ~41k
lines, ~87 public and ~108 internal declarations, ~708 calls into
reordered functions).

**Curves and surfaces.** Every constructor copies its data (control points,
knots, weights) into owned storage; objects are immutable after creation
(no setters, no lazy caches). Packed AoS only: no strides, no planar
layout; f32 or f64 per build (`QAWS_SCALAR_IS_FLOAT`), one precision per
binary on `meshing`. The only strided borrowed view is `qaws_field_view`
(differentiation). Swept, ruled, Coons, pipe, loft, Gordon, offset and
trimmed surfaces borrow their curves; composite curves take ownership of
theirs.

**Allocation.** `qaws_allocator` ({alloc(unsigned long), dealloc, user}: no
alignment, no size on free, 32-bit sizes on Win64) is honoured by Bezier,
Hermite, Catmull-Rom, B-spline, NURBS and trajectory curves, Bezier, B-spline
and NURBS surfaces, and the differentiation context. Everything else uses
malloc directly: rational Bezier, Yuksel, arc, polynomial, clothoid,
composite and subdivision curves; all derived operations (split, join,
offset, convert, reverse, fit, export); sampling (a table per call);
traversal; inspect; batches (worker scratch on the executor's threads);
repulsive (166 calls, scratch every step); most surfaces (loft and Gordon
fallbacks and trim queries allocate on every evaluation).

**Kernels.** `src/qaws/core/*.h` are header-only, allocation-free and
written for C, HLSL, GLSL and Halide (`qaws_platform.h`), with fixed arrays
(degree <= 16, derivative order <= 3, knot arrays of 34). The runtime uses
them only for Bezier, polynomial, the cubic spans and arcs; B-spline, NURBS,
rational Bezier and surfaces go through `internal/qaws_internal_basis.c`,
and their cores have no runtime caller and no test. `qaws_prepare.h` makes
caller-owned coefficient arrays for some families. SIMD batch kernels in
`batch/*.c` take raw borrowed arrays but are undeclared and unused.

**Backends, as the code stands (not as the docs say).** C: both precisions,
tested. HLSL: float only; compiled and run on a GPU once, for the repulsive
kernels (D3D11 harness, by hand, no threshold, not in CTest); the eval and
diff `.hlsl` examples are checked by hand only. GLSL: does not compile as
written (`qaws_scalar` undefined there, C casts, the struct macro). Halide:
B-spline, NURBS and surface cores need the knot span at pipeline-build time;
the surface normal block and the sampling/repulsive cores use C `if` on
expressions; no test builds or realizes anything. `docs/feature_per_backends.md`
and `docs/backends.md` claim every core runs on every backend: wrong, to be
replaced by the matrix of section 5.

**Defects found (fixed first, section 9).**
- Bezier curves overflow a stack buffer above degree 16 (no degree check);
  B-spline curves above 16 when D3 is asked; rational Bezier 3D at degree 48
  or more; the internal basis is limited to degree 31.
- `qaws_diff_ops.c:379` has a `static` grid inside `surface_solve`: the
  surface closest-point tangent and adjoint race when called concurrently.
- `qaws_rational_bezier_core.h`: D2 wrong when D1 is not asked; D2 of a
  weighted degree-1 curve returned as zero.
- `qaws_inline.h` says inline curves may be copied with memcpy: they hold
  absolute pointers into their buffer (and the buffer is 1-byte aligned).
- The arc-length reparameterized curve returns D1 with respect to its local
  parameter, unlike every other curve (untested).
- Surface flag values differ between the core and the runtime headers.
- `examples/eval_bspline.comp` searches spans in an uninitialised array.

**Repulsive Curves and Surfaces** (milestones R0-R6, `docs/repulsive_plan.md`)
are implemented: energy objects; explicit energy, gradient and
Hessian-vector product (curves, meshes, curves on a surface, control
points), validated against Mathematica; dense and Barnes-Hut/hierarchical
solvers with step and one-shot entry points, iteration and tolerance
limits, a termination reason and a report; caller-owned position buffers;
backend-neutral pair kernels (`core/qaws_repulsive_core.h`, HLSL parity
1.1e-4); solution JVP/VJP by the KKT system for curves (tested) and meshes
(JVP tested, **VJP untested**). Missing: cancellation (a dense step can take
a second), progress from `*_solve`, an iteration or time budget per step,
the allocator, a public evaluation of the full objective and of the metric,
constraint and projection operators, a public descent core, line-search and
GMRES limits as settings, and a validity report for solution derivatives.

## 1. Three layers

1. **Kernels** (`core/`): pure functions over fixed-size arrays, no
   pointers in the shader dialects, no allocation, deterministic; the same
   source for C, HLSL, GLSL and Halide. What a GPU or a Halide pipeline
   calls, and what every CPU path below calls too (one implementation per
   formula).
2. **Views and prepared evaluators** (new, C): non-owning descriptions of a
   curve or surface over host memory (`qaws_*_view`), and immutable prepared
   objects that add derived data (span tables, coefficients) without copying
   the host's arrays. Evaluation over ranges into caller buffers, no
   allocation, thread-safe. What Atelier's DAG and Alyaf use.
3. **Runtime objects** (existing `qaws_curve`, `qaws_surface`): owning,
   general, every operation. They become thin owners over a view and a
   prepared evaluator (owned storage + the same evaluation code), so owned
   and borrowed results agree by construction.

## 2. Borrowed control-point storage

**A typed strided array, the unit of borrowing:**

```c
typedef enum qaws_element_type { QAWS_ELEMENT_F32 = 1, QAWS_ELEMENT_F64 = 2 } qaws_element_type;

typedef struct qaws_array_view
{
	void const* data;               /* first component of the first element (packed), NULL when planar */
	void const* planes[4];          /* planar: each component's first element; unused when packed */
	qaws_element_type type;
	unsigned int components;        /* 1 (knots, weights), 2, 3, 4 (homogeneous) */
	unsigned int count;
	size_t stride;                  /* bytes between elements (packed: >= components * size; planar: between a plane's values) */
} qaws_array_view;
```

Packed AoS with any stride (interleaved host registers included), or
planar SoA (one pointer per component, as Atelier's registers or a Halide
buffer), f32 or f64 independently of the computing precision (elements are
converted as they are loaded, one at a time: no converted copy).
`qaws_array_view_validate` checks counts, strides (at least the element
size, a multiple of the scalar size), alignment of every pointer to its
scalar type, non-NULL pointers, and optionally the extent against a host
buffer size (`qaws_array_view_validate_extent(view, base, bytes)`).

**Curve and surface views** are plain structs of array views and scalars,
per family, e.g.

```c
typedef struct qaws_bspline_view
{
	qaws_array_view points;         /* count n, components = dimension (2 or 3), or 4 homogeneous */
	qaws_array_view knots;          /* n + degree + 1 values, or count 0: clamped uniform */
	qaws_array_view weights;        /* count 0: non-rational; else n values > 0 */
	unsigned int degree;
	unsigned int flags;             /* periodic, homogeneous points, ... */
} qaws_bspline_view;
```

and likewise Bezier, rational Bezier, Hermite, Catmull-Rom, polynomial,
arc, B-spline/NURBS surfaces and Bezier surfaces. A view is validated once
(`qaws_bspline_view_validate`: counts, knot monotonicity, weights > 0 when
asked, degree within limits) and evaluated directly (each call finds its
span), or prepared (section 3).

**Lifetime contract** (stated in the header, tested in section 8):
- qaws never writes to, frees, or keeps beyond a call any host memory a view
  points to, except as stated for prepared objects.
- A view is borrowed for the duration of each call that takes it. Between
  calls the host may change, move or free the arrays.
- A prepared object (section 3) borrows the view's arrays for its whole
  life: they must stay valid (not freed, not moved) until it is destroyed,
  and must not change while an evaluation with it is running. Which
  changes require re-preparation is a property of the family, queryable:
  `qaws_prepared_dependencies(p)` returns flags (POINTS, WEIGHTS, KNOTS,
  COUNTS); a change to anything not listed may be made between evaluations
  without re-preparing (e.g. a B-spline's span table depends on its knots
  only: control points may move freely).
- Owned runtime objects keep copying (their contract is unchanged);
  `qaws_curve_create_from_view` makes an owned copy only when asked.
  `qaws_curve_create_borrowed` makes a runtime object over a view without
  copying, with the prepared object's contract.

Versioning is the host's: qaws exposes what each prepared object depends
on, and Atelier maps that to its register versions.

## 3. Prepared evaluation

```c
qaws_status qaws_curve_prepare_size(size_t* out_bytes, size_t* out_alignment, qaws_curve_view const* view);
qaws_status qaws_curve_prepare(qaws_curve_prepared** out, qaws_curve_view const* view, qaws_allocator const* allocator);
qaws_status qaws_curve_prepare_into(qaws_curve_prepared** out, void* memory, size_t bytes, qaws_curve_view const* view);   /* no allocation */
void qaws_curve_prepared_destroy(qaws_curve_prepared* p);
```

A prepared object holds the view (pointers into host memory, never copies)
and derived data only: the distinct-knot span table and per-span knot
windows (B-spline, NURBS), power-basis coefficients (Hermite, Catmull-Rom,
trajectory), homogeneous weights where the view gives Cartesian points and
weights, arc-length tables when asked (a separate, optional prepared
measure). It is immutable after preparation and has no lazy caches, so any
number of threads may evaluate it at once. Scratch is the caller's: every
evaluation function either needs none (degree <= the kernel limit) or takes
a scratch buffer sized by `qaws_*_scratch_size`; no hidden global or
thread-local state anywhere (the `static` in `qaws_diff_ops.c` goes).
Invalidation is explicit: the host destroys and re-prepares when a listed
dependency changes; `qaws_curve_prepared_matches(p, view)` checks cheaply
(counts, pointers, a knot checksum when asked) for debug builds.

Evaluation over ranges, into caller buffers, with strides:

```c
qaws_status qaws_curve_prepared_evaluate(qaws_array_view_mut out_position, qaws_array_view_mut out_d1, qaws_array_view_mut out_d2,
	qaws_array_view_mut out_d3, qaws_curve_prepared const* p, qaws_array_view parameters, unsigned int begin, unsigned int end);
```

(outputs first, then the object, then inputs with strides, as the
parameter convention; any output may be absent). Same for surfaces with
(u, v) and the six partials and the normal. Deterministic: the same inputs
give the same bits whatever the range split.

## 4. Allocation and scheduling

**Allocation.** Every constructor, operation, sampling, traversal, batch,
mesh and repulsive entry point gets an allocator (the `_ex` pattern for
constructors; a field in descs; the context-last argument otherwise).
`qaws_allocator` gains, compatibly, an optional `alloc_aligned` and the
size on free (needed by arena hosts); sizes become `size_t`. Each API
documents what it allocates persistently (owned by a returned object) and
what it allocates as scratch (freed before returning), and scratch can be
passed in (`*_scratch_size` + a buffer) where an operation is meant to run
in a worker. The audit (section 0) is the work list; a test injects
failures at every allocation (section 8).

**Scheduling.** qaws starts no threads. Evaluation is ranges over
caller-owned buffers on immutable prepared objects (Atelier splits the
ranges over its own pool); `qaws_batch_executor` stays for operations whose
internal parallelism the caller can lend (batch queries, repulsive). For
iterative algorithms (repulsive solvers, remeshing, quads, CVT, fitting) a
common pattern: a `qaws_progress` struct in the desc (a cancel callback
polled at bounded intervals inside steps, a progress callback, a step
budget), step entry points where the algorithm has natural steps, and
`QAWS_STATUS_CANCELLED` leaving outputs untouched (or the last accepted
state, documented per operation).

## 5. Capability matrix (to be made authoritative)

`docs/capabilities.md` replaces `docs/feature_per_backends.md`: per family
(Bezier, rational Bezier, Hermite, Catmull-Rom, B-spline, NURBS,
polynomial, arc, clothoid, composite, trajectory, Yuksel, subdivision;
Bezier/B-spline/NURBS/swept/ruled/loft/... surfaces) and per operation
(position and derivatives, frames, curvature/torsion, arc length and
inverse, closest point, sampling/tessellation, derivatives with respect to
data), the columns CPU f32, CPU f64, HLSL, GLSL, Halide, derivative order,
preparation, allocation, limits and tolerances, and a **status** per cell:
tested (by which test), compiled only, written only, absent. Generated from
a table in the repository and checked by CI where a toolchain is present
(dxc for HLSL, glslang for GLSL, Halide when installed; the library itself
stays dependency-free, the checks are optional build steps).

Shader evaluation of B-splines, NURBS and surfaces without dynamic
allocation: a fixed shader ABI, a header per curve in a structured buffer
(degree, count, offsets into point/knot/weight buffers, flags) and the
arrays in buffers, a span search reading the knot buffer (bounded binary
search, no 34-knot array), the de Boor kernel on a window loaded from the
buffers (degree <= the compile-time limit), rational by homogeneous
coordinates. Halide: the span as a runtime `Expr` (a bounded unrolled
binary search over a knot `Func`), the window gathered by index.

## 6. Kernels

Reusable, object-free kernels (C and, where the formula allows, shaders),
with documented tolerances and determinism:

| Kernel | Exists | Work |
|---|---|---|
| position and derivatives (all families, to D3) | runtime, cores to D2 | cores to D3; runtime routed through the cores |
| tangent, normal, curvature, torsion | inspect (absolute 1e-12 thresholds), `qaws_*_geometry_eval` | scale-relative degeneracy, one kernel |
| parallel-transport frames | private to pipe surfaces | public kernel (double reflection, Wang et al. 2008) and a prepared frame table |
| arc length and inversion | two implementations that disagree (5-point Gauss, 256-entry table; and the 8-piece CDF sampler) | one, error-controlled (adaptive Gauss-Kronrod per span) and an inverse by safeguarded Newton, as a prepared measure |
| curve and surface closest point | 32 samples + Gauss-Newton (can miss the global minimum); 8x8 grid + Newton | certified by Bezier subdivision and bounds (curves), patch subdivision with bounds (surfaces); batched |
| sampling and adaptive tessellation | uniform/arc-length/adaptive; tessellation quadtree with an absolute "curvature" threshold | chord and angle tolerances, deterministic, into caller buffers with counts first |
| surface evaluation and partials | runtime to D2, cores to D2 | prepared surfaces, D3 via jets where needed |
| sweeps, ruled, loft | swept and pipe by finite differences; loft/Gordon fallbacks allocate per evaluation | analytic derivatives, prepared, allocation-free |
| parameterization | chord-length in fitting; Catmull-Rom alpha shape only | uniform/chord/centripetal for interpolation and fitting, documented |

## 7. Repulsive Curves and Surfaces

Already implemented (section 0). Work: the common progress/cancellation
pattern inside steps (bounded intervals, a dense step included), a
progress callback and intermediate positions from `*_solve`, a step
budget, the allocator for persistent solver state and per-step scratch, a
public evaluation of the full objective (energy plus length, area, volume
and obstacle terms) with its gradient, settings for the line search and
GMRES limits, a validity report on solution derivatives (residual,
conditioning, as `qaws_diff_ops`), and a test for the mesh VJP. Interfaces
are documented as implemented only where tests back them: nothing is
claimed from a header alone.

## 8. Tests

1. Owned versus borrowed parity: every family, same bits.
2. Lifetime: the host mutates points between evaluations (allowed, seen),
   knots without re-preparing (detected by `matches` in debug), frees after
   destroy only (checked with a guarded allocator).
3. f32/f64: data of either type into either precision; tolerances per pair.
4. Packed, strided and planar layouts give the same results.
5. Concurrent evaluation of one prepared object from many host threads (C11
   threads or the platform's), results identical to serial.
6. Allocation-failure injection: an allocator failing at the k-th call, for
   every k, through every constructor and operation: a status, no leak, no
   crash.
7. Prepared-data invalidation and `prepared_dependencies`.
8. CPU versus shader core: the cores compiled as C against the runtime
   (always), HLSL/GLSL/Halide when their toolchains are present.
9. Derivatives against numerical references (central differences, and
   Mathematica fixtures as in `tests/reference/`).
10. Stress: high degree (to the limit, and refused past it), repeated knots,
    zero-length spans, degenerate weights, closed and periodic curves.

## 9. Order of work

1. The defects of section 0 (overflows, the static, the rational core,
   the inline claim, reparam D1, flags, the GLSL example), with tests.
2. Array views and the lifetime contract; B-spline/NURBS views (curves,
   then surfaces), direct evaluation through the cores; parity tests 1-4.
3. Prepared objects for those, ranges, scratch; tests 5, 7.
4. The other families as views; the runtime objects re-based on views and
   prepared evaluators (owned = borrowed by construction).
5. The allocator through every entry point; test 6.
6. The capability table, the shader ABI for B-splines/NURBS/surfaces, the
   GLSL fixes, Halide runtime spans; test 8.
7. The kernels of section 6.
8. Repulsive: section 7.

## 10. Decisions for the user

1. **Base branch.** Merge `parameter_convention` into `meshing` first, so
   the new API is written once in the 3.0.0 order with f32 and f64 in one
   binary (recommended: views of either element type then serve both
   computing precisions, and nothing is written twice); or build on
   `meshing` in the 2.0.0 order and convert later.
2. **Mixed precision.** f32 host data computed in f64 (and the reverse) is
   the default here, conversion as elements are loaded. Is a pure-f32 path
   (f32 data, f32 compute) also required for speed on the CPU?
3. **Shader toolchains in CI.** dxc and glslang (and Halide) as optional
   build steps on machines that have them, or kept to manual checks?

## 11. Status (2026-10-19)

Decisions taken: work on `meshing` in the 2.0.0 order; mixed precision yes
(f32 or f64 host data into either computing precision); no shader CI, the
checks stay manual.

| Item | State | Commits, tests |
|---|---|---|
| Defects of section 0 | fixed: degree limits, the static grid, the rational core, the inline copy, reparam D1 (and D3's coefficient, the length-matched curve's chain rule), surface flags, the GLSL example | d1ec231 (140), fee7393 (141), 4d72562 (144) |
| Array views, curve views, lifetime contract | done (`qaws_view.h`, `docs/api/view.md`) | 91ee140 (142) |
| Prepared curves, ranges, borrowed `qaws_curve`, owned copies | done | 91ee140 (142) |
| Surface views, prepared surfaces, borrowed `qaws_surface` with jets | done | 6978734 (143) |
| Capability matrix | `docs/capabilities.md` | d6b4640, cef91cc |
| GLSL compiles; shader checks | the cores and examples compile (`tests/shaders/check_shaders.ps1`) | 40d6c76 |
| CPU/GPU validation | C cores = runtime (144); GPU = C cores, nine families, by hand (`examples/eval_gpu_parity.c`) | cef91cc |
| Repulsive: cancellation, progress, objective | done for the four solvers | bd8a376 (145) |
| Hermite, Catmull-Rom and polynomial views | done (`qaws_curve_view.family`) | 2593741 (142) |
| The allocator through every constructor | every curve and surface family; loft and Gordon evaluation scratch from it beyond 32 sections | 31e5f0b, 41df199 (146) |
| Surface intersection at its default tolerance | fixed (Gauss-Newton seeds); fillets work | f830321 (48) |

Remaining, in order:

- Views of the remaining families (arc, clothoid, Yuksel, trajectory,
  subdivision, composite); owned objects are not rebuilt on views.
- The allocator through operations, sampling, traversal, inspection,
  batches and repulsive (the constructors are done).
- The shader ABI for B-splines, NURBS and surfaces read from buffers, and
  Halide's runtime span.
- The kernels of section 6.
- Repulsive: the allocator, line-search and GMRES settings, a validity
  report for solution derivatives, and the mesh VJP test.
