# API Reference: qaws_view.h

Borrowed geometry: curves and surfaces evaluated from host memory without
copying it. A view describes arrays the host owns; a prepared object adds
the derived data evaluation needs (knots in the computing precision, the
span table); a borrowed curve or surface wraps a prepared object as an
ordinary `qaws_curve` or `qaws_surface`.

Precision-bearing names exist as `name_f32` and `name_f64`; this page uses
the unsuffixed names (see [precision.md](precision.md)). View data itself
may be f32 or f64 whatever `qaws_scalar` is: values are converted one at a
time as they are read, with no converted copy.

See also [views_and_batches.md](../views_and_batches.md) for the concepts,
[alloc.md](alloc.md) for allocators, [eval.md](eval.md) for the result
types and evaluation flags.

---

## Array views

```c
typedef enum qaws_element_type {
    QAWS_ELEMENT_F32 = 1,
    QAWS_ELEMENT_F64 = 2
} qaws_element_type;

typedef struct qaws_array_view {
    void const* data;          /* packed: element i, component c at data + i * stride + c * size */
    void const* planes[4];     /* planar (data NULL): at planes[c] + i * stride */
    qaws_element_type type;
    unsigned int components;   /* 1 to 4 */
    unsigned int count;
    size_t stride;             /* bytes; 0: tight (packed: components * size; planar: size) */
} qaws_array_view;

typedef struct qaws_array_view_mut {   /* the same, writable: outputs */
    void* data;
    void* planes[4];
    qaws_element_type type;
    unsigned int components;
    unsigned int count;
    size_t stride;
} qaws_array_view_mut;
```

`size` is 4 for f32 and 8 for f64. A view describes a packed array,
interleaved records (a stride, with the pointer at the member), or planar
arrays (one pointer per component). Values are read bytewise: no alignment
is required (doubles 4-aligned on 32-bit x86, packed host structs with any
stride).

```c
qaws_status qaws_array_view_validate(qaws_array_view const* view);
```

Returns `QAWS_STATUS_INVALID_ARGUMENT` unless the type is known, there are
1 to 4 components, and, when `count > 0`, the pointers are non-NULL and the
stride is at least an element (packed) or a value (planar).

```c
qaws_status qaws_array_view_check_extent(qaws_array_view const* view, void const* base, size_t bytes);
```

Returns `QAWS_STATUS_OUT_OF_RANGE` unless every byte the view reads lies in
`[base, base + bytes)`: a check against a host buffer's extent.

## Curve views

```c
#define QAWS_VIEW_MAX_DEGREE 31

typedef enum qaws_curve_view_family {
    QAWS_VIEW_SPLINE = 0,      /* B-spline, NURBS or Bezier (the default) */
    QAWS_VIEW_HERMITE,
    QAWS_VIEW_CATMULL_ROM,
    QAWS_VIEW_POLYNOMIAL,
    QAWS_VIEW_ARC,
    QAWS_VIEW_TRAJECTORY,
    QAWS_VIEW_YUKSEL,
    QAWS_VIEW_SUBDIVISION
} qaws_curve_view_family;

typedef struct qaws_curve_view {
    qaws_dimension dimension;           /* 2 or 3 */
    unsigned int degree;                /* splines: 1 to QAWS_VIEW_MAX_DEGREE; other families: theirs */
    qaws_array_view points;             /* n points: dimension components, or (splines) dimension + 1 homogeneous */
    qaws_array_view weights;            /* splines: count 0 none; else n weights > 0 (with Cartesian points) */
    qaws_array_view knots;              /* splines: count 0 clamped uniform; else n + degree + 1 values */
    qaws_curve_view_family family;      /* zero: QAWS_VIEW_SPLINE */
    qaws_array_view tangents;           /* Hermite: n tangents; trajectory: none or n velocities */
    qaws_parameterization parameterization;   /* Catmull-Rom (DEFAULT: centripetal) */
    int closed;                         /* Catmull-Rom, trajectory, Yuksel, subdivision */
    qaws_array_view axes_u, axes_v;     /* arc, 3D: n plane axes (3 components) */
    unsigned int mode;                  /* Yuksel: qaws_yuksel_mode; subdivision: qaws_subdivision_scheme */
    unsigned int levels;                /* subdivision: refinement levels (0: 6) */
} qaws_curve_view;

qaws_status qaws_curve_view_validate(qaws_curve_view const* view);
```

`qaws_curve_view_validate` checks the arrays and their values: counts,
layouts, knots non-decreasing over a non-empty domain, weights (or
homogeneous w) positive and finite. Every function taking a view runs this
check.

### Splines

- **Points.** n points of `dimension` components, or `dimension + 1`
  homogeneous components (w x, w y, [w z,] w).
- **Weights.** Count 0: none. Otherwise n positive weights, with Cartesian
  points.
- **Knots.** Count 0: clamped uniform 0, 1, ..., n - degree, as
  `qaws_curve_create_bspline` makes them (with n = degree + 1 this is a
  Bezier curve on [0, 1]). Otherwise n + degree + 1 non-decreasing values;
  the domain is `[knots[degree], knots[n]]`.
- **Degree.** 1 to `QAWS_VIEW_MAX_DEGREE` (31).

### Other families

`family` selects what the arrays describe. Each evaluates as the owned
family it names.

| Family | Arrays | Domain |
|---|---|---|
| `QAWS_VIEW_HERMITE` | `points` and `tangents` (n each, `dimension` components); degree 3 (`qaws_curve_create_hermite`) | n - 1 spans, span i on [i, i + 1] |
| `QAWS_VIEW_CATMULL_ROM` | `points`, `parameterization`, `closed` (`qaws_curve_create_catmull_rom`) | n - 3 spans (closed: n), span i on [i, i + 1] |
| `QAWS_VIEW_POLYNOMIAL` | `points`: the coefficients c_0 .. c_degree of sum c_i t^i (count degree + 1) | `knots`: none ([0, 1]) or 2 values, the range |
| `QAWS_VIEW_ARC` | `points`: n centers; `knots`: n elements of 3 components (radius, angle_start, angle_end); 3D: `axes_u`, `axes_v`; degree 2 | arc-length parameter |
| `QAWS_VIEW_TRAJECTORY` | `points`: n key positions; `knots`: n key times (1 component, strictly increasing, spans of at least 1e-10); `tangents`: none (velocities from the neighbours, as the owned curve) or n key velocities; `closed`; degree 3 | the key times |
| `QAWS_VIEW_YUKSEL` | `points`, `mode` (`qaws_yuksel_mode`), `closed`; degree 2 | n - 1 spans (closed: n), span i on [i, i + 1] |
| `QAWS_VIEW_SUBDIVISION` | `points`, `mode` (`qaws_subdivision_scheme`), `levels` (0: 6), `closed`; degree 2, 3 or 4 (the scheme's) | unit spans over the refined polygon |

Span coefficients of these families are computed at each evaluation from
the borrowed values, so points, tangents and axes may be rewritten in place
like spline points. Arc radii and angles and trajectory times are knots:
the span bounds come from them (`QAWS_DEPENDS_KNOTS`). A subdivision view
keeps no refined polygon: each evaluation refines only the window of
control points its span depends on.

## Prepared curves

```c
typedef struct qaws_curve_prepared qaws_curve_prepared;

#define QAWS_DEPENDS_POINTS  1u   /* the points' values */
#define QAWS_DEPENDS_WEIGHTS 2u   /* the weights' values */
#define QAWS_DEPENDS_KNOTS   4u   /* the knots' values */
#define QAWS_DEPENDS_LAYOUT  8u   /* pointers, strides, counts, types, the degree */
```

### qaws_curve_prepare

```c
qaws_status qaws_curve_prepare(qaws_curve_prepared** out_prepared, qaws_curve_view const* view, qaws_allocator const* allocator);
```

Validates the view and builds the prepared object, its memory from
`allocator` (NULL: the default order, see [alloc.md](alloc.md)). The
object borrows the view's arrays for its whole life. Destroy with
`qaws_curve_prepared_destroy`.

### qaws_curve_prepared_size, qaws_curve_prepare_into

```c
size_t qaws_curve_prepared_size(qaws_curve_view const* view);
qaws_status qaws_curve_prepare_into(qaws_curve_prepared** out_prepared, void* memory, qaws_curve_view const* view, size_t bytes);
```

`qaws_curve_prepared_size` returns the bytes `qaws_curve_prepare_into`
needs for this view (alignment 16), or 0 when the view is invalid.
`qaws_curve_prepare_into` builds the object in caller memory, 16-byte
aligned and at least that many bytes, with no allocation at all. The object
lives there until the caller reuses the memory;
`qaws_curve_prepared_destroy` is a no-op for it and may be called.

### qaws_curve_prepared_destroy

```c
void qaws_curve_prepared_destroy(qaws_curve_prepared* prepared);
```

### qaws_curve_prepared_dependencies

```c
unsigned int qaws_curve_prepared_dependencies(qaws_curve_prepared const* prepared);
```

`QAWS_DEPENDS_*` flags: what the host must not change without preparing
again. It returns `QAWS_DEPENDS_LAYOUT`, plus `QAWS_DEPENDS_KNOTS` when the
view has knots. The points' and weights' values are never depended on.

### qaws_curve_prepared_get_range

```c
qaws_range qaws_curve_prepared_get_range(qaws_curve_prepared const* prepared);
```

The parameter domain.

### qaws_curve_prepared_evaluate_2d / _3d

```c
qaws_status qaws_curve_prepared_evaluate_2d(qaws_eval_result_2d* out_result, qaws_curve_prepared const* prepared,
    qaws_scalar parameter, unsigned int eval_flags);
qaws_status qaws_curve_prepared_evaluate_3d(qaws_eval_result_3d* out_result, qaws_curve_prepared const* prepared,
    qaws_scalar parameter, unsigned int eval_flags);
```

One parameter, clamped to the domain. Flags and results as
`qaws_curve_evaluate_*` (`QAWS_EVAL_FLAG_*`). Thread-safe, no allocation.

### qaws_curve_prepared_evaluate_range

```c
typedef struct qaws_curve_outputs {
    qaws_array_view_mut position, d1, d2, d3;
} qaws_curve_outputs;

qaws_status qaws_curve_prepared_evaluate_range(qaws_curve_prepared const* prepared, qaws_array_view const* parameters,
    unsigned int begin, unsigned int end, qaws_curve_outputs const* outputs);
```

Evaluates parameters `begin` to `end - 1` of a one-component view, writing
element i of each output that has a pointer (`data` or `planes[0]`). Each
output has `dimension` components, at least `end` elements, f32 or f64.
Disjoint ranges write disjoint elements, so a host's workers may evaluate
ranges of one call at once. No allocation. Returns the first failing
status, with the outputs written up to it.

## Curves over views

### qaws_curve_create_borrowed

```c
qaws_status qaws_curve_create_borrowed(qaws_curve** out_curve, qaws_curve_view const* view, qaws_allocator const* allocator);
```

A `qaws_curve` over the view (kind `QAWS_CURVE_KIND_VIEW`) without copying
its points, weights or knots: a prepared object inside, from the allocator.
It has the lifetime contract of a prepared object for as long as the curve
lives. Destroy with `qaws_curve_destroy`.

- Works with every operation that evaluates curves: sampling, arc length,
  closest points, frames.
- Extract, split and reverse return owned curves, made through an owned
  copy.
- Operations that read a family's own arrays refuse it: control point and
  knot inspection, transforms, derivatives with respect to control points.
  The host holds those arrays.

### qaws_curve_create_from_view

```c
qaws_status qaws_curve_create_from_view(qaws_curve** out_curve, qaws_curve_view const* view, qaws_allocator const* allocator);
```

An owned copy. Spline views give a B-spline curve (polynomial views) or a
NURBS curve (weights or homogeneous points, divided by w), with the view's
knots or the clamped uniform ones. The other families give their own owned
curve. Values are converted to `qaws_scalar`; the view is borrowed for the
call only.

## Surface views

```c
typedef struct qaws_surface_view {
    unsigned int u_degree, v_degree;    /* 1 to QAWS_VIEW_MAX_DEGREE */
    unsigned int u_count, v_count;      /* at least degree + 1 each */
    qaws_array_view points;             /* u_count * v_count points: 3 components, or 4 homogeneous (w x, w y, w z, w) */
    qaws_array_view weights;            /* count 0: none; else u_count * v_count weights > 0 (with Cartesian points) */
    qaws_array_view u_knots, v_knots;   /* count 0: clamped uniform over [0, 1]; else count + degree + 1 */
} qaws_surface_view;

qaws_status qaws_surface_view_validate(qaws_surface_view const* view);
```

A tensor-product B-spline or NURBS surface. Point (i, j), with
`i < u_count` and `j < v_count`, is element `i * v_count + j` of `points`
(row major, as `qaws_surface_create_bspline`). Default knots are clamped
uniform over [0, 1], as the owned surfaces make them.

## Prepared surfaces

```c
typedef struct qaws_surface_prepared qaws_surface_prepared;

qaws_status qaws_surface_prepare(qaws_surface_prepared** out_prepared, qaws_surface_view const* view, qaws_allocator const* allocator);
size_t qaws_surface_prepared_size(qaws_surface_view const* view);
qaws_status qaws_surface_prepare_into(qaws_surface_prepared** out_prepared, void* memory, qaws_surface_view const* view, size_t bytes);
void qaws_surface_prepared_destroy(qaws_surface_prepared* prepared);
unsigned int qaws_surface_prepared_dependencies(qaws_surface_prepared const* prepared);
qaws_range qaws_surface_prepared_get_u_range(qaws_surface_prepared const* prepared);
qaws_range qaws_surface_prepared_get_v_range(qaws_surface_prepared const* prepared);
```

As for curves. `qaws_surface_prepared_dependencies` returns
`QAWS_DEPENDS_LAYOUT`, plus `QAWS_DEPENDS_KNOTS` when either knot array is
given.

```c
qaws_status qaws_surface_prepared_evaluate(qaws_surface_eval_result* out_result, qaws_surface_prepared const* prepared,
    qaws_scalar u, qaws_scalar v, unsigned int eval_flags);
```

`(u, v)` clamped to the domain; flags (`QAWS_SURFACE_EVAL_*`) and results
as `qaws_surface_evaluate` (a normal also gives du and dv). Thread-safe, no
allocation.

```c
typedef struct qaws_surface_outputs {
    qaws_array_view_mut position, du, dv, normal, duu, duv, dvv;
} qaws_surface_outputs;

qaws_status qaws_surface_prepared_evaluate_range(qaws_surface_prepared const* prepared, qaws_array_view const* parameters,
    unsigned int begin, unsigned int end, qaws_surface_outputs const* outputs);
```

Parameters `begin` to `end - 1` of a two-component (u, v) view; outputs of
3 components each; otherwise as `qaws_curve_prepared_evaluate_range`.

## Surfaces over views

```c
qaws_status qaws_surface_create_borrowed(qaws_surface** out_surface, qaws_surface_view const* view, qaws_allocator const* allocator);
qaws_status qaws_surface_create_from_view(qaws_surface** out_surface, qaws_surface_view const* view, qaws_allocator const* allocator);
```

`qaws_surface_create_borrowed` gives a `qaws_surface` over the view (kind
`QAWS_SURFACE_KIND_VIEW`) with the contract of a prepared object for its
life. It supports evaluation, sampling, area and jets to third order (so
closest points and curvatures); it has no derivatives with respect to the
host's arrays. `qaws_surface_create_from_view` makes an owned B-spline or
NURBS surface.

## Example

```c
qaws_curve_view v = { 0 };
v.dimension = QAWS_DIMENSION_3D;
v.degree = 3;
v.points.data = &records[0].position; v.points.type = QAWS_ELEMENT_F32;
v.points.components = 3; v.points.count = n; v.points.stride = sizeof(records[0]);
v.knots.data = knots; v.knots.type = QAWS_ELEMENT_F64; v.knots.components = 1; v.knots.count = n + 4;

qaws_curve_prepared* p = NULL;
qaws_curve_prepare(&p, &v, NULL);          /* or qaws_curve_prepare_into(&p, memory, &v, bytes) */
qaws_curve_prepared_evaluate_3d(&r, p, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
qaws_curve_prepared_evaluate_range(p, &params, begin, end, &outputs);
qaws_curve_prepared_destroy(p);
```

## Lifetime contract

- qaws never writes to, frees or keeps host memory a view points to,
  except as stated here for prepared objects and borrowed curves and
  surfaces.
- A function taking a view borrows it for the call only.
- A prepared object or a borrowed curve or surface borrows the view's
  arrays for its whole life: they must stay valid (not freed, not moved)
  until it is destroyed, and must not change while an evaluation with it
  runs.
- Between evaluations the host may change what the object does not depend
  on (`qaws_*_prepared_dependencies`). The points' and weights' values may
  be rewritten in place (an animation, a solver writing control points)
  without preparing again. A change to anything depended on requires
  destroying and preparing again.
- The prepared object's own memory comes from the given allocator, or from
  caller memory (`_prepare_into`). Evaluation never allocates.
- A prepared object is immutable: any number of threads may evaluate it.
