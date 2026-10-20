# API Reference: qaws_view.h

Curves and surfaces evaluated from host memory without copying it (tests 142
and 143).

---

## Array views

```c
typedef struct qaws_array_view {
    void const* data;          /* packed: element i, component c at data + i * stride + c * size */
    void const* planes[4];     /* planar (data NULL): at planes[c] + i * stride */
    qaws_element_type type;    /* QAWS_ELEMENT_F32 or QAWS_ELEMENT_F64 */
    unsigned int components;   /* 1 to 4 */
    unsigned int count;
    size_t stride;             /* bytes; 0: tight */
} qaws_array_view;             /* qaws_array_view_mut: the same, writable (outputs) */

qaws_status qaws_array_view_validate(qaws_array_view const* view);
qaws_status qaws_array_view_check_extent(qaws_array_view const* view, void const* base, size_t bytes);
```

A view describes one of these layouts:

- a packed array;
- interleaved records: a stride, with the pointer at the member;
- planar arrays: one pointer per component.

Data may be f32 or f64 whatever `qaws_scalar` is. Values are converted one
at a time as they are read, with no converted copy. `_check_extent`
returns OUT_OF_RANGE when the view reads past a host buffer.

## Curves

```c
qaws_curve_view v = { 0 };
v.dimension = QAWS_DIMENSION_3D;
v.degree = 3;
v.points.data = &records[0].position; v.points.type = QAWS_ELEMENT_F32;
v.points.components = 3; v.points.count = n; v.points.stride = sizeof(records[0]);
v.knots.data = knots; v.knots.type = QAWS_ELEMENT_F64; v.knots.components = 1; v.knots.count = n + 4;
/* weights: n positive values with Cartesian points, or points of dimension + 1 homogeneous components */

qaws_curve_prepared* p = NULL;
qaws_curve_prepare(&v, allocator, &p);                       /* or qaws_curve_prepare_into(&v, memory, bytes, &p) */
qaws_curve_prepared_evaluate_3d(p, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
qaws_curve_prepared_evaluate_range(p, &params, begin, end, &outputs);   /* f32/f64 outputs, caller buffers */
qaws_curve_prepared_destroy(p);

qaws_curve* c = NULL;
qaws_curve_create_borrowed(&v, allocator, &c);               /* a qaws_curve, kind QAWS_CURVE_KIND_VIEW */
qaws_curve_create_from_view(&v, allocator, &c);              /* an owned B-spline or NURBS copy */
```

- **Knots.** With a count of 0, the knots are clamped uniform 0, 1, ...,
  n - degree, as `qaws_curve_create_bspline` makes them. With n =
  degree + 1 this is a Bezier curve on [0, 1].
- **Degree.** From 1 to `QAWS_VIEW_MAX_DEGREE` (31).
- **Evaluation.** It is the owned B-spline's and NURBS's to the bit (same
  span search, basis, sums and quotient rule; test 142). Derivatives above
  the degree are valid zeros for polynomial curves and exact for rational
  ones.
- **Prepared object.** It owns only derived data: the knots in
  `qaws_scalar` and the span bounds, in one 16-byte-aligned block. The
  block comes from the allocator, or from caller memory with
  `_prepare_into` (`qaws_curve_prepared_size` bytes, no allocation at all).
- **Ranges.** `_evaluate_range` writes element i of each present output
  (position, d1, d2, d3) for parameters begin to end - 1. Disjoint ranges
  write disjoint elements, so a host's workers may split one call (test
  142: four threads, bit-identical to serial).
- **Borrowed curve.** It works with every operation that evaluates curves:
  sampling, arc length, closest points, frames. Extract, split and reverse
  return owned curves made through an owned copy. Control-point and knot
  inspection, transforms, and derivatives with respect to control points
  refuse it: the host holds those arrays.

### Other families

`view.family` selects what the arrays describe (zero: `QAWS_VIEW_SPLINE`, above). Each evaluates as the owned family it names, bit for bit (test 142):

| Family | Arrays | Domain |
|---|---|---|
| `QAWS_VIEW_HERMITE` | `points` and `tangents` (n each), degree 3 | [0, n - 1], span i on [i, i + 1] |
| `QAWS_VIEW_CATMULL_ROM` | `points`, `parameterization`, `closed`, degree 3 | [0, n - 3] (closed: [0, n]) |
| `QAWS_VIEW_POLYNOMIAL` | `points`: the coefficients c_0 .. c_degree of sum c_i t^i | `knots`: none ([0, 1]) or 2 values |
| `QAWS_VIEW_ARC` | `points`: n centers; `knots`: n elements of 3 components (radius, angle_start, angle_end); 3D: `axes_u`, `axes_v` (n x 3); degree 2 | arc length, [0, total length] |
| `QAWS_VIEW_TRAJECTORY` | `points`: n key positions; `knots`: n key times (strictly increasing, spans of at least 1e-10); `tangents`: none (velocities from the neighbours, as the owned curve) or n key velocities; `closed`; degree 3 | the key times |
| `QAWS_VIEW_YUKSEL` | `points`, `mode` (`qaws_yuksel_mode`), `closed`; degree 2 | [0, n - 1] (closed: [0, n]) |
| `QAWS_VIEW_SUBDIVISION` | `points`, `mode` (`qaws_subdivision_scheme`), `levels` (0: 6), `closed`; degree 2, 3, 4 (Chaikin, Lane-Riesenfeld 3, 4) | one unit span per refined edge |

Their span coefficients are computed at each evaluation from the borrowed values, by the owned curves' routines. Nothing depending on the points, tangents or axes is cached, so these may be rewritten in place like spline points. Arc radii and angles and trajectory times are knots: the span bounds come from them (`QAWS_DEPENDS_KNOTS`). A subdivision view keeps no refined polygon: each evaluation refines only the window of control points its span depends on (a few points per level), with the owned steps' arithmetic. Owned copies (`qaws_curve_create_from_view`, and extract, split and reverse of a borrowed curve) make the owned curve of the same family.

## Surfaces

```c
qaws_surface_view s = { 0 };
s.u_degree = 3; s.v_degree = 2; s.u_count = 6; s.v_count = 5;
s.points = ...;                 /* u_count * v_count points, row major: (i, j) at i * v_count + j */
/* weights, u_knots, v_knots as for curves; default knots clamped uniform over [0, 1], as the owned surfaces */

qaws_surface_prepare(&s, allocator, &sp);
qaws_surface_prepared_evaluate(sp, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
qaws_surface_prepared_evaluate_range(sp, &uv, begin, end, &surface_outputs);   /* position, du, dv, normal, duu, duv, dvv */
qaws_surface_create_borrowed(&s, allocator, &surface);       /* kind QAWS_SURFACE_KIND_VIEW */
qaws_surface_create_from_view(&s, allocator, &surface);
```

Evaluation matches the owned B-spline and NURBS surfaces bit for bit (test
143). A borrowed surface also provides the primal jet to third order, so
closest points, curvatures and offsets work on it. It has no derivatives
with respect to the host's arrays.

## Lifetime contract

- qaws never writes to, frees or keeps host memory that a view points to,
  except as stated below.
- A function taking a view borrows it for the call only.
- A prepared object or a borrowed curve or surface borrows the view's
  arrays for its whole life:
  - the arrays stay valid (not freed, not moved) until it is destroyed;
  - they do not change while an evaluation runs.
- Between evaluations the host may change what the object does not depend
  on. `qaws_*_prepared_dependencies` returns `QAWS_DEPENDS_LAYOUT` (pointers,
  strides, counts, types, degrees), plus `QAWS_DEPENDS_KNOTS` when knots
  were given.
- The points' and weights' values are never depended on. They may be
  rewritten in place without preparing again: an animation frame, or a
  solver writing control points (tests 142 and 143).
- Changing anything that is depended on requires destroying the object and
  preparing again.
- A prepared object is immutable: any number of threads may evaluate it.
- Evaluation never allocates.
