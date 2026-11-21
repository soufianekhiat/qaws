# qaws_interpolating.h: interpolating splines

Smooth 2D splines solved from points alone. The tangents (and, for the curvature-matching methods, the shapes) come
from a global solve over all the points, so moving one point moves the whole curve a little. Each name exists as
`name_f32` and `name_f64` ([precision.md](precision.md)).

The result is a composite of Bezier segments ([composite.md](composite.md)), parameter domain [0, segment_count].
Consecutive duplicate points are dropped first (and, for a closed curve, trailing points equal to the first).

## qaws_interpolating_method

| Value | Segments | Continuity | Where point i lies |
|---|---|---|---|
| `QAWS_INTERPOLATING_CUBIC_C2` (0) | cubic, equal parameter intervals; an open curve has zero second derivative at its ends | C2 at every point | parameter i |
| `QAWS_INTERPOLATING_HOBBY` | cubic, with Hobby's tangent angles and velocities: close to circular arcs where the points allow | G1 at every point | parameter i |
| `QAWS_INTERPOLATING_KAPPA` | quadratic, one per interior point, each with its curvature maximum at its point | G2 where segments join (G1 where the curve inflects at a join: there only the size of the curvature matches) | inside a segment, not at an integer parameter; segments join between points |
| `QAWS_INTERPOLATING_CUBIC_G2` | cubic, tangent directions and lengths solved so the curvature matches at every point | G2 | parameter i |
| `QAWS_INTERPOLATING_S9` | segments whose curvature is a function of arc length, set by their end angles, the angles solved so the curvature matches at every point; output as cubics, several per point | G2 to the accuracy of the cubic approximation; rounder than a single cubic allows, low curvature variation | at a join between cubics, not at the integer parameter matching its index |

An S9 solve uses about 145 KB of stack.

## qaws_interpolating_spline_desc

```c
typedef struct qaws_interpolating_spline_desc {
    qaws_vec2 const* points;
    unsigned int point_count;
    qaws_interpolating_method method;
    int closed;                        /* non-zero: the last point joins the first */
} qaws_interpolating_spline_desc;
```

| Member | Meaning |
|---|---|
| `points` | the points to pass through, finite |
| `point_count` | at least 1; after duplicates are dropped, at least 2 distinct points open and 3 closed |
| `method` | see above |
| `closed` | non-zero: the last point joins the first |

## qaws_curve_create_interpolating_spline

```c
qaws_status qaws_curve_create_interpolating_spline(
    qaws_curve** out_curve,
    qaws_interpolating_spline_desc const* desc);
```

Solves the spline and creates the composite curve. An open curve of two distinct points is the straight segment for
every method. The points are read during the call only. Destroy the result with `qaws_curve_destroy`
([curve.md](curve.md)); it owns its segments.

| Status | When |
|---|---|
| `QAWS_STATUS_OK` | created |
| `QAWS_STATUS_INVALID_ARGUMENT` | `out_curve`, `desc` or `points` NULL; `point_count` 0; a non-finite coordinate; an unknown method |
| `QAWS_STATUS_DEGENERATE_CURVE` | too few distinct points remain (under 2 open, under 3 closed) |
| `QAWS_STATUS_NUMERICAL_FAILURE` | the solve gave no finite result |
| `QAWS_STATUS_ALLOCATION_FAILURE` | allocation failed |

`*out_curve` is set to NULL once the pointers are checked. There is no allocator-taking form: memory comes from the
thread or default allocator ([alloc.md](alloc.md)).
