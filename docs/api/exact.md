# Exact / certified geometry (`qaws_exact.h`)

Certified predicates on doubles, exact rational curves and surfaces prepared
from the float objects, certified winding numbers and certified closest
points. C only; compiled with `QAWS_ENABLE_EXACT=1`, which the Sharpmake build
always defines (the CMake option of the same name, default ON, mirrors it and
defines it for users too). The model, the
guarantees and the statuses are explained in
[the exact kernel concept page](../exact_kernel.md).

Every name exists as `name_f32` and `name_f64` (see [precision.md](precision.md));
the exact kernel itself works on `double` in both, only the float objects it
prepares from (`qaws_curve`, `qaws_surface`) carry the precision.

Statuses specific to this header (see [status.md](status.md)):
`QAWS_STATUS_EXACT_UNSUPPORTED`, `QAWS_STATUS_EXACT_RANGE_EXCEEDED`,
`QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE`, `QAWS_STATUS_CERTIFICATION_FAILED`.
Non-finite inputs and NULL required pointers return
`QAWS_STATUS_INVALID_ARGUMENT`.

## Predicates

```c
typedef enum qaws_exact_sign
{
	QAWS_EXACT_NEGATIVE = -1,
	QAWS_EXACT_ZERO = 0,
	QAWS_EXACT_POSITIVE = 1
} qaws_exact_sign;

typedef enum qaws_exact_path
{
	QAWS_EXACT_PATH_FILTER = 0,   /* the f64 evaluation and its error bound */
	QAWS_EXACT_PATH_EXACT         /* exact integer evaluation */
} qaws_exact_path;
```

```c
qaws_status qaws_exact_orient2d(qaws_exact_sign* out_sign, qaws_exact_path* out_path,
	double const a[2], double const b[2], double const c[2]);
```

Sign of `(a - c) x (b - c)`: positive when a, b, c turn counterclockwise.

```c
qaws_status qaws_exact_orient3d(qaws_exact_sign* out_sign, qaws_exact_path* out_path,
	double const a[3], double const b[3], double const c[3], double const d[3]);
```

Sign of the volume `(a - d) . ((b - d) x (c - d))`.

```c
qaws_status qaws_exact_compare_ratio(qaws_exact_sign* out_sign, qaws_exact_path* out_path,
	double a, double b, double c, double d);
```

Sign of `a / b - c / d`. `b` and `d` must be non-zero
(`QAWS_STATUS_INVALID_ARGUMENT` otherwise).

Common to the three:

- `out_sign` is required; `out_path` may be NULL and tells which path decided.
- The f64 filter uses Shewchuk's error bounds and decides only a sign it
  proves; it never decides a zero. Otherwise the inputs are decoded exactly
  and the polynomial is evaluated in multi-limb integers. The filter is also
  skipped near underflow or when a term overflows.
- Points sharing a coordinate (all equal in x, say) give a zero column: an
  exact zero decided at once, reported as `QAWS_EXACT_PATH_EXACT`.
- `QAWS_STATUS_EXACT_RANGE_EXCEEDED` when the inputs span more binary
  exponents than the integer budget covers (about 2^900 between the smallest
  and the largest magnitude); `QAWS_STATUS_INVALID_ARGUMENT` for a non-finite
  input.

## Exact curves

### Types

```c
typedef struct qaws_exact_curve qaws_exact_curve;   /* opaque */

typedef enum qaws_numeric_quality
{
	QAWS_NUMERIC_APPROXIMATE = 0,
	QAWS_NUMERIC_CERTIFIED,
	QAWS_NUMERIC_EXACT_RATIONAL,
	QAWS_NUMERIC_UNSUPPORTED
} qaws_numeric_quality;

typedef enum qaws_exact_flag
{
	QAWS_EXACT_FLAG_NONE = 0,
	QAWS_EXACT_FLAG_INPUT_QUANTIZED = 1 << 0,   /* inputs were rounded onto the lattices */
	QAWS_EXACT_FLAG_PREP_QUANTIZED = 1 << 1     /* an approximate preparation was frozen */
} qaws_exact_flag;

typedef struct qaws_exact_desc
{
	int space_exp2;              /* lattice step 2^space_exp2 world units */
	unsigned int coord_bits;     /* usable coordinate magnitude bits (default 26) */
	unsigned int param_bits;     /* parameter lattice bits (default 24) */
	unsigned int weight_bits;    /* weight bits (default 24) */
} qaws_exact_desc;

typedef struct qaws_exact_report
{
	qaws_numeric_quality quality;
	unsigned int flags;                      /* qaws_exact_flag bits */
	unsigned int storage_bits;               /* widest homogeneous control value */
	double max_position_quantization_error;  /* world units */
	double max_weight_quantization_error;    /* relative */
	double parameter_quantization_error;     /* last evaluation */
} qaws_exact_report;
```

The preparation functions set `quality` to `QAWS_NUMERIC_EXACT_RATIONAL`;
the other `qaws_numeric_quality` values are part of the vocabulary but are
not produced by this header.

`qaws_exact_desc` fixes the exact space:

- coordinates are rounded (nearest, ties to even) onto `x = i 2^space_exp2`,
  `|i| < 2^coord_bits`;
- weights are scaled by a common power of two to `weight_bits`;
- knots and parameters are rounded onto `t = T 2^-shift`, the shift chosen
  from the domain and `param_bits`.

Valid ranges: `coord_bits` 1..32, `param_bits` 1..56, `weight_bits` 1..30
(`QAWS_STATUS_INVALID_ARGUMENT` otherwise). Objects used together (a winding
loop, a batch) must share `space_exp2`.

### Functions

```c
void qaws_exact_desc_default(qaws_exact_desc* desc);
```

Defaults: `space_exp2 = -20` (2^-20 world units), 26 coordinate bits, 24
parameter and weight bits.

```c
qaws_status qaws_exact_curve_prepare(qaws_exact_curve** out_curve, qaws_exact_report* out_report,
	qaws_curve const* curve, qaws_exact_desc const* desc);
void qaws_exact_curve_destroy(qaws_exact_curve* curve);
```

Prepares `curve` once into one integer homogeneous Bezier per span. `desc`
may be NULL (defaults); `out_report` may be NULL. On success `*out_curve` is
owned by the caller and freed with `qaws_exact_curve_destroy` (NULL is
accepted); it does not reference `curve` afterwards. On failure
`*out_curve` is NULL.

Supported families, 2D and 3D:

| family | preparation |
|---|---|
| Bezier, rational Bezier | one span |
| B-spline, NURBS | one span per non-empty knot interval (exact blossoms, de Boor triangles with the span ends as arguments, cleared to one denominator) |
| cubic Hermite | integer cubic Bezier per unit span, common factor 3 |
| uniform Catmull-Rom | integer cubic Bezier per unit span, common factor 6 |
| chordal / centripetal Catmull-Rom | the runtime's own per-segment cubics frozen and taken exactly; report flag `QAWS_EXACT_FLAG_PREP_QUANTIZED` |
| polynomial | dyadic coefficients taken exactly; only the domain ends are quantized |
| composite | each segment prepared exactly, segment `i` moved to `[i, i + 1]`; bounds off the lattice are reported as a parameter quantization |

Other families (arcs, clothoids, Yuksel curves, ...) return
`QAWS_STATUS_EXACT_UNSUPPORTED`. Degree above 16, more than 256 knots, a
dimension other than 2 or 3, or a coordinate past `2^coord_bits` lattice
steps return `QAWS_STATUS_EXACT_RANGE_EXCEEDED`; so does any integer that
outgrows the 2048-bit budget. `QAWS_STATUS_ALLOCATION_FAILURE` when memory
runs out.

The report: `flags` has `QAWS_EXACT_FLAG_INPUT_QUANTIZED` when any position,
weight or knot was rounded; `storage_bits` is the widest stored integer;
the quantization errors are the largest rounding of a position (world
units), of a weight (relative) and of a knot / domain end.

```c
qaws_status qaws_exact_curve_evaluate(double* out, qaws_exact_report* out_report,
	qaws_exact_curve const* curve, double t, unsigned int order);
```

Position and derivatives at `t`: `out` receives `order + 1` points of the
curve's dimension, `C, C', ..., C^(order)`, `order <= 3`. `t` is rounded onto
the parameter lattice and clamped to the domain; span selection is
`a <= T < b`, the last span at the domain end. Evaluation is division-free
De Casteljau on integer derivative polygons with the exact quotient rule;
each value is converted to double once, correctly rounded. `out_report` (may
be NULL) receives `quality` and `parameter_quantization_error` (only those
two fields are written). `QAWS_STATUS_INVALID_ARGUMENT` for `order > 3` or a
non-finite `t`; `QAWS_STATUS_EXACT_RANGE_EXCEEDED` when the derivative
outgrows the integer budget (high degrees reach it before `C'''`).

```c
unsigned int qaws_exact_curve_span_count(qaws_exact_curve const* curve);
```

Number of spans, one exact rational Bezier each (0 for NULL).

```c
qaws_status qaws_exact_curve_span_bezier(unsigned int* out_degree, double* out_t0, double* out_t1,
	double* out_points, double* out_weights, qaws_exact_curve const* curve, unsigned int s);
```

Span `s` as a rational Bezier, rounded once: degree, parameter interval
`[t0, t1]`, `degree + 1` control points of the curve's dimension (world
units) and `degree + 1` weights relative to the largest. Every output may be
NULL. `QAWS_STATUS_INVALID_ARGUMENT` for `s >= span_count`.

## Certified winding number

```c
qaws_status qaws_exact_winding_2d(int* out_winding, qaws_exact_curve const* const* pieces,
	double const p[2], unsigned int count);
```

Winding number of the closed 2D loop `pieces[0..count)` around the double
point `p`: the signed crossings of the ray `x > p_x`, decided on the exact
integer Bernstein coefficients of `Y - p_y W` and `X - p_x W` with exact
halving subdivision.

- Each piece's end must equal the next one's start exactly (the last piece
  joins the first); a loop that does not close, a NULL or non-2D piece
  returns `QAWS_STATUS_INVALID_ARGUMENT`.
- Pieces in different exact spaces: `QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE`.
- `QAWS_STATUS_CERTIFICATION_FAILED` when `p` lies on the curve, or so close
  in parameter that 200 halvings do not separate it: the winding number is
  undefined there and never guessed.
- `p` needing more than 1500 bits below the lattice step:
  `QAWS_STATUS_EXACT_RANGE_EXCEEDED`.

## Exact surfaces

```c
typedef struct qaws_exact_surface qaws_exact_surface;   /* opaque */

qaws_status qaws_exact_surface_prepare(qaws_exact_surface** out_surface, qaws_exact_report* out_report,
	qaws_surface const* surface, qaws_exact_desc const* desc);
void qaws_exact_surface_destroy(qaws_exact_surface* surface);
```

Bezier, B-spline and NURBS surfaces become one integer homogeneous Bezier
patch per non-empty knot rectangle, on the same lattices as the curves
(coordinates, weights, knots). Ownership and `desc` / `out_report` handling
as for curves. Other surface families return
`QAWS_STATUS_EXACT_UNSUPPORTED`; degree 0 or above 16, more than 256 knots
or more than 64 control points along u or v return
`QAWS_STATUS_EXACT_RANGE_EXCEEDED`.

```c
qaws_status qaws_exact_surface_evaluate(double* out, double* out_normal, qaws_exact_surface const* surface,
	double u, double v, unsigned int order);
```

Values at `(u, v)` (rounded onto the parameter lattices), exact and rounded
once, as vec3 entries: order 0 gives `S` (3 doubles); order 1 adds `Su, Sv`
(9); order 2 adds `Suu, Suv, Svv` (18). They come from integer derivative
nets and the exact bivariate quotient rule. `out_normal` (may be NULL)
receives the unnormalized `Su x Sv`, exact and rounded once (its unit length
needs a square root and is not offered). `QAWS_STATUS_INVALID_ARGUMENT` for
`order > 2` or a non-finite parameter.

```c
void qaws_exact_surface_patch_count(unsigned int* out_u_count, unsigned int* out_v_count,
	qaws_exact_surface const* surface);
```

Patches along u and v (0 for NULL). Either output may be NULL.

```c
qaws_status qaws_exact_surface_patch_bezier(unsigned int* out_p, unsigned int* out_q, double out_rect[4],
	double* out_points, double* out_weights, qaws_exact_surface const* surface, unsigned int iu, unsigned int iv);
```

Patch `(iu, iv)` rounded once: degrees `p`, `q`; parameter rectangle
`{u0, u1, v0, v1}`; `(p + 1)(q + 1)` vec3 control points, row-major with the
u index major; weights relative to the largest. Every output may be NULL.
`QAWS_STATUS_INVALID_ARGUMENT` for an index out of range.

## Batches

### Descriptors and statistics

```c
#define QAWS_EXACT_BATCH_SELF 1u

typedef struct qaws_exact_batch_desc
{
	qaws_exact_curve const* const* curves;
	unsigned int curve_count;
	unsigned int const* families;          /* optional, one id per curve */
	unsigned int flags;                    /* QAWS_EXACT_BATCH_* */
	qaws_batch_executor const* executor;   /* optional: runs independent work in parallel */
} qaws_exact_batch_desc;

typedef struct qaws_exact_surface_batch_desc
{
	qaws_exact_curve const* const* curves;
	unsigned int curve_count;
	qaws_exact_surface const* const* surfaces;
	unsigned int surface_count;
	qaws_batch_executor const* executor;
} qaws_exact_surface_batch_desc;

typedef struct qaws_exact_ssi_batch_desc
{
	qaws_exact_surface const* const* surfaces;
	unsigned int surface_count;
	unsigned int const* families;          /* optional, one id per surface */
	unsigned int min_depth;                /* subdivision depth before certifying */
	qaws_batch_executor const* executor;
} qaws_exact_ssi_batch_desc;

typedef struct qaws_exact_batch_stats
{
	unsigned int span_count;        /* spans over all curves */
	unsigned int cell_count;        /* grid cells */
	unsigned int candidate_count;   /* span pairs whose boxes overlap (closest: spans or patches solved) */
	unsigned int uncertified_count; /* items that failed certification */
	unsigned int hit_count;
	unsigned int patch_count;       /* surface patches */
} qaws_exact_batch_stats;
```

`QAWS_EXACT_BATCH_SELF` asks an intersection batch to add each curve's
self-intersections; `families` skips pairs with the same id.
`qaws_exact_surface_batch_desc` describes N curves against M surfaces.
`executor` is the one of the float batches ([batch_executor.md](batch_executor.md)).
The kernel has no mutable globals, so results are identical whatever the
executor does.

### Certified closest points on curves

```c
#define QAWS_EXACT_CLOSEST_NONE 0xFFFFFFFFu

typedef struct qaws_exact_closest_point
{
	unsigned int curve;               /* QAWS_EXACT_CLOSEST_NONE without curves */
	double t_lo, t_hi;                /* the nearest point's parameter (lo == hi: exact) */
	double distance_lo, distance_hi;  /* its distance, world units */
	int certified;                    /* the nearest point is proven unique */
} qaws_exact_closest_point;

qaws_status qaws_exact_curve_batch_closest(qaws_exact_closest_point* out_points, qaws_exact_batch_stats* out_stats,
	double const* points, unsigned int point_count, qaws_exact_batch_desc const* desc);
```

For each of `point_count` query points (the curves' dimension, exact
doubles, world units), the nearest point over `desc->curves` (`families` and
`flags` are ignored).

- On a span `C = X / W` the squared distance to `p = P / Q` is stationary at
  the roots of the integer polynomial
  `g(s) = sum_c (Q X_c - P_c W)(X_c' W - X_c W')` (degree `3n - 2`). Its
  roots, isolated exactly, and the span ends are the candidates, each with
  an enclosure of its distance (knots and the start of a closed curve count
  once).
- Spans come from one grid over their sound control boxes, nearest first,
  while they can still win.
- `certified = 1` when the winner's distance enclosure lies below every
  other candidate's (overlapping root intervals are refined first).
  Otherwise the nearest candidate is still reported with `certified = 0`: an
  exact tie (the centre of a circle) or candidates closer than refinement
  can separate. An uncertified answer does not change the status.
- Parameter and distance enclosures are rounded outward.
- `out_stats` (may be NULL): `span_count`, `candidate_count` (spans solved),
  `hit_count`, `uncertified_count`.
- Every output entry starts as `curve = QAWS_EXACT_CLOSEST_NONE`; with no
  curves or no points the call returns `QAWS_STATUS_OK` with those entries.
- Curves of different dimensions, a NULL curve: `QAWS_STATUS_INVALID_ARGUMENT`;
  different exact spaces: `QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE`.
- Query points run in chunks of 4 on the executor.

### Certified closest points on surfaces

```c
typedef struct qaws_exact_surface_closest_point
{
	unsigned int surface;             /* QAWS_EXACT_CLOSEST_NONE without surfaces */
	double u_lo, u_hi, v_lo, v_hi;    /* the nearest point's parameters */
	double distance_lo, distance_hi;  /* its distance, world units */
	int certified;
} qaws_exact_surface_closest_point;

qaws_status qaws_exact_surface_batch_closest(qaws_exact_surface_closest_point* out_points, qaws_exact_batch_stats* out_stats,
	double const* points, unsigned int point_count, qaws_exact_ssi_batch_desc const* desc);
```

For each 3D query point, the nearest point over `desc->surfaces`
(`families` and `min_depth` are ignored). On a patch `S = X / W` the squared
distance to `p = P / Q` is stationary where

    G_u = sum_c (Q X_c - P_c W)(X_c,u W - X_c W_u) = 0
    G_v = sum_c (Q X_c - P_c W)(X_c,v W - X_c W_v) = 0

With `3w - 1 = 0` a certified 3D solver isolates the interior roots on boxes
kept by a branch and bound on the distance; the patch edges, exact curves,
add theirs and the corners (each edge and corner once over the surface).
Certification, tie handling, enclosures, initial `QAWS_EXACT_CLOSEST_NONE`,
statuses and executor chunks are as for curves. `out_stats`: `patch_count`,
`candidate_count` (patches solved), `hit_count`, `uncertified_count`.
Patches up to bicubic.

## Intersection result types

The header also declares the result vocabulary of certified intersections,
Booleans and ray casts. qaws declares no function that fills these types;
they are the shared format for intersection code built on the exact curves
and surfaces. Enclosures are rounded outward; the true value is inside.

```c
typedef enum qaws_exact_hit_kind
{
	QAWS_EXACT_HIT_CROSSING = 0,   /* (t_lo, t_hi) holds exactly one transversal intersection */
	QAWS_EXACT_HIT_POINT,          /* an intersection at t_lo == t_hi exactly (enclosed when not a double) */
	QAWS_EXACT_HIT_OVERLAP         /* the curve lies in the line / plane on [t_lo, t_hi] */
} qaws_exact_hit_kind;

typedef struct qaws_exact_hit
{
	qaws_exact_hit_kind kind;
	double t_lo, t_hi;   /* curve parameters */
} qaws_exact_hit;

typedef struct qaws_exact_pair
{
	qaws_exact_hit_kind kind;
	double a_lo, a_hi;
	double b_lo, b_hi;
	unsigned int multiplicity;  /* 1 transversal, even a touch, 0 not determined */
	int opposite;               /* OVERLAP only: b runs from b_hi to b_lo as a goes up */
} qaws_exact_pair;

typedef struct qaws_exact_batch_hit
{
	unsigned int curve_a;
	unsigned int curve_b;
	qaws_exact_pair pair;       /* a_* on curve_a, b_* on curve_b */
} qaws_exact_batch_hit;
```

`qaws_exact_pair`: POINT has both parameters exact (`lo == hi`); CROSSING
encloses one intersection; OVERLAP means the curves coincide over a's
`[a_lo, a_hi]` and b's `[b_lo, b_hi]`.

```c
typedef struct qaws_exact_piece
{
	unsigned int region;    /* 0: a, 1: b */
	double t0_lo, t0_hi;    /* where the piece starts along its traversal (lo == hi when exact) */
	double t1_lo, t1_hi;    /* where it ends */
	int reversed;           /* traversed against the curve's parameter */
} qaws_exact_piece;

typedef struct qaws_exact_loop
{
	unsigned int first, count;   /* pieces[first .. first + count - 1], each ending where the next starts */
} qaws_exact_loop;
```

A 2D Boolean result: pieces of input boundary curves (region 0 or 1)
between two parameters, grouped into loops.

```c
typedef struct qaws_exact_surface_hit
{
	qaws_exact_hit_kind kind;   /* POINT: u, v and t exact (lo == hi); CROSSING: enclosures */
	double u_lo, u_hi;
	double v_lo, v_hi;
	double t_lo, t_hi;          /* x = p0 + t (p1 - p0) */
} qaws_exact_surface_hit;

typedef struct qaws_exact_curve_surface_hit
{
	double t_lo, t_hi;   /* curve parameter */
	double u_lo, u_hi;
	double v_lo, v_hi;
} qaws_exact_curve_surface_hit;

typedef struct qaws_exact_curve_surface_batch_hit
{
	unsigned int curve;
	unsigned int surface;
	qaws_exact_curve_surface_hit hit;
} qaws_exact_curve_surface_batch_hit;
```

A line / surface hit and a curve / surface hit, with patch parameters and
the line or curve parameter enclosed.

```c
typedef struct qaws_exact_ssi_point
{
	double u1_lo, u1_hi, v1_lo, v1_hi;   /* on surface a */
	double u2_lo, u2_hi, v2_lo, v2_hi;   /* on surface b */
} qaws_exact_ssi_point;

typedef struct qaws_exact_ssi_branch
{
	unsigned int first, count;   /* points[first .. first + count - 1] */
	int closed;                  /* the last point joins the first */
} qaws_exact_ssi_branch;

typedef struct qaws_exact_ssi_batch_branch
{
	unsigned int surface_a;
	unsigned int surface_b;
	qaws_exact_ssi_branch branch;   /* first indexes the shared point buffer */
} qaws_exact_ssi_batch_branch;
```

A surface / surface intersection curve as branches of enclosed parameter
pairs.

```c
typedef struct qaws_exact_ray_hit
{
	unsigned int surface;           /* QAWS_EXACT_CLOSEST_NONE: no hit */
	qaws_exact_surface_hit hit;     /* u, v and t enclosures (t along p0 -> p1) */
	int certified;                  /* the first hit (or the miss) is proven */
} qaws_exact_ray_hit;
```

The first hit of a ray `x = p0 + t (p1 - p0)`, `t >= 0`, over a set of exact
surfaces.
