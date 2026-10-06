# API Reference: qaws_diff*.h

Differentiation of curves, surfaces, geometric quantities, implicit operations, geometry-building operations, integral functionals and inverse-CDF sampling.

Headers: `qaws_diff_types.h`, `qaws_diff.h`, `qaws_diff_geometry.h`, `qaws_diff_ops.h`, `qaws_diff_map.h`, `qaws_diff_functionals.h`, `qaws_diff_sampling.h`.

- [Types](#types)
- [Context and reports](#context-and-reports)
- [Views and parameter keys](#views-and-parameter-keys)
- [Curves](#curves)
- [Surfaces](#surfaces)
- [Geometry](#geometry)
- [Implicit operations](#implicit-operations)
- [Differential maps](#differential-maps)
- [Integral functionals](#integral-functionals)
- [Inverse-CDF sampling](#inverse-cdf-sampling)

---

## Overview

Every differentiable object exposes typed parameter fields (control points, weights, knots, radius, ...). Tangents and adjoints are stored in views that mirror those fields, so no flat parameter vector is packed by the library. Flattening, bounds and optimizer variables belong to the caller.

- **Tangent**: forward propagation, `y' = J x'`.
- **Adjoint**: reverse propagation, `x_bar += J^T y_bar`.
- **Tangent2**: second directional derivative along the same tangent.

Derivatives with respect to the evaluation coordinate (`t`, or `u`/`v`) and derivatives with respect to parameters are independent axes. Results are returned as spatial jets: the tangent jet holds the derivative of every spatial derivative along the tangent direction. No numerical differentiation is returned through this contract.

General conventions:

- `ctx` may be NULL everywhere: first order, scatter accumulation, `malloc` scratch, no report.
- Tangent outputs are overwritten. Adjoint outputs are accumulated (`+=`); clear them first with `qaws_diff_views_clear` when needed.
- Objects are never modified, so calls on the same object may run concurrently as long as they do not share an adjoint view or a report (`ctx->report` is written).
- Batch entry points are the primary API; single-sample functions are thin wrappers around them.

---

## Types

All types in this section are declared in `qaws_diff_types.h`.

### qaws_diff_field

```c
typedef enum qaws_diff_field
{
	QAWS_FIELD_NONE = 0,
	QAWS_FIELD_CONTROL_POINTS = 1,
	QAWS_FIELD_WEIGHTS = 2,
	QAWS_FIELD_KNOTS = 3,
	QAWS_FIELD_U_KNOTS = 4,
	QAWS_FIELD_V_KNOTS = 5,
	QAWS_FIELD_POINTS = 6,           /* interpolated points (Hermite, Catmull-Rom, ...) */
	QAWS_FIELD_DERIVATIVES = 7,      /* Hermite tangents */
	QAWS_FIELD_COEFFICIENTS = 8,     /* monomial coefficients */
	QAWS_FIELD_KEY_TIMES = 9,
	QAWS_FIELD_CENTER = 10,
	QAWS_FIELD_RADIUS = 11,
	QAWS_FIELD_RADIUS_B = 12,
	QAWS_FIELD_ANGLE_START = 13,
	QAWS_FIELD_ANGLE_END = 14,
	QAWS_FIELD_AXIS_U = 15,
	QAWS_FIELD_AXIS_V = 16,
	QAWS_FIELD_OFFSET_DISTANCE = 17,
	QAWS_FIELD_SCALE = 18,
	QAWS_FIELD_CURVATURE = 19,
	QAWS_FIELD_CURVATURE_RATE = 20,
	QAWS_FIELD_DIRECTION = 21,
	QAWS_FIELD_PARAMETER = 22,       /* an operation input parameter (split parameter, ...) */
	QAWS_FIELD_COUNT
} qaws_diff_field;
```

Identifies one parameter field of an object. Numeric values are stable (they are part of `qaws_param_key`) and the list is append-only.

Scalar fields (one component per element): `WEIGHTS`, `KNOTS`, `U_KNOTS`, `V_KNOTS`, `KEY_TIMES`, `RADIUS`, `RADIUS_B`, `ANGLE_START`, `ANGLE_END`, `OFFSET_DISTANCE`, `SCALE`, `CURVATURE`, `CURVATURE_RATE`, `PARAMETER`. Every other field has as many components as the object dimension (2 or 3). Views whose `components` do not match are rejected with `QAWS_STATUS_INVALID_ARGUMENT` by the evaluation functions.

### qaws_value_type

```c
typedef enum qaws_value_type
{
	QAWS_VALUE_SCALAR = 1,
	QAWS_VALUE_VEC2 = 2,
	QAWS_VALUE_VEC3 = 3
} qaws_value_type;
```

Element type of a field.

### qaws_value_domain

```c
typedef enum qaws_value_domain
{
	QAWS_DOMAIN_GENERIC = 0,
	QAWS_DOMAIN_POSITION,
	QAWS_DOMAIN_DIRECTION,      /* tangent vector / free vector */
	QAWS_DOMAIN_LENGTH,
	QAWS_DOMAIN_ANGLE,
	QAWS_DOMAIN_WEIGHT,
	QAWS_DOMAIN_TIME,
	QAWS_DOMAIN_PARAMETRIC,     /* curve/surface parameter coordinate */
	QAWS_DOMAIN_ARC_LENGTH,
	QAWS_DOMAIN_CURVATURE
} qaws_value_domain;
```

Semantic meaning of a field's values.

### qaws_value_constraint

```c
typedef enum qaws_value_constraint
{
	QAWS_CONSTRAINT_NONE = 0,
	QAWS_CONSTRAINT_POSITIVE,       /* > 0 */
	QAWS_CONSTRAINT_NON_NEGATIVE,   /* >= 0 */
	QAWS_CONSTRAINT_MONOTONIC,      /* elements non-decreasing (knot vectors) */
	QAWS_CONSTRAINT_UNIT_INTERVAL   /* in [0, 1] */
} qaws_value_constraint;
```

Constraint the field values must satisfy for the object to stay valid. Enforcing it during optimization is the caller's job.

### qaws_diff_class

```c
typedef enum qaws_diff_class
{
	QAWS_DIFF_SMOOTH = 0,
	QAWS_DIFF_PIECEWISE_SMOOTH,  /* smooth inside pieces (knot spans, segments) */
	QAWS_DIFF_ACTIVE_SET,        /* smooth while a discrete choice is unchanged */
	QAWS_DIFF_SUBGRADIENT,       /* non-smooth (max/min); a subgradient is returned */
	QAWS_DIFF_DISCRETE,          /* only derivatives under frozen topology exist */
	QAWS_DIFF_UNSUPPORTED
} qaws_diff_class;
```

Differentiability class of an object, field or result, ordered from best to worst.

### qaws_diff_capability

```c
typedef enum qaws_diff_capability
{
	QAWS_CAP_TANGENT    = 1 << 0,
	QAWS_CAP_ADJOINT    = 1 << 1,
	QAWS_CAP_TANGENT2   = 1 << 2,
	QAWS_CAP_DIRECT_HVP = 1 << 3,
	QAWS_CAP_LOCAL_SUPPORT = 1 << 4,  /* qaws_*_local_support is available */
	QAWS_CAP_LINEAR     = 1 << 5      /* output is linear in all fields (weights in support are exact) */
} qaws_diff_capability;
```

Capability bits reported per object (`qaws_curve_get_diff_capabilities`, `qaws_surface_get_diff_capabilities`) and per field (`qaws_field_desc::capabilities`).

### qaws_diff_freeze

```c
typedef enum qaws_diff_freeze
{
	QAWS_FREEZE_NONE       = 0,
	QAWS_FREEZE_SPAN       = 1 << 0,   /* knot span / segment selection */
	QAWS_FREEZE_ACTIVE_SET = 1 << 1,   /* closest primitive, branch, extremum */
	QAWS_FREEZE_TOPOLOGY   = 1 << 2,   /* trim loops, boolean topology, mesh connectivity */
	QAWS_FREEZE_SAMPLE_SET = 1 << 3    /* sample counts / adaptive sampling decisions */
} qaws_diff_freeze;
```

Discrete states a derivative may rely on. The caller states which ones it accepts in `qaws_diff_context::frozen`; the report lists the ones actually used in `qaws_diff_report::frozen_used`.

### qaws_diff_validity

```c
typedef enum qaws_diff_validity
{
	QAWS_DIFF_VALID = 0,
	QAWS_DIFF_VALID_LOCALLY,     /* valid while the frozen discrete state is unchanged */
	QAWS_DIFF_AT_BOUNDARY,       /* evaluated exactly on a piece boundary: one-sided */
	QAWS_DIFF_ILL_CONDITIONED,   /* implicit solve near singular */
	QAWS_DIFF_AMBIGUOUS,         /* several competing solutions */
	QAWS_DIFF_INVALID
} qaws_diff_validity;
```

Runtime validity of a derivative, ordered from best to worst.

### qaws_diff_accumulation

```c
typedef enum qaws_diff_accumulation
{
	QAWS_ACCUMULATE_SCATTER = 0, /* each sample adds into the parameter adjoints */
	QAWS_ACCUMULATE_GATHER,      /* each parameter gathers its samples (deterministic order) */
	QAWS_ACCUMULATE_TILED        /* tile-local accumulation, one flush per tile */
} qaws_diff_accumulation;
```

Strategy used by batch adjoint functions to accumulate into parameter views. All three produce the same sum up to floating-point ordering:

- `SCATTER` -- each sample adds its contributions directly into the views.
- `GATHER` -- contributions are bucketed per parameter element, then each element sums its bucket in sample order. Deterministic and free of write conflicts.
- `TILED` -- samples are processed in tiles of `qaws_diff_context::tile_size` (0 = 256); each tile accumulates into a local buffer and flushes the touched range once.

`GATHER` and `TILED` fall back to `SCATTER` when a view set has more than 16 fields. They allocate scratch through `qaws_diff_context::allocator`.

### qaws_coordinate_kind

```c
typedef enum qaws_coordinate_kind
{
	QAWS_COORDINATE_PARAMETRIC = 0,
	QAWS_COORDINATE_NORMALIZED,
	QAWS_COORDINATE_ARC_LENGTH,
	QAWS_COORDINATE_TIME,
	QAWS_COORDINATE_WORLD_DISTANCE
} qaws_coordinate_kind;
```

Kind of the evaluation coordinate of an object. See `qaws_curve_get_coordinate_kind`.

### qaws_field_desc

```c
typedef struct qaws_field_desc
{
	qaws_diff_field field;
	qaws_value_type value_type;
	unsigned int count;               /* number of elements */
	qaws_value_domain domain;
	qaws_value_constraint constraint;
	qaws_scalar default_lower;        /* -inf/+inf encoded as -/+QAWS_DIFF_UNBOUNDED */
	qaws_scalar default_upper;
	qaws_diff_class diff_class;
	unsigned int capabilities;        /* qaws_diff_capability bits for this field */
} qaws_field_desc;

#define QAWS_DIFF_UNBOUNDED ((qaws_scalar)1e30)
```

Schema entry for one field, returned by `qaws_curve_describe_fields` and `qaws_surface_describe_fields`.

- `field`, `value_type`, `count` -- identity, element type and number of elements.
- `domain`, `constraint` -- meaning of the values and the constraint they must satisfy.
- `default_lower`, `default_upper` -- suggested bounds; `-QAWS_DIFF_UNBOUNDED` / `+QAWS_DIFF_UNBOUNDED` mean unbounded.
- `diff_class`, `capabilities` -- how this field can be differentiated.

### qaws_param_key

```c
#define QAWS_PARAM_KEY_MAX_DEPTH 4

typedef struct qaws_param_key
{
	unsigned short child[QAWS_PARAM_KEY_MAX_DEPTH];
	unsigned char depth;      /* number of valid child indices */
	unsigned char field;      /* qaws_diff_field */
	unsigned char component;  /* 0..2 */
	unsigned char components; /* components of the field; 1 = scalar (no component in text) */
	unsigned int element;
} qaws_param_key;
```

Stable structural identity of one scalar parameter: child path, field, element, component. Text form: `child[0]/control_points/4/z`.

- `child`, `depth` -- path through derived objects (up to 4 levels).
- `field` -- a `qaws_diff_field` value.
- `element` -- element index within the field.
- `component`, `components` -- component index and the field's component count; when `components` is 1 the text form has no component.

### qaws_field_view

```c
typedef struct qaws_field_view
{
	qaws_diff_field field;
	qaws_scalar* data;
	unsigned int count;               /* elements */
	unsigned int components;          /* 1, 2 or 3 */
	unsigned int stride;              /* scalars between elements; 0 = components */
	unsigned char const* active;      /* optional per-element activity, NULL = all active */
	unsigned int component_mask;      /* bit c set = component active; 0 = all */
} qaws_field_view;
```

Addresses one field of tangent or adjoint storage. Element `i`, component `c` is at `data[i * stride + c]`.

- `active` -- optional per-element mask; inactive elements read as zero and are not written.
- `component_mask` -- optional per-component mask; 0 means all components are active.
- A view with `data == NULL` is treated as absent.

### qaws_diff_views

```c
typedef struct qaws_diff_views
{
	qaws_field_view* fields;
	unsigned int field_count;
	struct qaws_diff_views* children;  /* for derived objects, indexed like child keys */
	unsigned int child_count;
} qaws_diff_views;
```

Set of field views for one object. Missing fields are treated as zero tangents or inactive adjoints. `children[i]` holds the views of child `i` of a derived object (see `qaws_diff_child`).

### qaws_diff_child

```c
typedef struct qaws_diff_child
{
	qaws_curve const* curve;
	struct qaws_surface const* surface;
} qaws_diff_child;
```

One object a derived object (offset, extrusion, ruled or swept surface, composite curve, ...) chains its rules into. Child `i` is addressed as `child[i]` in parameter keys and `views->children[i]` in tangent and adjoint storage. Exactly one pointer is set.

### qaws_curve_jet_2d / 3d

```c
#define QAWS_CURVE_JET_ORDER 3

typedef struct qaws_curve_jet_2d
{
	qaws_vec2 d[QAWS_CURVE_JET_ORDER + 1];
	unsigned int channels;
} qaws_curve_jet_2d;

typedef struct qaws_curve_jet_3d
{
	qaws_vec3 d[QAWS_CURVE_JET_ORDER + 1];
	unsigned int channels;
} qaws_curve_jet_3d;
```

Curve spatial jet: `d[k] = d^k C / dt^k` for `k = 0..3`. `channels` uses `QAWS_EVAL_FLAG_POSITION`, `QAWS_EVAL_FLAG_D1`, `QAWS_EVAL_FLAG_D2`, `QAWS_EVAL_FLAG_D3` (bit `k` = order `k`) and marks which entries are meaningful.

### qaws_surface_jet_channel

```c
typedef enum qaws_surface_jet_channel
{
	QAWS_SJET_P   = 1 << 0,
	QAWS_SJET_U   = 1 << 1,
	QAWS_SJET_V   = 1 << 2,
	QAWS_SJET_UU  = 1 << 3,
	QAWS_SJET_UV  = 1 << 4,
	QAWS_SJET_VV  = 1 << 5,
	QAWS_SJET_UUU = 1 << 6,
	QAWS_SJET_UUV = 1 << 7,
	QAWS_SJET_UVV = 1 << 8,
	QAWS_SJET_VVV = 1 << 9,
	QAWS_SJET_ORDER1 = QAWS_SJET_P | QAWS_SJET_U | QAWS_SJET_V,
	QAWS_SJET_ORDER2 = QAWS_SJET_ORDER1 | QAWS_SJET_UU | QAWS_SJET_UV | QAWS_SJET_VV,
	QAWS_SJET_ORDER3 = QAWS_SJET_ORDER2 | QAWS_SJET_UUU | QAWS_SJET_UUV | QAWS_SJET_UVV | QAWS_SJET_VVV
} qaws_surface_jet_channel;
```

Surface jet channels: every partial derivative up to total order 3. `ORDER1`, `ORDER2`, `ORDER3` select all partials up to that order.

### qaws_surface_jet

```c
#define QAWS_SURFACE_JET_COUNT 10

typedef struct qaws_surface_jet
{
	qaws_vec3 d[QAWS_SURFACE_JET_COUNT];
	unsigned int channels;
} qaws_surface_jet;
```

`d[i]` is the channel with bit `i`: `d[0]=S`, `d[1]=Su`, `d[2]=Sv`, `d[3]=Suu`, `d[4]=Suv`, `d[5]=Svv`, `d[6]=Suuu`, `d[7]=Suuv`, `d[8]=Suvv`, `d[9]=Svvv`. `channels` is a mask of `qaws_surface_jet_channel` bits.

### qaws_support_kind

```c
typedef enum qaws_support_kind
{
	QAWS_SUPPORT_LOCAL = 0,   /* a few contiguous elements per field */
	QAWS_SUPPORT_WIDE,        /* many elements, still bounded */
	QAWS_SUPPORT_GLOBAL,      /* every element may contribute */
	QAWS_SUPPORT_IMPLICIT     /* defined through a solve */
} qaws_support_kind;
```

How far the influence of the parameters on one evaluation extends.

### qaws_support_range / qaws_local_support

```c
#define QAWS_DIFF_MAX_ORDER 5
#define QAWS_DIFF_MAX_SUPPORT 32
#define QAWS_DIFF_MAX_RANGES 4

typedef struct qaws_support_range
{
	qaws_diff_field field;
	unsigned int first;
	unsigned int count;
} qaws_support_range;

typedef struct qaws_local_support
{
	qaws_support_kind kind;
	unsigned int range_count;
	qaws_support_range ranges[QAWS_DIFF_MAX_RANGES];
	qaws_diff_field weight_field;      /* QAWS_FIELD_WEIGHTS for rational curves, else NONE */
	int has_weights;
	unsigned int order;
	qaws_scalar weights[QAWS_DIFF_MAX_RANGES][QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
} qaws_local_support;
```

Which curve parameters influence one evaluation. Each range `r` covers elements `first .. first + count - 1` of `field`. Weights are given for derivative orders `0..order`.

- For linear families (`QAWS_CAP_LINEAR`) the weights are exact: `d^k C / dt^k = sum_j weights[r][k][j] * field_r[first_r + j]`.
- For rational families (`weight_field != QAWS_FIELD_NONE`) the weights are the homogeneous basis: `A^(k) = sum_j weights * w_j * P_j` and `W^(k) = sum_j weights * w_j`, with `C = A / W`.

### qaws_surface_support

```c
typedef struct qaws_surface_support
{
	qaws_support_kind kind;
	qaws_diff_field field;
	unsigned int u_first, u_count;
	unsigned int v_first, v_count;
	unsigned int u_stride, v_stride;
	int on_boundary;              /* (u,v) lies exactly on an interior knot line */
	qaws_diff_field weight_field;      /* QAWS_FIELD_WEIGHTS for rational surfaces, else NONE */
	int has_weights;
	unsigned int order;
	qaws_scalar u_weights[QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
	qaws_scalar v_weights[QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
} qaws_surface_support;
```

Tensor-product surface support. The weight of element `(i, j)` for the partial `d^(a+b) S / du^a dv^b` is `u_weights[a][i] * v_weights[b][j]`. Elements are addressed as `(u_first + i) * u_stride + (v_first + j) * v_stride`. `on_boundary` is set when `(u, v)` lies exactly on an interior knot line.

---

## Context and reports

### qaws_diff_report

```c
typedef struct qaws_diff_report
{
	qaws_diff_class diff_class;       /* worst class touched */
	qaws_diff_validity validity;      /* worst runtime validity */
	unsigned int frozen_used;         /* qaws_diff_freeze bits the result relies on */
	qaws_scalar condition_number;     /* largest seen, 1 when not applicable */
	qaws_scalar residual;             /* largest implicit residual */
	qaws_scalar branch_gap;           /* smallest gap to a competing solution */
	unsigned int worst_index;         /* sample index of the worst validity */
	unsigned int evaluation_count;
} qaws_diff_report;
```

Optional diagnostics filled through `qaws_diff_context::report`. Reports accumulate across calls: each recorded evaluation keeps the worst class and validity, ORs `frozen_used`, keeps the largest condition number and residual and the smallest branch gap, and increments `evaluation_count`. `worst_index` is the sample index at which the worst validity was first seen.

### qaws_diff_context

```c
typedef struct qaws_diff_context
{
	unsigned int order;                     /* 1 = first order, 2 = also second */
	unsigned int frozen;                    /* qaws_diff_freeze bits the caller accepts */
	qaws_diff_accumulation accumulation;
	unsigned int tile_size;                 /* samples per tile, 0 = default */
	qaws_allocator const* allocator;        /* scratch allocations, NULL = malloc */
	qaws_diff_report* report;               /* optional output */
} qaws_diff_context;
```

Options shared by all differentiation calls. Passing NULL in place of a context is equivalent to an initialized context with no report.

### qaws_diff_context_init

```c
void qaws_diff_context_init(qaws_diff_context* ctx);
```

Zeroes the context and sets `order = 1`, `frozen = QAWS_FREEZE_NONE`, `accumulation = QAWS_ACCUMULATE_SCATTER`. `tile_size`, `allocator` and `report` are left at zero / NULL. Does nothing if `ctx` is NULL.

### qaws_diff_report_reset

```c
void qaws_diff_report_reset(qaws_diff_report* report);
```

Resets a report to its best state: `diff_class = QAWS_DIFF_SMOOTH`, `validity = QAWS_DIFF_VALID`, `condition_number = 1`, `branch_gap = QAWS_DIFF_UNBOUNDED`, all other fields zero. Call it before a sequence of calls whose diagnostics should be collected together. Does nothing if `report` is NULL.

---

## Views and parameter keys

### qaws_field_view_make

```c
qaws_field_view qaws_field_view_make(
	qaws_diff_field field,
	qaws_scalar* data,
	unsigned int count,
	unsigned int components);
```

Returns a view over caller storage: element `i`, component `c` at `data[i * components + c]`. `stride` is set to `components`; `active` and `component_mask` are cleared (everything active).

### qaws_diff_views_find

```c
qaws_field_view* qaws_diff_views_find(
	qaws_diff_views const* views,
	qaws_diff_field field);
```

Returns the first view of `field` with non-NULL `data`, or NULL when the field is absent (inactive) or `views` is NULL. Child views are not searched.

### qaws_diff_views_clear

```c
void qaws_diff_views_clear(qaws_diff_views* views);
```

Zeroes every active element and active component of every view, recursively including children. Use before an adjoint call when the views should hold only that call's result.

### qaws_param_key_make

```c
qaws_param_key qaws_param_key_make(
	qaws_diff_field field,
	unsigned int element,
	unsigned int component,
	unsigned int components);
```

Builds a key with an empty child path.

### qaws_param_key_prepend_child

```c
qaws_status qaws_param_key_prepend_child(
	qaws_param_key* key,
	unsigned int child_index);
```

Prefixes the key with a child index (used by derived objects). Existing child indices shift by one.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `key` is NULL or `child_index > 0xFFFF`; `QAWS_STATUS_OUT_OF_RANGE` if the key already has `QAWS_PARAM_KEY_MAX_DEPTH` children.

### qaws_param_key_compare

```c
int qaws_param_key_compare(
	qaws_param_key const* a,
	qaws_param_key const* b);
```

Total order on keys. Returns -1, 0 or 1. Compares depth, then child indices, then field, element and component. `components` is not compared.

### qaws_param_key_to_string

```c
unsigned int qaws_param_key_to_string(
	qaws_param_key const* key,
	char* buffer,
	unsigned int capacity);
```

Formats the key, e.g. `child[0]/control_points/4/z`. Components 0..2 print as `x`, `y`, `z`; no component is printed when `components` is 1. Writes at most `capacity` bytes (always terminated when `capacity > 0`; `buffer` may be NULL). Returns the length the full string needs, excluding the terminator, so the call can be used to size a buffer.

### qaws_param_key_parse

```c
qaws_status qaws_param_key_parse(
	char const* text,
	qaws_param_key* out_key);
```

Parses the text form produced by `qaws_param_key_to_string`. Accepts zero or more `child[i]` segments, a field name, an element index, and an optional component (`x`/`y`/`z` or `0`/`1`/`2`). When a component is present, `components` is set to 3; otherwise to 1.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on malformed text or an unknown field name; `QAWS_STATUS_OUT_OF_RANGE` for more than `QAWS_PARAM_KEY_MAX_DEPTH` children or a child index above `0xFFFF`.

### qaws_diff_field_name / qaws_diff_field_from_name

```c
char const* qaws_diff_field_name(qaws_diff_field field);
qaws_diff_field qaws_diff_field_from_name(char const* name);
```

Convert between a field and its lower-case name (`"control_points"`, `"weights"`, `"knots"`, `"u_knots"`, ..., `"parameter"`). `qaws_diff_field_name` returns `"unknown"` for out-of-range values. `qaws_diff_field_from_name` returns `QAWS_FIELD_NONE` for NULL or unknown names.

---

## Curves

Declared in `qaws_diff.h`.

### qaws_curve_get_diff_capabilities

```c
unsigned int qaws_curve_get_diff_capabilities(qaws_curve const* curve);
```

Returns the `qaws_diff_capability` bits of the curve family, or 0 if the curve is NULL or not differentiable.

### qaws_curve_get_diff_class

```c
qaws_diff_class qaws_curve_get_diff_class(qaws_curve const* curve);
```

Returns the differentiability class of the curve, or `QAWS_DIFF_UNSUPPORTED` if the curve is NULL or not differentiable.

### qaws_curve_get_coordinate_kind

```c
qaws_coordinate_kind qaws_curve_get_coordinate_kind(qaws_curve const* curve);
```

Returns the kind of the evaluation coordinate `t`: `QAWS_COORDINATE_ARC_LENGTH` for reparameterized curves, `QAWS_COORDINATE_PARAMETRIC` otherwise.

### qaws_curve_describe_fields

```c
qaws_status qaws_curve_describe_fields(
	qaws_curve const* curve,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count);
```

Writes the schema of the curve's differentiable fields. `*out_count` receives the number of fields. Pass `out_fields = NULL` to query the count only.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_BUFFER_TOO_SMALL` if `capacity` is smaller than the field count; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not differentiable; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `curve` or `out_count`.

### qaws_curve_read_field

```c
qaws_status qaws_curve_read_field(
	qaws_curve const* curve,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count);
```

Copies the primal values of one field (`count * components` scalars, element-major). `*out_scalar_count` receives the scalar count. Pass `out_values = NULL` to query the size only.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_BUFFER_TOO_SMALL` if `capacity` is too small; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not differentiable; the family's error status if the field does not exist.

### qaws_curve_clone_with_fields

```c
qaws_status qaws_curve_clone_with_fields(
	qaws_curve const* curve,
	qaws_diff_views const* values,
	qaws_curve** out_curve);
```

Creates a new curve of the same family with fields replaced by the values in `values`; absent fields are copied unchanged. Typical use: apply an optimizer step. The caller owns the returned curve.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_UNSUPPORTED_OPERATION` for derived curves and families without rebuild support; `QAWS_STATUS_INVALID_ARGUMENT` on NULL arguments or a view whose count or component count does not match the field.

### qaws_curve_diff_children

```c
qaws_status qaws_curve_diff_children(
	qaws_curve const* curve,
	qaws_diff_child* out_children,
	unsigned int capacity,
	unsigned int* out_count);
```

Lists the objects a derived curve chains into (composite segments, ...). Child `i` matches `views->children[i]`. Curves without children return `QAWS_STATUS_OK` with `*out_count = 0`. Pass `out_children = NULL` to query the count only.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_BUFFER_TOO_SMALL` if `capacity` is too small; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not differentiable.

### qaws_curve_local_support

```c
qaws_status qaws_curve_local_support(
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int order,
	qaws_local_support* out_support);
```

Returns the parameters influencing the evaluation at `t`. For `QAWS_CAP_LINEAR` curves the weights are exact basis derivatives for orders `0..order`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_OUT_OF_RANGE` if `order > QAWS_DIFF_MAX_ORDER`; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the family has no local support (no `QAWS_CAP_LOCAL_SUPPORT`).

### qaws_curve_build_support_index

```c
qaws_status qaws_curve_build_support_index(
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count);
```

Builds a support index for gather-style adjoint accumulation (CPU or GPU): for every element `e` of `field`, the samples influenced by it are `out_samples[out_offsets[e] .. out_offsets[e+1])`, in sample order. `out_offsets` needs `element_count + 1` entries.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_BUFFER_TOO_SMALL` when `offset_capacity` is below `element_count + 1`, or when `out_samples` is NULL or `sample_capacity` is too small (then `*out_entry_count` holds the required size); `QAWS_STATUS_UNSUPPORTED_OPERATION` if the family has no local support.

### qaws_curve_eval_batch_tangent_2d / 3d

```c
qaws_status qaws_curve_eval_batch_tangent_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent);

qaws_status qaws_curve_eval_batch_tangent_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent);
```

Evaluates `count` samples and their tangent along a combined coordinate and parameter direction.

**Parameters:**
- `t` -- `count` evaluation parameters.
- `t_tangent` -- `count` coordinate tangents, or NULL (coordinate held fixed).
- `channels` -- `QAWS_EVAL_FLAG_POSITION/D1/D2/D3` bits to compute.
- `param_tangent` -- parameter tangent views, or NULL (parameters held fixed).
- `out_primal` -- `count` primal jets, or NULL.
- `out_tangent` -- `count` tangent jets (required).

Outputs are overwritten. The curve dimension must match the function suffix.

**Report:** for each sample with a non-zero coordinate tangent, the curve's class is recorded; `QAWS_DIFF_AT_BOUNDARY` is recorded when the sample lies exactly on a span boundary, and `QAWS_FREEZE_SPAN` for piecewise-smooth curves.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `curve`, `t` or `out_tangent`, or a view with the wrong component count; `QAWS_STATUS_INVALID_DIMENSION` on dimension mismatch; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not differentiable.

### qaws_curve_eval_batch_tangent2_2d / 3d

```c
qaws_status qaws_curve_eval_batch_tangent2_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent,
	qaws_curve_jet_2d* out_tangent2);

qaws_status qaws_curve_eval_batch_tangent2_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent,
	qaws_curve_jet_3d* out_tangent2);
```

Same as the tangent functions, and additionally writes the second directional derivative along the same direction into `out_tangent2` (required).

**Returns:** as above; additionally `QAWS_STATUS_INVALID_ARGUMENT` if `out_tangent2` is NULL and `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve lacks `QAWS_CAP_TANGENT2`.

### qaws_curve_eval_batch_adjoint_2d / 3d

```c
qaws_status qaws_curve_eval_batch_adjoint_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

qaws_status qaws_curve_eval_batch_adjoint_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);
```

Pulls back `count` jet adjoints to parameter and coordinate adjoints.

**Parameters:**
- `out_adjoint` -- `count` jet adjoints. Only the channels in `channels` of each jet are read.
- `param_adjoint` -- parameter adjoint views, accumulated (`+=`), or NULL.
- `t_adjoint` -- `count` coordinate adjoints, accumulated (`+=`), or NULL.

Parameter accumulation follows `ctx->accumulation`. Knot adjoints are produced when a `QAWS_FIELD_KNOTS` view is present and the family supports knot derivatives.

**Report:** as for the tangent functions, with the coordinate adjoint playing the role of the coordinate tangent.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `curve`, `t` or `out_adjoint`, or a view with the wrong component count; `QAWS_STATUS_INVALID_DIMENSION`; `QAWS_STATUS_UNSUPPORTED_OPERATION`; `QAWS_STATUS_ALLOCATION_FAILURE` if scratch allocation fails.

### qaws_curve_eval_tangent_2d / 3d

```c
qaws_status qaws_curve_eval_tangent_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_scalar t_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal,
	qaws_curve_jet_2d* out_tangent);

qaws_status qaws_curve_eval_tangent_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_scalar t_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal,
	qaws_curve_jet_3d* out_tangent);
```

Single-sample wrapper around the batch tangent functions. Pass `t_tangent = 0` to hold the coordinate fixed.

### qaws_curve_eval_adjoint_2d / 3d

```c
qaws_status qaws_curve_eval_adjoint_2d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);

qaws_status qaws_curve_eval_adjoint_3d(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* t_adjoint);
```

Single-sample wrapper around the batch adjoint functions. Accumulates (`+=`) into `param_adjoint` and `*t_adjoint` (either may be NULL).

---

## Surfaces

Declared in `qaws_diff.h`. Channels are `qaws_surface_jet_channel` bits (all partials up to total order 3). `(u, v)` are clamped to the surface domain like `qaws_surface_evaluate`. Surface fields are always 3D.

### qaws_surface_get_diff_capabilities

```c
unsigned int qaws_surface_get_diff_capabilities(qaws_surface const* surface);
```

Returns the `qaws_diff_capability` bits of the surface family, or 0 if not differentiable.

### qaws_surface_get_diff_class

```c
qaws_diff_class qaws_surface_get_diff_class(qaws_surface const* surface);
```

Returns the differentiability class of the surface, or `QAWS_DIFF_UNSUPPORTED`.

### qaws_surface_describe_fields

```c
qaws_status qaws_surface_describe_fields(
	qaws_surface const* surface,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count);
```

Same contract as `qaws_curve_describe_fields`.

### qaws_surface_read_field

```c
qaws_status qaws_surface_read_field(
	qaws_surface const* surface,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count);
```

Same contract as `qaws_curve_read_field`.

### qaws_surface_diff_children

```c
qaws_status qaws_surface_diff_children(
	qaws_surface const* surface,
	qaws_diff_child* out_children,
	unsigned int capacity,
	unsigned int* out_count);
```

Lists the objects the surface rules chain into; child `i` matches `views->children[i]`. Same contract as `qaws_curve_diff_children`.

### qaws_surface_local_support

```c
qaws_status qaws_surface_local_support(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int order,
	qaws_surface_support* out_support);
```

Returns the tensor-product support of the evaluation at `(u, v)` (clamped to the domain), with weights for orders `0..order` in each direction.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_OUT_OF_RANGE` if `order > QAWS_DIFF_MAX_ORDER`; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the family has no local support.

### qaws_surface_build_support_index

```c
qaws_status qaws_surface_build_support_index(
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count);
```

Surface counterpart of `qaws_curve_build_support_index`, with samples given as `(u[i], v[i])`. Same buffer and status contract.

### qaws_surface_eval_jet

```c
qaws_status qaws_surface_eval_jet(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int channels,
	qaws_surface_jet* out_jet);
```

Evaluates the primal spatial jet analytically (no finite differences). Only the requested channels are filled; `out_jet->channels` is set to `channels` masked to `QAWS_SJET_ORDER3`.

### qaws_surface_eval_batch_tangent

```c
qaws_status qaws_surface_eval_batch_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	qaws_scalar const* u_tangent,
	qaws_scalar const* v_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent);
```

Evaluates `count` samples and their tangent along a combined coordinate and parameter direction.

**Parameters:**
- `u`, `v` -- `count` coordinates each.
- `u_tangent`, `v_tangent` -- coordinate tangents, or NULL (held fixed).
- `param_tangent` -- parameter tangent views, or NULL.
- `out_primal` -- primal jets, or NULL.
- `out_tangent` -- tangent jets (required). Overwritten.

**Report:** for each sample with a coordinate tangent, the surface's class is recorded; `QAWS_DIFF_AT_BOUNDARY` when `(u, v)` lies exactly on an interior knot line, and `QAWS_FREEZE_SPAN` for piecewise-smooth surfaces.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `surface`, `u`, `v` or `out_tangent`, or a view with the wrong component count; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the surface is not differentiable.

### qaws_surface_eval_batch_tangent2

```c
qaws_status qaws_surface_eval_batch_tangent2(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	qaws_scalar const* u_tangent,
	qaws_scalar const* v_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent,
	qaws_surface_jet* out_tangent2);
```

Same as `qaws_surface_eval_batch_tangent`, and additionally writes the second directional derivative into `out_tangent2` (required).

**Returns:** as above; additionally `QAWS_STATUS_INVALID_ARGUMENT` if `out_tangent2` is NULL and `QAWS_STATUS_UNSUPPORTED_OPERATION` if the surface lacks `QAWS_CAP_TANGENT2`.

### qaws_surface_eval_batch_adjoint

```c
qaws_status qaws_surface_eval_batch_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	unsigned int channels,
	qaws_surface_jet const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint,
	qaws_scalar* v_adjoint);
```

Pulls back `count` jet adjoints. Only the channels in `channels` of each jet are read. Accumulates (`+=`) into `param_adjoint`, `u_adjoint` and `v_adjoint` (each may be NULL). Parameter accumulation follows `ctx->accumulation`. Knot adjoints are produced when `QAWS_FIELD_U_KNOTS` / `QAWS_FIELD_V_KNOTS` views are present and the family supports them.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `surface`, `u`, `v` or `out_adjoint`, or a view with the wrong component count; `QAWS_STATUS_UNSUPPORTED_OPERATION`; `QAWS_STATUS_ALLOCATION_FAILURE`.

### qaws_surface_eval_tangent

```c
qaws_status qaws_surface_eval_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar u_tangent,
	qaws_scalar v_tangent,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent);
```

Single-sample wrapper around `qaws_surface_eval_batch_tangent`.

### qaws_surface_eval_adjoint

```c
qaws_status qaws_surface_eval_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int channels,
	qaws_surface_jet const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint,
	qaws_scalar* v_adjoint);
```

Single-sample wrapper around `qaws_surface_eval_batch_adjoint`. Accumulates (`+=`).

---

## Geometry

Declared in `qaws_diff_geometry.h`. Local differential geometry as closed-form functions of a spatial jet, so the rules compose with the curve and surface rules:

```
parameters --eval_tangent--> jet tangent --geometry_eval--> quantity tangent
quantity adjoint --geometry_adjoint--> jet adjoint --eval_adjoint--> parameters
```

Validity (written to `*out_validity`):

- `QAWS_DIFF_VALID` -- all quantities are differentiable here.
- `QAWS_DIFF_ILL_CONDITIONED` -- a quantity is undefined or not differentiable (curve: zero curvature leaves normal, binormal and torsion undefined; surface: umbilic point for principal curvatures). Those quantities and their derivatives are returned as zero.
- `QAWS_DIFF_INVALID` -- singular jet (zero speed, degenerate normal): everything is returned as zero.

### qaws_curve_geometry_2d

```c
typedef struct qaws_curve_geometry_2d
{
	qaws_vec2 tangent;       /* unit tangent */
	qaws_vec2 normal;        /* tangent rotated by +90 degrees */
	qaws_scalar speed;       /* |C'| */
	qaws_scalar curvature;   /* signed: (x'y'' - y'x'') / |C'|^3 */
} qaws_curve_geometry_2d;
```

### qaws_curve_geometry_3d

```c
typedef struct qaws_curve_geometry_3d
{
	qaws_vec3 tangent;       /* T */
	qaws_vec3 normal;        /* N = B x T */
	qaws_vec3 binormal;      /* B = (C' x C'') / |C' x C''| */
	qaws_scalar speed;
	qaws_scalar curvature;   /* |C' x C''| / |C'|^3 */
	qaws_scalar torsion;     /* (C' x C'') . C''' / |C' x C''|^2 */
} qaws_curve_geometry_3d;
```

Frenet frame, speed, curvature and torsion.

### qaws_surface_geometry

```c
typedef struct qaws_surface_geometry
{
	qaws_vec3 normal;        /* (Su x Sv) / |Su x Sv| */
	qaws_scalar E, F, G;     /* first fundamental form */
	qaws_scalar L, M, N;     /* second fundamental form */
	qaws_scalar gaussian;    /* (LN - M^2) / (EG - F^2) */
	qaws_scalar mean;        /* (EN + GL - 2FM) / (2 (EG - F^2)) */
	qaws_scalar kappa1;      /* H + sqrt(H^2 - K) */
	qaws_scalar kappa2;      /* H - sqrt(H^2 - K) */
} qaws_surface_geometry;
```

Unit normal, fundamental forms, Gaussian and mean curvature, principal curvatures.

### qaws_curve_geometry_eval_2d / 3d

```c
qaws_status qaws_curve_geometry_eval_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_jet_2d const* tangent,
	qaws_curve_jet_2d const* tangent2,
	qaws_curve_geometry_2d* out_value,
	qaws_curve_geometry_2d* out_tangent,
	qaws_curve_geometry_2d* out_tangent2,
	qaws_diff_validity* out_validity);

qaws_status qaws_curve_geometry_eval_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_jet_3d const* tangent,
	qaws_curve_jet_3d const* tangent2,
	qaws_curve_geometry_3d* out_value,
	qaws_curve_geometry_3d* out_tangent,
	qaws_curve_geometry_3d* out_tangent2,
	qaws_diff_validity* out_validity);
```

Forward evaluation of the curve quantities and, optionally, their tangent and second tangent. Required primal channels: D1, D2, and D3 for torsion in 3D. `tangent`, `tangent2`, `out_tangent`, `out_tangent2` and `out_validity` may be NULL; a tangent output needs a tangent input, a tangent2 output needs both. Outputs are overwritten.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `primal` or `out_value` is NULL or a requested output lacks its inputs.

### qaws_surface_geometry_eval

```c
qaws_status qaws_surface_geometry_eval(
	qaws_surface_jet const* primal,
	qaws_surface_jet const* tangent,
	qaws_surface_jet const* tangent2,
	qaws_surface_geometry* out_value,
	qaws_surface_geometry* out_tangent,
	qaws_surface_geometry* out_tangent2,
	qaws_diff_validity* out_validity);
```

Surface counterpart. Required primal channels: U, V, UU, UV, VV. Same optional-argument and status contract as the curve functions.

### qaws_curve_geometry_adjoint_2d / 3d

```c
qaws_status qaws_curve_geometry_adjoint_2d(
	qaws_curve_jet_2d const* primal,
	qaws_curve_geometry_2d const* adjoint,
	qaws_curve_jet_2d* inout_jet_adjoint,
	qaws_diff_validity* out_validity);

qaws_status qaws_curve_geometry_adjoint_3d(
	qaws_curve_jet_3d const* primal,
	qaws_curve_geometry_3d const* adjoint,
	qaws_curve_jet_3d* inout_jet_adjoint,
	qaws_diff_validity* out_validity);
```

Reverse propagation: accumulates (`+=`) into `inout_jet_adjoint` the pullback of the quantity adjoint. Only the jet channels the quantities depend on are touched (D1, D2, plus D3 in 3D), and `inout_jet_adjoint->channels` is extended with them. Feed the result to `qaws_curve_eval_adjoint_*`. `out_validity` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `primal`, `adjoint` or `inout_jet_adjoint` is NULL.

### qaws_surface_geometry_adjoint

```c
qaws_status qaws_surface_geometry_adjoint(
	qaws_surface_jet const* primal,
	qaws_surface_geometry const* adjoint,
	qaws_surface_jet* inout_jet_adjoint,
	qaws_diff_validity* out_validity);
```

Surface counterpart. Accumulates (`+=`) into channels U, V, UU, UV, VV and adds them to `inout_jet_adjoint->channels`.

---

## Implicit operations

Declared in `qaws_diff_ops.h`. Operations defined by a solve. The solution `x*` of `g(x, theta, q) = 0` is differentiated with the implicit function theorem, never by differentiating solver iterations:

```
d(x*) / d(theta) = -(dg / dx)^-1 (dg / dtheta)
```

Every call solves the primal problem itself, then either writes the tangent of the solution (outputs overwritten) or accumulates (`+=`) the pullback of an output adjoint into parameter and input adjoints. Views passed as NULL are inactive.

Every call reports through `ctx->report` (when set):

- `condition_number` -- conditioning of `dg/dx` relative to the metric.
- `residual` -- relative `|g|` at the returned solution.
- `branch_gap` -- distance gap to the best competing local solution.
- `validity` -- `QAWS_DIFF_ILL_CONDITIONED` near singular `dg/dx`; `QAWS_DIFF_AMBIGUOUS` when a competing solution is within 0.5% of the distance; `QAWS_DIFF_VALID_LOCALLY` (with `QAWS_FREEZE_ACTIVE_SET`) when the solution sits on the domain boundary and the boundary coordinate is held fixed.

When `dg/dx` is singular, coordinate tangents and the implicit part of adjoints are returned as zero.

### Closest point on a curve

Optimality: `C'(t) . (C(t) - q) = 0`. Works for 2D curves too; the `z` of the query and of all outputs is then zero.

#### qaws_curve_closest_point

```c
typedef struct qaws_curve_closest_point
{
	qaws_scalar t;
	qaws_vec3 position;
	qaws_scalar distance;
} qaws_curve_closest_point;
```

Solution (or its tangent / adjoint): curve parameter, point on the curve, distance to the query.

#### qaws_curve_closest_point_tangent

```c
qaws_status qaws_curve_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_closest_point* out_value,
	qaws_curve_closest_point* out_tangent);
```

Solves the closest point to `query` and returns its tangent for a query tangent and/or a parameter tangent (either may be NULL). `out_value` is required; `out_tangent` may be NULL to solve only.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `curve` or `out_value` is NULL.

#### qaws_curve_closest_point_adjoint

```c
qaws_status qaws_curve_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 query,
	qaws_curve_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint);
```

Solves the closest point and accumulates (`+=`) the pullback of an output adjoint (`t`, `position`, `distance`) into `param_adjoint` and `query_adjoint` (either may be NULL).

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `curve` or `adjoint` is NULL.

### Closest point on a surface

Optimality: `Su . (S - q) = 0` and `Sv . (S - q) = 0`.

#### qaws_surface_closest_point

```c
typedef struct qaws_surface_closest_point
{
	qaws_scalar u, v;
	qaws_vec3 position;
	qaws_scalar distance;
} qaws_surface_closest_point;
```

#### qaws_surface_closest_point_tangent

```c
qaws_status qaws_surface_closest_point_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_vec3 const* query_tangent,
	qaws_diff_views const* param_tangent,
	qaws_surface_closest_point* out_value,
	qaws_surface_closest_point* out_tangent);
```

Surface counterpart of `qaws_curve_closest_point_tangent`.

#### qaws_surface_closest_point_adjoint

```c
qaws_status qaws_surface_closest_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_vec3 query,
	qaws_surface_closest_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* query_adjoint);
```

Surface counterpart of `qaws_curve_closest_point_adjoint`. Accumulates (`+=`).

### Roots seeded by the discrete finders

The operations below differentiate one solution. The caller passes a seed, typically from `qaws_curve_find_intersections_*`, `qaws_curve_find_plane_intersections`, `qaws_surface_find_curve_intersections`, `qaws_curve_find_extrema` or `qaws_curve_find_inflection_points`. Seeds are clamped to the parameter range, refined by Newton and differentiated by the implicit function theorem. Which solution exists, and how many, is frozen (`QAWS_FREEZE_ACTIVE_SET`). Tangential configurations report `QAWS_DIFF_ILL_CONDITIONED` and return zero coordinate tangents.

#### qaws_curve_pair_point

```c
typedef struct qaws_curve_pair_point
{
	qaws_scalar t_a, t_b;
	qaws_vec3 position_a, position_b;
	qaws_scalar distance;
} qaws_curve_pair_point;
```

Closest approach of two curves (an intersection when `distance = 0`). Stationarity of `|A(ta) - B(tb)|^2 / 2`: `A' . (A - B) = 0`, `B' . (A - B) = 0`. A parameter on its range end with the gradient pointing outward is held fixed.

#### qaws_curve_pair_point_tangent

```c
qaws_status qaws_curve_pair_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_scalar t_a_seed,
	qaws_scalar t_b_seed,
	qaws_diff_views const* tangent_a,
	qaws_diff_views const* tangent_b,
	qaws_curve_pair_point* out_value,
	qaws_curve_pair_point* out_tangent);
```

Refines the pair point from `(t_a_seed, t_b_seed)` and returns its tangent. `tangent_a` and `tangent_b` select the parameter tangents of A and of B (either may be NULL). `out_tangent` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if a curve or `out_value` is NULL.

#### qaws_curve_pair_point_adjoint

```c
qaws_status qaws_curve_pair_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_scalar t_a_seed,
	qaws_scalar t_b_seed,
	qaws_curve_pair_point const* adjoint,
	qaws_diff_views* adjoint_a,
	qaws_diff_views* adjoint_b);
```

Refines the pair point and accumulates (`+=`) the pullback of `adjoint` into the parameter adjoints of A and B (either may be NULL).

#### qaws_curve_plane_point

```c
typedef struct qaws_curve_plane_point
{
	qaws_scalar t;
	qaws_vec3 position;
} qaws_curve_plane_point;
```

Curve crossing a plane: `n . (C(t) - p) = 0`. The plane (`qaws_plane`, from `qaws_inspect.h`) is an input with its own tangent and adjoint (point and normal).

#### qaws_curve_plane_point_tangent

```c
qaws_status qaws_curve_plane_point_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_plane const* plane,
	qaws_scalar t_seed,
	qaws_plane const* plane_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_plane_point* out_value,
	qaws_curve_plane_point* out_tangent);
```

Refines the crossing from `t_seed` and returns its tangent for a plane tangent and/or a parameter tangent (either may be NULL). `out_tangent` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `curve`, `plane` or `out_value` is NULL.

#### qaws_curve_plane_point_adjoint

```c
qaws_status qaws_curve_plane_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_plane const* plane,
	qaws_scalar t_seed,
	qaws_curve_plane_point const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_plane* plane_adjoint);
```

Accumulates (`+=`) the pullback into `param_adjoint` and into the `point` and `normal` of `plane_adjoint` (either may be NULL).

#### qaws_surface_curve_point

```c
typedef struct qaws_surface_curve_point
{
	qaws_scalar u, v, t;
	qaws_vec3 position;
} qaws_surface_curve_point;
```

Curve piercing a surface: `S(u, v) - C(t) = 0`.

#### qaws_surface_curve_point_tangent

```c
qaws_status qaws_surface_curve_point_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_curve const* curve,
	qaws_scalar u_seed,
	qaws_scalar v_seed,
	qaws_scalar t_seed,
	qaws_diff_views const* surface_tangent,
	qaws_diff_views const* curve_tangent,
	qaws_surface_curve_point* out_value,
	qaws_surface_curve_point* out_tangent);
```

Refines the piercing point from `(u_seed, v_seed, t_seed)` and returns its tangent for surface and/or curve parameter tangents (either may be NULL). `out_tangent` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `surface`, `curve` or `out_value`; `QAWS_STATUS_INVALID_DIMENSION` if the curve is not 3D.

#### qaws_surface_curve_point_adjoint

```c
qaws_status qaws_surface_curve_point_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_curve const* curve,
	qaws_scalar u_seed,
	qaws_scalar v_seed,
	qaws_scalar t_seed,
	qaws_surface_curve_point const* adjoint,
	qaws_diff_views* surface_adjoint,
	qaws_diff_views* curve_adjoint);
```

Accumulates (`+=`) the pullback into the surface and curve parameter adjoints (either may be NULL). Same statuses as the tangent function.

#### qaws_curve_extremum

```c
typedef struct qaws_curve_extremum
{
	qaws_scalar t;
	qaws_vec3 position;
	qaws_scalar value;
} qaws_curve_extremum;
```

Extremum of the height `e . C(t)` along a direction: `e . C'(t) = 0`, `value = e . C(t)`. The direction is an input with its own tangent and adjoint. For 2D curves the direction is projected to the xy plane.

#### qaws_curve_extremum_tangent

```c
qaws_status qaws_curve_extremum_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 direction,
	qaws_scalar t_seed,
	qaws_vec3 const* direction_tangent,
	qaws_diff_views const* param_tangent,
	qaws_curve_extremum* out_value,
	qaws_curve_extremum* out_tangent);
```

Refines the extremum from `t_seed` and returns its tangent for a direction tangent and/or a parameter tangent (either may be NULL). `out_tangent` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `curve` or `out_value` is NULL.

#### qaws_curve_extremum_adjoint

```c
qaws_status qaws_curve_extremum_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_vec3 direction,
	qaws_scalar t_seed,
	qaws_curve_extremum const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_vec3* direction_adjoint);
```

Accumulates (`+=`) the pullback into `param_adjoint` and `direction_adjoint` (either may be NULL).

#### qaws_curve_inflection

```c
typedef struct qaws_curve_inflection
{
	qaws_scalar t;
	qaws_vec3 position;
} qaws_curve_inflection;
```

Inflection of the xy projection: `x'y'' - y'x'' = 0`.

#### qaws_curve_inflection_tangent

```c
qaws_status qaws_curve_inflection_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t_seed,
	qaws_diff_views const* param_tangent,
	qaws_curve_inflection* out_value,
	qaws_curve_inflection* out_tangent);
```

Refines the inflection from `t_seed` and returns its tangent for a parameter tangent. `out_tangent` may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` if `curve` or `out_value` is NULL.

#### qaws_curve_inflection_adjoint

```c
qaws_status qaws_curve_inflection_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t_seed,
	qaws_curve_inflection const* adjoint,
	qaws_diff_views* param_adjoint);
```

Accumulates (`+=`) the pullback into `param_adjoint` (may be NULL).

---

## Differential maps

Declared in `qaws_diff_map.h`. An operation that builds new geometry (split, join, conversion, degree change, fitting) returns, next to its output objects, a map from the parameters of its inputs to the parameters of its outputs:

- `qaws_diff_map_tangent`: `out = J in` (outputs are overwritten).
- `qaws_diff_map_adjoint`: `in += J^T out`.

Inputs and outputs are numbered objects, each addressed through its own `qaws_diff_views` (NULL entries are inactive). Each operation documents its numbering. Maps are either materialized sparse Jacobians (exact, evaluated at the operation's inputs) or operators that apply `J` without forming it.

### qaws_diff_map / qaws_diff_map_kind

```c
typedef struct qaws_diff_map qaws_diff_map;

typedef enum qaws_diff_map_kind
{
	QAWS_DIFF_MAP_LINEAR_SPARSE = 0,
	QAWS_DIFF_MAP_OPERATOR
} qaws_diff_map_kind;
```

Opaque map handle and its representation. All operations in this header currently return `QAWS_DIFF_MAP_LINEAR_SPARSE` maps.

### qaws_diff_map_entry

```c
#define QAWS_DIFF_MAP_ALL_COMPONENTS 0xFFu

typedef struct qaws_diff_map_entry
{
	unsigned int out_object, out_element;
	unsigned int in_object, in_element;
	unsigned char out_field, in_field;         /* qaws_diff_field */
	unsigned char out_component, in_component; /* or QAWS_DIFF_MAP_ALL_COMPONENTS */
	qaws_scalar weight;
} qaws_diff_map_entry;
```

One sparse Jacobian entry: `out[out_object].out_field[out_element][out_component] += weight * in[in_object].in_field[in_element][in_component]`. When both components are `QAWS_DIFF_MAP_ALL_COMPONENTS`, the same weight applies to every component.

### qaws_diff_map_get_kind

```c
qaws_diff_map_kind qaws_diff_map_get_kind(qaws_diff_map const* map);
```

Returns the map representation (`QAWS_DIFF_MAP_LINEAR_SPARSE` for NULL).

### qaws_diff_map_input_count / qaws_diff_map_output_count

```c
unsigned int qaws_diff_map_input_count(qaws_diff_map const* map);
unsigned int qaws_diff_map_output_count(qaws_diff_map const* map);
```

Return the number of input and output objects (0 for NULL).

### qaws_diff_map_tangent

```c
qaws_status qaws_diff_map_tangent(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* in_tangents,
	unsigned int in_count,
	qaws_diff_views* const* out_tangents,
	unsigned int out_count);
```

Applies `out = J in`. Every non-NULL output view set is cleared first, then filled. Entries whose input or output object is beyond the given counts, NULL, or missing the field are skipped.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `map` or NULL arrays with non-zero counts; `QAWS_STATUS_UNSUPPORTED_OPERATION` when a tangent is given for an input the map does not differentiate (for example the split parameter of a non-Bezier curve).

### qaws_diff_map_adjoint

```c
qaws_status qaws_diff_map_adjoint(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* out_adjoints,
	unsigned int out_count,
	qaws_diff_views* const* in_adjoints,
	unsigned int in_count);
```

Applies `in += J^T out`. Input adjoints are accumulated, not cleared.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `map` or NULL arrays with non-zero counts.

### qaws_diff_map_get_entries

```c
qaws_status qaws_diff_map_get_entries(
	qaws_diff_map const* map,
	qaws_diff_map_entry* out_entries,
	unsigned int capacity,
	unsigned int* out_count);
```

Copies the sparse entries of a `LINEAR_SPARSE` map. `*out_count` receives the entry count. Pass `out_entries = NULL` to query the count only.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_BUFFER_TOO_SMALL` if `capacity` is too small; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `map` or `out_count`.

### qaws_diff_map_destroy

```c
void qaws_diff_map_destroy(qaws_diff_map* map);
```

Frees a map. NULL is ignored.

### Operations returning maps

Polynomial families are exact linear maps of their fields; NURBS inputs are linear in homogeneous coordinates and their maps include the projection. Knots of the inputs are held fixed. The split parameter is differentiable for Bezier curves; for other families a split parameter tangent is refused.

All operations accept `out_map = NULL`, in which case only the geometry is built. The caller owns the returned curves and map. Each returns `QAWS_STATUS_INVALID_ARGUMENT` on NULL inputs or output curve pointers, the status of the underlying geometry operation on failure, and `QAWS_STATUS_ALLOCATION_FAILURE` if the map cannot be allocated.

#### qaws_curve_split_diff

```c
qaws_status qaws_curve_split_diff(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_curve** out_left,
	qaws_curve** out_right,
	qaws_diff_map** out_map);
```

Splits `curve` at `parameter`. Inputs: 0 = curve, 1 = split parameter (`QAWS_FIELD_PARAMETER`, element 0). Outputs: 0 = left, 1 = right.

#### qaws_curve_join_diff

```c
qaws_status qaws_curve_join_diff(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_curve** out_joined,
	qaws_diff_map** out_map);
```

Joins two curves. Inputs: 0 = `curve_a`, 1 = `curve_b`. Outputs: 0 = joined.

#### qaws_curve_convert_hermite_to_bezier_diff

```c
qaws_status qaws_curve_convert_hermite_to_bezier_diff(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_curve** out_bezier,
	qaws_diff_map** out_map);
```

Converts one Hermite span to a Bezier curve. Inputs: 0 = curve. Outputs: 0 = converted curve.

#### qaws_curve_convert_bezier_to_bspline_diff

```c
qaws_status qaws_curve_convert_bezier_to_bspline_diff(
	qaws_curve const* curve,
	qaws_curve** out_bspline,
	qaws_diff_map** out_map);
```

Inputs: 0 = curve. Outputs: 0 = converted curve.

#### qaws_curve_convert_bspline_to_nurbs_diff

```c
qaws_status qaws_curve_convert_bspline_to_nurbs_diff(
	qaws_curve const* curve,
	qaws_curve** out_nurbs,
	qaws_diff_map** out_map);
```

Inputs: 0 = curve. Outputs: 0 = converted curve.

#### qaws_curve_elevate_degree_diff

```c
qaws_status qaws_curve_elevate_degree_diff(
	qaws_curve const* curve,
	qaws_curve** out_elevated,
	qaws_diff_map** out_map);
```

Inputs: 0 = curve. Outputs: 0 = elevated curve.

#### qaws_curve_reduce_degree_diff

```c
qaws_status qaws_curve_reduce_degree_diff(
	qaws_curve const* curve,
	qaws_curve** out_reduced,
	qaws_diff_map** out_map);
```

Inputs: 0 = curve. Outputs: 0 = reduced curve.

#### qaws_curve_fit_bspline_diff

```c
struct qaws_bspline_fit_desc;

qaws_status qaws_curve_fit_bspline_diff(
	struct qaws_bspline_fit_desc const* desc,
	qaws_curve** out_curve,
	qaws_diff_map** out_map);
```

Least-squares B-spline fit (`qaws_curve_fit_bspline`, descriptor in `qaws_export.h`) with the exact Jacobian of its normal equations.

- Inputs: 0 = data points (`QAWS_FIELD_POINTS`, element = data point); 1 = sample parameters (`QAWS_FIELD_PARAMETER`, element = data point), present only when `desc->parameters` is given.
- Outputs: 0 = fitted curve (`QAWS_FIELD_CONTROL_POINTS`, `QAWS_FIELD_KNOTS`).

The derivative follows every smooth stage of the fit: chord-length parameters (when `desc->parameters` is NULL they are functions of the data points), averaged knots, basis functions and the solve. Knot span membership of each sample is frozen (`QAWS_FREEZE_SPAN`).

**Returns:** `QAWS_STATUS_OK`; the status of `qaws_curve_fit_bspline`; `QAWS_STATUS_UNSUPPORTED_OPERATION` if `degree + 1 > QAWS_DIFF_MAX_SUPPORT`; `QAWS_STATUS_NUMERICAL_FAILURE` if the normal equations are not positive definite; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `desc` or `out_curve`.

---

### qaws_curve_offset_3d_diff

```c
qaws_status qaws_curve_offset_3d_diff(
	qaws_curve const* curve,
	qaws_scalar distance,
	int direction_mode,
	qaws_vec3 const* direction,
	qaws_surface const* surface,
	unsigned int sample_count,
	qaws_curve** out_curve,
	qaws_diff_map** out_map);
```

Runs `qaws_curve_offset_3d` and returns its exact map. The offset samples `C(t_i) + distance * n(t_i)` at fixed parameters (the sample set is frozen, `QAWS_FREEZE_SAMPLE_SET`) and fits a B-spline through them; the fit is linear in the samples, so each map column is the fit of the sample tangents.

- inputs: 0 = curve (every differentiable field, including knots), 1 = distance (`QAWS_FIELD_PARAMETER`, element 0), 2 = direction (`QAWS_FIELD_DIRECTION`, element 0, three components; `direction_mode` 0 only, otherwise frozen), 3 = surface (every differentiable field but knots; `direction_mode` 2 only).
- outputs: 0 = offset curve (`QAWS_FIELD_CONTROL_POINTS`).

`direction_mode` 0 (constant direction, normalized), 1 (Frenet normal, differentiated through `qaws_curve_geometry_eval_3d`) and 2 (normal of `surface` at the closest point of each sample) are supported. In mode 2 the closest point `(u, v)` is differentiated through its optimality conditions `S_u . (S - q) = S_v . (S - q) = 0` at the points the offset itself finds, so `(u', v') = -H^-1 g'`, and the normal through `normalize(S_u x S_v)` (oriented like the evaluated normal); the curve moves the query `q`, the surface fields move `S`, `S_u`, `S_v`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `curve` or `out_curve`, or a NULL `surface` in mode 2; `QAWS_STATUS_UNSUPPORTED_OPERATION` for an unknown mode, in mode 1 for curves with straight pieces where the Frenet normal is undefined, and in mode 2 where the closest-point Hessian is singular; errors of `qaws_curve_offset_3d`. On failure `*out_curve` is NULL.

---

## Integral functionals

Declared in `qaws_diff_functionals.h`. Integral functionals of curves and surfaces with exact derivatives of their quadrature:

- value, tangent, tangent2 along one parameter direction;
- gradient, accumulated (`+=`) into parameter views;
- Hessian-vector product (`+=`), direct for families linear in their fields (`QAWS_CAP_LINEAR`); other families return `QAWS_STATUS_UNSUPPORTED_OPERATION` and the caller composes Hessian-vector products from tangents and adjoints.

Curves integrate span by span with a composite Gauss-Legendre rule: every span is split into 8 pieces of `quadrature` points (0 = 6; clamped to 2..8), so the length (a square root of the speed) is exact to about 1e-12. Surfaces integrate a grid of cells with Gauss points (`quadrature` = cells per direction, 0 = 8; 4 x 4 points per cell).

### qaws_curve_functional

```c
typedef enum qaws_curve_functional
{
	QAWS_FUNCTIONAL_LENGTH = 0,            /* integral of |C'| dt */
	QAWS_FUNCTIONAL_BENDING,               /* integral of |C''|^2 dt */
	QAWS_FUNCTIONAL_CURVATURE_SQUARED      /* integral of kappa^2 ds */
} qaws_curve_functional;
```

### qaws_surface_functional

```c
typedef enum qaws_surface_functional
{
	QAWS_FUNCTIONAL_AREA = 0,              /* integral of |Su x Sv| du dv */
	QAWS_FUNCTIONAL_THIN_PLATE,            /* integral of |Suu|^2 + 2|Suv|^2 + |Svv|^2 */
	QAWS_FUNCTIONAL_WILLMORE               /* integral of H^2 dA */
} qaws_surface_functional;
```

### qaws_curve_functional_eval

```c
qaws_status qaws_curve_functional_eval(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_scalar* out_value,
	qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2);
```

Evaluates the functional and its first and second derivatives along the parameter direction `direction` (NULL = zero direction). Each output may be NULL; outputs are overwritten. A `QAWS_FIELD_KNOTS` direction is differentiated when the curve's knot field has `QAWS_CAP_TANGENT`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `curve` or an unknown functional; errors from the underlying curve evaluation.

### qaws_curve_functional_gradient

```c
qaws_status qaws_curve_functional_gradient(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views* gradient,
	qaws_scalar* out_value);
```

Accumulates (`+=`) the gradient of the functional into `gradient` and writes the value to `*out_value` (may be NULL). Knot gradients are produced under the same condition as `qaws_curve_functional_eval`.

### qaws_curve_functional_hvp

```c
qaws_status qaws_curve_functional_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_curve_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);
```

Accumulates (`+=`) the Hessian-vector product `H * direction` into `out_hv`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not `QAWS_CAP_LINEAR`, or if `direction` or `out_hv` contains a knot view (knots enter non-linearly).

### qaws_surface_functional_eval

```c
qaws_status qaws_surface_functional_eval(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_scalar* out_value,
	qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2);
```

Surface counterpart of `qaws_curve_functional_eval`. The quadrature runs on a cell grid over the domain; `U_KNOTS` / `V_KNOTS` directions change the integrand at fixed nodes, and the domain end knots also move every node and scale the weights.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL `surface` or an unknown functional.

### qaws_surface_functional_gradient

```c
qaws_status qaws_surface_functional_gradient(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views* gradient,
	qaws_scalar* out_value);
```

Accumulates (`+=`) the gradient into `gradient` and writes the value to `*out_value` (may be NULL). Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` if `gradient` contains a knot view.

### qaws_surface_functional_hvp

```c
qaws_status qaws_surface_functional_hvp(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_surface_functional functional,
	unsigned int quadrature,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);
```

Accumulates (`+=`) the Hessian-vector product into `out_hv`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_UNSUPPORTED_OPERATION` if the surface is not `QAWS_CAP_LINEAR` or a knot view is present in `direction` or `out_hv`.

---

## Inverse-CDF sampling

`qaws_diff_sampling.h`: samples at prescribed values of a measure along the curve, with exact first and second order derivatives, forward and backward. The measure `M(t) = int m dt`, `m = rho |C'|`, acts as a CDF:

| kind | `rho` | use |
|---|---|---|
| `QAWS_MEASURE_ARC_LENGTH` | 1 | constant-speed sampling |
| `QAWS_MEASURE_CURVATURE` | `sqrt(floor^2 + kappa^2)` | samples gather where the curve bends; `floor > 0` sets the share of plain arc length and keeps `m` smooth |
| `QAWS_MEASURE_DENSITY` | `density(C(t))` | a positive field in space with gradient and Hessian (importance maps, distance fields) |

A target asks for the point at `sigma = distance + fraction * M_total`: `fraction` 0 is an absolute measure, `distance` 0 a normalized one (`i / (n - 1)` gives `n` equidistributed samples), `distance = -d, fraction = 1` a measure from the end. The parameter `t` solves `M(t) = sigma`; `M` uses a composite Gauss-Legendre rule (every span and the last partial span split into 8 pieces of `quadrature` points, 0 = 8), and every derivative differentiates that same relation:

```
t'  = (sigma' - M'(t)) / m(t)
t'' = (sigma'' - M''(t) - 2 m_e(t) t' - m_t(t) t'^2) / m(t)
p'  = dC + C' t'
p'' = d2C + 2 dC' t' + C'' t'^2 + C' t''
```

`m_e` is the rate of the integrand along the parameter direction at fixed `t`, `m_t` its `t`-derivative (the curvature measure reads `C'''` for it). Knots are parameters of the tangent and the adjoint: a knot moves the basis at fixed `t`, and the span boundaries it bounds move the quadrature nodes `t = a (1 - s) + b s` and weights `W = (b - a) c`, so the cumulative measures differentiate as `(W m)' = W m' + W' m`, `(W m)'' = W m'' + 2 W' m'`. The partial span `[a_k, t]` moves with its start knot only (`t` is held for `M'`), which also covers the moving start of the domain. The HVP takes knot views by polarizing this second order forward pass. Samples where `m` vanishes are reported as `QAWS_DIFF_ILL_CONDITIONED` with zero derivatives. 2D curves are lifted (z = 0).

### qaws_sample_measure_desc

```c
typedef enum qaws_sample_measure
{
	QAWS_MEASURE_ARC_LENGTH = 0,
	QAWS_MEASURE_CURVATURE,
	QAWS_MEASURE_DENSITY
} qaws_sample_measure;

typedef qaws_scalar (*qaws_density_fn)(qaws_vec3 position, void* user_data, qaws_vec3* out_gradient,
	qaws_scalar* out_hessian);   /* hessian: xx, xy, xz, yy, yz, zz */

typedef struct qaws_sample_measure_desc
{
	qaws_sample_measure kind;
	qaws_scalar curvature_floor;     /* QAWS_MEASURE_CURVATURE, > 0 */
	qaws_density_fn density;         /* QAWS_MEASURE_DENSITY */
	void* density_user_data;
} qaws_sample_measure_desc;
```

A NULL measure is arc length. An invalid measure (non-positive floor, missing density) returns `QAWS_STATUS_INVALID_ARGUMENT`.

### qaws_cdf_target / qaws_cdf_sample

```c
typedef struct qaws_cdf_target
{
	qaws_scalar distance;    /* measure from the start */
	qaws_scalar fraction;    /* plus this fraction of the total measure */
} qaws_cdf_target;

typedef struct qaws_cdf_sample
{
	qaws_scalar t;           /* curve parameter */
	qaws_vec3 position;
} qaws_cdf_sample;
```

### qaws_curve_cdf_sample_tangent

```c
qaws_status qaws_curve_cdf_sample_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	qaws_scalar const* distance_tangent,
	qaws_scalar const* distance_tangent2,
	unsigned int count,
	unsigned int quadrature,
	qaws_diff_views const* param_tangent,
	qaws_cdf_sample* out_value,
	qaws_cdf_sample* out_tangent,
	qaws_cdf_sample* out_tangent2,
	qaws_scalar* out_total);
```

Solves every target and writes the samples (`out_value`), their first (`out_tangent`) and second (`out_tangent2`) directional derivatives along `param_tangent` and the per-target distance rates `distance_tangent` and `distance_tangent2` (first and second derivative of the target distance along the same path), and the total measure. Any of the rates, `param_tangent` and the outputs may be NULL.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on a NULL curve, NULL targets or an invalid measure; `QAWS_STATUS_UNSUPPORTED_OPERATION` for a knot view or when second order is asked of a family without `QAWS_CAP_TANGENT2`.

### qaws_curve_cdf_sample_adjoint

```c
qaws_status qaws_curve_cdf_sample_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* distance_adjoint);
```

Accumulates (`+=`) the pullback of the sample adjoints `(t_bar, p_bar)`: `J^T p_bar + lambda (fraction grad M_total - grad M(t))` into `param_adjoint`, `lambda = (t_bar + p_bar . C') / m` into `distance_adjoint`. The measure gradients of the whole batch are one quadrature pass with suffix-summed weights, pulling `dm/d(C, C', C'')` back through the curve jets.

### qaws_curve_cdf_sample_hvp

```c
qaws_status qaws_curve_cdf_sample_hvp(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_sample_measure_desc const* measure,
	qaws_cdf_target const* targets,
	unsigned int count,
	unsigned int quadrature,
	qaws_cdf_sample const* adjoint,
	qaws_diff_views const* direction,
	qaws_diff_views* out_hv);
```

Accumulates (`+=`) `H * direction`, `H` the Hessian with respect to the parameters of `sum_i (t_bar_i t_i + p_bar_i . p_i)`. For families linear in their fields (`QAWS_CAP_LINEAR`) it differentiates the adjoint along the direction (forward over reverse); the integrand's Hessian enters by polarization of its dual second derivative. When `direction` or `out_hv` holds knots, or the family is rational, it polarizes the exact second order forward pass instead: `e_j^T H d = (q(d + e_j) - q(d) - q(e_j)) / 2` with `q(x)` the second derivative of the objective along `x`, about two forward passes per parameter of `out_hv`.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_UNSUPPORTED_OPERATION` for views with children on the polarized path.

### qaws_traversal_cdf_targets

```c
qaws_status qaws_traversal_cdf_targets(
	qaws_traversal const* traversal,
	qaws_scalar const* inputs,
	unsigned int count,
	qaws_cdf_target* out_targets,
	qaws_scalar* out_rate,
	qaws_scalar* out_rate2);
```

Turns traversal inputs (time through easing and motion profile, or arc length) into arc-length targets after the wrap mode, written as `distance + fraction * L_total` so that the derivatives of the total length are carried: loop gives fraction `-k`, ping-pong `(s, -2k)` going out and `(-s, 2k + 2)` coming back, clamp `(0, 1)` past the end. `out_rate` / `out_rate2` (may be NULL) receive the first and second derivatives of the target distance with respect to the input: the profile speed and acceleration through the easing. Passed (scaled by the input rate) as `distance_tangent` / `distance_tangent2` of `qaws_curve_cdf_sample_tangent` with a NULL measure, they give exact first and second time derivatives of the traversal samples; the adjoint and the HVP of `qaws_curve_cdf_sample_*` apply unchanged. The samples match `qaws_traversal_evaluate_*` to the accuracy of the traversal's table.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` on NULL arguments; `QAWS_STATUS_UNSUPPORTED_OPERATION` for parameter-mode and multi-curve traversals.

### Surfaces: qaws_surface_cdf_sample_*

```c
#define QAWS_MEASURE_AREA QAWS_MEASURE_ARC_LENGTH

typedef struct qaws_surface_cdf_sample
{
	qaws_scalar u;
	qaws_scalar v;
	qaws_vec3 position;
} qaws_surface_cdf_sample;

qaws_status qaws_surface_cdf_sample_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi, qaws_scalar const* xi_tangent, unsigned int count,
	unsigned int cells, unsigned int quadrature, qaws_diff_views const* param_tangent,
	qaws_surface_cdf_sample* out_value, qaws_surface_cdf_sample* out_tangent,
	qaws_surface_cdf_sample* out_tangent2, qaws_scalar* out_total);

qaws_status qaws_surface_cdf_sample_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi, unsigned int count, unsigned int cells, unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint, qaws_diff_views* param_adjoint, qaws_scalar* xi_adjoint);

qaws_status qaws_surface_cdf_sample_hvp(
	qaws_diff_context const* ctx, qaws_surface const* surface, qaws_sample_measure_desc const* measure,
	qaws_scalar const* xi, unsigned int count, unsigned int cells, unsigned int quadrature,
	qaws_surface_cdf_sample const* adjoint, qaws_diff_views const* direction, qaws_diff_views* out_hv);
```

Warps points `xi` of the unit square (stratified, blue noise, low discrepancy...) onto a surface by the inverse CDFs of the measure `w = rho(S) |S_u x S_v|` (`QAWS_MEASURE_AREA`, `QAWS_MEASURE_DENSITY`, or `QAWS_MEASURE_CURVATURE` with `rho = sqrt(floor^2 + k1^2 + k2^2)` from the principal curvatures, `k1^2 + k2^2 = 4 H^2 - 2 K`; its rates read the third and, at second order, the fourth derivatives of the patch): `u` solves the marginal `A(u) = xi_u A_total`, `v` the conditional `B(v; u) = xi_v B(v1; u)`, so samples follow the measure and keep the stratification of `xi`. The integrals use a grid of `cells x cells` cells (0 = 6) with `quadrature` Gauss points per cell and direction (0 = 8).

- **tangent**: the forward pass solves the discrete equations in dual numbers (the unknown and the quadrature nodes moving with it carry first and second order rates; one dual Newton step per order), so the derivatives are exact for the discretization, cross terms included. `xi_tangent` (2 per point) differentiates with respect to the points.
- **adjoint**: pulls back through the conditional equation (`mu = v_bar / G_v`) then the marginal one (`lambda = u_bar / F_u`); `xi_adjoint` receives `lambda A_total` and `mu B(v1; u)`. The full cells of the marginal integral are pulled back once for the batch.
- **hvp**: for families linear in their fields (`QAWS_CAP_LINEAR`), forward over reverse: the backward pass runs once in dual numbers along the direction, every coefficient (`mu`, `lambda`, cell weights, moving nodes) carrying its rate, and each node pulls back `J^T (C' dw/dy + C H_w ydot)` plus the moving-node term (the same pullback shifted to the next jet channels by `u'`, `v'`); the slopes `F_u`, `G_u`, `G_v` get their rates by polarizing second order dual passes. Other families (rational patches) polarize the second order forward pass over the parameters of `out_hv` (`e_j^T H d = (q(d + e_j) - q(d) - q(e_j)) / 2`, about two passes per parameter). Views with children are refused.

**Returns:** `QAWS_STATUS_OK`; `QAWS_STATUS_INVALID_ARGUMENT` for a NULL surface or points, a missing density or a non-positive curvature floor; `QAWS_STATUS_UNSUPPORTED_OPERATION` without `QAWS_CAP_TANGENT2` or for knot views. Points where the measure vanishes are reported as `QAWS_DIFF_ILL_CONDITIONED` with zero derivatives.
