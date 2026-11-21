# qaws_types.h: core types

The vector types, enums, descriptor structs and opaque handles shared by the curve API, plus the fixed-size inline
curve. `qaws_scalar` and `QAWS_SCALAR_IS_FLOAT` come from `qaws_platform.h`; every type here exists once per precision
(see [precision.md](precision.md)). The allocator types this header pulls in (`qaws_allocator.h`) are documented in
[alloc.md](alloc.md).

## Thread safety

- `qaws_curve` objects are immutable after creation: any number of threads may evaluate, sample, inspect or traverse
  the same curve without synchronization.
- `qaws_traversal` objects are safe for concurrent reads (evaluate, map queries); `qaws_traversal_advance_*` and
  `qaws_traversal_reset` mutate state and need external synchronization.
- Creation and destruction are not thread-safe on the same object: do not destroy a curve another thread is using.

## Vectors

```c
typedef struct qaws_vec2 { qaws_scalar x; qaws_scalar y; } qaws_vec2;
typedef struct qaws_vec3 { qaws_scalar x; qaws_scalar y; qaws_scalar z; } qaws_vec3;
```

Descriptors that take `void const*` point arrays expect `dimension` scalars per point, packed: an array of `qaws_vec2`
or `qaws_vec3` has that layout.

## Enums

### qaws_dimension

```c
typedef enum qaws_dimension { QAWS_DIMENSION_2D = 2, QAWS_DIMENSION_3D = 3 } qaws_dimension;
```

### qaws_curve_kind

What `qaws_curve_get_kind` reports ([inspect.md](inspect.md)).

| Value | Created by |
|---|---|
| `QAWS_CURVE_KIND_INVALID` (0) | - |
| `QAWS_CURVE_KIND_BEZIER` | [bezier.md](bezier.md) |
| `QAWS_CURVE_KIND_HERMITE` | [hermite.md](hermite.md) |
| `QAWS_CURVE_KIND_CATMULL_ROM` | [catmull_rom.md](catmull_rom.md) |
| `QAWS_CURVE_KIND_BSPLINE` | [bspline.md](bspline.md) |
| `QAWS_CURVE_KIND_NURBS` | [nurbs.md](nurbs.md) |
| `QAWS_CURVE_KIND_TRAJECTORY` | [trajectory.md](trajectory.md) |
| `QAWS_CURVE_KIND_YUKSEL` | [yuksel.md](yuksel.md) |
| `QAWS_CURVE_KIND_RATIONAL_BEZIER` | [rational_bezier.md](rational_bezier.md) |
| `QAWS_CURVE_KIND_COMPOSITE` | [composite.md](composite.md), [interpolating.md](interpolating.md) |
| `QAWS_CURVE_KIND_ARC` | [arc.md](arc.md) |
| `QAWS_CURVE_KIND_POLYNOMIAL` | [polynomial.md](polynomial.md) |
| `QAWS_CURVE_KIND_CLOTHOID` | [clothoid.md](clothoid.md) |
| `QAWS_CURVE_KIND_SUBDIVISION` | [subdivision.md](subdivision.md) |
| `QAWS_CURVE_KIND_REPARAMETERIZED` | a curve with a remapped parameter |
| `QAWS_CURVE_KIND_VIEW` | a borrowed curve evaluated from host arrays it does not own ([view.md](view.md)) |

### qaws_parameterization

`QAWS_PARAMETERIZATION_DEFAULT` (0), `QAWS_PARAMETERIZATION_UNIFORM`, `QAWS_PARAMETERIZATION_CHORDAL`,
`QAWS_PARAMETERIZATION_CENTRIPETAL`. Knot spacing for point-interpolating curves (Catmull-Rom: `DEFAULT` is
centripetal).

### qaws_continuity

`QAWS_CONTINUITY_NONE` (0), `QAWS_CONTINUITY_C0`, `QAWS_CONTINUITY_C1`, `QAWS_CONTINUITY_C2`, `QAWS_CONTINUITY_C3`.

### qaws_traversal_mode

`QAWS_TRAVERSAL_MODE_PARAMETER` (0), `QAWS_TRAVERSAL_MODE_ARC_LENGTH`, `QAWS_TRAVERSAL_MODE_TIME`. How an input value
is read: as a curve parameter, a distance along the curve, or a time.

### qaws_motion_profile

| Value | Meaning |
|---|---|
| `QAWS_MOTION_PROFILE_NONE` (0) | no speed control |
| `QAWS_MOTION_PROFILE_CONSTANT_SPEED` | uniform speed |
| `QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION` | constant acceleration |
| `QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED` | accelerate, cruise, decelerate |
| `QAWS_MOTION_PROFILE_S_CURVE` | jerk-limited ramps |
| `QAWS_MOTION_PROFILE_CUSTOM` | speed from a `qaws_speed_fn` |

### qaws_easing

`QAWS_EASING_LINEAR` (0), `QAWS_EASING_QUAD_IN`, `QAWS_EASING_QUAD_OUT`, `QAWS_EASING_QUAD_IN_OUT`,
`QAWS_EASING_CUBIC_IN`, `QAWS_EASING_CUBIC_OUT`, `QAWS_EASING_CUBIC_IN_OUT`, `QAWS_EASING_SINE_IN`,
`QAWS_EASING_SINE_OUT`, `QAWS_EASING_SINE_IN_OUT`. Applied to traversal input values.

### qaws_wrap_mode

| Value | Meaning |
|---|---|
| `QAWS_WRAP_CLAMP` (0) | clamp to the domain |
| `QAWS_WRAP_LOOP` | wrap around to the start |
| `QAWS_WRAP_PING_PONG` | reverse direction at the ends |

### qaws_eval_flags

Bitmask, combine with `|`: `QAWS_EVAL_FLAG_NONE` = 0, `QAWS_EVAL_FLAG_POSITION` = 1 << 0, `QAWS_EVAL_FLAG_D1` = 1 << 1,
`QAWS_EVAL_FLAG_D2` = 1 << 2, `QAWS_EVAL_FLAG_D3` = 1 << 3.

## Callback

```c
typedef qaws_scalar (*qaws_speed_fn)(qaws_scalar distance, void* user_data);
```

Returns the speed at a distance along the curve; used with `QAWS_MOTION_PROFILE_CUSTOM`.

## Structs

### qaws_range

```c
typedef struct qaws_range { qaws_scalar min_value; qaws_scalar max_value; } qaws_range;
```

A parameter interval.

### qaws_eval_result_2d, qaws_eval_result_3d

```c
typedef struct qaws_eval_result_2d {
    qaws_vec2 position;
    qaws_vec2 d1;
    qaws_vec2 d2;
    qaws_vec2 d3;
    unsigned int valid_flags;   /* QAWS_EVAL_FLAG_* bits of the fields written */
} qaws_eval_result_2d;
```

`qaws_eval_result_3d` has the same members, each a `qaws_vec3`. Filled by the evaluation calls ([eval.md](eval.md)).

### qaws_sampling_desc

```c
typedef struct qaws_sampling_desc {
    qaws_traversal_mode traversal_mode;   /* parameter or arc-length spacing */
    unsigned int sample_count;
    qaws_scalar error_tolerance;          /* for adaptive sampling */
    int use_adaptive_sampling;            /* 0: uniform, non-zero: adaptive */
} qaws_sampling_desc;
```

See [sampling.md](sampling.md).

### qaws_traversal_desc

```c
typedef struct qaws_traversal_desc {
    qaws_traversal_mode traversal_mode;
    qaws_motion_profile motion_profile;
    qaws_scalar speed;
    qaws_scalar acceleration;
    qaws_scalar max_speed;
    qaws_scalar start_time;
    qaws_scalar end_time;
    int clamp_to_domain;
    qaws_easing easing;
    qaws_wrap_mode wrap_mode;
    qaws_scalar jerk;                    /* max jerk for the S-curve profile */
    qaws_speed_fn custom_speed_fn;       /* for QAWS_MOTION_PROFILE_CUSTOM */
    void* custom_speed_fn_user_data;     /* passed to custom_speed_fn */
} qaws_traversal_desc;
```

See [traversal.md](traversal.md) for how each member is used.

## Opaque handles

```c
typedef struct qaws_curve qaws_curve;
typedef struct qaws_traversal qaws_traversal;
```

Created and destroyed through the API only ([curve.md](curve.md), [traversal.md](traversal.md)).

## Inline curve

```c
#define QAWS_INLINE_STORAGE_BYTES 512

typedef union qaws_curve_inline {
    char storage[QAWS_INLINE_STORAGE_BYTES];   /* opaque: laid out by the init functions */
    double align_double;
    void* align_pointer;
} qaws_curve_inline;

qaws_curve* qaws_curve_inline_get_curve(qaws_curve_inline* ic);
```

A curve held in caller storage with no heap allocation, for single-span curves up to degree 7 in 3D (8 control
points x 3 scalars, plus the curve's own fields, within 512 bytes). Initialize it with the `qaws_curve_init_*_inline`
functions ([inline.md](inline.md)); never call `qaws_curve_destroy` on it.

`qaws_curve_inline_get_curve` returns the usable curve pointer. The curve points into its own storage, so after a copy
(memcpy or assignment) call it on the copy once, before other threads use the copy: it re-points the copy at itself.
A pointer taken from the original keeps reading the original.
