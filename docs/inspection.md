# Inspection and Analysis

Qaws provides functions to query curve properties and perform geometric analysis.

## Generic inspection

These work on any curve regardless of family:

```c
qaws_curve_kind kind = qaws_curve_get_kind(curve);
qaws_dimension  dim  = qaws_curve_get_dimension(curve);
unsigned int  degree  = qaws_curve_get_degree(curve);
unsigned int  spans   = qaws_curve_get_span_count(curve);
qaws_range    range   = qaws_curve_get_parameter_range(curve);
int           closed  = qaws_curve_is_closed(curve);
int         periodic  = qaws_curve_is_periodic(curve);
int         rational  = qaws_curve_is_rational(curve);
qaws_continuity cont  = qaws_curve_get_continuity(curve);
```

| Property | Description |
|---|---|
| `kind` | Curve family (BEZIER, HERMITE, CATMULL_ROM, BSPLINE, NURBS, TRAJECTORY, YUKSEL) |
| `dimension` | 2D or 3D |
| `degree` | Polynomial degree of the curve |
| `span_count` | Number of piecewise segments |
| `parameter_range` | Valid parameter domain [min, max] |
| `is_closed` | 1 if the curve forms a closed loop |
| `is_periodic` | 1 if the curve wraps periodically |
| `is_rational` | 1 for NURBS (weighted control points) |
| `continuity` | C0, C1, C2, or C3 at span boundaries |

## Arc length

Compute the distance along the curve between two parameters:

```c
qaws_scalar length;
qaws_curve_compute_arc_length(&length,
    curve, range.min_value, range.max_value);
```

Internally uses Gauss-Legendre quadrature for accuracy.

## Bounding box

Compute an axis-aligned bounding box by sampling:

```c
qaws_vec2 bbmin, bbmax;
qaws_curve_compute_bounds_2d(&bbmin, &bbmax, curve);

qaws_vec3 bbmin3, bbmax3;
qaws_curve_compute_bounds_3d(&bbmin3, &bbmax3, curve);
```

## Closest point

Find the curve parameter nearest to a given point:

```c
qaws_scalar closest_t;

// 2D
qaws_vec2 target2 = { 1.5f, 0.5f };
qaws_curve_find_closest_parameter_2d(&closest_t, curve, target2);

// 3D
qaws_vec3 target3 = { 1.5f, 0.5f, 0.0f };
qaws_curve_find_closest_parameter_3d(&closest_t, curve, target3);

// Evaluate at closest_t to get the actual closest point on the curve
qaws_eval_result_2d result;
qaws_curve_evaluate_2d(&result, curve,
    closest_t, QAWS_EVAL_FLAG_POSITION);
```

Uses coarse sampling followed by Newton refinement.

## Span continuity

Query continuity at a specific span boundary:

```c
qaws_continuity cont;
qaws_curve_get_span_continuity(&cont, curve, 0);
// boundary_index is 0 to span_count-2
```

## Derived geometric helpers

Convenience functions that compute geometric quantities from derivatives:

### Tangent

```c
qaws_vec2 tangent;
qaws_curve_compute_tangent_2d(&tangent, curve, t);
// tangent is a unit vector in the direction of travel
```

### Normal (2D)

```c
qaws_vec2 normal;
qaws_curve_compute_normal_2d(&normal, curve, t);
// Unit normal perpendicular to the tangent (rotated 90 degrees CCW)
```

### Frenet frame (3D)

Compute the full Frenet frame at a point on a 3D curve:

```c
qaws_vec3 tangent, normal, binormal;
qaws_curve_compute_frenet_frame_3d(&tangent, &normal,
    &binormal, curve, t);
// tangent  -- unit tangent vector (direction of travel)
// normal   -- unit principal normal (points toward center of curvature)
// binormal -- tangent x normal, completes the right-handed frame
```

### Curvature

```c
qaws_scalar kappa;
qaws_curve_compute_curvature_2d(&kappa, curve, t);
// 2D: signed curvature (positive = CCW)
// 3D: unsigned curvature |d1 x d2| / |d1|^3
```

### Torsion (3D only)

```c
qaws_scalar tau;
qaws_curve_compute_torsion_3d(&tau, curve, t);
// Rate of change of the osculating plane
```

### Speed

```c
qaws_scalar speed;
qaws_curve_compute_speed(&speed, curve, t);
// |d1| -- magnitude of the first derivative
```

## Inflection points

Find parameter values where the curvature changes sign (inflection points):

```c
qaws_scalar inflections[32];
unsigned int inflection_count;
qaws_curve_find_inflection_points(inflections,
    &inflection_count, curve, 32);
// inflections[0..inflection_count-1] are the parameter values
```

## Extrema detection

Find parameter values where the curve reaches a minimum or maximum along a given axis:

```c
qaws_scalar extrema[32];
unsigned int extrema_count;
unsigned int axis = 0; // 0 = x, 1 = y, 2 = z
qaws_curve_find_extrema(extrema, &extrema_count,
    curve, 32, axis);
// extrema[0..extrema_count-1] are the parameter values
```

## Curvature comb

Sample curvature comb data along the curve for visualization. Each sample contains the position on the curve, the signed curvature value, and the normal direction:

```c
// 2D
qaws_curvature_sample_2d samples_2d[100];
qaws_curve_compute_curvature_comb_2d(// number of samples
    samples_2d,
    curve,        100, // output buffer
    100);       // buffer capacity

// Each sample has:
//   .position  -- point on the curve
//   .curvature -- signed curvature at that point
//   .normal    -- unit normal at that point

// 3D
qaws_curvature_sample_3d samples_3d[100];
qaws_curve_compute_curvature_comb_3d(// number of samples
    samples_3d,
    curve,        100, // output buffer
    100);       // buffer capacity
```

The sample structs are defined as:

```c
typedef struct qaws_curvature_sample_2d {
    qaws_vec2   position;
    qaws_scalar curvature;
    qaws_vec2   normal;
} qaws_curvature_sample_2d;

typedef struct qaws_curvature_sample_3d {
    qaws_vec3   position;
    qaws_scalar curvature;
    qaws_vec3   normal;
} qaws_curvature_sample_3d;
```

## Winding number

Compute the winding number of a closed 2D curve around a point. The winding number counts how many times the curve winds around the point (positive for CCW, negative for CW):

```c
int winding;
qaws_vec2 point = { 0.5f, 0.5f };
qaws_curve_compute_winding_number_2d(&winding, curve, point);
// winding == 0 means the point is outside the curve
// winding != 0 means the point is inside
```

## Curve-curve intersection

Find all intersection points between two curves:

```c
// 2D
qaws_intersection_2d hits_2d[32];
unsigned int hit_count;
qaws_curve_find_intersections_2d(hits_2d, &hit_count,
    curve_a, curve_b, 32);

// Each intersection has:
//   .parameter_a -- parameter on curve_a
//   .parameter_b -- parameter on curve_b
//   .position    -- intersection point in space

// 3D
qaws_intersection_3d hits_3d[32];
qaws_curve_find_intersections_3d(hits_3d, &hit_count,
    curve_a, curve_b, 32);
```

The intersection result structs are defined as:

```c
typedef struct qaws_intersection_2d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec2   position;
} qaws_intersection_2d;

typedef struct qaws_intersection_3d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec3   position;
} qaws_intersection_3d;
```

## Self-intersection detection

Find points where a curve crosses itself:

```c
// 2D
qaws_intersection_2d self_hits_2d[32];
unsigned int self_count;
qaws_curve_find_self_intersections_2d(self_hits_2d,
    &self_count, curve, 32);

// parameter_a and parameter_b are the two distinct parameter values
// where the curve passes through the same position

// 3D
qaws_intersection_3d self_hits_3d[32];
qaws_curve_find_self_intersections_3d(self_hits_3d,
    &self_count, curve, 32);
```

## Family-specific inspection

### Bezier control points

```c
qaws_vec2 cps[10];
unsigned int count;
qaws_status s = qaws_bezier_get_control_points(cps,
    &count, curve, 10);
// Returns QAWS_STATUS_UNSUPPORTED_OPERATION if not a Bezier curve
```

### B-Spline / NURBS knot vector

```c
qaws_scalar knots[64];
unsigned int knot_count;
qaws_bspline_get_knots(knots, &knot_count, curve, 64);
// Works for both B-Spline and NURBS curves
```

### NURBS weights

```c
qaws_scalar weights[32];
unsigned int weight_count;
qaws_nurbs_get_weights(weights, &weight_count, curve, 32);
// Returns QAWS_STATUS_UNSUPPORTED_OPERATION if not NURBS
```
