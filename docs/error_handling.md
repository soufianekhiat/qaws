# Error Handling

All Qaws functions that can fail return a `qaws_status` code. Check it before using the output.

## Status codes

| Code | Meaning |
|---|---|
| `QAWS_STATUS_OK` | Success |
| `QAWS_STATUS_INVALID_ARGUMENT` | NULL pointer or invalid input |
| `QAWS_STATUS_INVALID_DIMENSION` | Dimension mismatch (e.g., evaluating 3D on a 2D curve) |
| `QAWS_STATUS_INVALID_DEGREE` | Unsupported degree for the curve family |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | Wrong number of control points for degree |
| `QAWS_STATUS_INVALID_KNOT_COUNT` | Knot vector length doesn't match control points + degree + 1 |
| `QAWS_STATUS_INVALID_KNOT_VECTOR` | Knot vector is not non-decreasing |
| `QAWS_STATUS_INVALID_WEIGHT_COUNT` | Weight count doesn't match control point count |
| `QAWS_STATUS_INVALID_PARAMETER_RANGE` | Parameter range is degenerate (min >= max) |
| `QAWS_STATUS_OUT_OF_RANGE` | Parameter outside curve domain |
| `QAWS_STATUS_UNSUPPORTED_OPERATION` | Operation not available for this curve family |
| `QAWS_STATUS_DEGENERATE_CURVE` | Curve is degenerate (e.g., zero-length) |
| `QAWS_STATUS_NUMERICAL_FAILURE` | Numerical solver failed to converge |
| `QAWS_STATUS_BUFFER_TOO_SMALL` | Caller's buffer is too small for the result |
| `QAWS_STATUS_ALLOCATION_FAILURE` | malloc failed |
| `QAWS_STATUS_INTERNAL_ERROR` | Bug in the library |

## Converting to string

```c
qaws_status s = qaws_curve_create_bezier(&curve, &desc);
if (s != QAWS_STATUS_OK) {
    printf("Error: %s\n", qaws_status_to_string(s));
}
```

## Pattern

```c
qaws_curve* curve = NULL;
qaws_status s = qaws_curve_create_bezier(&curve, &desc);
if (s != QAWS_STATUS_OK) {
    // handle error -- curve is NULL
    return;
}
// use curve...
qaws_curve_destroy(curve);
```

On error, output pointers (`out_curve`, `out_result`, etc.) are not modified or are set to NULL/zero.

## Family-specific errors

Calling a family-specific inspection function on the wrong curve type returns `QAWS_STATUS_UNSUPPORTED_OPERATION`:

```c
// curve is a Hermite curve
qaws_scalar knots[10];
unsigned int count;
qaws_status s = qaws_bspline_get_knots(knots, &count, curve, 10);
// s == QAWS_STATUS_UNSUPPORTED_OPERATION
```
