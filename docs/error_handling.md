# Error Handling

Every Qaws function that can fail returns a `qaws_status` (`qaws_status.h`,
precision-free). Check it before using any output.

## Status codes

| Code | Meaning |
|---|---|
| `QAWS_STATUS_OK` | Success |
| `QAWS_STATUS_INVALID_ARGUMENT` | NULL pointer, non-finite value or otherwise invalid input |
| `QAWS_STATUS_INVALID_DIMENSION` | Dimension mismatch (e.g. a `_3d` call on a 2D curve) |
| `QAWS_STATUS_INVALID_DEGREE` | Degree not supported by the family |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | Wrong number of control points for the degree |
| `QAWS_STATUS_INVALID_KNOT_COUNT` | Knot vector length does not match |
| `QAWS_STATUS_INVALID_KNOT_VECTOR` | Knots not non-decreasing, or an empty domain |
| `QAWS_STATUS_INVALID_WEIGHT_COUNT` | Weight count does not match the control points |
| `QAWS_STATUS_INVALID_PARAMETER_RANGE` | Degenerate parameter range |
| `QAWS_STATUS_OUT_OF_RANGE` | Parameter or extent outside what the call accepts |
| `QAWS_STATUS_UNSUPPORTED_OPERATION` | Operation not available for this kind of object |
| `QAWS_STATUS_DEGENERATE_CURVE` | Degenerate geometry (no extent, vanishing speed, too few distinct points) |
| `QAWS_STATUS_NUMERICAL_FAILURE` | A solve gave no finite or converged result |
| `QAWS_STATUS_BUFFER_TOO_SMALL` | Caller's buffer too small (many calls report the needed count) |
| `QAWS_STATUS_ALLOCATION_FAILURE` | The allocator returned NULL |
| `QAWS_STATUS_INTERNAL_ERROR` | A bug in the library |
| `QAWS_STATUS_EXACT_UNSUPPORTED` | Exact kernel: operation outside the rational domain |
| `QAWS_STATUS_EXACT_RANGE_EXCEEDED` | Exact kernel: integer budget exceeded |
| `QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE` | Exact kernel: exact spaces that do not match |
| `QAWS_STATUS_CERTIFICATION_FAILED` | Exact kernel: a sign could not be certified |
| `QAWS_STATUS_NOT_CONVERGED` | A result defined at convergence asked for before it |
| `QAWS_STATUS_INVALID_TOPOLOGY` | A mesh edit refused by its validity checks (the mesh unchanged); no public function documents it |
| `QAWS_STATUS_CANCELLED` | Stopped by a cancel callback; no public function takes one |

What a given function returns, and when, is stated in its header and in its
API page ([api/status.md](api/status.md), [api/index.md](api/index.md)).

## Converting to a string

```c
char const* qaws_status_to_string(qaws_status status);
```

Returns the code's name (`"QAWS_STATUS_OK"`, ...), `"QAWS_STATUS_UNKNOWN"`
for a value outside the enum.

## Pattern

```c
qaws_curve* curve = NULL;
qaws_status s = qaws_curve_create_bezier(&curve, &desc);
if (s != QAWS_STATUS_OK) {
    fprintf(stderr, "create failed: %s\n", qaws_status_to_string(s));
    return;   /* nothing to destroy */
}
/* use curve ... */
qaws_curve_destroy(curve);
```

- On failure a creator returns no object: there is nothing to destroy, and
  inputs whose ownership would have moved (a composite's segments) stay with
  the caller.
- Other outputs are written only as the function documents. Calls that fill
  caller buffers and return `QAWS_STATUS_BUFFER_TOO_SMALL` usually report the
  needed count first (tessellation, adaptive sampling, knot and control point
  queries), so the call can be repeated with a larger buffer.
- Batch calls return the first failing status; `qaws_curve_evaluate_batch_*`
  marks failed entries with `valid_flags == 0`.

## Family-specific errors

A family-specific function called on another kind returns
`QAWS_STATUS_UNSUPPORTED_OPERATION`:

```c
/* curve is a Hermite curve */
qaws_scalar knots[10];
unsigned int count;
qaws_status s = qaws_bspline_get_knots(knots, &count, curve, 10);
/* s == QAWS_STATUS_UNSUPPORTED_OPERATION */
```

Borrowed curves and surfaces (`QAWS_CURVE_KIND_VIEW`,
`QAWS_SURFACE_KIND_VIEW`) also refuse operations that read a family's own
arrays (control point and knot inspection, transforms, derivatives with
respect to control points).
