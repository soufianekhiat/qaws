# qaws_status.h: status codes

Every fallible function returns a `qaws_status`. The header is precision-free (one declaration for both precisions).
See [../error_handling.md](../error_handling.md) for how outputs behave on failure.

## qaws_status

| Value | Meaning |
|---|---|
| `QAWS_STATUS_OK` (0) | success |
| `QAWS_STATUS_INVALID_ARGUMENT` | a NULL, out-of-range or otherwise unusable argument |
| `QAWS_STATUS_INVALID_DIMENSION` | dimension not 2D or 3D, or dimensions that do not match |
| `QAWS_STATUS_INVALID_DEGREE` | degree outside what the family accepts |
| `QAWS_STATUS_INVALID_CONTROL_POINT_COUNT` | wrong number of control points (or of a paired array) |
| `QAWS_STATUS_INVALID_KNOT_COUNT` | knot count not control_point_count + degree + 1, or no knots |
| `QAWS_STATUS_INVALID_KNOT_VECTOR` | knots decreasing, or no non-empty span |
| `QAWS_STATUS_INVALID_WEIGHT_COUNT` | wrong weight count, or a weight not positive |
| `QAWS_STATUS_INVALID_PARAMETER_RANGE` | an empty or reversed parameter range, or times not increasing |
| `QAWS_STATUS_OUT_OF_RANGE` | a parameter or index outside its domain |
| `QAWS_STATUS_UNSUPPORTED_OPERATION` | the operation is not available for this input or setting |
| `QAWS_STATUS_DEGENERATE_CURVE` | the input describes a degenerate curve (zero length, too few distinct points) |
| `QAWS_STATUS_NUMERICAL_FAILURE` | a solve or iteration gave no usable result |
| `QAWS_STATUS_BUFFER_TOO_SMALL` | a caller buffer (or fixed storage) is too small |
| `QAWS_STATUS_ALLOCATION_FAILURE` | an allocation failed |
| `QAWS_STATUS_INTERNAL_ERROR` | an internal invariant failed |
| `QAWS_STATUS_EXACT_UNSUPPORTED` | exact kernel: the operation is outside the rational domain |
| `QAWS_STATUS_EXACT_RANGE_EXCEEDED` | exact kernel: the integer budget was exceeded |
| `QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE` | exact kernel: exact spaces that do not match |
| `QAWS_STATUS_CERTIFICATION_FAILED` | exact kernel: a sign could not be certified |
| `QAWS_STATUS_NOT_CONVERGED` | a result defined at convergence asked for before it (repulsive solution derivatives) |
| `QAWS_STATUS_INVALID_TOPOLOGY` | a mesh edit refused by its validity checks; the mesh is unchanged |
| `QAWS_STATUS_CANCELLED` | stopped by the host's cancel callback; outputs are as the operation documents |

## qaws_status_to_string

```c
char const* qaws_status_to_string(qaws_status status);
```

Returns the enumerator's name as a static string (for example `"QAWS_STATUS_INVALID_DEGREE"`), or
`"QAWS_STATUS_UNKNOWN"` for a value outside the enum. Never NULL; do not free it.
