# API Reference: qaws_status.h

Error codes returned by all fallible Qaws functions.

---

## qaws_status

```c
typedef enum qaws_status {
    QAWS_STATUS_OK = 0,
    QAWS_STATUS_INVALID_ARGUMENT,
    QAWS_STATUS_INVALID_DIMENSION,
    QAWS_STATUS_INVALID_DEGREE,
    QAWS_STATUS_INVALID_CONTROL_POINT_COUNT,
    QAWS_STATUS_INVALID_KNOT_COUNT,
    QAWS_STATUS_INVALID_KNOT_VECTOR,
    QAWS_STATUS_INVALID_WEIGHT_COUNT,
    QAWS_STATUS_INVALID_PARAMETER_RANGE,
    QAWS_STATUS_OUT_OF_RANGE,
    QAWS_STATUS_UNSUPPORTED_OPERATION,
    QAWS_STATUS_DEGENERATE_CURVE,
    QAWS_STATUS_NUMERICAL_FAILURE,
    QAWS_STATUS_BUFFER_TOO_SMALL,
    QAWS_STATUS_ALLOCATION_FAILURE,
    QAWS_STATUS_INTERNAL_ERROR,
    QAWS_STATUS_EXACT_UNSUPPORTED,          /* exact kernel: outside the rational domain */
    QAWS_STATUS_EXACT_RANGE_EXCEEDED,       /* exact kernel: integer budget exceeded */
    QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE,   /* exact kernel: spaces that do not match */
    QAWS_STATUS_CERTIFICATION_FAILED,       /* exact kernel: a sign not certified */
    QAWS_STATUS_NOT_CONVERGED,              /* a result defined at convergence asked for before it */
    QAWS_STATUS_INVALID_TOPOLOGY,           /* a mesh edit refused by its validity checks */
    QAWS_STATUS_CANCELLED                   /* stopped by the host's cancel callback */
} qaws_status;
```

---

## qaws_status_to_string

```c
char const* qaws_status_to_string(qaws_status status);
```

Returns a human-readable string for the given status code. Never returns NULL.

**Example:**
```c
printf("Error: %s\n", qaws_status_to_string(QAWS_STATUS_INVALID_DEGREE));
// prints: "Error: QAWS_STATUS_INVALID_DEGREE"
```
