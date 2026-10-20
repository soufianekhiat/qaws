# API Reference: qaws_traversal.h

Traversal creation, evaluation, and parameter mapping.

---

## qaws_traversal_create

```c
qaws_status qaws_traversal_create(
    qaws_traversal** out_traversal,
    qaws_curve const* curve,
    qaws_traversal_desc const* desc);
```

Creates a traversal object that wraps a curve with a motion profile.

**Parameters:**
- `curve` -- The curve to traverse. Must remain alive for the traversal's lifetime.
- `desc` -- Traversal configuration.
- `out_traversal` -- Receives the created traversal.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The traversal does **not** own the curve.
- Internally builds an arc-length lookup table for distance-based operations.

---

## qaws_traversal_destroy

```c
void qaws_traversal_destroy(qaws_traversal* traversal);
```

Destroys a traversal object. May be NULL (no-op).

---

## qaws_traversal_evaluate_2d

```c
qaws_status qaws_traversal_evaluate_2d(
    qaws_eval_result_2d* out_result,
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags);
```

Evaluates the traversal at the given input value.

**Parameters:**
- `input_value` -- Interpreted based on traversal mode: parameter, distance, or time.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*`.
- `out_result` -- Receives position and derivatives.

---

## qaws_traversal_evaluate_3d

```c
qaws_status qaws_traversal_evaluate_3d(
    qaws_eval_result_3d* out_result,
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags);
```

Same as `evaluate_2d` but for 3D curves.

---

## qaws_traversal_map_time_to_parameter

```c
qaws_status qaws_traversal_map_time_to_parameter(
    qaws_scalar* out_parameter,
    qaws_traversal const* traversal,
    qaws_scalar time_value);
```

Converts a time value to a curve parameter using the traversal's motion profile.

---

## qaws_traversal_map_distance_to_parameter

```c
qaws_status qaws_traversal_map_distance_to_parameter(
    qaws_scalar* out_parameter,
    qaws_traversal const* traversal,
    qaws_scalar distance_value);
```

Converts a distance-along-curve value to a curve parameter.

---

## qaws_traversal_map_parameter_to_distance

```c
qaws_status qaws_traversal_map_parameter_to_distance(
    qaws_scalar* out_distance,
    qaws_traversal const* traversal,
    qaws_scalar parameter);
```

Converts a curve parameter to a distance-along-curve value.

---

## qaws_traversal_create_multi

```c
qaws_status qaws_traversal_create_multi(
    qaws_traversal** out_traversal,
    qaws_curve const* const* curves,
    unsigned int curve_count,
    qaws_traversal_desc const* desc);
```

Creates a traversal that spans multiple curves joined end-to-end.

**Parameters:**
- `curves` -- Array of `curve_count` curve pointers. All must remain alive for the traversal's lifetime.
- `curve_count` -- Number of curves in the array.
- `desc` -- Traversal configuration.
- `out_traversal` -- Receives the created traversal.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The traversal does **not** own the curves.
- The domain is the concatenation of the individual curve domains.

---

## qaws_traversal_advance_2d

```c
qaws_status qaws_traversal_advance_2d(
    qaws_eval_result_2d* out_result,
    qaws_traversal* traversal,
    qaws_scalar delta_time,
    unsigned int eval_flags);
```

Advances the traversal by a time step and evaluates the curve at the new position. Mutates internal traversal state.

**Parameters:**
- `traversal` -- The traversal (mutable; not thread-safe with concurrent calls).
- `delta_time` -- Time increment to advance by.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*`.
- `out_result` -- Receives position and derivatives at the new time.

---

## qaws_traversal_advance_3d

```c
qaws_status qaws_traversal_advance_3d(
    qaws_eval_result_3d* out_result,
    qaws_traversal* traversal,
    qaws_scalar delta_time,
    unsigned int eval_flags);
```

Same as `advance_2d` but for 3D curves.

---

## qaws_traversal_reset

```c
qaws_status qaws_traversal_reset(
    qaws_traversal* traversal);
```

Resets the traversal to its initial state (time = start_time). Mutates internal state and requires external synchronization if used from multiple threads.
