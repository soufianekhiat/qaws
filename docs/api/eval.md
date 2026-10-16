# API Reference: qaws_eval.h

Curve evaluation functions.

---

## qaws_curve_evaluate_2d

```c
qaws_status qaws_curve_evaluate_2d(
    qaws_eval_result_2d* out_result,
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags);
```

Evaluates a 2D curve at a global parameter.

**Parameters:**
- `curve` -- The curve. Must have dimension 2D.
- `parameter` -- Parameter value within the curve's domain.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*` specifying what to compute.
- `out_result` -- Filled with position and/or derivatives. Zeroed before filling.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_evaluate_3d

```c
qaws_status qaws_curve_evaluate_3d(
    qaws_eval_result_3d* out_result,
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags);
```

Same as `evaluate_2d` but for 3D curves.

---

## qaws_curve_evaluate_span_2d

```c
qaws_status qaws_curve_evaluate_span_2d(
    qaws_eval_result_2d* out_result,
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags);
```

Evaluates a 2D curve at a specific span with a local parameter.

**Parameters:**
- `span_index` -- 0-based span index (0 to span_count-1).
- `local_parameter` -- Parameter within [0, 1] for that span.
- Other parameters same as `evaluate_2d`.

---

## qaws_curve_evaluate_span_3d

```c
qaws_status qaws_curve_evaluate_span_3d(
    qaws_eval_result_3d* out_result,
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags);
```

Same as `evaluate_span_2d` but for 3D curves.

---

## qaws_curve_evaluate_batch_2d

```c
qaws_status qaws_curve_evaluate_batch_2d(
    qaws_eval_result_2d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);
```

Evaluates a 2D curve at multiple parameters in one call.

**Parameters:**
- `curve` -- The curve. Must have dimension 2D.
- `parameters` -- Array of `count` parameter values.
- `count` -- Number of parameters to evaluate.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*` specifying what to compute.
- `out_results` -- Caller-allocated array of `count` result structs.

**Returns:** `QAWS_STATUS_OK` if all evaluations succeed. On partial failure, evaluates as many as possible and returns the first error status; `out_results[i].valid_flags == 0` indicates that entry failed.

---

## qaws_curve_evaluate_batch_3d

```c
qaws_status qaws_curve_evaluate_batch_3d(
    qaws_eval_result_3d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);
```

Same as `evaluate_batch_2d` but for 3D curves.

---

## qaws_surface_evaluate_batch

```c
qaws_status qaws_surface_evaluate_batch(
    qaws_surface_eval_result* out_results,
    qaws_surface const* surface,
    qaws_scalar const* u_params,
    qaws_scalar const* v_params,
    unsigned int count,
    unsigned int eval_flags);
```

Evaluates a surface at multiple (u, v) parameter pairs in one call.

**Parameters:**
- `surface` -- The surface to evaluate.
- `u_params` -- Array of `count` u-parameter values.
- `v_params` -- Array of `count` v-parameter values.
- `count` -- Number of parameter pairs to evaluate.
- `eval_flags` -- Bitmask of `QAWS_SURFACE_EVAL_*` specifying what to compute.
- `out_results` -- Caller-allocated array of `count` `qaws_surface_eval_result` structs.

**Returns:** `QAWS_STATUS_OK` if all evaluations succeed.
