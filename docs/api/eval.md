# API Reference: qaws_eval.h

Point evaluation of curves (single, per span, batched) and batched surface evaluation. All functions are
thread-safe on immutable curves and surfaces.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)).

Related: result structs and `QAWS_EVAL_FLAG_*` in [types.md](types.md); `qaws_surface_eval_result` and
`QAWS_SURFACE_EVAL_*` in [surface.md](surface.md); concepts in [../evaluation.md](../evaluation.md).

`eval_flags` for curves is a bitmask of `QAWS_EVAL_FLAG_POSITION`, `QAWS_EVAL_FLAG_D1`, `QAWS_EVAL_FLAG_D2`,
`QAWS_EVAL_FLAG_D3` (`QAWS_EVAL_FLAG_NONE` = 0). The result's `valid_flags` says what was actually filled.

---

## Curve evaluation at a global parameter

```c
qaws_status qaws_curve_evaluate_2d(
    qaws_eval_result_2d* out_result,
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags);

qaws_status qaws_curve_evaluate_3d(
    qaws_eval_result_3d* out_result,
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags);
```

Finds the span holding `parameter` and evaluates it there.

- `curve` must have the matching dimension (2D or 3D).
- Returns `QAWS_STATUS_OK`, `QAWS_STATUS_INVALID_ARGUMENT` (null `curve` or `out_result`),
  `QAWS_STATUS_INVALID_DIMENSION` (dimension mismatch), or a status from the curve family's evaluator.

---

## Curve evaluation on one span

```c
qaws_status qaws_curve_evaluate_span_2d(
    qaws_eval_result_2d* out_result,
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags);

qaws_status qaws_curve_evaluate_span_3d(
    qaws_eval_result_3d* out_result,
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags);
```

Evaluates span `span_index` (0 to span_count - 1) at a local parameter in [0, 1]; values outside are clamped.

- Returns `QAWS_STATUS_OK`, `QAWS_STATUS_INVALID_ARGUMENT`, `QAWS_STATUS_INVALID_DIMENSION`, or
  `QAWS_STATUS_OUT_OF_RANGE` (`span_index >= span_count`).

---

## Batched curve evaluation

```c
qaws_status qaws_curve_evaluate_batch_2d(
    qaws_eval_result_2d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);

qaws_status qaws_curve_evaluate_batch_3d(
    qaws_eval_result_3d* out_results,
    qaws_curve const* curve,
    qaws_scalar const* parameters,
    unsigned int count,
    unsigned int eval_flags);
```

Evaluates `curve` at `count` global parameters.

- `parameters`: `count` values. `out_results`: caller-allocated array of `count` structs; each entry is zeroed
  before it is filled.
- Returns `QAWS_STATUS_OK` if all evaluations succeed. On partial failure every entry is still attempted and the
  first error status is returned; `out_results[i].valid_flags == 0` marks a failed entry.
- `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer, `QAWS_STATUS_INVALID_DIMENSION` for a dimension mismatch.
- B-splines of degree 1 to 5 compute the basis functions of several parameters at once, one per SIMD lane (2 or 4
  in f64, 4 or 8 in f32, SSE2 or AVX2; NEON f32), and reuse a span while the parameters stay inside it: each result
  equals the one-point evaluation bit for bit. Sorted parameters are the fastest.

---

## Batched surface evaluation

```c
qaws_status qaws_surface_evaluate_batch(
    qaws_surface_eval_result* out_results,
    qaws_surface const* surface,
    qaws_scalar const* u_params,
    qaws_scalar const* v_params,
    unsigned int count,
    unsigned int eval_flags);
```

Evaluates `surface` at `count` (u, v) pairs, calling `qaws_surface_evaluate` for each.

- `u_params`, `v_params`: `count` values each. `out_results`: caller-allocated array of `count` structs, each
  zeroed before it is filled.
- `eval_flags`: bitmask of `QAWS_SURFACE_EVAL_*` (`POSITION`, `DU`, `DV`, `NORMAL`, `DUU`, `DUV`, `DVV`).
- Returns `QAWS_STATUS_OK` if all succeed, otherwise the first error status (all entries attempted);
  `QAWS_STATUS_INVALID_ARGUMENT` for a null pointer.
