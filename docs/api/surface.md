# API Reference: qaws_surface.h, qaws_surface_types.h

The generic surface object: its kinds, evaluation flags and results, and the
functions shared by every surface family. Surfaces are made by the family
headers ([surface_families.md](surface_families.md),
[surface_trim.md](surface_trim.md), [surface_tspline.md](surface_tspline.md),
[surface_subdiv.md](surface_subdiv.md)) and by views ([view.md](view.md)).

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)).

Thread safety: a `qaws_surface` is immutable after creation. Any number of
threads may evaluate or sample the same surface at once.

---

## Types

```c
typedef struct qaws_surface qaws_surface;   /* opaque */
```

### qaws_surface_kind

| Member | Made by |
|---|---|
| `QAWS_SURFACE_KIND_INVALID` (0) | returned for a NULL surface |
| `QAWS_SURFACE_KIND_BEZIER` | `qaws_surface_create_bezier` |
| `QAWS_SURFACE_KIND_BSPLINE` | `qaws_surface_create_bspline` |
| `QAWS_SURFACE_KIND_NURBS` | `qaws_surface_create_nurbs` |
| `QAWS_SURFACE_KIND_SWEPT` | `qaws_surface_create_swept` |
| `QAWS_SURFACE_KIND_RULED` | `qaws_surface_create_ruled` |
| `QAWS_SURFACE_KIND_COONS` | `qaws_surface_create_coons` |
| `QAWS_SURFACE_KIND_EXTRUSION` | `qaws_surface_create_extrusion` |
| `QAWS_SURFACE_KIND_REVOLUTION` | `qaws_surface_create_revolution` |
| `QAWS_SURFACE_KIND_PIPE` | `qaws_surface_create_pipe` |
| `QAWS_SURFACE_KIND_BILINEAR` | `qaws_surface_create_bilinear` |
| `QAWS_SURFACE_KIND_BIQUADRATIC` | `qaws_surface_create_biquadratic` |
| `QAWS_SURFACE_KIND_LOFT` | `qaws_surface_create_loft` |
| `QAWS_SURFACE_KIND_GORDON` | `qaws_surface_create_gordon` |
| `QAWS_SURFACE_KIND_OFFSET` | `qaws_surface_create_offset` |
| `QAWS_SURFACE_KIND_TRIMMED` | `qaws_surface_create_trimmed` |
| `QAWS_SURFACE_KIND_TSPLINE` | `qaws_surface_create_tspline` |
| `QAWS_SURFACE_KIND_SUBDIV` | `qaws_surface_create_subdiv` |
| `QAWS_SURFACE_KIND_FILLET` | a fillet surface kind; no public header in this release creates one |
| `QAWS_SURFACE_KIND_VIEW` | a borrowed surface evaluated from host arrays it does not own ([view.md](view.md)) |

### qaws_surface_eval_flags

Bit flags passed as `eval_flags` to `qaws_surface_evaluate`.

| Flag | Value | Result field |
|---|---|---|
| `QAWS_SURFACE_EVAL_NONE` | 0 | - |
| `QAWS_SURFACE_EVAL_POSITION` | 1 << 0 | `position` |
| `QAWS_SURFACE_EVAL_DU` | 1 << 1 | `du`, dS/du |
| `QAWS_SURFACE_EVAL_DV` | 1 << 2 | `dv`, dS/dv |
| `QAWS_SURFACE_EVAL_NORMAL` | 1 << 3 | `normal`, unit (du x dv) |
| `QAWS_SURFACE_EVAL_DUU` | 1 << 4 | `duu`, d^2S/du^2 |
| `QAWS_SURFACE_EVAL_DUV` | 1 << 5 | `duv`, d^2S/dudv |
| `QAWS_SURFACE_EVAL_DVV` | 1 << 6 | `dvv`, d^2S/dv^2 |

### qaws_surface_eval_result

```c
typedef struct qaws_surface_eval_result {
	qaws_vec3 position;
	qaws_vec3 du;       /* dS/du */
	qaws_vec3 dv;       /* dS/dv */
	qaws_vec3 normal;   /* unit normal */
	qaws_vec3 duu;      /* d^2S/du^2 */
	qaws_vec3 duv;      /* d^2S/dudv */
	qaws_vec3 dvv;      /* d^2S/dv^2 */
	unsigned int valid_flags;   /* QAWS_SURFACE_EVAL_* bits of the fields filled */
} qaws_surface_eval_result;
```

---

## Properties

```c
qaws_surface_kind qaws_surface_get_kind(qaws_surface const* surface);
unsigned int      qaws_surface_get_u_degree(qaws_surface const* surface);
unsigned int      qaws_surface_get_v_degree(qaws_surface const* surface);
qaws_range        qaws_surface_get_u_range(qaws_surface const* surface);
qaws_range        qaws_surface_get_v_range(qaws_surface const* surface);
int               qaws_surface_is_rational(qaws_surface const* surface);
```

On a NULL surface these return `QAWS_SURFACE_KIND_INVALID`, 0, the range
{0, 0} and 0. `qaws_surface_is_rational` returns 1 for a rational surface
(NURBS, weighted T-spline and so on), else 0.

## Evaluation

```c
qaws_status qaws_surface_evaluate(qaws_surface_eval_result* out_result, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int eval_flags);
```

Evaluates the quantities named by `eval_flags` at (u, v). The parameters are
clamped to the surface's u and v ranges. The result is zeroed first; asking
for `QAWS_SURFACE_EVAL_NORMAL` also computes `du` and `dv`. Returns
`QAWS_STATUS_INVALID_ARGUMENT` on a NULL argument and
`QAWS_STATUS_UNSUPPORTED_OPERATION` for a surface that cannot be evaluated.

```c
qaws_status qaws_surface_get_normal(qaws_vec3* out_normal, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v);
```

The unit normal at (u, v): `qaws_surface_evaluate` with
`QAWS_SURFACE_EVAL_NORMAL`.

## Sampling

```c
qaws_status qaws_surface_sample(qaws_vec3* out_positions, unsigned int* out_position_count,
	qaws_surface const* surface, unsigned int u_count, unsigned int v_count,
	unsigned int position_capacity);
```

Positions on a uniform `u_count x v_count` grid spanning both ranges,
row-major with v varying fastest (`out_positions[i * v_count + j]`).
`out_positions` must hold `u_count * v_count` elements.

- `QAWS_STATUS_INVALID_ARGUMENT`: a NULL pointer, `u_count` or `v_count`
  below 2, or a grid whose count overflows `unsigned int`.
- `QAWS_STATUS_BUFFER_TOO_SMALL`: `u_count * v_count > position_capacity`.

## Lifetime

```c
void qaws_surface_destroy(qaws_surface* surface);
```

Frees a surface made by any `qaws_surface_create_*` function. NULL is
accepted. Surfaces that borrow curves or other surfaces do not free them;
destroy the borrowing surface first.
