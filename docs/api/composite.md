# API Reference: qaws_composite.h

Composite / piecewise multi-segment curve creation.

---

## qaws_composite_desc

```c
typedef struct qaws_composite_desc {
    qaws_dimension dimension;
    qaws_curve** segments;
    unsigned int segment_count;
} qaws_composite_desc;
```

| Field | Description |
|---|---|
| `dimension` | `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D` |
| `segments` | Array of segment curve pointers. Ownership transfers to the composite. |
| `segment_count` | Number of segments (>= 1) |

---

## qaws_curve_create_composite

```c
qaws_status qaws_curve_create_composite(
    qaws_curve** out_curve,
    qaws_composite_desc const* desc);
```

Chains multiple heterogeneous curves end-to-end into a single composite curve.

Each segment occupies one unit in the composite parameter space: segment i maps to [i, i+1]. The total parameter range is [0, segment_count].

**Parameters:**
- `desc` -- Composite descriptor. Must not be NULL.
- `out_curve` -- Receives the created composite curve. Set to NULL on error.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- Ownership: the create function takes ownership of every curve pointer in the `segments` array. The caller must NOT destroy them afterwards; they will be destroyed when the composite curve is destroyed.
- Segments can be of different curve families (Bezier, Hermite, NURBS, etc.) as long as they share the same dimension.
- Parameter range is [0, segment_count].
