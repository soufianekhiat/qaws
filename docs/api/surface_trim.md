# API Reference: qaws_surface_trim.h

Trimmed surfaces: a base surface restricted to a region bounded by trim
loops in its (u, v) parameter space. The result is a `qaws_surface` of kind
`QAWS_SURFACE_KIND_TRIMMED`, evaluated and destroyed with the generic
functions of [surface.md](surface.md).

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)).

---

## Types

```c
typedef struct qaws_trim_loop {
	qaws_curve const* const* curves;   /* 2D curves forming a closed loop (borrowed) */
	unsigned int curve_count;
	int is_outer;                      /* 1 = outer boundary, 0 = hole */
} qaws_trim_loop;

typedef struct qaws_surface_trim_desc {
	qaws_surface const* base;          /* base surface (borrowed) */
	qaws_trim_loop const* loops;       /* trim loops */
	unsigned int loop_count;
} qaws_surface_trim_desc;
```

- A trim loop is a closed sequence of 2D curves in parameter space.
- Points inside an outer loop are kept; points inside a hole are removed.
- With no outer loop, the full [0,1] x [0,1] domain is the outer boundary.
- The base surface and the trim curves are borrowed: they must outlive the
  trimmed surface.

## Functions

```c
qaws_status qaws_surface_create_trimmed(qaws_surface** out_surface, qaws_surface_trim_desc const* desc);
```

Creates the trimmed surface. Free it with `qaws_surface_destroy`; the base
and curves are not freed.

```c
qaws_status qaws_surface_trim_contains(int* out_inside, qaws_surface const* trimmed_surface,
	qaws_scalar u, qaws_scalar v);
```

Whether (u, v) lies in the kept region. Returns `QAWS_STATUS_OK` and sets
`*out_inside` to 1 (inside) or 0 (outside).

Trimmed surfaces are also the pieces produced by the surface shatter
([shatter.md](shatter.md)).
