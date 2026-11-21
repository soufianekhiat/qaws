# API Reference: qaws_shatter.h

Seams and shatter, in the sense of Houdini's Boolean "Seam" and "Shatter",
for curved surfaces in 3D and for 2D regions of curves.

- **Seams**: the curves (3D) or points (2D) where the operands cross.
- **Shatter**: operand A cut along its crossings with B into pieces, each
  told inside or outside B when B bounds a solid or region.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)). All functions are thread-safe on
immutable inputs. Concepts: [../shatter.md](../shatter.md).

---

## Seam kinds

```c
#define QAWS_SEAM_A_B 1u   /* an A element crosses a B element */
#define QAWS_SEAM_A_A 2u   /* two elements of A cross */
#define QAWS_SEAM_B_B 4u   /* two elements of B cross */
```

Bits for the `kinds` fields and the `kind` of each result.

---

## Surfaces (3D)

Operands A and B are lists of surfaces, trimmed or not: the faces of a body
or a sheet.

How it works: each surface is laid on a grid of `resolution`^2 cells over
its (u, v) domain, two triangles a cell. Triangles of surfaces that may meet
are paired by a grid of their boxes, each pair's crossing found, and the
crossings chained into curves by the triangle edges they cross. Every point
is then moved onto both true surfaces by Newton's method on
S1(u1, v1) = S2(u2, v2) (the least change of the four parameters, while the
gap shrinks). A point of a trimmed surface outside its trim is left out,
breaking the curve there.

Degenerate cases: exact coincidences between the operands (a grid vertex on
the other surface, a border lying in it, coplanar faces) are decided one
way for every triangle pair by moving each operand's grids a little (1e-8 of
the extent, 1e-5 in f32) along a direction of its own. A seam no bigger than
a hundred such moves (a touch at a point) is left out. A surface lying on
B's boundary is cut along B's edges there, its part on the boundary told
inside or outside as that move decides. Tangential contacts and crossings
finer than the grid are not resolved: raise `resolution`.

### Types

```c
typedef struct qaws_surface_seam_desc {
	qaws_surface const* const* a;   /* operand A's surfaces */
	unsigned int a_count;
	qaws_surface const* const* b;   /* operand B's surfaces (may be none for A_A seams) */
	unsigned int b_count;
	unsigned int kinds;             /* QAWS_SEAM_*: the seams wanted (shatter: A_B always) */
	unsigned int resolution;        /* grid cells along each side of a surface's domain; 0: 64 */
	qaws_scalar tolerance;          /* a point's gap between the two surfaces; 0: 1e-9 (f32: 1e-5) of the extent */
	qaws_batch_executor const* executor;   /* optional: the triangle pairs tested in parallel */
	qaws_progress const* progress;         /* optional: cancel and progress reports */
} qaws_surface_seam_desc;

typedef struct qaws_seam_point {
	qaws_vec3 position;
	qaws_scalar u1, v1;             /* on surface1 */
	qaws_scalar u2, v2;             /* on surface2 */
} qaws_seam_point;

typedef struct qaws_seam_info {
	unsigned int kind;                 /* QAWS_SEAM_* */
	unsigned int operand1, surface1;   /* operand 0: A, 1: B; the surface's index in it */
	unsigned int operand2, surface2;
	unsigned int point_count;
	int closed;
} qaws_seam_info;

typedef struct qaws_surface_seams qaws_surface_seams;       /* opaque */
typedef struct qaws_surface_shatter qaws_surface_shatter;   /* opaque */
```

### Seams

```c
qaws_status  qaws_surface_seams_create(qaws_surface_seams** out_seams, qaws_surface_seam_desc const* desc);
void         qaws_surface_seams_destroy(qaws_surface_seams* seams);
unsigned int qaws_surface_seams_get_count(qaws_surface_seams const* seams);
qaws_status  qaws_surface_seams_get(qaws_seam_info* out_info, qaws_seam_point* out_points,
	qaws_surface_seams const* seams, unsigned int index);
```

`qaws_surface_seams_create` finds the seams of the kinds asked:
`QAWS_SEAM_A_B` where a surface of A crosses one of B, `QAWS_SEAM_A_A`
(`QAWS_SEAM_B_B`) where two surfaces of A (of B) cross, not along a border
they share. A seam runs between ends (a surface's border, three surfaces
meeting) or is closed (its first point not repeated). surface1 is the lower
operand, then the lower index. Returns `QAWS_STATUS_INVALID_ARGUMENT` for no
surface, a NULL surface, or `kinds` 0.

`qaws_surface_seams_get` fills seam `index`'s description and writes its
`info.point_count` points to `out_points`, which may be NULL to read the
description only. Size the buffer from a first call with NULL.

### Shatter

```c
qaws_status  qaws_surface_shatter_create(qaws_surface_shatter** out_shatter, qaws_surface_seam_desc const* desc,
	int b_is_solid);
void         qaws_surface_shatter_destroy(qaws_surface_shatter* shatter);
unsigned int qaws_surface_shatter_get_count(qaws_surface_shatter const* shatter);
qaws_status  qaws_surface_shatter_get(qaws_surface const** out_piece, unsigned int* out_surface, int* out_inside,
	qaws_surface_shatter const* shatter, unsigned int index);
```

Every surface of A is cut along its A-B seams (`kinds` is taken as A_B)
into pieces. Each piece is a trimmed surface ([surface_trim.md](surface_trim.md))
over the A surface's base: its trim loops and the seams, flattened into
polygons in (u, v).

- `b_is_solid`: B's surfaces bound a solid with outward normals. Each piece
  is then told inside or outside B by the winding number of B's triangles
  at a point of it (inside from 1/2).
- `qaws_surface_shatter_create` returns `QAWS_STATUS_INVALID_ARGUMENT` as
  `qaws_surface_seams_create` does, or when B has no surface.
- `qaws_surface_shatter_get` gives piece `index`: the trimmed surface
  (owned by the shatter, do not destroy it), the index of the A surface it
  was cut from, and 1 inside B, 0 outside, -1 not told (B not a solid). Any
  output may be NULL.
- A piece is valid while the shatter and A's surfaces live.

---

## Regions (2D)

Operands are lists of paths ([path.md](path.md)). A is a region (closed
paths read with `a_rule`). B is a region (`b_rule`) or, with `b_is_region`
0, curves with no inside (open or closed paths), as a surface is in 3D.

How it works: every curve is flattened into chords within `flatness`, the
chords of both operands arranged, and every crossing moved onto both curves
by Newton's method. Collinear chords split one another where they overlap,
so boundaries running together (shared sides, identical regions) share
their edges and every operand's winding change across them counts.
Crossings closer than the flatness are not resolved.

### Types

```c
typedef struct qaws_region_seam_desc {
	qaws_path_2d const* a;
	unsigned int a_count;
	qaws_fill_rule a_rule;
	qaws_path_2d const* b;
	unsigned int b_count;
	qaws_fill_rule b_rule;
	int b_is_region;           /* 0: B's paths are curves with no inside */
	unsigned int kinds;        /* QAWS_SEAM_*: the seams wanted */
	qaws_scalar flatness;      /* the chords' deviation from the curves; 0: 1e-3 of the extent */
	qaws_scalar tolerance;     /* a crossing's gap between its curves; 0: 1e-12 (f32: 1e-6) of the extent */
	qaws_batch_executor const* executor;   /* optional: the chord pairs tested in parallel */
	qaws_progress const* progress;         /* optional: cancel and progress reports */
} qaws_region_seam_desc;

typedef struct qaws_region_seam_point {
	qaws_vec2 position;
	unsigned int kind;                         /* QAWS_SEAM_* */
	unsigned int operand1, path1, curve1;
	qaws_scalar t1;
	unsigned int operand2, path2, curve2;
	qaws_scalar t2;
} qaws_region_seam_point;

typedef struct qaws_region_shatter qaws_region_shatter;   /* opaque */
```

### Seams

```c
qaws_status qaws_region_seams_2d(qaws_region_seam_point* out_points, unsigned int* out_count,
	qaws_region_seam_desc const* desc, unsigned int capacity);
```

The points where a curve of A crosses one of B (`QAWS_SEAM_A_B`), or two
curves of one operand cross (`QAWS_SEAM_A_A`, `QAWS_SEAM_B_B`; not where a
path's curves join). Each carries its parameter on both curves; curve1 is
the lower operand, path and curve. A meeting where the curves run together
on either side (an overlap) is left out; one where they only touch (a
tangency, a corner on a side) is given.

- At most `capacity` points are written; `out_points` may be NULL to count.
  `*out_count` receives the number of points.
- `QAWS_STATUS_BUFFER_TOO_SMALL`: more than `capacity` points.
- `QAWS_STATUS_INVALID_ARGUMENT`: no path in A (or in B), a path without
  curves, `kinds` 0.

### Shatter

```c
qaws_status  qaws_region_shatter_create_2d(qaws_region_shatter** out_shatter, qaws_region_seam_desc const* desc);
void         qaws_region_shatter_destroy(qaws_region_shatter* shatter);
unsigned int qaws_region_shatter_get_count(qaws_region_shatter const* shatter);
qaws_status  qaws_region_shatter_get(qaws_path_2d const** out_paths, unsigned int* out_path_count, int* out_inside,
	qaws_region_shatter const* shatter, unsigned int index);
```

A's region cut along B's curves into pieces. A's own curves do not cut it
where its region lies on both sides (crossings of its paths inside it).

- Each piece is a region of closed paths, its outer path counter-clockwise
  and its holes clockwise. Its curves are exact pieces of A's and B's
  (`qaws_curve_extract`, reversed where walked backwards).
- When B is a region, each piece is told inside or outside B by winding
  numbers carried from face to face across the curves (each operand's fill
  rule on its own), not by sampling.
- `qaws_region_shatter_create_2d`: `QAWS_STATUS_INVALID_ARGUMENT` as
  `qaws_region_seams_2d` (`kinds` is ignored).
- `qaws_region_shatter_get` gives piece `index`: its paths (outer first,
  then holes) and their count, and 1 inside B, 0 outside, -1 not told. Any
  output may be NULL. The paths and their curves are owned by the shatter
  and live until it is destroyed.
