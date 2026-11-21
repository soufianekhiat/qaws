# API Reference: qaws_clip.h

The 2D clip types and the clip result: Boolean operations on regions of
curves in the style of Clipper2 (Clipper64 / ClipperD, on curves), their
output paths, nesting and vertices. The header declares the descriptor and
the result accessors; it declares no function that computes a result, so
the types are shared by code built on qaws that does.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)). Thread-safe on immutable curves. Paths
are `qaws_path_2d` ([path.md](path.md)) and may hold curves of any kind.

---

## Model

Subjects and clips are regions: closed paths read with the fill rule. Open
subjects are clipped by the closed regions. The result is new curves lying
exactly on the inputs (`qaws_curve_extract` pieces), in closed paths with
their nesting (Clipper2's PolyTree) and open paths.

The header describes the method: every input curve goes through one batched
intersection (crossings, touches, shared stretches); hits within tolerance
are one vertex; curves are cut at every vertex into edges, and edges that
run along each other are merged into one edge carrying the winding steps of
all of them. The edges form a planar graph whose faces get their subject and
clip winding numbers by walking across edges from the unbounded face. A face
is inside the result by the clip type on the fill rule of both:

| Clip type | Inside the result |
|---|---|
| intersection | in subject and in clip |
| union | in subject or in clip |
| difference | in subject and not in clip |
| xor | in exactly one |

The boundary between inside and outside faces, walked with the inside on
the left, gives the result: outer paths counter-clockwise, holes clockwise
(reversed with `QAWS_CLIP_REVERSE_SOLUTION`). Two parts touching at a point
come out as two paths. Nesting comes from the faces, not from
point-in-polygon tests.

## Types and flags

```c
typedef enum qaws_clip_type {
	QAWS_CLIP_NONE = 0,
	QAWS_CLIP_INTERSECTION,
	QAWS_CLIP_UNION,
	QAWS_CLIP_DIFFERENCE,
	QAWS_CLIP_XOR
} qaws_clip_type;

#define QAWS_CLIP_PRESERVE_COLLINEAR 1u   /* keep vertices between collinear line pieces */
#define QAWS_CLIP_REVERSE_SOLUTION   2u   /* outer paths clockwise, holes counter-clockwise */

#define QAWS_CLIP_NONE_INDEX 0xFFFFFFFFu
```

```c
/* Where an output vertex comes from: the input curve it lies on, and the
   second input curve at a crossing (operand 0 subject, 1 clip, 2 open
   subject; path and curve index into that operand's array). */
typedef struct qaws_clip_vertex {
	qaws_vec2 position;
	unsigned int operand_a, path_a, curve_a;
	qaws_scalar parameter_a;
	unsigned int operand_b, path_b, curve_b;   /* QAWS_CLIP_NONE_INDEX when on one curve only */
	qaws_scalar parameter_b;
	qaws_scalar z;                             /* from z_fn, else 0 */
} qaws_clip_vertex;

/* Called for every output vertex where two input curves meet (Clipper2's Z
   callback); returns the vertex's z. */
typedef qaws_scalar (*qaws_clip_z_fn)(void* user, qaws_clip_vertex const* vertex);

typedef struct qaws_clip_desc {
	qaws_path_2d const* subjects;
	unsigned int subject_count;
	qaws_path_2d const* open_subjects;
	unsigned int open_subject_count;
	qaws_path_2d const* clips;
	unsigned int clip_count;
	qaws_clip_type clip_type;
	qaws_fill_rule fill_rule;
	unsigned int flags;                    /* QAWS_CLIP_* */
	qaws_scalar tolerance;                 /* points this close are one vertex; 0 = 6.4e-9 of the extent (1e-4 in float) */
	qaws_batch_executor const* executor;   /* optional, for the batched intersection */
	qaws_clip_z_fn z_fn;                   /* optional */
	void* z_user;
} qaws_clip_desc;

typedef struct qaws_clip_result qaws_clip_result;   /* opaque */
```

## Result

```c
void         qaws_clip_result_destroy(qaws_clip_result* result);

unsigned int qaws_clip_result_get_path_count(qaws_clip_result const* result);
qaws_status  qaws_clip_result_get_path(qaws_path_2d* out_path, qaws_clip_result const* result, unsigned int index);
unsigned int qaws_clip_result_get_parent(qaws_clip_result const* result, unsigned int index);
int          qaws_clip_result_is_hole(qaws_clip_result const* result, unsigned int index);
unsigned int qaws_clip_result_get_depth(qaws_clip_result const* result, unsigned int index);

unsigned int qaws_clip_result_get_open_path_count(qaws_clip_result const* result);
qaws_status  qaws_clip_result_get_open_path(qaws_path_2d* out_path, qaws_clip_result const* result, unsigned int index);

qaws_status  qaws_clip_result_get_vertices(qaws_clip_vertex const** out_vertices, unsigned int* out_count,
	qaws_clip_result const* result, int open, unsigned int index);
```

- **Destroy** frees the result and every curve it owns. NULL is accepted.
- **Closed paths**: `qaws_clip_result_get_path` fills a view (`closed` = 1)
  whose curves stay valid until the result is destroyed.
- **Nesting** (PolyTree): `get_parent` is the enclosing path,
  `QAWS_CLIP_NONE_INDEX` at the top level; depth 0 is the top; holes are
  the paths at odd depth. A hole's parent is the outer path of the inside
  part around it; an outer path's parent is the hole it sits in.
- **Open paths**: pieces of the open subjects, read the same way
  (`closed` = 0).
- **Vertices**: the vertex at the start of each curve of a path; `open`
  selects an open path (1) or a closed one (0). Each records the input
  curves and parameters it lies on and, where two inputs meet, the z from
  `z_fn`. The array is owned by the result.
- Statuses: `QAWS_STATUS_INVALID_ARGUMENT` for a NULL pointer or an index
  out of range. The int and count accessors return 0 (and
  `get_parent` returns `QAWS_CLIP_NONE_INDEX`) for a NULL result or an index
  out of range.

**Output curves.** Consecutive pieces of one input curve are fused into one
curve, so an input the operation does not touch comes back whole. Runs of
line pieces become one polyline (a degree-1 B-spline), so a polygon result
is one closed polyline; unless `QAWS_CLIP_PRESERVE_COLLINEAR` is set,
collinear vertices are dropped. Curved pieces have their end control points
placed on the shared vertices, so consecutive curves meet exactly.
