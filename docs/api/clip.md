# API Reference: qaws_clip.h

The 2D clip types and the clip result: regions of curves, their nesting
and their vertices. The Boolean operations that make a result are not part
of qaws; the types stay here so that libraries built on qaws share them.

```c
typedef enum qaws_clip_type { QAWS_CLIP_NONE, QAWS_CLIP_INTERSECTION, QAWS_CLIP_UNION,
                              QAWS_CLIP_DIFFERENCE, QAWS_CLIP_XOR } qaws_clip_type;

#define QAWS_CLIP_PRESERVE_COLLINEAR 1u   /* keep corners between collinear line pieces */
#define QAWS_CLIP_REVERSE_SOLUTION   2u   /* outer paths clockwise, holes counter-clockwise */

typedef struct qaws_clip_desc {
	qaws_path_2d const* subjects;      unsigned int subject_count;
	qaws_path_2d const* open_subjects; unsigned int open_subject_count;
	qaws_path_2d const* clips;         unsigned int clip_count;
	qaws_clip_type clip_type;
	qaws_fill_rule fill_rule;          /* EVEN_ODD, NON_ZERO, POSITIVE, NEGATIVE */
	unsigned int flags;
	qaws_scalar tolerance;             /* points this close are one vertex; 0 = 6.4e-9 of the extent */
	qaws_batch_executor const* executor;
	qaws_clip_z_fn z_fn; void* z_user; /* Clipper2's Z callback */
} qaws_clip_desc;
```

The paths are `qaws_path_2d` (see [path.md](path.md)) and may hold curves of
any kind.

## Result

```c
unsigned int qaws_clip_result_get_path_count(r);
qaws_status  qaws_clip_result_get_path(qaws_path_2d* out, r, i);   /* view, valid until destroy */
unsigned int qaws_clip_result_get_parent(r, i);                    /* QAWS_CLIP_NONE_INDEX at the top */
int          qaws_clip_result_is_hole(r, i);
unsigned int qaws_clip_result_get_depth(r, i);
unsigned int qaws_clip_result_get_open_path_count(r);
qaws_status  qaws_clip_result_get_open_path(qaws_path_2d* out, r, i);
qaws_status  qaws_clip_result_get_vertices(qaws_clip_vertex const** v, unsigned int* n, r, open, i);
void         qaws_clip_result_destroy(r);
```

**Output curves.** Every output curve is a new curve that lies exactly on
an input (see `qaws_curve_extract`):

- Consecutive pieces of one input curve are fused into one curve, so an
  input the operation does not touch comes back whole.
- Runs of line pieces become one polyline, a degree-1 B-spline. A
  polygon result is therefore one closed polyline. Unless
  `QAWS_CLIP_PRESERVE_COLLINEAR` is set, collinear corners are dropped.
- Curved pieces have their end control points placed on the shared
  vertices, so consecutive curves meet exactly.

**Nesting** is Clipper2's PolyTree, taken from the faces rather than from
point-in-polygon tests:

- A hole's parent is the outer path of the inside part on its left.
- An outer path's parent is the hole around the outside part on its right.
- Depth 0 is the top level. Holes are at odd depths.

**Vertices.** `qaws_clip_result_get_vertices` lists a path's vertices in
order: every point where two inputs meet, every curve junction and every
polyline corner. Each vertex records the input curves and parameters it
lies on. Where two inputs meet, it carries the z returned by `z_fn`.

