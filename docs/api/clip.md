# API Reference: qaws_clip.h

Boolean operations on regions of curves. This is Clipper2's `Clipper64` and
`ClipperD`, working on curves as well as polygons.

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

qaws_status qaws_clip_execute(qaws_clip_result** out_result, qaws_clip_desc const* desc);
qaws_status qaws_clip_boolean(out_result, subjects, clips, subject_count, clip_count, qaws_clip_type, qaws_fill_rule);
```

The paths are `qaws_path_2d` (see [path.md](path.md)) and may hold curves of
any kind.

- Subjects and clips are regions, read with the fill rule.
- A region path that does not close is closed by a straight line.
- Open subjects are clipped by the regions.

## How it works

1. **Hits.** One batched intersection of every input curve, including
   self-intersections, gives crossings, touches and shared stretches (see
   [curve_batch.md](curve_batch.md#hit-kinds)).
2. **Vertices.** Hits and curve ends within the tolerance are one vertex.
   - Lines whose ends lie within the tolerance of each other's line are
     collinear: their crossing would land anywhere along them, so they
     overlap instead. An end within the tolerance of another line's inside
     cuts it there.
   - A vertex within the tolerance of a line cuts the line there (snap
     rounding), so no line passes a vertex that clustering moved.
   - A tangent contact is only known to within about the square root of
     the hits' precision. One that near a curve's end is moved onto the end,
     where the curves meet at a corner.
3. **Edges.**
   - Each curve is cut at its vertices.
   - Edges joining the same vertices that run along each other become one
     edge, which carries the winding steps of all of them. A shared polygon
     edge is one such case, and so is an arc on a circle.
   - Edges whose steps cancel are dropped. Spikes and doubled-back
     stretches are such edges.
4. **Graph.**
   - The half-edges are sorted around each vertex by tangent angle. When
     tangents tie, the curve that bends left first goes first.
     - A tie is within 1e-5 rad, the noise of a tangent contact's angles.
       It is decided by probing each curve 1e-2 of the shorter edge along,
       so that bending beats rounding.
     - A direction along −x can come out at −π or at π by rounding, so
       angles that close to −π count as π. Tangent edges there then tie,
       instead of sorting to opposite ends of the fan.
   - A face is a cycle of "turn clockwise at the next vertex".
   - Each counter-clockwise cycle bounds a face. A clockwise cycle is the
     outside of a connected part, and belongs to the smallest face of
     another part that winds around it.
5. **Windings.** Faces get their subject and clip winding numbers by walking
   from the unbounded face (0, 0) across edges. A curve's area is closed to
   its vertices with chords, so a sliver's orientation is not lost to the
   gaps between the curve's ends and the vertices.
6. **Truth tables.** A face is inside the result by Clipper2's tables, with
   *in* meaning the fill rule holds:

   | Clip type | Inside |
   |---|---|
   | intersection | in subject and in clip |
   | union | in subject or in clip |
   | difference | in subject and not in clip |
   | xor | in exactly one |

7. **Output.**
   - The half-edges with an inside face on their left are walked into
     loops. At a vertex, the walk takes the next such half-edge clockwise.
   - Two parts touching at a point therefore come out as two paths, as in
     Clipper2.
   - Outer paths run counter-clockwise and holes clockwise.
8. **Open subjects.**
   - They are cut where they meet the regions.
   - Each piece is kept by where its middle lies: inside the clip for an
     intersection, outside both regions for a union, and outside the clip
     otherwise.
   - Crossings between open subjects do not cut them, as in Clipper2.

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

## Numbers

Test 84 checks these against closed forms:

- Two squares, with every clip type.
- Squares sharing an edge (one rectangle), touching at a corner (two
  paths), and identical.
- A pentagram under each fill rule.
- A square with a hole and an island, at depths 0, 1 and 2.
- Two circles: the lens area, in NURBS pieces.
- A circle minus a square.
- A line clipped by a square.
- The z callback.

`qaws_better2d_showcase` (figure 2) runs 12 operations on polygons, NURBS
ellipses and cubic B-spline blobs in 12 ms.

## Batches

```c
qaws_status qaws_clip_execute_batch(qaws_clip_result** out_results, qaws_clip_desc const* descs, unsigned int count, qaws_batch_executor const* executor);
qaws_status qaws_clip64_execute_batch(qaws_clip64_result** out_results, qaws_clip64_desc const* descs, unsigned int count, qaws_batch_executor const* executor);
qaws_status qaws_offset_execute_batch(qaws_clip_result** out_results, qaws_offset_desc const* descs, unsigned int count, qaws_batch_executor const* executor);
qaws_status qaws_offset64_execute_batch(qaws_clip64_result** out_results, qaws_offset64_desc const* descs, unsigned int count, qaws_batch_executor const* executor);
```

These run many independent jobs in one call, in parallel through the
executor.
- Each job runs serially; its own `executor` field is not used, because
  the jobs are the parallel part.
- The results come back in order, identical to one call each, whatever
  order the executor runs them in.
- If a job fails, its status is returned and no result is kept.

Test 95 runs 48 jobs of each kind serially, backwards, and interleaved,
and checks every result vertex for vertex.
