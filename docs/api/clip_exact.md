# Exact Boolean operations on curves (`qaws_clip_exact.h`)

```c
typedef struct qaws_clip_exact_desc {
	qaws_path_2d const* subjects;  unsigned int subject_count;
	qaws_path_2d const* clips;     unsigned int clip_count;
	qaws_clip_type clip_type;      qaws_fill_rule fill_rule;
	unsigned int flags;                     /* QAWS_CLIP_REVERSE_SOLUTION, QAWS_CLIP_PRESERVE_COLLINEAR */
	qaws_exact_desc const* exact;           /* the lattice; NULL: qaws_exact_desc_default */
	qaws_batch_executor const* executor;
	qaws_clip_z_fn z_fn;  void* z_user;
} qaws_clip_exact_desc;

qaws_status qaws_clip_exact_execute(qaws_clip_exact_desc const* desc, qaws_clip_result** out_result,
	qaws_clip_exact_stats* out_stats);
```

`qaws_clip_exact_execute` runs the same arrangement as `qaws_clip`, but decides every topological
question on the exact kernel (`qaws_exact.h`), so no tolerance is involved. The result is an ordinary
`qaws_clip_result`: closed paths of new curves on the inputs, their nesting, and their vertices. You
read it with the `qaws_clip_result_*` functions.

## How it decides

The input curves are prepared onto one exact lattice. From there, each step is exact or certified.

1. **Hits.** One certified batch over all curves, self-intersections included:
   - crossings;
   - tangencies, with their multiplicity: the resultant is split into square-free factors (Yun's
     algorithm with modular gcds), and an even multiplicity means the curves touch;
   - overlaps, whose ends are exact rationals. Two spans on one algebraic curve map onto each other by
     a Möbius change of parameter, and the overlap is that image met with [0, 1].
2. **Vertices.** Hits are one vertex when that is proven:
   - exact rational points are compared exactly (segment crossings, dyadic parameters, overlap ends,
     curve ends);
   - an exact point is tested against the partner curve's implicit equation;
   - two algebraic points on one curve are tested by the gcd of the partners' implicit equations
     substituted into that curve, counting its roots in the common enclosure.
3. **Edges.** Each curve is cut at its vertices. Edges along an overlap are merged, using a union-find
   that tracks their relative direction, and they carry the sum of their winding steps. Edges whose
   steps cancel are dropped.
4. **Order round a vertex.** It comes from the exact sub-curve leaving the vertex, computed by integer
   de Casteljau at the rational parameter.
   - At an exact vertex, the first non-zero derivative gives exact signs. Curves leaving tangent are
     ordered by their curvature.
   - At an irrational vertex, the sub-curve over the parameter enclosure gives a cone of tangents
     (Sederberg and Wang). Neighbouring cones must be apart.
5. **Faces.** The cycles of the planar graph. Each part (connected component) lies in the innermost
   bounded face of the other parts that holds it. That face is found by a ray from an exact point of
   the part:
   - crossings are taken nearest first;
   - meeting another part from outside skips that part;
   - the first bounded face facing the start of the ray is the answer;
   - a ray that touches a curve, passes a vertex, or has two crossings too close to order is replaced
     by another one.
6. **Windings and output.** Windings come from the unbounded face across edges, then Clipper2's truth
   tables are applied on the fill rule. `qaws_clip`'s output stage does the rest: pieces of one source
   are fused, polylines joined, and nesting comes from the faces.

## What is not exact

- **Output geometry.** The output vertices and the ends of the new curves are the doubles nearest the
  exact points.
- **Cycle orientation.** It comes from areas in double. The topology does not depend on it.
- **Rays.** Their crossing distances are compared in double, with a margin. A ray that comes too close
  is given up, never guessed.

## When it refuses

`QAWS_STATUS_CERTIFICATION_FAILED` is returned when a decision cannot be proven:

- curves leave an irrational vertex tangent to each other;
- curves leave a vertex with the same curvature (a higher-order contact);
- two irrational vertices on one curve cannot be told apart by the gcd test;
- a contact passes through a singular point.

## Scope and statistics

Supported curves are those `qaws_exact_curve_prepare` handles. Polylines and composites are taken piece
by piece. Open subjects are not supported yet.

`qaws_clip_exact_stats` reports:

- the sources, hits, tangencies and overlaps;
- the vertices, edges, faces and parts;
- how many identity proofs ran;
- how many rays were given up.

## Tests and figure

**Test 97** (results agree with `qaws_clip`, the float engine):

- two squares in every clip type, with exact areas;
- squares sharing an edge, and squares touching at a corner;
- a parabola tangent to a square's side at an exact point, where curvature orders the vertex;
- three triangles with edges through (1/3, 1/3);
- a square with a hole and an island, nested to depths 0, 1 and 2;
- a pentagram under the even-odd and non-zero rules;
- crossing ellipses, with a duplicate of one.

**Test 85** runs all 195 records of Clipper2's `Polygons.txt` through the engine, on the integer
lattice:

- every record is certified;
- 121 records give `qaws_clip64`'s exact area;
- the other 74 contain slivers under one unit that `qaws_clip64` drops because they round away. Those
  agree with `qaws_clip` to within 5e-15.

**Figure 10** of `qaws_better2d_showcase` draws nine of these cases, one of them refused.

**Test 96** checks the exact contacts below the engine:

- dyadic and irrational tangencies;
- an inflection, of multiplicity 3;
- collinear overlaps, in the same and in opposite senses;
- overlapping parabola arcs.
