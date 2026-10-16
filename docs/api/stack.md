# API Reference: qaws_stack.h

2.5D shapes as stacks of parallel contours.

```c
typedef struct qaws_stack_level { qaws_scalar z; qaws_path_2d const* paths; unsigned int path_count; } qaws_stack_level;

typedef struct qaws_stack_2d {
	qaws_stack_level const* levels; unsigned int level_count;   /* increasing z, any spacing */
	qaws_fill_rule fill_rule;
	qaws_stack_interp interp;     /* QAWS_STACK_LINEAR, QAWS_STACK_CUBIC */
	qaws_stack_blend blend;       /* QAWS_STACK_BLEND_AUTO, _MATCH, _DISTANCE */
	unsigned int grid;            /* distance blend: marching squares cells across; 0 = 192 */
} qaws_stack_2d;

qaws_status qaws_stack_section_2d(qaws_clip_result** out, qaws_stack_2d const* stack, qaws_scalar z);
qaws_status qaws_stack_boolean_2d(qaws_stack_result** out, qaws_stack_2d const* a, qaws_stack_2d const* b,
	qaws_clip_type ct, qaws_batch_executor const* executor);
```

A stack is a list of levels at increasing heights, with any spacing. Each
level is a region: closed paths read with the fill rule. Between two levels
the shape is interpolated, and outside the first and last levels it is
empty. Examples:

- a cone: a circle at z = 0 and a point at the apex;
- a sphere: 64 circles;
- a height map: its contour regions { h ≥ z }.

## Sections

At a level's own height, the section is that level's region. Between
levels k and k + 1, with s = (z − z_k) / (z_k+1 − z_k), it is a blend:

**Matched contours** (`MATCH`). Contours are paired between the levels by
nearest centre and blended: C = (1 − s) C_k + s C_k+1.
- Contours of the same structure (kind, degree, knots, weights) blend
  their control points, homogeneous for rational kinds. This is exact: two
  concentric NURBS circles give the circle of the blended radius.
- A point contour (a path of one curve with no extent) takes the other
  contour's structure. A circle and a point therefore give a cone exactly:
  every section is a NURBS circle.
- Other pairs are resampled at 128 points of equal arc length. Their
  directions and starts are aligned, and the blend goes through a closed
  centripetal Catmull-Rom curve.

**Distance fields** (`DISTANCE`). The section is { (1 − s) d_k + s d_k+1 ≤ 0 },
where d is the signed distance to a level's region, negative inside.
- This handles contours that split, merge or appear between levels.
- The distances come from batched closest points on a grid. The sign comes
  from crossing counts with the fill rule.
- Marching squares traces the zero set, with saddles decided by the cell
  centre. The result is polylines.

**Automatic** (`AUTO`). This uses matching when the levels have the same
number of contours, each one curve, and distance fields otherwise.

**Cubic interpolation** (`CUBIC`) blends four levels with Hermite weights.
Their tangents are finite differences over the neighbouring levels, which
works for uneven spacing. Because both blend modes are linear in the
levels, the cubic weights serve either one. A sphere given as circles then
follows its curved profile rather than a faceted one.

## Booleans of stacks

The result has one level at every height where either stack has a level,
sorted, with each height once. Each level is the clip type applied to the
two sections at that height. A stack contributes nothing where it has no
level. The levels are independent jobs on the executor, and the results do
not depend on it.

```c
unsigned int qaws_stack_result_get_level_count(r);
qaws_scalar qaws_stack_result_get_z(r, i);
qaws_clip_result const* qaws_stack_result_get_level(r, i);
```

## Numbers

Test 89 checks:

- **A cone,** from a circle of radius 2 at z = 0 and a point at z = 4: the
  section at z = 1 has area 2.25π to 1e-9, and is a NURBS curve.
  - Above the apex the section is empty.
  - At a level, it is that level's region.
- **A sphere from 17 circles at uneven heights:** away from the poles, the
  cubic blend's radius error is under a quarter of the linear blend's.
- **An ellipse splitting into two circles:** the distance section at 0.9
  has 2 parts, with an area between the two levels'.
- **Cone × sphere:** 13 levels. At z = 1 the result is the closed-form
  lens of the two discs. The result is identical on a reversed executor.

Figure 5 of `qaws_better2d_showcase` draws these stacks.

## Surface meshes

```c
typedef struct qaws_stack_mesh_desc {
	qaws_stack_2d const* stack;
	qaws_scalar spacing;        /* target distance between contour samples; 0 = 1/64 of the extent */
	unsigned int sections;      /* extra sections between two levels; 0 = none */
	int caps;                   /* close the first and last sections */
	qaws_batch_executor const* executor;
} qaws_stack_mesh_desc;

qaws_status qaws_stack_mesh_3d(qaws_mesh_3d** out_mesh, qaws_stack_mesh_desc const* desc);
/* qaws_mesh_3d_get_vertex_count / get_vertex (qaws_vec3), get_triangle_count / get_triangle, destroy */
```

This builds the stack's surface as a 3D triangle mesh.

1. **Heights.** One section at every level, plus `sections` more between
   two levels so that a cubic profile is followed. The sections are
   computed in parallel on the executor.
2. **Samples.** Each contour is sampled on its curves, with its own count
   (its length over `spacing`, at least 8). Corners are always samples:
   where pieces meet at an angle or a polyline kinks (a turn over about
   25°). The arc length between corners is shared out evenly.
3. **Linking.** Between neighbouring heights, contours of the same
   orientation that overlap in plan are grouped:
   - **One to one:** a band.
   - **One to several (a split) or several to one (a merge):** the single
     contour's samples go to their nearest partner. Each run becomes a loop,
     closed across the neck by a bridge that the next loop shares. Each loop
     is banded to its partner. With three partners or more, the polygon of
     the cut points is filled.
   - **Anything else:** pairs by nearest centres.
4. **Bands.** Two contours with the same number of corners are stitched
   corner to corner, so creases run straight up the stack. The corner
   pairing is the rotation with the least total distance, since creases
   often come as close twins. Otherwise the upper contour starts at the
   sample nearest the lower one's start. In both cases the stitch is
   greedy, taking the shorter diagonal first, so the two sample counts
   need not match.
5. **Tips and ends.**
   - A contour left without a partner on one side closes in a fan to a
     tip when it is shrinking to a point: a cone's tip, a sphere's pole,
     a cone's tip inside a Boolean. The pace of shrinking comes from its
     partner on the other side. The tip sits where the size would reach
     zero, at most at the next height.
   - Otherwise it gets a flat cap at its own height, such as a cube's top
     face when the stack goes on above it, empty.
   - With `caps`, the first and last sections are closed by flat faces.
   - Caps are triangulated (`qaws_triangulate`) on the same samples. A
     cap edge that the triangulation merged across a straight run is split
     again through those samples, so caps meet bands edge for edge.

Triangles turn counter-clockwise seen from outside. The mesh is closed.

The levels of a stack Boolean (`qaws_stack_boolean_2d`) can be meshed
directly with `sections = 0`. Its contours merge, split and end level by
level.

Test 92 checks that every directed edge has its reverse exactly once, and
that the volume by the divergence theorem matches the closed form:

| Shape | Volume |
|---|---|
| cylinder | 2π |
| cone from a circle and a point | π/3 |
| tube with a hole | 3π |
| sphere from 17 circles with point poles, cubic | 4π/3, within 0.6% |
| box from two squares (corners kept) | 8 |
| two unit spheres united level by level (merge, then split) | the union's, within 3% |

The test also checks that an ellipse splitting into two circles stays
closed, and that an executor gives the same mesh. Figure 9 of
`qaws_better2d_showcase` renders nine such meshes on a 3x3 grid:

- cone, sphere, drifting tube;
- ellipse split, cone ∪ sphere, twisting star;
- cube − sphere, two spheres ∪, cube ∩ sphere.
