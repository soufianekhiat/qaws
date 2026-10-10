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

qaws_status qaws_stack_section_2d(qaws_stack_2d const* stack, qaws_scalar z, qaws_clip_result** out);
qaws_status qaws_stack_boolean_2d(qaws_clip_type ct, qaws_stack_2d const* a, qaws_stack_2d const* b,
	qaws_batch_executor const* executor, qaws_stack_result** out);
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
