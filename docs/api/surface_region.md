# API Reference: qaws_surface_region.h

Curves and regions on a parametric surface, traced by the surface's own
geometry (the plan: `docs/surface_regions_plan.md`). A geodesic is solved
on the surface. A geodesic circle is the set of points at one geodesic
distance from its centre. Both are kept as their images in the parameter
plane (u, v), which are generally curved. They are not projections of
planar shapes.

```c
qaws_surface_curve* c;
qaws_surface_curve_create_geodesic(surface, a_uv, b_uv, 0, &c);           /* shortest path a -> b */
qaws_surface_curve_create_geodesic_ray(surface, p_uv, dir_uv, len, 0, &c); /* from p along dir */
qaws_surface_curve_create_geodesic_arc(surface, centre_uv, r, a0, a1, 0, &c); /* geodesic circle arc */
qaws_surface_curve_create_uv(surface, uv_curve, &c);                       /* a curve drawn in (u, v) */

qaws_surface_curve_get_uv(c);              /* its (u, v) image: a qaws 2D curve */
qaws_surface_curve_get_length(c, &len);    /* surface length */
qaws_surface_curve_evaluate(c, t, &uv, &xyz);
qaws_surface_curve_destroy(c);

qaws_surface_get_metric(surface, u, v, g, christoffel);   /* E F G, G^k_ij */
qaws_surface_exp(surface, p_uv, dir_uv, len, 0, &q, &vel); /* the exponential map */
```

## The geometry

The surface S(u, v) carries the metric

```
ds^2 = E du^2 + 2 F du dv + G dv^2,   E = S_u.S_u, F = S_u.S_v, G = S_v.S_v,
```

and its geodesics solve, by arc length,

```
u'' = -(G1_11 u'^2 + 2 G1_12 u' v' + G1_22 v'^2),
v'' = -(G2_11 u'^2 + 2 G2_12 u' v' + G2_22 v'^2),
```

with `G^k_ij = g^-1 (S_u . S_ij, S_v . S_ij)` from the surface's first and
second derivatives. This works for every surface kind that evaluates
them.

## How curves are built

- **Integration**: the Dormand-Prince 5(4) pair with error control on the
  (u, v) position.
- **Representation**: a cubic Hermite spline in (u, v) through equal steps,
  with positions and velocities exact at the nodes. Geodesics step by arc
  length, circles by angle. The steps are halved until every span's
  midpoint is within the tolerance of the solution.
- **Geodesic between two points**: shooting from a. Newton runs on the
  leaving angle and the length. The angle's derivative comes from
  neighbouring shots; the length's is the end velocity. The start is the
  chord. When Newton fails from there, continuation walks the target along
  the chord. The result is locally shortest: on surfaces with several
  geodesics between the points, it is the one the chord leads to.
- **Geodesic circle**: theta -> exp_c(r (cos theta e1 + sin theta e2)).
  - e1 is along S_u, and e2 is a quarter turn from it about S_u x S_v.
  - Tangents come from a five-point stencil in theta.
  - A full turn closes the curve.
  - A circle that crosses itself (past the injectivity radius at its
    centre) is refused with `INVALID_ARGUMENT`.
- **Closed directions** (a cylinder's or a revolution's u): opposite edges
  of the domain that coincide in space are found at creation. Curves cross
  them continuously, with their (u, v) images in the covering plane (u past
  the domain where they went round).
- **Poles** (a whole edge of the domain mapped to one point, a sphere's
  ends): a geodesic circle centred on one has its rays start along the
  meridians u = u0 + P theta / (2 pi), a hair from the pole where the
  metric is regular. Its (u, v) image runs once round the period, against u
  at v0 so that it is counter-clockwise seen from outside.

## Accuracy and cost

- **Tolerance**: relative to the parameter domain's diagonal. It defaults
  to 1e-10 (1e-5 in float builds).
- **Test 101**, in f64:
  - against closed forms: a geodesic's length on a sphere matches the great
    circle to 3e-14, a geodesic circle's points and length 2 pi R sin(r/R)
    match to 1.4e-10, geodesics on a cylinder match the unrolled distances
    to 1e-15, and a flat non-affine patch gives straight geodesics to 4e-14;
  - against Mathematica (NDSolve) on a bicubic Bezier patch: the
    exponential map, a geodesic between two points and a geodesic circle
    agree to 3e-14.
- **Cost**: a geodesic across a sphere takes about 2 ms; a full circle with
  512 spans about 40 ms.

## Tests

Test 101 (`tests/101_surface_geodesic.c`), f64 and f32, with ground truth
from `tests/reference/101_surface_geodesic.wls`. Figure 16 in
`examples/better2d_showcase` draws geodesic polygons and circles on a
sphere, a bicubic patch and a torus, beside their parameter planes.

## Regions: Booleans and measures

```c
qaws_surface_curve const* edges[4] = { e0, e1, e2, e3 };
qaws_surface_path square = { edges, 4, 1 }, circle = { &c, 1, 1 };
qaws_surface_region* r;
qaws_surface_clip(surface, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, &square, 1, &circle, 1, &r);

qaws_surface_region_get_area(r, &area);          /* surface area */
qaws_surface_region_get_perimeter(r, &length);   /* boundary length on the surface */
qaws_surface_region_integrate(r, f, user, dim, out);   /* integral of f dA */
qaws_surface_region_contains(r, uv, &inside);    /* 1, 0, or -1 on the boundary */
qaws_surface_region_get_path(r, i, &uv_path);    /* the loops in (u, v) */
qaws_surface_region_create_trimmed(r, &trimmed); /* a trimmed surface */
qaws_surface_region_destroy(r);
```

- **Topology.** A patch maps (u, v) one to one onto the surface. So two
  curves cross on the surface exactly where their (u, v) images cross, and
  a point is inside a region exactly when its (u, v) is. `qaws_clip` on the
  images gives the regions' topology, under every clip type and fill rule.
  The curves are the true ones, traced intrinsically, so there is no
  projection error.
- **Area** is ∬ sqrt(EG − F²) du dv. Integrals ∬ f dA use the same
  construction, with any density f against the area element. Both go by
  Green's theorem: P(u, v) is the integral of the area element along u,
  with 8-point Gauss panels over each eighth of the domain. ∮ P dv is then
  taken round the loops, with 8-point Gauss on every half span of the
  boundary curves.
- **Perimeter**: ∫ sqrt(g(c', c')) dt along the loops.
- **Seams.** On surfaces closed in u or v (cylinders, revolutions, tori),
  the operands are brought into one period before the Boolean:
  - each loop is cut into its period cells (clipped by each cell under the
    fill rule of its own orientation), the pieces translated back by whole
    periods;
  - a loop going once round a closed direction (a latitude line, a pole
    circle) is capped along the domain's edge, and the cap's constant
    restored with the whole period;
  - each operand is united under its fill rule, then the Boolean runs
    under non-zero.
  Loops going round more than once, or round both directions of a torus,
  return `UNSUPPORTED_OPERATION`.
- **Poles**: the area element vanishes there, so areas and integrals over
  caps are ordinary; edges along the seam or on a pole are left out of the
  perimeter.

### Accuracy (test 102, f64)

| Check | Agreement |
|---|---|
| Sphere: geodesic circles' areas against 2 pi R² (1 − cos(r/R)) | 7e-11 |
| Sphere: perimeters against 2 pi R sin(r/R) | 1.1e-10 |
| Sphere: a geodesic triangle against Girard's R² (angles − pi) | 2.2e-10 |
| Sphere: union, intersection, difference, xor by inclusion and exclusion | 4e-13 |
| Flat sheared patch: a geodesic circle's area against pi r² | 4e-11 |
| Flat sheared patch: a union against the 2D engine on the plane's coordinates | 3e-11 |
| Cylinder: a geodesic circle against the unrolled pi r² | 5e-11 |
| Wavy bicubic patch: Gauss-Bonnet on a geodesic triangle (∬ K dA against the angle sum less pi) | 6e-12 |
| Sphere: a geodesic circle across the seam against the cap | 7e-11 |
| Sphere: circles centred on the poles, area and perimeter | 1e-15 |
| Sphere: caps and a band from latitude lines drawn in (u, v) | 9e-15 |
| Sphere: two circles across the seam by inclusion and exclusion | 2e-12 |
| Cylinder: a band between two loops round it, area and perimeter | 4e-15 |
| Cylinder: a circle across the seam against pi r² | 5e-11 |
| Torus: circles across either seam against the same circle away from it | 4e-13 |

In f32 the same checks hold to 6e-6.

A Boolean of two circles with their areas takes about 10 ms. Figure 17
shades Booleans on a patch and a sphere with their surface areas.
Figure 18 shows regions across the seams of a sphere, a cylinder and a
torus, and round the sphere's poles.

## Offsets

```c
qaws_surface_offset(surface, QAWS_FILL_NON_ZERO, paths, n, distance, 0, &r);
qaws_surface_curve_create_offset(curve, distance, 0, &offset_curve);
```

- **Offset curves.** Each point of the base is carried |distance| along
  the geodesic leaving it at a right angle: to the left of travel (about
  S_u x S_v) for a positive distance, to the right for a negative one.
  The points are fitted by a Hermite spline in (u, v), refined until the
  midpoints agree with the exponential map to the tolerance.
- **Regions.** The points within |d| of a region are reached either along
  a conormal geodesic from a smooth boundary point or from a corner. So the
  grown region (d > 0) is the region with:
  - the strips the outward conormal geodesics of length d sweep;
  - geodesic circles of radius d round the convex corners.
  The shrunk region (d < 0) is the region less the inward strips and
  circles round the reflex corners.
- **Corners.** A corner is where the boundary's pieces meet turning, with
  the turn measured in the metric. Only the corners on the side the offset
  moves to need a circle: elsewhere the two strips already overlap. Smooth
  joins need none. A piece ending at a cut along a seam is matched with the
  piece leaving the same point on the seam's other side, so a corner on a
  seam counts like any other. Leaving out the unneeded circles also keeps
  large offsets on open patches inside the patch: an outward circle at an
  inward offset's convex corner would leave the domain.
- **Boundary.** The strips follow the region's own boundary (its loops
  from `qaws_surface_clip`, inside on their left), not the input paths. So
  overlapping operands and holes are handled, and cuts along seams and
  poles are left out. The boundary is split into smooth pieces at its
  corners.
- **Open paths** give their tubes: strips on both sides, with circles at
  every vertex and round ends.
- **Topology.** The strips, circles and region are combined by
  `qaws_clip` in (u, v), cut into one period on closed surfaces, so the
  topology is exact. Folds of the strips past cusps (where the boundary
  curves more than 1/|d| towards the offset) fall inside the other strips
  and circles.

### Accuracy (test 103, f64)

| Check | Agreement |
|---|---|
| Sphere: a geodesic circle grown and shrunk against the circles of radius rho ± d (area, perimeter) | 1e-10 |
| Sphere: a geodesic triangle grown against A + L R sin(d/R) + (2 pi − A/R²) R² (1 − cos(d/R)) | 1.5e-10 |
| Sphere: a geodesic segment's tube against 2 l R sin(d/R) + 2 pi R² (1 − cos(d/R)) | 1.3e-10 |
| Sphere: the offset curve of a circle at the radius rho + d | 6e-11 |
| Sphere: two overlapping circles' union grown against the union of the larger circles | 9e-10 |
| Sphere: a geodesic triangle across the seam, a corner on it, grown against the closed form | 1.5e-10 |
| Flat sheared patch: a parallelogram grown, A + P d + pi d², and shrunk, A − P d + d² Σ cot(θ/2) | 3e-11 |
| Cylinder: a circle across the seam grown and shrunk against pi (rho ± d)² | 6e-11 |
| Cylinder: a band round it grown on both sides | 4e-15 |
| Wavy patch and torus: a geodesic circle's offset against the circle of radius rho ± d (Gauss lemma) | 2.5e-10 |

In f32 the same checks hold to 7e-5. Figure 19 shows offsets on a sphere,
a wavy patch and a torus.

## Sampling and meshing

```c
qaws_surface_sample_desc sd = { QAWS_REGION_SAMPLE_SOBOL, seed, NULL, NULL };
qaws_surface_region_sample(region, &sd, count, uv, points);

qaws_surface_mesh_desc md = { edge_length, 0, 0, 3 };
qaws_mesh_2d* mesh;
qaws_surface_region_triangulate(region, &md, &mesh);
```

- **Sampling.** The region's loops in (u, v) go through `qaws_region_sampler`
  (qaws_region_sample.h). Its density is the area element sqrt(EG − F²),
  times an optional density per unit area, so its mass-keeping map sends the
  unit square's patterns (random, Sobol, stratified, Halton, R2) to points
  uniform in surface area. The patterns keep their structure through the
  map.
  - **Blue noise** is sample elimination (Yuksel 2015) on five candidates a
    point, with distances measured in space.
  - Regions across seams and round poles sample like any other: the area
    element vanishes at a pole.
- **Meshing.** A constrained Delaunay mesh of the (u, v) loops
  (qaws_triangulate.h), refined by the size field L / sqrt(λmax(g)). λmax is
  the metric's largest eigenvalue, so no edge is longer than about L on the
  surface. The boundary tolerance is converted the same way. Botsch-Kobbelt
  remeshing rounds can follow.
  - On closed surfaces the cuts along seams are boundary edges of the
    (u, v) mesh; their vertices on the two sides are one point on the
    surface.
  - Near a pole the triangles grow long in u on the surface, since the
    size field is isotropic in (u, v).

### Accuracy (test 104, f64)

| Check | Agreement |
|---|---|
| Sphere, Sobol: mean cos(angle from the centre) over a circle across the seam, against (1 + cos(r/R)) / 2 | 2.5e-6 |
| Sphere, Sobol: z uniform over a band (Archimedes), mean and mean square | 2.3e-6 |
| Sphere, Sobol with the density 1 + z/R: the mean of z against its closed form | 2.2e-6 |
| Sphere, blue noise: 500 points' smallest distance against the packing radius | 1.28 times |
| Sphere: a circle across the seam meshed to 0.05: the triangles' area (relative), the longest edge | 1e-4, 0.0496 |
| Flat sheared patch: a parallelogram's mesh against its area | 4e-15 |
| Torus: a circle across both seams meshed to 0.03 (relative area) | 7.6e-5 |

The same checks pass in f32. Figure 20 shows the points and the meshes.

## Exact topology

```c
qaws_surface_clip_exact(surface, QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, a, na, b, nb, NULL, &r);
```

- **What is exact.** Every Boolean of the (u, v) images goes through
  `qaws_clip_exact` (qaws_clip_exact.h) on an integer lattice, so shared
  edges, touching corners, coincident loops and tangencies are decided
  exactly, not within a tolerance. This includes the cuts into period cells
  and the unions there on closed surfaces.
  - Curves drawn in (u, v) are taken exactly.
  - Geodesics are taken as their (u, v) Hermite images, which lie within the
    tolerance of the true curves.
- **The default lattice** is 2^-32 of the domain's extent plus a period, with
  48 parameter bits (2^-26 and 24 in float builds).
- **Moved, not cut.** A loop lying within one period cell is moved into the
  first cell without being cut.
- **Cuts at seams.** A curved boundary cut at a seam gives pieces whose ends
  are the doubles nearest the exact cuts. A piece's end a few ulps past its
  curve's start is that start on the kernel's parameter lattice; the kernel
  orders it there (it used to make a vertex just after the start, an edge of
  no length, and fail to certify). Curved and straight boundaries across
  seams certify.
- **Measures and queries** on the result are the same as for
  `qaws_surface_clip`.

### Checks (test 105)

- **Against the float path:** two geodesic circles' four Booleans, a circle
  moved by a period, and a circle on a torus give the same areas (1e-8).
- **Degenerate inputs drawn in (u, v):**
  - a region with itself: union and intersection are the region, one loop,
    and the difference is empty;
  - squares sharing an edge unite into one loop, with the area the sum to
    2e-16;
  - squares touching at a corner unite into two loops;
  - squares sharing part of an edge across the seam unite with the area the
    sum.
- **Circles across the seam:** five circles and their unions with shifted
  copies certify, with the float path's areas (1e-8).

Figure 21 shows a ring of sixteen squares sharing their sides united into
one band, checkerboards, and a lens.

## Medial axis

```c
qaws_surface_medial_desc md = { 0.04 };          /* resolution on the surface */
qaws_surface_medial_axis* m;
qaws_surface_region_medial_axis(region, &md, &m);
qaws_surface_medial_axis_get_segment(m, i, uv, radius);
qaws_surface_medial_axis_get_junction(m, i, &uv, &radius);
qaws_surface_medial_axis_get_max_inscribed(m, &centre, &radius);
```

- **What it is.** The axis is the set of points at the same geodesic
  distance from two or more nearest boundary points: the centres of the
  largest geodesic circles inside the region.
- **Boundary features.** The boundary is taken as its smooth pieces, split
  at corners, plus its corners. Cuts along seams and poles are not features,
  and their ends are not corners, since the boundary continues smoothly
  across them.
- **Labels.** Conormal geodesics are traced inward from samples of the
  pieces, with fans of geodesics from the corners. They label a mesh of the
  region, sized by the metric, with each vertex's nearest feature. Labels
  are then checked across edges by true distances.
- **Crossings.** A mesh edge whose ends have nearest boundary points far
  apart crosses the axis. The crossing is found on the edge where the two
  features' true geodesic distances agree: a piece's by Newton on its foot
  point, a corner's by shooting. A crossing nearer a third feature is off
  the axis and is dropped.
- **Junctions.** Within a triangle, two crossings are joined by a segment.
  Where three features meet, the junction (equally far from all three) is
  found by Newton.
- **Accuracy.** The axis is resolved to the resolution along it, and its
  points are on the true axis to the tolerance.

### Checks (test 106, f64)

| Check | Agreement |
|---|---|
| Sphere: a band between latitudes, the axis on the middle latitude at the half-width | 2.5e-11 |
| Sphere: a geodesic triangle's junction and largest circle against the spherical inradius | 1.3e-11 |
| Cylinder: a rectangle's axis, unrolled, on the plane rectangle's, radii the side distances | 3.5e-12 |
| Wavy patch: axis points inside the region shrunk by r − 1e-6, outside it shrunk by r + 1e-6 | 5 of 5 |
| Wavy patch: shrunk by the largest inscribed radius less 1e-6 a part left, plus 1e-6 nothing | yes |

In f32 the same checks hold to 3e-6. A test run takes about 14 s. Figure 22
shows axes, junctions and the largest inscribed circles on a sphere, a wavy
patch and a torus.
