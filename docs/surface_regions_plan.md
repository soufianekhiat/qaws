# Regions on surfaces: the 2D toolkit on a curved surface, intrinsically

Branch `better_2D`.

Goal: the features of the 2D work run on a parametric surface, with shapes
defined by the surface's own geometry rather than projected onto it:
- Booleans and areas;
- geodesic offsets;
- sampling and meshing;
- skeletons.

A square's edges are geodesics, and a circle is the set of points at one
geodesic distance from its centre. Both are traced on the surface, and
every measure (length, area, distance) uses the surface's metric.

Decisions (2026-10-13):
- **Surfaces.** One parametric surface first (NURBS, B-spline, Bezier,
  revolution, extrusion, ...), with the seams and poles of closed ones.
  Regions crossing the faces of a B-rep come later.
- **Features, in order.** Booleans and area, geodesic offsets, sampling and
  meshing, skeleton and medial axis.
- **Primitives.** Geodesic polygons, geodesic circles and arcs, and any
  qaws 2D curve drawn in the parameter plane (u, v).
- **Paths.** Like every other feature: a float engine in f64 and f32 to a
  tolerance relative to the surface's size, and a certified counterpart in
  the exact kernel's style.

---

## 1. Why the parameter plane carries the topology

A patch S(u, v) is one-to-one away from its seams and poles. So two curves
on the surface cross exactly where their (u, v) images cross, and a point is
inside a region exactly when its (u, v) is. Once the curves are the true
ones, traced intrinsically, the 2D engine (`qaws_clip`) gives the exact
topology of the surface regions in (u, v). Only lengths, areas and
distances need the metric:

```
g = [E F; F G],  E = S_u . S_u,  F = S_u . S_v,  G = S_v . S_v,
ds^2 = E du^2 + 2 F du dv + G dv^2,  dA = sqrt(E G - F^2) du dv.
```

This is not a projection. Projecting a planar square onto a surface gives
the shadow of planar edges. Here the edges are the surface's own shortest
paths, solved on the surface, so their (u, v) images are generally curved.

## 2. Primitives

1. **Geodesic segment** between two surface points. It solves the geodesic
   equation in (u, v),
   ```
   u'' + G1_11 u'^2 + 2 G1_12 u' v' + G1_22 v'^2 = 0,
   v'' + G2_11 u'^2 + 2 G2_12 u' v' + G2_22 v'^2 = 0,
   ```
   with the Christoffel symbols from S's first and second derivatives.
   - Two points: shooting on the initial direction, by Newton with the
     Jacobi field (the variational equation) for the derivative. The start
     comes from today's mesh geodesics (`qaws_surface_compute_geodesic_flip`)
     when the plain straight-line guess fails.
   - Integration: an embedded Runge-Kutta pair (Dormand-Prince 8(5,3))
     with error control.
   - Stored as a (u, v) curve: a cubic Hermite spline through the steps
     (positions and velocities), refined until the curve is within the
     tolerance of the solution.
   - Reported: length, and the endpoint directions.
   - Shortest between its ends, or one of several geodesics: the caller can
     ask for the one leaving in a given direction.
2. **Geodesic circle and arc**: the exponential map of a tangent circle,
   theta -> exp_p(r (cos theta e1 + sin theta e2)).
   - Each point is a geodesic from p. The Jacobi field gives d/dtheta,
     hence the Hermite tangents.
   - Adaptive in theta.
   - It is the true set at distance r while r is below the injectivity
     radius at p (no conjugate or cut points). Beyond that the curve
     crosses itself, and is detected and reported.
3. **Curves drawn in (u, v)**: any qaws 2D curve, taken as given. Their
   images on the surface are generally not geodesics.
4. **Paths and regions**: closed loops of these, under a fill rule, as in 2D.

## 3. Seams and poles

- **Periodic directions** (a cylinder's or a revolution's u):
  - a curve is traced continuously in the covering plane (u unwrapped);
  - for the engine, a region is cut at the seam into copies in one period
    [u0, u0 + P];
  - each copy is clipped to the period strip, and the copies are united.
  - A loop that winds round the seam (a band round a cylinder) gives two
    boundary curves in the strip, as on the cylinder itself.
- **Poles** (a whole edge of the domain mapped to one point, a sphere's or a
  cone's apex): regions through or round a pole use a local chart there,
  geodesic polar coordinates or the tangent plane, and are mapped back.
  Containment of a pole is decided by the winding round it.
- **Singular points elsewhere** (S_u x S_v = 0 inside the domain) are
  reported, not handled.

## 4. Booleans and measures (first)

- `qaws_clip` on the (u, v) curves, after the seam cut, under every clip
  type and fill rule. The pieces of result edges keep their source
  (geodesic, circle or drawn curve, and its parameter), as in 2D.
- **Output:**
  - (u, v) paths;
  - a trimmed surface per region (`qaws_surface_create_trimmed`);
  - the boundary as 3D curves: S(c(t)), evaluated exactly, and fitted to
    a B-spline within the tolerance on request.
- **Area** ∬ sqrt(EG - F^2) du dv, by Green's theorem:
  - the inner integral along u is done with Gauss rules;
  - this is the same construction the region sampler uses for a density.
- **Perimeter** ∫ |S_u u' + S_v v'| dt.
- **Point containment** by the winding in (u, v), with the seam and pole
  rules.

## 5. Geodesic offsets (second)

- **The offset of a boundary curve** at distance d: each point moved along
  the geodesic leaving it at a right angle in the surface's tangent plane,
  for a length d (exp along the conormal).
  - Tangents from the Jacobi field.
  - An offset geodesic is not a geodesic. On a curved surface it is traced
    point by point, adaptively.
- **Joins** at convex corners:
  - round: geodesic circle arcs about the corner;
  - miter: the two edges' geodesic continuations, cut at the miter limit;
  - square and bevel: geodesic segments.
- **Clean-up**: a union under the positive fill rule, as in 2D.
- **Negative distances** shrink the region. Self-crossings are resolved by
  the union.
- **Checks**: on a sphere, the offset of a geodesic circle of radius r by
  d is the circle of radius r + d. On a cylinder, offsets equal the
  unrolled plane's.

## 6. Sampling and meshing (third)

- **Points uniform in surface area**: the region sampler on the (u, v)
  region with the density sqrt(EG - F^2), times a user density if given.
  The density sampling exists already.
- **Triangulation** (`qaws_triangulate`) in (u, v):
  - boundary samples are on the true curves;
  - refinement uses the metric: lengths, angles and areas measured on the
    surface, so the triangles are well shaped on the surface, not in
    (u, v);
  - remeshing keeps vertices on the surface and boundary vertices on their
    curves.

## 7. Skeleton and medial axis (fourth)

- **The medial axis** of a surface region: the points with two or more
  nearest boundary points in geodesic distance.
- **The wavefront** is the geodesic offset of section 5 run in time. Events
  are where offset fronts meet, found by tracking their (u, v) images with
  the generalised intersection event of `qaws_skeleton`.
- **Splice** at reflex corners: a growing geodesic circle, as in the curved
  2D skeleton.
- **Mitered (straight) skeletons** have no straight lines to run on. Their
  analogue moves each edge by its geodesic offset, with corners on the
  meeting of the offset fronts. Defined and implemented after the medial
  axis.
- **This is the research-heavy part.** Vertices travel on curves with no
  closed form, and every event is a root of three ODE-defined fronts.

## 8. The certified path

- **Curves drawn in (u, v)**: their crossings and the region topology are
  exact with `qaws_clip_exact` (rational or algebraic vertices, certified
  order), as in 2D. Only the lift to the surface is evaluated.
- **Geodesics and circles**: validated ODE enclosures (interval Taylor
  models with a Lohner-style wrapping) give a tube in (u, v) round each
  edge.
  - Crossings are certified when the tubes cross transversally and
    separate elsewhere.
  - The topology is certified when every crossing is.
  - Otherwise a refusal names the near-tangency, as `qaws_clip_exact`
    refuses equal curvatures.
- **Areas and lengths** come with rigorous error bounds (interval
  quadrature over the tubes).

## 9. API sketch

```c
/* qaws_surface_region.h */
typedef struct qaws_surface_curve qaws_surface_curve;    /* a curve on a surface: (u, v) image + source */

qaws_status qaws_surface_curve_create_uv(qaws_surface const* s, qaws_curve const* uv, qaws_surface_curve** out);
qaws_status qaws_surface_curve_create_geodesic(qaws_surface const* s, qaws_vec2 a_uv, qaws_vec2 b_uv,
    qaws_surface_geodesic_desc const* opt, qaws_surface_curve** out);
qaws_status qaws_surface_curve_create_geodesic_arc(qaws_surface const* s, qaws_vec2 centre_uv,
    qaws_scalar radius, qaws_scalar angle0, qaws_scalar angle1, qaws_surface_curve** out);
qaws_status qaws_surface_curve_length(qaws_surface_curve const* c, qaws_scalar* out);
qaws_status qaws_surface_curve_evaluate(qaws_surface_curve const* c, qaws_scalar t, qaws_vec2* uv, qaws_vec3* xyz);

typedef struct qaws_surface_path { qaws_surface_curve const* const* curves; unsigned int curve_count; int closed; } qaws_surface_path;

qaws_status qaws_surface_clip(qaws_surface const* s, qaws_clip_type type, qaws_fill_rule rule,
    qaws_surface_path const* subjects, unsigned int ns, qaws_surface_path const* clips, unsigned int nc,
    qaws_surface_region** out);
qaws_status qaws_surface_region_area(qaws_surface_region const* r, qaws_scalar* out);
qaws_status qaws_surface_region_contains(qaws_surface_region const* r, qaws_vec2 uv, int* out);
qaws_status qaws_surface_region_trimmed(qaws_surface_region const* r, qaws_surface** out);
qaws_status qaws_surface_offset(qaws_surface const* s, qaws_surface_path const* paths, unsigned int n,
    qaws_scalar distance, qaws_join_type join, qaws_surface_region** out);
qaws_status qaws_surface_region_sample(...);   /* through qaws_region_sampler with the metric density */
qaws_status qaws_surface_region_triangulate(...);
```

## 10. Tests and ground truth

- **Plane (a flat NURBS patch, affine and then curved in (u, v))**: every
  result equals the 2D engine's on the plane's own coordinates. Geodesics
  are straight in space, and areas are planar.
- **Cylinder and cone (developable)**: unrolled, they are a plane. Booleans,
  offsets and areas equal the 2D results on the unrolled shapes,
  including bands round the seam.
- **Sphere** (closed forms):
  - a geodesic triangle's area is R^2 (angle sum - pi);
  - a geodesic circle of radius r has perimeter 2 pi R sin(r / R) and area
    2 pi R^2 (1 - cos(r / R));
  - two circles' lens has a closed form;
  - poles inside regions.
- **Torus** (Gauss-Bonnet): ∮ k_g + ∬ K + Σ turning angles = 2 pi chi for
  every region, a check that needs no closed form.
- **Free-form NURBS**: Mathematica (`tests/reference/*.wls`, NDSolve on the
  geodesic equation) for geodesics, circles, crossings and areas.
- **Float and certified agree**: the certified enclosures hold the float
  results, and in f32 to its tolerance.

## 11. Roadmap

1. **Metric and geodesics.**
   - Christoffel symbols for every surface kind, from its derivatives.
   - The geodesic integrator with error control and the Jacobi field.
   - Shooting with the mesh start.
   - The (u, v) Hermite representation.
   - Geodesic circles and arcs with injectivity checks.
   - Tests: plane, cylinder, sphere, Mathematica.
2. **Surface curves and paths**: the curve-on-surface type, evaluation,
   length, 3D lift and fit.
3. **Booleans, area, perimeter, containment** on open patches: `qaws_clip`
   in (u, v), the metric area by Green's theorem, trimmed surfaces out.
4. **Seams and poles**: the cut and union at seams, and pole charts.
   Tests on the cylinder, torus, sphere and revolutions.
5. **Geodesic offsets**: exp along conormals, joins, clean-up.
6. **Sampling and meshing** in the metric.
7. **The certified path**: drawn curves exact through `qaws_clip_exact`,
   geodesics through validated enclosures.
8. **Skeleton and medial axis** on surfaces.
9. **Across B-rep faces**: curves continued over shared edges, regions
   spanning faces.

Each item lands with tests in f64 and f32, docs, and a showcase figure.

## 12. Status

- Items 1 to 3: done (tests 101 and 102, figures 16 and 17).
- Item 4: done.
  - Loops are cut into period cells and united there, and loops round a
    closed direction are capped.
  - Circles can be centred on poles.
  - Cut pieces are cleaned before the Booleans: control points within
    rounding of the seam lines go onto them, and end spans of no length
    (left by cuts a rounding from a knot) are dropped.
  - Test 102 covers the sphere, cylinder and torus; figure 18.
- Item 5: done. The offset strips follow the region's own boundary, split
  at corners with seam cuts left out. Test 103 checks closed forms on the
  sphere and a sheared plane, and the Gauss lemma on a free-form patch and
  a torus; figure 19.
  - Found on the way: `qaws_path_compute_area_2d` compared 8- and 16-point
    Gauss rules on polynomial spans, where 16 points are exact. Where
    rounding beat the tolerance, the spans subdivided to depth 24. On the
    B-splines that `qaws_clip` returns, this took 98 s for a 512-span
    circle; it now takes 2 ms.
- Item 6: done. Sampling goes through qaws_region_sampler with the area
  element as density (every pattern, plus a density), and blue noise by
  sample elimination in space. Meshes come from the (u, v) Delaunay mesher
  with the size field L / sqrt(λmax(g)). Test 104; figure 20.
