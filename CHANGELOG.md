# Changelog

All notable changes to Qaws are documented here.

## [Unreleased]

### Changed (3.0.0, breaking)

- **The parameter convention of Alwan, Suwar and Abaad** (CONTRIBUTING.md):
  outputs first, then reports, in-out buffers, the object handle, input
  buffers (strides beside them), sizes, descriptors, scalars, and the
  context, allocator or executor last. Constructors put the new object
  first: `qaws_curve_create_bezier(&curve, &desc)`. Every public, internal,
  mesh-module and shader-core declaration follows it; `tools/qaws_api_tool.c`
  checks it (CTest `qaws_api_convention`). Callers reorder their arguments;
  a type-checking compiler finds most call sites, but two arguments of the
  same type swap silently, so review each call.
- **Both precisions in one binary.** Every scalar-bearing function and type
  exists as `name_f32` and `name_f64`. `QAWS_SCALAR_IS_FLOAT` no longer
  selects the library; it selects which precision a caller's unsuffixed
  names, `qaws_scalar` and the scalar macros mean (default 1, as before).
  `qaws_both.h` declares both in one translation unit; `qaws_f32` and
  `qaws_f64` are the scalar types. The allocator (`qaws_allocator.h`: the
  descriptor, the default, thread and scratch allocators, the arena) is
  precision-free and exists once: an allocator installed once serves both
  precisions. The amalgamation is compiled twice (through two generated
  wrappers), and `QAWS_UNITY_BUILD` builds one translation unit per
  precision.

### Changed

- Self-intersections and Booleans scale: pairs from a grid fitted to the triangles
  (`qaws_internal_grid_create_fit`; the volume rule made a surface's cells hold
  dozens of boxes), filtered read-only through `qaws_intersect_desc.executor` and
  cut in a fixed order (results independent of threads), neighbours rejected by
  their free corners' sides, edge records made only for faces that meet
  something, the check limited to pieces of different faces. A 522k-triangle
  resolution: 4.05 s to 0.77 s on one thread, 0.55 s on 8. The grid walk reads
  each cell's boxes once and compares cell indices instead of coordinates; the
  repair's intersection search uses the fitted grid.

### Fixed

- Mesh Booleans: a self-union of parts overlapping in a plane (copies facing
  the same way or folded back) kept every copy, leaving the result open; the
  lowest input face's pieces now stand for the overlap.
- Remeshing against a reference with input feature lines that are not its
  creases (cut isolines against the uncut mesh): the lines were snapped
  toward the reference's nearest crease (the border), splitting without end
  (13k faces where 2k were due). They are now kept as lines across the
  surface.
- Self-intersection resolution: points where a face's edge crosses another's
  plane were made before the pair was known to meet, splitting edges beside
  the other face; coplanar neighbours touching at a corner or along an edge
  were clipped against each other within the tolerance. A clean mesh now comes
  back unchanged, and a chain of Booleans no longer grows stray vertices (the
  CSG emblem of figure 27: 3784 faces instead of 7918). Boolean results have the
  edges inside flat regions turned toward Delaunay.

### Added

- Offsets of meshes (`qaws_mesh_3d_offset`, `qaws_mesh_offset.h`): the surface
  at a signed distance, shells and thickened sheets, by contouring the exact
  distance (dual contouring, shared with the watertight surfaces) or by
  moving the vertices along mitred directions with the self-intersections
  resolved by the self-union. Test 171, figure 39 (key D).
- Repulsive curves on a mesh (`qaws_repulsive_mesh_curves.h`): curve
  vertices as locations on a mesh's triangles, the tangent-point energy with
  each pair's distance the larger of the chord and the geodesic (heat
  method), steps walked across the triangles. Test 170, figures 37 and 38
  (keys B, C).
- The mesh connector (`qaws_mesh_3d_connect`): two border loops (on two
  meshes or one) joined by a bridge matched by arc length, Hermite-blended
  along both surfaces (tangent to both), rows and counts graded, quads or
  zipped triangles, optional smoothing; chains of borders (strips), the
  attributes and features carried, the bridge measured by the quality
  report; planar meshes (`qaws_mesh_2d_connect`); volume patches joined by
  prisms and hexahedra (`qaws_volume_mesh_connect`). Test 169, figures 35
  and 36 (keys z, A).
- Metric CVT: `QAWS_REMESH_CVT` takes `qaws_remesh_desc.metric` (it was
  UNSUPPORTED_OPERATION): an anisotropic CVT with each Voronoi cell
  measured in its seed's metric (Du and Wang 2005), density sqrt det M, the
  count from the metric area, Delaunay flips in the metric. On a square at
  0.02 by 0.2, edges 1.009 +- 0.076 in the metric (LOCAL 1.059 +- 0.122).
  Test 158.
- Remeshing from attributes: `qaws_mesh_3d_cut_isolines` (an attribute's
  isolines cut in as feature lines, kept through remeshing),
  `qaws_metric_field_create_intersection` (several attributes' metrics
  made one) and `qaws_metric_field_create_graded` (H-gradation). Test 168,
  figure 34 (key y).
- Decimation by quadric error metrics (`qaws_mesh_3d_decimate`): a target
  count or an error bound, border and feature lines held by constraint
  planes, feature corners kept, folds refused, attributes carried, the
  certified deviation from the input on request; creases found by angle,
  silhouettes kept from given views (planes through the view rays) or
  weighted from octahedral views (faces seen, outlines), per-vertex
  importance, a bound on the faces' turn. Test 167, figures 32 and 33
  (keys w, x).
- Deformation and parameterization (`qaws_mesh_3d_deform`,
  `qaws_mesh_3d_parameterize`): ARAP, Laplacian and biharmonic editing by
  handles in a region; LSCM and ARAP flattening with seams cut open, charts
  laid side by side, distortion reported (BFF reserved). The smoothing
  module's sparse Laplacian and factored systems shared. Test 166, figure 30
  (key u).
- Smoothing and fairing (`qaws_mesh_3d_smooth`, `qaws_mesh_3d_fair`):
  Laplacian, Taubin, implicit mean curvature and bi-Laplacian smoothing,
  cotangent or uniform weights, border, held vertices and feature corners
  kept, feature lines smoothed along themselves, the volume kept on request;
  fairing of free vertices as a membrane, thin plate or curvature-continuous
  surface. Test 165, figure 29 (key t).
- Geodesics on meshes (`qaws_mesh_geodesic`): distances by the heat method
  (Crane et al. 2013) on the intrinsic Delaunay triangulation (stretched and
  obtuse triangles handled), its matrices factored once by the sparse LDL^T so each
  query is a few solves; sources anywhere on the triangles; Neumann,
  Dirichlet or mean border conditions; paths down the distance shortened
  through the triangles they cross (unfolded, funnel algorithm). Test 164,
  figure 28 (key s).
- Mesh Booleans (`qaws_mesh_boolean`): union, intersection, difference and
  self-union of closed triangle meshes, on the resolved arrangement: patches
  bounded by the cuts, each sampled once in front of and behind it by the
  operands' winding numbers, kept and turned outward where they separate the
  result; coplanar overlaps kept once. Test 163, meshing figure 27. The repair's
  exact pair test skips shared corners (no exact zero ground out), and the
  resolution checks only the cut faces' pieces.
- Self-intersections resolved numerically (`qaws_mesh_resolve_intersections`):
  crossing triangles cut along their intersections into a conforming arrangement,
  points keyed by what made them (edge-plane, three planes, edge-edge) and snapped
  to the input's features within a tolerance, coplanar overlaps cut along each
  other's edges, cut faces triangulated again (constrained Delaunay), the result
  checked by exact predicates. Test 162, meshing figure 26.
- Watertight surfaces from any triangle input (`qaws_mesh_make_watertight`): the
  generalized winding number's inside (or within an offset) contoured on a grid by
  dual contouring with Hermite data from rays against the input (sharp features
  kept), holes first filled tangent to the surface (faired fills), overlaps
  united, sheets wrapped. Test 161, meshing figure 25.
- Repair completed: holes filled along creases (crease ends joined across the hole,
  corners recovered, each part faired on its own), long holes halved by short
  diagonals, unused vertices reported, small closed parts removed by volume,
  non-manifold edges and vertices split into manifold sheets
  (`QAWS_MESH_REPAIR_SPLIT_NON_MANIFOLD`), self-intersections found with exact
  predicates and marked (`QAWS_MESH_REPAIR_FIND_INTERSECTIONS`). Test 139,
  meshing figure 24.
- Quality reports (`qaws_quality_report.h`): angles, dihedral angles, shape (mean
  ratio), aspect ratio, scaled Jacobian, edge ratio and size per element of 2D,
  surface and volume meshes (inverted elements negative), and a summary
  (statistics, histogram, worst elements, threshold counts). Test 160, meshing
  figure 23.
- Spatial queries on meshes (`qaws_mesh_bvh.h`): a BVH over segments, triangles or
  tetrahedra (binned SAH), from arrays, soups, `qaws_mesh_3d`, `qaws_mesh_2d` or
  `qaws_volume_mesh`; closest points and point location, watertight rays (closest
  or any hit), generalized winding numbers (second-order far field), signed
  distances signed by them, refit; batches through the executor. Test 159,
  meshing figure 22.
- Single-file amalgamation: `buildsystem/amalgamate.py` collapses the tree into
  one `qaws.h` and one `qaws.c`, for consumers with no build system. Exposed as
  the CMake `amalgamate` target and the `QAWS_AMALGAMATION` /
  `QAWS_TEST_AMALGAMATION` options, and as the Sharpmake `QawsAmalgam` project.
- `QAWS_UNITY_BUILD` CMake option to compile the library as one translation unit.
- `examples/amalgam_standalone.c`, a worked no-build-system example.
- `docs/integration.md` covering integration without CMake, including Xcode/iOS
  notes on SIMD dispatch.
- `.github/workflows/permutations.yml` builds every combination of OS, scalar
  type, SIMD, and translation-unit layout, and compiles the amalgamation
  standalone on all three platforms.
- Boundary layers: `qaws_volume_mesh_grow_layers` (`qaws_boundary_layers.h`)
  grows layers of prisms and hexahedra from a surface's triangles and quads,
  inward or outward, at a first height and a growth ratio, thinned where
  they would meet the surface or each other; every corner's Jacobian
  positive (test 156).
- Curved elements: `qaws_mesh_generate_nodes_surface` / `_shell` (with
  `qaws_mesh_generate_node_counts`) give a generated mesh's nodes at order 2
  or 3 (P2, P3 triangles; Q2, Q3 quads) on the exact geometry: order - 1 per
  edge shared by both faces, a cubic triangle's centroid, a quad's
  (order - 1)^2 inside, each snapped onto its element's face. A torus at
  0.4 comes 50 to 140 times closer to the surface than its straight
  elements (test 157, showcase figure 20).
- Metric remeshing: `qaws_remesh_desc.metric` remeshes (LOCAL) to an
  anisotropic metric field (as Mmgs): every length, angle and shape it
  weighs measured in the metric, each vertex's taken in its tangent plane.
  Metric fields (`qaws_mesh_field.h`) are given per reference vertex
  (`qaws_metric_field_create_vertices`) or built from a solution's values
  for its linear interpolation to stay within an error
  (`qaws_metric_field_create_hessian`: the Hessian fitted on two rings,
  M = (2/9) |H| / error, Frey and Alauzet 2005). A tanh layer on a 64 x 64
  grid comes out in 133 vertices, interpolated within 1.3 times the error,
  its triangles six times longer along the layer than across (isotropic at
  the smallest size, the layer alone would take about 2300); a constant
  metric gives every edge within [1/sqrt 2, sqrt 2] of unit length (test
  158, showcase figure 21).
- Halide evaluates the curve and surface records (`qaws_curve_record_core.h`,
  `qaws_surface_record_core.h`) with the knot span found at run time, an
  `Expr` per element; `qaws_index` (an `int`, an `Expr` in Halide) for the
  indices the data decide. `examples/eval_record_halide.cpp` checks them
  against the C cores by hand (1.5e-15 relative).

### Changed

- LOCAL remeshing collapses along straight feature lines and borders: the
  feature edge a collapse carries on was refused as a chord over the crease
  (its midpoint on it), so a straight border never coarsened. A 64 x 64
  grid of the unit square remeshed to 0.2 now has 42 vertices (338 before);
  generated boxes come out coarser along their edges.

- Deduplicated `compute_normal` (16 identical copies) and `catmull_rom_blend`
  (2 identical copies) across the surface sources into
  `internal/qaws_internal_surface.{h,c}`. This removes the only symbol
  collisions that prevented single-translation-unit builds.
- CMake compiler settings moved into a shared `qaws_configure_target` function
  so the per-file and amalgamated libraries cannot drift apart.
- Shell regions: the trims' nearest-point and inside tests skip the runs of
  their polygons whose boxes cannot matter (results unchanged, bit for bit):
  trimmed shells about 2.5 times faster (test 108: 133 s to 45 s).
- Tests: a fast default tier and `QAWS_TEST_HEAVY=1` for the slow checks
  (docs/building.md); the runner reports each suite's time. The default run
  takes 360 s serially instead of 783 s.
- Tests: `TEST_BOUND(bound, n)`, a double bound that rounding sets in float
  (at least n roundings of the scalar); the mesh suites 114 and 122 to 139
  pass with `QAWS_SCALAR_IS_FLOAT` but for five quality checks against
  baselines measured in double (126 box, 134 torus, 135 blob, 136 can, 137
  density, stroke and attractor), which the f64 build fails too with its
  input rounded to float: the quad meshers are sensitive to input noise of
  1e-7 (open).
- Sharpmake: the test project builds again at /W4 /WX (Release f32 and f64):
  tests 93 and 94 read their data through an undefined `QAWS_TEST_DATA_DIR`
  (now defaulted in test_common.h), variables possibly used uninitialized,
  shadowed locals, 32-bit shifts widened; `QAWS_UNUSED` for the core
  headers; constant conditions (the precision branches) and the SDK's
  winbase.h warning disabled for the tests only.

## [1.0.0] - 2026-03-12

First public release.

### Curve families (13)

- Bezier (linear through arbitrary degree, 2D/3D)
- Hermite (cubic, 2D/3D)
- Catmull-Rom (uniform, chordal, centripetal parameterization, open/closed, 2D/3D)
- B-Spline (arbitrary degree, uniform/custom knots, open/closed, 2D/3D)
- NURBS (arbitrary degree, weighted, open/closed, 2D/3D)
- Trajectory (time-keyed cubic, optional velocity/acceleration, open/closed, 2D/3D)
- Yuksel C2 (Bezier, Circular, Elliptical, Hybrid modes, open/closed, 2D/3D)
- Rational Bezier (weighted De Casteljau, 2D/3D)
- Arc (circular/elliptical, multi-segment, 2D/3D)
- Polynomial (monomial form, Horner evaluation, 2D/3D)
- Clothoid (Euler spiral, linear curvature, 2D)
- Subdivision (Chaikin, Lane-Riesenfeld 3/4, open/closed, 2D/3D)
- Composite (heterogeneous multi-segment, 2D/3D)

### Surface types (5)

- Bezier patch (tensor product)
- B-Spline surface
- NURBS surface
- Swept surface (profile along path)
- Ruled surface (linear blend of two curves)

### Evaluation

- Position, D1, D2, D3 derivatives
- Per-span evaluation
- Batch evaluation (SIMD-accelerated: AVX2, SSE2, NEON)
- Surface evaluation with partial derivatives and normals

### Sampling

- Uniform parameter sampling
- Adaptive (error-tolerance based)
- Curvature-weighted
- Feature-preserving
- Streaming (callback-based)

### Traversal

- Parameter, arc-length, and timed modes
- Motion profiles: constant speed, constant acceleration, trapezoidal, S-curve, custom
- 10 easing functions (quad, cubic, sine -- in/out/in-out)
- Wrap modes: clamp, loop, ping-pong
- Multi-curve traversal

### Inspection

- Arc length, bounding box, closest point
- Tangent, normal, curvature, torsion, speed
- Frenet frame (3D)
- Inflection points, extrema
- Winding number (2D)
- Curvature comb (2D/3D)
- Curve-curve intersection, self-intersection

### Operations

- Split at parameter
- Join two curves
- Offset / parallel curve (2D, with self-intersection cleanup)
- Reverse
- Arc-length reparameterization

### Conversion

- Hermite -> Bezier, Catmull-Rom -> Bezier
- Bezier -> B-Spline, B-Spline -> NURBS
- Degree elevation, degree reduction (Bezier)

### Import / export

- SVG path data export (exact for Bezier/Hermite/Catmull-Rom, approximate for others)
- Polyline export (uniform / curvature-weighted, 2D/3D)
- Polyline import as Catmull-Rom or Trajectory
- B-spline fitting (least-squares)

### Multi-backend

- C (full runtime, SIMD batch eval, float/double)
- HLSL (22 core eval functions, SM 5.0+)
- GLSL (22 core eval functions, Vulkan/OpenGL)
- Halide (22 core eval functions, symbolic IR)
- Automatic backend detection via compiler macros
- Prepare functions for CPU-to-GPU coefficient upload

### Infrastructure

- Zero dependencies (C11 standard library only)
- Immutable, thread-safe curves and surfaces
- Custom allocator support (`_ex` variants)
- Stack-allocated inline curves (no malloc, degree <= 7)
- CMake install rules, pkg-config, CMake config export
- Sharpmake build system support
- MIT license
