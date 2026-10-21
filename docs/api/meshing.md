# API Reference: meshing

Meshes as editable objects, the surfaces they sample, and the operators
that build and rebuild them: remeshing, quadrangulation, generation from
exact geometry, tetrahedra, hexahedra, direction and frame fields. Data
enters and leaves as arrays; qaws reads and writes no files.

| Header | What |
|---|---|
| `qaws_mesh.h` | `qaws_mesh_3d`: polygon meshes in space, edits, features, attributes, footpoints, quality |
| `qaws_triangulate.h` | `qaws_mesh_2d`: planar meshes, the same store at z = 0 ([triangulate](triangulate.md)) |
| `qaws_mesh_reference.h` | `qaws_mesh_reference`: the smooth surface a mesh with normals implies (PN) |
| `qaws_mesh_field.h` | size fields; metric fields; cross fields on surfaces; frame fields in volumes |
| `qaws_remesh.h` | surface remeshing (LOCAL, CVT), isotropic or to a metric |
| `qaws_quadrangulate.h` | quad meshes (INSTANT) |
| `qaws_mesh_generate.h` | meshes of surfaces and shells (REMESH); curved elements (P2, P3, Q2, Q3 nodes) |
| `qaws_mesh_repair.h` | repair for farq: soups and broken meshes made closed and outward |
| `qaws_mesh_watertight.h` | watertight surfaces from any input: winding-number inside, dual contouring, offsets |
| `qaws_mesh_intersect.h` | self-intersections resolved numerically: a conforming arrangement of the input's triangles |
| `qaws_mesh_bvh.h` | spatial queries: closest points, point location, rays, winding numbers, signed distance |
| `qaws_quality_report.h` | quality reports: a measure per element of 2D, surface or volume meshes, histograms, worst elements, thresholds |
| `qaws_volume_mesh.h` | `qaws_volume_mesh`: tetrahedra, hexahedra, prisms, pyramids |
| `qaws_tetrahedralize.h` | tetrahedra filling a closed surface (DELAUNAY) |
| `qaws_hexahedralize.h` | hexahedra filling a closed surface (GRID) |
| `qaws_boundary_layers.h` | boundary layers: prisms and hexahedra grown from a surface |

Every operator has a `method` field. The methods listed as reserved return
`QAWS_STATUS_UNSUPPORTED_OPERATION` until they land. Operators never change
their input and return a new object; each fills a `qaws_mesh_report` (what
it did and reached; a target not reached is not an error). The design and
its decisions are in [meshing_plan.md](../meshing_plan.md).

## Surface meshes: qaws_mesh_3d

```c
typedef struct qaws_mesh_3d_desc {
	qaws_scalar const* positions; unsigned int vertex_count;     /* 3 per vertex */
	unsigned int const* face_vertices; unsigned int const* face_sizes; unsigned int face_count;   /* sizes NULL: triangles */
} qaws_mesh_3d_desc;

qaws_status qaws_mesh_3d_create(qaws_mesh_3d** out_mesh, qaws_mesh_3d_desc const* desc);
qaws_status qaws_mesh_3d_get_faces(unsigned int* out_face_vertices, unsigned int* out_face_sizes, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_get_face_adjacency(unsigned int* out_faces, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_split_edge(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position, qaws_mesh_location const* footpoint,
	unsigned int a, unsigned int b);
qaws_status qaws_mesh_3d_collapse_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position,
	qaws_mesh_location const* footpoint);
qaws_status qaws_mesh_3d_flip_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b);
qaws_status qaws_mesh_3d_split_face(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position,
	qaws_mesh_location const* footpoint, unsigned int face);
```

The storage is geogram's facet corners: a face is a run of corners, each
holding its vertex and the opposite corner across its edge (the corner
table). A corner is a half-edge, so the faces' adjacency, the dual graph,
is stored rather than derived (`qaws_mesh_3d_get_face_adjacency`,
`qaws_mesh_3d_create_dual`). Faces turn counter-clockwise seen from
outside; indices stay dense through edits.

- **Creation** never refuses non-manifold input: a face whose edge is
  already used in the same direction gets copies of that edge's vertices,
  and a vertex shared by several fans of faces gets one copy per extra fan
  (`qaws_mesh_3d_get_vertex_source` names the original).
- **Edits** check their own validity (`qaws_mesh_3d_can_collapse_edge`,
  `..._can_flip_edge`): a collapse needs the link condition, a flip a
  convex quad. Refused edits return `QAWS_STATUS_INVALID_TOPOLOGY` and
  leave the mesh unchanged.
- **Features**: feature edges (creases, the border, seams) and feature
  corners, set or detected by dihedral angle; edits keep them.
- **Attributes** per vertex, per corner (a vertex seen from one face: UVs,
  split normals) or per face, each with a rule (`LINEAR`, `NORMALIZED`,
  `NEAREST`) that edits and operators follow.
- **Footpoints**: each vertex's location on the reference an operator
  worked on, kept through every operator (`qaws_mesh_3d_get_footpoints`).
- **Measures**: `qaws_mesh_3d_measure_quality` (faces by size, angles,
  edges, aspect, valences), `qaws_mesh_3d_get_topology` (components,
  boundary loops, Euler characteristic, genus).

## References: the surface between the vertices

```c
qaws_status qaws_mesh_reference_create(qaws_mesh_reference** out, qaws_mesh_3d const* mesh, qaws_mesh_reference_desc const* desc);
qaws_status qaws_mesh_reference_evaluate(qaws_scalar* out_positions, qaws_scalar* out_normals, qaws_scalar* out_curvatures,
	qaws_mesh_reference const* r, qaws_mesh_location const* locations, unsigned int count);
qaws_status qaws_mesh_reference_project(qaws_mesh_location* out_locations, qaws_scalar* out_distances, qaws_mesh_reference const* r,
	qaws_scalar const* points, qaws_mesh_location const* hints, unsigned int count);
qaws_status qaws_mesh_reference_walk(qaws_mesh_location* out_to, qaws_mesh_reference const* r, qaws_scalar const* tangent,
	qaws_mesh_location from);
```

A triangle mesh with vertex normals implies a smooth surface. The reference
is that surface, a patch per triangle built edge curves first so
neighbours share their edge curve:

- `PN` (the default): curved PN triangles (Vlachos et al. 2001) with a
  quadratic normal field. Creases follow their own tangent line and keep
  split normals on each side, so they stay sharp without cracks; the border
  is a crease too.
- `LINEAR`: the flat triangles, for meshes that are the truth.
- `PHONG`, `NAGATA`, `G1`, `SUBDIVISION`, `MLS`: reserved.

A location is a triangle and two barycentric coordinates. Curvatures come
from the normal field's derivative (the shape operator). The PN patch, its
derivatives and curvatures are checked against Mathematica
(`tests/reference/122_pn_triangle.wls`).

## Fields

```c
qaws_status qaws_size_field_create_uniform(qaws_size_field** out, qaws_scalar length);
qaws_status qaws_size_field_create_curvature(qaws_size_field** out, qaws_mesh_reference const* r, qaws_curvature_size_desc const* desc);
qaws_status qaws_size_field_create_feature_size(qaws_size_field** out, qaws_mesh_reference const* r, qaws_feature_size_desc const* desc);
qaws_status qaws_mesh_3d_direction_field(qaws_scalar* out_directions, int* out_face_index, qaws_mesh_3d const* mesh,
	qaws_direction_field_desc const* desc);
qaws_status qaws_volume_mesh_frame_field(qaws_scalar* out_frames, unsigned char* out_singular, qaws_volume_mesh const* mesh,
	qaws_frame_field_desc const* desc);
```

- **Size fields** give the target edge length at each point: uniform; by
  curvature from chord and angle tolerances (clamped, graded so sizes grow
  at most by `gradation` times the distance); per reference vertex; a
  callback; the smaller of two; by local feature size: the width across
  thin parts and gaps over `cells_across` (default 3), so they get that
  many edges across. Medial balls by shrinking (Ma, Bae, Choi, Kim 2012) on
  both sides of samples of the reference (thickness inward, gaps outward),
  their contacts carried onto the patches; a face across a crease never
  counts as opposite (normals within 60 degrees of perpendicular), so box
  edges are not refined. Clamped to `[min_length, max_length]` (defaults a
  thousandth and a tenth of the reference's size), graded as the curvature
  field; curvature is not included: combine with the curvature field by
  min. A fin 0.02 thick at 3 cells: 0.00667 on its faces, in 0.5 s.
- **Cross fields** on surfaces (4-symmetric, Instant Meshes'
  multiresolution smoothing along features); the index round each face in
  quarter turns marks the singularities, summing to 4 times the Euler
  characteristic.
- **Frame fields** in tetrahedral volumes: an octahedral frame per cell from
  a global solve of its spherical-harmonic coefficients (Huang et al. 2011;
  Ray, Sokolov, Levy 2016), the border's normals held, each cell projected
  onto the nearest frame, then polished over the stored facet adjacency;
  singular edges reported per cell as a mask of its six edges.

### Smoothest direction fields

```c
typedef enum qaws_direction_field_method { QAWS_DIRECTION_FIELD_LOCAL = 0, QAWS_DIRECTION_FIELD_SMOOTHEST } qaws_direction_field_method;
/* appended to qaws_direction_field_desc */
qaws_direction_field_method method;   /* LOCAL: Instant Meshes, the default */
int curvature_aligned;                /* SMOOTHEST: drawn to the principal curvature directions */
```

- **SMOOTHEST** is the globally smoothest n-direction field of Knoppel,
  Crane, Pinkall and Schroder (2013), n = 2, 4 or 6 (LOCAL stays crosses
  only). Per vertex a complex number whose argument is n times the
  direction's angle; neighbours compared through their edge's transport
  (corner angles rescaled to sum to 2 pi round a vertex, pi on the border:
  the s = 0 connection, without the paper's holonomy correction). One
  sparse factor of the connection Laplacian (cotangent weights) serves it
  all: vertices along feature lines and the border held along their lines
  (corners free), the rest from one solve; a connected part without lines
  (a closed surface) takes its own smallest eigenvector (inverse power
  iteration, at most 20 steps). Face indices are in 1/n turns and sum to
  n chi, coarse meshes included (an octahedron, an icosahedron, a cube): 8,
  12 and 4 on a sphere for n = 4, 6, 2; the held vertices lie along their
  creases to rounding.
- **curvature_aligned** draws the field to the direction of largest
  curvature on the mesh's PN reference, as strongly as the surface bends one
  way more than the other (Knoppel's lambda_t = 0, the same factor): a torus
  comes out along its parallels without a singularity.
- **The factor** (`internal/qaws_internal_ldlt.h`, for other solvers to
  share) is a sparse LDL^T: rows ordered by geometric nested dissection,
  Davis's elimination tree and column counts, then supernodal left-looking
  with dense panels. A 50k-vertex sphere's field takes about 1.3 s
  (factor 0.8 s, solves 0.03 s each). Solves work in the factor's own
  workspace: one factor is not for concurrent threads.
- **INSTANT** quads now also try the global fields (the smoothest, kept on a
  tie, and the one drawn to curvature, kept when it has fewer
  singularities); where lines hold the field it is first relaxed at the
  lattice's scale (Ginzburg-Landau by MBO steps, a second factor shared by
  both), which parts the zeros the relaxed energy lets coincide. The can of
  figure 4 (all quads, 0.14) drops from 22 irregular vertices to 10, all
  on the caps (6 on the top, a +1/-1 pair of them near its rim; 4 on the
  bottom); the blob stays at 12 with 8% fewer quads (mean edge 0.115 for
  0.11), the torus at 0.
- Limits: SMOOTHEST is the relaxed optimum, so where symmetry lets zeros
  coincide they sit close together (a disc held by its border: its four
  cones near the centre; INSTANT parts them, the public field does not).
  The smoothest field of a fat torus has pairs of singularities (+1
  outside, -1 inside); ask for `curvature_aligned` to avoid them. Polygons
  are weighted by their centroid's dual length, not a polygon Laplacian.

## Remeshing and quads

```c
qaws_status qaws_mesh_remesh(qaws_mesh_3d** out, qaws_mesh_report* report, qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc);
qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d** out, qaws_mesh_report* report, qaws_mesh_3d const* mesh,
	qaws_quadrangulate_desc const* desc);
```

- **LOCAL** remeshing (Botsch and Kobbelt 2004, as pmp, geometry-central
  and Mmgs refine them): split, collapse, flip toward each vertex's ideal
  valence, relax along the reference (walked, never projected), toward a
  uniform length or a size field; optionally every angle within bounds
  (Wang, Yan et al. 2019) and every face within a distance of the
  reference (as Mmgs's hausd). Feature corners stay, feature lines slide.
  Attributes are transferred from the reference's mesh at the footpoints,
  seams kept per side. The default.
- **CVT** remeshing: a centroidal Voronoi tessellation restricted to the
  reference (Lloyd, then L-BFGS), creases locked: the most regular and
  uniform triangles, about five times slower. `FIELD_ALIGNED` is reserved.
- **INSTANT** quads (Jakob et al. 2015, with QuadriFlow's consistency,
  Huang et al. 2018): a dense remesh; a cross field, the best of Instant
  Meshes' local one and the globally smoothest ones (free, or drawn to the
  principal curvature directions; Knoppel et al. 2013); a lattice position
  field; the lattice steps between neighbours made to agree round every
  triangle by a min-cost flow, cones on lattice points; the mesh read off
  them, half cells paired into quads, doublets removed; then finishing
  (valence pairs cancelled by local quad operations, the quads relaxed by
  their angles). Quad-dominant, or all quads with `pure` at the target
  spacing. Irregular vertices stay at the field's
  singularities: a torus is a grid along its parallels, a box's faces are
  grids. Lattice points are unfolded across creases before their offsets
  are read, so quads follow sharp edges.
- **GLOBAL** quads: an integer-grid map on INSTANT's dense mesh, read off
  as QEx does (below); all quads, the best layouts in the library, INSTANT's
  mesh returned when no valid map comes out. `PAIRING`, `SPLIT` and
  `LAYOUT` are reserved.

### LOCAL: quality by default, bounded angles, a distance guard
qaws_remesh_desc d;
memset(&d, 0, sizeof d);                       /* LOCAL, defaults */
d.target_length = 0.11f;
d.min_angle = (qaws_scalar)(35 * M_PI / 180);  /* optional: Wang et al. 2019's bounds */
d.max_angle = (qaws_scalar)(86 * M_PI / 180);
d.tolerance = 2e-3f;                           /* optional: Mmgs's hausd */
qaws_mesh_remesh(mesh, &d, &out, &report);
By default LOCAL now does what pmp, geometry-central and Mmgs do:
relaxation to the ring's area centroid weighted by the size field (a
graded mesh keeps its grading; three sweeps in the last third of the
iterations with a size field, angle bounds or a tolerance, one otherwise,
since INSTANT's dense mesh keeps fewer singularities with one), collapses onto the edge's midpoint on the
reference (or on the crease) and refused when a corner would open past 120
degrees, flips toward each sector's ideal valence (a smooth crease side or
a straight border wants 3 faces of 60 degrees) refused when they lower the
smallest angle below 30, needles collapsed, caps flipped, Delaunay flips to
finish, and a stop once nothing changes. The new fields are appended to
`qaws_remesh_desc`; zero keeps the defaults.
- `min_angle`, `max_angle` (radians, either may be 0) run the
  angle-bounded pass after the loop: rounds over the worst faces, large
  angles cut by a flip, moved away (a pattern search over the face's
  corners) or split away (the longest edge split, the flip round the new
  vertex that settles best), small ones collapsed, flipped or moved away;
  an operation is kept only when the worst angle of its star gets better
  and no angle inside the bounds leaves them. At border, crease and corner
  vertices a sector of angle theta is exempt when no whole k has theta / k
  inside the bounds (a 20-degree wedge is; a cube's 90-degree sector is
  not). `INVALID_ARGUMENT` for a negative value, `min_angle` of 60 degrees
  or more, `max_angle` of 60 or less, or `min_angle >= max_angle`.
- `tolerance` (> 0), as Mmgs's hausd: the sizes shrink where the surface
  bends (an edge's chord gap at most 3/4 of it); no collapse, flip or move
  leaves a face whose centre or edge midpoints are farther than it from the
  reference (in the loop, unless a face it replaces was farther already);
  edges whose midpoint is farther are split whatever their length. After
  the loop the guard samples the edges' quarter points too: faces farther
  have their longest edge split, and the Delaunay flips and the angle pass
  are held to it. The guard is sampled: between samples a face may bulge
  farther (dense sampling of the blob at 2e-3 finds 2.38e-3). It wins over
  the angle bounds.

Measured (test 131, figure 10): the can of test 125 at 0.12 runs 33.9 to
108.0 degrees by default (27.2 to 117.4 before); the blob at 0.11 has its
mean edge 1.5% over the target (7% under before), 82.5% of its vertices of
valence 6 (78%), all Delaunay; with [35, 86] every angle not exempt is
inside on the icosphere, the blob, the adaptive bumped sphere, the can, a
cube at 0.08, 0.15, 0.25 and 0.6 and an open square, and the can fits
[30, 90]. With a tolerance every sample (centre, edge midpoints, quarter
points) is within it: on the blob 2e-3 takes 66% more vertices, 5e-4 and
2e-4 keep every angle within 9.8..145 degrees, and with [35, 86] too the
bounds hold. The angle pass costs 1.1 times the default remesh on the
blob, but more where it fights to the last corners: about 3 s on the cube
at 0.15 (16 times its default), 1.5 s at 0.6 (80 times a 20 ms remesh). A
tolerance costs 5 to 10 times the default at the same density (the guard
on every move). Figure 3 runs about 1.05 times as long as before. Not
reached: a strip narrower than the target (its border corners stay near 26
degrees), and the faces next to a 20-degree corner.
### Metric remeshing

```c
qaws_status qaws_metric_field_create_vertices(qaws_metric_field** out_field, qaws_mesh_reference const* reference, qaws_scalar const* metrics);
qaws_status qaws_metric_field_create_hessian(qaws_metric_field** out_field, qaws_mesh_reference const* reference, qaws_scalar const* values,
	qaws_hessian_metric_desc const* desc);

qaws_remesh_desc d;
memset(&d, 0, sizeof d);                       /* LOCAL */
d.reference = reference;                       /* the metric's */
d.metric = metric;
```

A metric field is a symmetric positive definite tensor M per reference
vertex (six values: xx, xy, xz, yy, yz, zz), interpolated linearly: an edge
e is of unit length where e^T M e = 1, so R diag(1/h1^2, 1/h2^2, 1/h3^2)
R^T asks for lengths h1, h2, h3 along R's columns. For adaptation to a
solution, `qaws_metric_field_create_hessian` takes one value per vertex
and an error: the Hessian H fitted by least squares on each vertex's two
rings in its tangent plane, M = (2/9) |H| / error (Frey and Alauzet 2005),
its sizes clamped to [min_length, max_length], so the linear interpolation
of the values stays within about the error.

With `metric` set, LOCAL measures everything in the metric (each vertex's
restricted to its tangent plane): the split and collapse lengths (4/3 and
4/5 of unit length), the guards' angles and shapes, the Delaunay flips and
the relaxation's weights (metric areas). A size field or the tolerance caps
its sizes; `target_length` is unused; 30 iterations by default (it stops
when nothing changes). CVT refuses a metric (UNSUPPORTED_OPERATION).

Test 158: a constant metric of 0.02 by 0.15 on the unit square gives 723
triangles (the metric area over a unit triangle: 770), every edge within
[1/sqrt 2, sqrt 2] of unit length, mean shape 0.97 in the metric; the
Hessian metric of tanh(12 (x - 1/2)) at an error of 0.01 on a 64 x 64 grid
(recovered to 8% along x), 133 vertices, 92% of the edges of unit length
(6% shorter), the interpolation within 0.0126, the triangles of the layer
6.2 times longer along it than across; an open cylinder at 0.08 around and
0.4 along, 99.6% of unit length. About 0.3 s each. Figure 21 (key l).

### CVT remeshing
`QAWS_REMESH_CVT` places the vertices at a centroidal Voronoi tessellation
of the reference (Du, Faber, Gunzburger 1999), on its restricted Voronoi
diagram (Yan, Levy et al. 2009), by Lloyd then L-BFGS (Liu et al. 2009):
the most uniform and regular isotropic triangles. It takes `reference`,
`target_length` and `size` as LOCAL does; `iterations` counts the L-BFGS
steps after 5 Lloyd steps, per level (0: 30).

1. LOCAL runs its 10 iterations: a valid start with features and
   footpoints, and no vertex left inside a part thinner than the size.
2. The vertex count is set from the size: the integral of
   2 / (sqrt(3) h^2) over the reference, one hexagon of side h per vertex
   (the shortest free edges collapsed or the longest split, to 1%).
3. Vertices on feature lines and the border are spread along them by a
   one-dimensional Lloyd (density h^-3), then locked. Corners stay.
4. The free vertices move on the CVT energy, density h^-4. The reference is
   sampled once into flat sub-triangles (up to 16 per edge, about the
   size, fewer on slivers); each cell is cut from them by its neighbours'
   bisectors, up to the radius of security, region by region, so no cell
   crosses a crease (a locked vertex's cell starts on every side of its
   creases). A vertex walks along the reference, at most half the size
   per step, never onto the other side of a thin part, and no move turns
   a face by more than 60 degrees or against the reference's normals.
5. The connectivity stays the working mesh's, kept extrinsic-Delaunay by
   flips (never across or along a crease): the restricted Delaunay
   triangulation on a well-sampled surface, a manifold whatever the seeds
   do.
Lloyd and L-BFGS remove a lattice's dislocations slowly, so when the
surface takes 400 vertices or more, steps 2 to 4 run first on a coarse
level: a quarter of the vertices at twice the size, every other vertex of
each feature line merged. Each triangle is then split in four, which keeps
every vertex's valence and gives the lines back LOCAL's count, and the full
level is relaxed from there. Where the split would fold faces (a part
curved tighter than twice the size) the coarse level is dropped.

Measured (test 132, against LOCAL on the same input):
| Case | CVT | LOCAL |
|---|---|---|
| Cube with creases at 0.25: angles | 35.5..100.1, creases exact | |
| Icosphere at 0.12: valence 6 | 95.2% (5 to 7: 100%) | 77.6% |
| Icosphere: angles, edge CV, mean edge | 40..95, 0.058, +0.8% | 36..105, 0.117, -6.6% |
| Blob at 0.11: irregular vertices | 5.4% | 22.4% |
| Cube with creases at 0.25: angles, crease edges | 36..97, 96, exact | 22..122, 96 |
| Adaptive bump: mean edge / size - 1 | 0.045 | 0.106 |
| Slab 0.02 thick at 0.05: angles, folds | 18..97, none | 18..91, none (16..129 before its cap flips) |
| Ellipsoid 0.3 thick at 0.06: angles, faces turned over | 39..99, none | 34..102, none |
| 2.4k vertices | about 2 s (75 energy evaluations) | about 5 times faster |

Features and attributes are as exact as LOCAL's: crease vertices on the
creases, corners kept, attributes transferred at the footpoints. Known
limits: cells are found breadth first over connected sub-triangles, so a
cell does not reach across empty space (thin parts are sampled per side);
a smooth part about as thin as the size, with no crease along its rim
(an ellipsoid 0.06 thick at 0.06), keeps poor triangles and sharp folds
there, as LOCAL does; near locked corners and borders angles reach 29
degrees; a surface curved tighter than twice the size loses the coarse
level (about 80% valence 6 instead of 95%) and pays for it in time.
### Quads by min-cost flow
INSTANT's consistency follows QuadriFlow (Huang, Zhou, Niessner, Shewchuk,
Guibas 2018):
- The rounded lattice steps are made to close round every dense triangle
  by a min-cost flow. Each face has two nodes, the rows of its closure sum
  in frames combed across the faces; a unit of an edge's step is a unit of
  flow between the rows it enters. The border is a sink, and a cone takes
  units at a price. A change costs how far it takes the cheaper end of its
  edge from its own lattice point, so extra rows run where the rounding
  was nearly tied. The solver is generic (`qaws_internal_flow.h`:
  successive shortest paths, Dijkstra on reduced costs, potentials moved
  only where a search went, deterministic).
- Each cone where the field turns by a quarter is put on a lattice point
  (its closure sum made even, since det(I - R) = 2) by unit paths to
  another odd cone or the border, then onto a vertex of its face where
  moving that vertex's lattice point there is cheap. A cone left off its
  vertices leaves a hole, filled from its centre: quads round a vertex of
  valence 3 or 5, or a fan whose triangles pair with their neighbours.
- Faces the changes turn over in the lattice are undone by flip removal:
  vertices moved by a lattice step, alone, with a neighbour or with the
  rest of their lattice point, while that leaves fewer turned over.
- Positions are refitted to the mended steps, and the mesh is read off the
  steps alone: a dense triangle spanning half a cell is that half cell,
  two sharing their diagonal are its quad, so cells sheared in space come
  out whole.
- `pure` runs at the target spacing. Only when a face other than a quad
  is left does it fall back to fields at twice the target with every face
  split in four.

Measured (test 134, figure 13). The integer steps turn a change of
1e-12 in the input into another mesh, so the blob and the can are each
run six times, scaled by 1 + k 1e-12, and the counts are means over the
six (before: e280d24 on the same inputs), with 131's LOCAL, 133's field
and 135's finishing in place:

- the blob quad-dominant at 0.16: no face but quads on all six (were 6 to
  8 triangles each), 8.0 irregular (8 on all six; were 8.3);
- the blob all quads at 0.11, at the target spacing: 8.8 irregular (8 to
  11; were 13.2, split from twice the target). The field has 8 singular
  faces on all six; each cone left off its vertices is a 3-5 pair (one on
  two of the six). A cone is placed at up to two steps of two, which flip
  removal mends: at one, 10.2;
- the can all quads at 0.14: 8.0 (8 on all six; were 17.2);
- a torus at 0.3: 640 quads, none irregular; a wavy B-spline patch in
  quads: none (were 4).

Figure 13's blob (figure 4's, its coordinates rounded otherwise) gives 8
and 10 (were 8 and 12). Consistency takes about 0.5 s on 50k dense faces.
A target about the shape's size is a knife edge: a disc of radius 1 at
0.82 to 1.20 comes out empty at half the spacings in either mode.
### Quad finishing
INSTANT's quads are finished before the attributes are transferred
(`src/qaws/mesh/qaws_quad_improve.c`, internal; no public change):
- **Ideal valence** per vertex: 4 inside; on a crease, the border or a
  corner, the quads its corner angles hold (an n-gon's corner counts as
  n - 2 of n straight angles), plus one on the border. A box's corners want
  3, its edges 4, a straight border 3, a square corner 2. E sums the
  squared deviations.
- **Operations** on the corner store (`qaws_mesh_core.h`): `hm_rotate` (the
  edge between two quads turned to the next or previous diagonal of their
  hexagon; turned back, the same faces at the same indices),
  `hm_collapse_diagonal` (a quad's opposite corners merged, the quad
  removed) and `hm_remove_doublet` (a vertex of valence 2 between two quads),
  each with its legality check: border and feature edges, corner vertices,
  non-manifold results and faces sharing more than the edge are refused;
  feature flags and corner attributes are kept.
- **Cancellation** round the vertices off their ideal, the largest drop of
  E first (the candidates collected once per pass, at most four passes). An
  operation is made only when E drops with no more vertices off their
  ideal, corners keep their valence, line vertices never go below theirs,
  and its patch, relaxed on a copy first, keeps its smallest scaled Jacobian
  (at least max(0.3, before - 0.1), no lower where it was poorer) with no
  face newly turned against the reference. Where no operation alone helps,
  a rotation and a second operation beside it, judged together against the
  quads as they were before the rotation. A quad turned inside out (scaled
  Jacobian at or under 0; the extraction may leave a few) has its diagonals
  collapsed whatever E says, made only when its patch's smallest scaled
  Jacobian rises.
- **Relaxation** on the reference, vertex by vertex while their
  neighbourhoods move (at most 30 sweeps): a free vertex goes to the better
  of Zhou and Shimada's angle target and a Newton step down a shape energy
  (its corners' angles against the regular polygons', its edges' lengths
  against the target length, weighted 0.1), walked; a line vertex to the
  middle of its line neighbours along the creases. No move lowers the
  scaled Jacobian round the vertex (or adds a quad at or under 0 there),
  shortens its shortest edge below half the target length (or below what it
  was), or takes a quad's shortest to longest edge under 0.3 (or lower where
  it was already). Then the quads
  under 0.5, worst first, their free vertices searched in eight directions
  under the same guards.

Measured (test 135, figure 14), all quads: the can at 0.14 from 22 to 20
irregular vertices, SJ mean 0.966 to 0.974, angle deviation 5.8 to 4.7
degrees; the blob at 0.11 keeps 12, angle deviation 5.7 to 4.8; the wavy
patch keeps 4, angle deviation 4.6 to 3.7; test 125's torus at 0.3 from 12
to 11, SJ min -0.55 to 0.74. Quads with an edge under 0.2 of their longest
and the shortest edge stay at the plain relaxation's or better; the mean
edge-length deviation rises from 0.124 to 0.146 on the can and from 0.105
to 0.118 on the blob (the angles bought with the lengths). About 20% more
time (medians of alternating runs: can 353 to 414 ms, blob 668 to 793 ms,
on a loaded machine). The box keeps its 8 corners of valence 3, the tori
stay grids. Limits: single operations that lower E exist but each would add
vertices off their ideal (a 6 split into two 5s and the like) and is
refused on purpose, so the gains come from pairs; what is left are 3-5
pairs two edges apart and lone cones; the can's rims keep slivers with
three crease vertices in a row (SJ near 0 or below), which only a collapse
along the crease would remove.

### Guiding the mesh
What a modelling tool's user steers with (Instant Meshes' strokes and
painting). All of it is arrays appended to the descs; zeroed fields keep
today's results exactly (tests 125, 126, 133 to 135 unchanged to the byte).

```c
qaws_quadrangulate_desc d = {0};
d.target_length = 0.12f;
d.pure = 1;
d.size = painted;                    /* qaws_size_field_create_vertices on d.reference */
d.guide_directions = dirs;           /* 3 per input vertex, zero: none */
d.guide_weights = weights;           /* > 0 pull (1: as much as the neighbours), < 0 or inf held, 0 none */
d.guide_edges = loop;                /* pairs of input vertices */
d.guide_edge_count = n_loop;
d.pinned = pins;                     /* input vertices */
d.pinned_count = n_pins;
d.cone_vertices = corners;           /* input vertices made irregular */
d.cone_indices = turns;              /* quarter turns: +1 valence 3, -1 valence 5 */
d.cone_count = n_cones;              /* a closed surface: they sum to 4 chi */
d.attractors = points;               /* 3 each: free singularities drawn there */
d.attractor_strengths = strengths;   /* SMOOTHEST in cones, < 0 pushes */
d.attractor_count = n_points;
d.attractor_radius = 0.12f;          /* near the spacing */
qaws_mesh_quadrangulate(&quads, &report, mesh, &d);
```

| Input | Field | What honours it |
|---|---|---|
| Density paint | `qaws_quadrangulate_desc.size` | the dense remesh graded by a third of it; the lattice's spacing per vertex (positions, steps, refit at each pair's mean spacing); faces whose spacings differ by over 2% may keep a dislocation (the extra row a graded lattice needs) at the price of three unit changes, so they are not sheared away; with `pure` its pentagons go through the twice-and-split path |
| Directions | `guide_directions`, `guide_weights` in `qaws_direction_field_desc` and `qaws_quadrangulate_desc` | LOCAL: MF_PREFER (pull) or MF_ALONG (hold) on level 0, carried up every level; SMOOTHEST: a held value, or a penalty w d_v \|u - g\|^2 (d_v the vertex's edge weights) in the one factor; quads: carried to the dense mesh through each vertex's footpoint triangle (pulls fade across it, holds cover the dense vertices nearest a held input vertex), every candidate field honours it |
| Edge loops | `guide_edges`, `guide_edge_count` | the dense vertices within half a dense edge of a segment moved onto it, held along it and on one lattice line (MF_GUIDE), the flow changing their steps along the line only; their quad vertices placed at their lattice points on it and held by the relaxation and the finishing; not creases, never flagged as features |
| Pinned vertices | `pinned`, `pinned_count` (also `qaws_remesh_desc`) | LOCAL and CVT keep them as corners (never moved, collapsed or placed); the lattice holds a point there (MF_PINNED); the quad vertex is the input vertex exactly, flagged a corner |
| Placed cones | `cone_vertices`, `cone_indices`, `cone_count` in `qaws_direction_field_desc` (SMOOTHEST; LOCAL: UNSUPPORTED_OPERATION) and `qaws_quadrangulate_desc` | SMOOTHEST: the smoothest field for a connection adjusted to the cones (Crane, Desbrun, Schroder 2010; Knoppel et al. 2013, 6.1): the least change of every edge's transport (min-norm, one sparse LDL^T of the faces' D D^T) that turns each face by its angle share of the cones round it instead of its own excess, so the field winds at the cones and nowhere else; edges between two held vertices keep theirs. Poincare-Hopf: a closed component without held vertices needs the indices to sum to n chi (8 quarter turns on a sphere), else INVALID_ARGUMENT; with lines or a border the remainder is left free (spread evenly over the faces that cannot close it, the solve places it). Quads: the dense remesh keeps each cone vertex, only the SMOOTHEST candidate with the cones is used, the lattice points round the cone's star are set on its vertex; with a border the flow keeps the star's steps (MF_CONE) and the vertex stays pinned, so the cone's quad vertex is the input vertex; on a closed surface the cones must absorb the lattice's mismatch (at a hundred times a free cone's price), so they may move a step and the nearest quad vertex of valence 4 - index within the spacing is moved onto the cone; either way it is held by the relaxation and the finishing like a corner, not flagged one |
| Attractors | `attractors`, `attractor_strengths`, `attractor_count`, `attractor_radius` (0: a tenth of the bounding diagonal) in both descs | a Gaussian exp(-d^2 / r^2) round each point. SMOOTHEST: the strength in cones (1: a +1 cone's turn) spread by area times the Gaussian, added to the adjusted connection's turns and taken back from every face by its area: a free cone forms there; negative pushes them off. LOCAL: the smoothing's edge weights times exp(-strength Gaussian) on every level. Quads: both candidates |

Measured (test 137, figure 16, a latitude-longitude unit sphere of 1106
vertices):
- density 0.08 above the equator, 0.16 below, all quads: mean edges 0.081
  and 0.137 away from the transition (without the dislocations: 0.088 and
  0.110), every quad's scaled Jacobian above 0.68, 33 irregular vertices;
- a stroke round the equator at 30 degrees to it: the field on the stroke
  within 0.01 degrees (LOCAL, pulled or held) and 0.07 (SMOOTHEST, the mesh
  normals against the sphere's); the quads' edges there 8.5 degrees off its
  crosses on average, 17.0 without it (14 irregular against 8);
- a parallel as guide edges at 0.12: every segment within 0.003 of a quad
  edge (0.06 without), all quads, 8 irregular;
- pins: the quad, LOCAL and CVT vertices at the input positions to the bit,
  flagged corners.
- cones placed: eight +1 at a cube's corner directions on the sphere: the
  SMOOTHEST field's 8 singular faces each round its cone vertex; the quads
  at 0.15 (pure) a cube-sphere of 592 quads with exactly 8 irregular
  vertices, of valence 3, at the input vertices to the bit, not corners
  (SJ above 0.56 INSTANT, 0.69 GLOBAL by the map); seven cones refused
  (INVALID_ARGUMENT); a +1 / -1 pair on a 24 x 24 square: 2 singular faces
  (the corners take the border's turn), quads at 0.12 with exactly two
  irregular inside, valence 3 and 5 on the vertices (INSTANT, GLOBAL);
- attractors (reach 0.3 on the field): SMOOTHEST's nearest cone to a point
  0.435 free, 0.074 drawn (strength 1), a free cone pushed 1.68 off
  (strength -1), 8 singular faces each time; LOCAL 0.137 free, 0.074 drawn
  (3); quads at 0.15 with strength 1, reach 0.1: the nearest irregular
  vertex 0.042 from the point (0.415 without), still 8 irregular.

Limits: a pull and a hold are honoured by the field, not by the lattice: a
lattice turned against the sphere's natural one ends rows near the stroke,
so the quads there are a few degrees off and add 3-5 pairs. Density
transitions are dislocations (pentagons, or 3-5 pairs after the split), not
smooth gradings; the finishing's length energy and guards use the smallest
painted size. Guide edges are matched to dense vertices by a brute-force
search (every vertex against every segment); crossing loops get no pinned
crossing; a loop through a cone of the field cannot stay one lattice line.
Directions need the reference made from the input.
Placed cones honour the field on the cones' vertices (in one face round
each); on a closed surface the lattice cannot hold every cone at a fixed
spacing, so a cone may be a step off its vertex before its quad vertex is
moved onto it (the quads there stretched by up to a step); with a border
the border takes the mismatch. On higher genus the transport's turns round
the handles are the least change's, not chosen (Crane et al.'s
generators). Index sums are checked only on closed components without held
vertices. An attractor reaching far beyond the lattice's spacing bends the
field away from the surface's parallel transport over that reach, which the
lattice pays for in 3-5 pairs (reach 0.3 at spacing 0.15 on the sphere:
16 irregular instead of 8): keep the reach near the spacing.
### Quads by integer-grid maps (GLOBAL)
`QAWS_QUADS_GLOBAL` (`src/qaws/mesh/qaws_quad_global.c`) follows
Mixed-Integer Quadrangulation and integer-grid maps (Bommes, Zimmer,
Kobbelt 2009; Bommes et al. 2013) with QEx's extraction (Ebke, Bommes,
Chen, Kobbelt 2013), on INSTANT's dense mesh and fields:
- **Integers, first try**: the lattice steps INSTANT's min-cost flow made
  consistent. Each dense vertex keeps two unknowns, its place in its
  lattice cell in its own frame; every edge's transition is a quarter turn
  and a whole shift, so no cut graph is needed. A singular face whose cone
  sits on a vertex's lattice point makes that vertex the cone (held on the
  point), else it is split in three round its cone at the face's centroid;
  feature corners are cones on their own lattice points, charted from
  their neighbours, the corner's place carried round its triangles from
  the fixed point of the turn round it (so neighbouring triangles agree on
  it); edges between line vertices on different creases by a corner have
  their steps read again in their faces' plane.
- **Integers, when those give no map**: quantized on the map. The cones
  (singular faces, corners going all the way round) are the nodes of a
  spanning tree over the faces whose paths cross no crease or border; the
  unknowns are whole steps along a node's path to its parent (two per
  node, steps that keep both ends on lattice points) and, for a singular
  face, along its path to the nearest border edge of each direction (two
  steps along the border). Every face between stays closed. The map's
  energy over those unknowns is a small dense quadratic: a step changes
  the triangles' translations along its path only, and the places respond
  by one back substitution each on the base map's factor. It is rounded
  greedily as MIQ does (the closest to whole first with those within a
  tenth, the rest solved again), then improved by whole steps from the
  rounding or from the flow's integers, whichever is lower. Each cone then
  goes onto the free vertex nearest it in the map with the cones split
  (that vertex's lattice point shifted onto it), as MIQ keeps its cones on
  vertices.
- **Places**: one sparse LDL^T solve of sum A (|grad u - X/h|^2 +
  |grad v - Y/h|^2) over the triangles (the field turned into each
  triangle's plane), line vertices held on their whole isoline, corners
  and cones on their lattice points. Folds left are undone vertex by
  vertex (each moved to the point of its kernel farthest from its
  triangles' far sides, a ring wider every 8 rounds, at most 64). A free
  vertex wound twice round (a cone beside it a full turn short) is no map.
- **Extraction**: the lattice points of each triangle are the vertices,
  registered in every triangle round them; from each, its unit isolines are
  traced across the triangles through the transitions, each nudged a
  millionth to its left so it never runs along an edge or through a
  vertex (a border or crease isoline is traced from the side it bounds);
  a face's next edge is a quarter turn left at each point. All faces must
  close in fours, else INSTANT's mesh is returned. Vertices are placed at
  their barycentric points and projected on the reference (creases for
  line points, corners exactly), then INSTANT's relaxation and finishing.
- **Report**: `iterations` 1 when the map made the mesh, 0 when INSTANT's
  was returned.
Measured (test 136, figure 15; GLOBAL against INSTANT, all quads): torus
at 0.3 640 quads both, 0 irregular, SJ min 0.994 (0.985), AD 3.13 (3.60);
creased box at 0.1 376 quads, AD 0.00 (0.04); blob at 0.11 1397 quads
(INSTANT 1276, from twice the target split), 8 irregular both, SJ min
0.815 (0.772), AD 2.74 (3.27); can at 0.14 766 quads both, 8 irregular,
SJ min 0.853 (0.801), AD 2.23 (2.23), ED 0.129 (0.138); wavy patch at 0.1
498 quads, 2 irregular, AD 2.97 (3.22). About INSTANT's time plus 10%;
the quantization, when it runs, adds about two more solves.
Over the five shapes at 0.75, 1.25, 1.6 and 2.2 times these targets the
map makes all 20 (was 14). Those the flow's integers alone missed: can
at 1.25 512 quads, SJ min 0.860 (INSTANT 0.178), AD 2.41 (2.78); can at
1.6 SJ min 0.854 (0.223), AD 2.71 (2.81); wavy patch at 0.75 2 irregular
(3), AD 2.41 (2.39); at 1.25 2 irregular (3), SJ min 0.793 (0.751), AD
3.59 (4.08); box at 2.2 AD 0.00 (0.03); box at 0.75 SJ min 0.029 (0.000),
2 irregular (2), AD 2.50 (1.45).

## Measures

```c
qaws_status qaws_mesh_3d_measure_shape(qaws_shape_quality* out, unsigned int* out_ideal, qaws_mesh_3d const* mesh);
qaws_status qaws_mesh_3d_measure_deviation(qaws_mesh_deviation* out, qaws_mesh_reference const* reference,
	qaws_mesh_3d const* mesh, qaws_deviation_desc const* desc);
```

`qaws_mesh_measure.h` reports a mesh the way remeshing papers do (Hu et
al. 2017 and Wang et al. 2019 for triangles, QuadWild for quads).
`qaws_mesh_quality` and the operators' reports stay as they were.

- **Shape**, in one pass over the faces:
  - triangles: Q = 2 sqrt(3) r / h_max (1 when equilateral), min and mean;
    the smallest, mean smallest and largest angle; the shares of corners
    under 30 and over 90 degrees;
  - quads: each quad's smallest corner scaled Jacobian, along the unit
    Newell normal (min and mean over quads), and AD, the mean |angle - 90|
    over quad corners;
  - ED, the edges' standard deviation over their mean length.
- **Ideal valence**, from the corners rather than one number for all:
  - a smooth interior vertex (off the border and the features) wants 6, or
    4 when quads are the majority; V567 is the share of those at 5 to 7;
  - a border, crease or corner vertex wants as many faces as its corners'
    angles hold ideal corners (60 degrees per triangle, 90 per quad), plus
    one edge on the border. It is regular within a half face (a tie at a
    half rounds to its own count, one face at the least); `irregular`
    counts the vertices off that rounded ideal.
  - So a box's corners (3 quads) and edges (4), and a quad grid's straight
    border (3 edges) and corners (2), are regular. A creased quad cube
    counts 0 irregular vertices, where `qaws_mesh_quality` counts 8.
    `out_ideal` gives each vertex's rounded ideal valence.
- **Deviation, certified both ways**: lower and upper bounds of the sup of
  the distance from the mesh to the reference, and from the reference to
  the mesh, their Hausdorff interval, the area-weighted RMS (projected, not
  certified), and the point where each lower bound is reached. One branch
  and bound serves both directions; the mesh enters as its LINEAR
  reference, so polygons count as their corner fans.
  - A triangle of the side measured is matched affinely to the other
    side. Each corner goes to its closest point there. A triangle whose
    corners land on different triangles is cut along a flat unfolding of
    the other side round a shared vertex.
  - Over each piece, the difference of the two cubic patches is a cubic
    Bezier triangle. Its control points (blossoms) bound it by the convex
    hull property, and they are refined by subdivision while above the
    goal. The distance's Lipschitz bound from the corners is taken too,
    over their flat triangle and, on a curved side, widened by the farthest
    its sub-patch's net strays from that triangle's: it still shrinks where
    the closest place on the other side jumps across a concave edge, which
    the matched patches do not. A zero-area triangle's wedge (a polygon
    with a straight corner, fanned) maps nothing and is left to its
    neighbours.
  - The largest bound is split into four until it is within `precision`
    (a share of the reference's box diagonal, default 1e-4) of the largest
    distance found. A face or patch splits `budget` times (default 256) past
    the splits that bring it down to the other side's triangle size.
    `converged` is 0 when a budget ran out first; the bounds are then wider
    but still hold.
- **What is certified**:
  - the upper bounds hold unconditionally: any point of a patch is a point
    of the surface, so a poor match only widens them;
  - the reference-to-mesh lower bound is an exact distance to the flat
    triangles (Ericson's closest point, the grid's rings searched outward);
  - the mesh-to-reference lower bound is a distance found by projection:
    exact on LINEAR references, and on PN up to the patch projection
    (the reference's Gauss-Newton, then a descent that runs along a
    patch's side where the gradient leads out, and on into the neighbour
    across it);
  - both are widened by a few units in the last place, and stored outward
    in float builds.

On figure 3's 3752-triangle blob both directions converge in 0.2 to 0.5 s;
a 10.7k-face remesh takes 0.5 to 0.8 s (Release). Showcase figure 9 tabulates
figures 3 and 4's meshes: the baseline the remeshing items quote.
## Meshes of exact geometry

```c
qaws_status qaws_mesh_generate_surface(qaws_mesh_3d** out, qaws_mesh_report* report, qaws_surface const* surface,
	qaws_generate_desc const* desc);
qaws_status qaws_mesh_generate_shell(qaws_mesh_3d** out, qaws_mesh_report* report, qaws_shell_atlas const* shell,
	qaws_generate_desc const* desc);
```

`REMESH` meshes a surface (a welded parameter grid, so seams and poles
close) or a shell (its faces triangulated conformally) densely, makes the
faces' shared edges and the free sides features, then remeshes
(`QAWS_ELEMENTS_TRIANGLES`) or quadrangulates (`QUADS`, `QUAD_DOMINANT`) its
PN reference at the target size, uniform or by chord and angle tolerances.
Every vertex is then snapped onto its face's exact surface by Gauss-Newton:
the output lies on the NURBS (or any `qaws_surface`) to rounding, its
footpoints are the face and (u, v), and a face attribute `face` names each
element's face. Border and crease vertices are held on their domain's side,
corners on its corners. `PARAMETRIC`, `RESTRICTED_DELAUNAY` and
`STRUCTURED` are reserved.

### Curved elements

```c
qaws_status qaws_mesh_generate_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_mesh_3d const* mesh,
	unsigned int order);
qaws_status qaws_mesh_generate_nodes_surface(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_surface const* surface, unsigned int order);
qaws_status qaws_mesh_generate_nodes_shell(qaws_scalar* out_positions, qaws_mesh_location* out_locations, unsigned int* out_element_nodes,
	qaws_mesh_3d const* mesh, qaws_shell_atlas const* shell, unsigned int order);
```

A generated mesh's elements at order 2 or 3 (P2, P3 triangles; Q2, Q3
quads) for a solver: the nodes on the exact geometry. The node list is the
mesh's vertices, then order - 1 per edge (shared by the faces either
side), then the inner nodes: a cubic triangle's centroid, a quad's
(order - 1)^2 (none in a polygon). Each element lists its corners, its
edges' nodes from corner k toward corner k + 1, then its inner nodes
(i fastest, along its sides 0-1 and 0-3). A new node is its element's
straight point snapped onto the element's face from its corners' (u, v);
nodes of feature and border edges are held on the domain's side as the
vertices are. The torus of test 157 at 0.4: the elements' distance to the
surface between their nodes 4.2e-2 (P1), 8.8e-4 (P2), 7.4e-4 (P3);
5.6e-2, 4.0e-4, 7.5e-4 in quads. Interpolation's distance to the surface is
O(h^4) at order 2 and 3 alike (an even order gains one by symmetry), so
order 3 pays off in the nodes' places along the surface, not in that
distance. Figure 20 (key k).

## Repair for farq

```c
qaws_status qaws_mesh_repair(qaws_mesh_3d** out_mesh, qaws_mesh_repair_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_mesh_repair_desc const* desc);
```

`qaws_mesh_repair.h` turns messy input (triangle and polygon soups, scans,
meshes exported from games and modelling tools) into a closed,
consistently oriented triangle surface: the input farq's exact Booleans
take. The split (decision 9 of [meshing_plan.md](../meshing_plan.md)):
farq is exact and tolerance-free and accepts anything with a defined
inside (every directed edge cancels: self-intersecting, non-manifold and
nested input included); qaws owns every tolerance and quality judgement.
So the repair closes and orients, and leaves to farq what needs no
tolerance.

The policies, in order. Tolerances are world lengths; 0 is the default,
relative to D, the bounding diagonal of the vertices the faces use; a
negative value turns the policy off.

| Policy | Field (default) | What |
|---|---|---|
| validation | | an index out of range, a face under 3 vertices, a position not finite: `INVALID_ARGUMENT`, `report.invalid_face` / `invalid_vertex` naming the first; a face with a vertex twice dropped; polygons ear-clipped in their plane (the output is triangles) |
| welding | `weld_tolerance` (1e-6 D; off: identical points only) | points that close merged transitively (union-find over a grid's pairs), each group at its mean; faces collapsed by it dropped |
| degenerate faces | `degenerate_tolerance` (the weld tolerance) | faces lower than it (twice the area over the longest edge): a needle's short edge collapsed; a cap (a corner on its long edge) removed and the faces across the long edge split at that corner, so no hole opens |
| duplicates | | faces on the same three vertices kept once, the majority's way; a pair turned opposite ways (a membrane that cancels) removed both |
| gaps | `gap_tolerance` (1e-3 D) | boundary vertices that close to another boundary vertex merged at their mean, or to a boundary edge: its face split at the vertex (the T-junction resolved); never two vertices already joined by an edge; repeated until the borders meet |
| junk | `min_component_area`, `min_component_faces` (0: none) | parts (faces sharing a vertex) below either removed |
| orientation | | each part made consistent by turning faces across its manifold edges; a part that cannot be (a Moebius strip) counted in `non_orientable` and left, its border open |
| holes | `hole_perimeter` (all) | boundary loops up to that perimeter filled: the triangulation largest in its smallest dihedral angle, then least in area, by dynamic programming (Liepa 2003) up to 500 vertices, a fan from the centroid beyond; refined toward the surrounding edge lengths (centroid splits, Delaunay flips) and faired (the scale-dependent bi-Laplacian, the border and its ring held, one LDL^T), unless `flags` has `QAWS_MESH_REPAIR_FLAT_FILLS` (refined, not faired) or `QAWS_MESH_REPAIR_COARSE_FILLS` (the minimal triangulation only) |
| outward | `nested_shells` (`KEEP`) | each part turned to a positive signed volume, except a part inside another closed part (the generalized winding number at its vertices): at an odd depth a cavity, turned inward, or removed with `QAWS_NESTED_SHELLS_REMOVE` |

Left for farq, and reported: non-manifold edges, self-intersections (not
looked for: farq's self-union resolves them exactly), the borders of holes
not filled. `report.closed` says every directed edge cancels; the report
also counts what each policy did (welded, needles, caps, duplicates,
stitched, holes filled and left, fill faces, faces turned back, nested
shells, boundary, non-manifold and open edges).

The output carries a face attribute `source_face`: the given face each
triangle comes from, -1 for fills (exact up to 2^24 in the f32 build). Its
creation splits non-manifold edges and vertices into copies at the same
positions (as every `qaws_mesh_3d`), so the indices from
`qaws_mesh_3d_get_faces` stop cancelling there. farq should read the
positions: its exact point keys merge the copies and the directed edges
cancel again; or map the indices through `qaws_mesh_3d_get_vertex_source`
first.

Test 139: a cube soup (36 points) to the closed cube of 8; a sphere soup
off by 1e-7 back to its 642 vertices; three holes in a sphere filled to
1.4% of the area removed (flat fills 8.6% under), fill edges within 12% of
the surrounding ones, the smallest fill angle 25.9 degrees; a box in two
halves at different spacings 3e-4 apart stitched closed without fills; half
of a sphere's faces turned and turned back; a cavity turned inward, an
island in it outward; a cap, needles and a zero-area face removed without
holes; a Moebius strip reported; two cubes on one edge closed with the edge
left non-manifold. Figure 18 (key i): a broken torus soup (holes, a gap,
three faces in ten turned, duplicates) repaired in 0.02 s.

### Repair, completed

- **Holes along creases.** Before a hole is filled, every crease edge that runs into it is followed (normals parting by more than 40°):
  - two crease ends facing each other across the hole (each within 37° of the chord) are joined by a cubic Hermite path;
  - three or more crease lines that meet at a point (a corner) are carried straight to it.

  The hole is split along these paths and each part is filled, refined and faired on its own. Path vertices are never flipped or faired, so a box's lost edge or corner comes back sharp.
- **Long holes.** Borders past 500 vertices are halved by their shortest diagonal, between vertices a third to a half of the loop apart, recursively, instead of fanned.
- **Floating vertices.** Given vertices that no face uses are dropped and counted (`unused_vertices`).
- **Small parts by volume.** `min_component_volume` removes closed parts holding less volume, after orientation.
- **Non-manifold edges and vertices.** With `QAWS_MESH_REPAIR_SPLIT_NON_MANIFOLD`, each sheet gets its own vertices:
  - the faces round an edge are sorted by angle and paired across the wedge inside the solid;
  - each vertex gets one copy per fan of faces (`split_vertices`).

  Two cubes on one edge, or on one corner, come out as two closed cubes.
- **Self-intersections.** With `QAWS_MESH_REPAIR_FIND_INTERSECTIONS`:
  - every pair of faces crossing beyond their shared corners and edge is found, using exact `orient3d` and `orient2d` (coplanar overlaps and folds included);
  - the pairs are counted (`self_intersections`, `intersecting_faces`) and marked in a face attribute `intersecting`.

  Resolving them belongs to the numerical Booleans (decision 10) or to farq.

Test 139 checks:
- two cubes on an edge are split into 16 vertices, closed, with a volume of 2; on a corner, one vertex is split;
- 2 floating vertices are reported, and a 0.01 cube is removed by volume;
- two overlapping cubes give 18 crossing pairs; a clean 5,120-face sphere gives 0; two faces folded flat give 1;
- a box with a hole across an edge, or across a corner, is filled with every vertex on the box (within 1e-14);
- a 600-vertex hole is filled closed.

Figure 24 (key o): a box with a hole round its corner, filled along its creases, and two crossing cubes with their intersecting faces.

## Watertight surfaces

`qaws_mesh_make_watertight` (`qaws_mesh_watertight.h`) makes a new closed, manifold, outward surface from any triangle input:
- soups;
- scans with holes;
- overlapping or self-intersecting parts;
- open sheets.

It follows the input rather than keeping its faces. To keep the input's faces, use `qaws_mesh_repair`.

```c
qaws_watertight_desc d = { 64, 0 };      /* resolution: cells along the longest side; offset */
qaws_mesh_make_watertight(&mesh, &report, positions, faces, NULL, vertex_count, face_count, &d);
```

**Holes first.** Holes are filled by `qaws_mesh_repair`'s fills before anything else: Liepa's triangulation, refined, then faired by the bi-Laplacian with the border's ring held. The caps are tangent to the surface round them and continue its curvature. On the test sphere they come within 0.017 of it, against 0.063 for the winding number's own ½-level caps, which are close to a minimal surface and leave a ridge at the rim. `QAWS_WATERTIGHT_WINDING_CAPS` keeps those caps. With an offset, nothing is filled first, since a sheet's border is not a hole. Input without holes is used as given.

**Inside.** A point is inside where the generalized winding number (through `qaws_mesh_bvh`) is at least ½: holes close smoothly, overlaps are united and inner faces drop out. With an offset, points within it of the input count as inside too, so an open sheet or a thin part becomes a closed shell (an envelope).

**Which grid points are evaluated:**
- points near the input, each one;
- every other point elsewhere;
- a 2×2×2 block of cells with no near point, whose corners agree by more than 0.1, takes their side;
- any other block has all its points evaluated, because the inside can turn far from the input, across a hole.

On the test sphere, about a third of the grid is evaluated.

**Surface.** Dual contouring (Ju et al. 2002) puts one vertex per crossed cell, at the least squares point of the crossings' tangent planes. Rank-truncated, it stays within the cell (or falls back to the mass point), so sharp edges and corners are kept.

Each crossing comes from one of two sources:
- a ray along the grid edge onto the input, which gives the input triangle's own point and normal;
- across a hole or on an offset, bisection, with the normal from the winding number's gradient or from the closest point.

The quads are split along their shorter diagonal. The result then goes through `qaws_mesh_repair`, which splits non-manifold places and drops slivers.

Detail below the cell size is lost, and the face count follows the grid; remesh the result to the size you want.

**Tests and figure.** Test 161:
- a sphere soup with three holes: closed, within 0.006 of the sphere away from the holes, volume 4.12;
- two crossing cubes: their exact union, every vertex on its surface (5e-17), the volume to 1e-6, no face crossing another;
- an open sheet wrapped at 0.05: closed, within 0.0033 of the offset;
- refusals.

Figure 25 (key p) shows the holed sphere and its watertight surface, the crossing cubes' union, and a wavy sheet wrapped at 0.06.

## Self-intersections resolved

`qaws_mesh_resolve_intersections` (`qaws_mesh_intersect.h`) cuts every pair of triangles that meet beyond their shared corners and edges along their intersection. The result's triangles meet only at shared vertices and edges: the curves where the input crossed itself become edges held by the faces of both sides (a conforming arrangement). Nothing is removed; deciding what is inside is left to Booleans and winding numbers.

```c
qaws_intersect_desc d = { 0 };          /* tolerance; 0: 1e-9 of the box's diagonal */
qaws_mesh_resolve_intersections(&mesh, &report, positions, faces, NULL, vertex_count, face_count, &d);
```

It is numerical, not exact (farq makes exact arrangements), and robust by construction:
- **Snapping.** Distances under the tolerance count as zero, so coplanar and nearly coplanar pairs are handled together (each overlap is cut along the other's edges), and a new point is snapped to the input feature near it: a vertex before an edge before a face.
- **Keyed points.** Every new point is keyed by what made it: an edge through a triangle's plane, three triangles' planes (where three cuts meet), or two edges of coplanar triangles. It is computed once, so every face holding it takes the same vertex.
- **Carriers.** Each input edge and each cut (the line where two faces cross) keeps its points sorted along it. Faces settle until no carrier gains a point: points lying on a piece are added to it, and pieces crossing inside a face are cut at their keyed point. Both sides of a cut therefore split it identically.
- **Triangulation.** Each cut face is triangulated again in its plane, constrained Delaunay through its cuts and the points on its edges. Points on an edge sit exactly on it, and orientations are exact.

The result is checked by exact predicates (`qaws_mesh_repair`'s test). `report.remaining` counts pairs still crossing, and `report.unresolved` counts cuts the triangulation could not keep. The face attribute `"source"` gives each triangle's input face. Where an intersection curve holds four faces, the mesh gives its vertices copies (`qaws_mesh_3d_get_vertex_source` gives the shared vertex).

**Tests and figure.** Test 162 covers five inputs: two cubes crossing at an angle, three cubes crossing (cuts meeting at triple points), two cubes sharing two planes, the same pair 1e-12 off and tilted, and two spheres crossing. In every case:
- no faces cross;
- nothing is unresolved;
- the directed edges cancel (vertices matched by position);
- the area is unchanged to 1e-12.

A cube alone comes back unchanged, and invalid input is refused. Figure 26 (key q) shows two spheres crossing, their faces outside the other sphere kept (the union) and inside (the lens), and three cubes resolved.

## Spatial queries: qaws_mesh_bvh

`qaws_mesh_bvh.h` answers point and ray queries on any mesh, broken or not.
A BVH holds simplices of one kind, in 2D or 3D:

- segments: curves, or a 2D mesh's border;
- triangles: surfaces, or a 2D mesh's faces;
- tetrahedra: volumes.

Each simplex belongs to an input element, which is what a hit reports:
- a polygon is fanned from its first vertex;
- a cell is split into tetrahedra from its first vertex (a hexahedron 6, a prism 3, a pyramid 2).

```c
qaws_mesh_bvh* bvh = NULL;
qaws_mesh_bvh_create_mesh_3d(&bvh, mesh, NULL);            /* or _create (arrays, soups), _create_mesh_2d, _create_volume */
qaws_mesh_bvh_closest(hits, bvh, points, n, 0, NULL);      /* nearest point; inside a 2D face or a cell: located, distance 0 */
qaws_mesh_bvh_ray(hits, bvh, origins, dirs, n, 0, 0, NULL);/* first hit; any_hit 1 for shadow rays */
qaws_mesh_bvh_winding_number(w, bvh, points, n, NULL);     /* 1 inside, 0 outside, between for open input */
qaws_mesh_bvh_signed_distance(d, NULL, bvh, points, n, NULL);
qaws_mesh_bvh_refit(bvh, moved_positions);                 /* same simplices, new positions */
```

A `qaws_mesh_hit` holds:
- the element and the simplex;
- the simplex's vertices with barycentric weights, for interpolating attributes;
- the point;
- the distance (for rays, the parameter t).

**Build.** Binned SAH with 8 bins per axis, as tinyBVH's default, and leaves of up to 4 simplices. Positions are held in double whatever the caller's precision. Building 199,200 triangles takes 0.15 s.

**Rays.** Triangles in 3D use the watertight test of Woop, Benthin and Wald (2013): a ray through a shared edge or vertex never slips between triangles. Segments in 2D decide each vertex's side by the same rule in both segments that share it.

**Winding numbers.** These are generalized winding numbers (Jacobson et al. 2013): solid angles over 4π in 3D, angles over 2π in 2D.
- Far nodes use their expansion to second order (Barill et al. 2018). The expansion uses each node's area vector, first moment and second moment, the second moment including each triangle's own spread.
- `winding_accuracy` (β, default 2) decides when a node is far: when the point lies beyond β times the node's radius.
- Worst error against the exact sum on a closed sphere: 2.2e-3 at the default, 1.3e-5 at β = 6.

**Signed distance.** The closest distance, negative where the winding number is at least 1/2. Its sign therefore holds on soups, with holes and with overlaps.

**Threads.** Queries run in chunks of 64 through the caller's `qaws_batch_executor`, and the results do not depend on it.

**Status.** `UNSUPPORTED_OPERATION` for rays or winding numbers except on triangles in 3D or segments in 2D. `INVALID_ARGUMENT` for:
- a dimension other than 2 or 3;
- tetrahedra in 2D;
- an index out of range;
- a position that is not finite.

**Tests and figure.** Test 159:
- closest points and rays equal a brute-force loop;
- winding numbers within 5e-3 of the exact sum;
- 1089 rays through a warped grid's vertices and edges, none missed;
- 2D location and border queries;
- points located in a block of hexahedra;
- refit equal to a new build;
- an executor running chunks backwards gives the same results;
- 100k closest points within 0.1 of a 200k-triangle sphere take 0.5 s.

Figure 22 (key m) is a torus soup with three holes:
- ray traced, with camera and shadow rays;
- its winding number on a plane;
- its signed distance on the same plane.

## Quality reports

`qaws_quality_report.h` measures every element of a 2D, surface or volume mesh, then summarizes those values in the forms a solver checks.

```c
qaws_scalar* v = malloc(sizeof(qaws_scalar) * face_count);
qaws_mesh_3d_measure_elements(v, mesh, QAWS_MEASURE_MIN_ANGLE);   /* or _mesh_2d_, _volume_mesh_ */
qaws_quality_summary_desc d = { QAWS_MEASURE_MIN_ANGLE, 1, 30, 24, { 0, 60 }, 10 };
unsigned int hist[24], worst[10];
qaws_quality_summary s;
qaws_quality_summarize(&s, hist, worst, v, face_count, &d);       /* s.failing: corners under 30 degrees */
```

Measures (1 is ideal where the measure is a ratio):

| Measure | Faces | Cells |
|---|---|---|
| `MIN_ANGLE`, `MAX_ANGLE` | corner angles, in degrees | dihedral angles between facets |
| `SHAPE` | triangles: mean ratio 4√3 A / Σl² | tetrahedra: 12 (3V)^(2/3) / Σl² |
| `ASPECT_RATIO` | triangles: R / 2r | tetrahedra: R / 3r |
| `SCALED_JACOBIAN` | smallest corner (triangles scaled by sin 60°, quads by the square) | smallest corner (scaled by the regular tetrahedron, cube, prism, and equal-edged pyramid; a pyramid's base corners) |
| `EDGE_RATIO` | longest over shortest | the same |
| `SIZE` | area | volume |

NaN marks a measure that does not apply to an element: SHAPE and ASPECT_RATIO on quads, polygons and non-tetrahedral cells, and SCALED_JACOBIAN on faces of more than four sides.

**Inverted elements.** A 2D mesh is measured in its plane, and cells by their orientation. So turned elements show as a negative SHAPE, SCALED_JACOBIAN or SIZE, and their angles pass 180°. A surface face's corners are measured along its Newell normal.

**Summary.** `qaws_quality_summarize` takes any array of values and returns:
- min, max, mean and median;
- a histogram over a given span, or over the values' own span when none is given;
- the worst elements in order;
- the count past a threshold.

The measure in the summary's description says which way is worse: higher for MAX_ANGLE, ASPECT_RATIO and EDGE_RATIO, lower for the others. NaN values are skipped and counted.

**Tests and figure.** Test 160 checks:
- the closed forms for regular and right triangles, the square, the regular tetrahedron, the cube, the prism and the pyramid, and inverted 2D faces;
- that the smallest angle and Jacobian agree with `qaws_mesh_3d_measure_shape` on a wavy quad grid;
- that the smallest scaled Jacobian agrees with `qaws_volume_mesh_measure_quality` on a sheared block;
- the summary's statistics, histogram, worst order and thresholds.

Figure 23 (key n) shows a jittered 2D mesh by shape, the blob by smallest angle, and its tetrahedra by smallest dihedral, each with its histogram.

## Volume meshes: qaws_volume_mesh

```c
typedef struct qaws_volume_mesh_desc {
	qaws_scalar const* positions; unsigned int vertex_count;
	unsigned int const* cell_vertices; qaws_cell_type const* cell_types; unsigned int cell_count;   /* types NULL: tetrahedra */
} qaws_volume_mesh_desc;

qaws_status qaws_volume_mesh_create(qaws_volume_mesh** out_mesh, qaws_volume_mesh_desc const* desc);
qaws_status qaws_volume_mesh_get_cells(unsigned int* out_cell_vertices, qaws_cell_type* out_types, qaws_volume_mesh const* mesh);
qaws_status qaws_volume_mesh_get_cell_adjacency(unsigned int* out_cells, qaws_volume_mesh const* mesh);
qaws_status qaws_volume_mesh_get_boundary(qaws_mesh_3d** out_surface, qaws_volume_mesh const* mesh);
qaws_status qaws_volume_mesh_measure_quality(qaws_volume_quality* out_quality, qaws_volume_mesh const* mesh);
```

Cells follow geogram's pattern: runs of corners, and per cell facet the
facet across it (stored; none on the border). Local numbering is VTK's,
with positive orientation:

| Cell | Vertices | Facets (outward) |
|---|---|---|
| tetrahedron | `0 1 2 3`, det(v1 − v0, v2 − v0, v3 − v0) > 0 | `1 2 3`, `0 3 2`, `0 1 3`, `0 2 1` (facet i opposite vertex i) |
| hexahedron | `0 1 2 3` the base (normal toward the top), `4 5 6 7` above | `0 3 2 1`, `4 5 6 7`, `0 1 5 4`, `1 2 6 5`, `2 3 7 6`, `3 0 4 7` |
| prism | `0 1 2` the base, `3 4 5` above | `0 2 1`, `3 4 5`, `0 1 4 3`, `1 2 5 4`, `2 0 3 5` |
| pyramid | `0 1 2 3` the base, `4` the apex | `0 3 2 1`, `0 1 4`, `1 2 4`, `2 3 4`, `3 0 4` |

The border comes back as a `qaws_mesh_3d` turned outward, with a vertex
attribute `vertex` (the volume's vertex) and a face attribute `cell`.
Quality: dihedral angles and circumradius over shortest edge for
tetrahedra, the scaled Jacobian for hexahedra, edges and volume.

## Tetrahedra

```c
qaws_status qaws_volume_mesh_tetrahedralize(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report,
	qaws_mesh_3d const* boundary, qaws_tetrahedralize_desc const* desc);
```

`DELAUNAY`: the Delaunay tetrahedralization of the border's vertices,
exact (an integer lattice of 2^30 steps, Bowyer-Watson with orientation and
insphere predicates in 384-bit integers behind a floating filter, each
cavity grown until its border faces the new point, so degenerate input
never leaves flat cells). The border is recovered in it: coplanar neighbours
flipped to the Delaunay choice in their plane, missing edges split at their
midpoints, missing triangles at their circumcentres. The inside is flooded
from the border, filled from a body-centred cubic lattice, refined at the
circumcentres of cells past the radius-edge bound with the border as walls,
then improved by 3-2 and 2-3 flips and smoothing on the worst dihedral
angle. The border's triangles are kept (split only where recovery needed
it), and the volume equals the surface's to rounding. `border_only` keeps
the border's vertices alone. `ENVELOPE` (triangle soups, fTetWild) and
`STUFFING` are reserved.

## Hexahedra

```c
qaws_status qaws_volume_mesh_hexahedralize(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report,
	qaws_mesh_3d const* boundary, qaws_hexahedralize_desc const* desc);
```

`GRID` (Maréchal 2009; Gao, Shen, Panozzo 2019): a grid over the box in
whole cells per axis, the cubes inside kept and made a manifold set, then a
layer of hexahedra padding the staircase to the surface. Each border vertex
of the staircase goes onto the surface's PN reference, onto a crease where
the staircase turns once and onto a corner where it turns twice, so boxes
keep their edges and corners. Smoothing follows, each move kept when it
raises the worst scaled Jacobian.

`FRAME_FIELD` (an integer-grid map, after CubeCover): the inside
tetrahedralized at half the cell size, its frame field combed into one
gauge, then a map fitted by three Poisson solves so that its gradient
follows the frames at the cell size, each border patch (facets whose normal
follows one axis) held on a whole-number plane of that axis. The hexahedra
are the unit cubes of the integer lattice in the map's image, their
vertices mapped back: a turned box comes out as a grid aligned with it, a
bent or twisted bar as a bent or twisted grid. Fields with singular edges
(a sphere, a round torus) or that turn round a hole give
`QAWS_STATUS_UNSUPPORTED_OPERATION`. `POLYCUBE` is reserved.

## Boundary layers

```c
qaws_status qaws_volume_mesh_grow_layers(qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report,
	qaws_mesh_3d const* surface, qaws_boundary_layers_desc const* desc);
```

Thin cells grown from a surface in `layer_count` layers, the first
`first_height` thick and each the one before times `growth` (0: 1.2), as
flow solvers want along walls. Triangles become stacks of prisms, quads
stacks of hexahedra; a closed surface turned outward grows inward unless
`outward` is set, an open surface grows too. Vertex v of layer k is
`k * vertex_count + v`, so the last layer is an offset copy of the surface
(the border of what a tetrahedralization of the rest would fill).

Each vertex moves along one direction: its faces' normals weighted by
their angles, turned until every face round it sees it, its height
stretched so the layers stay as thick across the faces (a box's corner
moves along its diagonal by sqrt 3 times the height: the layers of a box
are boxes). A stack is thinned where it would reach the surface (at most
0.4 of the distance along its direction) or its neighbours' (directions
converging in a concave region), neighbours kept within 1.25 of each
other, and halved where a cell would turn over: every corner's Jacobian
is positive. `report.moves` counts the stacks thinned, `iterations` the
halving rounds. Test 156.

## Example: a NURBS torus, its surface, tetrahedra and hexahedra

```c
qaws_generate_desc g = { 0 };
qaws_tetrahedralize_desc td = { 0 };
qaws_hexahedralize_desc hd = { 0 };
qaws_mesh_3d* surface = NULL;
qaws_volume_mesh *tets = NULL, *hexes = NULL;
qaws_volume_quality q;

g.target_length = 0.2;                                   /* triangles on the torus, exact to rounding */
qaws_mesh_generate_surface(&surface, NULL, torus, &g);
qaws_volume_mesh_tetrahedralize(&tets, NULL, surface, &td);
qaws_volume_mesh_measure_quality(&q, tets);              /* q.min_dihedral, q.max_radius_edge, q.volume */
hd.cell_size = 0.22;
qaws_volume_mesh_hexahedralize(&hexes, NULL, surface, &hd);
```

The showcase `qaws_meshing_showcase` draws figures 1 to 18: the corner
store, references, remeshing, quads, generation, tetrahedra, hexahedra,
hexahedra along frame fields, measures, angle-bounded LOCAL, CVT,
smoothest cross fields, quads by min-cost flow, quad finishing, GLOBAL
quads, guidance, local feature size and repair
(`examples/meshing_showcase/`); an argument listing figure keys draws
those alone (1 to 9, then a to i for figures 10 to 18: `5c` draws 5 and
12).

## Known limits

- Remeshing: the angle bounds are not reached on parts thinner than the
  target or next to sharp corners (exempt sectors aside), and on very
  coarse meshes the angle pass can end worse than the default; the
  distance guard is sampled (a face may bulge past it between samples);
  CVT keeps the crease spacing LOCAL gave it. To a metric: the metric is
  taken as given (no gradation of its sizes), the Hessian of a border
  vertex is fitted on one side, and a metric bent across a crease is
  measured in each vertex's tangent plane only.
- Quads: a cone that no vertex of its face can reach cheaply stays off
  its vertices (a 3-5 pair beside it); positions are refitted by sweeps,
  not a global solve; the field's singularities are those of a smooth field
  (or one drawn to curvature), not placed for the fewest irregular
  vertices; the pipeline is sensitive to tiny changes of the input (a
  scale of 1e-12 can change the quads; rounding the input to float, 1e-7,
  moved a painted density's fine spacing by 22% and a torus from 275 to
  230 quads, so five quality checks of 125 to 137 fail in single
  precision). GLOBAL keeps
  the flow's integers whenever they give a map, even a poorer one than its
  own quantization would (the box at 0.75 times its target: 2 irregular
  vertices, SJ min 0.03, where the quantized integers give 0 and 0.49;
  the two are not both run, to keep the time); on the blob the quantized
  integers, lower in the map's energy, come out with a larger angle
  deviation (3.5 against 2.7), so they are only the second try. The
  quantization moves cones only: the isovalues of creases (parallel
  creases a step too close or too far) and a torus's periods stay the
  flow's, a corner whose turn fixes no lattice point stops it, and it
  gives up past 48 cones. Paths to the border exist for singular faces,
  not for corners. INSTANT's mesh is still returned when neither set of
  integers leaves a fold-free map that winds once round every free vertex.
- Boundary layers: the distance ahead is found against every face (no
  hierarchy: slow past 10^5 faces); thinning is per stack, so a thin wall
  meshed coarsely across its edges (a vertex at a rim moving diagonally
  above a side's next vertex) is thinned by the rim's spacing, not the
  wall's thickness; no layers on a vertex no direction serves (a spike).
- Generation: vertices on a trimmed face's trim curve are snapped onto the
  surface, not the curve. Curved elements: nodes by the closest point of
  the straight element's (not equally spaced along the curve; no check of
  the curved element's Jacobian); a node on an edge between two faces lies
  on the face of the element that numbered it; meshes on a mesh reference
  (remeshed, not generated) get none; volume elements are not curved.
- Tetrahedra: improvement by flips and smoothing only (no edge removal or
  sliver exudation); no protection of sharp angles between border
  triangles, where recovery may stop (`QAWS_STATUS_NOT_CONVERGED`).
- Hexahedra: GRID is a uniform grid (no octree), features kept only where
  the staircase turns near them; FRAME_FIELD needs a field without singular
  edges.
- Frame fields: the singular edges are those of the smoothest field, not
  constrained to be meshable (Liu, Bommes 2023).
- Repair: welding is transitive (a chain of points each within the
  tolerance merges wider); a gap is closed by moving the vertices onto the
  other border, not by bridging strips between borders farther apart than
  the tolerance (those become holes, filled); loops past 500 vertices are
  filled by a fan from their centroid (folded where the loop is far from
  star-shaped round it); fills follow the loop and the ring round it, not
  features (a hole across a crease is filled smooth); nesting is judged at
  up to 16 vertices of a part (nested when all are inside), so a part
  crossing another stays outward unless those samples all fall inside
  (farq's self-union sorts crossing parts); a part with no consistent
  orientation is left as it is.
