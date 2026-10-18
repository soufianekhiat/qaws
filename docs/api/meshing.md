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
| `qaws_mesh_field.h` | size fields; cross fields on surfaces; frame fields in volumes |
| `qaws_remesh.h` | surface remeshing (LOCAL) |
| `qaws_quadrangulate.h` | quad meshes (INSTANT) |
| `qaws_mesh_generate.h` | meshes of surfaces and shells (REMESH) |
| `qaws_volume_mesh.h` | `qaws_volume_mesh`: tetrahedra, hexahedra, prisms, pyramids |
| `qaws_tetrahedralize.h` | tetrahedra filling a closed surface (DELAUNAY) |
| `qaws_hexahedralize.h` | hexahedra filling a closed surface (GRID) |

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

qaws_status qaws_mesh_3d_create(qaws_mesh_3d_desc const* desc, qaws_mesh_3d** out_mesh);
qaws_status qaws_mesh_3d_get_faces(qaws_mesh_3d const* mesh, unsigned int* out_face_vertices, unsigned int* out_face_sizes);
qaws_status qaws_mesh_3d_get_face_adjacency(qaws_mesh_3d const* mesh, unsigned int* out_faces);
qaws_status qaws_mesh_3d_split_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position,
	qaws_mesh_location const* footpoint, unsigned int* out_vertex);
qaws_status qaws_mesh_3d_collapse_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position,
	qaws_mesh_location const* footpoint);
qaws_status qaws_mesh_3d_flip_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b);
qaws_status qaws_mesh_3d_split_face(qaws_mesh_3d* mesh, unsigned int face, qaws_scalar const* position,
	qaws_mesh_location const* footpoint, unsigned int* out_vertex);
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
qaws_status qaws_mesh_reference_create(qaws_mesh_3d const* mesh, qaws_mesh_reference_desc const* desc, qaws_mesh_reference** out);
qaws_status qaws_mesh_reference_evaluate(qaws_mesh_reference const* r, qaws_mesh_location const* locations, unsigned int count,
	qaws_scalar* out_positions, qaws_scalar* out_normals, qaws_scalar* out_curvatures);
qaws_status qaws_mesh_reference_project(qaws_mesh_reference const* r, qaws_scalar const* points, unsigned int count,
	qaws_mesh_location const* hints, qaws_mesh_location* out_locations, qaws_scalar* out_distances);
qaws_status qaws_mesh_reference_walk(qaws_mesh_reference const* r, qaws_mesh_location from, qaws_scalar const* tangent,
	qaws_mesh_location* out_to);
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
qaws_status qaws_size_field_create_uniform(qaws_scalar length, qaws_size_field** out);
qaws_status qaws_size_field_create_curvature(qaws_mesh_reference const* r, qaws_curvature_size_desc const* desc, qaws_size_field** out);
qaws_status qaws_mesh_3d_direction_field(qaws_mesh_3d const* mesh, qaws_direction_field_desc const* desc, qaws_scalar* out_directions,
	int* out_face_index);
qaws_status qaws_volume_mesh_frame_field(qaws_volume_mesh const* mesh, qaws_frame_field_desc const* desc, qaws_scalar* out_frames,
	unsigned char* out_singular);
```

- **Size fields** give the target edge length at each point: uniform; by
  curvature from chord and angle tolerances (clamped, graded so sizes grow
  at most by `gradation` times the distance); per reference vertex; a
  callback; the smaller of two.
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
qaws_status qaws_mesh_remesh(qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc, qaws_mesh_3d** out, qaws_mesh_report* report);
qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_3d** out,
	qaws_mesh_report* report);
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
  are read, so quads follow sharp edges. `PAIRING`, `SPLIT`, `GLOBAL`
  (integer-grid maps) and `LAYOUT` are reserved.

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

## Measures

```c
qaws_status qaws_mesh_3d_measure_shape(qaws_mesh_3d const* mesh, qaws_shape_quality* out, unsigned int* out_ideal);
qaws_status qaws_mesh_3d_measure_deviation(qaws_mesh_3d const* mesh, qaws_mesh_reference const* reference,
	qaws_deviation_desc const* desc, qaws_mesh_deviation* out);
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
qaws_status qaws_mesh_generate_surface(qaws_surface const* surface, qaws_generate_desc const* desc, qaws_mesh_3d** out,
	qaws_mesh_report* report);
qaws_status qaws_mesh_generate_shell(qaws_shell_atlas const* shell, qaws_generate_desc const* desc, qaws_mesh_3d** out,
	qaws_mesh_report* report);
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

## Volume meshes: qaws_volume_mesh

```c
typedef struct qaws_volume_mesh_desc {
	qaws_scalar const* positions; unsigned int vertex_count;
	unsigned int const* cell_vertices; qaws_cell_type const* cell_types; unsigned int cell_count;   /* types NULL: tetrahedra */
} qaws_volume_mesh_desc;

qaws_status qaws_volume_mesh_create(qaws_volume_mesh_desc const* desc, qaws_volume_mesh** out_mesh);
qaws_status qaws_volume_mesh_get_cells(qaws_volume_mesh const* mesh, unsigned int* out_cell_vertices, qaws_cell_type* out_types);
qaws_status qaws_volume_mesh_get_cell_adjacency(qaws_volume_mesh const* mesh, unsigned int* out_cells);
qaws_status qaws_volume_mesh_get_boundary(qaws_volume_mesh const* mesh, qaws_mesh_3d** out_surface);
qaws_status qaws_volume_mesh_measure_quality(qaws_volume_mesh const* mesh, qaws_volume_quality* out_quality);
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
qaws_status qaws_volume_mesh_tetrahedralize(qaws_mesh_3d const* boundary, qaws_tetrahedralize_desc const* desc,
	qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report);
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
qaws_status qaws_volume_mesh_hexahedralize(qaws_mesh_3d const* boundary, qaws_hexahedralize_desc const* desc,
	qaws_volume_mesh** out_mesh, qaws_mesh_report* out_report);
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

## Example: a NURBS torus, its surface, tetrahedra and hexahedra

```c
qaws_generate_desc g = { 0 };
qaws_tetrahedralize_desc td = { 0 };
qaws_hexahedralize_desc hd = { 0 };
qaws_mesh_3d* surface = NULL;
qaws_volume_mesh *tets = NULL, *hexes = NULL;
qaws_volume_quality q;

g.target_length = 0.2;                                   /* triangles on the torus, exact to rounding */
qaws_mesh_generate_surface(torus, &g, &surface, NULL);
qaws_volume_mesh_tetrahedralize(surface, &td, &tets, NULL);
qaws_volume_mesh_measure_quality(tets, &q);              /* q.min_dihedral, q.max_radius_edge, q.volume */
hd.cell_size = 0.22;
qaws_volume_mesh_hexahedralize(surface, &hd, &hexes, NULL);
```

The showcase `qaws_meshing_showcase` draws figures 1 to 14: the corner
store, references, remeshing, quads, generation, tetrahedra, hexahedra,
hexahedra along frame fields, measures, angle-bounded LOCAL, CVT,
smoothest cross fields, quads by min-cost flow and quad finishing
(`examples/meshing_showcase/`); an argument listing figure keys draws
those alone (1 to 9, then a to e for figures 10 to 14: `5c` draws 5 and
12).

## Known limits

- Remeshing: the angle bounds are not reached on parts thinner than the
  target or next to sharp corners (exempt sectors aside), and on very
  coarse meshes the angle pass can end worse than the default; the
  distance guard is sampled (a face may bulge past it between samples);
  CVT keeps the crease spacing LOCAL gave it.
- Quads: a cone that no vertex of its face can reach cheaply stays off
  its vertices (a 3-5 pair beside it); positions are refitted by sweeps,
  not a global solve; the field's singularities are those of a smooth field
  (or one drawn to curvature), not placed for the fewest irregular
  vertices; the pipeline is sensitive to tiny changes of the input (a
  scale of 1e-12 can change the quads, not their quality).
- Generation: vertices on a trimmed face's trim curve are snapped onto the
  surface, not the curve.
- Tetrahedra: improvement by flips and smoothing only (no edge removal or
  sliver exudation); no protection of sharp angles between border
  triangles, where recovery may stop (`QAWS_STATUS_NOT_CONVERGED`).
- Hexahedra: GRID is a uniform grid (no octree), features kept only where
  the staircase turns near them; FRAME_FIELD needs a field without singular
  edges.
- Frame fields: the singular edges are those of the smoothest field, not
  constrained to be meshable (Liu, Bommes 2023).
