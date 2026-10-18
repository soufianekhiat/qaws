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
- **Frame fields** in tetrahedral volumes: an octahedral frame per cell,
  the border's normals held, smoothed over the stored facet adjacency;
  singular edges reported per cell as a mask of its six edges.

## Remeshing and quads

```c
qaws_status qaws_mesh_remesh(qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc, qaws_mesh_3d** out, qaws_mesh_report* report);
qaws_status qaws_mesh_quadrangulate(qaws_mesh_3d const* mesh, qaws_quadrangulate_desc const* desc, qaws_mesh_3d** out,
	qaws_mesh_report* report);
```

- **LOCAL** remeshing (Botsch and Kobbelt 2004): split, collapse, flip,
  relax along the reference (walked, never projected), toward a uniform
  length or a size field. Feature corners stay, feature lines slide, no
  operation leaves a poor triangle poorer. Attributes are transferred from
  the reference's mesh at the footpoints, seams kept per side. `CVT` and
  `FIELD_ALIGNED` are reserved.
- **INSTANT** quads (Jakob et al. 2015): a dense remesh, a cross field and
  a lattice position field, the mesh read off them; quad-dominant, or all
  quads with `pure`. Lattice points are unfolded across creases before
  their offsets are read, so quads follow sharp edges. `PAIRING`, `SPLIT`,
  `GLOBAL` (QuadriFlow's consistency, integer-grid maps) and `LAYOUT` are
  reserved.

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
raises the worst scaled Jacobian. `POLYCUBE` and `FRAME_FIELD` (the frame
field's integer-grid map) are reserved; the frame field itself is
`qaws_volume_mesh_frame_field`.

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

The showcase `qaws_meshing_showcase` draws figures 1 to 7: the corner
store, references, remeshing, quads, generation, tetrahedra and hexahedra
(`examples/meshing_showcase/`).

## Known limits

- Quads: the read-off keeps the position field's local inconsistencies (a
  few extra singularities); QuadriFlow's global consistency is `GLOBAL`.
- Generation: vertices on a trimmed face's trim curve are snapped onto the
  surface, not the curve.
- Tetrahedra: improvement by flips and smoothing only (no edge removal or
  sliver exudation); no protection of sharp angles between border
  triangles, where recovery may stop (`QAWS_STATUS_NOT_CONVERGED`).
- Hexahedra: a uniform grid (no octree), features kept only where the
  staircase turns near them.
- Frame fields: local smoothing; on curved shapes the singular edges are
  more numerous and less regular than a global optimization gives.
