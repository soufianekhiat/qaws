# API Reference: qaws_surface_subdiv.h

Subdivision surfaces: the limit surface of a control mesh refined by
Catmull-Clark or Loop subdivision. The result is a `qaws_surface` of kind
`QAWS_SURFACE_KIND_SUBDIV` (see [surface.md](surface.md)); its triangle mesh
can also be read directly for rendering.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)).

---

## Model

- **Catmull-Clark**: quad-dominant meshes, based on cubic B-spline
  refinement; C2 away from extraordinary vertices.
- **Loop**: triangle meshes, based on quartic box-spline refinement; C2 away
  from extraordinary vertices.
- Faces are triangles or quads. Mixed meshes use Catmull-Clark, with
  triangles treated as degenerate quads.
- The mesh is subdivided `subdivision_level` times at creation and the
  result stored on a regular grid. Evaluation maps (u, v) in [0,1] x [0,1] to
  that grid by bilinear interpolation.
- The control mesh is copied: the input arrays need not outlive the call.

## Types

```c
typedef enum qaws_subdiv_scheme {
	QAWS_SUBDIV_CATMULL_CLARK = 0,
	QAWS_SUBDIV_LOOP
} qaws_subdiv_scheme;

typedef struct qaws_surface_subdiv_desc {
	qaws_vec3 const* vertices;          /* vertex positions */
	unsigned int vertex_count;
	unsigned int const* face_indices;   /* flat array of vertex indices, face after face */
	unsigned int const* face_sizes;     /* vertices per face (3 or 4) */
	unsigned int face_count;
	qaws_subdiv_scheme scheme;
	unsigned int subdivision_level;     /* subdivision iterations; 0 = 3 */
} qaws_surface_subdiv_desc;
```

## Functions

```c
qaws_status qaws_surface_create_subdiv(qaws_surface** out_surface, qaws_surface_subdiv_desc const* desc);
```

Subdivides the control mesh and creates the surface. Free it with
`qaws_surface_destroy`.

```c
qaws_status qaws_surface_subdiv_get_mesh(
	qaws_vec3 const** out_vertices,
	qaws_vec3 const** out_normals,
	unsigned int const** out_triangles,
	unsigned int* out_vertex_count,
	unsigned int* out_triangle_count,
	qaws_surface const* surface);
```

The subdivided triangle mesh. The pointers point into the surface's own data
and stay valid until the surface is destroyed; do not free them.
`out_triangles` holds `*out_triangle_count` counter-clockwise index triples
into the vertex array; `out_normals` holds one smooth normal per vertex.
