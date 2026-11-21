# API Reference: qaws_brep.h

B-rep (boundary representation) shells: the topology of a solid or sheet
body as faces bounded by oriented edges between vertices.

- Each face is a (possibly trimmed) surface with oriented boundary edges.
- An edge is shared by two faces (manifold) or is a boundary edge.
- Vertices are shared between edges.
- The topology is half-edge style: each edge has a face on its left
  (used with orientation +1) and on its right (orientation -1).
- Surfaces and curves are borrowed: they must outlive the shell.

Names are shown without precision suffix; each exists as `_f32` and `_f64`
(see [precision.md](precision.md)). Concepts: [../brep.md](../brep.md).

---

## Types

```c
typedef struct qaws_brep_shell qaws_brep_shell;   /* opaque */

typedef unsigned int qaws_brep_face_id;
typedef unsigned int qaws_brep_edge_id;
typedef unsigned int qaws_brep_vertex_id;
typedef unsigned int qaws_brep_halfedge_id;

#define QAWS_BREP_INVALID_ID ((unsigned int)0xFFFFFFFF)

typedef struct qaws_brep_face_desc {
	qaws_surface const* surface;   /* the face geometry (borrowed) */
	int is_reversed;               /* 1 = normal points inward */
} qaws_brep_face_desc;

typedef struct qaws_brep_edge_desc {
	qaws_brep_vertex_id v_start;
	qaws_brep_vertex_id v_end;
	qaws_curve const* curve_3d;    /* 3D edge curve (borrowed, may be NULL) */
} qaws_brep_edge_desc;
```

Ids are indices in creation order, starting at 0. Half-edge `2e` is edge
`e`'s use on its left face (orientation +1, from `v_start` to `v_end`) and
`2e + 1` its use on its right face (orientation -1, from `v_end` to
`v_start`).

### qaws_brep_validation_flags

| Flag | Value | Meaning |
|---|---|---|
| `QAWS_BREP_VALID` | 0 | no problem found |
| `QAWS_BREP_NOT_MANIFOLD` | 1 << 0 | an edge shared by more than 2 faces |
| `QAWS_BREP_NOT_CLOSED` | 1 << 1 | boundary edges exist |
| `QAWS_BREP_NOT_ORIENTABLE` | 1 << 2 | inconsistent face orientations |
| `QAWS_BREP_HAS_DEGENERATE` | 1 << 3 | a zero-area face or zero-length edge |

---

## Shell lifetime

```c
qaws_status qaws_brep_create(qaws_brep_shell** out_shell);
void        qaws_brep_destroy(qaws_brep_shell* shell);
```

Creates an empty shell / frees it. Destroying does not free the borrowed
surfaces and curves.

## Building

```c
qaws_status qaws_brep_add_vertex(qaws_brep_vertex_id* out_id, qaws_brep_shell* shell, qaws_vec3 position);

qaws_status qaws_brep_add_edge(qaws_brep_edge_id* out_id, qaws_brep_shell* shell,
	qaws_brep_edge_desc const* desc);

qaws_status qaws_brep_add_face(
	qaws_brep_face_id* out_id,
	qaws_brep_shell* shell,
	qaws_brep_edge_id const* edge_ids,
	int const* edge_orientations,   /* +1 or -1 per edge */
	unsigned int edge_count,
	qaws_brep_face_desc const* desc);
```

- `qaws_brep_add_edge`: the edge's vertices must exist
  (`QAWS_STATUS_INVALID_ARGUMENT` otherwise).
- `qaws_brep_add_face`: the face's loop is the listed edges, each used with
  orientation +1 (the face is on the edge's left) or -1 (on its right). The
  edge and orientation arrays are copied. Returns
  `QAWS_STATUS_INVALID_ARGUMENT` for an unknown edge, an orientation other
  than +1 or -1, or an edge whose side is already taken by another face (so
  an edge never has more than two faces); the shell is then unchanged.
- All three return `QAWS_STATUS_ALLOCATION_FAILURE` when memory runs out.

## Topology queries

```c
unsigned int qaws_brep_get_vertex_count(qaws_brep_shell const* shell);
unsigned int qaws_brep_get_edge_count(qaws_brep_shell const* shell);
unsigned int qaws_brep_get_face_count(qaws_brep_shell const* shell);

qaws_status qaws_brep_get_vertex_position(qaws_vec3* out_position, qaws_brep_shell const* shell,
	qaws_brep_vertex_id id);
qaws_status qaws_brep_get_face_surface(qaws_surface const** out_surface, qaws_brep_shell const* shell,
	qaws_brep_face_id id);
qaws_status qaws_brep_get_edge_vertices(qaws_brep_vertex_id* out_v_start, qaws_brep_vertex_id* out_v_end,
	qaws_brep_shell const* shell, qaws_brep_edge_id id);
qaws_status qaws_brep_get_edge_faces(qaws_brep_face_id* out_face_left, qaws_brep_face_id* out_face_right,
	qaws_brep_shell const* shell, qaws_brep_edge_id id);
qaws_status qaws_brep_get_face_edges(qaws_brep_edge_id* out_edges, unsigned int* out_count,
	qaws_brep_shell const* shell, unsigned int capacity, qaws_brep_face_id id);
```

- `qaws_brep_get_edge_faces`: a side with no face is `QAWS_BREP_INVALID_ID`.
- `qaws_brep_get_face_edges`: the face's edges in the order given to
  `qaws_brep_add_face`; `*out_count` is always set, and `out_edges` may be
  NULL to count only.
- Statuses: `QAWS_STATUS_INVALID_ARGUMENT` on a NULL pointer,
  `QAWS_STATUS_OUT_OF_RANGE` for an unknown id,
  `QAWS_STATUS_BUFFER_TOO_SMALL` when the list does not fit `capacity`.

## Half-edges

```c
typedef struct qaws_brep_halfedge_info {
	qaws_brep_edge_id edge;         /* h / 2 */
	qaws_brep_face_id face;
	qaws_brep_halfedge_id next;     /* the next in the face's edge list */
	qaws_brep_halfedge_id prev;
	qaws_brep_halfedge_id twin;
	qaws_brep_vertex_id origin;     /* v_start for 2e, v_end for 2e + 1 */
	qaws_brep_vertex_id target;
} qaws_brep_halfedge_info;

qaws_status qaws_brep_get_face_halfedge(qaws_brep_halfedge_id* out_halfedge, qaws_brep_shell const* shell,
	qaws_brep_face_id id);
qaws_status qaws_brep_get_halfedge(qaws_brep_halfedge_info* out_info, qaws_brep_shell const* shell,
	qaws_brep_halfedge_id h);
```

- `qaws_brep_get_face_halfedge`: the half-edge of the face's first edge, as given to `qaws_brep_add_face`;
  `QAWS_BREP_INVALID_ID` for a face without edges.
- `qaws_brep_get_halfedge`: a half-edge no face uses has `face`, `next`, `prev` and `twin` `QAWS_BREP_INVALID_ID`;
  `twin` is `h ^ 1` when a face uses the other side, else `QAWS_BREP_INVALID_ID` (a boundary edge).
- `next` and `prev` follow the face's edge list as given to `qaws_brep_add_face` (cyclically), which is its loop when
  the edges were given in loop order (the shell does not check that they chain). Walking a loop: `next` from the face's
  half-edge until it comes back. Around a vertex of chained loops: `twin` then `next`.
- Statuses: `QAWS_STATUS_INVALID_ARGUMENT` on a NULL pointer, `QAWS_STATUS_OUT_OF_RANGE` for a face past the count or
  `h >= 2 *` the edge count.

## Adjacency queries

```c
qaws_status qaws_brep_get_vertex_edges(qaws_brep_edge_id* out_edges, unsigned int* out_count,
	qaws_brep_shell const* shell, unsigned int capacity, qaws_brep_vertex_id id);
qaws_status qaws_brep_get_vertex_faces(qaws_brep_face_id* out_faces, unsigned int* out_count,
	qaws_brep_shell const* shell, unsigned int capacity, qaws_brep_vertex_id id);
```

The edges that start or end at the vertex, and the faces on those edges
(each listed once). The output array may be NULL to count only; for
`qaws_brep_get_vertex_faces` that count does not remove repeated faces, so
it is an upper bound to size the buffer with. Statuses as above.

## Validation

```c
qaws_status qaws_brep_validate(unsigned int* out_flags, qaws_brep_shell const* shell);
qaws_status qaws_brep_euler_characteristic(int* out_chi, qaws_brep_shell const* shell);
```

- `qaws_brep_validate` sets `*out_flags` to a combination of
  `qaws_brep_validation_flags` (`QAWS_BREP_VALID` when none apply).
- `qaws_brep_euler_characteristic` gives V - E + F; 2 for a closed genus-0
  shell.

## Measures

```c
qaws_status qaws_brep_compute_volume(qaws_scalar* out_volume, qaws_brep_shell const* shell);
qaws_status qaws_brep_compute_surface_area(qaws_scalar* out_area, qaws_brep_shell const* shell);
```

- **Volume** is meaningful for closed shells only. It uses the divergence
  theorem, V = 1/3 of the integral of x . n dA over the faces, sampled on a
  midpoint grid over each face's (u, v) domain, with the normal flipped on
  reversed faces; the absolute value is returned.
- **Surface area** is the sum of the faces' areas
  (`qaws_surface_compute_area`, [inspect.md](inspect.md)); an error
  from a face is returned.
- Faces without a surface are skipped by both.
