# Meshes from host buffers

`qaws_mesh_view.h`. A host keeps its geometry in its own buffers: a vertex
buffer (f32 or f64, positions interleaved with other attributes at any
stride, or planar) and an index buffer (16- or 32-bit, possibly padded).
These functions read them directly: no staging array on the host's side,
no allocation to update positions frame after frame.

Derived data (a mesh's half-edges, a BVH's tree) is built once and owned.
Positions are read from the view at each use, converted as they are read.
qaws never writes, frees or keeps the host's memory, except a borrowed
vertex motion, which reads its two views at every evaluation for its whole
life (see the lifetime contract in `view.md`).

```c
qaws_mesh_view v = { 0 };
v.positions.data = vertices;                  /* struct { float pos[3]; float uv[2]; } */
v.positions.type = QAWS_ELEMENT_F32;
v.positions.components = 3;
v.positions.count = vertex_count;
v.positions.stride = sizeof *vertices;
v.indices = indices;                          /* uint16_t, 3 per triangle */
v.index_type = QAWS_INDEX_U16;
v.triangle_count = triangle_count;

qaws_mesh_bvh_create_view(&bvh, &v, NULL);    /* once */
/* each frame, after the host wrote its vertex buffer: */
qaws_mesh_bvh_refit_view(bvh, &v.positions);  /* no allocation, the tree kept */
```

| Function | Does |
|---|---|
| `qaws_mesh_view_validate` | positions a valid view of 3 components, all finite; a known index type; a stride of at least 3 indices (or 0, tight); every index below the vertex count |
| `qaws_mesh_3d_create_from_view` | an owned `qaws_mesh_3d` of the triangles (one pass) |
| `qaws_mesh_3d_set_positions_from_view` | a frame's positions into a mesh of as many vertices; no allocation; a value not finite leaves it as it was |
| `qaws_mesh_bvh_create_view` | a BVH over the triangles (elements: triangle indices) |
| `qaws_mesh_bvh_refit_view` | a frame's positions into the BVH (in place: `qaws_mesh_bvh_refit` no longer allocates either) |
| `qaws_motion_create_vertices_borrowed` | each vertex on its chord at a constant rate from one positions view (t = 0) to another (t = 1): the last frame's and this one's |

A borrowed vertex motion is a straight vertex motion, so contact detection
(`qaws_mesh_3d_first_contact`, `qaws_mesh_3d_first_self_contact`) gives it
exact times. The host may rewrite either frame between calls: the motion
reads them as they are.

Test 185: an interleaved f32 buffer with 16-bit indices and planar f64
positions with 32-bit indices padded to 4 per triangle make the same box
(volume 1); an index past the vertices, a position not finite and 2
components are refused; a frame's positions moved by (1, 2, 3) reach the
mesh and the BVH (the closest point from (5, 0.5, 0.5) at 4, then
sqrt(17.5)), a refit of another count refused with the BVH unchanged; a
borrowed motion of a box's top rising from 1 to 3 meets a ceiling at 2 at
0.5 exactly, and at 0.25 exactly once the host rewrites its second frame
(the top to 5).
