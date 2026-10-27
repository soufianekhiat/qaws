# Meshes for hosts: plan

The curves and surfaces plan (`docs/host_compatibility_plan.md`) is done.
This one brings meshes to the same rules where a host (Atelier) holds the
geometry: its vertex and index buffers, f32 or f64, interleaved or planar,
16- or 32-bit indices, changing every frame.

## 1. What borrowing means for meshes

A curve view is evaluated straight from host memory. A mesh operation
builds derived structures (half-edges, a BVH, prisms) that the host does
not have, so the rule here is the "prepared" one of the curves plan:

- derived data is built once and owned (connectivity, the BVH tree);
- positions are read from the host's buffers at each use, converted as
  they are read, with no staging array on the host's side and no
  allocation for an update;
- the host's memory is never written, freed or kept past the call, except
  by objects documented to borrow it for their life.

Heavy operations that make new meshes (Booleans, remeshing, sweeps) take a
`qaws_mesh_3d`; making one from host buffers is one conversion pass.

## 2. Work

1. **Mesh views.** `qaws_mesh_view`: positions as a `qaws_array_view` (3
   components, f32 or f64, any stride or planar), indices as 16- or 32-bit
   integers with a stride, a triangle count; `qaws_mesh_view_validate`.
   `qaws_mesh_3d_create_from_view` (an owned mesh, the indices checked) and
   `qaws_mesh_3d_set_positions_from_view` (a frame's positions into an
   existing mesh: no allocation).
2. **BVH from views.** `qaws_mesh_bvh_create_view` and
   `qaws_mesh_bvh_refit_view`: a skinned or simulated mesh refit every
   frame from its vertex buffer, queries (closest, rays, winding, signed
   distance) as before.
3. **Borrowed vertex motions.** `qaws_motion_create_vertices_borrowed`:
   each vertex on its chord from one positions view to another (the last
   frame's and this one's), read at each evaluation; the motion keeps the
   two views for its life. Contact detection with it gets the exact times
   of straight vertex motions without copying either frame.
4. **Scratch.** The temporaries of contact detection (none through the
   scratch allocator today), of the implicit and explicit sweeps (partly)
   through `qaws_internal_scratch_malloc`, results and objects through the
   persistent one, as the rest of the library.
5. **An executor for contact detection.** The per-triangle prism tests of a
   slab through `qaws_sweep_desc.executor` (the same results on any number
   of threads).

Tests stay in f64 and light (a few hundred triangles); the f32 data paths
are covered by views of f32 buffers in an f64 build.

## 3. Status

| Item | State | Commits, tests |
|---|---|---|
| Mesh views | `qaws_mesh_view`, validation, `qaws_mesh_3d_create_from_view`, `qaws_mesh_3d_set_positions_from_view` (`docs/api/mesh_view.md`) | 185 |
| BVH from views | `qaws_mesh_bvh_create_view`, `qaws_mesh_bvh_refit_view`; `qaws_mesh_bvh_refit` in place (it allocated a new array per call) | 185, 159 |
| Borrowed vertex motions | `qaws_motion_create_vertices_borrowed`: the two frames read at each evaluation, exact contact times; a rewritten frame seen at the next call | 185 |
| Scratch | contact detection, the explicit and implicit sweeps: every temporary through the scratch allocator (lists first grown there), sizes checked; a leak found by the scratch arena (a sweep of a mesh without features kept a block per call) fixed: suites 180, 181, 183 to 185 hold nothing with `QAWS_TEST_SCRATCH=1` (180 held 6, 181 3, 184 1) | 180, 181, 183-185 |
| Executor for contact | the moving mesh's triangles in 64 shares through `qaws_sweep_desc.executor`, each with its own candidates (the per-obstacle stamp array dropped: the BVH's box query gives each element once), the first contact the earliest share's, pairs sorted: the same reports bit for bit | 183 |

Nothing remains of this plan.
