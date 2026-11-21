# B-rep Shells and Trimmed Surfaces

This page covers the two structures qaws uses to build bodies out of surfaces:

- A **B-rep shell** (`qaws_brep.h`) records which faces, edges and vertices make up a body and how they connect.
- A **trimmed surface** (`qaws_surface_trim.h`) restricts one surface to a region of its parameter space.

Signatures are in [api/brep.md](api/brep.md) and [api/surface_trim.md](api/surface_trim.md).

> Names are given without their precision suffix. Each one exists as `name_f32` and `name_f64`; see
> [precision.md](precision.md).

## The shell

A `qaws_brep_shell` represents a solid or a sheet as a set of faces:

- A **vertex** is a 3D position.
- An **edge** joins a start vertex and an end vertex. It can optionally carry a 3D curve (`curve_3d`).
- A **face** is a surface, optionally trimmed, with a flag `is_reversed` that marks a normal pointing inward. It also
  holds an ordered list of its edges, each with an orientation of +1 or -1.

All ids are dense indices in order of creation, and `QAWS_BREP_INVALID_ID` marks "none".

### How faces and edges connect

Each edge has two sides:

- A face that uses the edge with orientation +1 sits on its left.
- A face that uses it with orientation -1 sits on its right.

`qaws_brep_add_face` refuses (`INVALID_ARGUMENT`) a face whose edge side is already taken, so an edge never has more
than two faces. An edge with one face is a boundary edge. These two sides are what the header calls the half-edges of
the edge. No query takes a half-edge id.

### Building a shell

Build a shell bottom-up:

1. `qaws_brep_create` makes an empty shell.
2. Add the vertices.
3. Add the edges, which reference existing vertices.
4. Add the faces, which reference existing edges.

The shell copies its own records (positions, edge lists, orientations). It **borrows** the surfaces and curves: they
must outlive the shell, and `qaws_brep_destroy` does not free them. Shell memory comes from the allocator chain in
[memory.md](memory.md). There is no `_ex` variant.

### Queries

- **Counts:** the number of vertices, edges and faces.
- **One entity:** a vertex position, a face's surface, an edge's vertices, an edge's left and right faces.
- **Adjacency lists:** a face's edges, the edges at a vertex, the faces at a vertex (without duplicates).

List queries follow the library's buffer convention. They return the count through `out_count`, and if the output
array is NULL they only count. When the array is smaller than the list, they return `BUFFER_TOO_SMALL`. An id that is
out of range returns `OUT_OF_RANGE`.

### Validation and measures

`qaws_brep_validate` returns a set of `qaws_brep_validation_flags`:

- `NOT_CLOSED`: some edge lacks a face on one side.
- `NOT_MANIFOLD`: some edge has the same face on both sides.
- `NOT_ORIENTABLE`: two faces traverse a shared edge in the same direction.
- `HAS_DEGENERATE`: an edge is shorter than `QAWS_EPSILON`, or a face has fewer than 3 edges.

The measures:

- `qaws_brep_euler_characteristic` returns V - E + F, which is 2 for a closed genus-0 shell.
- `qaws_brep_compute_surface_area` sums `qaws_surface_compute_area` over the faces.
- `qaws_brep_compute_volume` applies the divergence theorem on a 20 x 20 midpoint grid over each face's parameter
  domain. A reversed face counts with its normal flipped. The function returns the absolute value. It is meaningful
  only for closed shells, and its accuracy is that of the grid. It integrates over each face's full domain, so a trimmed
  face counts as its whole square.

## Trimmed surfaces

`qaws_surface_create_trimmed` wraps a base surface with trim loops in (u, v):

- A `qaws_trim_loop` is a closed sequence of 2D curves.
- An outer loop (`is_outer = 1`) keeps the points inside it. A hole loop removes the points inside it.
- With no outer loop, the whole domain is kept.

The trimmed surface copies the loop arrays but **borrows** the base surface and the trim curves. Both must outlive it.
The `_ex` variant takes an allocator (see [memory.md](memory.md)).

### Domain and evaluation

The trimmed surface's domain is `[0, 1] x [0, 1]`:

- Evaluation clamps (u, v) to that square and evaluates the base there.
- Evaluation does not test the trim. Every function that evaluates a surface sees the whole base over the unit square.
- A base whose domain is not the unit square is reached only where the two overlap.

### Testing a point against the trim

`qaws_surface_trim_contains(u, v)` tests a point against the trim:

1. Each loop is sampled at 100 points per curve into a polygon.
2. The point passes a loop when the polygon's winding number around it is non-zero.
3. The point is inside when it is inside every outer loop and inside no hole.

The test is approximate near the loops, to the resolution of that sampling. The function allocates nothing and is
safe to call from several threads.

## How the two fit together

A B-rep face takes any `qaws_surface`, and a trimmed surface is the usual face geometry. The pieces that
`qaws_surface_shatter_create` returns are trimmed surfaces over the cut surface's base, with polygonal loops (see
[shatter.md](shatter.md)), so they can be used as faces directly. Shatter also reads the trims of its trimmed input
surfaces: seam points outside a trim are dropped.
