# Curves and regions across shell faces

`qaws_shell_region.h` carries the surface toolkit (`qaws_surface_region.h`)
over the shared edges of a shell: a cube of patches, a sphere of
hemispheres, a network of NURBS faces.

```c
qaws_shell_atlas* atlas;
qaws_shell_atlas_create(surfaces, n, 0, &atlas);        /* or _create_from_brep */

qaws_shell_curve* c;
qaws_shell_curve_create_geodesic_circle(atlas, centre, r, 0, &c);
qaws_shell_region* region;
qaws_shell_region_create(atlas, &c, 1, &region);
qaws_shell_region_get_area(region, &area);
qaws_shell_region_clip(a, b, QAWS_CLIP_UNION, &u);
qaws_shell_exp(atlas, p, direction, length, &q, &velocity, &crossings);
```

## The atlas

- **Shared edges.** Two faces share an edge where a side of one face's
  parameter domain coincides in space with a side of another's (sampled
  along the side, in either direction).
- **Not edges.** A side along a closed direction is a seam, and a side that
  is one point is a pole; neither counts as an edge.
- **Orientation.** The faces must be oriented alike, with each shared edge
  running opposite ways round the two faces. Otherwise creation returns
  `INVALID_ARGUMENT`.
- **Crossing.** A point crossing an edge goes to the neighbour's side point
  at the same place in space, found by Gauss-Newton along that side.

## Geodesics over edges

- **Integration.** On each face the geodesic is integrated as in
  `qaws_surface_region.h`. Where it leaves the domain, the exit is found by
  bisection on the length.
- **Continuation.** The direction keeps its component along the edge and
  its component across it, the latter turned into the next face's tangent
  plane. This is the straightest continuation: across a smooth join it is
  the geodesic itself, and across a crease (a cube's edge) it unfolds the
  two faces about the edge.
- **Vertices.** A geodesic reaching a face's corner, a vertex of the shell,
  is not continued and returns `OUT_OF_RANGE`.

## Circles and regions

- **Circles.** A geodesic circle is the rays of the exponential map across
  faces. The angles where it changes face are found by bisection. Each
  face's run is a (u, v) Hermite curve through ray ends, with derivatives
  by differences in the angle, refined until its midpoints agree with the
  rays.
- **Regions.** A region is a `qaws_surface_region` on each face. On a face,
  the runs are joined along the face's sides, counter-clockwise from where
  one leaves to where the next comes in (round corners on open faces; along
  a closed direction's side as one circuit). A face without runs is inside
  or outside as its neighbour's side says.
- **Measures and Booleans.** Areas and Booleans go face by face. The
  perimeter leaves out the stretches along shared sides.

## Checks (test 107, f64)

| Check | Agreement |
|---|---|
| Cube of six flat faces: shared edges found | 12 |
| Cube: a straightest ray over an edge against the unfolded line | 2e-12 |
| Cube: circles over one edge and over two by a vertex against pi r² (area) and 2 pi r | 1e-11 |
| Cube: two circles over edges, union + intersection − both | 8e-14 |
| Two hemispheres: a circle over the equator against the cap and 2 pi sin r | 8e-12 |
| Two hemispheres: a great circle over the equator twice, back at its start | 2e-12 |

In f32 the same checks hold to 3e-5. Figure 23 shows a fan of straightest
geodesics over a cube's edges, circles over edges and their union, and
circles across two hemispheres.

## Limits

- Edges are whole sides of the faces' domains. Faces trimmed along other
  curves, and edges joining parts of sides, are not covered.
- A circle round a vertex of a cube (a cone point), or reaching one, is not
  supported: its rays stop at the vertex.
