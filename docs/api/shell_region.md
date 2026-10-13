# Curves and regions across shell faces

`qaws_shell_region.h` carries the surface toolkit (`qaws_surface_region.h`)
over the shared edges of a shell: a cube of patches, a sphere of
hemispheres, a network of NURBS faces.

```c
qaws_shell_atlas* atlas;
qaws_shell_atlas_create(surfaces, n, 0, &atlas);        /* or _create_from_brep, _create_trimmed */

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
| Cube: a geodesic disk round a corner against the unfolded cone of angle 3 pi/2 | 5e-11 |
| Saddle (five squares round a vertex, angle 5 pi/2): a disk through it against pi r² + (pi/4)(r − rho)² | 3e-11 |
| Monkey saddle (six squares, angle 3 pi): a disk through the vertex against pi r² + (pi/2)(r − rho)² | 5e-11 |
| Monkey saddle of hypars (curved faces): the disk at two tolerances; warped either way | 1.4e-7; 0 |
| Monkey saddle of bicubic patches (curved faces and edges): two disks through the vertex, union + intersection − both; against the f64 disks in f32 | 1e-9; 1.3e-5 |

In f32 the same checks hold to 3e-5. Figure 23 shows a fan of straightest
geodesics over a cube's edges, circles over edges and their union, and
circles across two hemispheres, a disk round a cube corner and a disk
through a saddle vertex. Figure 24 shows disks through saddle vertices: five
and six faces round a vertex, flat and bent into bicubic patches, rays
parting at the vertex, and on the bent monkey saddle two disks through the
vertex under union, intersection, difference and exclusive-or, the surfaces
shaded and drawn back to front.

## Trimmed faces

`qaws_shell_atlas_create_trimmed(faces, n, 0, &atlas)` takes each face as a
surface and, optionally, loops of (u, v) curves bounding it
(`qaws_shell_face_desc`). Each curve of a loop is one side of the face.

- **Orientation.** Loops are turned to keep the face on their left: a loop
  inside an even number of others runs counter-clockwise, its holes
  clockwise. A loop that does not close is refused (`INVALID_ARGUMENT`); a
  trim on a surface closed both ways (a torus) is not covered
  (`UNSUPPORTED_OPERATION`).
- **Closed surfaces.** On a surface closed in u (or v) a loop is given in
  continuous (u, v), and may cross the seam (a hole on it) or wrap round
  once (a cylinder's rim: a line from (0, v) to (1, v)). It closes up to
  whole periods. Inside is then by the crossings of a ray across the closed
  direction over the copies of the point, and each loop is turned so that
  the point just left of it is inside. Exits and region walks follow the
  copy of the trim the geodesic meets.
- **Edges.** A trim curve meets any other face's side (a trim curve or a
  side of an untrimmed domain) where every point of each lies within the
  tolerance of the other, whatever the two parametrisations. A plate's hole
  meets the top of a tube, a square's circular hole meets a disk of another
  surface.
- **Inside.** A point is inside a trimmed face by its loops' polygons (128
  points per curve) away from the curves; within the polygons' sag of a
  curve, by which side of the curve itself it lies on.
- **Crossing.** Where a geodesic leaves a trimmed face, its exit is solved
  with the trim curve by Newton on (length, curve parameter); the crossing
  then takes the curve's own tangent and outward conormal.
- **Corners.** A join of two trim curves is a vertex of the shell where the
  boundary turns or the neighbour changes; a smooth join with one
  neighbour (a circle given as one curve) is not.
- **T-junctions.** A side running through another side's end (one long
  side meeting two faces side by side) is split there, so every edge joins
  whole sides: a trim curve into stretches of itself, and an open untrimmed
  face's side by turning the face into one trimmed by its domain. The T
  point is a vertex (its sides' neighbours differ). This holds for
  `qaws_shell_atlas_create` as well.
- **Faces without runs.** A face no boundary crosses is inside or outside as
  a neighbour's region is just inside their shared side, at five places
  along it (a region only touching the side does not count).
- **Regions.** A region's runs on a trimmed face are closed along the trim
  loops (stretches of the trim curves between them); a hole no run touches
  is taken out where the region holds it.

| Check (test 108, f64) | Agreement |
|---|---|
| Plane split by a circle (two trimmed faces of two surfaces): a ray over the disk and out against the straight line | 1e-15 |
| The same: a disk across the circle against pi r² and 2 pi r; each face's part against the lens | 3e-12; 6e-12; 2e-12 |
| Plate with a hole and a tube: three straightest rays into the tube against the helices of the unrolled cylinder | 7e-15 |
| Cube of planes trimmed to its faces: circles over one and two edges against pi r² and 2 pi r | 5e-11 |
| The same: a disk round a corner against the untrimmed cube's | 2e-12 |
| T-junction (a 2 × 1 face meeting two unit faces, flat and as walls): three disks against pi r² and 2 pi r | 1e-10 |
| The same: a disk on a wall just touching the edge below it stays on its wall | 4e-11 |
| Cylinder trimmed by its rims and a hole on the seam, the hole's disk of the same surface: three disks against pi r² and 2 pi r (unrolled) | 2e-9 |
| The same: a ray over the hole's disk and the seam against the unrolled line | 2e-15 |

In f32 the same checks hold to 3e-5. Figure 25 shows the split plane, the
plate and tube (helices through the hole, a disk over the rim, a union with
a disk inside the tube) and the trimmed cube.

## Limits

- A T-junction along a side of a face closed in u or v (a side that cannot
  be split) is not covered.
- Trims on surfaces closed both ways are not covered.
- Rays that hit a vertex exactly stop there; disks round vertices go
  through `qaws_shell_region_create_geodesic_disk`.

## Disks round vertices

`qaws_shell_region_create_geodesic_disk(atlas, centre, r, 0, &region)` gives
the points within r along the shell, vertices included.

- **Vertex directions.** Rays passing either side of a shell vertex part, so
  the circle's end point jumps there. Those jumps are found among 256 rays
  and pinned down by bisection on the angle.
- **Sectors.** The disk is the union of sectors between the vertex
  directions, stopping 1e-10 rad either side of each (1e-4 in float) and
  none wider than a quarter turn. Each sector is a ray out, the circle's arc
  and a ray back.
- **Behind a convex vertex.** Behind a cube's corner (a cone of angle
  3 pi/2), the sectors on the two sides overlap where the disk is reached
  both ways round, and the union is the disk.
- **Behind a saddle vertex.** Round a vertex with more than 2 pi of angle
  (five squares round it, 5 pi/2), the rays either side of the vertex
  direction part and leave a gap behind it. Every point within r − rho of
  the vertex (rho the vertex's distance from the centre) is within r of the
  centre, so the disk about the vertex of radius r − rho is added: the union
  of the corner sectors of each face at the vertex. Where geodesics through
  the vertex are shortest this fills the gap exactly; behind a convex vertex
  it adds nothing.
- **Runs by a corner.** A ray passing a vertex by a hair may leave a run a
  rounding long (left out) or, from its rounded start, overshoot the domain
  by a rounding at its end (shortened by as much).
- **Curved edges.** Faces bent into bicubic patches have curved edges, and
  rays from a vertex along them fold. A ray grazing a side, past it by a
  rounding while heading along it, crosses at the side it is nearest. A
  ray ending exactly on a shared side counts for either face. Arc fits stop
  refining when doubling no longer gains, since a fold leaves a kink no
  sampling removes. In float the vertex sectors stop 1e-4 rad inside the
  sides: a ray hugging a side meets the boundary walked along it nearly
  tangent, and the float clipper can lose a piece there.
- **Booleans checked.** Each face's Boolean must keep its bounds (a union no
  smaller than either part, an intersection no larger, a difference no
  larger than the first, an exclusive-or between the difference and the
  sum). Out of them, the face's Boolean is redone with exact topology
  (`qaws_surface_clip_exact`).
- **Chains.** A sector's runs on the centre's face meet inside the face, so
  regions take each loop's consecutive runs on one face as one chain (one
  composite curve).
