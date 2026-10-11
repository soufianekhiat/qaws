# API Reference: qaws_skeleton.h

The straight skeleton of a region: the boundary moves inwards, each edge
parallel to itself at its own speed (its weight), and the skeleton is the
trace of the moving boundary's vertices. Each input edge sweeps one face, and
the time a point is swept is the roof's height there.

```c
qaws_skeleton_desc d = { paths, path_count, QAWS_FILL_NON_ZERO };
d.weight_fn = speed; d.weight_user = user;               /* optional: weighted skeleton */
qaws_skeleton* s;
qaws_skeleton_create(&d, &s);

qaws_skeleton_get_node(s, i, &node);       /* position and time */
qaws_skeleton_get_arc(s, i, &arc);         /* node_a -> node_b, between face_a and face_b */
qaws_skeleton_get_face(s, i, &face);       /* the input edge (nodes 0, 1), then the skeleton back */
qaws_skeleton_get_max_time(s);             /* the roof's top */
qaws_skeleton_offset(s, t, &result);       /* the wavefront at time t, as a qaws_clip result */
qaws_skeleton_destroy(s);
```

## How it works

1. **Region.** The paths are united under the fill rule (`qaws_clip`), so
   outer boundaries run counter-clockwise and holes clockwise, with the
   interior on the left. The curves must be straight (degree 1 Béziers or
   B-splines, polylines).
2. **Planes.** Each edge on the line n · p = d (n the inward unit normal)
   lies at time t on n · p = d + w t: a plane over (x, y, t).
   - A wavefront vertex is two planes meeting. It moves at a constant
     velocity, the solution of n_a · v = w_a, n_b · v = w_b.
   - Every event is three planes meeting at a point, so one 3 x 3 solve
     covers edge, split and vertex events.
3. **Events.** Each new vertex queues the times it meets its two
   neighbours. A reflex vertex also queues the time it reaches every
   plane it is closing on. A heap orders them by time.
4. **The generalised intersection event** (Kelly and Wonka, *Interactive
   architectural modeling with procedural extrusions*, TOG 2011; their
   implementation is twak/campskeleton). At a popped point and time,
   every live vertex there and every wavefront edge through it is taken at
   once:
   - Runs of consecutive vertices there are chains. An edge through the
     point is a chain of its own, which it splits.
   - The edges inside a chain end at the point. A loop entirely there
     closes.
   - The chains, sorted by the angle of their first edge round the point,
     are joined: each chain's last edge to the next chain's first, by a new
     vertex.

   Edge events, split events, vertex events and any number of them
   coinciding are all this one case.
5. **Antiparallel edges.** A new vertex between antiparallel edges sits
   where they lie on one line, as when a corridor of constant width closes
   all at once. It runs along that line to the nearer end, or to both
   when they are as near, and the skeleton gets that stretch. If other
   vertices are already where it lands, they make one event with it.
6. **Faces.** Each plane's arcs are walked from its edge's end back to
   its start.
7. **Offsets.** `qaws_skeleton_offset` reruns the events up to t. It
   returns the loops of the vertices' positions at t, united with the
   positive fill rule.

## Precision

- **Tolerance.** Events within the tolerance (default 1e-9 of the
  extent, 1e-6 in float builds for the input's rounding) are one.
  Arithmetic is in double precision throughout.
- **Weights.** Weights must be > 0. Consecutive collinear edges with
  different weights are ambiguous (Biedl, Held, Huber, Kaaser and
  Palfrader, CGTA 2015), and their vertex moves at the average speed.
- **Cost.** Each reflex vertex is checked against every live plane, and
  each event scans the live vertices, so the cost is quadratic. A random
  star-shaped polygon takes 0.5 ms with 64 vertices, 80 ms with 1,024 and
  2 s with 4,096.

## Tests

Test 99 (`tests/99_skeleton.c`) checks every shape three ways:

- **Distance.** Both ends of every arc lie at weighted distance w t from
  both faces' edges.
- **Tiling.** The faces' areas sum to the region's (the paths united by
  `qaws_clip`).
- **Wavefront.** At five times, the wavefront's area (a separate run of
  the events) equals the faces' unswept parts: each face clipped by its
  edge's line moved to w t.

Closed forms:

- the 2 x 2 square: 5 nodes, 4 arcs, top at 1;
- the 3 x 1 rectangle: a ridge of 2 at time 0.5;
- regular 3-, 10-, 17- and 24-gons: one centre node at the inradius;
- a 4 x 4 square with a centred 2 x 2 hole: every event at time 0.5;
- a 3 x 1 rectangle with its short sides at speed 2: a ridge of 1.

Shapes:

- an L, a comb, a plus (every event degenerate), a floor plan, an
  orthogonal spiral of equal-width corridors;
- a 7-point star and random star-shaped polygons, with and without holes;
- twelve unions of random rectangles on an integer grid.

On a further 1,000 random grid unions, outside the suite, every face set
tiled its region.

It runs in f64 and f32.
