# API Reference: qaws_skeleton.h

Straight and curved skeletons of a region. The boundary moves inwards, each
edge at its own speed (its weight), and the skeleton is the trace of the
moving boundary's vertices. Each edge sweeps one face, and the time a point
is swept is the roof's height there.

- **Straight lines** move parallel to themselves: the straight skeleton.
- **Circular arcs** shrink or grow about their centres, and the vertices
  between them travel on conics: the mitered skeleton of circular arc
  polygons.
- **The splice angle sigma** controls reflex vertices. sigma = 2 pi (the
  default) keeps them mitered; sigma = pi rounds every one at once, which
  gives the medial axis and the classical offsets, for polygons as well as
  for arcs.

```c
qaws_skeleton_desc d = { paths, path_count, QAWS_FILL_NON_ZERO };
d.weight_fn = speed; d.weight_user = user;   /* optional: weighted (lines and arcs) */
d.splice_angle = 3.14159265358979;           /* optional: the medial axis */
qaws_skeleton* s;
qaws_skeleton_create(&s, &d);

qaws_skeleton_get_node(&node, s, i);         /* position and time */
qaws_skeleton_get_arc(&arc, s, i);           /* node_a -> node_b, faces, smooth, boundary */
qaws_skeleton_evaluate_arc(&p, &t, s, i, u); /* exact point and time along it */
qaws_skeleton_get_arc_curve(&curve, s, i);   /* a line, a circular arc, or the exact conic (NURBS) */
qaws_skeleton_get_face(&face, s, i);         /* the edge, then the skeleton back; its support */
qaws_skeleton_get_max_time(s);               /* the roof's top */
qaws_skeleton_offset(&result, s, t);         /* the wavefront at time t: lines and arcs */
qaws_skeleton_destroy(s);
```

## How it works

1. **Region.** The paths are united under the fill rule (`qaws_clip`), so
   outer boundaries run counter-clockwise and holes clockwise, with the
   interior on the left. Lines (degree 1 Béziers or B-splines, polylines)
   and circular arcs (`qaws_curve_create_arc`) are taken exactly. Any other
   curve (Béziers, NURBS such as ellipses, B-splines, ...) becomes a G1 arc
   spline:
   - each span is fitted with an equal-tangent biarc, which is halved until
     its distance to the curve is within `curve_tolerance` (default 1e-4 of
     the extent);
   - consecutive arcs on one circle merge;
   - the joints are smooth and carry no skeleton.

   The skeleton is exact for the spline. Where curvature peaks (a petal's
   tip), the spline's piecewise-constant curvature ends the branch in a
   short segment between the centres of adjacent arcs, instead of at one
   centre of curvature.
2. **Supports in space-time.** Each edge's carrier is a surface over
   (x, y, t):
   - a line n · p = d at time t is on n · p = d + w t: a plane;
   - an arc |p − m| = r is on |p − m| = r − s w t (s = +1 when the region
     is inside its circle, so it shrinks; −1 when outside, so it grows): a
     cone.

   Every event is three supports meeting at a point. The difference of two
   cones of one speed is a plane, so each event is a linear system and at
   most one quadratic. Cones of different speeds differ by a quadric: x
   and y are then quadratic in t and the event is a quartic, its real
   roots by bisection between its derivative's.
3. **Corners.** A corner between two lines, or between two supports that
   stay tangent (a smooth joint), moves at a constant velocity. Any other
   corner is one of the two crossings of its supports (its branch), moving
   on a conic.
4. **Events**, in a heap by time:
   - a corner meets a support (its neighbour's, or any support when it can
     split an edge);
   - two supports touch (squeeze);
   - an arc's two ends meet round its circle (bubble), and it becomes a loop
     of its own;
   - a circle shrinks to its centre;
   - a reflex corner's angle reaches sigma, or its supports turn tangent: it
     splices, opening a circle that grows from its point (Weiß, Jüttler and
     Aurenhammer, *Mitered offsets and straight skeletons for circular arc
     polygons*, EuroCG 2018 / IJCGA 2020).
5. **The generalised intersection event** (Kelly and Wonka, *Interactive
   architectural modeling with procedural extrusions*, TOG 2011; their
   implementation is twak/campskeleton). At a popped point and time,
   every live corner there and every edge through it is taken at once:
   - Runs of consecutive corners there are chains. An edge through the
     point is a chain of its own, which it splits.
   - The edges inside a chain end at the point. A loop entirely there
     closes.
   - The chains, sorted by the angle of their first edge round the point,
     are joined: each chain's last edge to the next chain's first, by a new
     corner.

   Edge, split and vertex events, and any number of them coinciding, are
   all this one case.
6. **Antiparallel edges.** Coincident opposite supports (lines, or
   concentric circles of one radius) leave a corner that runs along them
   to the nearer end. A corridor of constant width closing all at once
   therefore collapses along its middle.
7. **Near tangencies.** These measures handle them:
   - A corner with a circle among its supports counts as at an event's
     point when both its supports pass through the point. Its position
     there moves like a square root of the time, so position alone is
     ill-conditioned. Both crossings of its supports pass through the
     point, though, so it counts only when its own crossing is the nearer
     one, unless the supports touch then or are about to. Just after a
     touching the two crossings part, and each is where it is.
   - Events among circles within the tolerances of one another are one
     event, because their times are near-double roots and good to only
     half the digits.
   - A corner born at its supports' tangency is classified just after it.
   - A loop an event leaves enclosing nothing, which turns inside out at
     once, is a sliver between two supports just past their touching. It
     closes there, and its two faces meet along it.
   - A circle opened at a splice is checked against every corner that can
     split an edge, not only its neighbours.
   - A corner between two lines starts exactly on both. An event's point
     is good only to the tolerances (among circles, to the square root of
     one), and its error would ride along the corner's whole path.
   - A tangency where a corner's supports start to cross does not end it
     (a corner born just before, its event within the tolerances of the
     touching).
   - A corner next to a smooth joint between a line and a circle is also
     met from the joint's side. From its own side the meeting is at the
     joint's tangency, a double root.
   - Supports touch only where both have edges. The corners of a circle
     opened a moment ago are within the tolerances of each other, and are
     not an event.
   - Corners an event made a moment ago near its point are its result, and
     a later event there leaves them. An event that finds a corner but
     nothing for it to meet takes them in after all.
8. **Faces.** Each support's arcs are walked from its edge's end back to its
   start, turning so the face stays on the left. A face is entered along
   a smooth path (a slit) and left back along it. A circle opened and
   closed at one instant (in a sliver) has an empty face: its one node.
9. **Offsets.** `qaws_skeleton_offset` reruns the events up to t. It
   returns the loops of lines and arcs at t, united with the positive
   fill rule.

## Pure skeleton

Each arc carries two flags:
- `smooth`: it is a smooth joint's path. Such paths are part of the
  faces' boundaries, not of the skeleton proper.
- `boundary`: it has an end on the input (time 0).

The arcs with neither flag are the pure skeleton, clear of the boundary.

## Precision

- **Tolerance.** Events within the tolerance are one. It defaults to 1e-9
  of the extent. In float builds it is 1e-6 for the vertices' rounding,
  except when free-form curves are fitted: the fitted arcs are exact
  doubles, so their events keep the double default.
  Arithmetic is in double precision throughout.
- **Float input.** A float vertex misses a float arc by its rounding. In
  float builds each input edge is first made exact in double:
  - a line goes through its start and the next edge's;
  - an arc's centre moves onto the bisector of its two ends, so its circle
    passes through both.

  The meetings at the vertices are then exact. Near tangencies resolve on
  the double's scale, whatever the float tolerance.
- **Weights.** Weights apply to lines and arcs and must be > 0. A circle
  opened at a reflex vertex moves at the faster of its two edges' speeds,
  so at sigma = pi the result is the weighted medial axis: each edge's
  distance divided by its weight, a vertex's by the larger of its two.
  - The opened circle stays tangent to the edge of its own speed; a slower
    edge crosses it (a conic corner, the crossing behind the vertex along
    the edge coming in, ahead along the one going out).
  - A smooth joint whose edges part (a line next to a faster arc) leaves a
    gap: it opens a circle at once.
  - Supports of different speeds born touching take the crossing where the
    edge coming in passes from the unswept side of the one going out to
    its swept side.
  - Arcs between circles of different speeds are Cartesian ovals, not
    conics: `qaws_skeleton_get_arc_curve` fits them with rational
    quadratics (16 pieces at most).
  - `weight_fn` is called with the input's path and curve (an edge of the
    united region is matched back to the input curve nearest its middle)
    and the span: a Bézier's whole curve, a B-spline's knot span, an arc
    curve's segment.

  Consecutive collinear lines with different weights are ambiguous
  (Biedl, Held, Huber, Kaaser and Palfrader, CGTA 2015), and their vertex
  moves at the average speed.
- **Cost.**
  - A splitting corner finds its meetings with every live support ahead,
    keeps the 8 earliest and queues one at a time: the queue holds a few
    events per corner, not one per corner and support. A convex corner
    whose edges are lines or growing circles' arcs (ahead of their
    tangents) splits nothing first and looks for none.
  - Touchings are kept the same way: each support keeps its 8 earliest
    touchings with the others and queues one at a time. A support without
    edges is dropped from them, since it never has edges again.
  - An event finds the corners at its point, and the edges through it,
    through a grid. Each corner and edge is entered where it passes until
    it has moved about a cell, then entered again:
    - a straight path by its ends;
    - a conic by samples, with the largest step between them to spare
      (near its supports' tangency, where it speeds up, it is looked at
      by every event);
    - an arc by its boxes at the samples, with its ends' steps and the
      circle's growth to spare.
  - Collecting a splitting corner's meetings, and a new circle's
    touchings, still scans the supports, so those parts stay quadratic.
  - A random star-shaped polygon with spiky sides takes 4 ms with 256
    vertices, 25 ms with 1,024, 0.28 s with 4,096 and 1.0 s with 8,192
    (mitered). Its medial axis (sigma = pi) takes 7 ms, 59 ms, 0.58 s and
    2.7 s.
- **Robustness.**
  - Polygons are robust: 1,000 random unions of grid rectangles, whose
    events are almost all degenerate, all tile their regions. Their
    medial axis (sigma = pi) is weaker. Random star-shaped polygons,
    nearly round or with spiky sides, all build from 48 to 8,192
    vertices (685 tried).
  - Curved regions pass every structured case in the tests. Random
    star-shaped regions whose edges are random bulging arcs, with or
    without a round hole, were run 400 times at each sigma:
    - all 800 runs complete;
    - each completed run tiles its region;
    - its arcs lie on their supports to within 3e-6 of the size: nodes at
      events near tangencies are good to half the digits;
  - **Checks and retries.** A run fails when a loop never closes, when
    events keep remaking one another, or when the faces do not close up
    (a face's walk does not come back to its edge, or an arc between two
    faces is not on both). With the default tolerance, a failed run is
    redone with the tolerance 8 times finer, 8 times coarser, then 64
    times finer, since the near-coincidences that break a run depend on
    where the tolerance cuts. In float builds, whose input is made exact,
    it is redone 8 times coarser, then at the double default (1e-3 times
    finer) and 8 times finer still. An explicit `tolerance` is run once.
  - In float builds the same random regions (400 at each sigma) build in
    793 of 800 runs, and each run that builds tiles its region. Nodes are
    floats, so arcs lie on their supports to about 1e-6; in ten runs some
    are off by 1e-4 to 4e-3. The worst is a sliver closed up at one
    instant, its seam on the face of a circle opened and closed at once.

## Tests

Test 99 (`tests/99_skeleton.c`) runs in f64 and f32.

**Polygons.** Every shape is checked three ways:
- **Distance.** Both ends of every arc lie at weighted distance w t from
  both faces' edges.
- **Tiling.** The faces' areas sum to the region's.
- **Wavefront.** At five times, the wavefront's area (a separate run of
  the events) equals the faces' unswept parts.

Closed forms:
- the square and the 3 x 1 rectangle (its ridge);
- regular polygons: one centre node at the inradius;
- a square with a centred square hole: every event at one time;
- a weighted rectangle.

Shapes: an L, a comb, a plus, a floor plan, an orthogonal spiral, stars,
random star-shaped polygons with holes, and unions of grid rectangles.

**Curved.** These checks apply:
- **Supports.** Arcs sampled along their length (exactly, with
  `qaws_skeleton_evaluate_arc`) lie on both faces' supports at their
  times.
- **Medial axis.** At sigma = pi, every point off the smooth paths is at
  distance t from the input boundary, by qaws's closest-point queries.
- **Tiling.** The faces, sampled along their edges and arcs, tile the
  region.
- **Wavefront.** Where nothing reflex is kept mitered, the wavefront's area
  equals the round inward offset (`qaws_offset`, round joins). For a disc
  with an off-centre hole it is checked against the closed form instead.

Curved shapes:
- a disc, a stadium, a rounded rectangle (smooth joints), a half disc, a
  lens;
- a disc with an off-centre hole (squeeze, then an elliptic skeleton);
- a square with a round bite;
- a square with a round bump, mitered (it pinches off as a bubble) and at
  sigma = pi;
- the medial axes of an L, a comb and a star;
- weighted: a rectangle of four lines, its bottom at weight 2 (the top at
  2/3), a disc at weight 2 (the top at 1/2), the L's weighted medial axis
  (sides at 2 and 1.5: every point at the least of the edges' distances
  over their weights, a vertex's over the larger of its two), a stadium
  whose arcs move at 2 and 0.5 (the top at 1), the bump's weighted medial
  axis (its arc at 0.5: a three-circle event with collinear centres).

**Random curved regions** (double builds), at both sigma: 40 star-shaped
regions of 5 to 18 sides, each a line or an arc bulging either way, two
in five with a round hole.
- All 80 runs build.
- Every run that builds tiles its region, and its arcs lie on their
  supports to 1e-5 of the size.
- At sigma = pi, the wavefront equals the round offset.

**Free-form**, at both sigma: an ellipse (NURBS), an ellipse with a round
hole, and a blob of four cubics.
- Arcs lie exactly on their supports.
- The faces tile the region.
- At sigma = pi, every skeleton point is within twice the fit tolerance of
  its time, measured against the true curves, not the spline.
