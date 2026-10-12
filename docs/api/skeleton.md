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
d.weight_fn = speed; d.weight_user = user;   /* optional: weighted (lines) */
d.splice_angle = 3.14159265358979;           /* optional: the medial axis */
qaws_skeleton* s;
qaws_skeleton_create(&d, &s);

qaws_skeleton_get_node(s, i, &node);         /* position and time */
qaws_skeleton_get_arc(s, i, &arc);           /* node_a -> node_b, faces, smooth, boundary */
qaws_skeleton_evaluate_arc(s, i, u, &p, &t); /* exact point and time along it */
qaws_skeleton_get_arc_curve(s, i, &curve);   /* a line, a circular arc, or the exact conic (NURBS) */
qaws_skeleton_get_face(s, i, &face);         /* the edge, then the skeleton back; its support */
qaws_skeleton_get_max_time(s);               /* the roof's top */
qaws_skeleton_offset(s, t, &result);         /* the wavefront at time t: lines and arcs */
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
   - an arc |p − m| = r is on |p − m| = r − s t (s = +1 when the region is
     inside its circle, so it shrinks; −1 when outside, so it grows): a
     cone.

   Every event is three supports meeting at a point. The difference of two
   cones is a plane, so each event is a linear system and at most one
   quadratic.
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
8. **Faces.** Each support's arcs are walked from its edge's end back to its
   start, turning so the face stays on the left. A face is entered along
   a smooth path (a slit) and left back along it.
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
- **Weights.** Weights apply to lines and must be > 0; arcs move at speed 1.
  Consecutive collinear lines with different weights are ambiguous
  (Biedl, Held, Huber, Kaaser and Palfrader, CGTA 2015), and their vertex
  moves at the average speed.
- **Cost.** Each splitting corner is checked against every live support,
  and each event scans the live corners, so the cost is quadratic. A
  random star-shaped polygon takes 0.5 ms with 64 vertices, 90 ms with
  1,024 and 2.3 s with 4,096.
- **Robustness.**
  - Polygons are robust: 1,000 random unions of grid rectangles, whose
    events are almost all degenerate, all tile their regions.
  - Curved regions pass every structured case in the tests. Random
    star-shaped regions whose edges are random bulging arcs, with or
    without a round hole, were run 400 times at each sigma:
    - 789 of the 800 runs complete;
    - each completed run tiles its region;
    - its arcs lie on their supports to within 3e-6 of the size: nodes at
      events near tangencies are good to half the digits;
    - the other 11 stop with QAWS_STATUS_NUMERICAL_FAILURE.
  - **Checks and retries.** A run fails when a loop never closes, when
    events keep remaking one another, or when the faces do not close up
    (a face's walk does not come back to its edge, or an arc between two
    faces is not on both). With the default tolerance, a failed run is
    redone with the tolerance 8 times finer, 8 times coarser, then 64
    times finer, since the near-coincidences that break a run depend on
    where the tolerance cuts. In float builds only the coarser one is
    tried: below the vertices' rounding, their meetings come apart. An
    explicit `tolerance` is run once.
  - In float builds the vertices are rounded off the curves' meetings,
    and about six in ten such random regions fail.

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
- the medial axes of an L, a comb and a star.

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
