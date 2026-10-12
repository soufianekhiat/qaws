# API Reference: qaws_mean_value.h

Mean value coordinates over regions bounded by curves of any kind: concave,
with holes. They are generalized barycentric coordinates whose cage is the
boundary itself, with the boundary curves' control points as handles.

```c
qaws_mean_value_desc d = { paths, path_count };      /* closed loops; holes by nesting */
d.executor = &executor;                               /* optional: deform in parallel */
qaws_mean_value* m;
qaws_mean_value_create(&d, &m);

unsigned int n = qaws_mean_value_get_handle_count(m);
qaws_mean_value_get_handle(m, i, &p);                 /* the handles at rest */
qaws_mean_value_coordinates(m, x, weights);           /* n weights: sum 1, sum w_i P_i = x */
qaws_mean_value_deform(m, moved, points, count, out); /* x -> sum w_i(x) P_i' */
qaws_mean_value_interpolate(m, x, f, user, dim, out); /* the interpolant of data on the boundary */
qaws_mean_value_destroy(m);
```

## What they are

The transfinite mean value interpolant (Ju, Schaefer and Warren 2005;
Dyken and Floater 2009) of a function f on the boundary is

```
g(x) = ∮ f(c) k ds / ∮ k ds,     k = (c - x) · n / |c - x|^3,
```

over every boundary loop, n the outward normal. Seen from x, it averages f
over the directions, each boundary point weighted by 1 / distance. A ray
from x that crosses the boundary several times meets it with alternating
signs, which `(c - x) · n` gives, so concave regions and holes need nothing
more. g interpolates f on the boundary and reproduces linear functions
exactly.

Each boundary curve is its control points times its basis functions:
- Bernstein for Beziers;
- B-spline for B-splines;
- rational for weighted kinds.

So a handle's coordinate is the interpolant of its basis function:

```
λ_i(x) = ∮ B_i k ds / ∮ k ds.
```

- They sum to 1, and `Σ λ_i P_i = x`.
- On a polygon they are Floater's mean value coordinates (2003), and
  Hormann and Floater's for arbitrary polygons (2006).
- Moving the handles moves the boundary curves. `x -> Σ λ_i(x) P_i'`
  carries the inside along: a cage deformation whose cage is the region's
  own curves.
- On the boundary the coordinates are the curves' basis functions, so
  boundary points follow the moved curves exactly.

## Handles

Handles are listed in loop order, then curve order. A point shared by
consecutive curves, or by a loop's two ends, is one handle (within 1e-12 of
the extent, 1e-6 in float builds).

| Curve | Handles |
|---|---|
| Lines, polylines (degree 1) | the vertices |
| Bezier, rational Bezier | the control points |
| B-spline, NURBS | the control points (each span by Bezier extraction) |
| Arc | the control points of rational quadratic pieces of at most 90 degrees |
| Others (Hermite, Catmull-Rom, clothoid, ...) | those of a cubic B-spline fitted to the curve (exact for cubic Hermite and Catmull-Rom) |

## How it is computed

1. **Pieces.** Every curve becomes rational Bezier pieces. Their
   homogeneous control points are combinations of the handles: span by
   span, a B-spline's basis in Bernstein form (`a_ji`), times the weights.
   On a piece of degree d with Bernstein weights `w_j`, handle i's basis is
   `Σ_j a_ji B_j / W`, `W = Σ_j w_j B_j`. So the d + 1 integrals
   `J_j = ∫ B_j k / W du` give every handle's coordinate on the piece.
2. **Lines** (degree 1, equal weights) have the integrals in closed form:
   `J_0 = tan(α/2) / (r_0 w)` and `J_1 = tan(α/2) / (r_1 w)`, α the signed
   angle the segment subtends at x.
3. **Curved pieces** are integrated by adaptive Gauss-Kronrod (7-15).
   - The kernel peaks at the piece's point nearest x, with a width of the
     distance. The intervals start there and double outwards from that
     width.
   - Near the boundary, `c - x` is good only to about 1e-16 of the extent
     over the distance, and so is the kernel. The tolerance is no finer
     than that.
4. **Loops** are turned so that outside loops run counter-clockwise and
   holes clockwise. A loop's nesting depth (the winding of the others round
   it) decides which it is, whatever its orientation.
5. **On the boundary** (within 1e-13 of the extent, 1e-7 in float builds)
   the coordinates are the basis functions there, and the interpolant is f
   there.

Outside the region the formula still holds, but its denominator can vanish.
Where it does, `NUMERICAL_FAILURE`.

## Accuracy and cost

- **Tolerance.** `tolerance` is the relative accuracy of the integrals. It
  defaults to 1e-10 (1e-5 in float builds).
- **Test 100**, in f64:
  - against Mathematica: the twelve coordinates of a region of four cubic
    Beziers at four points (one near a corner, one near an edge), and the
    interpolants over a NURBS ellipse and over a disc with a square hole,
    agree to 1e-15;
  - on polygons, against Floater's closed form: 3e-16;
  - over every kind, the reproduction `Σ λ_i P_i = x`: 1e-15 of the size;
  - towards the boundary, at 1e-8 the interpolant is within 4e-8 of f there.
- **Cost.** One point takes about 2 µs on a polygon of 100 vertices,
  19 µs on 1,000, and 34 µs on a region of four cubic Beziers. That is
  linear in the boundary's pieces. `deform` runs its points on the
  executor.

## Tests

Test 100 (`tests/100_mean_value.c`), f64 and f32, with ground truth from
`tests/reference/100_mean_value.wls`:
- **Mathematica:** the blob's coordinates, the ellipse's and the disc's
  interpolants.
- **Polygons:** an L with a hole against Floater's closed form. The same
  polygon with straight quadratic Bezier edges, integrated rather than in
  closed form, gives the same coordinates. The interpolant of piecewise
  linear boundary data equals the coordinates' combination. A convex
  hexagon's coordinates are positive.
- **Every kind** (Beziers, a NURBS ellipse with a round hole, a lens of
  rational Beziers, a closed cubic B-spline, an arc with a Hermite cubic):
  the coordinates sum to 1 and reproduce the point, at Sobol points over
  each region.
- **Boundary:** the coordinates are the curve's Bernstein values on it,
  and the interpolant tends to f from inside.
- **Orientation:** either orientation of either loop changes nothing.
- **Deformation:** handles unmoved, translated or under an affine map move
  the points the same way, and an executor gives the serial points.
