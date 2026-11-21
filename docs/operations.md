# Curve Operations, Conversion, Import, Export and Fitting

This page explains what the curve operations produce, how accurate each one is, and who owns the results. For the
signatures, see the API pages: [operations](api/operations.md), [convert](api/convert.md),
[import](api/import.md), [export](api/export.md) and [interpolating](api/interpolating.md).

> Names are given without their precision suffix. Each one exists as `name_f32` and `name_f64`; see
> [precision.md](precision.md).

## Three kinds of result

Every operation falls into one of three groups. The group tells you how far the result can be from the true geometry.

| Group | Result | Examples |
|---|---|---|
| Exact | lies on the source to rounding | split, extract, join, reverse, family conversions, degree elevation, transform under a similarity |
| Within a stated tolerance | the error bound is a parameter or a documented constant | `qaws_curve_offset_cubic_2d`, `qaws_curve_merge_chain`, extract of a Yuksel curve |
| Sampled | the accuracy depends on a sample count | `qaws_curve_offset_2d`, `qaws_curve_offset_3d`, SVG export of non-Bezier families, polyline and OBJ export |

The error of a sampled result is not measured. Raise the sample count, or use an operation from one of the first two
groups when the error must be bounded.

## Exact operations

### Split and extract

`qaws_curve_extract(t0, t1)` returns the piece of a curve between two parameters as a new curve that lies exactly on
the source. With `t0 > t1` the piece is reversed. The piece has the simplest kind that represents it exactly:

- Bezier, rational Bezier, B-spline, NURBS, polynomial, arc and clothoid pieces keep their kind. B-spline and NURBS
  pieces keep the source's parameters.
- Hermite, Catmull-Rom, trajectory and subdivision pieces become one cubic Bezier (a single span) or a cubic B-spline
  whose knots are the source's parameters.
- A composite gives the piece of one segment, or a composite of pieces.
- A reparameterized curve gives the piece of its source.
- Yuksel curves have no polynomial form. Their pieces are cubic B-splines fitted to 1e-9 of the piece's extent (1e-5
  in float builds).

Parameters within rounding of the domain are clamped onto it.

`qaws_curve_split` cuts at one interior parameter and returns two curves. Bezier, B-spline, NURBS and trajectory
curves are cut in their own family (De Casteljau, knot insertion). A Hermite curve cut at an integer parameter stays a
Hermite curve. Every other cut goes through extract. A parameter outside the domain returns `OUT_OF_RANGE`, and a
parameter on either end returns `INVALID_ARGUMENT`.

`qaws_curve_join` joins two curves of the same kind and dimension. It supports Hermite, B-spline, NURBS, Catmull-Rom
and trajectory curves and returns `UNSUPPORTED_OPERATION` for the other kinds. To chain curves of different kinds, use
a composite curve or `qaws_curve_merge_chain`.

### Family conversions and degree

`qaws_curve_convert_*` re-express a curve in another family without changing its shape:

- one Hermite or Catmull-Rom span as a Bezier curve;
- a Bezier curve as a B-spline;
- a B-spline as a NURBS curve.

`qaws_curve_elevate_degree` is exact. `qaws_curve_reduce_degree` (Bezier only, degree n to n - 1) is an approximation
whenever the curve is not of the lower degree already.

### 2D transforms and constructors

`qaws_curve_transform_2d` (in [path](api/path.md)) applies an affine map. Bezier, rational Bezier, B-spline, NURBS
and polynomial curves keep their kind exactly. Arcs and clothoids keep their kind under a similarity. Under other maps
an arc becomes an exact NURBS, and a clothoid becomes a cubic B-spline fitted to 1e-9 of its extent. The other kinds
are first converted to their exact extract form.

`qaws_curve_create_ellipse_2d` builds an exact ellipse as a rational quadratic NURBS made of four quarter arcs.

## Operations within a tolerance

`qaws_curve_offset_cubic_2d` returns the parallel curve as a composite of cubic Bezier segments. The `tolerance`
argument is the largest distance error, in the curve's units. The steps are:

1. Each span is made cubic.
2. Each cubic is split at its inflections and curvature maxima.
3. Each piece is tried as a single offset cubic, then as circular arcs, and is halved when neither fits.

A cusp becomes an arc about it, and a corner between spans becomes a line. Where the curvature radius is below
`|distance|`, the loops are kept. A curve with no extent returns `DEGENERATE_CURVE`.

`qaws_curve_merge_chain` fits a chain of C0-joined segments with fewer B-splines of `target_degree` (3 by default).
The result deviates from the chain by at most `tolerance`.

## Sampled operations

- `qaws_curve_offset_2d` samples the offset at 512 parameters. It removes the backward loops at cusps, and at higher
  `trim` levels it also removes points closer than `|distance|` and splits the remaining self-intersections into
  separate closed loops. The results are Catmull-Rom curves through the kept samples. Use this function when you need
  clean, non-self-intersecting boundaries. Use the cubic offset when you need an error bound.
- `qaws_curve_offset_3d` offsets along a constant direction, the Frenet normal or a surface normal. The result is a
  curve through `sample_count` offset samples (256 by default).
- `qaws_curve_fillet_2d` and `qaws_curve_chamfer_2d` insert circular arcs or straight lines at the sharp corners of a
  2D composite curve and return a new composite.

## Arc-length wrappers

`qaws_curve_reparameterize_arc_length` and `qaws_curve_match_arc_length` return wrapper curves of kind
`REPARAMETERIZED`. A wrapper maps its parameter to the source through a lookup table of `table_resolution` entries
(256 by default). The wrapper does not own its source: the source must outlive it.

## Fitting and interpolation

- `qaws_curve_fit_bspline` fits a least-squares B-spline with `control_point_count` control points to `m` data
  points. The end points are interpolated. The data parameters are either given, or derived from the points by
  chord length (`DEFAULT`, `CHORDAL`), equal steps (`UNIFORM`) or the square roots of the chord lengths
  (`CENTRIPETAL`).
- `qaws_curve_create_interpolating_spline` builds a smooth 2D spline through points by a global solve, so moving one
  point moves the whole curve a little. The result is a composite of Bezier segments on `[0, segment_count]`.

| Method | Segments | Continuity | Point i at parameter i |
|---|---|---|---|
| `CUBIC_C2` | cubic | C2; an open curve has zero second derivative at its ends | yes |
| `HOBBY` | cubic | G1, close to circular arcs | yes |
| `KAPPA` | quadratic, one per interior point | G2, or G1 at an inflecting join | no: each point lies inside its segment, at that segment's curvature maximum |
| `CUBIC_G2` | cubic | G2 | yes |
| `S9` | cubics approximating arc-length-curvature segments | G2 to the approximation's accuracy | no |

Consecutive duplicate points are dropped first. An open curve of two points gives the straight segment. A closed
curve needs three distinct points, and `KAPPA` needs three points when open as well. An `S9` solve uses about
145 KB of stack.

- `qaws_polyline_import_catmull_rom` and `qaws_polyline_import_trajectory` build interpolating curves from raw point
  arrays. The Catmull-Rom spacing is centripetal by default, and the trajectory has uniform timing.

## Export

- `qaws_curve_export_svg_path` writes SVG path data for a 2D curve. Bezier curves map to `M`/`C`/`Q`/`L` commands
  exactly, and Hermite and Catmull-Rom spans are converted to exact cubics. All other families are sampled
  (`sample_count`, 64 by default) and fitted with cubics. Numbers are written with 6 significant digits. If
  `*out_length >= capacity`, the string was truncated.
- `qaws_polyline_export_2d` and `qaws_polyline_export_3d` sample positions uniformly in the parameter, or weighted by
  curvature. Curvature sampling may write fewer than `sample_count` points.
- `qaws_surface_export_obj` writes a `u_samples x v_samples` grid of a surface as Wavefront OBJ triangles, with
  normals optionally included.

## Ownership of results

- Every `qaws_curve**` output is a new curve owned by the caller and released with `qaws_curve_destroy`.
  `qaws_curve_offset_2d` and `qaws_curve_merge_chain` write up to `curve_capacity` curve pointers into the caller's
  array, and each of those curves is owned by the caller.
- Curves derived from one source (split, extract, reverse, join, conversions, fillet, chamfer) are allocated with the
  source curve's allocator. All other results come from the allocator chain described in [memory.md](memory.md).
- Wrapper curves (arc-length reparameterization) borrow their source.
- Export writes into caller buffers and allocates no result.
- On failure the outputs are `NULL` or zero; see [error_handling.md](error_handling.md).

All of these functions read their inputs only. They are safe to call from several threads on shared curves (see
[threading.md](threading.md)).
