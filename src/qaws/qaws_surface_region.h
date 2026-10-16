#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SURFACE_REGION_H_F32
#define QAWS_SURFACE_REGION_H_F32
#define QAWS_SURFACE_REGION_H_BODY
#endif
#else
#ifndef QAWS_SURFACE_REGION_H_F64
#define QAWS_SURFACE_REGION_H_F64
#define QAWS_SURFACE_REGION_H_BODY
#endif
#endif
#ifdef QAWS_SURFACE_REGION_H_BODY
#undef QAWS_SURFACE_REGION_H_BODY
#define QAWS_SURFACE_REGION_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_path.h"
#include "qaws_clip.h"

/*
 * Curves and regions on a parametric surface, defined by the surface's own
 * geometry (docs/surface_regions_plan.md).
 *
 * A surface S(u, v) carries the metric
 *
 *     ds^2 = E du^2 + 2 F du dv + G dv^2,
 *     E = S_u . S_u,  F = S_u . S_v,  G = S_v . S_v,
 *
 * and its geodesics (locally shortest paths) solve, by arc length,
 *
 *     u'' + G1_11 u'^2 + 2 G1_12 u' v' + G1_22 v'^2 = 0,
 *     v'' + G2_11 u'^2 + 2 G2_12 u' v' + G2_22 v'^2 = 0,
 *
 * the Christoffel symbols G^k_ij from S's first and second derivatives. A
 * curve on the surface is kept as its image in the parameter plane (a qaws
 * 2D curve with x = u, y = v), traced intrinsically: a geodesic between two
 * points is solved on the surface, a geodesic circle is the exponential map
 * of a tangent circle (every point at one geodesic distance from the
 * centre). Their (u, v) images are generally curved; they are not
 * projections of planar shapes.
 *
 * Closed directions (a cylinder's or a revolution's u, found by comparing
 * the surface's opposite edges) are followed continuously: a curve's (u, v)
 * image lives in the covering plane, u beyond the domain where it went
 * round. Poles (S_u x S_v = 0) are not crossed.
 *
 * Accuracy: tolerance is relative to the parameter domain's diagonal; a
 * curve's (u, v) image is within it of the true curve (0 = 1e-10, 1e-5 in
 * float builds). Arithmetic is in double precision throughout.
 */

typedef struct qaws_surface_curve qaws_surface_curve;

#ifndef QAWS_ENUM_QAWS_SURFACE_CURVE_KIND
#define QAWS_ENUM_QAWS_SURFACE_CURVE_KIND
typedef enum qaws_surface_curve_kind
{
	QAWS_SURFACE_CURVE_UV = 0,          /* a curve drawn in (u, v) */
	QAWS_SURFACE_CURVE_GEODESIC,        /* a geodesic segment (or ray) */
	QAWS_SURFACE_CURVE_GEODESIC_ARC,    /* an arc of a geodesic circle */
	QAWS_SURFACE_CURVE_GEODESIC_OFFSET  /* a curve's geodesic offset */
} qaws_surface_curve_kind;
#endif

/* The metric and Christoffel symbols at (u, v): g = { E, F, G },
   christoffel = { G1_11, G1_12, G1_22, G2_11, G2_12, G2_22 } (either may be
   NULL). NUMERICAL_FAILURE where the surface is singular (S_u x S_v = 0). */
qaws_status qaws_surface_get_metric(qaws_surface const* surface, qaws_scalar u, qaws_scalar v,
	qaws_scalar* out_g, qaws_scalar* out_christoffel);

/* The exponential map: the point reached from p along the geodesic leaving
   with the (u, v) velocity direction, after the surface length `length`;
   and the unit (u, v) velocity there (may be NULL). */
qaws_status qaws_surface_exp(qaws_surface const* surface, qaws_vec2 p, qaws_vec2 direction, qaws_scalar length,
	qaws_scalar tolerance, qaws_vec2* out_point, qaws_vec2* out_velocity);

/* A curve drawn in (u, v): the 2D curve is copied. */
qaws_status qaws_surface_curve_create_uv(qaws_surface const* surface, qaws_curve const* uv_curve,
	qaws_surface_curve** out_curve);

/* The geodesic from a to b (parameter points): shooting from a, by Newton
   on the leaving direction and the length, from the chord (and by
   continuation along it when that fails). Locally shortest; on surfaces
   with several geodesics between the points it is the one the chord leads
   to. */
qaws_status qaws_surface_curve_create_geodesic(qaws_surface const* surface, qaws_vec2 a, qaws_vec2 b,
	qaws_scalar tolerance, qaws_surface_curve** out_curve);

/* The geodesic leaving p in the (u, v) direction, for the surface length
   `length`. */
qaws_status qaws_surface_curve_create_geodesic_ray(qaws_surface const* surface, qaws_vec2 p, qaws_vec2 direction,
	qaws_scalar length, qaws_scalar tolerance, qaws_surface_curve** out_curve);

/* An arc of the geodesic circle of radius r about centre: the points
   exp_centre(r (cos a e1 + sin a e2)) for a from angle0 to angle1, e1 along
   S_u and e2 a quarter turn from it about S_u x S_v (angles in radians; a
   full turn closes the circle). Beyond the injectivity radius the circle
   crosses itself: INVALID_ARGUMENT. */
qaws_status qaws_surface_curve_create_geodesic_arc(qaws_surface const* surface, qaws_vec2 centre, qaws_scalar radius,
	qaws_scalar angle0, qaws_scalar angle1, qaws_scalar tolerance, qaws_surface_curve** out_curve);

/* The geodesic offset of a curve: each point carried the surface distance
   |distance| along the geodesic leaving the base square to it, to the left
   of travel (about S_u x S_v) for a positive distance, to the right for a
   negative one. The base must be smooth (G1); where it curves more than
   1/|distance| towards the offset, the offset has cusps and turns back. */
qaws_status qaws_surface_curve_create_offset(qaws_surface_curve const* base, qaws_scalar distance, qaws_scalar tolerance,
	qaws_surface_curve** out_curve);

/* A validated enclosure of the geodesic leaving p in the (u, v) direction
   for the surface length `length`, on a polynomial surface (Bezier,
   B-spline, bilinear, biquadratic): the length cut into steps, each with a
   box in (u, v) proven to hold the true geodesic over it, and a box proven
   to hold its end and unit velocity. Interval arithmetic with outward
   rounding throughout: the surface's derivatives up to the third bounded
   from its control nets (B-splines cut into Bezier patches), each step's
   box by Picard iteration, its end by a second-order Taylor form with the
   remainder over that box. width: the end box's side to aim for, relative
   to the domain's diagonal (0 = 1e-6; the steps shrink to reach it, about
   ten times the steps for a hundred times narrower).
   UNSUPPORTED_OPERATION for other surfaces; OUT_OF_RANGE when a box leaves
   the domain; NUMERICAL_FAILURE when a step cannot be validated (the
   metric degenerates). */
typedef struct qaws_geodesic_enclosure qaws_geodesic_enclosure;

qaws_status qaws_surface_geodesic_enclose(qaws_surface const* surface, qaws_vec2 p, qaws_vec2 direction, qaws_scalar length,
	qaws_scalar width, qaws_geodesic_enclosure** out_enclosure);
void qaws_geodesic_enclosure_destroy(qaws_geodesic_enclosure* enclosure);
unsigned int qaws_geodesic_enclosure_get_step_count(qaws_geodesic_enclosure const* enclosure);
/* Step i: its stretch of length [t0, t1] and the (u, v) box holding the
   geodesic over it. */
qaws_status qaws_geodesic_enclosure_get_step(qaws_geodesic_enclosure const* enclosure, unsigned int index, qaws_scalar* out_t0,
	qaws_scalar* out_t1, qaws_vec2* out_lo, qaws_vec2* out_hi);
/* The end: boxes holding the point and the unit (u, v) velocity there. */
qaws_status qaws_geodesic_enclosure_get_end(qaws_geodesic_enclosure const* enclosure, qaws_vec2* out_point_lo, qaws_vec2* out_point_hi,
	qaws_vec2* out_velocity_lo, qaws_vec2* out_velocity_hi);
/* The crossings of two enclosed geodesics of one surface, certified: none
   where no step boxes meet; exactly one in each stretch where they meet,
   when there the two are transversal (the cross product of their velocity
   boxes has one sign, so no second crossing) and each passes from one side
   of the other to the other (Miranda on A(s) - B(t), seen across each
   curve's mean direction, so one crossing). Their boxes (the hull of the
   meeting boxes' common parts), up to capacity; the count. CERTIFICATION_
   FAILED when a stretch cannot be told (a tangency, an overlap, boxes too
   wide). */
qaws_status qaws_geodesic_enclosure_cross(qaws_geodesic_enclosure const* a, qaws_geodesic_enclosure const* b, unsigned int capacity,
	qaws_vec2* out_lo, qaws_vec2* out_hi, unsigned int* out_count);

/* A geodesic segment: from start, leaving in the (u, v) direction, for the
   surface length length. */
typedef struct qaws_geodesic_segment
{
	qaws_vec2 start;
	qaws_vec2 direction;
	qaws_scalar length;
} qaws_geodesic_segment;

/* The crossings of geodesic segments of a polynomial surface certified
   against their float (u, v) images: each segment enclosed
   (qaws_surface_geodesic_enclose, width as there), every pair's crossings
   certified (qaws_geodesic_enclosure_cross), and the images
   (qaws_surface_curve_create_geodesic_ray) crossing exactly as often, each
   crossing inside its certified box.
   Corners: segments share an end only where declared, never because their
   ends are near. next[i] = j declares a polygon corner: segment i ends
   where segment j starts (next NULL, or QAWS_GEODESIC_NO_CORNER: none);
   two segments whose starts are the same point exactly share it too. A
   shot segment's true end is not exactly the next start, so a declared
   corner is a contract that changes the geometry: there the two tubes must
   meet (else CERTIFICATION_FAILED), the stretch where they meet is proven
   by transversality to hold at most one meeting, and the corner is taken
   as simple (the segments cut at their meeting, or joined across the gap
   when they do not meet). Inside that stretch's box crossings are neither
   counted nor checked. Undeclared ends are crossed like any other stretch:
   a segment ending too near another to tell whether they cross fails.
   Then a Boolean of regions bounded by these segments, taken exactly on
   their images (qaws_surface_clip_exact), has the true geodesics' topology
   away from the declared corners. out_crossings: the crossings certified
   (may be NULL). CERTIFICATION_FAILED when a pair cannot be told apart, a
   declared corner is not met, or the images disagree. */
#define QAWS_GEODESIC_NO_CORNER 0xffffffffu
qaws_status qaws_surface_geodesics_certify(qaws_surface const* surface, qaws_geodesic_segment const* segments, unsigned int count,
	unsigned int const* next, qaws_scalar width, unsigned int* out_crossings);

/* Rigorous bounds on the surface area of a geodesic polygon of a polynomial
   surface: segments in order round it, each enclosed (width as in
   qaws_surface_geodesic_enclose). The polygon is the closed curve of the
   segments, each true end joined to the next start by a straight (u, v)
   segment (the gap, a hair long for shot corners); its region is where
   that curve winds round a nonzero number of times. A grid of resolution
   cells (0: 256) along the longer side of the boxes round the curve; each
   free component of cells (clear of every box) is classified by its
   winding number, proven on a polyline through the boxes' common parts
   (homotopic to the curve within the boxes), each side test in intervals;
   inside cells count in both bounds, the boxes' cells in the upper, each
   by the interval of |Su x Sv| over it. CERTIFICATION_FAILED when a
   winding number cannot be proven. */
qaws_status qaws_surface_geodesic_polygon_area(qaws_surface const* surface, qaws_geodesic_segment const* segments, unsigned int count,
	qaws_scalar width, unsigned int resolution, qaws_scalar* out_lo, qaws_scalar* out_hi);

void qaws_surface_curve_destroy(qaws_surface_curve* curve);

qaws_surface_curve_kind qaws_surface_curve_get_kind(qaws_surface_curve const* curve);
qaws_surface const* qaws_surface_curve_get_surface(qaws_surface_curve const* curve);

/* The (u, v) image: a 2D curve owned by the surface curve (parameter range
   from qaws_curve_get_parameter_range). */
qaws_curve const* qaws_surface_curve_get_uv(qaws_surface_curve const* curve);

/* The surface length (geodesics: exact to the tolerance; others by Gauss
   quadrature of |S_u u' + S_v v'|). */
qaws_status qaws_surface_curve_get_length(qaws_surface_curve const* curve, qaws_scalar* out_length);

/* The point at the (u, v) image's parameter t, in (u, v) and on the surface. */
qaws_status qaws_surface_curve_evaluate(qaws_surface_curve const* curve, qaws_scalar t, qaws_vec2* out_uv,
	qaws_vec3* out_point);

/* ------------------------------------------------------------------------ */
/*  Regions                                                                 */
/* ------------------------------------------------------------------------ */

/*
 * A region on the surface: closed paths of surface curves under a fill rule,
 * or a Boolean of two such sets. A patch maps (u, v) one to one onto the
 * surface, so curves cross on the surface exactly where their (u, v) images
 * cross, and a point is inside exactly when its (u, v) is: the 2D engine
 * (qaws_clip) on the images gives the regions' topology, the curves being
 * the true ones. Lengths and areas use the metric:
 *
 *     area = double integral of sqrt(E G - F^2) du dv,
 *
 * by Green's theorem along the boundary (the integral along u of the area
 * element, by Gauss rules, carried round the loops).
 */

typedef struct qaws_surface_path
{
	qaws_surface_curve const* const* curves;   /* consecutive curves, end to start */
	unsigned int curve_count;
	int closed;
} qaws_surface_path;

typedef struct qaws_surface_region qaws_surface_region;

/* The subjects under the fill rule, combined with the clips by the clip
   type (QAWS_CLIP_NONE or no clips: the subjects' union under the rule).
   Every curve must lie on `surface`. */
qaws_status qaws_surface_clip(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	qaws_surface_region** out_region);

void qaws_surface_region_destroy(qaws_surface_region* region);

/* The same with exact topology: every Boolean of the (u, v) images (the
   cuts into period cells and the unions there too) through qaws_clip_exact
   on the lattice `exact` (NULL: 2^-32 of the domain's extent with a period
   either side; 2^-26 in float builds). Each decision is exact or certified
   for the curves as given: drawn curves exactly, geodesics as their (u, v)
   Hermite images (within the tolerance of the true ones). Measures and
   queries on the region are those of qaws_surface_clip's.
   QAWS_STATUS_CERTIFICATION_FAILED when a decision cannot be proven. */
struct qaws_exact_desc;
qaws_status qaws_surface_clip_exact(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	struct qaws_exact_desc const* exact, qaws_surface_region** out_region);

/* The boundary loops in (u, v): outside loops counter-clockwise, holes
   clockwise (views valid while the region lives). */
unsigned int qaws_surface_region_get_path_count(qaws_surface_region const* region);
qaws_status qaws_surface_region_get_path(qaws_surface_region const* region, unsigned int index, qaws_path_2d* out_uv_path);
int qaws_surface_region_is_hole(qaws_surface_region const* region, unsigned int index);

/* The region within the surface distance |distance| of the region bounded
   by the closed paths under the fill rule (distance > 0, grown) or the part
   of it farther than |distance| from its boundary (distance < 0, shrunk),
   with round joins: geodesic offsets of the boundary curves, geodesic
   circles round the corners. Open paths give their tubes of radius
   |distance| with round ends. Exact in the topology like qaws_surface_clip;
   each boundary curve must be smooth (corners between curves). */
qaws_status qaws_surface_offset(qaws_surface const* surface, qaws_fill_rule fill_rule,
	qaws_surface_path const* paths, unsigned int path_count, qaws_scalar distance, qaws_scalar tolerance,
	qaws_surface_region** out_region);

/* The surface area, and the boundary's surface length. */
qaws_status qaws_surface_region_get_area(qaws_surface_region const* region, qaws_scalar* out_area);
qaws_status qaws_surface_region_get_perimeter(qaws_surface_region const* region, qaws_scalar* out_length);

/* The integral over the region of dim values f(user, u, v, point) against
   the surface's area element (a density, a moment, a curvature). */
typedef void (*qaws_surface_region_fn)(void* user, qaws_scalar u, qaws_scalar v, qaws_vec3 point, qaws_scalar* out_values);
qaws_status qaws_surface_region_integrate(qaws_surface_region const* region, qaws_surface_region_fn fn, void* user,
	unsigned int dim, qaws_scalar* out_values);

/* 1 inside, 0 outside, -1 on the boundary (within the tolerance), at the
   parameter point uv. */
qaws_status qaws_surface_region_contains(qaws_surface_region const* region, qaws_vec2 uv, int* out_inside);

qaws_surface const* qaws_surface_region_get_surface(qaws_surface_region const* region);

/* ------------------------------------------------------------------------ */
/*  Geodesic queries                                                        */
/* ------------------------------------------------------------------------ */

/* The geodesic distance from a to b: the length of the geodesic between
   them (as qaws_surface_curve_create_geodesic: locally shortest, the one
   the chord leads to). */
qaws_status qaws_surface_distance(qaws_surface const* surface, qaws_vec2 a, qaws_vec2 b, qaws_scalar tolerance,
	qaws_scalar* out_distance);

/* The point of a curve nearest p along the surface: its parameter (of the
   (u, v) image), its (u, v) and the geodesic distance (any may be NULL).
   Inside the curve it is where the geodesic to p leaves the curve square
   (a root of that angle, bracketed round the nearest of 64 samples); the
   curve's ends are taken too. */
qaws_status qaws_surface_curve_closest_point(qaws_surface_curve const* curve, qaws_vec2 p, qaws_scalar tolerance,
	qaws_scalar* out_parameter, qaws_vec2* out_point, qaws_scalar* out_distance);

typedef struct qaws_surface_crossing
{
	qaws_scalar parameter_a;    /* of curve a's (u, v) image */
	qaws_scalar parameter_b;
	qaws_vec2 point;            /* in a's (u, v) */
} qaws_surface_crossing;

/* Where two curves of one surface cross: their (u, v) images met (on a
   closed surface, b's period copies too), up to capacity; the count of all. */
qaws_status qaws_surface_curve_intersect(qaws_surface_curve const* a, qaws_surface_curve const* b, qaws_surface_crossing* out_crossings,
	unsigned int capacity, unsigned int* out_count);

/* The closest points of two curves and their geodesic distance: a crossing
   when they cross (distance 0), else the pair each nearest the other
   (closest points taken in turn from the nearest of 64 x 64 samples). Any
   output may be NULL. */
qaws_status qaws_surface_curve_closest_points(qaws_surface_curve const* a, qaws_surface_curve const* b, qaws_scalar tolerance,
	qaws_scalar* out_parameter_a, qaws_scalar* out_parameter_b, qaws_vec2* out_point_a, qaws_vec2* out_point_b,
	qaws_scalar* out_distance);

/* A point's signed geodesic distance to a region: to its boundary (cuts
   along seams left out), negative inside; the boundary's nearest point
   (may be NULL). */
qaws_status qaws_surface_region_distance(qaws_surface_region const* region, qaws_vec2 p, qaws_scalar tolerance,
	qaws_scalar* out_distance, qaws_vec2* out_nearest);

/* The shortest path from p0 to p1 avoiding a region (NULL: none, the
   shortest geodesic, searched from a distance field as the shells' are): the
   geodesics touching the region's boundary joined along it, as
   qaws_shell_curve_create_shortest_path on the surface trimmed to the
   region's outside (its domain's edges, or the rims of a closed direction,
   bounding it). One curve drawn in (u, v), its pieces end to end.
   UNSUPPORTED_OPERATION when the region's boundary is cut along a seam. */
qaws_status qaws_surface_curve_create_shortest_path(qaws_surface const* surface, qaws_surface_region const* obstacles, qaws_vec2 p0,
	qaws_vec2 p1, qaws_scalar tolerance, qaws_surface_curve** out_curve);

/* The region as a trimmed surface (its loops borrowed: the region must
   outlive it). */
qaws_status qaws_surface_region_create_trimmed(qaws_surface_region const* region, qaws_surface** out_surface);

/* ------------------------------------------------------------------------ */
/*  Sampling and meshing                                                    */
/* ------------------------------------------------------------------------ */

/* Points over the region uniform in surface area, or in proportion to a
   density per unit area: the region's (u, v) loops through qaws_region_
   sampler with the area element sqrt(E G - F^2) (times the density) as its
   density, so the patterns of qaws_region_sample.h (QAWS_REGION_SAMPLE_*:
   random, Sobol, stratified, Halton, R2) keep their structure through a
   mass-keeping map. QAWS_REGION_SAMPLE_BLUE_NOISE is sample elimination
   (Yuksel 2015) on the surface's own distances (in space), uniform in area
   only. The BOX flag is not supported. */
typedef struct qaws_surface_sample_desc
{
	unsigned int kind;                  /* QAWS_REGION_SAMPLE_* */
	unsigned long long seed;
	qaws_scalar (*density)(void* user, qaws_scalar u, qaws_scalar v);  /* optional, >= 0, per unit surface area */
	void* density_user;
} qaws_surface_sample_desc;

/* count points; out_uv and out_points (either may be NULL) take them in
   (u, v) and on the surface. desc NULL: random, seed 0. */
qaws_status qaws_surface_region_sample(qaws_surface_region const* region, qaws_surface_sample_desc const* desc,
	unsigned int count, qaws_vec2* out_uv, qaws_vec3* out_points);

/* A triangle mesh of the region in (u, v) (qaws_triangulate.h) sized by the
   metric: edges at most about edge_length on the surface, through the size
   field edge_length / sqrt(the metric's largest eigenvalue). Its vertices
   on the surface: qaws_surface_evaluate. On closed surfaces the cuts along
   seams are boundary edges of the (u, v) mesh, their vertices one point on
   the surface. */
typedef struct qaws_surface_mesh_desc
{
	qaws_scalar edge_length;            /* surface length; 0 = 1/30 of the region's size */
	qaws_scalar tolerance;              /* boundary deviation on the surface; 0 = edge_length / 100 */
	qaws_scalar min_angle;              /* degrees in (u, v), 0 = 25 */
	unsigned int remesh_iterations;     /* Botsch-Kobbelt rounds after refinement, 0 = none */
} qaws_surface_mesh_desc;

struct qaws_mesh_2d;
qaws_status qaws_surface_region_triangulate(qaws_surface_region const* region, qaws_surface_mesh_desc const* desc,
	struct qaws_mesh_2d** out_mesh);

/* ------------------------------------------------------------------------ */
/*  Medial axis                                                             */
/* ------------------------------------------------------------------------ */

/* The medial axis of the region on the surface: the points at the same
   geodesic distance from two or more nearest boundary points, the centres
   of the largest geodesic circles inside.
   - The boundary is taken as its smooth pieces (split at corners; cuts
     along seams and poles left out) and their ends.
   - Conormal geodesics traced inward from samples of the pieces, and
     fans of geodesics from the ends, label a mesh of the region (sized by
     the metric) with each vertex's nearest feature.
   - A mesh edge whose ends have nearest boundary points far apart crosses
     the axis there: found on the edge where the two features' distances
     agree, each a true geodesic distance (a piece's by Newton on its
     foot point, an end's by shooting).
   - Within a triangle two crossings are joined by a segment; where three
     features meet, the junction (equally far from all three) by Newton.
   The axis is resolved to the resolution (the mesh's edge length on the
   surface) along it; its points are on the true axis within the
   tolerance. */
typedef struct qaws_surface_medial_axis qaws_surface_medial_axis;

typedef struct qaws_surface_medial_desc
{
	qaws_scalar resolution;     /* edge length on the surface; 0 = 1/40 of sqrt(area) */
} qaws_surface_medial_desc;

qaws_status qaws_surface_region_medial_axis(qaws_surface_region const* region, qaws_surface_medial_desc const* desc,
	qaws_surface_medial_axis** out_axis);
void qaws_surface_medial_axis_destroy(qaws_surface_medial_axis* axis);

/* The axis as segments in (u, v), with the radii (distances to the
   boundary) at their ends. */
unsigned int qaws_surface_medial_axis_get_segment_count(qaws_surface_medial_axis const* axis);
qaws_status qaws_surface_medial_axis_get_segment(qaws_surface_medial_axis const* axis, unsigned int index,
	qaws_vec2 out_uv[2], qaws_scalar out_radius[2]);

/* The junctions: points equally far from three boundary features. */
unsigned int qaws_surface_medial_axis_get_junction_count(qaws_surface_medial_axis const* axis);
qaws_status qaws_surface_medial_axis_get_junction(qaws_surface_medial_axis const* axis, unsigned int index,
	qaws_vec2* out_uv, qaws_scalar* out_radius);

/* The largest geodesic circle inside the region found on the axis. */
qaws_status qaws_surface_medial_axis_get_max_inscribed(qaws_surface_medial_axis const* axis, qaws_vec2* out_centre,
	qaws_scalar* out_radius);

/* The mitered (straight) skeleton of a region with corners: each edge (the
   boundary between two corners) moves as its line, the edge carried on
   past both ends along the geodesics of its tangents, at unit speed; the
   fronts meet at mitres. A point is swept first by the edge nearest it
   among those whose wedge holds it (at a convex corner nearer the edge's
   line than the neighbour's, at a reflex one farther); the skeleton is
   where that edge changes, found on a mesh as for the medial axis (labels,
   crossings by bisection on the mesh edges, junctions by Newton), with
   the time (the radius) at each point. The same accessors as the medial
   axis. Exact for convex regions and reflex corners whose fronts are not
   split (no split events: a front cut off by another counts on beyond the
   cut). UNSUPPORTED_OPERATION for a loop without a corner (its medial
   axis is the answer) or one cut along a seam or a pole. */
qaws_status qaws_surface_region_mitered_skeleton(qaws_surface_region const* region, qaws_surface_medial_desc const* desc,
	qaws_surface_medial_axis** out_axis);

#endif /* QAWS_SURFACE_REGION_H */
