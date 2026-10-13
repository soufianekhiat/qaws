#ifndef QAWS_SURFACE_REGION_H
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

typedef enum qaws_surface_curve_kind
{
	QAWS_SURFACE_CURVE_UV = 0,          /* a curve drawn in (u, v) */
	QAWS_SURFACE_CURVE_GEODESIC,        /* a geodesic segment (or ray) */
	QAWS_SURFACE_CURVE_GEODESIC_ARC     /* an arc of a geodesic circle */
} qaws_surface_curve_kind;

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

/* The boundary loops in (u, v): outside loops counter-clockwise, holes
   clockwise (views valid while the region lives). */
unsigned int qaws_surface_region_get_path_count(qaws_surface_region const* region);
qaws_status qaws_surface_region_get_path(qaws_surface_region const* region, unsigned int index, qaws_path_2d* out_uv_path);
int qaws_surface_region_is_hole(qaws_surface_region const* region, unsigned int index);

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

/* The region as a trimmed surface (its loops borrowed: the region must
   outlive it). */
qaws_status qaws_surface_region_create_trimmed(qaws_surface_region const* region, qaws_surface** out_surface);

#endif /* QAWS_SURFACE_REGION_H */
