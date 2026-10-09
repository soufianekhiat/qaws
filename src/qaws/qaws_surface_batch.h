#ifndef QAWS_SURFACE_BATCH_H
#define QAWS_SURFACE_BATCH_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_surface_intersect.h"
#include "qaws_curve_batch.h"

/*
 * Intersections of N 3D curves with M surfaces at once.
 *
 * Curves are flattened into chord segments and surfaces into (u, v)
 * patches, each with a box inflated by its measured deviation. One grid
 * over every segment and patch gives the segment / patch pairs that may
 * meet; each is seeded by intersecting the segment with the patch's two
 * triangles and refined by Newton on S(u, v) = C(t). The cost grows with the
 * pieces and their nearby pairs, not with the curve / surface pairs.
 *
 * Tangential contacts are not reported (singular Jacobian). Thread-safe on
 * immutable curves and surfaces.
 */

typedef struct qaws_curve_surface_batch_desc
{
	qaws_curve const* const* curves;        /* 3D */
	unsigned int curve_count;
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	qaws_scalar flatness;                   /* deviation bound; 0 = 2^-9 of the scene extent */
} qaws_curve_surface_batch_desc;

/* sorted by (curve, surface, t) */
typedef struct qaws_curve_surface_batch_hit
{
	unsigned int curve;
	unsigned int surface;
	qaws_scalar t, u, v;
	qaws_vec3 position;
} qaws_curve_surface_batch_hit;

typedef struct qaws_surface_batch_stats
{
	unsigned int segment_count;     /* curve chord segments */
	unsigned int patch_count;       /* surface patches */
	unsigned int cell_count;        /* grid cells */
	unsigned int candidate_count;   /* piece pairs whose boxes overlap */
	unsigned int newton_count;      /* candidates refined */
	unsigned int hit_count;
} qaws_surface_batch_stats;

/* Writes min(count, capacity) hits; *out_count is the total. out_stats may be NULL. */
qaws_status qaws_curve_surface_batch_find_intersections(
	qaws_curve_surface_batch_desc const* desc,
	qaws_curve_surface_batch_hit* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_surface_batch_stats* out_stats);

/*
 * Intersection curves of N surfaces at once.
 *
 * Every surface is flattened into patches; one grid over all patches gives
 * the patch pairs of different surfaces (and different families, when
 * `families` is given) that may meet. Their triangles are intersected, and
 * the end points of each small segment are moved onto both true surfaces
 * by a minimum-norm Newton on S_a(u1, v1) = S_b(u2, v2). The segments of
 * each surface pair are chained into polylines through end points closer
 * than a few flatness bounds.
 *
 * Tangential contacts and branches closer than the flatness bound are not
 * resolved. Thread-safe on immutable surfaces.
 */
typedef struct qaws_surface_batch_desc
{
	qaws_surface const* const* surfaces;
	unsigned int surface_count;
	unsigned int const* families;   /* optional, one id per surface */
	qaws_scalar flatness;           /* deviation bound; 0 = 2^-9 of the scene extent */
} qaws_surface_batch_desc;

/* points[first .. first + count - 1] of the point buffer; closed when the
   last point joins the first. Sorted by (surface_a, surface_b). */
typedef struct qaws_surface_batch_curve
{
	unsigned int surface_a;         /* surface_a < surface_b */
	unsigned int surface_b;
	unsigned int first, count;
	int closed;
} qaws_surface_batch_curve;

/* Writes min(n, capacity) curves and points (curves whose points do not fit
   are not written); *out_curve_count / *out_point_count are the totals. */
qaws_status qaws_surface_batch_find_intersections(
	qaws_surface_batch_desc const* desc,
	qaws_surface_batch_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* out_points,
	unsigned int point_capacity,
	unsigned int* out_point_count,
	qaws_surface_batch_stats* out_stats);


/*
 * Prepared surface sets: the flattening done once, for surfaces queried
 * again and again (a terrain against changing planes, curves or other
 * surfaces). A set keeps pointers to its surfaces: they must outlive it
 * and stay unchanged.
 *
 * qaws_surface_set_find_intersections(set, NULL, ...) intersects a set with
 * itself by its desc (families) as the one-shot call; (set, other, ...)
 * every surface of set with every surface of other: surface_a indexes
 * set, surface_b other. qaws_curve_set_find_surface_intersections takes a
 * prepared set of 3D curves against a prepared surface set.
 */
typedef struct qaws_surface_set qaws_surface_set;

qaws_status qaws_surface_set_create(qaws_surface_batch_desc const* desc, qaws_surface_set** out_set);
void qaws_surface_set_destroy(qaws_surface_set* set);
unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set);

qaws_status qaws_surface_set_find_intersections(
	qaws_surface_set const* set,
	qaws_surface_set const* other,
	qaws_surface_batch_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* out_points,
	unsigned int point_capacity,
	unsigned int* out_point_count,
	qaws_surface_batch_stats* out_stats);

qaws_status qaws_curve_set_find_surface_intersections(
	qaws_curve_set const* curves,
	qaws_surface_set const* surfaces,
	qaws_curve_surface_batch_hit* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_surface_batch_stats* out_stats);

#endif /* QAWS_SURFACE_BATCH_H */
