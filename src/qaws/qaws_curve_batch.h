#ifndef QAWS_CURVE_BATCH_H
#define QAWS_CURVE_BATCH_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_batch_executor.h"

/*
 * Intersections of N curves at once.
 *
 * Every curve is flattened once into chord segments whose boxes enclose the
 * curve (inflated by the measured chord deviation). One uniform grid over
 * all segments of all curves gives the candidate segment pairs, and only
 * those pairs are refined by Newton on the true curves. The cost grows with
 * the total segment count plus the number of nearby segment pairs, not with
 * the number of curve pairs.
 *
 * Families: when `families` is given, two curves with the same family id are
 * never intersected with each other. For a heightfield, give the contour
 * lines one family and the gradient lines another: contours of different
 * levels never cross, and gradient lines meet only at critical points.
 *
 * Thread-safe on immutable curves.
 */

#define QAWS_CURVE_BATCH_SELF 1u   /* also intersect each curve with itself */

typedef struct qaws_curve_batch_desc
{
	qaws_curve const* const* curves;
	unsigned int curve_count;
	unsigned int const* families;   /* optional, one id per curve */
	unsigned int flags;             /* QAWS_CURVE_BATCH_* */
	qaws_scalar flatness;           /* chord deviation bound; 0 = 2^-10 of the scene extent */
	qaws_batch_executor const* executor;   /* optional: runs the work in parallel */
} qaws_curve_batch_desc;

/* Hit kinds. A CROSSING passes from one side to the other (transversal, or
   tangent with an odd contact such as at an inflection); a TOUCH meets
   without crossing (an even tangency); an OVERLAP is a stretch the two
   curves share: [parameter_a, parameter_a_end] on curve_a (increasing) and
   [parameter_b, parameter_b_end] on curve_b (decreasing when they run
   opposite ways), from `position` to the end point. 3D hits are CROSSING or
   OVERLAP. */
#define QAWS_CURVE_HIT_CROSSING 0u
#define QAWS_CURVE_HIT_TOUCH    1u
#define QAWS_CURVE_HIT_OVERLAP  2u

/* curve_a < curve_b, or curve_a == curve_b with parameter_a < parameter_b
   for a self-intersection. Sorted by (curve_a, curve_b, parameter_a). */
typedef struct qaws_curve_batch_hit_2d
{
	unsigned int curve_a;
	unsigned int curve_b;
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec2 position;
	unsigned int kind;              /* QAWS_CURVE_HIT_* */
	qaws_scalar parameter_a_end;    /* OVERLAP: other end of the shared stretch; else parameter_a */
	qaws_scalar parameter_b_end;    /* OVERLAP: other end on curve_b; else parameter_b */
} qaws_curve_batch_hit_2d;

typedef struct qaws_curve_batch_hit_3d
{
	unsigned int curve_a;
	unsigned int curve_b;
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec3 position;
	unsigned int kind;              /* QAWS_CURVE_HIT_* */
	qaws_scalar parameter_a_end;    /* OVERLAP: other end of the shared stretch; else parameter_a */
	qaws_scalar parameter_b_end;    /* OVERLAP: other end on curve_b; else parameter_b */
} qaws_curve_batch_hit_3d;

typedef struct qaws_curve_batch_stats
{
	unsigned int segment_count;     /* chord segments over all curves */
	unsigned int cell_count;        /* grid cells */
	unsigned int candidate_count;   /* segment pairs whose boxes overlap */
	unsigned int newton_count;      /* candidates close enough to refine */
	unsigned int hit_count;         /* distinct intersections */
} qaws_curve_batch_stats;

/* Writes min(count, capacity) hits; *out_count is the total found.
   out_stats may be NULL. */
qaws_status qaws_curve_batch_find_intersections_2d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_2d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats);

qaws_status qaws_curve_batch_find_intersections_3d(
	qaws_curve_batch_desc const* desc,
	qaws_curve_batch_hit_3d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats);

/*
 * Level crossings: where N curves cross the level sets h = L of a scalar
 * field, for every level at once, without building the level curves. For
 * a heightfield, the gradient lines against the contour levels: the
 * crossings of every contour with every gradient line.
 *
 * Each curve is flattened once; along each segment the field values at its
 * ends bracket the levels it crosses (a segment is split where h is not
 * close to linear along it, so a level crossed twice inside it is not
 * missed); each bracket is solved for h(C(t)) = L by Newton with
 * bisection safeguards (Newton needs the gradient; without it, secant).
 * A curve touching a level without crossing it is not reported.
 */
typedef qaws_scalar (*qaws_scalar_field_fn)(
	void* user,
	qaws_scalar const* point,       /* x, y (, z) */
	qaws_scalar* gradient);         /* NULL, or receives dh/dx, dh/dy (, dh/dz) */

typedef struct qaws_level_crossing_desc
{
	qaws_curve const* const* curves;    /* all 2D or all 3D */
	unsigned int curve_count;
	qaws_scalar_field_fn field;
	void* user;
	qaws_scalar const* levels;          /* strictly ascending */
	unsigned int level_count;
	qaws_scalar flatness;               /* geometric, as for the batch; 0 = 2^-10 of the extent */
	qaws_batch_executor const* executor;   /* optional: runs the work in parallel */
} qaws_level_crossing_desc;

/* sorted by (curve, parameter) */
typedef struct qaws_level_crossing
{
	unsigned int curve;
	unsigned int level;
	qaws_scalar parameter;
	qaws_vec3 position;                 /* z = 0 for 2D */
} qaws_level_crossing;

/* Writes min(count, capacity) crossings; *out_count is the total. */
qaws_status qaws_curve_batch_find_level_crossings(
	qaws_level_crossing_desc const* desc,
	qaws_level_crossing* out_crossings,
	unsigned int capacity,
	unsigned int* out_count);

/*
 * Closest points: for every query point, the nearest point over N curves
 * at once. One grid over the flattened segments of all curves is searched
 * ring by ring outward from each point, a segment kept while its chord
 * distance minus its inflation can beat the best chord distance plus
 * inflation seen; the survivors are refined by Newton on
 * (C(t) - p) . C'(t) = 0 (clamped to the curve's domain) and the nearest
 * true point wins (the lower curve index on a tie).
 */
#define QAWS_CURVE_BATCH_NONE 0xFFFFFFFFu

typedef struct qaws_closest_desc
{
	qaws_curve const* const* curves;    /* all 2D or all 3D */
	unsigned int curve_count;
	qaws_scalar const* points;          /* point_count points of 2 or 3 scalars */
	unsigned int point_count;
	qaws_scalar max_distance;           /* 0 = no limit */
	qaws_scalar flatness;               /* 0 = 2^-10 of the extent */
	qaws_batch_executor const* executor;   /* optional: runs the work in parallel */
} qaws_closest_desc;

typedef struct qaws_closest_point
{
	unsigned int curve;                 /* QAWS_CURVE_BATCH_NONE: no curve within max_distance */
	qaws_scalar parameter;
	qaws_scalar distance;
	qaws_vec3 position;                 /* z = 0 for 2D */
} qaws_closest_point;

/* out_points receives point_count results. out_stats may be NULL
   (candidate_count: segments refined, hit_count: points with a curve). */
qaws_status qaws_curve_batch_find_closest(
	qaws_closest_desc const* desc,
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats);

/*
 * Prepared sets: the flattening (the costly part: curve evaluations) done
 * once, for curves queried again and again, e.g. fixed contour lines
 * against gradient lines that change. A set keeps pointers to its curves:
 * they must outlive it and stay unchanged.
 *
 * qaws_curve_set_find_intersections_*(set, NULL, ...) intersects a set
 * with itself, by its desc (families, flags) as the one-shot call;
 * (set, other, ...) intersects every curve of set with every curve of
 * other: curve_a indexes set, curve_b other.
 */
typedef struct qaws_curve_set qaws_curve_set;

qaws_status qaws_curve_set_create(qaws_curve_batch_desc const* desc, qaws_curve_set** out_set);
void qaws_curve_set_destroy(qaws_curve_set* set);
unsigned int qaws_curve_set_get_segment_count(qaws_curve_set const* set);

qaws_status qaws_curve_set_find_intersections_2d(
	qaws_curve_set const* set,
	qaws_curve_set const* other,
	qaws_curve_batch_hit_2d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats);

/* closest points of a prepared set (its own flatness) */
qaws_status qaws_curve_set_find_closest(
	qaws_curve_set const* set,
	qaws_scalar const* points,
	unsigned int point_count,
	qaws_scalar max_distance,
	qaws_closest_point* out_points,
	qaws_curve_batch_stats* out_stats);

qaws_status qaws_curve_set_find_intersections_3d(
	qaws_curve_set const* set,
	qaws_curve_set const* other,
	qaws_curve_batch_hit_3d* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_curve_batch_stats* out_stats);

#endif /* QAWS_CURVE_BATCH_H */
