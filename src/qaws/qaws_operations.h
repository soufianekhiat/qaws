#ifndef QAWS_OPERATIONS_H
#define QAWS_OPERATIONS_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"

/* Thread-safe: operates on immutable input curves. */

qaws_status qaws_curve_split(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_curve** out_left,
	qaws_curve** out_right);

/*
 * The piece of a curve between parameters t0 and t1, as a new curve that lies
 * exactly on the source (t0 > t1: the piece reversed). The new curve is of the
 * simplest kind that represents the piece exactly:
 *
 *   Bezier, rational Bezier, B-spline, NURBS, polynomial, arc, clothoid:
 *     the same kind (B-spline and NURBS keep the source parameters);
 *   Hermite, Catmull-Rom, trajectory, subdivision:
 *     a cubic Bezier (one span) or cubic B-spline whose knots are the
 *     source parameters;
 *   composite: the piece of one segment, or a composite of pieces;
 *   reparameterized: the piece of its source.
 *
 * Yuksel curves have no polynomial form; their pieces are cubic B-splines
 * fitted to 1e-9 (1e-5 in float builds) of the piece's extent.
 * Parameters within rounding of the domain are clamped onto it.
 */
qaws_status qaws_curve_extract(
	qaws_curve const* curve,
	qaws_scalar t0,
	qaws_scalar t1,
	qaws_curve** out_curve);

qaws_status qaws_curve_join(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_curve** out_joined);

/*
 * Compute the offset (parallel) curve at signed distance.
 * Positive distance offsets to the left of the travel direction.
 * 2D only.
 *
 * At cusps (where curvature radius < |distance|), the raw offset
 * reverses direction. The backward loops are trimmed at their
 * self-intersection points so the result is a clean boundary.
 *
 * trim levels:
 *   0  Raw offset — cusp loop removal only.
 *   1  + Distance-based trim: removes points closer than |distance|
 *      to the curve (inside the swept boundary).
 *   2  + Self-intersection cleanup: detects remaining polyline
 *      self-intersections, extracts each loop as a separate
 *      closed curve (closed at the crossing point), and
 *      discards the tails. Produces clean, non-self-intersecting
 *      closed loops.
 */
qaws_status qaws_curve_offset_2d(
	qaws_curve const* curve,
	qaws_scalar distance,
	int trim,
	qaws_curve** out_curves,
	unsigned int curve_capacity,
	unsigned int* out_count);

/*
 * Create an arc-length reparameterized wrapper curve.
 * The resulting curve has parameter domain [0, total_arc_length]
 * and maps uniformly to arc length along the source curve.
 *
 * The source curve must outlive the returned wrapper (non-owning reference).
 * table_resolution controls the arc-length lookup table size (0 = default 256).
 */
qaws_status qaws_curve_reparameterize_arc_length(
	qaws_curve const* curve,
	unsigned int table_resolution,
	qaws_curve** out_curve);

/*
 * Offset a 3D curve by a distance along a direction field.
 *
 * direction_mode:
 *   0 = constant direction (uses the direction parameter)
 *   1 = curve normal (Frenet frame normal)
 *   2 = surface normal (requires surface parameter)
 *
 * For mode 0, direction must be non-NULL.
 * For mode 2, surface must be non-NULL.
 * sample_count controls the output resolution (0 = default 256).
 */
qaws_status qaws_curve_offset_3d(
	qaws_curve const* curve,
	qaws_scalar distance,
	int direction_mode,
	qaws_vec3 const* direction,
	qaws_surface const* surface,
	unsigned int sample_count,
	qaws_curve** out_curve);

/*
 * Merge a chain of curve segments into fewer, higher-quality curves.
 * The input curves should share endpoints (C0 continuity).
 *
 * target_degree: degree of the output B-spline (0 = default 3).
 * tolerance: maximum allowed deviation from original chain.
 */
qaws_status qaws_curve_merge_chain(
	qaws_curve const* const* curves,
	unsigned int curve_count,
	unsigned int target_degree,
	qaws_scalar tolerance,
	qaws_curve** out_curves,
	unsigned int curve_capacity,
	unsigned int* out_count);

/*
 * Reparameterize curve_b so both curves traverse the same total arc length.
 * The result has parameter domain [0, arc_length(curve_a)].
 * curve_b must outlive the returned wrapper (non-owning reference).
 * table_resolution controls lookup table size (0 = default 256).
 */
qaws_status qaws_curve_match_arc_length(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	unsigned int table_resolution,
	qaws_curve** out_curve);

/*
 * Insert circular arc fillets at sharp corners of a 2D composite curve.
 * radius: fillet arc radius at each corner.
 * The output is a new composite curve with smooth transitions.
 */
qaws_status qaws_curve_fillet_2d(
	qaws_curve const* curve,
	qaws_scalar radius,
	qaws_curve** out_curve);

/*
 * Insert straight-line chamfers at sharp corners of a 2D composite curve.
 * distance: chamfer setback distance from each corner.
 * The output is a new composite curve with chamfered corners.
 */
qaws_status qaws_curve_chamfer_2d(
	qaws_curve const* curve,
	qaws_scalar distance,
	qaws_curve** out_curve);

#endif /* QAWS_OPERATIONS_H */
