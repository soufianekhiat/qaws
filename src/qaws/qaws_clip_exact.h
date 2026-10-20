#ifndef QAWS_CLIP_EXACT_H
#define QAWS_CLIP_EXACT_H

#include "qaws_clip.h"
#include "qaws_exact.h"

/*
 * Exact Boolean operations on regions of curves: qaws_clip's arrangement
 * built on the exact kernel (qaws_exact.h), so its topology does not depend
 * on a tolerance.
 *
 * The input curves are prepared onto one exact lattice (qaws_exact_desc);
 * from there every topological decision is exact or certified:
 *
 *   - hits: the certified batch (crossings, tangencies with their
 *     multiplicity, overlaps with exact ends, self-intersections);
 *   - vertices: hits at the same point are one vertex when proven so: exact
 *     rational points compared exactly, an exact point on a curve by its
 *     implicit equation, two algebraic points on one curve by the gcd of
 *     their defining polynomials;
 *   - edges: each curve cut at its vertices; stretches that overlap are one
 *     edge carrying the winding steps of all of them;
 *   - order round a vertex: from the exact sub-curve leaving it (integer de
 *     Casteljau at the rational parameter): the first non-zero derivative,
 *     then the curvature for curves leaving tangent; at an irrational
 *     vertex, the tangent cones of the sub-curve over its enclosure
 *     (Sederberg and Wang), which must be apart;
 *   - faces: the cycles of the planar graph; each part lies in the face of
 *     the others found by a ray from an exact point of it (crossings
 *     certified, a new ray when it passes too close to a vertex);
 *   - windings from the unbounded face across edges, Clipper2's truth
 *     tables on the fill rule.
 *
 * The result is a qaws_clip_result (paths of new curves on the inputs, their
 * nesting, the vertices), read with the qaws_clip_result_* functions; the
 * vertices and the new curves' ends are the doubles nearest the exact
 * points. Cycle orientations come from areas in double (the topology does
 * not change with them).
 *
 * QAWS_STATUS_CERTIFICATION_FAILED when a decision cannot be proven: two
 * irrational vertices on one curve that the gcd cannot tell apart, curves
 * leaving an irrational vertex tangent to each other, a contact through a
 * singular point, or curves leaving a vertex with the same curvature.
 * Supported curves: those of qaws_exact_curve_prepare; polylines and
 * composites are taken piece by piece. Open subjects are not supported yet.
 */

typedef struct qaws_clip_exact_desc
{
	qaws_path_2d const* subjects;
	unsigned int subject_count;
	qaws_path_2d const* clips;
	unsigned int clip_count;
	qaws_clip_type clip_type;
	qaws_fill_rule fill_rule;
	unsigned int flags;                     /* QAWS_CLIP_REVERSE_SOLUTION, QAWS_CLIP_PRESERVE_COLLINEAR */
	qaws_exact_desc const* exact;           /* the lattice; NULL: qaws_exact_desc_default */
	qaws_batch_executor const* executor;    /* optional, for the batched intersection */
	qaws_clip_z_fn z_fn;                    /* optional */
	void* z_user;
} qaws_clip_exact_desc;

typedef struct qaws_clip_exact_stats
{
	unsigned int source_count;      /* curves after polylines and composites are split */
	unsigned int hit_count;         /* from the batch */
	unsigned int tangency_count;    /* hits of even multiplicity */
	unsigned int overlap_count;
	unsigned int vertex_count;
	unsigned int edge_count;        /* after overlaps are merged and cancelled ones dropped */
	unsigned int face_count;        /* bounded faces + the unbounded one */
	unsigned int component_count;
	unsigned int identity_tests;    /* vertices proven equal or apart by an exact test */
	unsigned int ray_retries;       /* rays given up for passing too close to a vertex */
} qaws_clip_exact_stats;

qaws_status qaws_clip_exact_execute(qaws_clip_result** out_result, qaws_clip_exact_stats* out_stats, qaws_clip_exact_desc const* desc);

#endif /* QAWS_CLIP_EXACT_H */
