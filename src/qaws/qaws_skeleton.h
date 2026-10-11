#ifndef QAWS_SKELETON_H
#define QAWS_SKELETON_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"

/*
 * The straight skeleton of a region (Aichholzer and Aurenhammer): the
 * boundary moves inwards, every edge parallel to itself at its own speed
 * (its weight), and the skeleton is the trace of the moving boundary's
 * vertices. Each input edge sweeps one face; the time a point is swept is
 * the roof's height there.
 *
 * Space-time form: the edge on the line n . p = d (n the inward unit
 * normal) at time t is on n . p = d + w t, a plane over (x, y, t). A
 * wavefront vertex is two planes meeting, a line in space-time, and every
 * event is three planes meeting at a point: the same 3 x 3 solve for edge,
 * split and vertex events.
 *
 * Events follow Kelly and Wonka's generalised intersection event
 * (procedural extrusions, TOG 2011; twak/campskeleton). Every wavefront
 * vertex and edge at the event's point and time is taken at once. Runs of
 * consecutive vertices there form chains, and an edge through the point
 * is a chain of its own. The edges inside a chain end there. The chains,
 * sorted by angle round the point, are joined last edge to the next
 * chain's first. Coincident edge, split and vertex events, and several
 * vertices meeting, are one case with no special handling.
 *
 * Antiparallel edges that meet along a stretch (a corridor of constant
 * width closing all at once) collapse along their common line at the
 * event's time; the skeleton gets that stretch.
 *
 * The region is the paths united under the fill rule (qaws_clip), so outer
 * boundaries run counter-clockwise and holes clockwise. Its curves must be
 * straight (degree 1 Beziers or B-splines, polylines). Weights must be > 0;
 * consecutive collinear edges with different weights are ambiguous
 * (Biedl et al., CGTA 2015), and their vertex moves at the average speed.
 *
 * Thread-safe on an immutable skeleton.
 */

typedef struct qaws_skeleton qaws_skeleton;

/* The weight (speed) of the input edge on span `span` of curve `curve` of
   path `path`; a span is a Bezier's whole curve or a B-spline's knot span. */
typedef qaws_scalar (*qaws_skeleton_weight_fn)(void* user, unsigned int path, unsigned int curve, unsigned int span);

typedef struct qaws_skeleton_desc
{
	qaws_path_2d const* paths;      /* closed paths */
	unsigned int path_count;
	qaws_fill_rule fill_rule;
	qaws_skeleton_weight_fn weight_fn;   /* optional: NULL = every weight 1 */
	void* weight_user;
	qaws_scalar tolerance;          /* events this close are one; 0 = 1e-9 of the extent (1e-6 in float) */
} qaws_skeleton_desc;

/* A skeleton node: where wavefront vertices start (the input's corners, at
   time 0) or meet. */
typedef struct qaws_skeleton_node
{
	qaws_vec2 position;
	qaws_scalar time;
} qaws_skeleton_node;

/* A skeleton arc: one wavefront vertex's path from node_a to node_b, between
   the faces of its two edges (or a collapse along antiparallel edges). */
typedef struct qaws_skeleton_arc
{
	unsigned int node_a, node_b;
	unsigned int face_a, face_b;
} qaws_skeleton_arc;

/* A face: the part of the region its input edge sweeps. Its boundary is
   the input edge from node 0 to node 1, then the skeleton round to node 0
   (counter-clockwise). The source names the input edge. */
typedef struct qaws_skeleton_face
{
	unsigned int const* nodes;
	unsigned int node_count;
	unsigned int path, curve, span;
	qaws_scalar weight;
} qaws_skeleton_face;

qaws_status qaws_skeleton_create(qaws_skeleton_desc const* desc, qaws_skeleton** out_skeleton);
void qaws_skeleton_destroy(qaws_skeleton* skeleton);

unsigned int qaws_skeleton_get_node_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_node(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_node* out_node);
unsigned int qaws_skeleton_get_arc_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_arc(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_arc* out_arc);
unsigned int qaws_skeleton_get_face_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_face(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_face* out_face);

/* The time the region is swept out (the roof's top). */
qaws_scalar qaws_skeleton_get_max_time(qaws_skeleton const* skeleton);

/* The wavefront at time t >= 0 (the mitered inward offset by t, weighted):
   its loops as polygons, outer counter-clockwise, holes clockwise, as a
   qaws_clip result. Empty past the max time. */
qaws_status qaws_skeleton_offset(qaws_skeleton const* skeleton, qaws_scalar time, qaws_clip_result** out_result);

#endif /* QAWS_SKELETON_H */
