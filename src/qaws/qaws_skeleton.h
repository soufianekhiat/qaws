#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SKELETON_H_F32
#define QAWS_SKELETON_H_F32
#define QAWS_SKELETON_H_BODY
#endif
#else
#ifndef QAWS_SKELETON_H_F64
#define QAWS_SKELETON_H_F64
#define QAWS_SKELETON_H_BODY
#endif
#endif
#ifdef QAWS_SKELETON_H_BODY
#undef QAWS_SKELETON_H_BODY
#define QAWS_SKELETON_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"

/*
 * Straight and curved skeletons of a region: the boundary moves inwards,
 * every edge at its own speed (its weight), and the skeleton is the trace of
 * the moving boundary's vertices. Each edge sweeps one face; the time a
 * point is swept is the roof's height there.
 *
 * Edges are straight lines or circular arcs. A line moves parallel to
 * itself (the straight skeleton, Aichholzer and Aurenhammer). An arc
 * shrinks or grows about its centre (the mitered offsets and skeletons of
 * circular arc polygons, Weiss, Juttler and Aurenhammer, EuroCG 2018 /
 * IJCGA 2020): its vertices travel on conics. A reflex vertex whose edges
 * turn tangent, or whose angle reaches the splice angle sigma, opens a
 * circle that grows from its point. sigma = 2 pi (the default) keeps the
 * mitered (straight) skeleton; sigma = pi opens every reflex vertex at
 * once: the medial axis and the classical offsets.
 *
 * Space-time form: over (x, y, t) a line is the plane n . p = d + w t (n its
 * inward normal) and an arc the cone |p - m| = r -+ t. Every event is three
 * of these meeting at a point; the difference of two cones is a plane, so
 * each is a linear system and at most one quadratic. Two supports touching
 * (squeeze), an arc whose ends meet round its circle (bubble) and a circle
 * shrinking to its centre are the other events.
 *
 * Coincident events follow Kelly and Wonka's generalised intersection event
 * (procedural extrusions, TOG 2011; twak/campskeleton): every wavefront
 * vertex and edge at the event's point and time is taken at once. Runs of
 * consecutive vertices there form chains, and an edge through the point is
 * a chain of its own. The edges inside a chain end. The chains, sorted by
 * angle round the point, are joined last edge to the next chain's first.
 * Edges that meet along a stretch (a corridor of constant width closing all
 * at once) collapse along it.
 *
 * The region is the paths united under the fill rule (qaws_clip), so outer
 * boundaries run counter-clockwise and holes clockwise. Lines (degree 1
 * Beziers or B-splines, polylines) and circular arcs (qaws_curve_create_arc)
 * are taken exactly; any other curve (Beziers, NURBS, ellipses, ...) as a G1
 * arc spline: equal-tangent biarcs, halved until within the curve tolerance,
 * whose smooth joints carry no skeleton. The skeleton is then exact for
 * that spline. Line weights must be > 0 (arcs move at speed
 * 1); consecutive collinear lines with different weights are ambiguous
 * (Biedl et al., CGTA 2015), and their vertex moves at the average speed.
 *
 * Thread-safe on an immutable skeleton.
 */

typedef struct qaws_skeleton qaws_skeleton;

/* The weight (speed) of the input line on span `span` of curve `curve` of
   path `path`; a span is a Bezier's whole curve or a B-spline's knot span. */
typedef qaws_scalar (*qaws_skeleton_weight_fn)(void* user, unsigned int path, unsigned int curve, unsigned int span);

typedef struct qaws_skeleton_desc
{
	qaws_path_2d const* paths;      /* closed paths */
	unsigned int path_count;
	qaws_fill_rule fill_rule;
	qaws_skeleton_weight_fn weight_fn;   /* optional: NULL = every weight 1 */
	void* weight_user;
	qaws_scalar tolerance;          /* events this close are one; 0 = 1e-9 of the extent (in float 1e-6, for the vertices' rounding, unless free-form curves are fitted) */
	qaws_scalar splice_angle;       /* sigma in [pi, 2 pi]; 0 = 2 pi (mitered); pi: the medial axis */
	qaws_scalar curve_tolerance;    /* free-form curves as G1 arc splines within this; 0 = 1e-4 of the extent */
} qaws_skeleton_desc;

/* A skeleton node: where wavefront vertices start (the input's corners, at
   time 0) or meet. */
typedef struct qaws_skeleton_node
{
	qaws_vec2 position;
	qaws_scalar time;
} qaws_skeleton_node;

/* A skeleton arc: one wavefront vertex's path from node_a to node_b (a
   segment, or a conic between curved edges), between the faces of its two
   edges; or a collapse along coincident edges. smooth: the path of a
   vertex where the edges meet tangentially (part of the faces' boundaries,
   not of the skeleton proper). boundary: an end on the input (time 0); the
   arcs with neither flag are the pure skeleton, clear of the boundary. */
typedef struct qaws_skeleton_arc
{
	unsigned int node_a, node_b;
	unsigned int face_a, face_b;
	int smooth;
	int boundary;
} qaws_skeleton_arc;

/* A face: the part of the region one support sweeps. The input's edges come
   first, then the circles opened by splices (their faces start at the
   spliced vertex). The boundary runs from node 0 to node 1 along the edge,
   then through the skeleton back to node 0 (counter-clockwise). The
   support: a line n . p = offset (n the inward normal) moving at weight,
   or a circle (center, radius at time 0, side +1 when it shrinks). */
typedef struct qaws_skeleton_face
{
	unsigned int const* nodes;
	unsigned int node_count;
	qaws_curve const* edge;         /* the input edge, owned by the skeleton; NULL for a splice */
	unsigned int path, curve, span; /* the input edge's source */
	qaws_scalar weight;
	int circle;
	qaws_vec2 normal;
	qaws_scalar offset;
	qaws_vec2 center;
	qaws_scalar radius;
	int side;
} qaws_skeleton_face;

qaws_status qaws_skeleton_create(qaws_skeleton_desc const* desc, qaws_skeleton** out_skeleton);
void qaws_skeleton_destroy(qaws_skeleton* skeleton);

unsigned int qaws_skeleton_get_node_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_node(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_node* out_node);
unsigned int qaws_skeleton_get_arc_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_arc(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_arc* out_arc);
unsigned int qaws_skeleton_get_face_count(qaws_skeleton const* skeleton);
qaws_status qaws_skeleton_get_face(qaws_skeleton const* skeleton, unsigned int index, qaws_skeleton_face* out_face);

/* The point of an arc at u in [0, 1] (node_a to node_b), exactly, and the
   time it is swept. */
qaws_status qaws_skeleton_evaluate_arc(qaws_skeleton const* skeleton, unsigned int index, qaws_scalar u,
	qaws_vec2* out_point, qaws_scalar* out_time);

/* The arc as a new curve (destroy it): a line, a circular arc, or a degree
   2 rational NURBS (the exact conic, in pieces when it turns far). */
qaws_status qaws_skeleton_get_arc_curve(qaws_skeleton const* skeleton, unsigned int index, qaws_curve** out_curve);

/* The time the region is swept out (the roof's top). */
qaws_scalar qaws_skeleton_get_max_time(qaws_skeleton const* skeleton);

/* The wavefront at time t >= 0 (the inward offset by t, mitered or round
   by the splice angle, weighted): its loops of lines and arcs, outer
   counter-clockwise, holes clockwise, as a qaws_clip result. Empty past
   the max time. */
qaws_status qaws_skeleton_offset(qaws_skeleton const* skeleton, qaws_scalar time, qaws_clip_result** out_result);

#endif /* QAWS_SKELETON_H */
