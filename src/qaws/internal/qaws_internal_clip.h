#ifndef QAWS_INTERNAL_CLIP_H
#define QAWS_INTERNAL_CLIP_H

#include "../qaws_clip.h"
#include "../qaws_clip64.h"

/* Moves everything of src (paths, open paths, curves, vertices) to the end
   of dst and destroys src; parents keep pointing at their own paths. */
qaws_status qaws_internal_clip_result_append(qaws_clip_result* dst, qaws_clip_result* src);
qaws_status qaws_internal_clip64_result_append(qaws_clip64_result* dst, qaws_clip64_result* src);

/* An empty result. */
qaws_status qaws_internal_clip_result_empty(qaws_clip_result** out);
qaws_status qaws_internal_clip64_result_empty(qaws_clip64_result** out);

/*
 * An arrangement decided elsewhere (the exact engine), turned into a
 * qaws_clip_result by qaws_clip's output stage: loops of the half-edges
 * between inside and outside faces, pieces of one source fused, nesting from
 * the faces.
 *
 * Half-edge 2 e runs along edge e from t0 to t1 (from v0 to v1), 2 e + 1
 * back; its face is the one on its left. Each vertex lists its outgoing
 * half-edges counter-clockwise in out[first .. first + count).
 */
typedef struct qaws_internal_arr_source
{
	qaws_curve const* curve;
	unsigned int operand, path, index;   /* where it comes from (operand 0 subject, 1 clip) */
	double base, scale;                  /* the input curve's parameter = base + scale t */
} qaws_internal_arr_source;

typedef struct qaws_internal_arr_vertex
{
	double p[2];
	unsigned int sa, sb;                 /* sources meeting there (sb QAWS_CLIP_NONE_INDEX: one) */
	double ta, tb;
	unsigned int first, count;           /* its outgoing half-edges in out */
} qaws_internal_arr_vertex;

typedef struct qaws_internal_arr_edge
{
	unsigned int src;
	double t0, t1;                       /* along src, t0 < t1 */
	unsigned int v0, v1;
	int ds, dc;                          /* subject / clip winding steps from its right to its left */
} qaws_internal_arr_edge;

typedef struct qaws_internal_arrangement
{
	qaws_clip_desc const* desc;          /* flags, z callback */
	qaws_internal_arr_source const* sources;
	unsigned int source_count;
	qaws_internal_arr_vertex const* vertices;
	unsigned int vertex_count;
	qaws_internal_arr_edge const* edges;
	unsigned int edge_count;
	unsigned int const* out;             /* outgoing half-edges per vertex, counter-clockwise */
	unsigned int const* half_face;       /* per half-edge */
	unsigned int face_count;
	int const* inside;                   /* per face */
	double tolerance, extent;            /* collinear corners dropped within tolerance; area accuracy */
} qaws_internal_arrangement;

qaws_status qaws_internal_clip_output_arrangement(qaws_internal_arrangement const* a, qaws_clip_result* result);

#endif /* QAWS_INTERNAL_CLIP_H */
