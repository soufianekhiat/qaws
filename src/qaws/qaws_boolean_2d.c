/*
 * qaws_boolean_2d: two closed curves through qaws_clip, each result path
 * returned as one curve (a composite when it has several pieces).
 */

#include "qaws_boolean_2d.h"
#include "qaws_clip.h"
#include "qaws_curve.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_composite.h"
#include <stdlib.h>
#include <string.h>

/* an owned copy of a curve: its whole exact piece */
static qaws_status bo_copy(qaws_curve const* c, qaws_curve** out)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	return qaws_curve_extract(c, r.min_value, r.max_value, out);
}

static qaws_status bo_path_curve(qaws_path_2d const* p, qaws_curve** out)
{
	qaws_curve** seg;
	qaws_composite_desc d;
	unsigned int i, made = 0;
	qaws_status s = QAWS_STATUS_OK;
	if (p->curve_count == 1)
		return bo_copy(p->curves[0], out);
	seg = (qaws_curve**)malloc(sizeof(qaws_curve*) * p->curve_count);
	if (!seg)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < p->curve_count && s == QAWS_STATUS_OK; i++)
		if ((s = bo_copy(p->curves[i], &seg[i])) == QAWS_STATUS_OK)
			made++;
	if (s == QAWS_STATUS_OK)
	{
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.segments = seg;
		d.segment_count = p->curve_count;
		s = qaws_curve_create_composite(&d, out);
	}
	if (s != QAWS_STATUS_OK)
		for (i = 0; i < made; i++)
			qaws_curve_destroy(seg[i]);
	free(seg);
	return s;
}

qaws_status qaws_boolean_2d(
	qaws_curve const* region_a,
	qaws_curve const* region_b,
	qaws_boolean_op operation,
	qaws_curve** out_boundaries,
	unsigned int boundary_capacity,
	unsigned int* out_count)
{
	qaws_path_2d pa, pb;
	qaws_clip_result* r = NULL;
	qaws_clip_type ct;
	qaws_status s;
	unsigned int i, n;

	if (!region_a || !region_b || !out_boundaries || !out_count || boundary_capacity == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(region_a) != QAWS_DIMENSION_2D ||
		qaws_curve_get_dimension(region_b) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;
	*out_count = 0;

	switch (operation)
	{
	case QAWS_BOOLEAN_UNION: ct = QAWS_CLIP_UNION; break;
	case QAWS_BOOLEAN_INTERSECTION: ct = QAWS_CLIP_INTERSECTION; break;
	case QAWS_BOOLEAN_DIFFERENCE: ct = QAWS_CLIP_DIFFERENCE; break;
	default: return QAWS_STATUS_INVALID_ARGUMENT;
	}

	pa.curves = &region_a; pa.curve_count = 1; pa.closed = 1;
	pb.curves = &region_b; pb.curve_count = 1; pb.closed = 1;
	s = qaws_clip_boolean(ct, QAWS_FILL_NON_ZERO, &pa, 1, &pb, 1, &r);
	if (s != QAWS_STATUS_OK)
		return s;
	n = qaws_clip_result_get_path_count(r);
	for (i = 0; i < n && i < boundary_capacity && s == QAWS_STATUS_OK; i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(r, i, &p);
		s = bo_path_curve(&p, &out_boundaries[i]);
		if (s == QAWS_STATUS_OK)
			(*out_count)++;
	}
	qaws_clip_result_destroy(r);
	if (s == QAWS_STATUS_OK && n > boundary_capacity)
		s = QAWS_STATUS_BUFFER_TOO_SMALL;
	return s;
}
