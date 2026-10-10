#include "qaws_boolean_2d.h"
#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "qaws_bezier.h"
#include "qaws_composite.h"
#include "qaws_curve.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_fit.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Maximum intersection points between two curves */
#define BOOL_MAX_ISECT 256
#define BOOL_SAMPLES_PER_ARC 48

/* ========================================================================== */
/*  Internal: descriptor for a kept arc (not yet sampled)                     */
/* ========================================================================== */

typedef struct {
	qaws_curve const* source;
	qaws_scalar t_start;
	qaws_scalar t_end;
	/* cached geometric endpoints */
	qaws_vec2 geo_start;
	qaws_vec2 geo_end;
} bool_arc;

static qaws_vec2 eval_pt(qaws_curve const* c, qaws_scalar t)
{
	qaws_eval_result_2d er;
	qaws_vec2 r = {QAWS_ZERO, QAWS_ZERO};
	memset(&er, 0, sizeof(er));
	qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &er);
	r = er.position;
	return r;
}

static qaws_scalar dist2_vec2(qaws_vec2 a, qaws_vec2 b)
{
	qaws_scalar dx = a.x - b.x;
	qaws_scalar dy = a.y - b.y;
	return dx * dx + dy * dy;
}

static void arc_cache_endpoints(bool_arc* a)
{
	a->geo_start = eval_pt(a->source, a->t_start);
	a->geo_end   = eval_pt(a->source, a->t_end);
}

/* ========================================================================== */
/*  Internal: sample a sub-arc of a curve between two parameters              */
/* ========================================================================== */

static qaws_status sample_arc_as_bspline(
	qaws_curve const* curve,
	qaws_scalar t_start,
	qaws_scalar t_end,
	unsigned int n_samples,
	qaws_curve** out_curve)
{
	qaws_scalar* params = NULL;
	qaws_scalar* coords = NULL;
	unsigned int i;
	unsigned int n_cp;
	qaws_status status;

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples);
	coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)n_samples * 2);
	if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

	for (i = 0; i < n_samples; i++)
	{
		qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)(n_samples - 1);
		qaws_scalar t = t_start + frac * (t_end - t_start);
		qaws_eval_result_2d er;

		memset(&er, 0, sizeof(er));
		qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
		params[i] = frac;
		coords[i * 2 + 0] = er.position.x;
		coords[i * 2 + 1] = er.position.y;
	}

	n_cp = n_samples / 4;
	if (n_cp < 4) n_cp = 4;
	if (n_cp > n_samples) n_cp = n_samples;

	status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, params, coords, n_samples, n_cp, out_curve);
	free(params);
	free(coords);
	return status;
}

/* ========================================================================== */
/*  Internal: evaluate midpoint of a segment for winding classification       */
/* ========================================================================== */

static qaws_vec2 segment_midpoint(
	qaws_curve const* curve,
	qaws_scalar t_start,
	qaws_scalar t_end)
{
	qaws_eval_result_2d er;
	qaws_scalar t_mid = (t_start + t_end) * QAWS_LITERAL(0.5);
	qaws_vec2 result = {QAWS_ZERO, QAWS_ZERO};

	memset(&er, 0, sizeof(er));
	qaws_curve_evaluate_2d(curve, t_mid, QAWS_EVAL_FLAG_POSITION, &er);
	result = er.position;
	return result;
}

/* ========================================================================== */
/*  Internal: check if a point is "inside" a closed curve (winding != 0)      */
/* ========================================================================== */

static int point_inside_curve(qaws_curve const* curve, qaws_vec2 point)
{
	int winding = 0;
	qaws_curve_compute_winding_number_2d(curve, point, &winding);
	return winding != 0;
}

/* ========================================================================== */
/*  Internal: classify a segment as kept or discarded for a boolean op        */
/* ========================================================================== */

static int segment_kept(
	qaws_boolean_op op,
	int from_a,  /* 1 if segment from curve A, 0 if from curve B */
	int inside_other)  /* 1 if midpoint is inside the other curve */
{
	switch (op)
	{
	case QAWS_BOOLEAN_UNION:
		/* Keep segments that are NOT inside the other curve */
		return !inside_other;

	case QAWS_BOOLEAN_INTERSECTION:
		/* Keep segments that ARE inside the other curve */
		return inside_other;

	case QAWS_BOOLEAN_DIFFERENCE:
		if (from_a)
			return !inside_other; /* A segments: keep those outside B */
		else
			return inside_other;  /* B segments: keep those inside A (reversed) */

	default:
		return 0;
	}
}

/* ========================================================================== */
/*  Internal: sample arcs into one continuous polyline and fit B-spline       */
/* ========================================================================== */

static qaws_status arcs_to_bspline(
	bool_arc const* arcs,
	unsigned int arc_count,
	unsigned int samples_per_arc,
	qaws_curve** out_curve)
{
	unsigned int total_samples = arc_count * samples_per_arc;
	qaws_scalar* params = NULL;
	qaws_scalar* coords = NULL;
	unsigned int idx = 0;
	unsigned int ai, si;
	unsigned int n_cp;
	qaws_status status;

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)total_samples);
	coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)total_samples * 2);
	if (!params || !coords) { free(params); free(coords); return QAWS_STATUS_ALLOCATION_FAILURE; }

	for (ai = 0; ai < arc_count; ai++)
	{
		/* Skip duplicate point at junction (except first arc) */
		unsigned int start_si = (ai == 0) ? 0 : 1;
		for (si = start_si; si < samples_per_arc; si++)
		{
			qaws_scalar frac = (qaws_scalar)si / (qaws_scalar)(samples_per_arc - 1);
			qaws_scalar t = arcs[ai].t_start + frac * (arcs[ai].t_end - arcs[ai].t_start);
			qaws_eval_result_2d er;

			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(arcs[ai].source, t, QAWS_EVAL_FLAG_POSITION, &er);
			params[idx] = (qaws_scalar)idx;
			coords[idx * 2 + 0] = er.position.x;
			coords[idx * 2 + 1] = er.position.y;
			idx++;
		}
	}

	/* Normalize params to [0,1] */
	if (idx > 1)
	{
		unsigned int pi;
		for (pi = 0; pi < idx; pi++)
			params[pi] = (qaws_scalar)pi / (qaws_scalar)(idx - 1);
	}

	n_cp = idx / 3;
	if (n_cp < 4) n_cp = 4;
	if (n_cp > idx) n_cp = idx;

	status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3, params, coords, idx, n_cp, out_curve);
	free(params);
	free(coords);
	return status;
}

/* ========================================================================== */
/*  Internal: reorder arcs into a connected chain by endpoint matching        */
/* ========================================================================== */

static void reorder_arcs(bool_arc* arcs, unsigned int count)
{
	unsigned int i, j;
	qaws_scalar best_dist;
	unsigned int best_j;
	int best_reversed;

	for (i = 0; i < count - 1; i++)
	{
		qaws_vec2 cur_end = arcs[i].geo_end;
		best_dist = QAWS_LITERAL(1e30);
		best_j = i + 1;
		best_reversed = 0;

		for (j = i + 1; j < count; j++)
		{
			qaws_scalar d_start = dist2_vec2(cur_end, arcs[j].geo_start);
			qaws_scalar d_end   = dist2_vec2(cur_end, arcs[j].geo_end);

			if (d_start < best_dist)
			{
				best_dist = d_start;
				best_j = j;
				best_reversed = 0;
			}
			if (d_end < best_dist)
			{
				best_dist = d_end;
				best_j = j;
				best_reversed = 1;
			}
		}

		/* Swap best into position i+1 */
		if (best_j != i + 1)
		{
			bool_arc tmp = arcs[i + 1];
			arcs[i + 1] = arcs[best_j];
			arcs[best_j] = tmp;
		}

		/* Reverse if needed (swap t_start/t_end and geo endpoints) */
		if (best_reversed)
		{
			qaws_scalar tmp_t = arcs[i + 1].t_start;
			qaws_vec2 tmp_v = arcs[i + 1].geo_start;
			arcs[i + 1].t_start = arcs[i + 1].t_end;
			arcs[i + 1].t_end = tmp_t;
			arcs[i + 1].geo_start = arcs[i + 1].geo_end;
			arcs[i + 1].geo_end = tmp_v;
		}
	}
}

/* ========================================================================== */
/*  Internal: merge first and last kept arcs from the same curve              */
/*  (wrap-around for closed curves where domain start ≈ domain end)          */
/* ========================================================================== */

static void merge_wraparound_arcs(
	bool_arc* arcs, unsigned int* count,
	qaws_curve const* curve, qaws_range range)
{
	unsigned int n = *count;
	unsigned int first_idx = (unsigned int)-1;
	unsigned int last_idx  = (unsigned int)-1;
	unsigned int i;
	qaws_scalar domain_tol;
	bool_arc merged;

	if (n < 2) return;

	domain_tol = (range.max_value - range.min_value) * QAWS_LITERAL(0.001);

	/* Find first and last arcs from this curve */
	for (i = 0; i < n; i++)
	{
		if (arcs[i].source == curve)
		{
			if (first_idx == (unsigned int)-1) first_idx = i;
			last_idx = i;
		}
	}

	if (first_idx == last_idx || first_idx == (unsigned int)-1) return;

	/* Check if first arc starts near domain min and last arc ends near domain max */
	if (QAWS_FABS(arcs[first_idx].t_start - range.min_value) < domain_tol &&
		QAWS_FABS(arcs[last_idx].t_end - range.max_value) < domain_tol)
	{
		/* Merge: last arc's range + first arc's range form a single wrap-around arc.
		   We'll sample last arc then first arc as one continuous piece. */
		/* Replace first_idx with merged, remove last_idx */
		merged.source = curve;
		merged.t_start = arcs[last_idx].t_start;
		merged.t_end   = arcs[first_idx].t_end;
		/* For the merged arc, we need special handling since t_start > t_end
		   (wraps through domain boundary). Mark with a flag by using the
		   domain max as a breakpoint. We'll handle this by storing t_end
		   as t_end + (domain_max - domain_min) to indicate wrap-around. */
		/* Actually, simpler: geo endpoints of merged arc */
		merged.geo_start = arcs[last_idx].geo_start;
		merged.geo_end   = arcs[first_idx].geo_end;

		arcs[first_idx] = merged;

		/* Remove last_idx */
		for (i = last_idx; i + 1 < n; i++)
			arcs[i] = arcs[i + 1];
		(*count)--;
	}
}

/* ========================================================================== */
/*  Internal: sample a possibly wrap-around arc                               */
/* ========================================================================== */

static unsigned int sample_arc_points(
	bool_arc const* arc,
	qaws_range range,
	qaws_scalar* coords,  /* output: x,y pairs */
	unsigned int max_samples)
{
	unsigned int idx = 0;
	unsigned int i;

	if (arc->t_start <= arc->t_end)
	{
		/* Normal arc */
		for (i = 0; i < max_samples && idx < max_samples; i++)
		{
			qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)(max_samples - 1);
			qaws_scalar t = arc->t_start + frac * (arc->t_end - arc->t_start);
			qaws_eval_result_2d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(arc->source, t, QAWS_EVAL_FLAG_POSITION, &er);
			coords[idx * 2 + 0] = er.position.x;
			coords[idx * 2 + 1] = er.position.y;
			idx++;
		}
	}
	else
	{
		/* Wrap-around: sample from t_start to range.max, then range.min to t_end */
		unsigned int half = max_samples / 2;
		qaws_scalar span1 = range.max_value - arc->t_start;
		qaws_scalar span2 = arc->t_end - range.min_value;
		qaws_scalar total_span = span1 + span2;
		unsigned int n1, n2;

		if (total_span < QAWS_LITERAL(1e-10)) return 0;
		n1 = (unsigned int)(half * (span1 / total_span));
		if (n1 < 2) n1 = 2;
		n2 = max_samples - n1;
		if (n2 < 2) n2 = 2;

		/* First part: t_start → range.max */
		for (i = 0; i < n1 && idx < max_samples; i++)
		{
			qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)(n1 - 1);
			qaws_scalar t = arc->t_start + frac * (range.max_value - arc->t_start);
			qaws_eval_result_2d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(arc->source, t, QAWS_EVAL_FLAG_POSITION, &er);
			coords[idx * 2 + 0] = er.position.x;
			coords[idx * 2 + 1] = er.position.y;
			idx++;
		}
		/* Second part: range.min → t_end (skip first to avoid duplicate) */
		for (i = 1; i < n2 && idx < max_samples; i++)
		{
			qaws_scalar frac = (qaws_scalar)i / (qaws_scalar)(n2 - 1);
			qaws_scalar t = range.min_value + frac * (arc->t_end - range.min_value);
			qaws_eval_result_2d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_2d(arc->source, t, QAWS_EVAL_FLAG_POSITION, &er);
			coords[idx * 2 + 0] = er.position.x;
			coords[idx * 2 + 1] = er.position.y;
			idx++;
		}
	}
	return idx;
}

/* ========================================================================== */
/*  Internal: sample ordered arcs into a composite curve (preserves corners)  */
/* ========================================================================== */

static qaws_status ordered_arcs_to_composite(
	bool_arc const* arcs,
	unsigned int arc_count,
	qaws_range range_a,
	qaws_range range_b,
	unsigned int samples_per_arc,
	qaws_curve** out_curve)
{
	qaws_curve** segments = NULL;
	unsigned int seg_count = 0;
	unsigned int ai;
	qaws_status status = QAWS_STATUS_OK;

	segments = (qaws_curve**)calloc((size_t)arc_count, sizeof(qaws_curve*));
	if (!segments) return QAWS_STATUS_ALLOCATION_FAILURE;

	for (ai = 0; ai < arc_count; ai++)
	{
		qaws_range range = (arcs[ai].source == arcs[0].source) ? range_a : range_b;
		qaws_scalar* coords = NULL;
		qaws_scalar* params = NULL;
		unsigned int n_pts, n_cp, pi;

		coords = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)samples_per_arc * 2);
		params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)samples_per_arc);
		if (!coords || !params)
		{
			free(coords); free(params);
			status = QAWS_STATUS_ALLOCATION_FAILURE;
			goto cleanup_composite;
		}

		n_pts = sample_arc_points(&arcs[ai], range, coords, samples_per_arc);
		if (n_pts < 4)
		{
			free(coords); free(params);
			status = QAWS_STATUS_DEGENERATE_CURVE;
			goto cleanup_composite;
		}

		for (pi = 0; pi < n_pts; pi++)
			params[pi] = (qaws_scalar)pi / (qaws_scalar)(n_pts - 1);

		n_cp = n_pts / 4;
		if (n_cp < 4) n_cp = 4;
		if (n_cp > n_pts) n_cp = n_pts;

		status = qaws_internal_fit_bspline(QAWS_DIMENSION_2D, 3,
			params, coords, n_pts, n_cp, &segments[seg_count]);
		free(params); free(coords);

		if (status != QAWS_STATUS_OK) goto cleanup_composite;
		seg_count++;
	}

	if (seg_count == 1)
	{
		*out_curve = segments[0];
		free(segments);
		return QAWS_STATUS_OK;
	}

	/* Create composite curve - sharp corners are preserved at junctions */
	{
		qaws_composite_desc cdesc;
		memset(&cdesc, 0, sizeof(cdesc));
		cdesc.dimension = QAWS_DIMENSION_2D;
		cdesc.segments = segments;
		cdesc.segment_count = seg_count;
		status = qaws_curve_create_composite(&cdesc, out_curve);
	}

	if (status != QAWS_STATUS_OK) goto cleanup_composite;
	free(segments);
	return QAWS_STATUS_OK;

cleanup_composite:
	{
		unsigned int ci;
		for (ci = 0; ci < seg_count; ci++)
			qaws_curve_destroy(segments[ci]);
		free(segments);
	}
	return status;
}

/* ========================================================================== */
/*  Public API                                                                */
/* ========================================================================== */

qaws_status qaws_boolean_2d(
	qaws_curve const* region_a,
	qaws_curve const* region_b,
	qaws_boolean_op operation,
	qaws_curve** out_boundaries,
	unsigned int boundary_capacity,
	unsigned int* out_count)
{
	qaws_intersection_2d isects[BOOL_MAX_ISECT];
	unsigned int isect_count = 0;
	qaws_status status;
	qaws_range range_a, range_b;
	qaws_scalar* t_a = NULL;
	qaws_scalar* t_b = NULL;
	bool_arc* arcs = NULL;
	unsigned int arc_count = 0;
	unsigned int max_arcs;
	unsigned int i, j;

	if (!region_a || !region_b || !out_boundaries || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (boundary_capacity == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_curve_get_dimension(region_a) != QAWS_DIMENSION_2D ||
		qaws_curve_get_dimension(region_b) != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	*out_count = 0;

	/* Find intersections between the two boundaries */
	status = qaws_curve_find_intersections_2d(
		region_a, region_b, isects, BOOL_MAX_ISECT, &isect_count);
	if (status != QAWS_STATUS_OK)
		return status;
	/* the count is the total found; only the buffer's worth was written */
	if (isect_count > BOOL_MAX_ISECT)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	range_a = qaws_curve_get_parameter_range(region_a);
	range_b = qaws_curve_get_parameter_range(region_b);

	/* No intersections: one may be entirely inside the other */
	if (isect_count == 0)
	{
		qaws_eval_result_2d er_a, er_b;
		int a_in_b, b_in_a;

		memset(&er_a, 0, sizeof(er_a));
		memset(&er_b, 0, sizeof(er_b));
		qaws_curve_evaluate_2d(region_a, range_a.min_value, QAWS_EVAL_FLAG_POSITION, &er_a);
		qaws_curve_evaluate_2d(region_b, range_b.min_value, QAWS_EVAL_FLAG_POSITION, &er_b);

		a_in_b = point_inside_curve(region_b, er_a.position);
		b_in_a = point_inside_curve(region_a, er_b.position);

		switch (operation)
		{
		case QAWS_BOOLEAN_UNION:
			if (a_in_b)
			{
				status = sample_arc_as_bspline(region_b, range_b.min_value, range_b.max_value, 128, &out_boundaries[0]);
				if (status == QAWS_STATUS_OK) *out_count = 1;
				return status;
			}
			else if (b_in_a)
			{
				status = sample_arc_as_bspline(region_a, range_a.min_value, range_a.max_value, 128, &out_boundaries[0]);
				if (status == QAWS_STATUS_OK) *out_count = 1;
				return status;
			}
			else
			{
				status = sample_arc_as_bspline(region_a, range_a.min_value, range_a.max_value, 128, &out_boundaries[0]);
				if (status != QAWS_STATUS_OK) return status;
				if (boundary_capacity > 1)
				{
					status = sample_arc_as_bspline(region_b, range_b.min_value, range_b.max_value, 128, &out_boundaries[1]);
					if (status == QAWS_STATUS_OK) *out_count = 2;
					else { qaws_curve_destroy(out_boundaries[0]); *out_count = 0; return status; }
				}
				else *out_count = 1;
				return QAWS_STATUS_OK;
			}

		case QAWS_BOOLEAN_INTERSECTION:
			if (a_in_b)
			{
				status = sample_arc_as_bspline(region_a, range_a.min_value, range_a.max_value, 128, &out_boundaries[0]);
				if (status == QAWS_STATUS_OK) *out_count = 1;
				return status;
			}
			else if (b_in_a)
			{
				status = sample_arc_as_bspline(region_b, range_b.min_value, range_b.max_value, 128, &out_boundaries[0]);
				if (status == QAWS_STATUS_OK) *out_count = 1;
				return status;
			}
			else
			{
				*out_count = 0;
				return QAWS_STATUS_OK;
			}

		case QAWS_BOOLEAN_DIFFERENCE:
			if (a_in_b)
			{
				*out_count = 0;
				return QAWS_STATUS_OK;
			}
			else if (b_in_a)
			{
				status = sample_arc_as_bspline(region_a, range_a.min_value, range_a.max_value, 128, &out_boundaries[0]);
				if (status != QAWS_STATUS_OK) return status;
				if (boundary_capacity > 1)
				{
					status = sample_arc_as_bspline(region_b, range_b.min_value, range_b.max_value, 128, &out_boundaries[1]);
					if (status == QAWS_STATUS_OK) *out_count = 2;
					else { qaws_curve_destroy(out_boundaries[0]); *out_count = 0; return status; }
				}
				else *out_count = 1;
				return QAWS_STATUS_OK;
			}
			else
			{
				status = sample_arc_as_bspline(region_a, range_a.min_value, range_a.max_value, 128, &out_boundaries[0]);
				if (status == QAWS_STATUS_OK) *out_count = 1;
				return status;
			}
		}
		return QAWS_STATUS_OK;
	}

	/* Sort intersection parameters on each curve */
	t_a = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(isect_count + 2));
	t_b = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(isect_count + 2));
	if (!t_a || !t_b) { free(t_a); free(t_b); return QAWS_STATUS_ALLOCATION_FAILURE; }

	/* Include domain endpoints */
	t_a[0] = range_a.min_value;
	t_b[0] = range_b.min_value;
	for (i = 0; i < isect_count; i++)
	{
		t_a[i + 1] = isects[i].parameter_a;
		t_b[i + 1] = isects[i].parameter_b;
	}
	t_a[isect_count + 1] = range_a.max_value;
	t_b[isect_count + 1] = range_b.max_value;

	/* Simple insertion sort for t_a */
	for (i = 1; i < isect_count + 2; i++)
	{
		qaws_scalar key = t_a[i];
		j = i;
		while (j > 0 && t_a[j - 1] > key)
		{
			t_a[j] = t_a[j - 1];
			j--;
		}
		t_a[j] = key;
	}
	/* Sort t_b */
	for (i = 1; i < isect_count + 2; i++)
	{
		qaws_scalar key = t_b[i];
		j = i;
		while (j > 0 && t_b[j - 1] > key)
		{
			t_b[j] = t_b[j - 1];
			j--;
		}
		t_b[j] = key;
	}

	/* Collect kept arc descriptors (not yet sampled) */
	max_arcs = (isect_count + 1) * 2;
	arcs = (bool_arc*)calloc((size_t)max_arcs, sizeof(bool_arc));
	if (!arcs) { free(t_a); free(t_b); return QAWS_STATUS_ALLOCATION_FAILURE; }

	/* Process segments from curve A */
	for (i = 0; i + 1 < isect_count + 2; i++)
	{
		qaws_scalar ta0 = t_a[i];
		qaws_scalar ta1 = t_a[i + 1];
		qaws_vec2 mid;
		int inside;

		if (ta1 - ta0 < QAWS_LITERAL(1e-10)) continue;

		mid = segment_midpoint(region_a, ta0, ta1);
		inside = point_inside_curve(region_b, mid);

		if (segment_kept(operation, 1, inside) && arc_count < max_arcs)
		{
			arcs[arc_count].source = region_a;
			arcs[arc_count].t_start = ta0;
			arcs[arc_count].t_end = ta1;
			arc_cache_endpoints(&arcs[arc_count]);
			arc_count++;
		}
	}

	/* Process segments from curve B */
	for (i = 0; i + 1 < isect_count + 2; i++)
	{
		qaws_scalar tb0 = t_b[i];
		qaws_scalar tb1 = t_b[i + 1];
		qaws_vec2 mid;
		int inside;

		if (tb1 - tb0 < QAWS_LITERAL(1e-10)) continue;

		mid = segment_midpoint(region_b, tb0, tb1);
		inside = point_inside_curve(region_a, mid);

		if (segment_kept(operation, 0, inside) && arc_count < max_arcs)
		{
			arcs[arc_count].source = region_b;
			arcs[arc_count].t_start = tb0;
			arcs[arc_count].t_end = tb1;
			arc_cache_endpoints(&arcs[arc_count]);
			arc_count++;
		}
	}

	free(t_a);
	free(t_b);

	if (arc_count == 0)
	{
		free(arcs);
		*out_count = 0;
		return QAWS_STATUS_OK;
	}

	/* Merge wrap-around segments for closed curves */
	merge_wraparound_arcs(arcs, &arc_count, region_a, range_a);
	merge_wraparound_arcs(arcs, &arc_count, region_b, range_b);

	/* Reorder arcs to form a connected chain */
	if (arc_count > 1)
		reorder_arcs(arcs, arc_count);

	/* Build composite from individual arc B-splines (preserves sharp corners) */
	status = ordered_arcs_to_composite(arcs, arc_count, range_a, range_b,
		BOOL_SAMPLES_PER_ARC, &out_boundaries[0]);

	free(arcs);

	if (status == QAWS_STATUS_OK)
		*out_count = 1;
	return status;
}
