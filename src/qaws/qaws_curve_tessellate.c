/*
 * Adaptive curve tessellation (qaws_sampling.h): per span, depth-first halving on an explicit stack, left half first, so the
 * vertices come out in parameter order and the same every call. A segment is tested at 1/4, 1/2 and 3/4: the distance of
 * those points to its chord, the tangent's turning summed through them. Spans are evaluated from their own side, so a
 * kink at a span boundary does not count as turning.
 */

#include "qaws_sampling.h"
#include "internal/qaws_internal_types.h"
#include <math.h>
#include <string.h>

#define TC_DEPTH_DEFAULT 16
#define TC_DEPTH_MAX 30

typedef struct tc_point
{
	double u;                   /* span-local parameter */
	double p[3];
	double t[3];                /* unit tangent, 0 where D1 vanishes */
} tc_point;

static void tc_eval(qaws_curve const* c, unsigned int s, double u, tc_point* o)
{
	double l = 0;
	unsigned int k;
	memset(o, 0, sizeof *o);
	o->u = u;
	if (c->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d r;
		c->vtable->eval_span_2d(c, s, (qaws_scalar)u, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		o->p[0] = r.position.x; o->p[1] = r.position.y;
		o->t[0] = r.d1.x; o->t[1] = r.d1.y;
	}
	else
	{
		qaws_eval_result_3d r;
		c->vtable->eval_span_3d(c, s, (qaws_scalar)u, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		o->p[0] = r.position.x; o->p[1] = r.position.y; o->p[2] = r.position.z;
		o->t[0] = r.d1.x; o->t[1] = r.d1.y; o->t[2] = r.d1.z;
	}
	for (k = 0; k < 3; k++)
		l += o->t[k] * o->t[k];
	l = sqrt(l);
	for (k = 0; k < 3; k++)
		o->t[k] = l > 0 ? o->t[k] / l : 0;
}

/* the angle between two unit tangents; 0 when one is missing */
static double tc_turn(tc_point const* a, tc_point const* b)
{
	double d = a->t[0] * b->t[0] + a->t[1] * b->t[1] + a->t[2] * b->t[2];
	double x[3], s;
	if (a->t[0] == 0 && a->t[1] == 0 && a->t[2] == 0)
		return 0;
	if (b->t[0] == 0 && b->t[1] == 0 && b->t[2] == 0)
		return 0;
	x[0] = a->t[1] * b->t[2] - a->t[2] * b->t[1];
	x[1] = a->t[2] * b->t[0] - a->t[0] * b->t[2];
	x[2] = a->t[0] * b->t[1] - a->t[1] * b->t[0];
	s = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	return atan2(s, d);
}

/* the distance from x to the segment ab */
static double tc_chord(tc_point const* a, tc_point const* b, tc_point const* x)
{
	double ab[3], ax[3], l = 0, f = 0, d = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = b->p[k] - a->p[k];
		ax[k] = x->p[k] - a->p[k];
		l += ab[k] * ab[k];
		f += ab[k] * ax[k];
	}
	f = l > 0 ? f / l : 0;
	f = f < 0 ? 0 : f > 1 ? 1 : f;
	for (k = 0; k < 3; k++)
		d += (ax[k] - f * ab[k]) * (ax[k] - f * ab[k]);
	return sqrt(d);
}

static void tc_emit(qaws_curve const* c, unsigned int s, tc_point const* v, unsigned int n, qaws_scalar* params, void* points,
	unsigned int capacity)
{
	if (!params || n >= capacity)
		return;
	params[n] = (qaws_scalar)(c->span_boundaries[s] + v->u * ((double)c->span_boundaries[s + 1] - c->span_boundaries[s]));
	if (points && c->dimension == QAWS_DIMENSION_2D)
	{
		((qaws_vec2*)points)[n].x = (qaws_scalar)v->p[0];
		((qaws_vec2*)points)[n].y = (qaws_scalar)v->p[1];
	}
	else if (points)
	{
		((qaws_vec3*)points)[n].x = (qaws_scalar)v->p[0];
		((qaws_vec3*)points)[n].y = (qaws_scalar)v->p[1];
		((qaws_vec3*)points)[n].z = (qaws_scalar)v->p[2];
	}
}

static qaws_status tc_tessellate(qaws_curve const* c, qaws_dimension dim, qaws_curve_tessellation_desc const* desc,
	qaws_scalar* params, void* points, unsigned int capacity, unsigned int* out_count)
{
	struct { tc_point a, b; unsigned int depth; } stack[TC_DEPTH_MAX + 2];
	double chord, angle;
	unsigned int depth_max, n = 0, s, top;
	if (!c || !desc || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (c->dimension != dim)
		return QAWS_STATUS_INVALID_DIMENSION;
	chord = desc->chord_tolerance;
	angle = desc->angle_tolerance;
	if (!(chord >= 0) || !(angle >= 0) || !(chord > 0 || angle > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	depth_max = desc->max_depth ? desc->max_depth : TC_DEPTH_DEFAULT;
	depth_max = depth_max > TC_DEPTH_MAX ? TC_DEPTH_MAX : depth_max;
	for (s = 0; s < c->span_count; s++)
	{
		tc_eval(c, s, 0, &stack[0].a);
		tc_eval(c, s, 1, &stack[0].b);
		stack[0].depth = 0;
		if (s == 0)
			tc_emit(c, 0, &stack[0].a, n++, params, points, capacity);
		top = 1;
		while (top > 0)
		{
			tc_point a = stack[top - 1].a, b = stack[top - 1].b, q1, m = a, q3;
			unsigned int depth = stack[--top].depth;
			int ok = depth >= depth_max;
			if (!ok)
			{
				tc_eval(c, s, a.u + 0.25 * (b.u - a.u), &q1);
				tc_eval(c, s, 0.5 * (a.u + b.u), &m);
				tc_eval(c, s, a.u + 0.75 * (b.u - a.u), &q3);
				ok = (chord == 0 || (tc_chord(&a, &b, &q1) <= chord && tc_chord(&a, &b, &m) <= chord && tc_chord(&a, &b, &q3) <= chord))
					&& (angle == 0 || tc_turn(&a, &q1) + tc_turn(&q1, &m) + tc_turn(&m, &q3) + tc_turn(&q3, &b) <= angle);
			}
			if (ok)
			{
				tc_emit(c, s, &b, n++, params, points, capacity);
				continue;
			}
			stack[top].a = m; stack[top].b = b; stack[top].depth = depth + 1; top++;
			stack[top].a = a; stack[top].b = m; stack[top].depth = depth + 1; top++;
		}
	}
	*out_count = n;
	return params && n > capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

qaws_status qaws_curve_tessellate_2d(qaws_curve const* curve, qaws_curve_tessellation_desc const* desc, qaws_scalar* out_parameters,
	qaws_vec2* out_points, unsigned int capacity, unsigned int* out_count)
{
	return tc_tessellate(curve, QAWS_DIMENSION_2D, desc, out_parameters, out_points, capacity, out_count);
}

qaws_status qaws_curve_tessellate_3d(qaws_curve const* curve, qaws_curve_tessellation_desc const* desc, qaws_scalar* out_parameters,
	qaws_vec3* out_points, unsigned int capacity, unsigned int* out_count)
{
	return tc_tessellate(curve, QAWS_DIMENSION_3D, desc, out_parameters, out_points, capacity, out_count);
}
