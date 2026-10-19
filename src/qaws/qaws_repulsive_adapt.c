/*
 * Adapters between qaws curves and the repulsive energies' working
 * polylines (qaws_repulsive.h): sampling with provenance, and the exact
 * degree-1 refit.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_bspline.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void ra_eval(qaws_curve const* c, unsigned int dim, double t, double* p)
{
	if (dim == 2)
	{
		qaws_eval_result_2d r;
		qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = 0;
	}
	else
	{
		qaws_eval_result_3d r;
		qaws_curve_evaluate_3d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = r.position.z;
	}
}

static double ra_dist(double const* a, double const* b)
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/* whether the curve's ends coincide (within 1e-9 of its extent, 1e-5 in f32) */
static int ra_closed(qaws_curve const* c, unsigned int dim)
{
	qaws_range rg = qaws_curve_get_parameter_range(c);
	double a[3], b[3], p[3], ext = 0;
	unsigned int k;
	ra_eval(c, dim, (double)rg.min_value, a);
	ra_eval(c, dim, (double)rg.max_value, b);
	for (k = 1; k < 8; k++)
	{
		ra_eval(c, dim, (double)rg.min_value + ((double)rg.max_value - rg.min_value) * k / 8, p);
		ext = fmax(ext, ra_dist(a, p));
	}
	return ext > 0 && ra_dist(a, b) <= (sizeof(qaws_scalar) == 4 ? 1e-5 : 1e-9) * ext;
}

typedef struct ra_out
{
	unsigned int dim, count;
	double h, tol;
	qaws_curve const* c;
	unsigned int cap, curve;
	qaws_scalar* pos;
	unsigned int* src;
	qaws_scalar* srct;
	unsigned int base;
} ra_out;

static void ra_emit(ra_out* o, double t, double const* p)
{
	unsigned int v = o->base + o->count, q;
	if (v < o->cap)
	{
		if (o->pos) for (q = 0; q < o->dim; q++) o->pos[o->dim * v + q] = (qaws_scalar)p[q];
		if (o->src) o->src[v] = o->curve;
		if (o->srct) o->srct[v] = (qaws_scalar)t;
	}
	o->count++;
}

/* [t0, t1] split until short and flat; emits the starts of the pieces */
static void ra_split(ra_out* o, double t0, double const* p0, double t1, double const* p1, unsigned int depth)
{
	double tm = 0.5 * (t0 + t1), pm[3], d[3], e = 0, l = ra_dist(p0, p1), u = 0;
	unsigned int q;
	ra_eval(o->c, o->dim, tm, pm);
	/* the middle's distance to the chord */
	for (q = 0; q < 3; q++) d[q] = p1[q] - p0[q];
	if (l > 0)
	{
		for (q = 0; q < 3; q++) u += (pm[q] - p0[q]) * d[q] / (l * l);
		for (q = 0; q < 3; q++) { double w = pm[q] - p0[q] - u * d[q]; e += w * w; }
		e = sqrt(e);
	}
	else
		e = ra_dist(pm, p0);
	if (depth < 40 && (l > o->h || e > o->tol))
	{
		ra_split(o, t0, p0, tm, pm, depth + 1);
		ra_split(o, tm, pm, t1, p1, depth + 1);
	}
	else
		ra_emit(o, t0, p0);
}

qaws_status qaws_repulsive_polyline_from_curves(qaws_repulsive_polyline_desc const* desc,
	qaws_scalar* out_positions, unsigned int vertex_capacity, unsigned int* out_edges, unsigned int edge_capacity,
	unsigned int* out_source_curve, qaws_scalar* out_source_t, unsigned int* out_vertex_count, unsigned int* out_edge_count)
{
	unsigned int i, nv = 0, ne = 0, k;
	if (!desc || !desc->curves || !desc->curve_count || (desc->dimension != 2 && desc->dimension != 3) || !(desc->edge_length > 0) ||
		!out_vertex_count || !out_edge_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->curve_count; i++)
		if (!desc->curves[i] || (unsigned int)qaws_curve_get_dimension(desc->curves[i]) != desc->dimension)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->curve_count; i++)
	{
		qaws_curve const* c = desc->curves[i];
		qaws_range rg = qaws_curve_get_parameter_range(c);
		double a = (double)rg.min_value, b = (double)rg.max_value, pa[3], pb[3], pm[3];
		int closed = ra_closed(c, desc->dimension);
		ra_out o;
		memset(&o, 0, sizeof(o));
		o.dim = desc->dimension; o.h = (double)desc->edge_length;
		o.tol = desc->tolerance > 0 ? (double)desc->tolerance : o.h / 100;
		o.c = c; o.cap = vertex_capacity; o.curve = i;
		o.pos = out_positions; o.src = out_source_curve; o.srct = out_source_t; o.base = nv;
		ra_eval(c, o.dim, a, pa);
		ra_eval(c, o.dim, b, pb);
		/* (a loop split once first: its chord is a point) */
		ra_eval(c, o.dim, 0.5 * (a + b), pm);
		ra_split(&o, a, pa, 0.5 * (a + b), pm, 0);
		ra_split(&o, 0.5 * (a + b), pm, b, pb, 0);
		if (!closed)
			ra_emit(&o, b, pb);
		for (k = 0; k + 1 < o.count + (closed ? 1u : 0u); k++, ne++)
			if (ne < edge_capacity && out_edges)
			{
				out_edges[2 * ne] = nv + k;
				out_edges[2 * ne + 1] = nv + (k + 1) % o.count;
			}
		nv += o.count;
	}
	*out_vertex_count = nv;
	*out_edge_count = ne;
	return nv > vertex_capacity || ne > edge_capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_polyline_to_curve(unsigned int dimension, qaws_scalar const* positions, unsigned int vertex_count,
	unsigned int const* source_curve, qaws_scalar const* source_t, unsigned int curve, qaws_curve const* source,
	qaws_curve** out_curve)
{
	unsigned int i, m = 0, q, first = 0;
	int closed;
	qaws_scalar *cp, *kn;
	qaws_bspline_desc d;
	qaws_status st;
	if (!positions || !source_curve || !source_t || !source || !out_curve || (dimension != 2 && dimension != 3) ||
		(unsigned int)qaws_curve_get_dimension(source) != dimension)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	closed = ra_closed(source, dimension);
	for (i = 0; i < vertex_count; i++)
		if (source_curve[i] == curve) { if (!m) first = i; m++; }
	if (m < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	cp = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), dimension, (m + 1)));
	kn = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (m + 3)));
	if (!cp || !kn)
	{
		qaws_internal_free(cp); qaws_internal_free(kn);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the vertices in their order, knots at their parameters (clamped) */
	m = 0;
	for (i = 0; i < vertex_count; i++)
		if (source_curve[i] == curve)
		{
			for (q = 0; q < dimension; q++) cp[dimension * m + q] = positions[dimension * i + q];
			kn[m + 1] = source_t[i];
			m++;
		}
	if (closed)
	{
		for (q = 0; q < dimension; q++) cp[dimension * m + q] = positions[dimension * first + q];
		kn[m + 1] = qaws_curve_get_parameter_range(source).max_value;
		m++;
	}
	kn[0] = kn[1];
	kn[m + 1] = kn[m];
	memset(&d, 0, sizeof(d));
	d.dimension = (qaws_dimension)dimension; d.degree = 1; d.control_points = cp; d.control_point_count = m;
	d.knots = kn; d.knot_count = m + 2;
	st = qaws_curve_create_bspline(&d, out_curve);
	qaws_internal_free(cp); qaws_internal_free(kn);
	return st;
}
