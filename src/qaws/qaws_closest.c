/*
 * Certified closest points on curves with Bezier spans (internal/qaws_internal_closest.h): branch and bound over each span's
 * Bezier form, from its B-spline window by blossoming. A piece is dropped when a lower bound of its distance is no less
 * than the best distance found: the distance to the box of its control points, or to its fat line (the chord's segment
 * widened by the points' largest offset from it, second order in the piece's size, so only a few pieces survive near the
 * minimum). Pieces are halved by de Casteljau until their box is under the tolerance; a Newton polish ends it. Surfaces
 * the same over Bezier patches, with the box in each patch frame in place of the fat line.
 */

#include "internal/qaws_internal_closest.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_surface.h"
#include "qaws_surface.h"
#include "qaws_eval.h"
#include "qaws_platform.h"
#include <math.h>
#include <string.h>

#define CL_MAX (QAWS_INTERNAL_BASIS_MAX_DEGREE + 1)
#define CL_DEPTH 56
#define CL_WORK 4096

/* span s as its B-spline window and knots (see qaws_internal_view_span_window); 0 for families without Bezier spans */
static int cl_window(qaws_curve const* c, unsigned int s, unsigned int* p, double* win, double* uk)
{
	unsigned int dim = (unsigned int)c->dimension, h = dim + 1, d = c->degree, j, i;
	qaws_scalar const* pts = NULL;
	qaws_scalar const* w = NULL;
	if (c->kind == QAWS_CURVE_KIND_VIEW)
		return qaws_internal_view_span_window(c, s, p, win, uk);
	if (c->kind == QAWS_CURVE_KIND_BEZIER || c->kind == QAWS_CURVE_KIND_RATIONAL_BEZIER)
	{
		if (c->kind == QAWS_CURVE_KIND_BEZIER)
			pts = ((qaws_bezier_impl const*)c->impl)->control_points;
		for (j = 0; j <= d; j++)
			for (i = 0; i < h; i++)
				win[j * h + i] = pts ? (i < dim ? (double)pts[j * dim + i] : 1.0)
					: (double)((qaws_rational_bezier_impl const*)c->impl)->weighted_points[j * h + i];
		for (j = 0; j < 2 * d; j++)
			uk[j] = j < d ? 0.0 : 1.0;
	}
	else if (c->kind == QAWS_CURVE_KIND_BSPLINE || c->kind == QAWS_CURVE_KIND_NURBS)
	{
		qaws_scalar const* knots;
		unsigned int kc, n, k;
		if (c->kind == QAWS_CURVE_KIND_BSPLINE)
		{
			qaws_bspline_impl const* b = (qaws_bspline_impl const*)c->impl;
			pts = b->control_points; knots = b->knots; kc = b->knot_count; n = b->control_point_count;
		}
		else
		{
			qaws_nurbs_impl const* b = (qaws_nurbs_impl const*)c->impl;
			pts = b->control_points; knots = b->knots; kc = b->knot_count; n = b->control_point_count; w = b->weights;
		}
		k = qaws_internal_find_knot_span(knots, kc, d, n, (c->span_boundaries[s] + c->span_boundaries[s + 1]) / 2);
		for (j = 0; j <= d; j++)
		{
			double wj = w ? (double)w[k - d + j] : 1.0;
			for (i = 0; i < dim; i++)
				win[j * h + i] = wj * (double)pts[(k - d + j) * dim + i];
			win[j * h + dim] = wj;
		}
		for (j = 0; j < 2 * d; j++)
			uk[j] = (double)knots[k - d + 1 + j];
	}
	else
		return 0;
	*p = d;
	return 1;
}

/* Bezier point i of the span is the blossom at (u_k ^ (p - i), u_k+1 ^ i): de Boor with those arguments (uk[p - 1] = u_k) */
static void cl_bezier(double const* win, double const* uk, unsigned int p, unsigned int h, double* bez)
{
	double tmp[CL_MAX * 4];
	unsigned int i, r, j, c;
	for (i = 0; i <= p; i++)
	{
		memcpy(tmp, win, sizeof(double) * (p + 1) * h);
		for (r = 1; r <= p; r++)
		{
			double x = r <= p - i ? uk[p - 1] : uk[p];
			for (j = p; j >= r; j--)
			{
				double a = uk[j - 1], al = (x - a) / (uk[j + p - r] - a);
				for (c = 0; c < h; c++)
					tmp[j * h + c] = (1 - al) * tmp[(j - 1) * h + c] + al * tmp[j * h + c];
			}
		}
		memcpy(bez + i * h, tmp + p * h, sizeof(double) * h);
	}
}

/* de Casteljau at u in place: keeps the right part (or the left) */
static void cl_split(double* b, unsigned int p, unsigned int h, double u, int right)
{
	unsigned int r, j, c;
	for (r = 1; r <= p; r++)
	{
		if (right)
			for (j = 0; j + r <= p; j++)
				for (c = 0; c < h; c++)
					b[j * h + c] += u * (b[(j + 1) * h + c] - b[j * h + c]);
		else
			for (j = p; j >= r; j--)
				for (c = 0; c < h; c++)
					b[j * h + c] = b[(j - 1) * h + c] + u * (b[j * h + c] - b[(j - 1) * h + c]);
	}
}

static double cl_dist(double const* a, double const* b, unsigned int dim)
{
	double s = 0;
	unsigned int c;
	for (c = 0; c < dim; c++)
		s += (a[c] - b[c]) * (a[c] - b[c]);
	return sqrt(s);
}

/* the piece's points projected (q[i]); its lower bound of the distance to x, its box diagonal */
static double cl_lower(double const* qp, unsigned int p, unsigned int dim, double const* x, double* diag)
{
	double lo[3], hi[3], box = 0, dg = 0, u[3], len = 0, smin = 0, smax = 0, off = 0, sx = 0, px = 0, fat;
	unsigned int i, c;
	for (c = 0; c < dim; c++)
	{
		lo[c] = hi[c] = qp[c];
		for (i = 1; i <= p; i++)
		{
			lo[c] = qp[i * 3 + c] < lo[c] ? qp[i * 3 + c] : lo[c];
			hi[c] = qp[i * 3 + c] > hi[c] ? qp[i * 3 + c] : hi[c];
		}
		box += x[c] < lo[c] ? (lo[c] - x[c]) * (lo[c] - x[c]) : x[c] > hi[c] ? (x[c] - hi[c]) * (x[c] - hi[c]) : 0;
		dg += (hi[c] - lo[c]) * (hi[c] - lo[c]);
		u[c] = qp[p * 3 + c] - qp[c];
		len += u[c] * u[c];
	}
	*diag = sqrt(dg);
	box = sqrt(box);
	len = sqrt(len);
	if (!(len > 1e-300))
		return box;
	for (c = 0; c < dim; c++)
		u[c] /= len;
	/* every point within off of the chord line, between smin and smax along it */
	for (i = 0; i <= p; i++)
	{
		double s = 0, r = 0;
		for (c = 0; c < dim; c++)
			s += (qp[i * 3 + c] - qp[c]) * u[c];
		for (c = 0; c < dim; c++)
		{
			double e = qp[i * 3 + c] - qp[c] - s * u[c];
			r += e * e;
		}
		smin = s < smin ? s : smin;
		smax = s > smax ? s : smax;
		off = sqrt(r) > off ? sqrt(r) : off;
	}
	for (c = 0; c < dim; c++)
		sx += (x[c] - qp[c]) * u[c];
	for (c = 0; c < dim; c++)
	{
		double e = x[c] - qp[c] - sx * u[c];
		px += e * e;
	}
	sx = sx < smin ? smin - sx : sx > smax ? sx - smax : 0;
	fat = sqrt(sx * sx + px) - off;
	return fat > box ? fat : box;
}

static double cl_point_distance(qaws_curve const* c, qaws_scalar t, double const* x)
{
	double pos[3];
	if (c->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
		pos[0] = e.position.x; pos[1] = e.position.y;
	}
	else
	{
		qaws_eval_result_3d e;
		qaws_curve_evaluate_3d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
		pos[0] = e.position.x; pos[1] = e.position.y; pos[2] = e.position.z;
	}
	return cl_dist(pos, x, (unsigned int)c->dimension);
}

/* Newton on (C - x) . C' = 0 from t, kept in the range; the closer of the start and the end */
static qaws_scalar cl_polish(qaws_curve const* c, qaws_scalar t, double const* x, double* d)
{
	qaws_scalar t0 = t, lo = c->parameter_range.min_value, hi = c->parameter_range.max_value;
	unsigned int it, k, dim = (unsigned int)c->dimension;
	double dn;
	for (it = 0; it < 24; it++)
	{
		double v[3][3], f = 0, g = 0, step;
		if (dim == 2)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
			v[0][0] = e.position.x; v[0][1] = e.position.y; v[1][0] = e.d1.x; v[1][1] = e.d1.y; v[2][0] = e.d2.x; v[2][1] = e.d2.y;
		}
		else
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(c, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
			v[0][0] = e.position.x; v[0][1] = e.position.y; v[0][2] = e.position.z;
			v[1][0] = e.d1.x; v[1][1] = e.d1.y; v[1][2] = e.d1.z; v[2][0] = e.d2.x; v[2][1] = e.d2.y; v[2][2] = e.d2.z;
		}
		for (k = 0; k < dim; k++)
		{
			f += (v[0][k] - x[k]) * v[1][k];
			g += v[1][k] * v[1][k] + (v[0][k] - x[k]) * v[2][k];
		}
		if (!(g > 0))
			break;
		step = f / g;
		t = (qaws_scalar)((double)t - step);
		t = t < lo ? lo : t > hi ? hi : t;
		if (fabs(step) <= 1e-15 * ((double)hi - (double)lo))
			break;
	}
	dn = cl_point_distance(c, t, x);
	if (dn < *d)
	{
		*d = dn;
		return t;
	}
	return t0;
}

qaws_status qaws_internal_closest_certify(qaws_curve const* curve, double const* x, qaws_scalar* io_t, double* out_distance,
	double* out_gap)
{
	double win[CL_MAX * 4], uk[2 * CL_MAX], bez[CL_MAX * 4], sub[CL_MAX * 4], qp[CL_MAX * 3];
	double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 }, ext = 0, tol, best, open = 1e300;
	struct { double u0, u1, lb; unsigned int depth; } stack[CL_DEPTH + 4];
	unsigned int dim = (unsigned int)curve->dimension, h = dim + 1, s, p, i, c, top, work = 0;
	qaws_scalar best_t = *io_t;

	/* the extent, for the tolerance, and the families and weights this handles */
	for (s = 0; s < curve->span_count; s++)
	{
		if (!cl_window(curve, s, &p, win, uk))
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		for (i = 0; i <= p; i++)
		{
			if (!(win[i * h + dim] > 0))
				return QAWS_STATUS_UNSUPPORTED_OPERATION;
			for (c = 0; c < dim; c++)
			{
				double v = win[i * h + c] / win[i * h + dim];
				lo[c] = v < lo[c] ? v : lo[c];
				hi[c] = v > hi[c] ? v : hi[c];
			}
		}
	}
	for (c = 0; c < dim; c++)
		ext += (hi[c] - lo[c]) * (hi[c] - lo[c]);
	tol = sqrt(ext) * ((double)QAWS_EPSILON > 1e-9 ? (double)QAWS_EPSILON : 1e-9) + 1e-300;
	best = cl_point_distance(curve, best_t, x);

	for (s = 0; s < curve->span_count; s++)
	{
		double ta = curve->span_boundaries[s], tb = curve->span_boundaries[s + 1];
		double root_diag;
		cl_window(curve, s, &p, win, uk);
		/* the span lies in its window's hull too: most spans end here, before the blossom */
		for (i = 0; i <= p; i++)
			for (c = 0; c < dim; c++)
				qp[i * 3 + c] = win[i * h + c] / win[i * h + dim];
		if (cl_lower(qp, p, dim, x, &root_diag) >= best)
			continue;
		cl_bezier(win, uk, p, h, bez);
		stack[0].u0 = 0; stack[0].u1 = 1; stack[0].lb = 0; stack[0].depth = 0;
		top = 1;
		while (top > 0)
		{
			double u0 = stack[top - 1].u0, u1 = stack[top - 1].u1, lb = stack[top - 1].lb, diag, d0, d1, m;
			unsigned int depth = stack[--top].depth;
			if (lb >= best)
				continue;
			memcpy(sub, bez, sizeof(double) * (p + 1) * h);
			if (u0 > 0)
				cl_split(sub, p, h, u0, 1);
			if (u1 < 1)
				cl_split(sub, p, h, (u1 - u0) / (1 - u0), 0);
			for (i = 0; i <= p; i++)
				for (c = 0; c < dim; c++)
					qp[i * 3 + c] = sub[i * h + c] / sub[i * h + dim];
			lb = cl_lower(qp, p, dim, x, &diag);
			/* the ends are on the curve */
			d0 = cl_dist(qp, x, dim);
			d1 = cl_dist(qp + p * 3, x, dim);
			if (d0 < best)
			{
				best = d0;
				best_t = (qaws_scalar)(ta + u0 * (tb - ta));
			}
			if (d1 < best)
			{
				best = d1;
				best_t = (qaws_scalar)(ta + u1 * (tb - ta));
			}
			if (lb >= best)
				continue;
			if (diag <= tol || depth >= CL_DEPTH || ++work > CL_WORK)
			{
				open = lb < open ? lb : open;
				continue;
			}
			/* the half nearer x on top */
			m = 0.5 * (u0 + u1);
			stack[top].u0 = d0 < d1 ? m : u0; stack[top].u1 = d0 < d1 ? u1 : m; stack[top].lb = lb; stack[top].depth = depth + 1; top++;
			stack[top].u0 = d0 < d1 ? u0 : m; stack[top].u1 = d0 < d1 ? m : u1; stack[top].lb = lb; stack[top].depth = depth + 1; top++;
		}
	}
	best_t = cl_polish(curve, best_t, x, &best);
	*io_t = best_t;
	*out_distance = best;
	*out_gap = open < best ? best - open : 0;
	return QAWS_STATUS_OK;
}

/* ---- surfaces ------------------------------------------------------------ */

#define SF_MAX 16
#define SF_DEPTH 48
#define SF_WORK 8192

/* a surface's control net: degrees, counts, knots (u, v); the points through sf_point */
typedef struct sf_net
{
	qaws_surface const* s;
	unsigned int p[2], n[2];
	qaws_scalar const* k[2];
	qaws_scalar const* cp;              /* owned: 3 scalars per point, u-major */
	qaws_scalar const* w;
	qaws_scalar bez[2][2 * SF_MAX];     /* Bezier: the clamped knots of [0, 1] */
} sf_net;

static int sf_net_init(qaws_surface const* s, sf_net* m)
{
	unsigned int d, i;
	memset(m, 0, sizeof *m);
	m->s = s;
	if (qaws_internal_surface_view_net(s, m->p, m->n, &m->k[0], &m->k[1]))
		return m->p[0] < SF_MAX && m->p[1] < SF_MAX;
	if (s->kind == QAWS_SURFACE_KIND_BEZIER)
	{
		qaws_surface_bezier_impl const* b = (qaws_surface_bezier_impl const*)s->impl;
		m->n[0] = b->u_count; m->n[1] = b->v_count; m->cp = b->control_points;
		m->p[0] = b->u_count - 1; m->p[1] = b->v_count - 1;
		if (m->p[0] >= SF_MAX || m->p[1] >= SF_MAX)
			return 0;
		for (d = 0; d < 2; d++)
		{
			for (i = 0; i < 2 * (m->p[d] + 1); i++)
				m->bez[d][i] = i <= m->p[d] ? QAWS_ZERO : QAWS_ONE;
			m->k[d] = m->bez[d];
		}
		return 1;
	}
	if (s->kind == QAWS_SURFACE_KIND_BSPLINE)
	{
		qaws_surface_bspline_impl const* b = (qaws_surface_bspline_impl const*)s->impl;
		m->n[0] = b->u_count; m->n[1] = b->v_count; m->cp = b->control_points; m->k[0] = b->u_knots; m->k[1] = b->v_knots;
	}
	else if (s->kind == QAWS_SURFACE_KIND_NURBS)
	{
		qaws_surface_nurbs_impl const* b = (qaws_surface_nurbs_impl const*)s->impl;
		m->n[0] = b->u_count; m->n[1] = b->v_count; m->cp = b->control_points; m->k[0] = b->u_knots; m->k[1] = b->v_knots;
		m->w = b->weights;
	}
	else
		return 0;
	m->p[0] = s->u_degree;
	m->p[1] = s->v_degree;
	return m->p[0] < SF_MAX && m->p[1] < SF_MAX;
}

/* net point (i, j) homogeneous */
static void sf_point(sf_net const* m, unsigned int i, unsigned int j, double* h)
{
	unsigned int e = i * m->n[1] + j, c;
	double w;
	if (!m->cp)
	{
		qaws_internal_surface_view_point(m->s, i, j, h);
		return;
	}
	w = m->w ? (double)m->w[e] : 1.0;
	for (c = 0; c < 3; c++)
		h[c] = w * (double)m->cp[e * 3 + c];
	h[3] = w;
}

/* the patch over knot spans (ku, kv) as a Bezier net: rows (v) blossomed, then columns (u) */
static void sf_bezier(sf_net const* m, unsigned int ku, unsigned int kv, double* bez)
{
	double win[SF_MAX * SF_MAX * 4], col[SF_MAX * 4], out[SF_MAX * 4], uk[2 * SF_MAX], vk[2 * SF_MAX];
	unsigned int pu = m->p[0], pv = m->p[1], a, b;
	for (a = 0; a <= pu; a++)
		for (b = 0; b <= pv; b++)
			sf_point(m, ku - pu + a, kv - pv + b, win + (a * (pv + 1) + b) * 4);
	for (a = 0; a < 2 * pu; a++)
		uk[a] = (double)m->k[0][ku - pu + 1 + a];
	for (b = 0; b < 2 * pv; b++)
		vk[b] = (double)m->k[1][kv - pv + 1 + b];
	for (a = 0; a <= pu; a++)
		cl_bezier(win + a * (pv + 1) * 4, vk, pv, 4, bez + a * (pv + 1) * 4);
	for (b = 0; b <= pv; b++)
	{
		for (a = 0; a <= pu; a++)
			memcpy(col + a * 4, bez + (a * (pv + 1) + b) * 4, sizeof(double) * 4);
		cl_bezier(col, uk, pu, 4, out);
		for (a = 0; a <= pu; a++)
			memcpy(bez + (a * (pv + 1) + b) * 4, out + a * 4, sizeof(double) * 4);
	}
}

/* the piece [a0, a1] x [b0, b1] of a Bezier net */
static void sf_sub(double const* bez, unsigned int pu, unsigned int pv, double a0, double a1, double b0, double b1, double* sub)
{
	double col[SF_MAX * 4];
	unsigned int a, b;
	memcpy(sub, bez, sizeof(double) * (pu + 1) * (pv + 1) * 4);
	for (a = 0; a <= pu; a++)
	{
		if (b0 > 0)
			cl_split(sub + a * (pv + 1) * 4, pv, 4, b0, 1);
		if (b1 < 1)
			cl_split(sub + a * (pv + 1) * 4, pv, 4, (b1 - b0) / (1 - b0), 0);
	}
	for (b = 0; b <= pv; b++)
	{
		for (a = 0; a <= pu; a++)
			memcpy(col + a * 4, sub + (a * (pv + 1) + b) * 4, sizeof(double) * 4);
		if (a0 > 0)
			cl_split(col, pu, 4, a0, 1);
		if (a1 < 1)
			cl_split(col, pu, 4, (a1 - a0) / (1 - a0), 0);
		for (a = 0; a <= pu; a++)
			memcpy(sub + (a * (pv + 1) + b) * 4, col + a * 4, sizeof(double) * 4);
	}
}

static void sf_cross(double const* a, double const* b, double* r)
{
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
}

/* the distance from x to the hull of the projected net q: to its box, or to its prism in the patch frame (u and v along the corner
   edges, the normal across the corner diagonals: second order in the patch size); the box diagonal */
static double sf_lower(double const* q, unsigned int pu, unsigned int pv, double const* x, double* diag)
{
	unsigned int n = (pu + 1) * (pv + 1), i, c, k;
	double lo[3], hi[3], box = 0, dg = 0, e[3][3], f[3], l0, l2, smin[3], smax[3], obb = 0;
	double const* c00 = q;
	double const* cp0 = q + pu * (pv + 1) * 3;
	double const* c0q = q + pv * 3;
	double const* cpq = q + (pu * (pv + 1) + pv) * 3;
	for (c = 0; c < 3; c++)
	{
		lo[c] = hi[c] = q[c];
		for (i = 1; i < n; i++)
		{
			lo[c] = q[i * 3 + c] < lo[c] ? q[i * 3 + c] : lo[c];
			hi[c] = q[i * 3 + c] > hi[c] ? q[i * 3 + c] : hi[c];
		}
		box += x[c] < lo[c] ? (lo[c] - x[c]) * (lo[c] - x[c]) : x[c] > hi[c] ? (x[c] - hi[c]) * (x[c] - hi[c]) : 0;
		dg += (hi[c] - lo[c]) * (hi[c] - lo[c]);
		e[0][c] = cp0[c] - c00[c] + cpq[c] - c0q[c];
		f[c] = c0q[c] - c00[c] + cpq[c] - cp0[c];
	}
	*diag = sqrt(dg);
	box = sqrt(box);
	sf_cross(e[0], f, e[2]);
	l0 = sqrt(e[0][0] * e[0][0] + e[0][1] * e[0][1] + e[0][2] * e[0][2]);
	l2 = sqrt(e[2][0] * e[2][0] + e[2][1] * e[2][1] + e[2][2] * e[2][2]);
	if (!(l0 > 1e-300) || !(l2 > 1e-300))
		return box;
	for (c = 0; c < 3; c++)
	{
		e[0][c] /= l0;
		e[2][c] /= l2;
	}
	sf_cross(e[2], e[0], e[1]);
	/* in the plane the patch is a parallelogram (u and v along its corner edges), whose bounding rectangle would be a first
	   order larger: oblique coordinates (s along u, t along v), the normal apart */
	{
		double fn = f[0] * e[2][0] + f[1] * e[2][1] + f[2] * e[2][2], cu = 0, cv = 0, fl = 0, px, py, ps, pt, din, dn;
		for (c = 0; c < 3; c++)
			fl += (f[c] - fn * e[2][c]) * (f[c] - fn * e[2][c]);
		fl = sqrt(fl);
		for (c = 0; c < 3; c++)
		{
			cu += (f[c] - fn * e[2][c]) / fl * e[0][c];
			cv += (f[c] - fn * e[2][c]) / fl * e[1][c];
		}
		if (!(fabs(cv) > 1e-6))
			return box;
		smin[0] = smax[0] = smin[1] = smax[1] = smin[2] = smax[2] = 0;
		for (i = 0; i < n; i++)
		{
			double w[3], o[3];
			for (c = 0; c < 3; c++)
				w[c] = q[i * 3 + c] - c00[c];
			for (k = 0; k < 3; k++)
				o[k] = w[0] * e[k][0] + w[1] * e[k][1] + w[2] * e[k][2];
			o[1] /= cv;                  /* t */
			o[0] -= o[1] * cu;           /* s */
			for (k = 0; k < 3; k++)
			{
				smin[k] = o[k] < smin[k] ? o[k] : smin[k];
				smax[k] = o[k] > smax[k] ? o[k] : smax[k];
			}
		}
		px = (x[0] - c00[0]) * e[0][0] + (x[1] - c00[1]) * e[0][1] + (x[2] - c00[2]) * e[0][2];
		py = (x[0] - c00[0]) * e[1][0] + (x[1] - c00[1]) * e[1][1] + (x[2] - c00[2]) * e[1][2];
		dn = (x[0] - c00[0]) * e[2][0] + (x[1] - c00[1]) * e[2][1] + (x[2] - c00[2]) * e[2][2];
		dn = dn < smin[2] ? smin[2] - dn : dn > smax[2] ? dn - smax[2] : 0;
		pt = py / cv;
		ps = px - pt * cu;
		din = 0;
		if (ps < smin[0] || ps > smax[0] || pt < smin[1] || pt > smax[1])
		{
			/* the nearest of the four sides, in the plane */
			double cs[4] = { smin[0], smax[0], smax[0], smin[0] }, ct[4] = { smin[1], smin[1], smax[1], smax[1] };
			din = 1e300;
			for (k = 0; k < 4; k++)
			{
				double ax = cs[k] + ct[k] * cu, ay = ct[k] * cv;
				double bx = cs[(k + 1) % 4] + ct[(k + 1) % 4] * cu, by = ct[(k + 1) % 4] * cv;
				double dx = bx - ax, dy = by - ay, l = dx * dx + dy * dy, g = l > 0 ? ((px - ax) * dx + (py - ay) * dy) / l : 0;
				g = g < 0 ? 0 : g > 1 ? 1 : g;
				g = (px - ax - g * dx) * (px - ax - g * dx) + (py - ay - g * dy) * (py - ay - g * dy);
				din = g < din ? g : din;
			}
			din = sqrt(din);
		}
		obb = sqrt(din * din + dn * dn);
	}
	return obb > box ? obb : box;
}

static double sf_distance(qaws_surface const* s, qaws_scalar u, qaws_scalar v, double const* x)
{
	qaws_surface_eval_result r;
	memset(&r, 0, sizeof r);
	qaws_surface_evaluate(s, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
	return sqrt(((double)r.position.x - x[0]) * ((double)r.position.x - x[0]) + ((double)r.position.y - x[1]) * ((double)r.position.y - x[1])
		+ ((double)r.position.z - x[2]) * ((double)r.position.z - x[2]));
}

/* Newton on the gradient of the squared distance while its Hessian is positive definite, kept in the domain; kept if closer */
static void sf_polish(qaws_surface const* s, double const* x, qaws_scalar* io_u, qaws_scalar* io_v, double* d)
{
	qaws_range ur = s->u_range, vr = s->v_range;
	double u = *io_u, v = *io_v, dn;
	unsigned int it;
	for (it = 0; it < 30; it++)
	{
		qaws_surface_eval_result e;
		double r[3], gu, gv, huu, huv, hvv, det, su, sv;
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
			| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV, &e);
		r[0] = e.position.x - x[0]; r[1] = e.position.y - x[1]; r[2] = e.position.z - x[2];
		gu = r[0] * e.du.x + r[1] * e.du.y + r[2] * e.du.z;
		gv = r[0] * e.dv.x + r[1] * e.dv.y + r[2] * e.dv.z;
		huu = (double)e.du.x * e.du.x + (double)e.du.y * e.du.y + (double)e.du.z * e.du.z + r[0] * e.duu.x + r[1] * e.duu.y + r[2] * e.duu.z;
		huv = (double)e.du.x * e.dv.x + (double)e.du.y * e.dv.y + (double)e.du.z * e.dv.z + r[0] * e.duv.x + r[1] * e.duv.y + r[2] * e.duv.z;
		hvv = (double)e.dv.x * e.dv.x + (double)e.dv.y * e.dv.y + (double)e.dv.z * e.dv.z + r[0] * e.dvv.x + r[1] * e.dvv.y + r[2] * e.dvv.z;
		det = huu * hvv - huv * huv;
		if (!(huu > 0) || !(det > 0))
			break;
		su = (hvv * gu - huv * gv) / det;
		sv = (huu * gv - huv * gu) / det;
		u -= su;
		v -= sv;
		u = u < ur.min_value ? ur.min_value : u > ur.max_value ? ur.max_value : u;
		v = v < vr.min_value ? vr.min_value : v > vr.max_value ? vr.max_value : v;
		if (fabs(su) <= 1e-15 * ((double)ur.max_value - ur.min_value) && fabs(sv) <= 1e-15 * ((double)vr.max_value - vr.min_value))
			break;
	}
	dn = sf_distance(s, (qaws_scalar)u, (qaws_scalar)v, x);
	if (dn < *d)
	{
		*d = dn;
		*io_u = (qaws_scalar)u;
		*io_v = (qaws_scalar)v;
	}
}

qaws_status qaws_internal_surface_closest_certify(qaws_surface const* surface, double const* x, qaws_scalar* io_u, qaws_scalar* io_v,
	double* out_distance, double* out_gap)
{
	double bez[SF_MAX * SF_MAX * 4], sub[SF_MAX * SF_MAX * 4], qp[SF_MAX * SF_MAX * 3], h[4];
	double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 }, ext = 0, tol, best, open = 1e300;
	struct { double a0, a1, b0, b1, lb; unsigned int depth; } stack[3 * SF_DEPTH + 8];
	sf_net m;
	unsigned int i, j, c, ku, kv, top, work = 0, pu, pv, np;
	qaws_scalar bu = *io_u, bv = *io_v;

	if (!sf_net_init(surface, &m))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	pu = m.p[0];
	pv = m.p[1];
	np = (pu + 1) * (pv + 1);
	/* the extent, for the tolerance; positive weights */
	for (i = 0; i < m.n[0]; i++)
		for (j = 0; j < m.n[1]; j++)
		{
			sf_point(&m, i, j, h);
			if (!(h[3] > 0))
				return QAWS_STATUS_UNSUPPORTED_OPERATION;
			for (c = 0; c < 3; c++)
			{
				lo[c] = h[c] / h[3] < lo[c] ? h[c] / h[3] : lo[c];
				hi[c] = h[c] / h[3] > hi[c] ? h[c] / h[3] : hi[c];
			}
		}
	for (c = 0; c < 3; c++)
		ext += (hi[c] - lo[c]) * (hi[c] - lo[c]);
	tol = sqrt(ext) * ((double)QAWS_EPSILON > 1e-9 ? (double)QAWS_EPSILON : 1e-9) + 1e-300;
	best = sf_distance(surface, bu, bv, x);

	for (ku = pu; ku < m.n[0]; ku++)
	{
		double u0 = m.k[0][ku], u1 = m.k[0][ku + 1];
		if (!(u1 > u0))
			continue;
		for (kv = pv; kv < m.n[1]; kv++)
		{
			double v0 = m.k[1][kv], v1 = m.k[1][kv + 1], diag;
			if (!(v1 > v0))
				continue;
			/* the patch lies in its window hull: most patches end here, before the blossom */
			for (i = 0; i <= pu; i++)
				for (j = 0; j <= pv; j++)
				{
					sf_point(&m, ku - pu + i, kv - pv + j, h);
					for (c = 0; c < 3; c++)
						qp[(i * (pv + 1) + j) * 3 + c] = h[c] / h[3];
				}
			if (sf_lower(qp, pu, pv, x, &diag) >= best)
				continue;
			sf_bezier(&m, ku, kv, bez);
			stack[0].a0 = 0; stack[0].a1 = 1; stack[0].b0 = 0; stack[0].b1 = 1; stack[0].lb = 0; stack[0].depth = 0;
			top = 1;
			while (top > 0)
			{
				double a0 = stack[top - 1].a0, a1 = stack[top - 1].a1, b0 = stack[top - 1].b0, b1 = stack[top - 1].b1, lb = stack[top - 1].lb;
				double cd[4], am, bm;
				unsigned int depth = stack[--top].depth, order[4] = { 0, 1, 2, 3 }, q, r;
				if (lb >= best)
					continue;
				sf_sub(bez, pu, pv, a0, a1, b0, b1, sub);
				for (i = 0; i < np; i++)
					for (c = 0; c < 3; c++)
						qp[i * 3 + c] = sub[i * 4 + c] / sub[i * 4 + 3];
				lb = sf_lower(qp, pu, pv, x, &diag);
				/* the corners are on the surface: (a0, b0), (a1, b0), (a0, b1), (a1, b1) */
				cd[0] = cl_dist(qp, x, 3);
				cd[1] = cl_dist(qp + pu * (pv + 1) * 3, x, 3);
				cd[2] = cl_dist(qp + pv * 3, x, 3);
				cd[3] = cl_dist(qp + (pu * (pv + 1) + pv) * 3, x, 3);
				for (q = 0; q < 4; q++)
					if (cd[q] < best)
					{
						best = cd[q];
						bu = (qaws_scalar)(u0 + (q & 1 ? a1 : a0) * (u1 - u0));
						bv = (qaws_scalar)(v0 + (q & 2 ? b1 : b0) * (v1 - v0));
					}
				if (lb >= best)
					continue;
				if (diag <= tol || depth >= SF_DEPTH || ++work > SF_WORK)
				{
					open = lb < open ? lb : open;
					continue;
				}
				/* the quarters, the one at the nearest corner on top */
				for (q = 1; q < 4; q++)
					for (r = q; r > 0 && cd[order[r]] > cd[order[r - 1]]; r--)
					{
						unsigned int t = order[r];
						order[r] = order[r - 1];
						order[r - 1] = t;
					}
				am = 0.5 * (a0 + a1);
				bm = 0.5 * (b0 + b1);
				for (q = 0; q < 4; q++)
				{
					stack[top].a0 = order[q] & 1 ? am : a0; stack[top].a1 = order[q] & 1 ? a1 : am;
					stack[top].b0 = order[q] & 2 ? bm : b0; stack[top].b1 = order[q] & 2 ? b1 : bm;
					stack[top].lb = lb; stack[top].depth = depth + 1;
					top++;
				}
			}
		}
	}
	sf_polish(surface, x, &bu, &bv, &best);
	*io_u = bu;
	*io_v = bv;
	*out_distance = best;
	*out_gap = open < best ? best - open : 0;
	return QAWS_STATUS_OK;
}
