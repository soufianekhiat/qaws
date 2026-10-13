/*
 * Geodesic queries on a surface (qaws_surface_region.h): distances, a
 * curve's closest point, two curves' crossings and closest points, a
 * point's distance to a region.
 *
 * Each distance is the length of the geodesic between two points, aimed
 * by shooting (qaws_internal_sg_reach) from the chord. A curve's closest
 * point to p inside it is where the geodesic from the curve to p leaves it
 * square: the root of f(t) = g(c'(t), w(t)), w the geodesic's unit
 * starting velocity, bracketed round the nearest sample.
 */

#include "qaws_surface_region.h"
#include "qaws_surface.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "internal/qaws_internal_surface_region.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SQ_SAMPLES 64

static double sq_dist3(qaws_vec3 a, qaws_vec3 b)
{
	double dx = (double)a.x - b.x, dy = (double)a.y - b.y, dz = (double)a.z - b.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}

static qaws_vec3 sq_pos(qaws_internal_sg_ctx const* c, double const* uv)
{
	qaws_surface_eval_result r;
	memset(&r, 0, sizeof(r));
	qaws_internal_sg_eval(c, uv[0], uv[1], QAWS_SURFACE_EVAL_POSITION, &r);
	return r.position;
}

/* the geodesic from a to x: its length, and (may be NULL) its unit (u, v)
   velocity leaving a; th and L warm starts (L <= 0: from the chord) */
static int sq_geo(qaws_internal_sg_ctx const* c, double const* a, double const* x, double* th, double* L, double* w)
{
	double end[2], vel[2];
	if (!(*L > 0))
	{
		double d[2];
		d[0] = x[0] - a[0]; d[1] = x[1] - a[1];
		if (c->pu > 0) d[0] -= c->pu * floor(d[0] / c->pu + 0.5);
		if (c->pv > 0) d[1] -= c->pv * floor(d[1] / c->pv + 0.5);
		if (!(hypot(d[0], d[1]) > 0) || qaws_internal_sg_angle(c, a[0], a[1], d, th) != QAWS_STATUS_OK)
			return 0;
		*L = sq_dist3(sq_pos(c, a), sq_pos(c, x));
	}
	if (!qaws_internal_sg_reach(c, a, x, th, L))
		return 0;
	if (w)
	{
		if (qaws_internal_sg_shoot(c, a, *th, 0, end, vel) != QAWS_STATUS_OK)
			return 0;
		w[0] = vel[0]; w[1] = vel[1];
	}
	return 1;
}

qaws_status qaws_surface_distance(qaws_surface const* surface, qaws_vec2 a, qaws_vec2 b, qaws_scalar tolerance, qaws_scalar* out_distance)
{
	qaws_internal_sg_ctx c;
	double pa[2], pb[2], th = 0, L = 0;
	qaws_status st;
	if (!surface || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_internal_sg_init(&c, surface, (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	pa[0] = a.x; pa[1] = a.y; pb[0] = b.x; pb[1] = b.y;
	if (sq_dist3(sq_pos(&c, pa), sq_pos(&c, pb)) <= 1e-3 * c.tol)
	{
		*out_distance = 0;
		return QAWS_STATUS_OK;
	}
	if (!sq_geo(&c, pa, pb, &th, &L, NULL))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	*out_distance = (qaws_scalar)L;
	return QAWS_STATUS_OK;
}

/* the curve's point and (u, v) derivative at t */
static void sq_curve_at(qaws_curve const* cv, double t, double* p, double* d)
{
	qaws_eval_result_2d e;
	qaws_curve_evaluate_2d(cv, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
	p[0] = e.position.x; p[1] = e.position.y;
	if (d) { d[0] = e.d1.x; d[1] = e.d1.y; }
}

/* f(t) = g(c'(t), w(t)) and the distance at t; 0 when the geodesic fails */
static int sq_angle_at(qaws_internal_sg_ctx const* c, qaws_curve const* cv, double t, double const* x, double* th, double* L, double* f,
	double* dist)
{
	double p[2], d[2], w[2], g[3];
	sq_curve_at(cv, t, p, d);
	if (sq_dist3(sq_pos(c, p), sq_pos(c, x)) <= 1e-3 * c->tol)
	{
		*f = 0;
		*dist = 0;
		return 1;
	}
	if (!sq_geo(c, p, x, th, L, w) || qaws_internal_sg_first_form(c, p[0], p[1], g) != QAWS_STATUS_OK)
		return 0;
	*f = g[0] * d[0] * w[0] + g[1] * (d[0] * w[1] + d[1] * w[0]) + g[2] * d[1] * w[1];
	*dist = *L;
	return 1;
}

/* the closest point of the (u, v) curve cv to x: parameter, distance */
static int sq_closest(qaws_internal_sg_ctx const* c, qaws_curve const* cv, double const* x, double* out_t, double* out_d)
{
	qaws_range rg = qaws_curve_get_parameter_range(cv);
	double t0 = (double)rg.min_value, t1 = (double)rg.max_value, best = 1e300, bt = t0, e3[SQ_SAMPLES + 1];
	unsigned int i, k, cand[3], nc = 0;
	qaws_vec3 X = sq_pos(c, x);
	for (i = 0; i <= SQ_SAMPLES; i++)
	{
		double p[2];
		sq_curve_at(cv, t0 + (t1 - t0) * i / SQ_SAMPLES, p, NULL);
		e3[i] = sq_dist3(sq_pos(c, p), X);
	}
	/* the three nearest local minima of the samples (by the chord in space) */
	for (k = 0; k < 3; k++)
	{
		unsigned int bi = SQ_SAMPLES + 1;
		for (i = 0; i <= SQ_SAMPLES; i++)
		{
			unsigned int j;
			int used = 0;
			if ((i > 0 && e3[i - 1] < e3[i]) || (i < SQ_SAMPLES && e3[i + 1] < e3[i]))
				continue;
			for (j = 0; j < nc; j++) used |= cand[j] == i;
			if (!used && (bi > SQ_SAMPLES || e3[i] < e3[bi]))
				bi = i;
		}
		if (bi <= SQ_SAMPLES)
			cand[nc++] = bi;
	}
	for (k = 0; k < nc; k++)
	{
		/* the angle's root in the samples either side, by regula falsi
		   (Illinois); an end of the curve taken as it is */
		unsigned int ia = cand[k] > 0 ? cand[k] - 1 : 0, ib = cand[k] < SQ_SAMPLES ? cand[k] + 1 : SQ_SAMPLES, it;
		double a = t0 + (t1 - t0) * ia / SQ_SAMPLES, b = t0 + (t1 - t0) * ib / SQ_SAMPLES, fa, fb, da, db, tha = 0, La = 0, thb = 0, Lb = 0;
		int okA = sq_angle_at(c, cv, a, x, &tha, &La, &fa, &da), okB = sq_angle_at(c, cv, b, x, &thb, &Lb, &fb, &db), side = 0;
		if (okA && da < best) { best = da; bt = a; }
		if (okB && db < best) { best = db; bt = b; }
		if (!okA || !okB || (fa > 0) == (fb > 0))
		{
			/* no root between: the sample itself, then an end */
			double tm = t0 + (t1 - t0) * cand[k] / SQ_SAMPLES, fm, dm, th = 0, L = 0;
			if (sq_angle_at(c, cv, tm, x, &th, &L, &fm, &dm) && dm < best) { best = dm; bt = tm; }
			continue;
		}
		for (it = 0; it < 80 && fabs(b - a) > 1e-14 * (fabs(t1 - t0) + 1); it++)
		{
			double m = (a * fb - b * fa) / (fb - fa), fm, dm, th = tha, L = La;
			if (!(m > a && m < b)) m = 0.5 * (a + b);
			if (!sq_angle_at(c, cv, m, x, &th, &L, &fm, &dm))
				break;
			if (dm < best) { best = dm; bt = m; }
			if (fm == 0)
				break;
			if ((fm > 0) == (fa > 0))
			{
				a = m; fa = fm; tha = th; La = L;
				if (side == -1) fb *= 0.5;
				side = -1;
			}
			else
			{
				b = m; fb = fm;
				if (side == 1) fa *= 0.5;
				side = 1;
			}
		}
	}
	{
		/* the ends */
		double th = 0, L = 0, f, d;
		if (sq_angle_at(c, cv, t0, x, &th, &L, &f, &d) && d < best) { best = d; bt = t0; }
		th = 0; L = 0;
		if (sq_angle_at(c, cv, t1, x, &th, &L, &f, &d) && d < best) { best = d; bt = t1; }
	}
	if (!(best < 1e300))
		return 0;
	*out_t = bt;
	*out_d = best;
	return 1;
}

qaws_status qaws_surface_curve_closest_point(qaws_surface_curve const* curve, qaws_vec2 p, qaws_scalar tolerance,
	qaws_scalar* out_parameter, qaws_vec2* out_point, qaws_scalar* out_distance)
{
	qaws_internal_sg_ctx c;
	double x[2], t, d, q[2];
	qaws_curve const* cv;
	qaws_status st;
	if (!curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_internal_sg_init(&c, qaws_surface_curve_get_surface(curve), (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	cv = qaws_surface_curve_get_uv(curve);
	x[0] = p.x; x[1] = p.y;
	if (!sq_closest(&c, cv, x, &t, &d))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	sq_curve_at(cv, t, q, NULL);
	if (out_parameter) *out_parameter = (qaws_scalar)t;
	if (out_point) { out_point->x = (qaws_scalar)q[0]; out_point->y = (qaws_scalar)q[1]; }
	if (out_distance) *out_distance = (qaws_scalar)d;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_curve_intersect(qaws_surface_curve const* a, qaws_surface_curve const* b, qaws_surface_crossing* out_crossings,
	unsigned int capacity, unsigned int* out_count)
{
	qaws_internal_sg_ctx c;
	qaws_curve const *ca, *cb;
	unsigned int count = 0;
	int ku, kv, nu, nv;
	qaws_status st;
	if (!a || !b || !out_count || (capacity && !out_crossings))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_surface_curve_get_surface(a) != qaws_surface_curve_get_surface(b))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_internal_sg_init(&c, qaws_surface_curve_get_surface(a), 0)) != QAWS_STATUS_OK)
		return st;
	ca = qaws_surface_curve_get_uv(a);
	cb = qaws_surface_curve_get_uv(b);
	/* b's copies a period either way in each closed direction (the images
	   live in the covering plane) */
	nu = c.pu > 0 ? 2 : 0;
	nv = c.pv > 0 ? 2 : 0;
	for (ku = -nu; ku <= nu && st == QAWS_STATUS_OK; ku++)
		for (kv = -nv; kv <= nv && st == QAWS_STATUS_OK; kv++)
		{
			qaws_curve* moved = NULL;
			qaws_intersection_2d hits[64];
			unsigned int n = 0, i;
			qaws_scalar m[6];
			m[0] = 1; m[1] = 0; m[2] = (qaws_scalar)(ku * c.pu); m[3] = 0; m[4] = 1; m[5] = (qaws_scalar)(kv * c.pv);
			if ((ku || kv) && (st = qaws_curve_transform_2d(cb, m, &moved)) != QAWS_STATUS_OK)
				break;
			st = qaws_curve_find_intersections_2d(ca, moved ? moved : cb, hits, 64, &n);
			qaws_curve_destroy(moved);
			if (st == QAWS_STATUS_BUFFER_TOO_SMALL)
				st = QAWS_STATUS_OK;
			for (i = 0; i < n && i < 64 && st == QAWS_STATUS_OK; i++, count++)
				if (count < capacity)
				{
					out_crossings[count].parameter_a = hits[i].parameter_a;
					out_crossings[count].parameter_b = hits[i].parameter_b;
					out_crossings[count].point = hits[i].position;
				}
		}
	*out_count = count;
	return st;
}

qaws_status qaws_surface_curve_closest_points(qaws_surface_curve const* a, qaws_surface_curve const* b, qaws_scalar tolerance,
	qaws_scalar* out_parameter_a, qaws_scalar* out_parameter_b, qaws_vec2* out_point_a, qaws_vec2* out_point_b,
	qaws_scalar* out_distance)
{
	qaws_internal_sg_ctx c;
	qaws_curve const *ca, *cb;
	qaws_surface_crossing x;
	unsigned int n = 0, i, j, it;
	double ta, tb, d = 0, pa[2], pb[2];
	qaws_status st;
	if (!a || !b)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_surface_curve_intersect(a, b, &x, 1, &n)) != QAWS_STATUS_OK)
		return st;
	if ((st = qaws_internal_sg_init(&c, qaws_surface_curve_get_surface(a), (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	ca = qaws_surface_curve_get_uv(a);
	cb = qaws_surface_curve_get_uv(b);
	if (n)
	{
		ta = x.parameter_a; tb = x.parameter_b;
		sq_curve_at(ca, ta, pa, NULL);
		sq_curve_at(cb, tb, pb, NULL);
	}
	else
	{
		/* the nearest pair of samples (in space), then closest points in turn */
		qaws_range ra = qaws_curve_get_parameter_range(ca), rb = qaws_curve_get_parameter_range(cb);
		qaws_vec3 P[SQ_SAMPLES + 1], Q[SQ_SAMPLES + 1];
		double best = 1e300, prev = 1e300;
		ta = ra.min_value; tb = rb.min_value;
		for (i = 0; i <= SQ_SAMPLES; i++)
		{
			sq_curve_at(ca, ra.min_value + (ra.max_value - ra.min_value) * (double)i / SQ_SAMPLES, pa, NULL);
			sq_curve_at(cb, rb.min_value + (rb.max_value - rb.min_value) * (double)i / SQ_SAMPLES, pb, NULL);
			P[i] = sq_pos(&c, pa);
			Q[i] = sq_pos(&c, pb);
		}
		for (i = 0; i <= SQ_SAMPLES; i++)
			for (j = 0; j <= SQ_SAMPLES; j++)
				if (sq_dist3(P[i], Q[j]) < best)
				{
					best = sq_dist3(P[i], Q[j]);
					ta = ra.min_value + (ra.max_value - ra.min_value) * (double)i / SQ_SAMPLES;
					tb = rb.min_value + (rb.max_value - rb.min_value) * (double)j / SQ_SAMPLES;
				}
		for (it = 0; it < 60; it++)
		{
			double t2;
			sq_curve_at(ca, ta, pa, NULL);
			if (!sq_closest(&c, cb, pa, &tb, &d))
				return QAWS_STATUS_NUMERICAL_FAILURE;
			sq_curve_at(cb, tb, pb, NULL);
			t2 = ta;
			if (!sq_closest(&c, ca, pb, &ta, &d))
				return QAWS_STATUS_NUMERICAL_FAILURE;
			if (fabs(ta - t2) <= 1e-13 * (ra.max_value - ra.min_value) || !(d < prev - 1e-15 * (d + 1)))
				break;
			prev = d;
		}
		sq_curve_at(ca, ta, pa, NULL);
		sq_curve_at(cb, tb, pb, NULL);
	}
	if (out_parameter_a) *out_parameter_a = (qaws_scalar)ta;
	if (out_parameter_b) *out_parameter_b = (qaws_scalar)tb;
	if (out_point_a) { out_point_a->x = (qaws_scalar)pa[0]; out_point_a->y = (qaws_scalar)pa[1]; }
	if (out_point_b) { out_point_b->x = (qaws_scalar)pb[0]; out_point_b->y = (qaws_scalar)pb[1]; }
	if (out_distance) *out_distance = (qaws_scalar)(n ? 0 : d);
	return QAWS_STATUS_OK;
}

/* a region's cut along a seam: a curve on a closed direction's line */
static int sq_on_seam(qaws_internal_sg_ctx const* c, qaws_curve const* cv)
{
	qaws_range rg = qaws_curve_get_parameter_range(cv);
	unsigned int i, ax;
	for (ax = 0; ax < 2; ax++)
	{
		double per = ax ? c->pv : c->pu, base = ax ? c->v0 : c->u0, w0 = 0;
		int on = per > 0;
		for (i = 0; i <= 4 && on; i++)
		{
			double p[2], w;
			sq_curve_at(cv, rg.min_value + (rg.max_value - rg.min_value) * i / 4.0, p, NULL);
			w = (p[ax] - base) / per;
			if (i == 0) w0 = w;
			on = fabs(w - floor(w + 0.5)) <= 1e-9 && fabs(w - w0) <= 1e-9;
		}
		if (on)
			return 1;
	}
	return 0;
}

qaws_status qaws_surface_region_distance(qaws_surface_region const* region, qaws_vec2 p, qaws_scalar tolerance,
	qaws_scalar* out_distance, qaws_vec2* out_nearest)
{
	qaws_internal_sg_ctx c;
	unsigned int np, i, j;
	double x[2], best = 1e300, bq[2] = { 0, 0 };
	int inside = 0;
	qaws_status st;
	if (!region || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = qaws_internal_sg_init(&c, qaws_surface_region_get_surface(region), (double)tolerance)) != QAWS_STATUS_OK)
		return st;
	x[0] = p.x; x[1] = p.y;
	np = qaws_surface_region_get_path_count(region);
	{
		/* the boundary's curves by a lower bound on their distance (the chord
		   in space to 16 samples less the samples' spacing: a geodesic is no
		   shorter than its chord), the nearest first; past the best, done */
		qaws_curve const** cv = NULL;
		double* lb = NULL;
		unsigned int n = 0, cap = 0, k;
		qaws_vec3 X = sq_pos(&c, x);
		for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
		{
			qaws_path_2d path;
			if ((st = qaws_surface_region_get_path(region, i, &path)) != QAWS_STATUS_OK)
				break;
			for (j = 0; j < path.curve_count; j++)
			{
				qaws_range rg = qaws_curve_get_parameter_range(path.curves[j]);
				double lo = 1e300, gap = 0, q[2];
				qaws_vec3 prev = X;
				if (sq_on_seam(&c, path.curves[j]))
					continue;
				if (n == cap)
				{
					unsigned int nc = cap ? 2 * cap : 64;
					qaws_curve const** g1 = (qaws_curve const**)realloc((void*)cv, sizeof(*cv) * nc);
					double* g2;
					if (g1) cv = g1;
					g2 = g1 ? (double*)realloc(lb, sizeof(double) * nc) : NULL;
					if (!g1 || !g2)
					{
						st = QAWS_STATUS_ALLOCATION_FAILURE;
						break;
					}
					lb = g2;
					cap = nc;
				}
				for (k = 0; k <= 16; k++)
				{
					qaws_vec3 y;
					sq_curve_at(path.curves[j], (double)rg.min_value + ((double)rg.max_value - rg.min_value) * k / 16.0, q, NULL);
					y = sq_pos(&c, q);
					if (sq_dist3(y, X) < lo) lo = sq_dist3(y, X);
					if (k && sq_dist3(y, prev) > gap) gap = sq_dist3(y, prev);
					prev = y;
				}
				cv[n] = path.curves[j];
				lb[n] = lo - gap;
				n++;
			}
		}
		while (st == QAWS_STATUS_OK)
		{
			unsigned int bi = n;
			double t, d;
			for (k = 0; k < n; k++)
				if (cv[k] && (bi == n || lb[k] < lb[bi]))
					bi = k;
			if (bi == n || lb[bi] > best)
				break;
			if (sq_closest(&c, cv[bi], x, &t, &d) && d < best)
			{
				best = d;
				sq_curve_at(cv[bi], t, bq, NULL);
			}
			cv[bi] = NULL;
		}
		free((void*)cv);
		free(lb);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	if (!(best < 1e300))
		return QAWS_STATUS_NUMERICAL_FAILURE;
	if ((st = qaws_surface_region_contains(region, p, &inside)) != QAWS_STATUS_OK)
		return st;
	*out_distance = (qaws_scalar)(inside == 1 ? -best : inside == -1 ? 0 : best);
	if (out_nearest) { out_nearest->x = (qaws_scalar)bq[0]; out_nearest->y = (qaws_scalar)bq[1]; }
	return QAWS_STATUS_OK;
}
