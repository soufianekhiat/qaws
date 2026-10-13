/*
 * Regions on surfaces (qaws_surface_region.h): Booleans through the 2D
 * engine on the curves' (u, v) images, areas and integrals by Green's
 * theorem with the surface's area element.
 *
 * For a density h(u, v) (the area element sqrt(E G - F^2) times f), with
 * P(u, v) = integral from u_ref to u of h(s, v) ds,
 *
 *     double integral of h du dv = sum over the loops of the line integral of P dv,
 *
 * outside loops counter-clockwise and holes clockwise. P is computed by
 * Gauss-Legendre panels along u (8 points per eighth of the domain), the
 * loop integral by 8-point Gauss on every span of the boundary curves.
 */
#include "qaws_surface_region.h"
#include "qaws_surface.h"
#include "qaws_surface_trim.h"
#include "qaws_curve.h"
#include "qaws_bezier.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface_region.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

struct qaws_surface_region
{
	qaws_surface const* s;
	qaws_internal_sg_ctx ctx;
	qaws_clip_result* res;
};

static double const sr_x8[4] = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
static double const sr_w8[4] = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

/* the paths' (u, v) images */
static qaws_path_2d* sr_paths(qaws_surface const* s, qaws_surface_path const* p, unsigned int n, qaws_curve const*** keep, qaws_status* st)
{
	qaws_path_2d* out;
	qaws_curve const** cv;
	unsigned int i, k, total = 0, at = 0;
	*st = QAWS_STATUS_OK;
	for (i = 0; i < n; i++) total += p[i].curve_count;
	out = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (n + 1));
	cv = (qaws_curve const**)malloc(sizeof(qaws_curve const*) * (total + 1));
	if (!out || !cv)
	{
		free(out); free(cv);
		*st = QAWS_STATUS_ALLOCATION_FAILURE;
		return NULL;
	}
	for (i = 0; i < n; i++)
	{
		out[i].curves = cv + at;
		out[i].curve_count = p[i].curve_count;
		out[i].closed = p[i].closed;
		for (k = 0; k < p[i].curve_count; k++)
		{
			if (!p[i].curves[k] || qaws_surface_curve_get_surface(p[i].curves[k]) != s)
				*st = QAWS_STATUS_INVALID_ARGUMENT;
			cv[at++] = p[i].curves[k] ? qaws_surface_curve_get_uv(p[i].curves[k]) : NULL;
		}
	}
	*keep = cv;
	return out;
}

/* ------------------------------------------------------------------------ */
/*  Closed directions: the regions cut into one period                      */
/* ------------------------------------------------------------------------ */

/* curves made on the way, destroyed after the clip */
typedef struct sr_pool
{
	qaws_curve** c;
	unsigned int n, cap;
	qaws_curve const** views;   /* path curve arrays */
	unsigned int nv, capv;
	qaws_path_2d* paths;
	unsigned int np, capp;
	qaws_clip_result** res;
	unsigned int nr, capr;
} sr_pool;

static int sr_grow(void** p, unsigned int* cap, unsigned int need, size_t size)
{
	unsigned int n;
	void* q;
	if (need <= *cap)
		return 1;
	n = *cap ? *cap : 16;
	while (n < need) n *= 2;
	q = realloc(*p, size * n);
	if (!q)
		return 0;
	*p = q;
	*cap = n;
	return 1;
}

static qaws_status sr_keep(sr_pool* pl, qaws_curve* c)
{
	if (!c)
		return QAWS_STATUS_NUMERICAL_FAILURE;
	if (!sr_grow((void**)&pl->c, &pl->cap, pl->n + 1, sizeof(qaws_curve*)))
	{
		qaws_curve_destroy(c);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	pl->c[pl->n++] = c;
	return QAWS_STATUS_OK;
}

static void sr_pool_free(sr_pool* pl)
{
	unsigned int i;
	for (i = 0; i < pl->n; i++) qaws_curve_destroy(pl->c[i]);
	for (i = 0; i < pl->nr; i++) qaws_clip_result_destroy(pl->res[i]);
	free(pl->c); free(pl->views); free(pl->paths); free(pl->res);
	memset(pl, 0, sizeof(*pl));
}

static qaws_status sr_line(sr_pool* pl, double x0, double y0, double x1, double y1, qaws_curve const** out)
{
	qaws_bezier_desc d;
	qaws_scalar cp[4];
	qaws_curve* c = NULL;
	qaws_status st;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
	st = qaws_curve_create_bezier(&d, &c);
	if (st != QAWS_STATUS_OK)
		return st;
	*out = c;
	return sr_keep(pl, c);
}

/* a path of n curves (copied views) */
static qaws_status sr_add_path(sr_pool* pl, qaws_curve const* const* curves, unsigned int n)
{
	unsigned int i;
	if (!sr_grow((void**)&pl->views, &pl->capv, pl->nv + n, sizeof(qaws_curve const*)) ||
		!sr_grow((void**)&pl->paths, &pl->capp, pl->np + 1, sizeof(qaws_path_2d)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n; i++) pl->views[pl->nv + i] = curves[i];
	/* the views move when the array grows: offsets until the end */
	pl->paths[pl->np].curves = (qaws_curve const* const*)(size_t)pl->nv;
	pl->paths[pl->np].curve_count = n;
	pl->paths[pl->np].closed = 1;
	pl->nv += n;
	pl->np++;
	return QAWS_STATUS_OK;
}

static void sr_fix_views(sr_pool* pl)
{
	unsigned int i;
	for (i = 0; i < pl->np; i++)
		pl->paths[i].curves = pl->views + (size_t)pl->paths[i].curves;
}

static double sr_signed_area(qaws_path_2d const* p)
{
	qaws_scalar a = 0;
	qaws_path_compute_area_2d(p, &a);
	return (double)a;
}

/* the loop's part in each period cell, translated into the first, its
   orientation kept */
static qaws_status sr_cut_loop(qaws_internal_sg_ctx const* c, sr_pool* pl, qaws_path_2d const* loop)
{
	qaws_vec2 lo, hi;
	double area = sr_signed_area(loop);
	int k0 = 0, k1 = 0, m0 = 0, m1 = 0, k, m;
	qaws_status st = qaws_path_compute_bounds_2d(loop, &lo, &hi);
	if (st != QAWS_STATUS_OK)
		return st;
	if (area == 0)
		return QAWS_STATUS_OK;
	if (c->pu > 0) { k0 = (int)floor(((double)lo.x - c->u0) / c->pu); k1 = (int)floor(((double)hi.x - c->u0) / c->pu - 1e-12); }
	if (c->pv > 0) { m0 = (int)floor(((double)lo.y - c->v0) / c->pv); m1 = (int)floor(((double)hi.y - c->v0) / c->pv - 1e-12); }
	for (k = k0; k <= k1; k++)
		for (m = m0; m <= m1; m++)
		{
			qaws_curve const* rc[4];
			qaws_path_2d cell;
			qaws_clip_result* res = NULL;
			double big = 0.25 * (c->diag + (double)hi.x - lo.x + (double)hi.y - lo.y);   /* past the loop, not so far as to hurt float clips */
			double ua = c->pu > 0 ? c->u0 + k * c->pu : (double)lo.x - big, ub = c->pu > 0 ? c->u0 + (k + 1) * c->pu : (double)hi.x + big;
			double va = c->pv > 0 ? c->v0 + m * c->pv : (double)lo.y - big, vb = c->pv > 0 ? c->v0 + (m + 1) * c->pv : (double)hi.y + big;
			unsigned int np, i, j;
			/* the cell runs the way the loop does: the rule applies to both */
			if (area > 0)
			{
				if ((st = sr_line(pl, ua, va, ub, va, &rc[0])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ub, va, ub, vb, &rc[1])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ub, vb, ua, vb, &rc[2])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ua, vb, ua, va, &rc[3])) != QAWS_STATUS_OK) return st;
			}
			else
			{
				if ((st = sr_line(pl, ua, va, ua, vb, &rc[0])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ua, vb, ub, vb, &rc[1])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ub, vb, ub, va, &rc[2])) != QAWS_STATUS_OK) return st;
				if ((st = sr_line(pl, ub, va, ua, va, &rc[3])) != QAWS_STATUS_OK) return st;
			}
			cell.curves = rc; cell.curve_count = 4; cell.closed = 1;
			st = qaws_clip_boolean(QAWS_CLIP_INTERSECTION, area > 0 ? QAWS_FILL_POSITIVE : QAWS_FILL_NEGATIVE, loop, 1, &cell, 1, &res);
			if (st != QAWS_STATUS_OK)
				return st;
			if (!sr_grow((void**)&pl->res, &pl->capr, pl->nr + 1, sizeof(qaws_clip_result*)))
			{
				qaws_clip_result_destroy(res);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			pl->res[pl->nr++] = res;
			np = qaws_clip_result_get_path_count(res);
			for (i = 0; i < np; i++)
			{
				qaws_path_2d p;
				qaws_curve const* moved[512];
				qaws_scalar t[6];
				if ((st = qaws_clip_result_get_path(res, i, &p)) != QAWS_STATUS_OK)
					return st;
				if (p.curve_count > 512)
					return QAWS_STATUS_UNSUPPORTED_OPERATION;
				t[0] = 1; t[1] = 0; t[2] = (qaws_scalar)(c->pu > 0 ? -k * c->pu : 0);
				t[3] = 0; t[4] = 1; t[5] = (qaws_scalar)(c->pv > 0 ? -m * c->pv : 0);
				for (j = 0; j < p.curve_count; j++)
				{
					qaws_curve *mv = NULL, *rv = NULL;
					unsigned int src = area > 0 ? j : p.curve_count - 1 - j;
					if ((st = qaws_curve_transform_2d(p.curves[src], t, &mv)) != QAWS_STATUS_OK)
						return st;
					if (area < 0)
					{
						/* a negative loop stays negative: its pieces run back */
						st = qaws_curve_reverse(mv, &rv);
						qaws_curve_destroy(mv);
						if (st != QAWS_STATUS_OK)
							return st;
						mv = rv;
					}
					if ((st = sr_keep(pl, mv)) != QAWS_STATUS_OK)
						return st;
					moved[j] = mv;
				}
				if ((st = sr_add_path(pl, moved, p.curve_count)) != QAWS_STATUS_OK)
					return st;
			}
		}
	return QAWS_STATUS_OK;
}

/* an operand's paths as loops of one period: loops going once round a
   closed direction capped along the domain's edge (u loops down to v0, v
   loops left to u0), their constant restored by whole periods */
static qaws_status sr_period(qaws_internal_sg_ctx const* c, qaws_path_2d const* in, unsigned int n, sr_pool* pl)
{
	unsigned int i;
	int wu = 0, wv = 0;
	qaws_status st = QAWS_STATUS_OK;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
	{
		qaws_path_2d const* p = &in[i];
		qaws_eval_result_2d a, b;
		qaws_range r0, r1;
		double du, dv;
		int ku = 0, kv = 0;
		if (p->curve_count == 0)
			continue;
		r0 = qaws_curve_get_parameter_range(p->curves[0]);
		r1 = qaws_curve_get_parameter_range(p->curves[p->curve_count - 1]);
		qaws_curve_evaluate_2d(p->curves[0], r0.min_value, QAWS_EVAL_FLAG_POSITION, &a);
		qaws_curve_evaluate_2d(p->curves[p->curve_count - 1], r1.max_value, QAWS_EVAL_FLAG_POSITION, &b);
		du = (double)b.position.x - a.position.x;
		dv = (double)b.position.y - a.position.y;
		if (c->pu > 0) ku = (int)floor(du / c->pu + 0.5);
		if (c->pv > 0) kv = (int)floor(dv / c->pv + 0.5);
		if (ku == 0 && kv == 0)
			st = sr_cut_loop(c, pl, p);
		else if ((ku == 1 || ku == -1) && kv == 0 && c->pv <= 0)
		{
			/* round u once: capped down to v0 */
			qaws_curve const* cv[1024];
			qaws_path_2d q;
			unsigned int j;
			if (p->curve_count > 1020)
				return QAWS_STATUS_UNSUPPORTED_OPERATION;
			for (j = 0; j < p->curve_count; j++) cv[j] = p->curves[j];
			if ((st = sr_line(pl, b.position.x, b.position.y, b.position.x, c->v0, &cv[j])) != QAWS_STATUS_OK) return st;
			if ((st = sr_line(pl, b.position.x, c->v0, a.position.x, c->v0, &cv[j + 1])) != QAWS_STATUS_OK) return st;
			if ((st = sr_line(pl, a.position.x, c->v0, a.position.x, a.position.y, &cv[j + 2])) != QAWS_STATUS_OK) return st;
			q.curves = cv; q.curve_count = p->curve_count + 3; q.closed = 1;
			st = sr_cut_loop(c, pl, &q);
			wu += ku;
		}
		else if ((kv == 1 || kv == -1) && ku == 0 && c->pu <= 0)
		{
			/* round v once: capped left to u0 */
			qaws_curve const* cv[1024];
			qaws_path_2d q;
			unsigned int j;
			if (p->curve_count > 1020)
				return QAWS_STATUS_UNSUPPORTED_OPERATION;
			for (j = 0; j < p->curve_count; j++) cv[j] = p->curves[j];
			if ((st = sr_line(pl, b.position.x, b.position.y, c->u0, b.position.y, &cv[j])) != QAWS_STATUS_OK) return st;
			if ((st = sr_line(pl, c->u0, b.position.y, c->u0, a.position.y, &cv[j + 1])) != QAWS_STATUS_OK) return st;
			if ((st = sr_line(pl, c->u0, a.position.y, a.position.x, a.position.y, &cv[j + 2])) != QAWS_STATUS_OK) return st;
			q.curves = cv; q.curve_count = p->curve_count + 3; q.closed = 1;
			st = sr_cut_loop(c, pl, &q);
			wv += kv;
		}
		else
			return QAWS_STATUS_UNSUPPORTED_OPERATION;      /* round several times, or round both ways */
	}
	/* the constant: a loop round u capped down counts -1 below it, while its
	   inside is above (its left) */
	for (i = 0; st == QAWS_STATUS_OK && (int)i < (wu > 0 ? wu : 0) + (wv < 0 ? -wv : 0); i++)
	{
		qaws_curve const* rc[4];
		double ua = c->u0, ub = c->pu > 0 ? c->u0 + c->pu : c->u1, va = c->v0, vb = c->pv > 0 ? c->v0 + c->pv : c->v1;
		if ((st = sr_line(pl, ua, va, ub, va, &rc[0])) != QAWS_STATUS_OK) return st;
		if ((st = sr_line(pl, ub, va, ub, vb, &rc[1])) != QAWS_STATUS_OK) return st;
		if ((st = sr_line(pl, ub, vb, ua, vb, &rc[2])) != QAWS_STATUS_OK) return st;
		if ((st = sr_line(pl, ua, vb, ua, va, &rc[3])) != QAWS_STATUS_OK) return st;
		st = sr_add_path(pl, rc, 4);
	}
	return st;
}

qaws_status qaws_surface_clip(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	qaws_surface_region** out_region)
{
	qaws_surface_region* r;
	qaws_path_2d *ps = NULL, *pc = NULL;
	qaws_curve const **ks = NULL, **kc = NULL;
	qaws_clip_desc d;
	qaws_status st, st2 = QAWS_STATUS_OK;
	if (!surface || !out_region || (!subjects && subject_count) || (!clips && clip_count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_region = NULL;
	r = (qaws_surface_region*)calloc(1, sizeof(*r));
	if (!r)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	r->s = surface;
	if ((st = qaws_internal_sg_init(&r->ctx, surface, 0)) != QAWS_STATUS_OK)
	{
		free(r);
		return st;
	}
	ps = sr_paths(surface, subjects, subject_count, &ks, &st);
	if (st == QAWS_STATUS_OK && clip_count)
		pc = sr_paths(surface, clips, clip_count, &kc, &st2);
	if (st == QAWS_STATUS_OK) st = st2;
	if (st == QAWS_STATUS_OK && r->ctx.pu <= 0 && r->ctx.pv <= 0)
	{
		memset(&d, 0, sizeof(d));
		d.subjects = ps; d.subject_count = subject_count;
		d.clips = pc; d.clip_count = clip_count;
		d.clip_type = clip_count && clip_type != QAWS_CLIP_NONE ? clip_type : QAWS_CLIP_UNION;
		d.fill_rule = fill_rule;
		st = qaws_clip_execute(&d, &r->res);
	}
	else if (st == QAWS_STATUS_OK)
	{
		/* each operand in one period under its rule, then the Boolean there */
		sr_pool pa, pb;
		qaws_clip_result *ua = NULL, *ub = NULL;
		memset(&pa, 0, sizeof(pa));
		memset(&pb, 0, sizeof(pb));
		st = sr_period(&r->ctx, ps, subject_count, &pa);
		if (st == QAWS_STATUS_OK && clip_count)
			st = sr_period(&r->ctx, pc, clip_count, &pb);
		sr_fix_views(&pa);
		sr_fix_views(&pb);
		if (st == QAWS_STATUS_OK)
			st = qaws_clip_boolean(QAWS_CLIP_UNION, fill_rule, pa.paths, pa.np, NULL, 0, &ua);
		if (st == QAWS_STATUS_OK && clip_count)
			st = qaws_clip_boolean(QAWS_CLIP_UNION, fill_rule, pb.paths, pb.np, NULL, 0, &ub);
		if (st == QAWS_STATUS_OK)
		{
			unsigned int na = qaws_clip_result_get_path_count(ua), nb = ub ? qaws_clip_result_get_path_count(ub) : 0, i;
			qaws_path_2d* qa = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (na + 1));
			qaws_path_2d* qb = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (nb + 1));
			if (!qa || !qb)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; i < na && st == QAWS_STATUS_OK; i++) st = qaws_clip_result_get_path(ua, i, &qa[i]);
			for (i = 0; i < nb && st == QAWS_STATUS_OK; i++) st = qaws_clip_result_get_path(ub, i, &qb[i]);
			if (st == QAWS_STATUS_OK)
			{
				memset(&d, 0, sizeof(d));
				d.subjects = qa; d.subject_count = na;
				d.clips = qb; d.clip_count = nb;
				d.clip_type = clip_count && clip_type != QAWS_CLIP_NONE ? clip_type : QAWS_CLIP_UNION;
				d.fill_rule = QAWS_FILL_NON_ZERO;
				st = qaws_clip_execute(&d, &r->res);
			}
			free(qa); free(qb);
		}
		qaws_clip_result_destroy(ua);
		qaws_clip_result_destroy(ub);
		sr_pool_free(&pa);
		sr_pool_free(&pb);
	}
	free(ps); free(pc); free(ks); free(kc);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_region_destroy(r);
		return st;
	}
	*out_region = r;
	return QAWS_STATUS_OK;
}

void qaws_surface_region_destroy(qaws_surface_region* region)
{
	if (!region)
		return;
	qaws_clip_result_destroy(region->res);
	free(region);
}

unsigned int qaws_surface_region_get_path_count(qaws_surface_region const* region)
{
	return region ? qaws_clip_result_get_path_count(region->res) : 0;
}

qaws_status qaws_surface_region_get_path(qaws_surface_region const* region, unsigned int index, qaws_path_2d* out_uv_path)
{
	if (!region || !out_uv_path)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return qaws_clip_result_get_path(region->res, index, out_uv_path);
}

int qaws_surface_region_is_hole(qaws_surface_region const* region, unsigned int index)
{
	return region ? qaws_clip_result_is_hole(region->res, index) : 0;
}

/* P(u, v): the integral of h along u from u_ref; dim values */
typedef struct sr_job
{
	qaws_surface_region const* r;
	qaws_surface_region_fn fn;
	void* user;
	unsigned int dim;
	double u_ref, panel;
} sr_job;

static qaws_status sr_density(sr_job const* jb, double u, double v, double* out)
{
	qaws_surface_eval_result e;
	double a;
	unsigned int k;
	qaws_status st = qaws_internal_sg_area_element(&jb->r->ctx, u, v, &a);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!jb->fn)
	{
		out[0] = a;
		return QAWS_STATUS_OK;
	}
	{
		qaws_scalar vals[16];
		st = qaws_internal_sg_eval(&jb->r->ctx, u, v, QAWS_SURFACE_EVAL_POSITION, &e);
		if (st != QAWS_STATUS_OK)
			return st;
		jb->fn(jb->user, (qaws_scalar)u, (qaws_scalar)v, e.position, vals);
		for (k = 0; k < jb->dim; k++)
			out[k] = (double)vals[k] * a;
	}
	return QAWS_STATUS_OK;
}

static qaws_status sr_potential(sr_job const* jb, double u, double v, double* out)
{
	double len = u - jb->u_ref, t[16];
	unsigned int m, p, i, k, d = jb->dim;
	for (k = 0; k < d; k++) out[k] = 0;
	if (len == 0)
		return QAWS_STATUS_OK;
	m = (unsigned int)ceil(fabs(len) / jb->panel);
	if (m < 1) m = 1;
	for (p = 0; p < m; p++)
	{
		double a = jb->u_ref + len * p / m, b = jb->u_ref + len * (p + 1) / m;
		for (i = 0; i < 4; i++)
		{
			unsigned int side;
			for (side = 0; side < 2; side++)
			{
				double s = 0.5 * (a + b) + (side ? 0.5 : -0.5) * (b - a) * sr_x8[i];
				qaws_status st = sr_density(jb, s, v, t);
				if (st != QAWS_STATUS_OK)
					return st;
				for (k = 0; k < d; k++)
					out[k] += 0.5 * (b - a) * sr_w8[i] * t[k];
			}
		}
	}
	return QAWS_STATUS_OK;
}

static qaws_status sr_green(sr_job const* jb, double* out)
{
	unsigned int np = qaws_clip_result_get_path_count(jb->r->res), pi, ci, k, d = jb->dim;
	for (k = 0; k < d; k++) out[k] = 0;
	for (pi = 0; pi < np; pi++)
	{
		qaws_path_2d path;
		qaws_status st = qaws_clip_result_get_path(jb->r->res, pi, &path);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ci = 0; ci < path.curve_count; ci++)
		{
			qaws_curve const* c = path.curves[ci];
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int spans = qaws_curve_get_span_count(c), sp, i, side;
			double t0 = (double)rg.min_value, t1 = (double)rg.max_value;
			if (spans < 1) spans = 1;
			for (sp = 0; sp < 2 * spans; sp++)
			{
				double a = t0 + (t1 - t0) * sp / (2 * spans), b = t0 + (t1 - t0) * (sp + 1) / (2 * spans);
				for (i = 0; i < 4; i++)
					for (side = 0; side < 2; side++)
					{
						double t = 0.5 * (a + b) + (side ? 0.5 : -0.5) * (b - a) * sr_x8[i], P[16];
						qaws_eval_result_2d e;
						st = qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
						if (st != QAWS_STATUS_OK)
							return st;
						st = sr_potential(jb, (double)e.position.x, (double)e.position.y, P);
						if (st != QAWS_STATUS_OK)
							return st;
						for (k = 0; k < d; k++)
							out[k] += 0.5 * (b - a) * sr_w8[i] * P[k] * (double)e.d1.y;
					}
			}
		}
	}
	return QAWS_STATUS_OK;
}

static void sr_job_init(sr_job* jb, qaws_surface_region const* r)
{
	memset(jb, 0, sizeof(*jb));
	jb->r = r;
	jb->dim = 1;
	jb->u_ref = r->ctx.u0;
	jb->panel = (r->ctx.u1 - r->ctx.u0) / 8;
}

qaws_status qaws_surface_region_get_area(qaws_surface_region const* region, qaws_scalar* out_area)
{
	sr_job jb;
	double a;
	qaws_status st;
	if (!region || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;
	sr_job_init(&jb, region);
	st = sr_green(&jb, &a);
	*out_area = (qaws_scalar)a;
	return st;
}

qaws_status qaws_surface_region_integrate(qaws_surface_region const* region, qaws_surface_region_fn fn, void* user,
	unsigned int dim, qaws_scalar* out_values)
{
	sr_job jb;
	double v[16];
	unsigned int k;
	qaws_status st;
	if (!region || !fn || !out_values || dim == 0 || dim > 16)
		return QAWS_STATUS_INVALID_ARGUMENT;
	sr_job_init(&jb, region);
	jb.fn = fn; jb.user = user; jb.dim = dim;
	st = sr_green(&jb, v);
	for (k = 0; k < dim; k++)
		out_values[k] = (qaws_scalar)v[k];
	return st;
}

/* a boundary piece along a seam line (a cut, not the region's edge) */
static int sr_on_seam(qaws_internal_sg_ctx const* c, qaws_curve const* cv, double ta, double tb)
{
	unsigned int k, j;
	for (j = 0; j < 2; j++)
	{
		double line = j ? c->v0 : c->u0, per = j ? c->pv : c->pu;
		int all = per > 0;
		for (k = 0; k <= 4 && all; k++)
		{
			qaws_eval_result_2d e;
			double x, t;
			qaws_curve_evaluate_2d(cv, (qaws_scalar)(ta + (tb - ta) * k / 4), QAWS_EVAL_FLAG_POSITION, &e);
			x = j ? (double)e.position.y : (double)e.position.x;
			t = fmod(x - line, per);
			if (t < 0) t += per;
			if (!(t <= 1e3 * c->tol || per - t <= 1e3 * c->tol))
				all = 0;
		}
		if (all)
			return 1;
	}
	return 0;
}

qaws_status qaws_surface_region_get_perimeter(qaws_surface_region const* region, qaws_scalar* out_length)
{
	unsigned int np, pi, ci;
	double total = 0;
	if (!region || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;
	np = qaws_clip_result_get_path_count(region->res);
	for (pi = 0; pi < np; pi++)
	{
		qaws_path_2d path;
		qaws_status st = qaws_clip_result_get_path(region->res, pi, &path);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ci = 0; ci < path.curve_count; ci++)
		{
			qaws_curve const* c = path.curves[ci];
			qaws_range rg = qaws_curve_get_parameter_range(c);
			unsigned int spans = qaws_curve_get_span_count(c), sp, i, side;
			double t0 = (double)rg.min_value, t1 = (double)rg.max_value;
			if (spans < 1) spans = 1;
			for (sp = 0; sp < 2 * spans; sp++)
			{
				double a = t0 + (t1 - t0) * sp / (2 * spans), b = t0 + (t1 - t0) * (sp + 1) / (2 * spans);
				if (sr_on_seam(&region->ctx, c, a, b))
					continue;          /* a cut along a seam, not the region's edge */
				for (i = 0; i < 4; i++)
					for (side = 0; side < 2; side++)
					{
						double t = 0.5 * (a + b) + (side ? 0.5 : -0.5) * (b - a) * sr_x8[i], g[3];
						qaws_eval_result_2d e;
						st = qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
						if (st != QAWS_STATUS_OK)
							return st;
						st = qaws_internal_sg_first_form(&region->ctx, (double)e.position.x, (double)e.position.y, g);
						if (st != QAWS_STATUS_OK)
							return st;
						total += 0.5 * (b - a) * sr_w8[i] * sqrt(g[0] * e.d1.x * e.d1.x + 2 * g[1] * e.d1.x * e.d1.y + g[2] * e.d1.y * e.d1.y);
					}
			}
		}
	}
	*out_length = (qaws_scalar)total;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_region_contains(qaws_surface_region const* region, qaws_vec2 uv, int* out_inside)
{
	unsigned int np, i;
	qaws_path_2d* paths;
	qaws_point_location loc = QAWS_POINT_OUTSIDE;
	qaws_status st = QAWS_STATUS_OK;
	if (!region || !out_inside)
		return QAWS_STATUS_INVALID_ARGUMENT;
	np = qaws_clip_result_get_path_count(region->res);
	*out_inside = 0;
	if (np == 0)
		return QAWS_STATUS_OK;
	paths = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * np);
	if (!paths)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
		st = qaws_clip_result_get_path(region->res, i, &paths[i]);
	if (region->ctx.pu > 0)
	{
		double t = fmod((double)uv.x - region->ctx.u0, region->ctx.pu);
		if (t < 0) t += region->ctx.pu;
		uv.x = (qaws_scalar)(region->ctx.u0 + t);
	}
	if (region->ctx.pv > 0)
	{
		double t = fmod((double)uv.y - region->ctx.v0, region->ctx.pv);
		if (t < 0) t += region->ctx.pv;
		uv.y = (qaws_scalar)(region->ctx.v0 + t);
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_region_locate_point_2d(paths, np, QAWS_FILL_NON_ZERO, uv, (qaws_scalar)region->ctx.tol, &loc);
	free(paths);
	*out_inside = loc == QAWS_POINT_ON ? -1 : (loc == QAWS_POINT_INSIDE ? 1 : 0);
	return st;
}

qaws_status qaws_surface_region_create_trimmed(qaws_surface_region const* region, qaws_surface** out_surface)
{
	unsigned int np, i;
	qaws_trim_loop* loops;
	qaws_surface_trim_desc d;
	qaws_status st = QAWS_STATUS_OK;
	if (!region || !out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	np = qaws_clip_result_get_path_count(region->res);
	loops = (qaws_trim_loop*)malloc(sizeof(qaws_trim_loop) * (np + 1));
	if (!loops)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
	{
		qaws_path_2d p;
		st = qaws_clip_result_get_path(region->res, i, &p);
		loops[i].curves = p.curves;
		loops[i].curve_count = p.curve_count;
		loops[i].is_outer = !qaws_clip_result_is_hole(region->res, i);
	}
	if (st == QAWS_STATUS_OK)
	{
		memset(&d, 0, sizeof(d));
		d.base = region->s;
		d.loops = loops;
		d.loop_count = np;
		st = qaws_surface_create_trimmed(&d, out_surface);
	}
	free(loops);
	return st;
}
