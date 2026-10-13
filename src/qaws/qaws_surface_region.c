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
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface_region.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
	if (st == QAWS_STATUS_OK)
	{
		memset(&d, 0, sizeof(d));
		d.subjects = ps; d.subject_count = subject_count;
		d.clips = pc; d.clip_count = clip_count;
		d.clip_type = clip_count && clip_type != QAWS_CLIP_NONE ? clip_type : QAWS_CLIP_UNION;
		d.fill_rule = fill_rule;
		st = qaws_clip_execute(&d, &r->res);
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
	double g[3], a;
	unsigned int k;
	qaws_status st = qaws_internal_sg_metric(&jb->r->ctx, u, v, g, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	a = sqrt(g[0] * g[2] - g[1] * g[1]);
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
				for (i = 0; i < 4; i++)
					for (side = 0; side < 2; side++)
					{
						double t = 0.5 * (a + b) + (side ? 0.5 : -0.5) * (b - a) * sr_x8[i], g[3];
						qaws_eval_result_2d e;
						st = qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
						if (st != QAWS_STATUS_OK)
							return st;
						st = qaws_internal_sg_metric(&region->ctx, (double)e.position.x, (double)e.position.y, g, NULL);
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
