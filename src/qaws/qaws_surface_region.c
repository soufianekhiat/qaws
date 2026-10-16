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
#include "qaws_bspline.h"
#include "qaws_composite.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_region_sample.h"
#include "qaws_triangulate.h"
#include "qaws_clip_exact.h"
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

#if QAWS_SCALAR_IS_FLOAT
#define SR_EXACT_BITS 26
#define SR_EXACT_PARAM_BITS 24
#else
#define SR_EXACT_BITS 32   /* the kernel's widest coordinates */
#define SR_EXACT_PARAM_BITS 48   /* long splines: their pieces' ends a rounding off the true cuts */
#endif

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
	st = qaws_curve_create_bezier(&c, &d);
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

/* a piece cut at a period cell, cleaned for the Boolean that follows (a
   curve that touches or runs along a seam a few ulps off it meets its
   copies across the seam nearly tangentially, which the clip resolves
   badly):
   - an end span of no length (a cut a rounding away from a knot leaves
     one, its tangent noise) dropped where the knot joining it to the next
     span is C0 (multiplicity >= the degree), the end point kept exactly;
   - control points within rounding of the cell's lines put on them.
   Lines: lines[2j] is 0 for u = const, 1 for v = const, lines[2j + 1] the
   value. Kinds other than Bezier and B-spline are kept. */
static qaws_status sr_snap_cell(qaws_curve const* cv, double const* lines, unsigned int nl, qaws_curve** out)
{
	qaws_curve_kind kind = qaws_curve_get_kind(cv);
	unsigned int n = 0, nk = 0, deg = qaws_curve_get_degree(cv), i, j, f = 0, b = 0, changed = 0;
	qaws_scalar probe[2], *cp, *kn = NULL;
	qaws_status st;
#if QAWS_SCALAR_IS_FLOAT
	double const near_line = 4e-6, tiny = 2e-5;   /* tiny: past what snapping can fold to a point */
#else
	double const near_line = 1e-13, tiny = 1e-12;
#endif
	*out = NULL;
	if (kind != QAWS_CURVE_KIND_BEZIER && kind != QAWS_CURVE_KIND_BSPLINE)
		return QAWS_STATUS_OK;
	qaws_curve_get_control_points(probe, &n, cv, 0);
	if (kind == QAWS_CURVE_KIND_BSPLINE)
		qaws_bspline_get_knots(probe, &nk, cv, 0);
	cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (2 * n + nk + 2));
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if ((st = qaws_curve_get_control_points(cp, &n, cv, n)) == QAWS_STATUS_OK && nk)
	{
		kn = cp + 2 * n;
		st = qaws_bspline_get_knots(kn, &nk, cv, nk);
	}
	/* f, b: control points dropped at the front and the back */
	if (st == QAWS_STATUS_OK && kn && deg > 1 && nk == n + deg + 1)
	{
		unsigned int e;
		for (e = 0; e < 2; e++)
		{
			/* the end span's control points: deg + 1 from the end; its far knot's
			   multiplicity */
			unsigned int m = 1, k, base = e ? n - 1 - deg - b : f;
			double kfar = (double)kn[e ? nk - deg - 2 - b : deg + 1 + f], span = 0;
			if (n - f - b < 2 * deg + 1)
				break;
			for (k = 1; k < nk; k++)
			{
				unsigned int idx = e ? nk - deg - 2 - b - k : deg + 1 + f + k;
				if (idx >= nk || (double)kn[idx] != kfar) break;
				m++;
			}
			if (m < deg)
				continue;
			for (k = 0; k <= deg; k++)
			{
				double dx = (double)cp[2 * (base + k)] - cp[2 * (e ? n - 1 - b : f)];
				double dy = (double)cp[2 * (base + k) + 1] - cp[2 * (e ? n - 1 - b : f) + 1];
				if (hypot(dx, dy) > span) span = hypot(dx, dy);
			}
			if (span <= tiny * (1 + fabs((double)cp[2 * base]) + fabs((double)cp[2 * base + 1])))
			{
				if (e) b += deg; else f += deg;
			}
		}
	}
	if (f || b)
	{
		/* the new ends where the old ones were; the clamped knots */
		qaws_scalar s0x = cp[0], s0y = cp[1], s1x = cp[2 * n - 2], s1y = cp[2 * n - 1];
		unsigned int nn = n - f - b, nnk = nn + deg + 1;
		qaws_scalar* k2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nnk + 1));
		if (!k2)
		{
			free(cp);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < nnk; i++)
			k2[i] = kn[i + f];
		for (i = 0; i <= deg; i++)
		{
			k2[i] = kn[f + deg];
			k2[nnk - 1 - i] = kn[nk - 1 - b - deg];
		}
		memmove(cp, cp + 2 * f, sizeof(qaws_scalar) * 2 * nn);
		n = nn;
		kn = cp + 2 * n;
		memcpy(kn, k2, sizeof(qaws_scalar) * nnk);
		nk = nnk;
		free(k2);
		cp[0] = s0x; cp[1] = s0y;
		cp[2 * n - 2] = s1x; cp[2 * n - 1] = s1y;
		changed = 1;
	}
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
		for (j = 0; j < nl; j++)
		{
			unsigned int k = (unsigned int)lines[2 * j];
			double x = (double)cp[2 * i + k], l = lines[2 * j + 1];
			if (x != l && fabs(x - l) <= near_line * (1 + fabs(l)))
			{
				cp[2 * i + k] = (qaws_scalar)l;
				changed = 1;
			}
		}
	if (st == QAWS_STATUS_OK && changed)
	{
		if (kind == QAWS_CURVE_KIND_BEZIER)
		{
			qaws_bezier_desc d;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.degree = deg; d.control_points = cp; d.control_point_count = n;
			st = qaws_curve_create_bezier(out, &d);
		}
		else
		{
			qaws_bspline_desc d;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.degree = deg; d.control_points = cp; d.control_point_count = n;
			d.knots = kn; d.knot_count = nk;
			st = qaws_curve_create_bspline(out, &d);
		}
	}
	free(cp);
	return st;
}

/* the Booleans' path: float (qaws_clip) or exact topology (qaws_clip_exact
   on its lattice) */
typedef struct sr_mode
{
	int exact;
	qaws_exact_desc ex;
} sr_mode;

static qaws_status sr_bool(sr_mode const* md, qaws_clip_type type, qaws_fill_rule rule, qaws_path_2d const* a, unsigned int na,
	qaws_path_2d const* b, unsigned int nb, qaws_clip_result** out)
{
	qaws_clip_exact_desc d;
	if (!md || !md->exact)
		return qaws_clip_boolean(type, rule, a, na, b, nb, out);
	memset(&d, 0, sizeof(d));
	d.subjects = a; d.subject_count = na;
	d.clips = b; d.clip_count = nb;
	d.clip_type = nb && type != QAWS_CLIP_NONE ? type : QAWS_CLIP_UNION;
	d.fill_rule = rule;
	d.exact = &md->ex;
	return qaws_clip_exact_execute(&d, out, NULL);
}

/* the loop's part in each period cell, translated into the first, its
   orientation kept */
static qaws_status sr_cut_loop(qaws_internal_sg_ctx const* c, sr_pool* pl, qaws_path_2d const* loop, sr_mode const* md)
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
	if (k0 == k1 && m0 == m1)
	{
		/* in one cell: the loop itself, moved into the first (no clip, so an
		   exact path takes no pieces back as input) */
		qaws_curve const* moved[1024];
		qaws_scalar t[6];
		unsigned int j;
		if (loop->curve_count > 1024)
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		t[0] = 1; t[1] = 0; t[2] = (qaws_scalar)(c->pu > 0 ? -k0 * c->pu : 0);
		t[3] = 0; t[4] = 1; t[5] = (qaws_scalar)(c->pv > 0 ? -m0 * c->pv : 0);
		for (j = 0; j < loop->curve_count; j++)
		{
			qaws_curve* mv = NULL;
			if ((st = qaws_curve_transform_2d(loop->curves[j], t, &mv)) != QAWS_STATUS_OK)
				return st;
			if ((st = sr_keep(pl, mv)) != QAWS_STATUS_OK)
				return st;
			moved[j] = mv;
		}
		return sr_add_path(pl, moved, loop->curve_count);
	}
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
			st = sr_bool(md, QAWS_CLIP_INTERSECTION, area > 0 ? QAWS_FILL_POSITIVE : QAWS_FILL_NEGATIVE, loop, 1, &cell, 1, &res);
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
					qaws_curve* snapped = NULL;
					double lines[8];
					unsigned int nl = 0;
					if (c->pu > 0) { lines[2 * nl] = 0; lines[2 * nl + 1] = ua; nl++; lines[2 * nl] = 0; lines[2 * nl + 1] = ub; nl++; }
					if (c->pv > 0) { lines[2 * nl] = 1; lines[2 * nl + 1] = va; nl++; lines[2 * nl] = 1; lines[2 * nl + 1] = vb; nl++; }
					if ((!md || !md->exact) && (st = sr_snap_cell(p.curves[src], lines, nl, &snapped)) != QAWS_STATUS_OK)
						return st;
					st = qaws_curve_transform_2d(snapped ? snapped : p.curves[src], t, &mv);
					qaws_curve_destroy(snapped);
					if (st != QAWS_STATUS_OK)
						return st;
					if (area < 0)
					{
						/* a negative loop stays negative: its pieces run back */
						st = qaws_curve_reverse(&rv, mv);
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
static qaws_status sr_period(qaws_internal_sg_ctx const* c, qaws_path_2d const* in, unsigned int n, sr_pool* pl, sr_mode const* md)
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
		qaws_curve_evaluate_2d(&a, p->curves[0], r0.min_value, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&b, p->curves[p->curve_count - 1], r1.max_value, QAWS_EVAL_FLAG_POSITION);
		du = (double)b.position.x - a.position.x;
		dv = (double)b.position.y - a.position.y;
		if (c->pu > 0) ku = (int)floor(du / c->pu + 0.5);
		if (c->pv > 0) kv = (int)floor(dv / c->pv + 0.5);
		if (ku == 0 && kv == 0)
			st = sr_cut_loop(c, pl, p, md);
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
			st = sr_cut_loop(c, pl, &q, md);
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
			st = sr_cut_loop(c, pl, &q, md);
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

static qaws_status sr_clip(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	sr_mode const* md, qaws_surface_region** out_region)
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
		st = sr_bool(md, d.clip_type, d.fill_rule, d.subjects, d.subject_count, d.clips, d.clip_count, &r->res);
	}
	else if (st == QAWS_STATUS_OK)
	{
		/* each operand in one period under its rule, then the Boolean there */
		sr_pool pa, pb;
		qaws_clip_result *ua = NULL, *ub = NULL;
		memset(&pa, 0, sizeof(pa));
		memset(&pb, 0, sizeof(pb));
		st = sr_period(&r->ctx, ps, subject_count, &pa, md);
		if (st == QAWS_STATUS_OK && clip_count)
			st = sr_period(&r->ctx, pc, clip_count, &pb, md);
		sr_fix_views(&pa);
		sr_fix_views(&pb);
		if (st == QAWS_STATUS_OK)
			st = sr_bool(md, QAWS_CLIP_UNION, fill_rule, pa.paths, pa.np, NULL, 0, &ua);
		if (st == QAWS_STATUS_OK && clip_count)
			st = sr_bool(md, QAWS_CLIP_UNION, fill_rule, pb.paths, pb.np, NULL, 0, &ub);
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
				st = sr_bool(md, d.clip_type, d.fill_rule, d.subjects, d.subject_count, d.clips, d.clip_count, &r->res);
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

qaws_status qaws_surface_clip(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	qaws_surface_region** out_region)
{
	return sr_clip(surface, clip_type, fill_rule, subjects, subject_count, clips, clip_count, NULL, out_region);
}

qaws_status qaws_surface_clip_exact(qaws_surface const* surface, qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_surface_path const* subjects, unsigned int subject_count,
	qaws_surface_path const* clips, unsigned int clip_count,
	qaws_exact_desc const* exact, qaws_surface_region** out_region)
{
	sr_mode md;
	md.exact = 1;
	if (exact)
		md.ex = *exact;
	else
	{
		/* the lattice over the domain and a period either side of it, fine
		   enough for the (u, v) images' tolerance */
		qaws_internal_sg_ctx c;
		qaws_status st;
		double ext;
		int e;
		if (!surface)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if ((st = qaws_internal_sg_init(&c, surface, 0)) != QAWS_STATUS_OK)
			return st;
		ext = fabs(c.u0) > fabs(c.u1) ? fabs(c.u0) : fabs(c.u1);
		if (fabs(c.v0) > ext) ext = fabs(c.v0);
		if (fabs(c.v1) > ext) ext = fabs(c.v1);
		ext += (c.pu > c.pv ? c.pu : c.pv) + 0.5 * c.diag;
		qaws_exact_desc_default(&md.ex);
		md.ex.coord_bits = SR_EXACT_BITS;
		md.ex.param_bits = SR_EXACT_PARAM_BITS;
		e = (int)ceil(log2(ext > 0 ? ext : 1));
		md.ex.space_exp2 = e - (int)md.ex.coord_bits;
	}
	return sr_clip(surface, clip_type, fill_rule, subjects, subject_count, clips, clip_count, &md, out_region);
}

qaws_surface const* qaws_surface_region_get_surface(qaws_surface_region const* region)
{
	return region ? region->s : NULL;
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
						st = qaws_curve_evaluate_2d(&e, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
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
			qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(ta + (tb - ta) * k / 4), QAWS_EVAL_FLAG_POSITION);
			x = j ? (double)e.position.y : (double)e.position.x;
			t = fmod(x - line, per);
			if (t < 0) t += per;
			/* the cuts are the cell's lines themselves: within rounding */
#if QAWS_SCALAR_IS_FLOAT
			if (!(t <= 1e-5 * (1 + fabs(line) + per) || per - t <= 1e-5 * (1 + fabs(line) + per)))
#else
			if (!(t <= 1e-12 * (1 + fabs(line) + per) || per - t <= 1e-12 * (1 + fabs(line) + per)))
#endif
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
						st = qaws_curve_evaluate_2d(&e, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
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
		st = qaws_surface_create_trimmed(out_surface, &d);
	}
	free(loops);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  Offsets                                                                 */
/* ------------------------------------------------------------------------ */

/* the curves made for an offset, and the corners already given a circle */
typedef struct sr_tube
{
	qaws_surface const* s;
	qaws_internal_sg_ctx const* c;
	double r, tol;
	sr_pool raw;
	qaws_surface_curve** disks;
	unsigned int nd, capd;
	double* centres;
	unsigned int nc, capc;
} sr_tube;

/* the strip swept by the conormal geodesics of length |sd| on the side of
   sd (> 0 left), counter-clockwise */
static qaws_status sr_strip(sr_tube* tb, qaws_curve const* base, double sd)
{
	sr_pool* pl = &tb->raw;
	qaws_curve *off = NULL, *r0 = NULL, *r1 = NULL, *x = NULL, *y = NULL;
	qaws_curve const* loop[4];
	qaws_status st = qaws_internal_sg_offset(tb->c, base, sd, &off, &r0, &r1);
	if (st != QAWS_STATUS_OK)
		return st;
	if ((st = sr_keep(pl, off)) != QAWS_STATUS_OK) { qaws_curve_destroy(r0); qaws_curve_destroy(r1); return st; }
	if ((st = sr_keep(pl, r0)) != QAWS_STATUS_OK) { qaws_curve_destroy(r1); return st; }
	if ((st = sr_keep(pl, r1)) != QAWS_STATUS_OK) return st;
	if (sd > 0)
	{
		/* along the base, out to the offset's end, back along it, in */
		if ((st = qaws_curve_reverse(&x, off)) != QAWS_STATUS_OK || (st = sr_keep(pl, x)) != QAWS_STATUS_OK) return st;
		if ((st = qaws_curve_reverse(&y, r0)) != QAWS_STATUS_OK || (st = sr_keep(pl, y)) != QAWS_STATUS_OK) return st;
		loop[0] = base; loop[1] = r1; loop[2] = x; loop[3] = y;
	}
	else
	{
		/* along the offset, in to the base's end, back along it, out */
		if ((st = qaws_curve_reverse(&x, r1)) != QAWS_STATUS_OK || (st = sr_keep(pl, x)) != QAWS_STATUS_OK) return st;
		if ((st = qaws_curve_reverse(&y, base)) != QAWS_STATUS_OK || (st = sr_keep(pl, y)) != QAWS_STATUS_OK) return st;
		loop[0] = off; loop[1] = x; loop[2] = y; loop[3] = r0;
	}
	return sr_add_path(pl, loop, 4);
}

/* a geodesic circle of radius r about the point (once per point: wrapped
   into the period and compared within the tolerance) */
static qaws_status sr_disk(sr_tube* tb, qaws_vec2 p)
{
	qaws_surface_curve* d = NULL;
	qaws_curve const* one[1];
	double u = (double)p.x, v = (double)p.y;
	unsigned int i;
	qaws_status st;
	if (tb->c->pu > 0) { u = fmod(u - tb->c->u0, tb->c->pu); if (u < 0) u += tb->c->pu; }
	if (tb->c->pv > 0) { v = fmod(v - tb->c->v0, tb->c->pv); if (v < 0) v += tb->c->pv; }
	for (i = 0; i < tb->nc; i++)
	{
		double du = fabs(tb->centres[2 * i] - u), dv = fabs(tb->centres[2 * i + 1] - v);
		if (tb->c->pu > 0 && du > 0.5 * tb->c->pu) du = tb->c->pu - du;
		if (tb->c->pv > 0 && dv > 0.5 * tb->c->pv) dv = tb->c->pv - dv;
		if (du <= tb->c->tol && dv <= tb->c->tol)
			return QAWS_STATUS_OK;
	}
	if (!sr_grow((void**)&tb->centres, &tb->capc, 2 * (tb->nc + 1), sizeof(double)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	tb->centres[2 * tb->nc] = u; tb->centres[2 * tb->nc + 1] = v; tb->nc++;
	st = qaws_surface_curve_create_geodesic_arc(tb->s, p, (qaws_scalar)tb->r, 0, (qaws_scalar)(2 * 3.14159265358979323846), (qaws_scalar)tb->tol, &d);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!sr_grow((void**)&tb->disks, &tb->capd, tb->nd + 1, sizeof(qaws_surface_curve*)))
	{
		qaws_surface_curve_destroy(d);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	tb->disks[tb->nd++] = d;
	one[0] = qaws_surface_curve_get_uv(d);
	return sr_add_path(&tb->raw, one, 1);
}

/* a piece of no length on the surface (along a pole) */
static int sr_on_pole(qaws_internal_sg_ctx const* c, qaws_curve const* cv, double ta, double tb)
{
	unsigned int k;
	for (k = 1; k <= 3; k++)
	{
		qaws_eval_result_2d e;
		double g[3], d2;
		qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(ta + (tb - ta) * k / 4), QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		if (qaws_internal_sg_first_form(c, (double)e.position.x, (double)e.position.y, g) != QAWS_STATUS_OK)
			return 0;
		d2 = (double)e.d1.x * e.d1.x + (double)e.d1.y * e.d1.y;
		if (!(g[0] * e.d1.x * e.d1.x + 2 * g[1] * e.d1.x * e.d1.y + g[2] * e.d1.y * e.d1.y <= 1e-16 * (g[0] + g[2]) * d2))
			return 0;
	}
	return 1;
}

/* a corner at t: the tangents' directions just before and after it apart */
static int sr_kink(qaws_curve const* cv, double t, qaws_range rg)
{
#if QAWS_SCALAR_IS_FLOAT
	double const eps = 1e-4, ang = 1e-2;
#else
	double const eps = 1e-9, ang = 1e-5;
#endif
	double h = eps * ((double)rg.max_value - rg.min_value), la, lb, cr, dt;
	qaws_eval_result_2d a, b;
	qaws_curve_evaluate_2d(&a, cv, (qaws_scalar)(t - h), QAWS_EVAL_FLAG_D1);
	qaws_curve_evaluate_2d(&b, cv, (qaws_scalar)(t + h), QAWS_EVAL_FLAG_D1);
	la = hypot((double)a.d1.x, (double)a.d1.y);
	lb = hypot((double)b.d1.x, (double)b.d1.y);
	if (!(la > 0) || !(lb > 0))
		return 1;
	cr = ((double)a.d1.x * b.d1.y - (double)a.d1.y * b.d1.x) / (la * lb);
	dt = ((double)a.d1.x * b.d1.x + (double)a.d1.y * b.d1.y) / (la * lb);
	return atan2(fabs(cr), dt) > ang;
}

/* a piece [lo, hi]'s end: its point, and its tangent from inside the piece
   (at a knot the curve's own derivative is the next piece's) */
static void sr_piece_end(qaws_curve const* cv, double lo, double hi, qaws_eval_result_2d* r)
{
	qaws_eval_result_2d d;
#if QAWS_SCALAR_IS_FLOAT
	double const eps = 1e-4;
#else
	double const eps = 1e-9;
#endif
	qaws_curve_evaluate_2d(r, cv, (qaws_scalar)hi, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
	qaws_curve_evaluate_2d(&d, cv, (qaws_scalar)(hi - eps * (hi - lo)), QAWS_EVAL_FLAG_D1);
	r->d1 = d.d1;
}

/* the parameters splitting the curve into smooth pieces: its ends and the
   knots of a B-spline where it turns a corner (malloc'ed) */
static qaws_status sr_breaks(qaws_curve const* cv, double** out, unsigned int* out_n)
{
	qaws_range rg = qaws_curve_get_parameter_range(cv);
	qaws_curve_kind kind = qaws_curve_get_kind(cv);
	double* br = NULL;
	unsigned int nb = 0, i;
	*out = NULL;
	*out_n = 0;
	br = (double*)malloc(sizeof(double) * 2);
	if (!br)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	br[nb++] = (double)rg.min_value;
	if (kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS)
	{
		unsigned int nk = 0, deg = qaws_curve_get_degree(cv);
		qaws_scalar probe[2];
		qaws_bspline_get_knots(probe, &nk, cv, 0);
		if (nk > 0)
		{
			qaws_scalar* kn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nk);
			double* nbr = (double*)realloc(br, sizeof(double) * (nk + 2));
			if (!kn || !nbr)
			{
				free(kn); free(nbr ? nbr : br);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			br = nbr;
			if (qaws_bspline_get_knots(kn, &nk, cv, nk) == QAWS_STATUS_OK)
				for (i = 0; i < nk;)
				{
					unsigned int m = 1;
					while (i + m < nk && kn[i + m] == kn[i]) m++;
					if (m >= deg && (double)kn[i] > br[nb - 1] && (double)kn[i] < (double)rg.max_value && sr_kink(cv, (double)kn[i], rg))
						br[nb++] = (double)kn[i];
					i += m;
				}
			free(kn);
		}
	}
	br[nb++] = (double)rg.max_value;
	*out = br;
	*out_n = nb;
	return QAWS_STATUS_OK;
}

/* the curve in smooth pieces (sr_breaks), less its cuts along seams and
   poles: each piece's strips on the sides asked (bit 1 left, bit 2
   right), circles at its ends */
static qaws_status sr_tube_curve(sr_tube* tb, qaws_curve const* cv, int sides, int disks)
{
	double* br = NULL;
	unsigned int nb = 0, i;
	qaws_status st = sr_breaks(cv, &br, &nb);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i + 1 < nb && st == QAWS_STATUS_OK; i++)
	{
		double ta = br[i], tb2 = br[i + 1];
		qaws_curve* piece = NULL;
		qaws_curve const* base = cv;
		qaws_eval_result_2d e0, e1;
		if (!(tb2 > ta) || sr_on_seam(tb->c, cv, ta, tb2) || sr_on_pole(tb->c, cv, ta, tb2))
			continue;
		if (nb > 2)
		{
			if ((st = qaws_curve_extract(&piece, cv, (qaws_scalar)ta, (qaws_scalar)tb2)) != QAWS_STATUS_OK || (st = sr_keep(&tb->raw, piece)) != QAWS_STATUS_OK)
				break;
			base = piece;
		}
		if (sides & 1) st = sr_strip(tb, base, tb->r);
		if (st == QAWS_STATUS_OK && (sides & 2)) st = sr_strip(tb, base, -tb->r);
		if (st == QAWS_STATUS_OK && disks)
		{
			qaws_range pr = qaws_curve_get_parameter_range(base);
			qaws_curve_evaluate_2d(&e0, base, pr.min_value, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&e1, base, pr.max_value, QAWS_EVAL_FLAG_POSITION);
			st = sr_disk(tb, e0.position);
			if (st == QAWS_STATUS_OK) st = sr_disk(tb, e1.position);
		}
	}
	free(br);
	return st;
}

/* the circles a region's corners need: where its boundary's smooth pieces
   meet turning left (convex, inside on the left) when it grows, right
   (reflex) when it shrinks; none at smooth joins. A piece ending at a cut
   along a seam goes on from the piece leaving the same point on the
   seam's other side (another loop, maybe): matched there by the point
   modulo the periods. */
typedef struct sr_end
{
	double p[2], d[2];          /* the point and the tangent there */
} sr_end;

static double sr_turn(qaws_internal_sg_ctx const* c, double const* at, double const* din, double const* dout, int* ok)
{
	double aa = 0, ab = 0;
	*ok = qaws_internal_sg_angle(c, at[0], at[1], din, &aa) == QAWS_STATUS_OK &&
		qaws_internal_sg_angle(c, at[0], at[1], dout, &ab) == QAWS_STATUS_OK;
	return remainder(ab - aa, 2 * 3.14159265358979323846);
}

static qaws_status sr_corner_disks(sr_tube* tb, qaws_clip_result const* res, int grow)
{
	typedef struct sr_piece { qaws_curve const* cv; double lo, hi; } sr_piece;
	sr_piece* seq = NULL;
	sr_end *arr = NULL, *dep = NULL;
	unsigned int ns, cap = 0, na = 0, capa = 0, nd = 0, capd = 0, i, j, k, np = qaws_clip_result_get_path_count(res);
	qaws_status st = QAWS_STATUS_OK;
	qaws_internal_sg_ctx const* c = tb->c;
#if QAWS_SCALAR_IS_FLOAT
	double const turn_min = 1e-3, same = 1e-5;
#else
	double const turn_min = 1e-6, same = 1e-9;
#endif
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
	{
		qaws_path_2d p;
		if ((st = qaws_clip_result_get_path(res, i, &p)) != QAWS_STATUS_OK)
			break;
		/* the loop's pieces in order, a cut as cv NULL */
		ns = 0;
		for (j = 0; j < p.curve_count && st == QAWS_STATUS_OK; j++)
		{
			qaws_curve const* cv = p.curves[j];
			double* br = NULL;
			unsigned int nb = 0;
			if ((st = sr_breaks(cv, &br, &nb)) != QAWS_STATUS_OK)
				break;
			for (k = 0; k + 1 < nb; k++)
			{
				if (!(br[k + 1] > br[k]))
					continue;
				if (!sr_grow((void**)&seq, &cap, ns + 1, sizeof(sr_piece)))
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				seq[ns].cv = sr_on_seam(c, cv, br[k], br[k + 1]) || sr_on_pole(c, cv, br[k], br[k + 1]) ? NULL : cv;
				seq[ns].lo = br[k]; seq[ns].hi = br[k + 1];
				ns++;
			}
			free(br);
		}
		for (j = 0; j < ns && st == QAWS_STATUS_OK; j++)
		{
			sr_piece const* cur = &seq[j];
			sr_piece const* prev = &seq[(j + ns - 1) % ns];
			sr_piece const* next = &seq[(j + 1) % ns];
			qaws_eval_result_2d a, b;
			if (!cur->cv)
				continue;
			qaws_curve_evaluate_2d(&b, cur->cv, (qaws_scalar)cur->lo, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			if (!prev->cv)
			{
				/* leaving a cut: a departure to match */
				if (!sr_grow((void**)&dep, &capd, nd + 1, sizeof(sr_end))) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
				dep[nd].p[0] = b.position.x; dep[nd].p[1] = b.position.y; dep[nd].d[0] = b.d1.x; dep[nd].d[1] = b.d1.y;
				nd++;
			}
			else
			{
				double da[2], db[2], at[2], turn;
				int ok;
				sr_piece_end(prev->cv, prev->lo, prev->hi, &a);
				da[0] = a.d1.x; da[1] = a.d1.y; db[0] = b.d1.x; db[1] = b.d1.y;
				at[0] = b.position.x; at[1] = b.position.y;
				turn = sr_turn(c, at, da, db, &ok);
				if (!ok || (grow && turn > turn_min) || (!grow && turn < -turn_min))
					st = sr_disk(tb, b.position);
			}
			if (st == QAWS_STATUS_OK && !next->cv)
			{
				/* arriving at a cut: an arrival to match */
				sr_piece_end(cur->cv, cur->lo, cur->hi, &a);
				if (!sr_grow((void**)&arr, &capa, na + 1, sizeof(sr_end))) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
				arr[na].p[0] = a.position.x; arr[na].p[1] = a.position.y; arr[na].d[0] = a.d1.x; arr[na].d[1] = a.d1.y;
				na++;
			}
		}
	}
	/* the arrivals at cuts, each with the departure from its point */
	for (i = 0; i < na && st == QAWS_STATUS_OK; i++)
		for (j = 0; j < nd; j++)
		{
			double du = arr[i].p[0] - dep[j].p[0], dv = arr[i].p[1] - dep[j].p[1], turn;
			int ok;
			if (c->pu > 0) du -= c->pu * floor(du / c->pu + 0.5);
			if (c->pv > 0) dv -= c->pv * floor(dv / c->pv + 0.5);
			if (fabs(du) > same * (c->u1 - c->u0) || fabs(dv) > same * (c->v1 - c->v0))
				continue;
			turn = sr_turn(c, dep[j].p, arr[i].d, dep[j].d, &ok);
			if (!ok || (grow && turn > turn_min) || (!grow && turn < -turn_min))
			{
				qaws_vec2 q;
				q.x = (qaws_scalar)dep[j].p[0]; q.y = (qaws_scalar)dep[j].p[1];
				st = sr_disk(tb, q);
			}
			break;
		}
	free(seq); free(arr); free(dep);
	return st;
}

qaws_status qaws_surface_offset(qaws_surface const* surface, qaws_fill_rule fill_rule,
	qaws_surface_path const* paths, unsigned int path_count, qaws_scalar distance, qaws_scalar tolerance,
	qaws_surface_region** out_region)
{
	double d = (double)distance, r = fabs(d);
	qaws_surface_region *R = NULL, *out = NULL;
	qaws_surface_path* closed = NULL;
	unsigned int nc = 0, i, k;
	sr_tube tb;
	sr_pool cut;
	qaws_clip_result* T = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!surface || !out_region || (!paths && path_count) || !(r == r))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_region = NULL;
	memset(&tb, 0, sizeof(tb));
	memset(&cut, 0, sizeof(cut));
	closed = (qaws_surface_path*)malloc(sizeof(qaws_surface_path) * (path_count + 1));
	out = (qaws_surface_region*)calloc(1, sizeof(*out));
	if (!closed || !out)
	{
		free(closed); free(out);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	out->s = surface;
	st = qaws_internal_sg_init(&out->ctx, surface, (double)tolerance);
	tb.s = surface; tb.c = &out->ctx; tb.r = r; tb.tol = (double)tolerance;
	for (i = 0; i < path_count && st == QAWS_STATUS_OK; i++)
	{
		for (k = 0; k < paths[i].curve_count; k++)
			if (!paths[i].curves[k] || qaws_surface_curve_get_surface(paths[i].curves[k]) != surface)
				st = QAWS_STATUS_INVALID_ARGUMENT;
		if (paths[i].closed && paths[i].curve_count)
			closed[nc++] = paths[i];
	}
	/* the region the closed paths bound */
	if (st == QAWS_STATUS_OK && nc)
		st = qaws_surface_clip(surface, QAWS_CLIP_NONE, fill_rule, closed, nc, NULL, 0, &R);
	/* the tube: its boundary (inside on the left of its loops) swept along
	   the conormal geodesics on the side the offset moves to, circles round
	   the corners; open paths on both sides, circles at every vertex */
	if (st == QAWS_STATUS_OK && r > 0)
	{
		unsigned int np = R ? qaws_clip_result_get_path_count(R->res) : 0;
		for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
		{
			qaws_path_2d p;
			if ((st = qaws_clip_result_get_path(R->res, i, &p)) != QAWS_STATUS_OK)
				break;
			for (k = 0; k < p.curve_count && st == QAWS_STATUS_OK; k++)
				st = sr_tube_curve(&tb, p.curves[k], d > 0 ? 2 : 1, 0);
		}
		if (st == QAWS_STATUS_OK && R)
			st = sr_corner_disks(&tb, R->res, d > 0);
		for (i = 0; i < path_count && st == QAWS_STATUS_OK; i++)
			if (!paths[i].closed)
				for (k = 0; k < paths[i].curve_count && st == QAWS_STATUS_OK; k++)
					st = sr_tube_curve(&tb, qaws_surface_curve_get_uv(paths[i].curves[k]), 3, 1);
	}
	sr_fix_views(&tb.raw);
	if (st == QAWS_STATUS_OK && tb.raw.np)
	{
		/* in one period on closed surfaces, then united */
		if (out->ctx.pu > 0 || out->ctx.pv > 0)
		{
			st = sr_period(&out->ctx, tb.raw.paths, tb.raw.np, &cut, NULL);
			sr_fix_views(&cut);
			if (st == QAWS_STATUS_OK)
				st = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, cut.paths, cut.np, NULL, 0, &T);
		}
		else
			st = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, tb.raw.paths, tb.raw.np, NULL, 0, &T);
	}
	if (st == QAWS_STATUS_OK)
	{
		/* the region with the tube (grown) or less it (shrunk) */
		unsigned int na = R ? qaws_clip_result_get_path_count(R->res) : 0, nb = T ? qaws_clip_result_get_path_count(T) : 0;
		qaws_path_2d* qa = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (na + 1));
		qaws_path_2d* qb = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (nb + 1));
		if (!qa || !qb)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < na && st == QAWS_STATUS_OK; i++) st = qaws_clip_result_get_path(R->res, i, &qa[i]);
		for (i = 0; i < nb && st == QAWS_STATUS_OK; i++) st = qaws_clip_result_get_path(T, i, &qb[i]);
		if (st == QAWS_STATUS_OK)
			st = qaws_clip_boolean(d < 0 && na ? QAWS_CLIP_DIFFERENCE : QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, qa, na, qb, nb, &out->res);
		free(qa); free(qb);
	}
	qaws_clip_result_destroy(T);
	sr_pool_free(&tb.raw);
	sr_pool_free(&cut);
	for (i = 0; i < tb.nd; i++) qaws_surface_curve_destroy(tb.disks[i]);
	free(tb.disks);
	free(tb.centres);
	qaws_surface_region_destroy(R);
	free(closed);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_region_destroy(out);
		return st;
	}
	*out_region = out;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Sampling and meshing                                                    */
/* ------------------------------------------------------------------------ */

/* the region's loops in (u, v) (views into its result) */
static qaws_status sr_loops(qaws_surface_region const* r, qaws_path_2d** out, unsigned int* n)
{
	unsigned int np = qaws_clip_result_get_path_count(r->res), i;
	qaws_status st = QAWS_STATUS_OK;
	*out = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (np + 1));
	*n = np;
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
		st = qaws_clip_result_get_path(r->res, i, &(*out)[i]);
	return st;
}

/* the sampler's density: the area element, times the caller's */
typedef struct sr_dens
{
	qaws_internal_sg_ctx const* c;
	qaws_scalar (*f)(void*, qaws_scalar, qaws_scalar);
	void* user;
} sr_dens;

static qaws_scalar sr_area_density(void* user, qaws_scalar u, qaws_scalar v)
{
	sr_dens const* d = (sr_dens const*)user;
	double a = 0;
	if (qaws_internal_sg_area_element(d->c, (double)u, (double)v, &a) != QAWS_STATUS_OK)
		return 0;
	if (d->f)
		a *= (double)d->f(d->user, u, v);
	return (qaws_scalar)(a > 0 ? a : 0);
}

/* sample elimination (Yuksel 2015) in space: of the m candidates the n
   kept, each removal the candidate most crowded by its neighbours within
   2 r_max, r_max the densest packing's radius for n points over the area */
typedef struct sr_cell_key
{
	long long x, y, z;
	unsigned int i;
} sr_cell_key;

static int sr_key_cmp(void const* a, void const* b)
{
	sr_cell_key const* p = (sr_cell_key const*)a;
	sr_cell_key const* q = (sr_cell_key const*)b;
	if (p->x != q->x) return p->x < q->x ? -1 : 1;
	if (p->y != q->y) return p->y < q->y ? -1 : 1;
	if (p->z != q->z) return p->z < q->z ? -1 : 1;
	return 0;
}

static void sr_heap_down(unsigned int* heap, unsigned int* pos, double const* w, unsigned int n, unsigned int i)
{
	for (;;)
	{
		unsigned int l = 2 * i + 1, r = l + 1, m = i, t;
		if (l < n && w[heap[l]] > w[heap[m]]) m = l;
		if (r < n && w[heap[r]] > w[heap[m]]) m = r;
		if (m == i)
			return;
		t = heap[i]; heap[i] = heap[m]; heap[m] = t;
		pos[heap[i]] = i; pos[heap[m]] = m;
		i = m;
	}
}

static qaws_status sr_eliminate(qaws_vec3 const* p, unsigned int m, unsigned int n, double area, unsigned char* keep)
{
	double rmax = sqrt(area / (2 * sqrt(3.0) * (n ? n : 1))), rr = 2 * rmax;
	sr_cell_key* keys = (sr_cell_key*)malloc(sizeof(sr_cell_key) * (m + 1));
	unsigned int *start = NULL, *nbr = NULL, *heap = NULL, *pos = NULL, nn = 0, cap = 0, i, left;
	double *w = NULL, *wd = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!keys)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < m; i++)
	{
		keys[i].x = (long long)floor((double)p[i].x / rr);
		keys[i].y = (long long)floor((double)p[i].y / rr);
		keys[i].z = (long long)floor((double)p[i].z / rr);
		keys[i].i = i;
		keep[i] = 1;
	}
	qsort(keys, m, sizeof(sr_cell_key), sr_key_cmp);
	/* the neighbour lists, as one array with starts */
	start = (unsigned int*)malloc(sizeof(unsigned int) * (m + 1));
	w = (double*)calloc(m + 1, sizeof(double));
	if (!start || !w)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	{
		unsigned int* order = (unsigned int*)malloc(sizeof(unsigned int) * (m + 1));
		if (!order) st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < m && st == QAWS_STATUS_OK; i++) order[keys[i].i] = i;
		for (i = 0; i < m && st == QAWS_STATUS_OK; i++)
		{
			sr_cell_key const* ki = &keys[order[i]];
			int dx, dy, dz;
			start[i] = nn;
			for (dx = -1; dx <= 1; dx++)
				for (dy = -1; dy <= 1; dy++)
					for (dz = -1; dz <= 1; dz++)
					{
						sr_cell_key probe;
						unsigned int lo = 0, hi = m;
						probe.x = ki->x + dx; probe.y = ki->y + dy; probe.z = ki->z + dz;
						/* the first key not below the cell */
						while (lo < hi)
						{
							unsigned int mid = lo + (hi - lo) / 2;
							if (sr_key_cmp(&keys[mid], &probe) < 0) lo = mid + 1; else hi = mid;
						}
						for (; lo < m && sr_key_cmp(&keys[lo], &probe) == 0; lo++)
						{
							unsigned int j = keys[lo].i;
							double d;
							if (j == i)
								continue;
							d = sqrt(((double)p[i].x - p[j].x) * ((double)p[i].x - p[j].x) + ((double)p[i].y - p[j].y) * ((double)p[i].y - p[j].y) +
								((double)p[i].z - p[j].z) * ((double)p[i].z - p[j].z));
							if (d >= rr)
								continue;
							if (!sr_grow((void**)&nbr, &cap, nn + 1, sizeof(unsigned int)))
							{
								st = QAWS_STATUS_ALLOCATION_FAILURE;
								break;
							}
							nbr[nn++] = j;
						}
					}
		}
		start[m] = nn;
		free(order);
	}
	free(keys);
	/* the weights: sum over the neighbours of (1 - d / 2 r_max)^8 */
	if (st == QAWS_STATUS_OK)
	{
		wd = (double*)malloc(sizeof(double) * (nn + 1));
		heap = (unsigned int*)malloc(sizeof(unsigned int) * (m + 1));
		pos = (unsigned int*)malloc(sizeof(unsigned int) * (m + 1));
		if (!wd || !heap || !pos)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < m && st == QAWS_STATUS_OK; i++)
	{
		unsigned int k;
		for (k = start[i]; k < start[i + 1]; k++)
		{
			unsigned int j = nbr[k];
			double d = sqrt(((double)p[i].x - p[j].x) * ((double)p[i].x - p[j].x) + ((double)p[i].y - p[j].y) * ((double)p[i].y - p[j].y) +
				((double)p[i].z - p[j].z) * ((double)p[i].z - p[j].z));
			wd[k] = pow(1 - d / rr, 8);
			w[i] += wd[k];
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; i < m; i++) { heap[i] = i; pos[i] = i; }
		for (i = m / 2 + 1; i-- > 0;) sr_heap_down(heap, pos, w, m, i);
		left = m;
		while (left > n)
		{
			unsigned int top = heap[0], k;
			keep[top] = 0;
			left--;
			heap[0] = heap[left];
			pos[heap[0]] = 0;
			sr_heap_down(heap, pos, w, left, 0);
			/* its neighbours less crowded: their weights fall, they sift down */
			for (k = start[top]; k < start[top + 1]; k++)
			{
				unsigned int j = nbr[k];
				if (!keep[j])
					continue;
				w[j] -= wd[k];
				sr_heap_down(heap, pos, w, left, pos[j]);
			}
		}
	}
	free(start); free(nbr); free(w); free(wd); free(heap); free(pos);
	return st;
}

qaws_status qaws_surface_region_sample(qaws_surface_region const* region, qaws_surface_sample_desc const* desc,
	unsigned int count, qaws_vec2* out_uv, qaws_vec3* out_points)
{
	unsigned int kind = desc ? desc->kind : QAWS_REGION_SAMPLE_RANDOM, np = 0, i;
	unsigned long long seed = desc ? desc->seed : 0;
	qaws_path_2d* loops = NULL;
	qaws_region_sampler_desc sd;
	qaws_region_sampler* sm = NULL;
	qaws_vec2* uv = NULL;
	sr_dens dens;
	qaws_status st;
	if (!region || (count && !out_uv && !out_points))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (count == 0)
		return QAWS_STATUS_OK;
	if ((kind & QAWS_REGION_SAMPLE_BOX) || (kind == QAWS_REGION_SAMPLE_BLUE_NOISE && desc && desc->density))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if ((st = sr_loops(region, &loops, &np)) != QAWS_STATUS_OK || np == 0)
	{
		free(loops);
		return st != QAWS_STATUS_OK ? st : QAWS_STATUS_INVALID_ARGUMENT;
	}
	dens.c = &region->ctx;
	dens.f = desc ? desc->density : NULL;
	dens.user = desc ? desc->density_user : NULL;
	memset(&sd, 0, sizeof(sd));
	sd.paths = loops; sd.path_count = np; sd.fill_rule = QAWS_FILL_NON_ZERO;
	sd.density = sr_area_density; sd.density_user = &dens;
	st = qaws_region_sampler_create(&sd, &sm);
	if (st == QAWS_STATUS_OK)
	{
		if (kind == QAWS_REGION_SAMPLE_BLUE_NOISE)
		{
			/* five candidates a point, uniform in area, thinned in space */
			unsigned int m = 5 * count, k = 0;
			qaws_vec2* cu = (qaws_vec2*)malloc(sizeof(qaws_vec2) * m);
			qaws_vec3* cp = (qaws_vec3*)malloc(sizeof(qaws_vec3) * m);
			unsigned char* keep = (unsigned char*)malloc(m);
			uv = (qaws_vec2*)malloc(sizeof(qaws_vec2) * count);
			if (!cu || !cp || !keep || !uv)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			if (st == QAWS_STATUS_OK)
				st = qaws_region_sampler_generate(sm, m, QAWS_REGION_SAMPLE_RANDOM, seed, cu);
			for (i = 0; i < m && st == QAWS_STATUS_OK; i++)
			{
				qaws_surface_eval_result e;
				st = qaws_internal_sg_eval(&region->ctx, (double)cu[i].x, (double)cu[i].y, QAWS_SURFACE_EVAL_POSITION, &e);
				cp[i] = e.position;
			}
			if (st == QAWS_STATUS_OK)
				st = sr_eliminate(cp, m, count, (double)qaws_region_sampler_get_mass(sm), keep);
			for (i = 0; i < m && st == QAWS_STATUS_OK; i++)
				if (keep[i] && k < count)
					uv[k++] = cu[i];
			free(cu); free(cp); free(keep);
		}
		else
		{
			uv = (qaws_vec2*)malloc(sizeof(qaws_vec2) * count);
			if (!uv)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			else
				st = qaws_region_sampler_generate(sm, count, kind, seed, uv);
		}
	}
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
	{
		if (out_uv)
			out_uv[i] = uv[i];
		if (out_points)
		{
			qaws_surface_eval_result e;
			st = qaws_internal_sg_eval(&region->ctx, (double)uv[i].x, (double)uv[i].y, QAWS_SURFACE_EVAL_POSITION, &e);
			out_points[i] = e.position;
		}
	}
	free(uv);
	qaws_region_sampler_destroy(sm);
	free(loops);
	return st;
}

/* the mesh's size field: the surface length L in (u, v) along the metric's
   stretchiest direction */
typedef struct sr_size
{
	qaws_internal_sg_ctx const* c;
	double L;
} sr_size;

static double sr_stretch(qaws_internal_sg_ctx const* c, double u, double v)
{
	double g[3];
	if (qaws_internal_sg_first_form(c, u, v, g) != QAWS_STATUS_OK)
		return 0;
	return sqrt(0.5 * (g[0] + g[2] + sqrt((g[0] - g[2]) * (g[0] - g[2]) + 4 * g[1] * g[1])));
}

static qaws_scalar sr_size_fn(void* user, qaws_vec2 p)
{
	sr_size const* s = (sr_size const*)user;
	double k = sr_stretch(s->c, (double)p.x, (double)p.y);
	return (qaws_scalar)(k > 0 ? s->L / k : s->L);
}

qaws_status qaws_surface_region_triangulate(qaws_surface_region const* region, qaws_surface_mesh_desc const* desc,
	struct qaws_mesh_2d** out_mesh)
{
	qaws_path_2d* loops = NULL;
	unsigned int np = 0, i, j;
	qaws_triangulate_desc td;
	qaws_vec2 lo, hi;
	sr_size sz;
	double L = desc ? (double)desc->edge_length : 0, tol = desc ? (double)desc->tolerance : 0, kmax = 0;
	qaws_status st;
	if (!region || !out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if ((st = sr_loops(region, &loops, &np)) != QAWS_STATUS_OK || np == 0)
	{
		free(loops);
		return st != QAWS_STATUS_OK ? st : QAWS_STATUS_INVALID_ARGUMENT;
	}
	if (!(L > 0))
	{
		qaws_scalar a = 0;
		qaws_surface_region_get_area(region, &a);
		L = sqrt((double)a) / 30;
	}
	if (!(tol > 0))
		tol = L / 100;
	/* the boundary tolerance in (u, v): at the stretchiest point of the box */
	{
		qaws_vec2 a, b;
		lo.x = lo.y = (qaws_scalar)1e30; hi.x = hi.y = (qaws_scalar)-1e30;
		for (i = 0; i < np; i++)
			if (qaws_path_compute_bounds_2d(&loops[i], &a, &b) == QAWS_STATUS_OK)
			{
				if (a.x < lo.x) lo.x = a.x;
				if (a.y < lo.y) lo.y = a.y;
				if (b.x > hi.x) hi.x = b.x;
				if (b.y > hi.y) hi.y = b.y;
			}
	}
	for (i = 0; i <= 16; i++)
		for (j = 0; j <= 16; j++)
		{
			double k = sr_stretch(&region->ctx, (double)lo.x + ((double)hi.x - lo.x) * i / 16, (double)lo.y + ((double)hi.y - lo.y) * j / 16);
			if (k > kmax) kmax = k;
		}
	sz.c = &region->ctx;
	sz.L = L;
	memset(&td, 0, sizeof(td));
	td.paths = loops; td.path_count = np; td.fill_rule = QAWS_FILL_NON_ZERO;
	td.tolerance = (qaws_scalar)(kmax > 0 ? tol / kmax : tol);
	td.min_angle = (qaws_scalar)(desc && desc->min_angle > 0 ? desc->min_angle : 25);
	td.size_fn = sr_size_fn; td.size_user = &sz;
	st = qaws_triangulate_region_2d(&td, out_mesh);
	if (st == QAWS_STATUS_OK && desc && desc->remesh_iterations)
		st = qaws_mesh_2d_remesh(*out_mesh, (qaws_scalar)L, sr_size_fn, &sz, desc->remesh_iterations);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_2d_destroy(*out_mesh);
		*out_mesh = NULL;
	}
	free(loops);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  Medial axis                                                             */
/* ------------------------------------------------------------------------ */

#if QAWS_SCALAR_IS_FLOAT
#define SM_FD_STEP 1e-4     /* the junction's difference step (of the diagonal), past the rounding */
#define SM_CORNER 1e-3     /* a turn this sharp (radians) is a corner */
#else
#define SM_FD_STEP 1e-7
#define SM_CORNER 1e-6
#endif

/* a boundary feature: a smooth piece of a loop's curve, or a corner (a
   point; it is the nearest boundary only inside a reflex corner) */
typedef struct sm_feature
{
	int point;
	qaws_curve const* cv;
	double lo, hi;
	double p[2];
} sm_feature;

/* a traced point: on the geodesic from feature f leaving at s (a piece's
   parameter, a corner's angle), at distance t */
typedef struct sm_trace
{
	double x[3];
	unsigned int f;
	double s, t;
	long long cell[3];
} sm_trace;

/* the nearest boundary of a point */
typedef struct sm_label
{
	unsigned int f;
	double s, t;
	double foot[3];
} sm_label;

typedef struct sm_seg
{
	double uv[4];
	double r[2];
} sm_seg;

struct qaws_surface_medial_axis
{
	sm_seg* seg;
	unsigned int nseg, capseg;
	double* junction;            /* u, v, radius */
	unsigned int njun, capjun;
	double best[3];              /* the largest inscribed circle: u, v, radius */
};

typedef struct sm_ctx
{
	qaws_internal_sg_ctx const* c;
	sm_feature* f;
	unsigned int nf, capf;
	sm_trace* tr;
	unsigned int ntr, captr;
	double cell;
} sm_ctx;

static int sm_trace_cmp(void const* a, void const* b)
{
	sm_trace const* p = (sm_trace const*)a;
	sm_trace const* q = (sm_trace const*)b;
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (p->cell[k] != q->cell[k])
			return p->cell[k] < q->cell[k] ? -1 : 1;
	return 0;
}

static qaws_status sm_point3(qaws_internal_sg_ctx const* c, double u, double v, double* x)
{
	qaws_surface_eval_result e;
	qaws_status st = qaws_internal_sg_eval(c, u, v, QAWS_SURFACE_EVAL_POSITION, &e);
	x[0] = e.position.x; x[1] = e.position.y; x[2] = e.position.z;
	return st;
}

static double sm_d3(double const* a, double const* b)
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

static qaws_status sm_add_trace(sm_ctx* m, double const* uv, unsigned int f, double s, double t)
{
	sm_trace* q;
	qaws_status st;
	if (!sr_grow((void**)&m->tr, &m->captr, m->ntr + 1, sizeof(sm_trace)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	q = &m->tr[m->ntr];
	if ((st = sm_point3(m->c, uv[0], uv[1], q->x)) != QAWS_STATUS_OK)
		return st;
	q->f = f; q->s = s; q->t = t;
	m->ntr++;
	return QAWS_STATUS_OK;
}

/* the geodesic from state y traced in steps h up to D, every point kept;
   it stops quietly where it leaves the domain */
static qaws_status sm_trace_ray(sm_ctx* m, double* y, unsigned int f, double s, double h, double D)
{
	double t = 0;
	qaws_status st = sm_add_trace(m, y, f, s, 0);
	while (st == QAWS_STATUS_OK && t < D)
	{
		if (qaws_internal_sg_integrate(m->c, y, h) != QAWS_STATUS_OK)
			return QAWS_STATUS_OK;
		t += h;
		if (!(y[1] >= m->c->v0 - 1e-9 && y[1] <= m->c->v1 + 1e-9) || (m->c->pu <= 0 && !(y[0] >= m->c->u0 - 1e-9 && y[0] <= m->c->u1 + 1e-9)))
			return QAWS_STATUS_OK;
		st = sm_add_trace(m, y, f, s, t);
	}
	return st;
}

/* the traced point minimising t + |x - p| among those near x */
static int sm_nearest(sm_ctx const* m, double const* x, sm_label* out)
{
	long long c0[3];
	int ring, found = 0, k;
	double best = 1e300;
	for (k = 0; k < 3; k++) c0[k] = (long long)floor(x[k] / m->cell);
	for (ring = 1; ring <= 8 && !found; ring *= 2)
	{
		long long dx, dy, dz;
		for (dx = -ring; dx <= ring; dx++)
			for (dy = -ring; dy <= ring; dy++)
				for (dz = -ring; dz <= ring; dz++)
				{
					sm_trace key;
					unsigned int lo = 0, hi = m->ntr;
					key.cell[0] = c0[0] + dx; key.cell[1] = c0[1] + dy; key.cell[2] = c0[2] + dz;
					while (lo < hi)
					{
						unsigned int mid = lo + (hi - lo) / 2;
						if (sm_trace_cmp(&m->tr[mid], &key) < 0) lo = mid + 1; else hi = mid;
					}
					for (; lo < m->ntr && sm_trace_cmp(&m->tr[lo], &key) == 0; lo++)
					{
						double d = m->tr[lo].t + sm_d3(x, m->tr[lo].x);
						if (d < best)
						{
							best = d;
							out->f = m->tr[lo].f;
							out->s = m->tr[lo].s;
							out->t = d;
							found = 1;
						}
					}
				}
	}
	return found;
}

/* the distance from x (u, v) to feature f, refined from (s, t): a piece by
   its foot (Newton on the conormal geodesic), a corner by the geodesic to
   it; 1e300 when it does not land */
static double sm_dist(sm_ctx const* m, unsigned int f, double const* x, double* s, double* t)
{
	sm_feature const* ft = &m->f[f];
	double s0 = *s, t0 = *t, best = 1e300, bs = 0, bt = 0, x3[3], q3[3];
	unsigned int k;
	int ok;
	if (!(t0 > 0))
		t0 = 1e-3 * m->c->diag;
	*t = t0;
	ok = ft->point ? qaws_internal_sg_reach(m->c, ft->p, x, s, t) : qaws_internal_sg_foot(m->c, ft->cv, 1, ft->lo, ft->hi, x, s, t);
	if (ok)
		return *t;
	/* from a fan of starts: along the piece, or round the corner; the
	   nearest that lands (the chord in space for the length) */
	if (sm_point3(m->c, x[0], x[1], x3) != QAWS_STATUS_OK)
		return 1e300;
	for (k = 0; k < 8; k++)
	{
		double ss, tt;
		if (ft->point)
		{
			ss = 2 * 3.14159265358979323846 * k / 8;
			sm_point3(m->c, ft->p[0], ft->p[1], q3);
			tt = sm_d3(x3, q3);
			ok = qaws_internal_sg_reach(m->c, ft->p, x, &ss, &tt);
		}
		else
		{
			qaws_eval_result_2d r;
			ss = ft->lo + (ft->hi - ft->lo) * (k + 0.5) / 8;
			qaws_curve_evaluate_2d(&r, ft->cv, (qaws_scalar)ss, QAWS_EVAL_FLAG_POSITION);
			sm_point3(m->c, r.position.x, r.position.y, q3);
			tt = sm_d3(x3, q3);
			ok = tt > 0 && qaws_internal_sg_foot(m->c, ft->cv, 1, ft->lo, ft->hi, x, &ss, &tt);
		}
		if (ok && tt < best)
		{
			best = tt; bs = ss; bt = tt;
		}
	}
	(void)s0;
	if (best < 1e300)
	{
		*s = bs; *t = bt;
	}
	return best;
}

static qaws_status sm_add_seg(qaws_surface_medial_axis* a, double const* p, double rp, double const* q, double rq)
{
	sm_seg* g;
	if (!sr_grow((void**)&a->seg, &a->capseg, a->nseg + 1, sizeof(sm_seg)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	g = &a->seg[a->nseg++];
	g->uv[0] = p[0]; g->uv[1] = p[1]; g->uv[2] = q[0]; g->uv[3] = q[1];
	g->r[0] = rp; g->r[1] = rq;
	return QAWS_STATUS_OK;
}

/* a crossing on a mesh edge: where the distances to the two ends' nearest
   features agree */
typedef struct sm_cross
{
	unsigned int a, b;              /* the edge's vertices (a < b), or none */
	int ok;
	double uv[2], r;
	unsigned int fa, fb;
	double sa, ta, sb, tb;          /* the features' foot data there */
	int dom;                        /* a third feature nearer: off the axis proper */
	unsigned int fc;
	double sc, tc;
} sm_cross;

qaws_status qaws_surface_region_medial_axis(qaws_surface_region const* region, qaws_surface_medial_desc const* desc,
	qaws_surface_medial_axis** out_axis)
{
	sm_ctx m;
	qaws_surface_medial_axis* ax = NULL;
	qaws_path_2d* loops = NULL;
	unsigned int np = 0, i, j, k;
	qaws_mesh_2d* mesh = NULL;
	qaws_surface_mesh_desc md;
	qaws_scalar area = 0;
	double L, h, D;
	sm_label* lab = NULL;
	double* vx3 = NULL;
	sm_cross* cr = NULL;
	unsigned int ncr = 0, capcr = 0, *hash = NULL, hsize = 0;
	qaws_status st;
	if (!region || !out_axis)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_axis = NULL;
	memset(&m, 0, sizeof(m));
	m.c = &region->ctx;
	qaws_surface_region_get_area(region, &area);
	L = desc && desc->resolution > 0 ? (double)desc->resolution : sqrt((double)area) / 40;
	h = 0.5 * L;
	D = 2 * sqrt((double)area) + L;
	m.cell = 2 * h;
	ax = (qaws_surface_medial_axis*)calloc(1, sizeof(*ax));
	if (!ax)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* the features: each loop's smooth pieces (cuts along seams and poles
	   left out), and its corners: where one piece ends and the next starts
	   turning (not across a cut, where the boundary goes on smoothly) */
	st = sr_loops(region, &loops, &np);
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
	{
		/* the loop's pieces in order, a cut marked by cv NULL */
		sm_feature* seq = NULL;
		unsigned int ns = 0, caps = 0, first = m.nf;
		for (j = 0; j < loops[i].curve_count && st == QAWS_STATUS_OK; j++)
		{
			qaws_curve const* cv = loops[i].curves[j];
			double* br = NULL;
			unsigned int nb = 0;
			if ((st = sr_breaks(cv, &br, &nb)) != QAWS_STATUS_OK)
				break;
			for (k = 0; k + 1 < nb && st == QAWS_STATUS_OK; k++)
			{
				int cut;
				if (!(br[k + 1] > br[k]))
					continue;
				cut = sr_on_seam(m.c, cv, br[k], br[k + 1]) || sr_on_pole(m.c, cv, br[k], br[k + 1]);
				if (!sr_grow((void**)&seq, &caps, ns + 1, sizeof(sm_feature)))
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				seq[ns].point = 0; seq[ns].cv = cut ? NULL : cv; seq[ns].lo = br[k]; seq[ns].hi = br[k + 1];
				ns++;
			}
			free(br);
		}
		for (j = 0; j < ns && st == QAWS_STATUS_OK; j++)
		{
			sm_feature const* cur = &seq[j];
			sm_feature const* prev = &seq[(j + ns - 1) % ns];
			if (!cur->cv)
				continue;
			if (!sr_grow((void**)&m.f, &m.capf, m.nf + 2, sizeof(sm_feature)))
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			m.f[m.nf++] = *cur;
			/* a corner at its start: the previous piece's end tangent and its
			   start tangent apart (as angles in the metric) */
			if (prev->cv)
			{
				qaws_eval_result_2d a, b;
				double da[2], db[2], aa = 0, ab = 0, turn;
				sr_piece_end(prev->cv, prev->lo, prev->hi, &a);
				qaws_curve_evaluate_2d(&b, cur->cv, (qaws_scalar)cur->lo, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
				da[0] = a.d1.x; da[1] = a.d1.y; db[0] = b.d1.x; db[1] = b.d1.y;
				qaws_internal_sg_angle(m.c, (double)b.position.x, (double)b.position.y, da, &aa);
				qaws_internal_sg_angle(m.c, (double)b.position.x, (double)b.position.y, db, &ab);
				turn = fabs(remainder(ab - aa, 2 * 3.14159265358979323846));
				if (turn > SM_CORNER)
				{
					m.f[m.nf].point = 1; m.f[m.nf].cv = NULL;
					m.f[m.nf].p[0] = b.position.x; m.f[m.nf].p[1] = b.position.y;
					m.nf++;
				}
			}
		}
		(void)first;
		free(seq);
	}
	/* the traces: conormal geodesics inward from every piece, fans from every
	   end */
	for (i = 0; i < m.nf && st == QAWS_STATUS_OK; i++)
	{
		sm_feature const* ft = &m.f[i];
		if (!ft->point)
		{
			/* samples about h apart on the surface */
			double len = 0, prev[3], cur[3];
			unsigned int n, q;
			for (q = 0; q <= 32; q++)
			{
				qaws_eval_result_2d r;
				qaws_curve_evaluate_2d(&r, ft->cv, (qaws_scalar)(ft->lo + (ft->hi - ft->lo) * q / 32), QAWS_EVAL_FLAG_POSITION);
				sm_point3(m.c, r.position.x, r.position.y, cur);
				if (q) len += sm_d3(prev, cur);
				prev[0] = cur[0]; prev[1] = cur[1]; prev[2] = cur[2];
			}
			n = (unsigned int)ceil(len / h);
			if (n < 2) n = 2;
			for (q = 0; q <= n && st == QAWS_STATUS_OK; q++)
			{
				double y[4], s = ft->lo + (ft->hi - ft->lo) * q / n;
				if (qaws_internal_sg_conormal(m.c, ft->cv, s, 1, y) != QAWS_STATUS_OK)
					continue;
				st = sm_trace_ray(&m, y, i, s, h, D);
			}
		}
		else
		{
			unsigned int na = (unsigned int)ceil(2 * 3.14159265358979323846 * D / h), q;
			if (na > 4096) na = 4096;
			for (q = 0; q < na && st == QAWS_STATUS_OK; q++)
			{
				double th = 2 * 3.14159265358979323846 * q / na, end[2], vel[2], y[4];
				if (qaws_internal_sg_shoot(m.c, ft->p, th, 0, end, vel) != QAWS_STATUS_OK)
					continue;
				y[0] = end[0]; y[1] = end[1]; y[2] = vel[0]; y[3] = vel[1];
				st = sm_trace_ray(&m, y, i, th, h, D);
			}
		}
	}
	for (i = 0; i < m.ntr; i++)
		for (k = 0; k < 3; k++)
			m.tr[i].cell[k] = (long long)floor(m.tr[i].x[k] / m.cell);
	if (st == QAWS_STATUS_OK)
		qsort(m.tr, m.ntr, sizeof(sm_trace), sm_trace_cmp);
	/* the mesh, its vertices labelled by their nearest feature */
	if (st == QAWS_STATUS_OK)
	{
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)L;
		st = qaws_surface_region_triangulate(region, &md, &mesh);
	}
	if (st == QAWS_STATUS_OK)
	{
		unsigned int nv = qaws_mesh_2d_get_vertex_count(mesh), nt = qaws_mesh_2d_get_triangle_count(mesh);
		lab = (sm_label*)malloc(sizeof(sm_label) * (nv + 1));
		vx3 = (double*)malloc(sizeof(double) * 3 * (nv + 1));
		hsize = 1;
		while (hsize < 4 * nt + 16) hsize <<= 1;
		hash = (unsigned int*)malloc(sizeof(unsigned int) * hsize);
		if (!lab || !vx3 || !hash)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < hsize && st == QAWS_STATUS_OK; i++) hash[i] = 0xffffffffu;
		for (i = 0; i < nv && st == QAWS_STATUS_OK; i++)
		{
			qaws_vec2 w = qaws_mesh_2d_get_vertex(mesh, i);
			sm_feature const* ft;
			st = sm_point3(m.c, w.x, w.y, vx3 + 3 * i);
			if (!sm_nearest(&m, vx3 + 3 * i, &lab[i]))
			{
				lab[i].f = 0xffffffffu;
				continue;
			}
			ft = &m.f[lab[i].f];
			if (ft->point)
				sm_point3(m.c, ft->p[0], ft->p[1], lab[i].foot);
			else
			{
				qaws_eval_result_2d r;
				qaws_curve_evaluate_2d(&r, ft->cv, (qaws_scalar)lab[i].s, QAWS_EVAL_FLAG_POSITION);
				sm_point3(m.c, r.position.x, r.position.y, lab[i].foot);
			}
		}
		/* the labels checked across the edges whose ends' nearest boundary
		   points are far apart: an end truly nearer the other end's feature
		   (the traced field is first order) takes it, so the crossings
		   bracket */
		for (i = 0; i < nt && st == QAWS_STATUS_OK; i++)
		{
			unsigned int v[3], e;
			qaws_mesh_2d_get_triangle(mesh, i, v);
			for (e = 0; e < 3; e++)
			{
				unsigned int a = v[e], b = v[(e + 1) % 3], z;
				if (lab[a].f == 0xffffffffu || lab[b].f == 0xffffffffu || lab[a].f == lab[b].f ||
					sm_d3(lab[a].foot, lab[b].foot) <= 3 * sm_d3(vx3 + 3 * a, vx3 + 3 * b))
					continue;
				for (z = 0; z < 2; z++)
				{
					/* end p labelled with its feature; the other's feature tried there */
					unsigned int p = z ? b : a, o = z ? a : b;
					qaws_vec2 w = qaws_mesh_2d_get_vertex(mesh, p);
					double x[2], s0 = lab[p].s, t0 = lab[p].t, s1 = lab[o].s, t1 = lab[o].t, d0, d1;
					x[0] = w.x; x[1] = w.y;
					d0 = sm_dist(&m, lab[p].f, x, &s0, &t0);
					d1 = sm_dist(&m, lab[o].f, x, &s1, &t1);
					if (d1 < d0 - m.c->tol && d1 < 1e299)
					{
						sm_feature const* ft = &m.f[lab[o].f];
						lab[p].f = lab[o].f; lab[p].s = s1; lab[p].t = d1;
						if (ft->point)
							sm_point3(m.c, ft->p[0], ft->p[1], lab[p].foot);
						else
						{
							qaws_eval_result_2d r;
							qaws_curve_evaluate_2d(&r, ft->cv, (qaws_scalar)s1, QAWS_EVAL_FLAG_POSITION);
							sm_point3(m.c, r.position.x, r.position.y, lab[p].foot);
						}
					}
					else if (d0 < 1e299)
					{
						lab[p].s = s0; lab[p].t = d0;
					}
				}
			}
		}
		/* the triangles: their edges crossing the axis (the ends' nearest
		   boundary points far apart), each crossing found once */
		for (i = 0; i < nt && st == QAWS_STATUS_OK; i++)
		{
			unsigned int v[3], e, ids[3], nc = 0;
			double cuv[3][2], crad[3];
			qaws_mesh_2d_get_triangle(mesh, i, v);
			for (e = 0; e < 3; e++)
			{
				unsigned int a = v[e] < v[(e + 1) % 3] ? v[e] : v[(e + 1) % 3], b = v[e] < v[(e + 1) % 3] ? v[(e + 1) % 3] : v[e];
				unsigned int slot = (unsigned int)((a * 2654435761u) ^ (b * 40503u)) & (hsize - 1), id = 0xffffffffu;
				while (hash[slot] != 0xffffffffu)
				{
					if (cr[hash[slot]].a == a && cr[hash[slot]].b == b) { id = hash[slot]; break; }
					slot = (slot + 1) & (hsize - 1);
				}
				if (id == 0xffffffffu)
				{
					/* a new edge: find its crossing */
					sm_cross* x;
					if (!sr_grow((void**)&cr, &capcr, ncr + 1, sizeof(sm_cross)))
					{
						st = QAWS_STATUS_ALLOCATION_FAILURE;
						break;
					}
					x = &cr[ncr];
					memset(x, 0, sizeof(*x));
					x->a = a; x->b = b;
					if (lab[a].f != 0xffffffffu && lab[b].f != 0xffffffffu &&
						sm_d3(lab[a].foot, lab[b].foot) > 3 * sm_d3(vx3 + 3 * a, vx3 + 3 * b))
					{
						qaws_vec2 pa = qaws_mesh_2d_get_vertex(mesh, a), pb = qaws_mesh_2d_get_vertex(mesh, b);
						double xa[2], xb[2], sA = lab[a].s, tA = lab[a].t, sB = lab[b].s, tB = lab[b].t, l0 = 0, l1 = 1, g0, g1;
						unsigned int fa = lab[a].f, fb = lab[b].f, it;
						double sA1, tA1, sB0, tB0;
						xa[0] = pa.x; xa[1] = pa.y; xb[0] = pb.x; xb[1] = pb.y;
						/* g = d_A - d_B at both ends, each feature's foot carried along */
						g0 = sm_dist(&m, fa, xa, &sA, &tA);
						sB0 = sB; tB0 = tB;
						g0 -= sm_dist(&m, fb, xa, &sB0, &tB0);
						sA1 = sA; tA1 = tA;
						g1 = sm_dist(&m, fa, xb, &sA1, &tA1);
						g1 -= sm_dist(&m, fb, xb, &sB, &tB);
						if (g0 <= 0 && g1 >= 0 && g1 - g0 < 1e299)
						{
							/* the Illinois method on the edge */
							double gl = g0, gr = g1, lm = 0, sAm = sA, tAm = tA, sBm = sB0, tBm = tB0;
							int side = 0;
							for (it = 0; it < 60; it++)
							{
								double xm[2], gm, da, db;
								lm = gr != gl ? (l0 * gr - l1 * gl) / (gr - gl) : 0.5 * (l0 + l1);
								if (!(lm > l0 && lm < l1)) lm = 0.5 * (l0 + l1);
								xm[0] = xa[0] + lm * (xb[0] - xa[0]); xm[1] = xa[1] + lm * (xb[1] - xa[1]);
								da = sm_dist(&m, fa, xm, &sAm, &tAm);
								db = sm_dist(&m, fb, xm, &sBm, &tBm);
								if (!(da < 1e299 && db < 1e299))
									break;
								gm = da - db;
								x->uv[0] = xm[0]; x->uv[1] = xm[1]; x->r = 0.5 * (da + db);
								/* kept only once found: the distances agree, or the bracket is down to the tolerance */
								if (fabs(gm) <= 1e-3 * m.c->tol || (l1 - l0) * hypot(xb[0] - xa[0], xb[1] - xa[1]) <= 0.05 * m.c->tol)
								{
									x->ok = 1;
									break;
								}
								if (gm < 0) { l0 = lm; gl = gm; if (side == -1) gr *= 0.5; side = -1; }
								else { l1 = lm; gr = gm; if (side == 1) gl *= 0.5; side = 1; }
							}
							x->fa = fa; x->fb = fb;
							x->sa = sAm; x->ta = tAm; x->sb = sBm; x->tb = tBm;
							if (x->ok)
							{
								/* on the axis proper only when no third feature is nearer:
								   then a candidate for the largest inscribed circle */
								double p3[3];
								sm_label l3;
								int dom = 0;
								if (sm_point3(m.c, x->uv[0], x->uv[1], p3) == QAWS_STATUS_OK && sm_nearest(&m, p3, &l3) && l3.f != fa && l3.f != fb)
								{
									double s3 = l3.s, t3 = l3.t, d3 = sm_dist(&m, l3.f, x->uv, &s3, &t3);
									if (d3 < x->r - 10 * m.c->tol)
									{
										dom = 1;
										x->fc = l3.f; x->sc = s3; x->tc = t3;
									}
								}
								x->dom = dom;
								if (!dom && x->r > ax->best[2]) { ax->best[0] = x->uv[0]; ax->best[1] = x->uv[1]; ax->best[2] = x->r; }
							}
						}
					}
					hash[slot] = ncr;
					id = ncr++;
				}
				ids[e] = id;
			}
			if (st != QAWS_STATUS_OK)
				break;
			/* the triangle's crossings on the axis proper, and the features met
			   (a dominated crossing's third feature among them) */
			{
				unsigned int fs[6], nfs = 0, q, w;
				double fsd[6][2];
				for (e = 0; e < 3; e++)
				{
					sm_cross const* x = &cr[ids[e]];
					unsigned int ff[3], nff = 0, z;
					double fd[3][2];
					if (!x->ok)
						continue;
					if (!x->dom)
					{
						cuv[nc][0] = x->uv[0]; cuv[nc][1] = x->uv[1]; crad[nc] = x->r;
						nc++;
					}
					ff[nff] = x->fa; fd[nff][0] = x->sa; fd[nff][1] = x->ta; nff++;
					ff[nff] = x->fb; fd[nff][0] = x->sb; fd[nff][1] = x->tb; nff++;
					if (x->dom) { ff[nff] = x->fc; fd[nff][0] = x->sc; fd[nff][1] = x->tc; nff++; }
					for (q = 0; q < nff; q++)
					{
						int seen = 0;
						for (z = 0; z < nfs; z++) if (fs[z] == ff[q]) seen = 1;
						if (!seen && nfs < 6) { fs[nfs] = ff[q]; fsd[nfs][0] = fd[q][0]; fsd[nfs][1] = fd[q][1]; nfs++; }
					}
				}
				/* the vertices' own features when their nearest boundary points are
				   all far apart: a junction inside, though its crossings did not
				   settle */
				{
					double emax = 0;
					int far = 1;
					for (q = 0; q < 3; q++)
					{
						double ed = sm_d3(vx3 + 3 * v[q], vx3 + 3 * v[(q + 1) % 3]);
						if (ed > emax) emax = ed;
					}
					for (q = 0; q < 3; q++)
						if (lab[v[q]].f == 0xffffffffu || sm_d3(lab[v[q]].foot, lab[v[(q + 1) % 3]].foot) <= 3 * emax)
							far = 0;
					if (far)
						for (q = 0; q < 3; q++)
						{
							unsigned int z;
							int seen = 0;
							for (z = 0; z < nfs; z++) if (fs[z] == lab[v[q]].f) seen = 1;
							if (!seen && nfs < 6) { fs[nfs] = lab[v[q]].f; fsd[nfs][0] = lab[v[q]].s; fsd[nfs][1] = lab[v[q]].t; nfs++; }
						}
				}
				if (nfs >= 3)
				{
					/* a junction: three features equally far, by Newton in (u, v)
					   from the triangle's centre; accepted by its residual */
					double p[2], r = 0, bestres = 1e300, bp[2] = { 0, 0 }, br = 0;
					unsigned int itn, a2;
					qaws_vec2 w0 = qaws_mesh_2d_get_vertex(mesh, v[0]), w1 = qaws_mesh_2d_get_vertex(mesh, v[1]), w2 = qaws_mesh_2d_get_vertex(mesh, v[2]);
					p[0] = ((double)w0.x + w1.x + w2.x) / 3;
					p[1] = ((double)w0.y + w1.y + w2.y) / 3;
					for (itn = 0; itn < 30; itn++)
					{
						double d[3], F[2], J[4], det, res;
						for (q = 0; q < 3; q++) d[q] = sm_dist(&m, fs[q], p, &fsd[q][0], &fsd[q][1]);
						if (!(d[0] < 1e299 && d[1] < 1e299 && d[2] < 1e299))
							break;
						F[0] = d[0] - d[1]; F[1] = d[0] - d[2];
						res = fabs(F[0]) + fabs(F[1]);
						if (res < bestres) { bestres = res; bp[0] = p[0]; bp[1] = p[1]; br = (d[0] + d[1] + d[2]) / 3; }
						if (res <= 1e-3 * m.c->tol)
							break;
						for (a2 = 0; a2 < 2; a2++)
						{
							double px[2], dd[3], sv[3][2], hs = SM_FD_STEP * m.c->diag;
							for (q = 0; q < 3; q++) { sv[q][0] = fsd[q][0]; sv[q][1] = fsd[q][1]; }
							px[0] = p[0] + (a2 == 0 ? hs : 0); px[1] = p[1] + (a2 == 1 ? hs : 0);
							for (q = 0; q < 3; q++) dd[q] = sm_dist(&m, fs[q], px, &sv[q][0], &sv[q][1]);
							J[a2] = ((dd[0] - dd[1]) - F[0]) / hs;
							J[2 + a2] = ((dd[0] - dd[2]) - F[1]) / hs;
						}
						det = J[0] * J[3] - J[1] * J[2];
						if (!(fabs(det) > 0))
							break;
						p[0] -= (J[3] * F[0] - J[1] * F[1]) / det;
						p[1] -= (J[0] * F[1] - J[2] * F[0]) / det;
					}
					r = br;
					/* equally far within the tolerance, and inside (a radius of 0 is the
					   boundary at a corner, not a junction) */
					if (bestres <= 10 * m.c->tol && br > 100 * m.c->tol)
					{
						/* found already from a neighbouring triangle: the same junction */
						unsigned int z, dup = 0;
						for (z = 0; z < ax->njun && !dup; z++)
							if (hypot(ax->junction[3 * z] - bp[0], ax->junction[3 * z + 1] - bp[1]) <= 100 * m.c->tol)
							{
								bp[0] = ax->junction[3 * z]; bp[1] = ax->junction[3 * z + 1]; r = ax->junction[3 * z + 2];
								dup = 1;
							}
						if (dup)
							;
						else if (!sr_grow((void**)&ax->junction, &ax->capjun, 3 * (ax->njun + 1), sizeof(double)))
							st = QAWS_STATUS_ALLOCATION_FAILURE;
						else
						{
							ax->junction[3 * ax->njun] = bp[0]; ax->junction[3 * ax->njun + 1] = bp[1]; ax->junction[3 * ax->njun + 2] = r;
							ax->njun++;
							if (r > ax->best[2]) { ax->best[0] = bp[0]; ax->best[1] = bp[1]; ax->best[2] = r; }
						}
						for (q = 0; q < nc && st == QAWS_STATUS_OK; q++)
							st = sm_add_seg(ax, bp, r, cuv[q], crad[q]);
						continue;
					}
				}
				if (nc == 2)
					st = sm_add_seg(ax, cuv[0], crad[0], cuv[1], crad[1]);
				(void)w;
			}
		}
	}
	free(lab); free(vx3); free(hash); free(cr);
	qaws_mesh_2d_destroy(mesh);
	free(m.f); free(m.tr);
	free(loops);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_medial_axis_destroy(ax);
		return st;
	}
	*out_axis = ax;
	return QAWS_STATUS_OK;
}

void qaws_surface_medial_axis_destroy(qaws_surface_medial_axis* axis)
{
	if (!axis)
		return;
	free(axis->seg);
	free(axis->junction);
	free(axis);
}

unsigned int qaws_surface_medial_axis_get_segment_count(qaws_surface_medial_axis const* axis)
{
	return axis ? axis->nseg : 0;
}

qaws_status qaws_surface_medial_axis_get_segment(qaws_surface_medial_axis const* axis, unsigned int index,
	qaws_vec2 out_uv[2], qaws_scalar out_radius[2])
{
	if (!axis || index >= axis->nseg || !out_uv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out_uv[0].x = (qaws_scalar)axis->seg[index].uv[0]; out_uv[0].y = (qaws_scalar)axis->seg[index].uv[1];
	out_uv[1].x = (qaws_scalar)axis->seg[index].uv[2]; out_uv[1].y = (qaws_scalar)axis->seg[index].uv[3];
	if (out_radius)
	{
		out_radius[0] = (qaws_scalar)axis->seg[index].r[0];
		out_radius[1] = (qaws_scalar)axis->seg[index].r[1];
	}
	return QAWS_STATUS_OK;
}

unsigned int qaws_surface_medial_axis_get_junction_count(qaws_surface_medial_axis const* axis)
{
	return axis ? axis->njun : 0;
}

qaws_status qaws_surface_medial_axis_get_junction(qaws_surface_medial_axis const* axis, unsigned int index,
	qaws_vec2* out_uv, qaws_scalar* out_radius)
{
	if (!axis || index >= axis->njun || !out_uv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out_uv->x = (qaws_scalar)axis->junction[3 * index]; out_uv->y = (qaws_scalar)axis->junction[3 * index + 1];
	if (out_radius) *out_radius = (qaws_scalar)axis->junction[3 * index + 2];
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_medial_axis_get_max_inscribed(qaws_surface_medial_axis const* axis, qaws_vec2* out_centre,
	qaws_scalar* out_radius)
{
	if (!axis || !out_centre || !out_radius)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out_centre->x = (qaws_scalar)axis->best[0]; out_centre->y = (qaws_scalar)axis->best[1];
	*out_radius = (qaws_scalar)axis->best[2];
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Mitered skeletons                                                       */
/* ------------------------------------------------------------------------ */

/* an edge: the boundary between two corners, carried on by the geodesics
   leaving its ends along its tangents (the line it moves as) */
typedef struct mk_edge
{
	qaws_curve* cv;
	double lo, hi;
	unsigned int prev, next;
	int convex[2];                  /* the corners at its start and end */
	double corner[5];               /* its start corner: (u, v), point */
} mk_edge;

/* the signed distance from x to edge e's line: a foot on its left (the
   region's side) from (s, t) or a fan of starts; -1e300 when none lands
   (x right of the line) */
static double mk_dist(qaws_internal_sg_ctx const* c, mk_edge const* e, double const* x, double* s, double* t)
{
	double best = -1e300, x3[3], q3[3], side;
	unsigned int k;
	qaws_scalar tp = 0;
	qaws_vec2 xp;
	qaws_eval_result_2d r;
	if (*t > 0 && qaws_internal_sg_foot(c, e->cv, 1, e->lo, e->hi, x, s, t))
		return *t;
	if (sm_point3(c, x[0], x[1], x3) != QAWS_STATUS_OK)
		return best;
	/* the side first, in (u, v) (the map keeps it): right of the line, none;
	   on it, 0; left, the foot from the nearest point */
	xp.x = (qaws_scalar)x[0]; xp.y = (qaws_scalar)x[1];
	if (qaws_curve_find_closest_parameter_2d(&tp, e->cv, xp) == QAWS_STATUS_OK &&
		qaws_curve_evaluate_2d(&r, e->cv, tp, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1) == QAWS_STATUS_OK && hypot((double)r.d1.x, (double)r.d1.y) > 0)
	{
		side = ((double)r.d1.x * (x[1] - r.position.y) - (double)r.d1.y * (x[0] - r.position.x)) / hypot((double)r.d1.x, (double)r.d1.y);
		if (side < -10 * c->tol)
			return best;
		*s = tp;
		if (side <= 10 * c->tol)
		{
			*t = 0;
			return 0;
		}
		sm_point3(c, r.position.x, r.position.y, q3);
		*t = sm_d3(x3, q3);
		if (qaws_internal_sg_foot(c, e->cv, 1, e->lo, e->hi, x, s, t))
			return *t;
	}
	for (k = 0; k < 8; k++)
	{
		double ss = e->lo + (e->hi - e->lo) * (k + 0.5) / 8, tt;
		qaws_curve_evaluate_2d(&r, e->cv, (qaws_scalar)ss, QAWS_EVAL_FLAG_POSITION);
		sm_point3(c, r.position.x, r.position.y, q3);
		tt = sm_d3(x3, q3);
		if (tt > 0 && qaws_internal_sg_foot(c, e->cv, 1, e->lo, e->hi, x, &ss, &tt) && (best < 0 || tt < best))
		{
			best = tt; *s = ss; *t = tt;
		}
	}
	return best;
}

/* the edge whose front sweeps x first, and when: the least distance among
   the edges whose wedge holds x (at a convex corner nearer its own line
   than its neighbour's, at a reflex one farther); 0xffffffff when none.
   st holds each edge's foot (s, t), carried from call to call.
   shortcut: the wedges alone, no split events (a front cut off by another
   still counts beyond the cut); exact for convex polygons and reflex
   corners whose fronts are not split, upgrade to a wavefront for others */
static unsigned int mk_label(qaws_internal_sg_ctx const* c, mk_edge const* e, unsigned int ne, double const* x, double* st, double* d, double* out_t)
{
	unsigned int k, best = 0xffffffffu;
	double bt = 1e300, tol = c->tol;
	for (k = 0; k < ne; k++)
		d[k] = mk_dist(c, &e[k], x, &st[2 * k], &st[2 * k + 1]);
	for (k = 0; k < ne; k++)
	{
		double dp = d[e[k].prev], dn = d[e[k].next];
		if (d[k] < -tol || d[k] >= bt)
			continue;
		if (e[k].convex[0] ? d[k] > dp + tol : d[k] < dp - tol)
			continue;
		if (e[k].convex[1] ? d[k] > dn + tol : d[k] < dn - tol)
			continue;
		best = k; bt = d[k];
	}
	*out_t = bt;
	return best;
}

/* the uv curve of the geodesic from p along d (reversed: arriving at p),
   as long as it stays in the domain up to len (halved till it does) */
static qaws_curve* mk_ray(qaws_surface const* s, qaws_vec2 p, qaws_vec2 d, double len, int reversed)
{
	unsigned int k;
	for (k = 0; k < 8; k++, len *= 0.5)
	{
		qaws_surface_curve* sc = NULL;
		qaws_curve *a = NULL, *b = NULL;
		qaws_range rg;
		if (qaws_surface_curve_create_geodesic_ray(s, p, d, (qaws_scalar)len, 0, &sc) != QAWS_STATUS_OK)
			continue;
		rg = qaws_curve_get_parameter_range(qaws_surface_curve_get_uv(sc));
		if (qaws_curve_extract(&a, qaws_surface_curve_get_uv(sc), rg.min_value, rg.max_value) == QAWS_STATUS_OK && reversed)
		{
			qaws_curve_reverse(&b, a);
			qaws_curve_destroy(a);
			a = b;
		}
		qaws_surface_curve_destroy(sc);
		if (a)
			return a;
	}
	return NULL;
}

qaws_status qaws_surface_region_mitered_skeleton(qaws_surface_region const* region, qaws_surface_medial_desc const* desc,
	qaws_surface_medial_axis** out_axis)
{
	qaws_internal_sg_ctx const* c;
	qaws_surface_medial_axis* ax = NULL;
	qaws_path_2d* loops = NULL;
	mk_edge* E = NULL;
	unsigned int np = 0, ne = 0, cape = 0, i, j, k, nv = 0, nt, *lab = NULL, *hash = NULL, hsize = 1, ncr = 0, capcr = 0;
	qaws_mesh_2d* mesh = NULL;
	qaws_surface_mesh_desc md;
	qaws_scalar area = 0;
	double L, D, *vst = NULL, *vt = NULL, *dd = NULL, *ws = NULL;
	sm_cross* cr = NULL;
	qaws_status st;
	if (!region || !out_axis)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_axis = NULL;
	c = &region->ctx;
	qaws_surface_region_get_area(region, &area);
	L = desc && desc->resolution > 0 ? (double)desc->resolution : sqrt((double)area) / 40;
	D = 2 * sqrt((double)area) + L;
	ax = (qaws_surface_medial_axis*)calloc(1, sizeof(*ax));
	if (!ax)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* the edges: each loop's pieces between corners, smooth joins merged */
	st = sr_loops(region, &loops, &np);
	for (i = 0; i < np && st == QAWS_STATUS_OK; i++)
	{
		sm_feature* seq = NULL;
		double* turn = NULL;
		unsigned int ns = 0, caps = 0, c0 = 0, first = ne, q;
		for (j = 0; j < loops[i].curve_count && st == QAWS_STATUS_OK; j++)
		{
			qaws_curve const* cv = loops[i].curves[j];
			double* br = NULL;
			unsigned int nb = 0;
			if ((st = sr_breaks(cv, &br, &nb)) != QAWS_STATUS_OK)
				break;
			for (k = 0; k + 1 < nb && st == QAWS_STATUS_OK; k++)
			{
				if (!(br[k + 1] > br[k]))
					continue;
				/* shortcut: cuts along seams or poles are refused, upgrade by
				   joining the edge across the cut */
				if (sr_on_seam(c, cv, br[k], br[k + 1]) || sr_on_pole(c, cv, br[k], br[k + 1]))
					st = QAWS_STATUS_UNSUPPORTED_OPERATION;
				else if (!sr_grow((void**)&seq, &caps, ns + 1, sizeof(sm_feature)))
					st = QAWS_STATUS_ALLOCATION_FAILURE;
				else
				{
					seq[ns].point = 0; seq[ns].cv = cv; seq[ns].lo = br[k]; seq[ns].hi = br[k + 1];
					ns++;
				}
			}
			free(br);
		}
		turn = (double*)malloc(sizeof(double) * (ns + 1));
		if (st == QAWS_STATUS_OK && !turn)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		/* the turn at each piece's start (left: a convex corner) */
		for (j = 0; j < ns && st == QAWS_STATUS_OK; j++)
		{
			sm_feature const* cur = &seq[j];
			sm_feature const* prev = &seq[(j + ns - 1) % ns];
			qaws_eval_result_2d a, b;
			double da[2], db[2], aa = 0, ab = 0;
			sr_piece_end(prev->cv, prev->lo, prev->hi, &a);
			qaws_curve_evaluate_2d(&b, cur->cv, (qaws_scalar)cur->lo, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			da[0] = a.d1.x; da[1] = a.d1.y; db[0] = b.d1.x; db[1] = b.d1.y;
			qaws_internal_sg_angle(c, (double)b.position.x, (double)b.position.y, da, &aa);
			qaws_internal_sg_angle(c, (double)b.position.x, (double)b.position.y, db, &ab);
			turn[j] = remainder(ab - aa, 2 * 3.14159265358979323846);
			if (fabs(turn[j]) > SM_CORNER && fabs(turn[c0]) <= SM_CORNER)
				c0 = j;
		}
		/* shortcut: a loop with no corner (a smooth boundary) is refused; its
		   mitered skeleton is its medial axis */
		if (st == QAWS_STATUS_OK && (ns == 0 || fabs(turn[c0]) <= SM_CORNER))
			st = QAWS_STATUS_UNSUPPORTED_OPERATION;
		for (j = 0; j < ns && st == QAWS_STATUS_OK;)
		{
			/* the pieces from corner c0 + j to the next */
			qaws_curve* part[66];
			unsigned int n = 0, j0 = j, z;
			qaws_composite_desc cd;
			qaws_eval_result_2d r0, r1;
			sm_feature const* f0 = &seq[(c0 + j) % ns];
			sm_feature const* f1;
			qaws_vec2 back;
			do
			{
				sm_feature const* f = &seq[(c0 + j) % ns];
				if (n < 64 && qaws_curve_extract(&part[n + 1], f->cv, (qaws_scalar)f->lo, (qaws_scalar)f->hi) == QAWS_STATUS_OK)
					n++;
				j++;
			} while (j < ns && fabs(turn[(c0 + j) % ns]) <= SM_CORNER);
			if (!sr_grow((void**)&E, &cape, ne + 1, sizeof(mk_edge)))
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				for (z = 1; z <= n; z++) qaws_curve_destroy(part[z]);
				break;
			}
			/* carried on past both ends along the geodesics of their tangents */
			f1 = &seq[(c0 + j - 1) % ns];
			qaws_curve_evaluate_2d(&r0, f0->cv, (qaws_scalar)f0->lo, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
			sr_piece_end(f1->cv, f1->lo, f1->hi, &r1);
			back.x = -r0.d1.x; back.y = -r0.d1.y;
			part[0] = mk_ray(region->s, r0.position, back, D, 1);
			part[n + 1] = mk_ray(region->s, r1.position, r1.d1, D, 0);
			memset(&cd, 0, sizeof(cd));
			cd.dimension = QAWS_DIMENSION_2D;
			cd.segments = part[0] ? part : part + 1;
			cd.segment_count = n + (part[0] ? 1 : 0) + (part[n + 1] ? 1 : 0);
			E[ne].cv = NULL;
			st = qaws_curve_create_composite(&E[ne].cv, &cd);
			if (st != QAWS_STATUS_OK)
				for (z = 0; z <= n + 1; z++) qaws_curve_destroy(part[z]);
			else
			{
				qaws_range rg = qaws_curve_get_parameter_range(E[ne].cv);
				E[ne].lo = rg.min_value; E[ne].hi = rg.max_value;
				E[ne].corner[0] = r0.position.x; E[ne].corner[1] = r0.position.y;
				sm_point3(c, E[ne].corner[0], E[ne].corner[1], E[ne].corner + 2);
				E[ne].convex[0] = turn[(c0 + j0) % ns] > 0;
				E[ne].convex[1] = turn[(c0 + j) % ns] > 0;
				ne++;
			}
		}
		for (q = first; q < ne; q++)
		{
			E[q].prev = q == first ? ne - 1 : q - 1;
			E[q].next = q + 1 == ne ? first : q + 1;
		}
		free(seq); free(turn);
	}
	/* the mesh, each vertex labelled with the edge whose front sweeps it first */
	if (st == QAWS_STATUS_OK)
	{
		memset(&md, 0, sizeof(md));
		md.edge_length = (qaws_scalar)L;
		st = qaws_surface_region_triangulate(region, &md, &mesh);
	}
	if (st == QAWS_STATUS_OK)
	{
		nv = qaws_mesh_2d_get_vertex_count(mesh);
		lab = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
		vt = (double*)malloc(sizeof(double) * (nv + 1));
		vst = (double*)calloc(2 * ne * (nv + 1), sizeof(double));
		dd = (double*)malloc(sizeof(double) * (ne + 1));
		ws = (double*)malloc(sizeof(double) * 2 * (ne + 1));
		if (!lab || !vt || !vst || !dd || !ws)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < nv && st == QAWS_STATUS_OK; i++)
		{
			qaws_vec2 w = qaws_mesh_2d_get_vertex(mesh, i);
			double x[2], *seed = vst + 2 * ne * i;
			/* the feet carried from the last vertex */
			if (i) memcpy(seed, seed - 2 * ne, sizeof(double) * 2 * ne);
			x[0] = w.x; x[1] = w.y;
			lab[i] = mk_label(c, E, ne, x, seed, dd, &vt[i]);
		}
	}
	/* the crossings: on each mesh edge between differently labelled
	   vertices, where the label changes (bisection), each found once */
	if (st == QAWS_STATUS_OK)
	{
		nt = qaws_mesh_2d_get_triangle_count(mesh);
		while (hsize < 4 * nt + 16) hsize <<= 1;
		hash = (unsigned int*)malloc(sizeof(unsigned int) * hsize);
		if (!hash)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < hsize && st == QAWS_STATUS_OK; i++) hash[i] = 0xffffffffu;
		for (i = 0; i < nt && st == QAWS_STATUS_OK; i++)
		{
			unsigned int v[3], e, ids[3], nc = 0, lv[3];
			double cuv[3][2], crad[3];
			qaws_mesh_2d_get_triangle(mesh, i, v);
			for (e = 0; e < 3; e++)
			{
				unsigned int a = v[e] < v[(e + 1) % 3] ? v[e] : v[(e + 1) % 3], b = v[e] < v[(e + 1) % 3] ? v[(e + 1) % 3] : v[e];
				unsigned int slot = (unsigned int)((a * 2654435761u) ^ (b * 40503u)) & (hsize - 1), id = 0xffffffffu;
				while (hash[slot] != 0xffffffffu)
				{
					if (cr[hash[slot]].a == a && cr[hash[slot]].b == b) { id = hash[slot]; break; }
					slot = (slot + 1) & (hsize - 1);
				}
				if (id == 0xffffffffu)
				{
					sm_cross* x;
					if (!sr_grow((void**)&cr, &capcr, ncr + 1, sizeof(sm_cross)))
					{
						st = QAWS_STATUS_ALLOCATION_FAILURE;
						break;
					}
					x = &cr[ncr];
					memset(x, 0, sizeof(*x));
					x->a = a; x->b = b;
					if (lab[a] != 0xffffffffu && lab[b] != 0xffffffffu && lab[a] != lab[b])
					{
						qaws_vec2 pa = qaws_mesh_2d_get_vertex(mesh, a), pb = qaws_mesh_2d_get_vertex(mesh, b);
						double l0 = 0, l1 = 1, len = hypot((double)pb.x - pa.x, (double)pb.y - pa.y), tm = 0;
						unsigned int it;
						memcpy(ws, vst + 2 * ne * a, sizeof(double) * 2 * ne);
						for (it = 0; it < 60 && (l1 - l0) * len > 0.05 * c->tol; it++)
						{
							double xm[2], lm = 0.5 * (l0 + l1);
							xm[0] = pa.x + lm * ((double)pb.x - pa.x); xm[1] = pa.y + lm * ((double)pb.y - pa.y);
							if (mk_label(c, E, ne, xm, ws, dd, &tm) == lab[a]) l0 = lm; else l1 = lm;
						}
						x->uv[0] = pa.x + 0.5 * (l0 + l1) * ((double)pb.x - pa.x);
						x->uv[1] = pa.y + 0.5 * (l0 + l1) * ((double)pb.y - pa.y);
						mk_label(c, E, ne, x->uv, ws, dd, &x->r);
						x->ok = x->r < 1e299;
						/* at the boundary (feet so near a line do not settle): the
						   skeleton meets it only at corners, so the corner within a
						   mesh edge, else nothing */
						if (x->ok && x->r < 1e-3 * L)
						{
							double p3[3], dc = 1e300;
							unsigned int q, kc = 0;
							sm_point3(c, x->uv[0], x->uv[1], p3);
							for (q = 0; q < ne; q++)
								if (sm_d3(p3, E[q].corner + 2) < dc) { dc = sm_d3(p3, E[q].corner + 2); kc = q; }
							x->ok = dc <= 2 * L;
							x->uv[0] = E[kc].corner[0]; x->uv[1] = E[kc].corner[1]; x->r = 0;
						}
						x->fa = lab[a]; x->fb = lab[b];
						if (x->ok && x->r > ax->best[2]) { ax->best[0] = x->uv[0]; ax->best[1] = x->uv[1]; ax->best[2] = x->r; }
					}
					hash[slot] = ncr;
					id = ncr++;
				}
				ids[e] = id;
			}
			if (st != QAWS_STATUS_OK)
				break;
			for (e = 0; e < 3; e++)
			{
				sm_cross const* x = &cr[ids[e]];
				lv[e] = lab[v[e]];
				if (x->ok)
				{
					cuv[nc][0] = x->uv[0]; cuv[nc][1] = x->uv[1]; crad[nc] = x->r;
					nc++;
				}
			}
			if (nc == 2)
				st = sm_add_seg(ax, cuv[0], crad[0], cuv[1], crad[1]);
			else if (nc == 3 && lv[0] != lv[1] && lv[1] != lv[2] && lv[0] != lv[2])
			{
				/* three fronts meet: equally far from the three lines, by Newton
				   from the crossings' centre (kept there if it does not settle) */
				double p[2], sv[6], bp[2], br = (crad[0] + crad[1] + crad[2]) / 3, bres = 1e300;
				double span = hypot(cuv[1][0] - cuv[0][0], cuv[1][1] - cuv[0][1]) + hypot(cuv[2][0] - cuv[0][0], cuv[2][1] - cuv[0][1]);
				unsigned int itn, a2, q, dup = 0;
				p[0] = bp[0] = (cuv[0][0] + cuv[1][0] + cuv[2][0]) / 3;
				p[1] = bp[1] = (cuv[0][1] + cuv[1][1] + cuv[2][1]) / 3;
				for (q = 0; q < 3; q++) { sv[2 * q] = vst[2 * ne * v[q] + 2 * lv[q]]; sv[2 * q + 1] = vst[2 * ne * v[q] + 2 * lv[q] + 1]; }
				for (itn = 0; itn < 20; itn++)
				{
					double d3[3], F[2], J[4], det, res;
					for (q = 0; q < 3; q++) d3[q] = mk_dist(c, &E[lv[q]], p, &sv[2 * q], &sv[2 * q + 1]);
					if (!(d3[0] > -1e299 && d3[1] > -1e299 && d3[2] > -1e299))
						break;
					F[0] = d3[0] - d3[1]; F[1] = d3[0] - d3[2];
					res = fabs(F[0]) + fabs(F[1]);
					/* (only near the triangle: a far root is another meeting) */
					if (res < bres && hypot(p[0] - cuv[0][0], p[1] - cuv[0][1]) <= 2 * span + 4 * c->tol)
					{
						bres = res; bp[0] = p[0]; bp[1] = p[1]; br = (d3[0] + d3[1] + d3[2]) / 3;
					}
					if (res <= 1e-3 * c->tol)
						break;
					for (a2 = 0; a2 < 2; a2++)
					{
						double px[2], e3[3], s2[6], hs = SM_FD_STEP * c->diag;
						memcpy(s2, sv, sizeof(s2));
						px[0] = p[0] + (a2 == 0 ? hs : 0); px[1] = p[1] + (a2 == 1 ? hs : 0);
						for (q = 0; q < 3; q++) e3[q] = mk_dist(c, &E[lv[q]], px, &s2[2 * q], &s2[2 * q + 1]);
						J[a2] = ((e3[0] - e3[1]) - F[0]) / hs;
						J[2 + a2] = ((e3[0] - e3[2]) - F[1]) / hs;
					}
					det = J[0] * J[3] - J[1] * J[2];
					if (!(fabs(det) > 0))
						break;
					p[0] -= (J[3] * F[0] - J[1] * F[1]) / det;
					p[1] -= (J[0] * F[1] - J[2] * F[0]) / det;
				}
				/* found already from a neighbouring triangle: the same junction */
				for (q = 0; q < ax->njun && !dup; q++)
					if (hypot(ax->junction[3 * q] - bp[0], ax->junction[3 * q + 1] - bp[1]) <= 100 * c->tol)
					{
						bp[0] = ax->junction[3 * q]; bp[1] = ax->junction[3 * q + 1]; br = ax->junction[3 * q + 2];
						dup = 1;
					}
				if (!dup)
				{
					if (!sr_grow((void**)&ax->junction, &ax->capjun, 3 * (ax->njun + 1), sizeof(double)))
						st = QAWS_STATUS_ALLOCATION_FAILURE;
					else
					{
						ax->junction[3 * ax->njun] = bp[0]; ax->junction[3 * ax->njun + 1] = bp[1]; ax->junction[3 * ax->njun + 2] = br;
						ax->njun++;
						if (br > ax->best[2]) { ax->best[0] = bp[0]; ax->best[1] = bp[1]; ax->best[2] = br; }
					}
				}
				for (q = 0; q < 3 && st == QAWS_STATUS_OK; q++)
					st = sm_add_seg(ax, bp, br, cuv[q], crad[q]);
			}
		}
	}
	for (i = 0; i < ne; i++) qaws_curve_destroy(E[i].cv);
	free(E); free(lab); free(vt); free(vst); free(dd); free(ws); free(hash); free(cr);
	qaws_mesh_2d_destroy(mesh);
	free(loops);
	if (st != QAWS_STATUS_OK)
	{
		qaws_surface_medial_axis_destroy(ax);
		return st;
	}
	*out_axis = ax;
	return QAWS_STATUS_OK;
}
