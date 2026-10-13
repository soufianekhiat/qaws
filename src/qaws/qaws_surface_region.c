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
	qaws_curve_get_control_points(cv, probe, 0, &n);
	if (kind == QAWS_CURVE_KIND_BSPLINE)
		qaws_bspline_get_knots(cv, probe, 0, &nk);
	cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (2 * n + nk + 2));
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if ((st = qaws_curve_get_control_points(cv, cp, n, &n)) == QAWS_STATUS_OK && nk)
	{
		kn = cp + 2 * n;
		st = qaws_bspline_get_knots(cv, kn, nk, &nk);
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
			st = qaws_curve_create_bezier(&d, out);
		}
		else
		{
			qaws_bspline_desc d;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.degree = deg; d.control_points = cp; d.control_point_count = n;
			d.knots = kn; d.knot_count = nk;
			st = qaws_curve_create_bspline(&d, out);
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
		qaws_curve_evaluate_2d(p->curves[0], r0.min_value, QAWS_EVAL_FLAG_POSITION, &a);
		qaws_curve_evaluate_2d(p->curves[p->curve_count - 1], r1.max_value, QAWS_EVAL_FLAG_POSITION, &b);
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
		if ((st = qaws_curve_reverse(off, &x)) != QAWS_STATUS_OK || (st = sr_keep(pl, x)) != QAWS_STATUS_OK) return st;
		if ((st = qaws_curve_reverse(r0, &y)) != QAWS_STATUS_OK || (st = sr_keep(pl, y)) != QAWS_STATUS_OK) return st;
		loop[0] = base; loop[1] = r1; loop[2] = x; loop[3] = y;
	}
	else
	{
		/* along the offset, in to the base's end, back along it, out */
		if ((st = qaws_curve_reverse(r1, &x)) != QAWS_STATUS_OK || (st = sr_keep(pl, x)) != QAWS_STATUS_OK) return st;
		if ((st = qaws_curve_reverse(base, &y)) != QAWS_STATUS_OK || (st = sr_keep(pl, y)) != QAWS_STATUS_OK) return st;
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
		qaws_curve_evaluate_2d(cv, (qaws_scalar)(ta + (tb - ta) * k / 4), QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e);
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
	qaws_curve_evaluate_2d(cv, (qaws_scalar)(t - h), QAWS_EVAL_FLAG_D1, &a);
	qaws_curve_evaluate_2d(cv, (qaws_scalar)(t + h), QAWS_EVAL_FLAG_D1, &b);
	la = hypot((double)a.d1.x, (double)a.d1.y);
	lb = hypot((double)b.d1.x, (double)b.d1.y);
	if (!(la > 0) || !(lb > 0))
		return 1;
	cr = ((double)a.d1.x * b.d1.y - (double)a.d1.y * b.d1.x) / (la * lb);
	dt = ((double)a.d1.x * b.d1.x + (double)a.d1.y * b.d1.y) / (la * lb);
	return atan2(fabs(cr), dt) > ang;
}

/* the curve in smooth pieces (split at the knots of a B-spline where it
   turns a corner), less its cuts along seams and poles: each piece's strips on
   the sides asked (bit 1 left, bit 2 right), circles at its ends */
static qaws_status sr_tube_curve(sr_tube* tb, qaws_curve const* cv, int sides)
{
	qaws_range rg = qaws_curve_get_parameter_range(cv);
	qaws_curve_kind kind = qaws_curve_get_kind(cv);
	double* br = NULL;
	unsigned int nb = 0, i;
	qaws_status st = QAWS_STATUS_OK;
	br = (double*)malloc(sizeof(double) * 2);
	if (!br)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	br[nb++] = (double)rg.min_value;
	if (kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS)
	{
		unsigned int nk = 0, deg = qaws_curve_get_degree(cv);
		qaws_scalar probe[2];
		qaws_bspline_get_knots(cv, probe, 0, &nk);
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
			if (qaws_bspline_get_knots(cv, kn, nk, &nk) == QAWS_STATUS_OK)
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
			if ((st = qaws_curve_extract(cv, (qaws_scalar)ta, (qaws_scalar)tb2, &piece)) != QAWS_STATUS_OK || (st = sr_keep(&tb->raw, piece)) != QAWS_STATUS_OK)
				break;
			base = piece;
		}
		if (sides & 1) st = sr_strip(tb, base, tb->r);
		if (st == QAWS_STATUS_OK && (sides & 2)) st = sr_strip(tb, base, -tb->r);
		if (st == QAWS_STATUS_OK)
		{
			qaws_range pr = qaws_curve_get_parameter_range(base);
			qaws_curve_evaluate_2d(base, pr.min_value, QAWS_EVAL_FLAG_POSITION, &e0);
			qaws_curve_evaluate_2d(base, pr.max_value, QAWS_EVAL_FLAG_POSITION, &e1);
			st = sr_disk(tb, e0.position);
			if (st == QAWS_STATUS_OK) st = sr_disk(tb, e1.position);
		}
	}
	free(br);
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
				st = sr_tube_curve(&tb, p.curves[k], d > 0 ? 2 : 1);
		}
		for (i = 0; i < path_count && st == QAWS_STATUS_OK; i++)
			if (!paths[i].closed)
				for (k = 0; k < paths[i].curve_count && st == QAWS_STATUS_OK; k++)
					st = sr_tube_curve(&tb, qaws_surface_curve_get_uv(paths[i].curves[k]), 3);
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
