/*
 * 2.5D stacks of parallel contours: sections by matched-contour blending or
 * by blended signed distances, and Booleans level by level.
 */

#include "qaws_stack.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_curve_batch.h"
#include "qaws_bezier.h"
#include "qaws_rational_bezier.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_catmull_rom.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_parallel.h"
#include "internal/qaws_internal_clip.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SK_NONE 0xFFFFFFFFu
#define SK_RESAMPLE 128

/* ======================================================================== */
/*  Weights of the levels round z                                           */
/* ======================================================================== */

/* levels idx[0..n) and their weights at z (LINEAR: 2, CUBIC: up to 4);
   returns 0 when z is outside the stack, -1 when z is exactly level *at */
static int sk_weights(qaws_stack_2d const* st, double z, unsigned int* idx, double* w, unsigned int* n, unsigned int* at)
{
	unsigned int L = st->level_count, k;
	double s, d;
	if (L == 0 || z < st->levels[0].z || z > st->levels[L - 1].z)
		return 0;
	for (k = 0; k < L; k++)
		if ((double)st->levels[k].z == z)
		{
			*at = k;
			return -1;
		}
	for (k = 0; k + 1 < L && !(z < st->levels[k + 1].z); k++)
		;
	d = (double)st->levels[k + 1].z - (double)st->levels[k].z;
	s = (z - (double)st->levels[k].z) / d;
	if (st->interp == QAWS_STACK_LINEAR || L < 3)
	{
		idx[0] = k; w[0] = 1 - s;
		idx[1] = k + 1; w[1] = s;
		*n = 2;
		return 1;
	}
	{
		/* Hermite with finite-difference tangents on uneven heights:
		   value = h00 Q_k + h10 d m_k + h01 Q_k+1 + h11 d m_k+1 */
		double h00 = 2 * s * s * s - 3 * s * s + 1, h10 = s * s * s - 2 * s * s + s;
		double h01 = -2 * s * s * s + 3 * s * s, h11 = s * s * s - s * s;
		double wq[4] = { 0, 0, 0, 0 };   /* levels k - 1, k, k + 1, k + 2 */
		unsigned int j;
		wq[1] += h00;
		wq[2] += h01;
		/* m_k = (Q_k+1 - Q_k-1) / (z_k+1 - z_k-1), or one-sided at the bottom */
		if (k > 0)
		{
			double a = d / ((double)st->levels[k + 1].z - (double)st->levels[k - 1].z);
			wq[2] += h10 * a;
			wq[0] -= h10 * a;
		}
		else
		{
			wq[2] += h10;
			wq[1] -= h10;
		}
		if (k + 2 < L)
		{
			double a = d / ((double)st->levels[k + 2].z - (double)st->levels[k].z);
			wq[3] += h11 * a;
			wq[1] -= h11 * a;
		}
		else
		{
			wq[2] += h11;
			wq[1] -= h11;
		}
		*n = 0;
		for (j = 0; j < 4; j++)
		{
			int lv = (int)k - 1 + (int)j;
			if (lv < 0 || lv >= (int)L || wq[j] == 0)
				continue;
			idx[*n] = (unsigned int)lv;
			w[*n] = wq[j];
			(*n)++;
		}
		return 1;
	}
}

/* ======================================================================== */
/*  Regions                                                                 */
/* ======================================================================== */

static qaws_status sk_union(qaws_path_2d const* paths, unsigned int n, qaws_fill_rule fr, qaws_clip_result** out)
{
	if (!n)
		return qaws_internal_clip_result_empty(out);
	return qaws_clip_boolean(QAWS_CLIP_UNION, fr, paths, n, NULL, 0, out);
}

static void sk_center(qaws_path_2d const* p, double* c, double* size)
{
	qaws_vec2 lo, hi;
	if (qaws_path_compute_bounds_2d(p, &lo, &hi) != QAWS_STATUS_OK)
	{
		c[0] = c[1] = 0;
		*size = 0;
		return;
	}
	c[0] = 0.5 * (lo.x + hi.x);
	c[1] = 0.5 * (lo.y + hi.y);
	*size = hypot(hi.x - lo.x, hi.y - lo.y);
}

/* ======================================================================== */
/*  MATCH: contours blended                                                 */
/* ======================================================================== */

typedef struct sk_struct
{
	qaws_curve_kind kind;
	unsigned int degree, n;
	qaws_scalar const* cp;
	qaws_scalar const* w;       /* NULL: polynomial */
	qaws_scalar const* knots;   /* NULL: Bezier kinds */
	unsigned int knot_count;
} sk_struct;

static int sk_get_struct(qaws_curve const* c, sk_struct* s)
{
	memset(s, 0, sizeof(*s));
	s->kind = c->kind;
	s->degree = c->degree;
	switch (c->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
	{
		qaws_bezier_impl const* i = (qaws_bezier_impl const*)c->impl;
		s->cp = i->control_points; s->n = i->control_point_count;
		return 1;
	}
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	{
		qaws_rational_bezier_impl const* i = (qaws_rational_bezier_impl const*)c->impl;
		s->cp = i->control_points; s->n = i->control_point_count; s->w = i->weights;
		return 1;
	}
	case QAWS_CURVE_KIND_BSPLINE:
	{
		qaws_bspline_impl const* i = (qaws_bspline_impl const*)c->impl;
		s->cp = i->control_points; s->n = i->control_point_count; s->knots = i->knots; s->knot_count = i->knot_count;
		return 1;
	}
	case QAWS_CURVE_KIND_NURBS:
	{
		qaws_nurbs_impl const* i = (qaws_nurbs_impl const*)c->impl;
		s->cp = i->control_points; s->n = i->control_point_count; s->w = i->weights;
		s->knots = i->knots; s->knot_count = i->knot_count;
		return 1;
	}
	default:
		return 0;
	}
}

static int sk_same_struct(sk_struct const* a, sk_struct const* b)
{
	unsigned int i;
	if (a->kind != b->kind || a->degree != b->degree || a->n != b->n || a->knot_count != b->knot_count)
		return 0;
	for (i = 0; i < a->knot_count; i++)
		if (fabs((double)a->knots[i] - (double)b->knots[i]) > 1e-12 * (fabs((double)a->knots[i]) + 1))
			return 0;
	for (i = 0; a->w && i < a->n; i++)
		if (fabs((double)a->w[i] - (double)b->w[i]) > 1e-12 * fabs((double)a->w[i]))
			return 0;
	return 1;
}

/* a point contour: a curve of (almost) no extent */
static int sk_is_point(qaws_path_2d const* p, double ref_size, double* pt)
{
	double c[2], size;
	sk_center(p, c, &size);
	pt[0] = c[0];
	pt[1] = c[1];
	return size <= 1e-9 * (ref_size + 1e-300);
}

/* the blend of same-structure contours (points take the template's
   structure): control points (homogeneous for rational kinds) */
static qaws_status sk_blend_struct(sk_struct const* tpl, qaws_curve const* const* cs, double const* pts, int const* is_pt,
	double const* w, unsigned int n, qaws_curve** out)
{
	unsigned int dim = 2, i, k;
	qaws_scalar* cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (tpl->n * dim + tpl->n));
	qaws_scalar* wt;
	qaws_status s;
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	wt = cp + tpl->n * dim;
	for (k = 0; k < tpl->n; k++)
	{
		double hx = 0, hy = 0, hw = 0;
		for (i = 0; i < n; i++)
		{
			double px, py, ww = tpl->w ? (double)tpl->w[k] : 1.0;
			if (is_pt[i])
			{
				px = pts[2 * i];
				py = pts[2 * i + 1];
			}
			else
			{
				sk_struct si;
				sk_get_struct(cs[i], &si);
				px = (double)si.cp[2 * k];
				py = (double)si.cp[2 * k + 1];
				if (si.w) ww = (double)si.w[k];
			}
			hx += w[i] * ww * px;
			hy += w[i] * ww * py;
			hw += w[i] * ww;
		}
		if (!(hw > 0))
		{
			free(cp);
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
		}
		cp[2 * k] = (qaws_scalar)(hx / hw);
		cp[2 * k + 1] = (qaws_scalar)(hy / hw);
		wt[k] = (qaws_scalar)hw;
	}
	switch (tpl->kind)
	{
	case QAWS_CURVE_KIND_BEZIER:
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = tpl->degree; d.control_points = cp; d.control_point_count = tpl->n;
		s = qaws_curve_create_bezier(&d, out);
		break;
	}
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = tpl->degree; d.control_points = cp; d.control_point_count = tpl->n;
		d.weights = wt; d.weight_count = tpl->n;
		s = qaws_curve_create_rational_bezier(&d, out);
		break;
	}
	case QAWS_CURVE_KIND_BSPLINE:
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = tpl->degree; d.control_points = cp; d.control_point_count = tpl->n;
		d.knots = tpl->knots; d.knot_count = tpl->knot_count;
		s = qaws_curve_create_bspline(&d, out);
		break;
	}
	default:
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = tpl->degree; d.control_points = cp; d.control_point_count = tpl->n;
		d.knots = tpl->knots; d.knot_count = tpl->knot_count; d.weights = wt; d.weight_count = tpl->n;
		s = qaws_curve_create_nurbs(&d, out);
		break;
	}
	}
	(void)k;
	free(cp);
	return s;
}

/* a closed curve resampled at m points of equal arc length (dense chords) */
static void sk_resample(qaws_curve const* c, unsigned int m, double* out)
{
	enum { DENSE = 4096 };
	static double const dummy = 0;
	double* acc = (double*)malloc(sizeof(double) * (DENSE + 1));
	double* xy = (double*)malloc(sizeof(double) * 2 * (DENSE + 1));
	qaws_range r = qaws_curve_get_parameter_range(c);
	unsigned int i, j = 0;
	(void)dummy;
	if (!acc || !xy) { free(acc); free(xy); memset(out, 0, sizeof(double) * 2 * m); return; }
	for (i = 0; i <= DENSE; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * i / DENSE), QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * i] = e.position.x;
		xy[2 * i + 1] = e.position.y;
		acc[i] = i ? acc[i - 1] + hypot(xy[2 * i] - xy[2 * i - 2], xy[2 * i + 1] - xy[2 * i - 1]) : 0;
	}
	for (i = 0; i < m; i++)
	{
		double target = acc[DENSE] * i / m, f;
		while (j + 1 < DENSE && acc[j + 1] < target) j++;
		f = acc[j + 1] > acc[j] ? (target - acc[j]) / (acc[j + 1] - acc[j]) : 0;
		out[2 * i] = xy[2 * j] + f * (xy[2 * j + 2] - xy[2 * j]);
		out[2 * i + 1] = xy[2 * j + 1] + f * (xy[2 * j + 3] - xy[2 * j + 1]);
	}
	free(acc);
	free(xy);
}

/* contours of different structures: resampled, aligned to the first's start
   and direction, blended, a closed centripetal Catmull-Rom through them */
static qaws_status sk_blend_resampled(qaws_curve const* const* cs, double const* pts, int const* is_pt, double const* w, unsigned int n,
	qaws_curve** out)
{
	unsigned int m = SK_RESAMPLE, i, k, ref = SK_NONE;
	double* all = (double*)malloc(sizeof(double) * 2 * m * (n + 1));
	double* blend = all + 2 * m * n;
	qaws_scalar* sp;
	qaws_catmull_rom_desc d;
	qaws_status s;
	if (!all)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n; i++)
	{
		double* p = &all[2 * m * i];
		if (is_pt[i])
		{
			for (k = 0; k < m; k++) { p[2 * k] = pts[2 * i]; p[2 * k + 1] = pts[2 * i + 1]; }
			continue;
		}
		sk_resample(cs[i], m, p);
		if (ref == SK_NONE)
		{
			ref = i;
			continue;
		}
		{
			/* same direction as the reference (signed areas), then the shift
			   with the smallest squared distance */
			double* rp = &all[2 * m * ref];
			double a0 = 0, a1 = 0, best = 1e300;
			unsigned int sh, bs = 0, q;
			for (k = 0; k < m; k++)
			{
				unsigned int k2 = (k + 1) % m;
				a0 += rp[2 * k] * rp[2 * k2 + 1] - rp[2 * k2] * rp[2 * k + 1];
				a1 += p[2 * k] * p[2 * k2 + 1] - p[2 * k2] * p[2 * k + 1];
			}
			if ((a0 < 0) != (a1 < 0))
				for (k = 1; k < m / 2 + (m & 1); k++)
				{
					unsigned int k2 = m - k;
					double tx = p[2 * k], ty = p[2 * k + 1];
					if (k2 <= k) break;
					p[2 * k] = p[2 * k2]; p[2 * k + 1] = p[2 * k2 + 1];
					p[2 * k2] = tx; p[2 * k2 + 1] = ty;
				}
			for (sh = 0; sh < m; sh++)
			{
				double e = 0;
				for (q = 0; q < m && e < best; q++)
				{
					double dx = p[2 * ((q + sh) % m)] - rp[2 * q], dy = p[2 * ((q + sh) % m) + 1] - rp[2 * q + 1];
					e += dx * dx + dy * dy;
				}
				if (e < best) { best = e; bs = sh; }
			}
			if (bs)
			{
				double* tmp = (double*)malloc(sizeof(double) * 2 * m);
				if (!tmp) { free(all); return QAWS_STATUS_ALLOCATION_FAILURE; }
				for (q = 0; q < m; q++) { tmp[2 * q] = p[2 * ((q + bs) % m)]; tmp[2 * q + 1] = p[2 * ((q + bs) % m) + 1]; }
				memcpy(p, tmp, sizeof(double) * 2 * m);
				free(tmp);
			}
		}
	}
	for (k = 0; k < 2 * m; k++)
	{
		double v = 0;
		for (i = 0; i < n; i++)
			v += w[i] * all[2 * m * i + k];
		blend[k] = v;
	}
	sp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * m);
	if (!sp) { free(all); return QAWS_STATUS_ALLOCATION_FAILURE; }
	for (k = 0; k < 2 * m; k++) sp[k] = (qaws_scalar)blend[k];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.control_points = sp; d.control_point_count = m;
	d.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL; d.closed = 1;
	s = qaws_curve_create_catmull_rom(&d, out);
	free(sp);
	free(all);
	return s;
}

/* MATCH: 1 when it applied (*out set), 0 when the levels do not pair up */
static qaws_status sk_match(qaws_stack_2d const* st, unsigned int const* idx, double const* w, unsigned int n,
	qaws_clip_result** out, int* applied)
{
	qaws_stack_level const* L0 = &st->levels[idx[0]];
	unsigned int m = L0->path_count, i, j, k;
	qaws_curve** blended;
	qaws_curve const** views;
	qaws_path_2d* paths;
	unsigned int* pairs;
	double ref_size = 0;
	qaws_status s = QAWS_STATUS_OK;
	*applied = 0;
	for (i = 0; i < n; i++)
	{
		qaws_stack_level const* L = &st->levels[idx[i]];
		if (L->path_count != m)
			return QAWS_STATUS_OK;
		for (j = 0; j < m; j++)
			if (L->paths[j].curve_count != 1)
				return QAWS_STATUS_OK;
	}
	pairs = (unsigned int*)malloc(sizeof(unsigned int) * (n * m + 1));
	blended = (qaws_curve**)calloc(m + 1, sizeof(qaws_curve*));
	views = (qaws_curve const**)calloc(m + 1, sizeof(qaws_curve*));
	paths = (qaws_path_2d*)calloc(m + 1, sizeof(qaws_path_2d));
	if (!pairs || !blended || !views || !paths)
	{
		free(pairs); free(blended); free((void*)views); free(paths);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			double c[2], sz;
			sk_center(&st->levels[idx[i]].paths[j], c, &sz);
			if (sz > ref_size) ref_size = sz;
		}
	/* pair every level's contours with the first's: greedy nearest centres */
	for (i = 0; i < n; i++)
	{
		unsigned char* used = (unsigned char*)calloc(m + 1, 1);
		for (j = 0; j < m; j++)
		{
			double cj[2], sj, best = 1e300;
			unsigned int bk = 0;
			sk_center(&L0->paths[j], cj, &sj);
			for (k = 0; k < m; k++)
			{
				double ck[2], sk, dd;
				if (used && used[k]) continue;
				sk_center(&st->levels[idx[i]].paths[k], ck, &sk);
				dd = hypot(ck[0] - cj[0], ck[1] - cj[1]);
				if (dd < best) { best = dd; bk = k; }
			}
			pairs[i * m + j] = bk;
			if (used) used[bk] = 1;
		}
		free(used);
	}
	for (j = 0; j < m && s == QAWS_STATUS_OK; j++)
	{
		qaws_curve const* cs[4];
		double pts[8];
		int is_pt[4];
		sk_struct tpl, si;
		int have_tpl = 0, same = 1;
		for (i = 0; i < n; i++)
		{
			qaws_path_2d const* p = &st->levels[idx[i]].paths[pairs[i * m + j]];
			cs[i] = p->curves[0];
			is_pt[i] = sk_is_point(p, ref_size, &pts[2 * i]);
			if (is_pt[i])
				continue;
			if (!sk_get_struct(cs[i], &si))
				same = 0;
			else if (!have_tpl)
			{
				tpl = si;
				have_tpl = 1;
			}
			else if (!sk_same_struct(&tpl, &si))
				same = 0;
		}
		if (!have_tpl)
		{
			/* points only: nothing to enclose */
			continue;
		}
		s = QAWS_STATUS_UNSUPPORTED_OPERATION;
		if (same)
			s = sk_blend_struct(&tpl, cs, pts, is_pt, w, n, &blended[j]);
		if (s == QAWS_STATUS_UNSUPPORTED_OPERATION)
			s = sk_blend_resampled(cs, pts, is_pt, w, n, &blended[j]);
	}
	if (s == QAWS_STATUS_OK)
	{
		unsigned int np = 0;
		for (j = 0; j < m; j++)
			if (blended[j])
			{
				views[np] = blended[j];
				paths[np].curves = &views[np];
				paths[np].curve_count = 1;
				paths[np].closed = 1;
				np++;
			}
		s = sk_union(paths, np, st->fill_rule, out);
		*applied = s == QAWS_STATUS_OK;
	}
	for (j = 0; j < m; j++) qaws_curve_destroy(blended[j]);
	free(pairs); free(blended); free((void*)views); free(paths);
	return s;
}

/* ======================================================================== */
/*  DISTANCE: blended signed distances, marching squares                    */
/* ======================================================================== */

typedef struct sk_field
{
	unsigned int nx, ny;
	double x0, y0, h;
	double* f;        /* (nx + 1) (ny + 1) nodes */
} sk_field;

/* inside flags of the nodes of one row y for a level: crossings of its
   flattened curves with the row, wound by the fill rule */
static qaws_status sk_level_signs(qaws_stack_level const* L, qaws_fill_rule fr, sk_field const* g, double flat, unsigned char* inside)
{
	qaws_flat_seg* segs = NULL;
	unsigned int ns = 0, cap = 0, i, j, c;
	qaws_status s = QAWS_STATUS_OK;
	for (i = 0; i < L->path_count && s == QAWS_STATUS_OK; i++)
		for (c = 0; c < L->paths[i].curve_count && s == QAWS_STATUS_OK; c++)
			s = qaws_internal_flatten_curve(L->paths[i].curves[c], 2, (qaws_scalar)flat, i, &segs, &ns, &cap);
	if (s != QAWS_STATUS_OK) { free(segs); return s; }
	for (j = 0; j <= g->ny; j++)
	{
		double y = g->y0 + g->h * j;
		/* crossings of the row (half-open in y), as (x, +1 / -1) */
		double* xs = (double*)malloc(sizeof(double) * 2 * (ns + 1));
		unsigned int nc = 0, a, b;
		int wsum = 0;
		if (!xs) { free(segs); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < ns; i++)
		{
			double y0 = segs[i].p0[1], y1 = segs[i].p1[1];
			if ((y0 <= y) != (y1 <= y))
			{
				double t = (y - y0) / (y1 - y0);
				xs[2 * nc] = segs[i].p0[0] + t * (segs[i].p1[0] - segs[i].p0[0]);
				xs[2 * nc + 1] = y1 > y0 ? 1 : -1;
				nc++;
			}
		}
		/* sort crossings by x (insertion: rows have few) */
		for (a = 1; a < nc; a++)
			for (b = a; b > 0 && xs[2 * b] < xs[2 * b - 2]; b--)
			{
				double t0 = xs[2 * b], t1 = xs[2 * b + 1];
				xs[2 * b] = xs[2 * b - 2]; xs[2 * b + 1] = xs[2 * b - 1];
				xs[2 * b - 2] = t0; xs[2 * b - 1] = t1;
			}
		/* the winding number of a node: the crossings of the ray to its
		   right, upward +1 and downward -1 (all minus those to its left) */
		{
			int total = 0;
			for (a = 0; a < nc; a++) total += (int)xs[2 * a + 1];
			a = 0;
			for (i = 0; i <= g->nx; i++)
			{
				double x = g->x0 + g->h * i;
				while (a < nc && xs[2 * a] < x) { wsum += (int)xs[2 * a + 1]; a++; }
				inside[j * (g->nx + 1) + i] = (unsigned char)qaws_fill_rule_inside(fr, total - wsum);
			}
		}
		free(xs);
	}
	free(segs);
	return QAWS_STATUS_OK;
}

static qaws_status sk_distance(qaws_stack_2d const* st, unsigned int const* idx, double const* w, unsigned int n, qaws_clip_result** out)
{
	sk_field g;
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, ext;
	unsigned int i, j, k, nn, G = st->grid ? st->grid : 192;
	qaws_scalar* pts = NULL;
	qaws_closest_point* cp = NULL;
	unsigned char* ins = NULL;
	qaws_status s = QAWS_STATUS_OK;
	memset(&g, 0, sizeof(g));
	for (i = 0; i < n; i++)
	{
		qaws_stack_level const* L = &st->levels[idx[i]];
		for (j = 0; j < L->path_count; j++)
		{
			qaws_vec2 a, b;
			if (L->paths[j].curve_count && qaws_path_compute_bounds_2d(&L->paths[j], &a, &b) == QAWS_STATUS_OK)
			{
				if (a.x < lo[0]) lo[0] = a.x;
				if (a.y < lo[1]) lo[1] = a.y;
				if (b.x > hi[0]) hi[0] = b.x;
				if (b.y > hi[1]) hi[1] = b.y;
			}
		}
	}
	if (!(hi[0] >= lo[0]))
		return qaws_internal_clip_result_empty(out);
	ext = (hi[0] - lo[0] > hi[1] - lo[1] ? hi[0] - lo[0] : hi[1] - lo[1]);
	if (!(ext > 0)) ext = 1;
	g.h = ext * 1.1 / G;
	g.x0 = 0.5 * (lo[0] + hi[0]) - 0.55 * ext - 0.5 * g.h;
	g.y0 = 0.5 * (lo[1] + hi[1]) - 0.55 * ext - 0.5 * g.h;
	g.nx = g.ny = G + 1;
	nn = (g.nx + 1) * (g.ny + 1);
	g.f = (double*)calloc(nn, sizeof(double));
	pts = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nn);
	cp = (qaws_closest_point*)malloc(sizeof(qaws_closest_point) * nn);
	ins = (unsigned char*)malloc(nn);
	if (!g.f || !pts || !cp || !ins)
	{
		s = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (j = 0; j <= g.ny; j++)
		for (i = 0; i <= g.nx; i++)
		{
			pts[2 * (j * (g.nx + 1) + i)] = (qaws_scalar)(g.x0 + g.h * i);
			pts[2 * (j * (g.nx + 1) + i) + 1] = (qaws_scalar)(g.y0 + g.h * j);
		}
	/* F = sum w_i d_i, d_i the signed distance to level i (negative inside) */
	for (k = 0; k < n && s == QAWS_STATUS_OK; k++)
	{
		qaws_stack_level const* L = &st->levels[idx[k]];
		qaws_curve const** curves;
		unsigned int nc = 0, q;
		qaws_curve_batch_desc bd;
		qaws_curve_set* set = NULL;
		for (q = 0; q < L->path_count; q++) nc += L->paths[q].curve_count;
		curves = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (nc + 1));
		if (!curves) { s = QAWS_STATUS_ALLOCATION_FAILURE; break; }
		nc = 0;
		for (q = 0; q < L->path_count; q++)
			for (j = 0; j < L->paths[q].curve_count; j++)
				curves[nc++] = L->paths[q].curves[j];
		if (nc == 0)
		{
			/* an empty level: far outside everywhere */
			for (q = 0; q < nn; q++) g.f[q] += w[k] * 4 * ext;
			free((void*)curves);
			continue;
		}
		memset(&bd, 0, sizeof(bd));
		bd.curves = curves;
		bd.curve_count = nc;
		s = qaws_curve_set_create(&bd, &set);
		if (s == QAWS_STATUS_OK) s = qaws_curve_set_find_closest(set, pts, nn, 0, cp, NULL);
		if (s == QAWS_STATUS_OK) s = sk_level_signs(L, st->fill_rule, &g, ext / 4096, ins);
		if (s == QAWS_STATUS_OK)
			for (q = 0; q < nn; q++)
			{
				double d = cp[q].curve == QAWS_CURVE_BATCH_NONE ? 4 * ext : (double)cp[q].distance;
				g.f[q] += w[k] * (ins[q] ? -d : d);
			}
		qaws_curve_set_destroy(set);
		free((void*)curves);
	}
	if (s != QAWS_STATUS_OK)
		goto done;
	/* marching squares: one closed polyline per component of { F <= 0 } */
	{
		unsigned int W = g.nx, Hc = g.ny;
		/* edge points: horizontal edges (i, j)-(i + 1, j) id j W + i, vertical edges after */
		unsigned int nh = W * (Hc + 1), ne = nh + (W + 1) * Hc, e;
		double* ep = (double*)malloc(sizeof(double) * 2 * ne);
		unsigned int* nxt = (unsigned int*)malloc(sizeof(unsigned int) * ne);
		unsigned char* seen = (unsigned char*)calloc(ne, 1);
		qaws_curve** loops = NULL;
		unsigned int nloops = 0, caploops = 0;
		if (!ep || !nxt || !seen)
		{
			free(ep); free(nxt); free(seen);
			s = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (e = 0; e < ne; e++) nxt[e] = SK_NONE;
#define SK_F(ii, jj) g.f[(jj) * (g.nx + 1) + (ii)]
		for (j = 0; j < Hc; j++)
			for (i = 0; i < W; i++)
			{
				/* corners: 0 (i, j), 1 (i + 1, j), 2 (i + 1, j + 1), 3 (i, j + 1); in = F <= 0 */
				double f0 = SK_F(i, j), f1 = SK_F(i + 1, j), f2 = SK_F(i + 1, j + 1), f3 = SK_F(i, j + 1);
				unsigned int eb = j * W + i, et = (j + 1) * W + i, el = nh + j * (W + 1) + i, er = nh + j * (W + 1) + i + 1;
				int c = (f0 <= 0) | ((f1 <= 0) << 1) | ((f2 <= 0) << 2) | ((f3 <= 0) << 3);
				/* links run with the inside on the left */
				unsigned int lk[4][2];
				unsigned int nl = 0, q;
				switch (c)
				{
				case 1: lk[nl][0] = el; lk[nl][1] = eb; nl++; break;
				case 2: lk[nl][0] = eb; lk[nl][1] = er; nl++; break;
				case 3: lk[nl][0] = el; lk[nl][1] = er; nl++; break;
				case 4: lk[nl][0] = er; lk[nl][1] = et; nl++; break;
				case 6: lk[nl][0] = eb; lk[nl][1] = et; nl++; break;
				case 7: lk[nl][0] = el; lk[nl][1] = et; nl++; break;
				case 8: lk[nl][0] = et; lk[nl][1] = el; nl++; break;
				case 9: lk[nl][0] = et; lk[nl][1] = eb; nl++; break;
				case 11: lk[nl][0] = et; lk[nl][1] = er; nl++; break;
				case 12: lk[nl][0] = er; lk[nl][1] = el; nl++; break;
				case 13: lk[nl][0] = er; lk[nl][1] = eb; nl++; break;
				case 14: lk[nl][0] = eb; lk[nl][1] = el; nl++; break;
				case 5: case 10:
				{
					/* saddle: the centre decides whether the inside corners connect */
					int centre_in = 0.25 * (f0 + f1 + f2 + f3) <= 0;
					if (c == 5)
					{
						if (centre_in) { lk[nl][0] = el; lk[nl][1] = et; nl++; lk[nl][0] = er; lk[nl][1] = eb; nl++; }
						else { lk[nl][0] = el; lk[nl][1] = eb; nl++; lk[nl][0] = er; lk[nl][1] = et; nl++; }
					}
					else
					{
						if (centre_in) { lk[nl][0] = eb; lk[nl][1] = el; nl++; lk[nl][0] = et; lk[nl][1] = er; nl++; }
						else { lk[nl][0] = eb; lk[nl][1] = er; nl++; lk[nl][0] = et; lk[nl][1] = el; nl++; }
					}
					break;
				}
				default: break;
				}
				for (q = 0; q < nl; q++)
					nxt[lk[q][0]] = lk[q][1];
			}
		/* edge points by linear interpolation of F */
		for (j = 0; j <= Hc; j++)
			for (i = 0; i < W; i++)
			{
				double a = SK_F(i, j), b = SK_F(i + 1, j), t = (a - b) != 0 ? a / (a - b) : 0.5;
				e = j * W + i;
				ep[2 * e] = g.x0 + g.h * (i + t);
				ep[2 * e + 1] = g.y0 + g.h * j;
			}
		for (j = 0; j < Hc; j++)
			for (i = 0; i <= W; i++)
			{
				double a = SK_F(i, j), b = SK_F(i, j + 1), t = (a - b) != 0 ? a / (a - b) : 0.5;
				e = nh + j * (W + 1) + i;
				ep[2 * e] = g.x0 + g.h * i;
				ep[2 * e + 1] = g.y0 + g.h * (j + t);
			}
#undef SK_F
		/* walk the links into loops */
		for (e = 0; e < ne && s == QAWS_STATUS_OK; e++)
		{
			unsigned int q = e, m = 0, capm = 0;
			qaws_scalar* lp = NULL;
			if (nxt[e] == SK_NONE || seen[e]) continue;
			while (q != SK_NONE && !seen[q])
			{
				seen[q] = 1;
				if (m == capm)
				{
					unsigned int nc2 = capm ? 2 * capm : 64;
					qaws_scalar* gp = (qaws_scalar*)realloc(lp, sizeof(qaws_scalar) * 2 * nc2);
					if (!gp) { s = QAWS_STATUS_ALLOCATION_FAILURE; break; }
					lp = gp; capm = nc2;
				}
				lp[2 * m] = (qaws_scalar)ep[2 * q];
				lp[2 * m + 1] = (qaws_scalar)ep[2 * q + 1];
				m++;
				q = nxt[q];
			}
			if (s == QAWS_STATUS_OK && m >= 3)
			{
				qaws_curve* c = NULL;
				if (nloops == caploops)
				{
					unsigned int nc2 = caploops ? 2 * caploops : 16;
					qaws_curve** gl = (qaws_curve**)realloc(loops, sizeof(qaws_curve*) * nc2);
					if (!gl) { free(lp); s = QAWS_STATUS_ALLOCATION_FAILURE; break; }
					loops = gl; caploops = nc2;
				}
				s = qaws_curve_create_polyline_2d(lp, m, 1, &c);
				if (s == QAWS_STATUS_OK) loops[nloops++] = c;
			}
			free(lp);
		}
		if (s == QAWS_STATUS_OK)
		{
			qaws_path_2d* paths = (qaws_path_2d*)calloc(nloops + 1, sizeof(qaws_path_2d));
			qaws_curve const** views = (qaws_curve const**)calloc(nloops + 1, sizeof(qaws_curve*));
			if (!paths || !views)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
			else
			{
				for (i = 0; i < nloops; i++)
				{
					views[i] = loops[i];
					paths[i].curves = &views[i];
					paths[i].curve_count = 1;
					paths[i].closed = 1;
				}
				s = sk_union(paths, nloops, QAWS_FILL_NON_ZERO, out);
			}
			free(paths);
			free((void*)views);
		}
		for (i = 0; i < nloops; i++) qaws_curve_destroy(loops[i]);
		free(loops);
		free(ep); free(nxt); free(seen);
	}
done:
	free(g.f); free(pts); free(cp); free(ins);
	return s;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

qaws_status qaws_stack_section_2d(qaws_stack_2d const* st, qaws_scalar z, qaws_clip_result** out)
{
	unsigned int idx[4], n = 0, at = 0;
	double w[4];
	int r;
	if (!st || !out || (st->level_count && !st->levels))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	r = sk_weights(st, (double)z, idx, w, &n, &at);
	if (r == 0)
		return qaws_internal_clip_result_empty(out);
	if (r < 0)
		return sk_union(st->levels[at].paths, st->levels[at].path_count, st->fill_rule, out);
	if (st->blend != QAWS_STACK_BLEND_DISTANCE)
	{
		int applied = 0;
		qaws_status s = sk_match(st, idx, w, n, out, &applied);
		if (s != QAWS_STATUS_OK || applied)
			return s;
		if (st->blend == QAWS_STACK_BLEND_MATCH)
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
	}
	return sk_distance(st, idx, w, n, out);
}

struct qaws_stack_result
{
	unsigned int count;
	double* z;
	qaws_clip_result** levels;
};

typedef struct sk_job
{
	qaws_clip_type ct;
	qaws_stack_2d const *a, *b;
	double const* z;
	qaws_clip_result** levels;
} sk_job;

static qaws_status sk_level_job(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sk_job* j = (sk_job*)ctx;
	unsigned int i;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		qaws_clip_result *ra = NULL, *rb = NULL;
		qaws_status s = qaws_stack_section_2d(j->a, (qaws_scalar)j->z[i], &ra);
		if (s == QAWS_STATUS_OK) s = qaws_stack_section_2d(j->b, (qaws_scalar)j->z[i], &rb);
		if (s == QAWS_STATUS_OK)
		{
			unsigned int na = qaws_clip_result_get_path_count(ra), nb = qaws_clip_result_get_path_count(rb), k;
			qaws_path_2d* pa = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (na + 1));
			qaws_path_2d* pb = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (nb + 1));
			if (!pa || !pb)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
			else
			{
				for (k = 0; k < na; k++) qaws_clip_result_get_path(ra, k, &pa[k]);
				for (k = 0; k < nb; k++) qaws_clip_result_get_path(rb, k, &pb[k]);
				s = qaws_clip_boolean(j->ct, QAWS_FILL_NON_ZERO, pa, na, pb, nb, &j->levels[i]);
			}
			free(pa);
			free(pb);
		}
		qaws_clip_result_destroy(ra);
		qaws_clip_result_destroy(rb);
		if (s != QAWS_STATUS_OK)
			return s;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_stack_boolean_2d(qaws_clip_type ct, qaws_stack_2d const* a, qaws_stack_2d const* b,
	qaws_batch_executor const* executor, qaws_stack_result** out)
{
	qaws_stack_result* r;
	unsigned int i, ia = 0, ib = 0, n = 0;
	sk_job job;
	qaws_status s;
	if (!a || !b || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	r = (qaws_stack_result*)calloc(1, sizeof(qaws_stack_result));
	if (!r)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	r->z = (double*)malloc(sizeof(double) * (a->level_count + b->level_count + 1));
	r->levels = (qaws_clip_result**)calloc(a->level_count + b->level_count + 1, sizeof(qaws_clip_result*));
	if (!r->z || !r->levels)
	{
		qaws_stack_result_destroy(r);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* every height of either stack, sorted, once */
	while (ia < a->level_count || ib < b->level_count)
	{
		double za = ia < a->level_count ? (double)a->levels[ia].z : 1e300;
		double zb = ib < b->level_count ? (double)b->levels[ib].z : 1e300;
		double zz = za < zb ? za : zb;
		if (za == zz) ia++;
		if (zb == zz) ib++;
		if (!n || r->z[n - 1] != zz)
			r->z[n++] = zz;
	}
	r->count = n;
	job.ct = ct; job.a = a; job.b = b; job.z = r->z; job.levels = r->levels;
	s = qaws_internal_parallel(executor, n, 1, sk_level_job, &job);
	if (s != QAWS_STATUS_OK)
	{
		qaws_stack_result_destroy(r);
		return s;
	}
	for (i = 0; i < n; i++)
		if (!r->levels[i])
		{
			qaws_stack_result_destroy(r);
			return QAWS_STATUS_INTERNAL_ERROR;
		}
	*out = r;
	return QAWS_STATUS_OK;
}

void qaws_stack_result_destroy(qaws_stack_result* r)
{
	unsigned int i;
	if (!r) return;
	for (i = 0; r->levels && i < r->count; i++)
		qaws_clip_result_destroy(r->levels[i]);
	free(r->levels);
	free(r->z);
	free(r);
}

unsigned int qaws_stack_result_get_level_count(qaws_stack_result const* r) { return r ? r->count : 0; }
qaws_scalar qaws_stack_result_get_z(qaws_stack_result const* r, unsigned int i) { return r && i < r->count ? (qaws_scalar)r->z[i] : 0; }
qaws_clip_result const* qaws_stack_result_get_level(qaws_stack_result const* r, unsigned int i) { return r && i < r->count ? r->levels[i] : NULL; }
