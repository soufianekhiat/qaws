#include "qaws_exact_curve.h"
#include "qaws_exact_poly.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_broadphase.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * Certified closest points.
 *
 * For the query point p = P / Q (lattice units, Q a power of two) and a
 * span C = X / W, the squared distance |X / W - p|^2 is stationary where
 *   g(s) = sum_c (Q X_c - P_c W)(X_c' W - X_c W') = 0,
 * an integer polynomial of degree 3n - 2 in the local parameter. Its roots
 * in (0, 1), isolated exactly, and the span ends are the candidates; each
 * candidate's squared distance E / V, E = sum_c (Q X_c - P_c W)^2 and
 * V = (Q W)^2, is enclosed by the ratios of their Bernstein coefficients
 * restricted together to the candidate's interval (V > 0 there). The
 * nearest candidate is certified when its enclosure lies below every
 * other one. A span whose roots do not isolate (a double root, or g = 0:
 * every point at the same distance) enters as one candidate covering the
 * whole span.
 *
 * Spans are taken from one grid over their sound control boxes, nearest
 * box first, while a box can still beat the best enclosure.
 */

#define CL_DEPTH_START 40
#define CL_DEPTH_MAX QAWS_EXACT_ROOTS_MAX_DEPTH
#define CL_MAX_DEG (QAWS_EXACT_ROOTS_MAX_DEGREE)

typedef struct cl_cand
{
	unsigned int curve, span;
	qaws_exact_root r;          /* the candidate's interval (or exact point) */
	int whole;                  /* the whole span: roots not isolated */
	double lo, hi;              /* squared distance enclosure, lattice units */
} cl_cand;

typedef struct cl_span
{
	unsigned int curve, index;
	double lb;                  /* box distance from the query point (squared) */
} cl_span;

typedef struct cl_query
{
	unsigned int dim;
	qaws_exact_int P[3], Q;     /* p = P / Q */
	double pd[3];               /* p, rounded */
} cl_query;

/* polynomials of one span for one query point */
typedef struct cl_work
{
	qaws_exact_poly A[3], W, Xd, Wd, g, E, V, t1, t2;
	qaws_exact_int bg[CL_MAX_DEG + 1], bE[CL_MAX_DEG + 1], bV[CL_MAX_DEG + 1];
	qaws_exact_int rE[CL_MAX_DEG + 1], rV[CL_MAX_DEG + 1];
	qaws_exact_root roots[CL_MAX_DEG];
} cl_work;

#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

/* derivative in the power basis */
static qaws_status cl_deriv(qaws_exact_poly const* a, qaws_exact_poly* d)
{
	unsigned int i;
	qaws_status st;
	qaws_exact_poly_zero(d, a->deg ? a->deg - 1 : 0);
	for (i = 1; i <= a->deg; i++)
		TRY(qaws_exact_int_mul_i64(&d->c[i - 1], &a->c[i], (int64_t)i));
	return QAWS_STATUS_OK;
}

/* g, E and V of span sp for the query q (Bernstein, common degree for E and V) */
static qaws_status cl_build(qaws_exact_span const* sp, cl_query const* q, cl_work* w, unsigned int* ng, unsigned int* nev)
{
	unsigned int D = q->dim + 1, n = sp->degree, c, i;
	qaws_status st;
	TRY(qaws_exact_poly_from_bernstein(sp->h, n, D, q->dim, &w->W));
	TRY(cl_deriv(&w->W, &w->Wd));
	qaws_exact_poly_zero(&w->g, 0);
	qaws_exact_poly_zero(&w->E, 0);
	for (c = 0; c < q->dim; c++)
	{
		qaws_exact_poly X;
		/* A_c = Q X_c - P_c W */
		TRY(qaws_exact_poly_from_bernstein(sp->h, n, D, c, &X));
		TRY(qaws_exact_poly_scale(&w->A[c], &X, &q->Q));
		TRY(qaws_exact_poly_scale(&w->t1, &w->W, &q->P[c]));
		TRY(qaws_exact_poly_acc(&w->A[c], &w->t1, -1));
		/* B_c = X_c' W - X_c W' */
		TRY(cl_deriv(&X, &w->Xd));
		TRY(qaws_exact_poly_mul(&w->t1, &w->Xd, &w->W));
		TRY(qaws_exact_poly_mul(&w->t2, &X, &w->Wd));
		TRY(qaws_exact_poly_acc(&w->t1, &w->t2, -1));
		TRY(qaws_exact_poly_mul(&w->t2, &w->A[c], &w->t1));
		TRY(qaws_exact_poly_pad(&w->g, w->t2.deg > w->g.deg ? w->t2.deg : w->g.deg));
		TRY(qaws_exact_poly_pad(&w->t2, w->g.deg));
		TRY(qaws_exact_poly_acc(&w->g, &w->t2, 1));
		TRY(qaws_exact_poly_mul(&w->t2, &w->A[c], &w->A[c]));
		TRY(qaws_exact_poly_pad(&w->E, 2 * n));
		TRY(qaws_exact_poly_pad(&w->t2, 2 * n));
		TRY(qaws_exact_poly_acc(&w->E, &w->t2, 1));
	}
	/* V = (Q W)^2 */
	TRY(qaws_exact_poly_scale(&w->t1, &w->W, &q->Q));
	TRY(qaws_exact_poly_mul(&w->V, &w->t1, &w->t1));
	TRY(qaws_exact_poly_pad(&w->V, 2 * n));
	qaws_exact_poly_trim(&w->g);
	if (w->g.deg > CL_MAX_DEG || 2 * n > CL_MAX_DEG)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	*ng = w->g.deg;
	*nev = 2 * n;
	if (!qaws_exact_poly_is_zero(&w->g))
		TRY(qaws_exact_poly_to_bernstein(&w->g, w->bg));
	TRY(qaws_exact_poly_to_bernstein(&w->E, w->bE));
	TRY(qaws_exact_poly_to_bernstein(&w->V, w->bV));
	for (i = 0; i <= 2 * n; i++)
		if (qaws_exact_int_sign(&w->bV[i]) <= 0)
			return QAWS_STATUS_EXACT_UNSUPPORTED;   /* a weight not positive */
	return QAWS_STATUS_OK;
}

/* E / V on the candidate's interval (or at its exact point) */
static qaws_status cl_enclose(cl_work* w, unsigned int nev, qaws_exact_root const* r, int whole, double* lo, double* hi)
{
	unsigned int i, at_end;
	qaws_status st;
	*lo = HUGE_VAL;
	*hi = -HUGE_VAL;
	if (whole)
	{
		memcpy(w->rE, w->bE, sizeof(qaws_exact_int) * (nev + 1));
		memcpy(w->rV, w->bV, sizeof(qaws_exact_int) * (nev + 1));
	}
	else
	{
		/* an exact point: its value is the first coefficient of the cell it starts (or the last of the one it ends) */
		at_end = r->exact && r->index == ((uint64_t)1 << r->depth);
		TRY(qaws_exact_bernstein_restrict_pair(w->bE, w->bV, nev, r->index - at_end, r->depth, w->rE, w->rV));
		if (r->exact)
		{
			unsigned int e = at_end ? nev : 0;
			int ex;
			qaws_exact_ratio_enclose(&w->rE[e], &w->rV[e], lo, hi, &ex);
			return QAWS_STATUS_OK;
		}
	}
	for (i = 0; i <= nev; i++)
	{
		double l, h;
		int ex;
		qaws_exact_ratio_enclose(&w->rE[i], &w->rV[i], &l, &h, &ex);
		if (l < *lo) *lo = l;
		if (h > *hi) *hi = h;
	}
	return QAWS_STATUS_OK;
}

typedef struct cl_list
{
	cl_cand* c;
	unsigned int n, cap;
} cl_list;

static qaws_status cl_push(cl_list* L, cl_cand const* c)
{
	if (L->n == L->cap)
	{
		unsigned int cap = L->cap ? L->cap * 2 : 32;
		cl_cand* g = (cl_cand*)realloc(L->c, cap * sizeof(cl_cand));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		L->c = g;
		L->cap = cap;
	}
	L->c[L->n++] = *c;
	return QAWS_STATUS_OK;
}

static int cl_closed(qaws_exact_curve const* cv)
{
	unsigned int D = (unsigned int)cv->dimension + 1, i;
	qaws_exact_span const* f = &cv->spans[0];
	qaws_exact_span const* l = &cv->spans[cv->span_count - 1];
	qaws_exact_int x, y;
	for (i = 0; i < D - 1; i++)
	{
		if (qaws_exact_int_mul(&x, &f->h[i], &l->h[l->degree * D + D - 1]) != QAWS_STATUS_OK ||
		    qaws_exact_int_mul(&y, &l->h[l->degree * D + i], &f->h[D - 1]) != QAWS_STATUS_OK)
			return 0;
		if (qaws_exact_int_cmp(&x, &y) != 0)
			return 0;
	}
	return 1;
}

/* the candidates of one span: its roots and its ends (a knot once, the start
   of a closed curve once) */
static qaws_status cl_span_candidates(qaws_exact_curve const* cv, unsigned int ci, unsigned int si, int closed, cl_query const* q, cl_work* w,
	cl_list* L)
{
	qaws_exact_span const* sp = &cv->spans[si];
	unsigned int ng = 0, nev = 0, nr = 0, k;
	cl_cand c;
	qaws_status st, iso = QAWS_STATUS_CERTIFICATION_FAILED;
	TRY(cl_build(sp, q, w, &ng, &nev));
	memset(&c, 0, sizeof(c));
	c.curve = ci;
	c.span = si;
	if (!qaws_exact_poly_is_zero(&w->g) && ng > 0)
		iso = qaws_exact_bernstein_isolate(w->bg, ng, w->roots, CL_MAX_DEG, &nr);
	else if (!qaws_exact_poly_is_zero(&w->g))
		iso = QAWS_STATUS_OK;   /* a non-zero constant: no stationary point */
	if (iso == QAWS_STATUS_ALLOCATION_FAILURE)
		return iso;
	if (iso != QAWS_STATUS_OK)
	{
		/* not isolated: the whole span is one candidate */
		c.whole = 1;
		TRY(cl_enclose(w, nev, NULL, 1, &c.lo, &c.hi));
		return cl_push(L, &c);
	}
	for (k = 0; k < nr; k++)
	{
		qaws_exact_root r = w->roots[k];
		if (r.exact && (r.index == 0 || r.index == ((uint64_t)1 << r.depth)))
			continue;   /* a span end, below */
		if (!r.exact)
		{
			st = qaws_exact_bernstein_refine(w->bg, ng, &r, CL_DEPTH_START);
			if (st == QAWS_STATUS_EXACT_RANGE_EXCEEDED)
				r = w->roots[k];
			else if (st != QAWS_STATUS_OK)
				return st;
		}
		c.r = r;
		TRY(cl_enclose(w, nev, &c.r, 0, &c.lo, &c.hi));
		TRY(cl_push(L, &c));
	}
	/* ends: s = 1 always, s = 0 for the first span of an open curve */
	for (k = 0; k < 2; k++)
	{
		if (k == 0 && (si > 0 || closed))
			continue;
		c.r.exact = 1;
		c.r.depth = 0;
		c.r.index = k;
		TRY(cl_enclose(w, nev, &c.r, 0, &c.lo, &c.hi));
		TRY(cl_push(L, &c));
	}
	return QAWS_STATUS_OK;
}

static int cl_cmp_span(void const* a, void const* b)
{
	double x = ((cl_span const*)a)->lb, y = ((cl_span const*)b)->lb;
	return x < y ? -1 : (x > y ? 1 : 0);
}

typedef struct cl_ring
{
	qaws_exact_batch_desc const* desc;
	cl_span* all;               /* every span, with its box */
	double (*box)[2][3];
	unsigned int* stamp;
	unsigned int query, dim;
	double const* p;
	cl_span* got;
	unsigned int ngot, capgot;
	double best_ub;             /* squared: the nearest farthest box corner seen */
	int failed;
} cl_ring;

static double cl_box_dist2(double const* lo, double const* hi, double const* p, unsigned int dim)
{
	double d = 0;
	unsigned int k;
	for (k = 0; k < dim; k++)
	{
		double e = 0;
		if (!(lo[k] <= hi[k]))
			return 0;   /* unbounded */
		if (p[k] < lo[k]) e = lo[k] - p[k];
		else if (p[k] > hi[k]) e = p[k] - hi[k];
		d += e * e;
	}
	return d;
}

static void cl_ring_visit(void* user, unsigned int item)
{
	cl_ring* x = (cl_ring*)user;
	if (x->stamp[item] == x->query)
		return;
	x->stamp[item] = x->query;
	if (x->ngot == x->capgot)
	{
		unsigned int cap = x->capgot ? x->capgot * 2 : 64;
		cl_span* g = (cl_span*)realloc(x->got, cap * sizeof(cl_span));
		if (!g)
		{
			x->failed = 1;
			return;
		}
		x->got = g;
		x->capgot = cap;
	}
	x->got[x->ngot] = x->all[item];
	x->got[x->ngot].lb = cl_box_dist2(x->box[item][0], x->box[item][1], x->p, x->dim);
	x->ngot++;
	/* a box holds curve points: its farthest corner bounds the nearest distance */
	if (x->box[item][0][0] <= x->box[item][1][0])
	{
		double fd = 0;
		unsigned int c;
		for (c = 0; c < x->dim; c++)
		{
			double a = fabs(x->p[c] - x->box[item][0][c]), b = fabs(x->p[c] - x->box[item][1][c]);
			fd += (a > b ? a : b) * (a > b ? a : b);
		}
		if (fd < x->best_ub)
			x->best_ub = fd;
	}
}

/* the query point as P / Q on the lattice: p 2^-space_exp2 */
static qaws_status cl_query_make(double const* p, unsigned int dim, int space_exp2, cl_query* q)
{
	int64_t m[3];
	int e[3], k = 0, shift[3];
	unsigned int c;
	qaws_status st;
	q->dim = dim;
	for (c = 0; c < dim; c++)
	{
		if (!qaws_exact_split_double(p[c], &m[c], &e[c]))
			return QAWS_STATUS_INVALID_ARGUMENT;
		e[c] -= space_exp2;
		if (m[c] != 0 && -e[c] > k)
			k = -e[c];
		q->pd[c] = ldexp(p[c], -space_exp2);
	}
	/* Q = 2^k, P_c = m_c 2^(e_c + k) */
	for (c = 0; c < dim; c++)
	{
		shift[c] = m[c] ? e[c] + k : 0;
		qaws_exact_int_from_i64(&q->P[c], m[c]);
		if (shift[c] > 0)
			TRY(qaws_exact_int_shl(&q->P[c], &q->P[c], (unsigned int)shift[c]));
	}
	qaws_exact_int_from_i64(&q->Q, 1);
	if (k > 0)
		TRY(qaws_exact_int_shl(&q->Q, &q->Q, (unsigned int)k));
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_curve_batch_closest(qaws_exact_batch_desc const* desc, double const* points, unsigned int point_count,
	qaws_exact_closest_point* out_points, qaws_exact_batch_stats* out_stats)
{
	qaws_exact_batch_stats stats;
	cl_ring R;
	cl_list L;
	cl_work* w = NULL;
	qaws_bp_box* boxes = NULL;
	qaws_bp_grid* grid = NULL;
	int* closed = NULL;
	unsigned int i, k, nspan = 0, dim;
	int space;
	double scale;
	qaws_status st = QAWS_STATUS_OK;
	memset(&stats, 0, sizeof(stats));
	memset(&R, 0, sizeof(R));
	memset(&L, 0, sizeof(L));
	if (out_stats)
		*out_stats = stats;
	if (!desc || (desc->curve_count && !desc->curves) || (point_count && (!points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < point_count; i++)
	{
		memset(&out_points[i], 0, sizeof(out_points[i]));
		out_points[i].curve = QAWS_EXACT_CLOSEST_NONE;
	}
	if (!desc->curve_count || !point_count)
		return QAWS_STATUS_OK;
	for (i = 0; i < desc->curve_count; i++)
	{
		qaws_exact_curve const* c = desc->curves[i];
		if (!c || (c->dimension != 2 && c->dimension != 3) || c->dimension != desc->curves[0]->dimension)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (c->space_exp2 != desc->curves[0]->space_exp2)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		nspan += c->span_count;
	}
	dim = (unsigned int)desc->curves[0]->dimension;
	space = desc->curves[0]->space_exp2;
	scale = ldexp(1.0, space);
	R.desc = desc;
	R.dim = dim;
	R.all = (cl_span*)malloc(nspan * sizeof(cl_span));
	R.box = (double (*)[2][3])malloc(nspan * sizeof(double) * 6);
	R.stamp = (unsigned int*)malloc(nspan * sizeof(unsigned int));
	boxes = (qaws_bp_box*)malloc(nspan * sizeof(qaws_bp_box));
	closed = (int*)malloc(desc->curve_count * sizeof(int));
	w = (cl_work*)malloc(sizeof(cl_work));
	if (!R.all || !R.box || !R.stamp || !boxes || !closed || !w)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	nspan = 0;
	for (i = 0; i < desc->curve_count; i++)
	{
		closed[i] = cl_closed(desc->curves[i]);
		for (k = 0; k < desc->curves[i]->span_count; k++, nspan++)
		{
			unsigned int c;
			R.all[nspan].curve = i;
			R.all[nspan].index = k;
			R.stamp[nspan] = ~0u;
			qaws_exact_span_box(&desc->curves[i]->spans[k], dim, R.box[nspan][0], R.box[nspan][1]);
			for (c = 0; c < 3; c++)
			{
				boxes[nspan].lo[c] = R.box[nspan][0][c];
				boxes[nspan].hi[c] = R.box[nspan][1][c];
			}
		}
	}
	stats.span_count = nspan;
	st = qaws_internal_grid_create(boxes, nspan, dim, &grid);
	for (i = 0; i < point_count && st == QAWS_STATUS_OK; i++)
	{
		cl_query q;
		double best_hi = HUGE_VAL;
		unsigned int r, s, win = 0, j;
		int sure;
		st = cl_query_make(points + i * dim, dim, space, &q);
		if (st != QAWS_STATUS_OK)
			break;
		/* spans from the rings, until a ring cannot beat the best box bound */
		R.query = i;
		R.p = q.pd;
		R.ngot = 0;
		R.best_ub = HUGE_VAL;
		for (r = 0;; r++)
		{
			double b = qaws_internal_grid_ring(grid, q.pd, r, cl_ring_visit, &R);
			if (R.failed)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			if (b == HUGE_VAL || b * b > R.best_ub)
				break;
		}
		if (st != QAWS_STATUS_OK)
			break;
		qsort(R.got, R.ngot, sizeof(cl_span), cl_cmp_span);
		/* the spans nearest box first, while a box can beat the best enclosure */
		L.n = 0;
		for (s = 0; s < R.ngot && st == QAWS_STATUS_OK; s++)
		{
			unsigned int c0 = L.n;
			if (R.got[s].lb > best_hi * (1 + 1e-12))
				break;
			stats.candidate_count++;
			st = cl_span_candidates(desc->curves[R.got[s].curve], R.got[s].curve, R.got[s].index, closed[R.got[s].curve], &q, w, &L);
			for (j = c0; j < L.n && st == QAWS_STATUS_OK; j++)
				if (L.c[j].hi < best_hi)
					best_hi = L.c[j].hi;
		}
		if (st != QAWS_STATUS_OK)
			break;
		if (!L.n)
			continue;
		/* the winner: lowest upper bound; certified when below every other lower bound,
		   refining the contenders' intervals while they overlap */
		for (;;)
		{
			int refined = 0;
			win = 0;
			for (j = 1; j < L.n; j++)
				if (L.c[j].hi < L.c[win].hi)
					win = j;
			sure = !L.c[win].whole;
			for (j = 0; j < L.n && sure; j++)
				if (j != win && L.c[j].lo <= L.c[win].hi)
					sure = 0;
			if (sure)
				break;
			/* shrink every overlapping isolated root interval that can still shrink */
			for (j = 0; j < L.n && st == QAWS_STATUS_OK; j++)
			{
				cl_cand* c = &L.c[j];
				unsigned int ng = 0, nev = 0;
				if (c->whole || c->r.exact || c->r.depth >= CL_DEPTH_MAX || (j != win && c->lo > L.c[win].hi))
					continue;
				st = cl_build(&desc->curves[c->curve]->spans[c->span], &q, w, &ng, &nev);
				if (st == QAWS_STATUS_OK)
					st = qaws_exact_bernstein_refine(w->bg, ng, &c->r, c->r.depth + 8 > CL_DEPTH_MAX ? CL_DEPTH_MAX : c->r.depth + 8);
				if (st == QAWS_STATUS_EXACT_RANGE_EXCEEDED)
				{
					st = QAWS_STATUS_OK;
					c->r.depth = CL_DEPTH_MAX;
					continue;
				}
				if (st == QAWS_STATUS_OK)
					st = cl_enclose(w, nev, &c->r, 0, &c->lo, &c->hi);
				refined = 1;
			}
			if (st != QAWS_STATUS_OK || !refined)
				break;
		}
		if (st != QAWS_STATUS_OK)
			break;
		{
			cl_cand const* c = &L.c[win];
			qaws_exact_closest_point* o = &out_points[i];
			qaws_exact_span const* sp = &desc->curves[c->curve]->spans[c->span];
			int shift = desc->curves[c->curve]->param_shift, ex;
			o->curve = c->curve;
			o->certified = sure;
			if (c->whole)
			{
				o->t_lo = ldexp((double)sp->a, -shift);
				o->t_hi = ldexp((double)sp->b, -shift);
			}
			else
			{
				st = qaws_exact_span_param_to_double(sp, shift, c->r.index, c->r.depth, &o->t_lo, &ex);
				if (st == QAWS_STATUS_OK && !ex) o->t_lo = nextafter(o->t_lo, -HUGE_VAL);
				if (st == QAWS_STATUS_OK)
					st = qaws_exact_span_param_to_double(sp, shift, c->r.index + (c->r.exact ? 0 : 1), c->r.depth, &o->t_hi, &ex);
				if (st == QAWS_STATUS_OK && !ex) o->t_hi = nextafter(o->t_hi, HUGE_VAL);
			}
			/* world units, rounded outward */
			o->distance_lo = nextafter(sqrt(c->lo > 0 ? c->lo : 0), -HUGE_VAL) * scale;
			o->distance_hi = nextafter(sqrt(c->hi), HUGE_VAL) * scale;
			if (o->distance_lo < 0)
				o->distance_lo = 0;
			stats.hit_count++;
			if (!sure)
				stats.uncertified_count++;
		}
	}
done:
	if (out_stats)
		*out_stats = stats;
	qaws_internal_grid_destroy(grid);
	free(R.all);
	free(R.box);
	free(R.stamp);
	free(R.got);
	free(L.c);
	free(boxes);
	free(closed);
	free(w);
	return st;
}
