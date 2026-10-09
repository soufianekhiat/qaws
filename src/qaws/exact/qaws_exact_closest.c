#include "qaws_exact_surface.h"
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

/* the candidates of one span: its roots and the ends asked for (a knot once,
   the start of a closed curve once) */
static qaws_status cl_span_candidates(qaws_exact_span const* sp, unsigned int ci, unsigned int si, int with_start, int with_end, cl_query const* q,
	cl_work* w, cl_list* L)
{
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
		if ((k == 0 && !with_start) || (k == 1 && !with_end))
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
			st = cl_span_candidates(&desc->curves[R.got[s].curve]->spans[R.got[s].index], R.got[s].curve, R.got[s].index,
				R.got[s].index == 0 && !closed[R.got[s].curve], 1, &q, w, &L);
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

/* ------------------------------------------------------------------ */
/*  Surfaces                                                           */
/* ------------------------------------------------------------------ */

/*
 * On a patch S = X / W (degrees p, q) the squared distance to p = P / Q is
 * stationary where
 *   G_u = sum_c (Q X_c - P_c W)(X_c,u W - X_c W_u) = 0,
 *   G_v = sum_c (Q X_c - P_c W)(X_c,v W - X_c W_v) = 0,
 * integer tensors of degree (3p, 3q) at most. With the third equation
 * 3w - 1 = 0 (a root never on a dyadic cut) the system goes to the
 * certified 3D solver; its root boxes inside the patch, the edges as exact
 * curves (each edge once over the surface: the first row and column own
 * their lower edges) and their candidates, are the candidates of the patch.
 * Distances are enclosed by E / V restricted to the candidate's box.
 */

#define SC_MAX 10               /* 3 * 3 + 1 coefficients per direction (bicubic) */

typedef struct t2
{
	unsigned int du, dv;        /* power-basis degrees */
	qaws_exact_int c[SC_MAX * SC_MAX];
} t2;

#define T2(t, i, j) ((t)->c[(i) * SC_MAX + (j)])

static void t2_zero(t2* t, unsigned int du, unsigned int dv)
{
	unsigned int i, j;
	t->du = du;
	t->dv = dv;
	for (i = 0; i <= du; i++)
		for (j = 0; j <= dv; j++)
			qaws_exact_int_zero(&T2(t, i, j));
}

/* component c of a homogeneous patch net (u index major, 4 per point) to the power basis */
static qaws_status t2_from_patch(qaws_exact_int const* h, unsigned int p, unsigned int q, unsigned int comp, t2* out)
{
	qaws_exact_poly a;
	qaws_exact_int col[SC_MAX];
	unsigned int i, j;
	qaws_status st;
	t2 tmp;
	t2_zero(&tmp, p, q);
	/* along u for every v index j: the net seen with stride (q + 1) 4 */
	for (j = 0; j <= q; j++)
	{
		TRY(qaws_exact_poly_from_bernstein(h + j * 4, p, (q + 1) * 4, comp, &a));
		for (i = 0; i <= p && i <= a.deg; i++)
			T2(&tmp, i, j) = a.c[i];
	}
	/* along v for every u power i */
	t2_zero(out, p, q);
	for (i = 0; i <= p; i++)
	{
		for (j = 0; j <= q; j++)
			col[j] = T2(&tmp, i, j);
		TRY(qaws_exact_poly_from_bernstein(col, q, 1, 0, &a));
		for (j = 0; j <= q && j <= a.deg; j++)
			T2(out, i, j) = a.c[j];
	}
	return QAWS_STATUS_OK;
}

static qaws_status t2_mul(t2* r, t2 const* a, t2 const* b)
{
	unsigned int i, j, k, l;
	qaws_exact_int t;
	qaws_status st;
	if (a->du + b->du >= SC_MAX || a->dv + b->dv >= SC_MAX)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	t2_zero(r, a->du + b->du, a->dv + b->dv);
	for (i = 0; i <= a->du; i++)
		for (j = 0; j <= a->dv; j++)
		{
			if (qaws_exact_int_is_zero(&T2(a, i, j)))
				continue;
			for (k = 0; k <= b->du; k++)
				for (l = 0; l <= b->dv; l++)
				{
					TRY(qaws_exact_int_mul(&t, &T2(a, i, j), &T2(b, k, l)));
					TRY(qaws_exact_int_add(&T2(r, i + k, j + l), &T2(r, i + k, j + l), &t));
				}
		}
	return QAWS_STATUS_OK;
}

/* r += sign a (r grows to cover a) */
static qaws_status t2_acc(t2* r, t2 const* a, int sign)
{
	unsigned int i, j, du = a->du > r->du ? a->du : r->du, dv = a->dv > r->dv ? a->dv : r->dv;
	qaws_status st;
	for (i = 0; i <= du; i++)
		for (j = 0; j <= dv; j++)
			if (i > r->du || j > r->dv)
				qaws_exact_int_zero(&T2(r, i, j));
	r->du = du;
	r->dv = dv;
	for (i = 0; i <= a->du; i++)
		for (j = 0; j <= a->dv; j++)
			TRY(sign > 0 ? qaws_exact_int_add(&T2(r, i, j), &T2(r, i, j), &T2(a, i, j)) : qaws_exact_int_sub(&T2(r, i, j), &T2(r, i, j), &T2(a, i, j)));
	return QAWS_STATUS_OK;
}

static qaws_status t2_scale(t2* r, t2 const* a, qaws_exact_int const* k)
{
	unsigned int i, j;
	qaws_status st;
	t2_zero(r, a->du, a->dv);
	for (i = 0; i <= a->du; i++)
		for (j = 0; j <= a->dv; j++)
			TRY(qaws_exact_int_mul(&T2(r, i, j), &T2(a, i, j), k));
	return QAWS_STATUS_OK;
}

/* d/du (dir 0) or d/dv (dir 1) */
static qaws_status t2_deriv(t2* r, t2 const* a, int dir)
{
	unsigned int i, j;
	qaws_status st;
	t2_zero(r, dir == 0 && a->du ? a->du - 1 : a->du, dir == 1 && a->dv ? a->dv - 1 : a->dv);
	for (i = 0; i <= a->du; i++)
		for (j = 0; j <= a->dv; j++)
		{
			unsigned int e = dir == 0 ? i : j;
			if (e == 0)
				continue;
			TRY(qaws_exact_int_mul_i64(&T2(r, dir == 0 ? i - 1 : i, dir == 1 ? j - 1 : j), &T2(a, i, j), (int64_t)e));
		}
	return QAWS_STATUS_OK;
}

/* Bernstein coefficients of degree (Du, Dv) >= (du, dv) into b[a (Dv + 1) + b]
   (one positive factor for the whole tensor) */
static qaws_status t2_to_bern(t2 const* t, unsigned int Du, unsigned int Dv, qaws_exact_int* b)
{
	qaws_exact_poly a;
	qaws_exact_int tmp[SC_MAX * SC_MAX], col[SC_MAX];
	unsigned int i, j;
	qaws_status st;
	for (j = 0; j <= Dv; j++)
	{
		qaws_exact_poly_zero(&a, Du);
		for (i = 0; i <= t->du && j <= t->dv; i++)
			a.c[i] = T2(t, i, j);
		TRY(qaws_exact_poly_to_bernstein(&a, col));
		for (i = 0; i <= Du; i++)
			tmp[i * SC_MAX + j] = col[i];
	}
	for (i = 0; i <= Du; i++)
	{
		qaws_exact_poly_zero(&a, Dv);
		for (j = 0; j <= Dv; j++)
			a.c[j] = tmp[i * SC_MAX + j];
		TRY(qaws_exact_poly_to_bernstein(&a, col));
		for (j = 0; j <= Dv; j++)
			b[i * (Dv + 1) + j] = col[j];
	}
	return QAWS_STATUS_OK;
}

/* Bernstein coefficients b[0..n] (stride s) restricted to [lo, hi] / 2^d, in
   place, every coefficient times one positive factor: two integer
   de Casteljau splits */
static qaws_status sc_restrict1(qaws_exact_int* b, unsigned int n, unsigned int s, uint64_t lo, uint64_t hi, int d)
{
	qaws_exact_int lv[SC_MAX], r[SC_MAX], D, A0, A1, B0, B1, den, t, u;
	unsigned int k, i;
	qaws_status st;
	qaws_exact_int_from_i64(&D, 1);
	TRY(qaws_exact_int_shl(&D, &D, (unsigned int)d));
	/* right part at lo / 2^d, weights (2^d - lo, lo): r_j = b_j^(n - j), scaled by (2^d)^(n - j) */
	qaws_exact_int_from_i64(&A1, (int64_t)lo);
	TRY(qaws_exact_int_sub(&A0, &D, &A1));
	for (i = 0; i <= n; i++)
		lv[i] = b[i * s];
	r[n] = lv[n];
	for (k = 1; k <= n; k++)
	{
		for (i = 0; i + k <= n; i++)
		{
			TRY(qaws_exact_int_mul(&t, &A0, &lv[i]));
			TRY(qaws_exact_int_mul(&u, &A1, &lv[i + 1]));
			TRY(qaws_exact_int_add(&lv[i], &t, &u));
		}
		r[n - k] = lv[n - k];
	}
	for (i = 0; i <= n; i++)
		TRY(qaws_exact_int_shl(&r[i], &r[i], (unsigned int)(d * (int)i)));
	/* left part at (hi - lo) / (2^d - lo), weights (2^d - hi, hi - lo):
	   l_k = b_0^(k), scaled by (2^d - lo)^k */
	qaws_exact_int_from_i64(&B1, (int64_t)(hi - lo));
	qaws_exact_int_from_i64(&t, (int64_t)hi);
	TRY(qaws_exact_int_sub(&B0, &D, &t));
	den = A0;
	for (i = 0; i <= n; i++)
		lv[i] = r[i];
	b[0] = lv[0];
	for (k = 1; k <= n; k++)
	{
		for (i = 0; i + k <= n; i++)
		{
			TRY(qaws_exact_int_mul(&t, &B0, &lv[i]));
			TRY(qaws_exact_int_mul(&u, &B1, &lv[i + 1]));
			TRY(qaws_exact_int_add(&lv[i], &t, &u));
		}
		b[k * s] = lv[0];
	}
	/* one factor: l_k times (2^d - lo)^(n - k) */
	for (k = 0; k < n; k++)
		for (i = k; i < n; i++)
			TRY(qaws_exact_int_mul(&b[k * s], &b[k * s], &den));
	return QAWS_STATUS_OK;
}

/* the tensor b (degree Du, Dv) restricted to the box, in place */
static qaws_status sc_restrict2(qaws_exact_int* b, unsigned int Du, unsigned int Dv, uint64_t const* lo, uint64_t const* hi, int const* dep)
{
	unsigned int i;
	qaws_status st;
	for (i = 0; i <= Dv; i++)
		TRY(sc_restrict1(b + i, Du, Dv + 1, lo[0], hi[0], dep[0]));
	for (i = 0; i <= Du; i++)
		TRY(sc_restrict1(b + i * (Dv + 1), Dv, 1, lo[1], hi[1], dep[1]));
	return QAWS_STATUS_OK;
}

typedef struct sc_cand2
{
	unsigned int surface;
	double ulo, uhi, vlo, vhi;
	double lo, hi;              /* squared distance, lattice units */
	int whole;
} sc_cand2;

typedef struct sc_list
{
	sc_cand2* c;
	unsigned int n, cap;
} sc_list;

static qaws_status sc_push(sc_list* L, sc_cand2 const* c)
{
	if (L->n == L->cap)
	{
		unsigned int cap = L->cap ? L->cap * 2 : 32;
		sc_cand2* g = (sc_cand2*)realloc(L->c, cap * sizeof(sc_cand2));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		L->c = g;
		L->cap = cap;
	}
	L->c[L->n++] = *c;
	return QAWS_STATUS_OK;
}

typedef struct sc_work
{
	t2 X[4], A, B, Gu, Gv, E, V, d1, d2, t;
	qaws_exact_int F[3 * 2 * SC_MAX * SC_MAX];
	qaws_exact_int bE[SC_MAX * SC_MAX * 4], bV[SC_MAX * SC_MAX * 4], rE[SC_MAX * SC_MAX * 4], rV[SC_MAX * SC_MAX * 4];
	qaws_exact_int gbu[SC_MAX * SC_MAX], gbv[SC_MAX * SC_MAX], lgu[SC_MAX * SC_MAX], lgv[SC_MAX * SC_MAX];
	qaws_exact_box3 boxes[64];
	qaws_exact_int edge[(SC_MAX) * 4];
	cl_work cw;
	cl_list cl;
} sc_work;

/* E / V enclosure over the box (lo, hi, dep per direction) of the patch */
static qaws_status sc_enclose(sc_work* w, unsigned int De, uint64_t const* lo, uint64_t const* hi, int const* dep, double* elo, double* ehi)
{
	unsigned int i, n = (De + 1) * (De + 1);
	qaws_status st;
	memcpy(w->rE, w->bE, sizeof(qaws_exact_int) * n);
	memcpy(w->rV, w->bV, sizeof(qaws_exact_int) * n);
	if (lo)
	{
		TRY(sc_restrict2(w->rE, De, De, lo, hi, dep));
		TRY(sc_restrict2(w->rV, De, De, lo, hi, dep));
	}
	*elo = HUGE_VAL;
	*ehi = -HUGE_VAL;
	for (i = 0; i < n; i++)
	{
		double l, h;
		int ex;
		if (qaws_exact_int_sign(&w->rV[i]) <= 0)
			return QAWS_STATUS_EXACT_UNSUPPORTED;
		qaws_exact_ratio_enclose(&w->rE[i], &w->rV[i], &l, &h, &ex);
		if (l < *elo) *elo = l;
		if (h > *ehi) *ehi = h;
	}
	return QAWS_STATUS_OK;
}

static double sc_param(int64_t a, int64_t b, int shift, uint64_t x, int dep, int up)
{
	/* t = (a + (b - a) x / 2^dep) 2^-shift, rounded outward */
	double t = ldexp((double)a, -shift) + ldexp((double)(b - a), -shift) * ldexp((double)x, -dep);
	return up ? nextafter(t, HUGE_VAL) : nextafter(t, -HUGE_VAL);
}

#define SC_LEAF_DEPTH 4      /* interior localisation: 16 x 16 leaves at most */
#define SC_MAX_LEAVES 32
#define SC_LEAF_BOXES 400

/* a lower bound of the squared distance over patch (iu, iv): E / V over its Bernstein coefficients */
static qaws_status sc_patch_bound(qaws_exact_surface const* S, unsigned int iu, unsigned int iv, cl_query const* q, sc_work* w, double* elo)
{
	qaws_exact_int const* h = S->patch[iu * S->nv + iv];
	unsigned int p = S->p, qd = S->q, c, De = 2 * (p > qd ? p : qd);
	double ehi;
	qaws_status st;
	if (De + 1 > SC_MAX)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	for (c = 0; c < 4; c++)
		TRY(t2_from_patch(h, p, qd, c, &w->X[c]));
	t2_zero(&w->E, 0, 0);
	for (c = 0; c < 3; c++)
	{
		TRY(t2_scale(&w->A, &w->X[c], &q->Q));
		TRY(t2_scale(&w->t, &w->X[3], &q->P[c]));
		TRY(t2_acc(&w->A, &w->t, -1));
		TRY(t2_mul(&w->t, &w->A, &w->A));
		TRY(t2_acc(&w->E, &w->t, 1));
	}
	TRY(t2_scale(&w->t, &w->X[3], &q->Q));
	TRY(t2_mul(&w->V, &w->t, &w->t));
	TRY(t2_to_bern(&w->E, De, De, w->bE));
	TRY(t2_to_bern(&w->V, De, De, w->bV));
	return sc_enclose(w, De, NULL, NULL, NULL, elo, &ehi);
}

/* the candidates of patch (iu, iv) of surface si */
static qaws_status sc_patch_candidates(qaws_exact_surface const* S, unsigned int si, unsigned int iu, unsigned int iv, cl_query const* q, double best_hi,
	sc_work* w, sc_list* L)
{
	qaws_exact_int const* h = S->patch[iu * S->nv + iv];
	unsigned int p = S->p, qd = S->q, c, nb = 0, k, De = 2 * (p > qd ? p : qd), Dg = 3 * (p > qd ? p : qd), a, b;
	unsigned int n3[3];
	sc_cand2 cand;
	qaws_status st, s3;
	double u0 = ldexp((double)S->ub[iu], -S->u_shift), u1 = ldexp((double)S->ub[iu + 1], -S->u_shift);
	double v0 = ldexp((double)S->vb[iv], -S->v_shift), v1 = ldexp((double)S->vb[iv + 1], -S->v_shift);
	if (Dg + 1 > SC_MAX)
		return QAWS_STATUS_EXACT_UNSUPPORTED;
	/* X_c and W in the power basis */
	for (c = 0; c < 4; c++)
		TRY(t2_from_patch(h, p, qd, c, &w->X[c]));
	t2_zero(&w->Gu, 0, 0);
	t2_zero(&w->Gv, 0, 0);
	t2_zero(&w->E, 0, 0);
	for (c = 0; c < 3; c++)
	{
		int dir;
		/* A = Q X_c - P_c W */
		TRY(t2_scale(&w->A, &w->X[c], &q->Q));
		TRY(t2_scale(&w->t, &w->X[3], &q->P[c]));
		TRY(t2_acc(&w->A, &w->t, -1));
		TRY(t2_mul(&w->t, &w->A, &w->A));
		TRY(t2_acc(&w->E, &w->t, 1));
		for (dir = 0; dir < 2; dir++)
		{
			/* B = X_c,d W - X_c W_d */
			TRY(t2_deriv(&w->d1, &w->X[c], dir));
			TRY(t2_mul(&w->B, &w->d1, &w->X[3]));
			TRY(t2_deriv(&w->d2, &w->X[3], dir));
			TRY(t2_mul(&w->t, &w->X[c], &w->d2));
			TRY(t2_acc(&w->B, &w->t, -1));
			TRY(t2_mul(&w->t, &w->A, &w->B));
			TRY(t2_acc(dir == 0 ? &w->Gu : &w->Gv, &w->t, 1));
		}
	}
	/* V = (Q W)^2 */
	TRY(t2_scale(&w->t, &w->X[3], &q->Q));
	TRY(t2_mul(&w->V, &w->t, &w->t));
	TRY(t2_to_bern(&w->E, De, De, w->bE));
	TRY(t2_to_bern(&w->V, De, De, w->bV));
	/* interior: branch and bound on E / V localises the boxes that can beat the
	   best distance (the box corners, exact surface points, tighten it); only
	   those leaves are solved, G_u = G_v = 0 with 3 w - 1 = 0 restricted to
	   each, degrees (Dg, Dg, 1) */
	TRY(t2_to_bern(&w->Gu, Dg, Dg, w->gbu));
	TRY(t2_to_bern(&w->Gv, Dg, Dg, w->gbv));
	memset(&cand, 0, sizeof(cand));
	cand.surface = si;
	{
		typedef struct lb3 { uint64_t lo[2]; int d; } lb3;
		lb3 stack[4 * SC_LEAF_DEPTH + 4], leaves[SC_MAX_LEAVES];
		unsigned int ns = 1, nl = 0, gs = (Dg + 1) * (Dg + 1), size = gs * 2, l, visited = 0;
		int whole = 0;
		double best = best_hi;
		stack[0].lo[0] = stack[0].lo[1] = 0;
		stack[0].d = 0;
		while (ns > 0 && !whole)
		{
			lb3 bx = stack[--ns];
			uint64_t lo[3], hi[3];
			int dep[3], q4;
			double elo, ehi;
			unsigned int corner[4], c4;
			lo[0] = bx.lo[0]; hi[0] = bx.lo[0] + 1; lo[1] = bx.lo[1]; hi[1] = bx.lo[1] + 1;
			dep[0] = dep[1] = bx.d;
			lo[2] = 0; hi[2] = 1; dep[2] = 0;
			TRY(sc_enclose(w, De, bx.d ? lo : NULL, hi, dep, &elo, &ehi));
			/* the corners are surface points: exact distances bound the best */
			corner[0] = 0; corner[1] = De; corner[2] = De * (De + 1); corner[3] = De * (De + 1) + De;
			for (c4 = 0; c4 < 4; c4++)
			{
				double cl, ch;
				int ex;
				qaws_exact_ratio_enclose(&w->rE[corner[c4]], &w->rV[corner[c4]], &cl, &ch, &ex);
				if (ch < best)
					best = ch;
			}
			if (elo > best || ++visited > SC_LEAF_BOXES)
			{
				if (visited > SC_LEAF_BOXES)
					whole = 1;
				continue;
			}
			if (bx.d >= SC_LEAF_DEPTH)
			{
				if (nl >= SC_MAX_LEAVES)
					whole = 1;
				else
					leaves[nl++] = bx;
				continue;
			}
			for (q4 = 0; q4 < 4; q4++)
			{
				lb3 ch;
				ch.d = bx.d + 1;
				ch.lo[0] = 2 * bx.lo[0] + (uint64_t)(q4 & 1);
				ch.lo[1] = 2 * bx.lo[1] + (uint64_t)(q4 >> 1);
				stack[ns++] = ch;
			}
		}
		/* leaves ending before the depth (kept above as nearer) cannot occur: every
		   box nearer than best reaches the depth or is a leaf */
		for (l = 0; l < nl && !whole; l++)
		{
			uint64_t llo[3], lhi[3];
			int ldep[3];
			llo[0] = leaves[l].lo[0]; lhi[0] = llo[0] + 1;
			llo[1] = leaves[l].lo[1]; lhi[1] = llo[1] + 1;
			ldep[0] = ldep[1] = leaves[l].d;
			memcpy(w->lgu, w->gbu, sizeof(qaws_exact_int) * gs);
			memcpy(w->lgv, w->gbv, sizeof(qaws_exact_int) * gs);
			if (leaves[l].d)
			{
				TRY(sc_restrict2(w->lgu, Dg, Dg, llo, lhi, ldep));
				TRY(sc_restrict2(w->lgv, Dg, Dg, llo, lhi, ldep));
			}
			for (a = 0; a <= Dg; a++)
				for (b = 0; b <= Dg; b++)
					for (k = 0; k < 2; k++)
					{
						w->F[(a * (Dg + 1) + b) * 2 + k] = w->lgu[a * (Dg + 1) + b];
						w->F[size + (a * (Dg + 1) + b) * 2 + k] = w->lgv[a * (Dg + 1) + b];
					}
			for (a = 0; a <= Dg; a++)
				for (b = 0; b <= Dg; b++)
				{
					qaws_exact_int_from_i64(&w->F[2 * size + (a * (Dg + 1) + b) * 2], -1);
					qaws_exact_int_from_i64(&w->F[2 * size + (a * (Dg + 1) + b) * 2 + 1], 2);
				}
			n3[0] = Dg;
			n3[1] = Dg;
			n3[2] = 1;
			{
				static unsigned int const target[3] = { 22, 22, 0 };
				s3 = qaws_exact_solve3_to(n3, w->F, w->boxes, 64, &nb, 200000, target);
			}
			if (s3 == QAWS_STATUS_ALLOCATION_FAILURE)
				return s3;
			if (s3 != QAWS_STATUS_OK)
			{
				whole = 1;
				break;
			}
			for (k = 0; k < nb; k++)
			{
				/* leaf-local y to patch-local x = (leaf + y) 2^-d */
				qaws_exact_box3 const* B = &w->boxes[k];
				uint64_t plo[3], phi[3];
				int pdep[3];
				unsigned int dd;
				for (dd = 0; dd < 2; dd++)
				{
					pdep[dd] = leaves[l].d + B->dep[dd];
					plo[dd] = (llo[dd] << B->dep[dd]) + B->lo[dd];
					phi[dd] = (llo[dd] << B->dep[dd]) + B->hi[dd];
				}
				plo[2] = 0; phi[2] = 1; pdep[2] = 0;
				TRY(sc_enclose(w, De, plo, phi, pdep, &cand.lo, &cand.hi));
				cand.whole = 0;
				cand.ulo = sc_param(S->ub[iu], S->ub[iu + 1], S->u_shift, plo[0], pdep[0], 0);
				cand.uhi = sc_param(S->ub[iu], S->ub[iu + 1], S->u_shift, phi[0], pdep[0], 1);
				cand.vlo = sc_param(S->vb[iv], S->vb[iv + 1], S->v_shift, plo[1], pdep[1], 0);
				cand.vhi = sc_param(S->vb[iv], S->vb[iv + 1], S->v_shift, phi[1], pdep[1], 1);
				TRY(sc_push(L, &cand));
			}
		}
		if (whole)
		{
			/* not localised or not isolated: the whole patch is one candidate */
			cand.whole = 1;
			cand.ulo = u0; cand.uhi = u1; cand.vlo = v0; cand.vhi = v1;
			TRY(sc_enclose(w, De, NULL, NULL, NULL, &cand.lo, &cand.hi));
			TRY(sc_push(L, &cand));
		}
	}

	/* edges as exact curves: v = v0 (first row), v = v1, u = u0 (first column), u = u1;
	   corners through the edges along u */
	for (k = 0; k < 4; k++)
	{
		int along_u = k < 2, high = k % 2;
		unsigned int deg = along_u ? p : qd, i, j, e;
		qaws_exact_span sp;
		if (!high && (along_u ? iv > 0 : iu > 0))
			continue;
		for (i = 0; i <= deg; i++)
		{
			unsigned int ia = along_u ? i : (high ? p : 0), ib = along_u ? (high ? qd : 0) : i;
			for (j = 0; j < 4; j++)
				w->edge[i * 4 + j] = h[(ia * (qd + 1) + ib) * 4 + j];
		}
		sp.degree = deg;
		sp.a = along_u ? S->ub[iu] : S->vb[iv];
		sp.b = along_u ? S->ub[iu + 1] : S->vb[iv + 1];
		sp.h = w->edge;
		w->cl.n = 0;
		st = cl_span_candidates(&sp, 0, 0, along_u && iu == 0, along_u, q, &w->cw, &w->cl);
		if (st != QAWS_STATUS_OK)
			return st;
		for (e = 0; e < w->cl.n; e++)
		{
			cl_cand const* cc = &w->cl.c[e];
			double tlo, thi, fixed = along_u ? (high ? v1 : v0) : (high ? u1 : u0);
			int shift = along_u ? S->u_shift : S->v_shift, ex;
			if (cc->whole)
			{
				tlo = ldexp((double)sp.a, -shift);
				thi = ldexp((double)sp.b, -shift);
			}
			else
			{
				TRY(qaws_exact_span_param_to_double(&sp, shift, cc->r.index, cc->r.depth, &tlo, &ex));
				if (!ex) tlo = nextafter(tlo, -HUGE_VAL);
				TRY(qaws_exact_span_param_to_double(&sp, shift, cc->r.index + (cc->r.exact ? 0 : 1), cc->r.depth, &thi, &ex));
				if (!ex) thi = nextafter(thi, HUGE_VAL);
			}
			cand.whole = cc->whole;
			cand.lo = cc->lo;
			cand.hi = cc->hi;
			if (along_u)
			{
				cand.ulo = tlo; cand.uhi = thi; cand.vlo = cand.vhi = fixed;
			}
			else
			{
				cand.vlo = tlo; cand.vhi = thi; cand.ulo = cand.uhi = fixed;
			}
			TRY(sc_push(L, &cand));
		}
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_exact_surface_batch_closest(qaws_exact_ssi_batch_desc const* desc, double const* points, unsigned int point_count,
	qaws_exact_surface_closest_point* out_points, qaws_exact_batch_stats* out_stats)
{
	qaws_exact_batch_stats stats;
	cl_ring R;
	sc_list L;
	sc_work* w = NULL;
	qaws_bp_box* boxes = NULL;
	qaws_bp_grid* grid = NULL;
	unsigned int i, k, npatch = 0, iu, iv;
	int space;
	double scale;
	qaws_status st = QAWS_STATUS_OK;
	memset(&stats, 0, sizeof(stats));
	memset(&R, 0, sizeof(R));
	memset(&L, 0, sizeof(L));
	if (out_stats)
		*out_stats = stats;
	if (!desc || (desc->surface_count && !desc->surfaces) || (point_count && (!points || !out_points)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < point_count; i++)
	{
		memset(&out_points[i], 0, sizeof(out_points[i]));
		out_points[i].surface = QAWS_EXACT_CLOSEST_NONE;
	}
	if (!desc->surface_count || !point_count)
		return QAWS_STATUS_OK;
	for (i = 0; i < desc->surface_count; i++)
	{
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->surfaces[i]->space_exp2 != desc->surfaces[0]->space_exp2)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		npatch += desc->surfaces[i]->nu * desc->surfaces[i]->nv;
	}
	space = desc->surfaces[0]->space_exp2;
	scale = ldexp(1.0, space);
	R.dim = 3;
	R.all = (cl_span*)malloc(npatch * sizeof(cl_span));
	R.box = (double (*)[2][3])malloc(npatch * sizeof(double) * 6);
	R.stamp = (unsigned int*)malloc(npatch * sizeof(unsigned int));
	boxes = (qaws_bp_box*)malloc(npatch * sizeof(qaws_bp_box));
	w = (sc_work*)calloc(1, sizeof(sc_work));
	if (!R.all || !R.box || !R.stamp || !boxes || !w)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* cl_span reused: curve = surface, index = iu nv + iv */
	npatch = 0;
	for (i = 0; i < desc->surface_count; i++)
		for (iu = 0; iu < desc->surfaces[i]->nu; iu++)
			for (iv = 0; iv < desc->surfaces[i]->nv; iv++, npatch++)
			{
				unsigned int c;
				R.all[npatch].curve = i;
				R.all[npatch].index = iu * desc->surfaces[i]->nv + iv;
				R.stamp[npatch] = ~0u;
				qaws_exact_patch_box(desc->surfaces[i], iu, iv, R.box[npatch][0], R.box[npatch][1]);
				for (c = 0; c < 3; c++)
				{
					boxes[npatch].lo[c] = R.box[npatch][0][c];
					boxes[npatch].hi[c] = R.box[npatch][1][c];
				}
			}
	stats.patch_count = npatch;
	st = qaws_internal_grid_create(boxes, npatch, 3, &grid);
	for (i = 0; i < point_count && st == QAWS_STATUS_OK; i++)
	{
		cl_query q;
		double best_hi = HUGE_VAL;
		unsigned int r, s, win = 0, j;
		int sure;
		st = cl_query_make(points + 3 * i, 3, space, &q);
		if (st != QAWS_STATUS_OK)
			break;
		R.query = i;
		R.p = q.pd;
		R.ngot = 0;
		R.best_ub = HUGE_VAL;
		for (r = 0;; r++)
		{
			double bnd = qaws_internal_grid_ring(grid, q.pd, r, cl_ring_visit, &R);
			if (R.failed)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			if (bnd == HUGE_VAL || bnd * bnd > R.best_ub)
				break;
		}
		if (st != QAWS_STATUS_OK)
			break;
		/* the reachable patches by the lower bound of their distance (tighter
		   than their boxes): the nearest is usually solved first */
		for (s = 0; s < R.ngot && st == QAWS_STATUS_OK; s++)
		{
			qaws_exact_surface const* S = desc->surfaces[R.got[s].curve];
			double elo = 0;
			if (R.got[s].lb > R.best_ub)
				continue;
			st = sc_patch_bound(S, R.got[s].index / S->nv, R.got[s].index % S->nv, &q, w, &elo);
			if (st == QAWS_STATUS_EXACT_UNSUPPORTED)
				st = QAWS_STATUS_OK;
			else if (st == QAWS_STATUS_OK && elo > R.got[s].lb)
				R.got[s].lb = elo;
		}
		if (st != QAWS_STATUS_OK)
			break;
		qsort(R.got, R.ngot, sizeof(cl_span), cl_cmp_span);
		L.n = 0;
		for (s = 0; s < R.ngot && st == QAWS_STATUS_OK; s++)
		{
			qaws_exact_surface const* S = desc->surfaces[R.got[s].curve];
			unsigned int c0 = L.n;
			if (R.got[s].lb > best_hi * (1 + 1e-12))
				break;
			stats.candidate_count++;
			st = sc_patch_candidates(S, R.got[s].curve, R.got[s].index / S->nv, R.got[s].index % S->nv, &q, best_hi, w, &L);
			for (j = c0; j < L.n && st == QAWS_STATUS_OK; j++)
				if (L.c[j].hi < best_hi)
					best_hi = L.c[j].hi;
		}
		if (st != QAWS_STATUS_OK)
			break;
		if (!L.n)
			continue;
		/* one point found twice (an interior root on a patch edge): keep one */
		for (j = 0; j < L.n; j++)
			for (k = j + 1; k < L.n; k++)
				if (L.c[j].surface == L.c[k].surface && !L.c[j].whole && !L.c[k].whole && L.c[j].ulo <= L.c[k].uhi && L.c[k].ulo <= L.c[j].uhi
					&& L.c[j].vlo <= L.c[k].vhi && L.c[k].vlo <= L.c[j].vhi)
				{
					L.c[k] = L.c[--L.n];
					k--;
				}
		win = 0;
		for (j = 1; j < L.n; j++)
			if (L.c[j].hi < L.c[win].hi)
				win = j;
		sure = !L.c[win].whole;
		for (j = 0; j < L.n && sure; j++)
			if (j != win && L.c[j].lo <= L.c[win].hi)
				sure = 0;
		{
			sc_cand2 const* c = &L.c[win];
			qaws_exact_surface_closest_point* o = &out_points[i];
			o->surface = c->surface;
			o->u_lo = c->ulo; o->u_hi = c->uhi;
			o->v_lo = c->vlo; o->v_hi = c->vhi;
			o->distance_lo = nextafter(sqrt(c->lo > 0 ? c->lo : 0), -HUGE_VAL) * scale;
			o->distance_hi = nextafter(sqrt(c->hi), HUGE_VAL) * scale;
			if (o->distance_lo < 0)
				o->distance_lo = 0;
			o->certified = sure;
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
	if (w)
		free(w->cl.c);
	free(w);
	return st;
}
