/*
 * Uniform sampling of curved regions (see qaws_region_sample.h).
 *
 *   1. region: the paths united under the fill rule (qaws_clip), so the
 *      boundary is simple: outer loops counter-clockwise, holes clockwise
 *   2. pieces: every curve cut at its span ends and where dx/dt changes
 *      sign; each piece is x-monotone, sigma = +1 when x grows along it
 *   3. slabs: the pieces' end abscissas sorted and merged; each piece spans
 *      a run of slabs, and A at every slab's left side is summed up
 *   4. table: the chord length (or the density along the chords) g(x) in
 *      panels per slab, 8-point Gauss each, kept as Legendre series
 *   5. map: u * total located among the panels, x by the panel's series
 *      inverted; y in the chords at x, the pieces' y sorted and paired by
 *      the winding
 */

#include "qaws_region_sample.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_path.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_parallel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define RS_AREA_REL 1e-6
#else
#define RS_AREA_REL 1e-13
#endif
#define RS_ROOT_SAMPLES 32      /* dx/dt samples per span for its sign changes */
#define RS_MERGE_REL 1e-12      /* abscissas this close (of the box) are one */

typedef struct rs_piece
{
	qaws_curve const* c;
	double t0, t1;              /* t0 < t1 */
	double x0, x1;              /* x at t0, t1 */
	int sigma;                  /* +1 when x grows with t */
	unsigned int s0, s1;        /* slabs [s0, s1) it spans */
} rs_piece;

struct qaws_region_sampler
{
	qaws_clip_result* region;   /* owns the curves */
	rs_piece* pieces;
	unsigned int npieces;
	double* xs;                 /* critical abscissas, nxs of them */
	unsigned int nxs;
	double* acum;               /* A at xs[i] */
	unsigned int* first;        /* slab s: active pieces active[first[s] .. first[s + 1]) */
	unsigned int* active;
	double area, tol, ext;
	double ylo, yhi;            /* the region's y extent */
	/* the density, and its table: slab k's panels pfirst[k] .. pfirst[k + 1],
	   uniform in s where x = x_k + w (3 s^2 - 2 s^3) (flat at the slab's ends,
	   where chords grow like square roots); pmass[p] the mass left of panel
	   p's start; pcoef[8 p ..] the integrand over panel p in Legendre
	   polynomials, from its 8 Gauss values (the map inverts this interpolant) */
	qaws_scalar (*rho)(void* user, qaws_scalar x, qaws_scalar y);
	void* rho_user;
	unsigned int* pfirst;
	double* pmass;
	double* pcoef;
	double mass, rhomax;
	qaws_batch_executor const* executor;
};

static qaws_status rs_eval(qaws_curve const* c, double t, double* p, double* d)
{
	qaws_eval_result_2d r;
	qaws_status s = qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | (d ? QAWS_EVAL_FLAG_D1 : 0));
	if (s != QAWS_STATUS_OK)
	{
		/* (zeros: no caller reads garbage) */
		p[0] = p[1] = 0;
		if (d) d[0] = d[1] = 0;
		return s;
	}
	p[0] = r.position.x; p[1] = r.position.y;
	if (d) { d[0] = r.d1.x; d[1] = r.d1.y; }
	return QAWS_STATUS_OK;
}

/* the parameter of a piece at abscissa x (x within the piece's range):
   Newton on x(t) - x kept in a shrinking bracket */
static double rs_param(rs_piece const* pc, double x)
{
	double lo = pc->t0, hi = pc->t1, t;
	unsigned int it;
	if (pc->x1 == pc->x0)
		return 0.5 * (lo + hi);
	t = pc->t0 + (pc->t1 - pc->t0) * (x - pc->x0) / (pc->x1 - pc->x0);
	if (!(t >= lo && t <= hi)) t = 0.5 * (lo + hi);
	for (it = 0; it < 80; it++)
	{
		double p[2], d[2], f, nt;
		if (rs_eval(pc->c, t, p, d) != QAWS_STATUS_OK)
			break;
		f = (p[0] - x) * pc->sigma;      /* grows with t */
		if (f > 0) hi = t; else if (f < 0) lo = t; else return t;
		nt = d[0] != 0 ? t - (p[0] - x) / d[0] : 0.5 * (lo + hi);
		if (!(nt > lo && nt < hi))
			nt = 0.5 * (lo + hi);
		if (fabs(nt - t) <= (pc->t1 - pc->t0) * 1e-16 || hi - lo <= (pc->t1 - pc->t0) * 1e-16)
			return nt;
		t = nt;
	}
	return t;
}

/* integral of y dx along the piece from abscissa xa to xb: [x y] - integral x dy */
static double rs_int(qaws_region_sampler const* s, rs_piece const* pc, double xa, double xb, double* yb)
{
	double ta = rs_param(pc, xa), tb = rs_param(pc, xb), pa[2], pb[2], a = 0;
	rs_eval(pc->c, ta, pa, NULL);
	rs_eval(pc->c, tb, pb, NULL);
	if (yb) *yb = pb[1];
	if (ta != tb)
		qaws_internal_curve_area_2d(pc->c, ta, tb, s->tol, &a);
	return pb[0] * pb[1] - pa[0] * pa[1] - a;
}

static int rs_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : (x > y ? 1 : 0);
}

static qaws_status rs_add_piece(rs_piece** pcs, unsigned int* n, unsigned int* cap, qaws_curve const* c, double t0, double t1)
{
	double p0[2], p1[2];
	rs_piece* pc;
	qaws_status st;
	if (!(t1 > t0))
		return QAWS_STATUS_OK;
	st = rs_eval(c, t0, p0, NULL);
	if (st == QAWS_STATUS_OK) st = rs_eval(c, t1, p1, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (p0[0] == p1[0])
		return QAWS_STATUS_OK;      /* vertical: no area to its side */
	if (*n == *cap)
	{
		unsigned int nc = *cap ? 2 * *cap : 64;
		rs_piece* g = (rs_piece*)realloc(*pcs, sizeof(rs_piece) * nc);
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		*pcs = g;
		*cap = nc;
	}
	pc = &(*pcs)[(*n)++];
	memset(pc, 0, sizeof(*pc));
	pc->c = c;
	pc->t0 = t0; pc->t1 = t1;
	pc->x0 = p0[0]; pc->x1 = p1[0];
	pc->sigma = p1[0] > p0[0] ? 1 : -1;
	return QAWS_STATUS_OK;
}

/* a curve cut where dx/dt changes sign, span by span */
static qaws_status rs_cut_curve(qaws_curve const* c, rs_piece** pcs, unsigned int* n, unsigned int* cap)
{
	unsigned int sp, k;
	qaws_status st = QAWS_STATUS_OK;
	unsigned int ns = c->span_count ? c->span_count : 1;
	for (sp = 0; sp < ns && st == QAWS_STATUS_OK; sp++)
	{
		double s0 = c->span_count ? c->span_boundaries[sp] : c->parameter_range.min_value;
		double s1 = c->span_count ? c->span_boundaries[sp + 1] : c->parameter_range.max_value;
		double prev_t = s0, start = s0, p[2], d[2], prev_d;
		if (!(s1 > s0))
			continue;
		st = rs_eval(c, s0, p, d);
		prev_d = d[0];
		for (k = 1; k <= RS_ROOT_SAMPLES && st == QAWS_STATUS_OK; k++)
		{
			double t = s0 + (s1 - s0) * k / RS_ROOT_SAMPLES;
			st = rs_eval(c, t, p, d);
			if (st != QAWS_STATUS_OK)
				break;
			if ((prev_d > 0 && d[0] < 0) || (prev_d < 0 && d[0] > 0))
			{
				/* the root of dx/dt between: bisection */
				double lo = prev_t, hi = t, flo = prev_d;
				unsigned int it;
				for (it = 0; it < 64; it++)
				{
					double m = 0.5 * (lo + hi), pm[2], dm[2];
					rs_eval(c, m, pm, dm);
					if ((dm[0] > 0) == (flo > 0)) { lo = m; flo = dm[0]; } else hi = m;
				}
				st = rs_add_piece(pcs, n, cap, c, start, 0.5 * (lo + hi));
				start = 0.5 * (lo + hi);
			}
			if (d[0] != 0)
				prev_d = d[0];
			prev_t = t;
		}
		if (st == QAWS_STATUS_OK)
			st = rs_add_piece(pcs, n, cap, c, start, s1);
	}
	return st;
}

/* the merged abscissa index of x (xs sorted, merged within eps) */
static unsigned int rs_xindex(double const* xs, unsigned int n, double x)
{
	unsigned int lo = 0, hi = n - 1;
	while (lo < hi)
	{
		unsigned int m = (lo + hi) / 2;
		if (xs[m] < x) lo = m + 1; else hi = m;
	}
	/* the nearer of lo - 1 and lo */
	if (lo > 0 && fabs(xs[lo - 1] - x) < fabs(xs[lo] - x))
		lo--;
	return lo;
}

static qaws_status rs_build_density(qaws_region_sampler* s);

void qaws_region_sampler_destroy(qaws_region_sampler* s)
{
	if (!s)
		return;
	qaws_clip_result_destroy(s->region);
	free(s->pieces);
	free(s->xs);
	free(s->acum);
	free(s->first);
	free(s->active);
	free(s->pfirst);
	free(s->pmass);
	free(s->pcoef);
	free(s);
}

qaws_status qaws_region_sampler_create(qaws_region_sampler** out, qaws_region_sampler_desc const* desc)
{
	qaws_region_sampler* s;
	qaws_status st;
	unsigned int i, k, np, cap = 0;
	double* raw = NULL;
	if (!desc || !out || (desc->path_count && !desc->paths))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	s = (qaws_region_sampler*)calloc(1, sizeof(*s));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* 1. the region, its boundary simple */
	st = qaws_clip_boolean(&s->region, desc->paths, NULL, desc->path_count, 0, QAWS_CLIP_UNION, desc->fill_rule);
	np = st == QAWS_STATUS_OK ? qaws_clip_result_get_path_count(s->region) : 0;
	/* 2. x-monotone pieces */
	for (i = 0; st == QAWS_STATUS_OK && i < np; i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(&p, s->region, i);
		for (k = 0; k < p.curve_count && st == QAWS_STATUS_OK; k++)
			st = rs_cut_curve(p.curves[k], &s->pieces, &s->npieces, &cap);
	}
	if (st == QAWS_STATUS_OK && s->npieces == 0)
		st = QAWS_STATUS_DEGENERATE_CURVE;
	/* 3. slabs */
	if (st == QAWS_STATUS_OK)
	{
		double lo = 1e300, hi = -1e300, eps;
		unsigned int nraw = 2 * s->npieces, m = 0;
		raw = (double*)malloc(sizeof(double) * nraw);
		if (!raw)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; st == QAWS_STATUS_OK && i < s->npieces; i++)
		{
			raw[2 * i] = s->pieces[i].x0;
			raw[2 * i + 1] = s->pieces[i].x1;
		}
		if (st == QAWS_STATUS_OK)
		{
			qsort(raw, nraw, sizeof(double), rs_cmp_double);
			lo = raw[0]; hi = raw[nraw - 1];
			s->ext = hi - lo;
			{
				qaws_vec2 blo, bhi;
				double h;
				s->ylo = 1e300; s->yhi = -1e300;
				for (i = 0; i < np; i++)
				{
					qaws_path_2d p;
					qaws_clip_result_get_path(&p, s->region, i);
					if (qaws_path_compute_bounds_2d(&blo, &bhi, &p) == QAWS_STATUS_OK)
					{
						if ((double)blo.y < s->ylo) s->ylo = (double)blo.y;
						if ((double)bhi.y > s->yhi) s->yhi = (double)bhi.y;
					}
				}
				h = s->yhi > s->ylo ? s->yhi - s->ylo : 1;
				s->tol = desc->tolerance > 0 ? (double)desc->tolerance : RS_AREA_REL * (s->ext > 0 ? s->ext : 1) * h;
				s->executor = desc->executor;
			}
			eps = RS_MERGE_REL * (s->ext > 0 ? s->ext : 1);
			for (i = 0; i < nraw; i++)
				if (m == 0 || raw[i] - raw[m - 1] > eps)
					raw[m++] = raw[i];
			s->xs = raw;
			raw = NULL;
			s->nxs = m;
			if (m < 2)
				st = QAWS_STATUS_DEGENERATE_CURVE;
		}
	}
	/* each piece's slabs, then every slab's active pieces */
	if (st == QAWS_STATUS_OK)
	{
		unsigned int nslab = s->nxs - 1, total = 0;
		s->first = (unsigned int*)calloc(nslab + 2, sizeof(unsigned int));
		s->acum = (double*)malloc(sizeof(double) * s->nxs);
		if (!s->first || !s->acum)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; st == QAWS_STATUS_OK && i < s->npieces; i++)
		{
			rs_piece* pc = &s->pieces[i];
			unsigned int a = rs_xindex(s->xs, s->nxs, pc->x0 < pc->x1 ? pc->x0 : pc->x1);
			unsigned int b = rs_xindex(s->xs, s->nxs, pc->x0 < pc->x1 ? pc->x1 : pc->x0);
			pc->s0 = a; pc->s1 = b;
			for (k = a; k < b; k++) s->first[k + 1]++;
			total += b - a;
		}
		if (st == QAWS_STATUS_OK)
		{
			unsigned int* fill = (unsigned int*)calloc(nslab + 1, sizeof(unsigned int));
			for (k = 0; k < nslab; k++) s->first[k + 1] += s->first[k];
			s->active = (unsigned int*)malloc(sizeof(unsigned int) * (total + 1));
			if (!fill || !s->active)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; st == QAWS_STATUS_OK && i < s->npieces; i++)
				for (k = s->pieces[i].s0; k < s->pieces[i].s1; k++)
					s->active[s->first[k] + fill[k]++] = i;
			free(fill);
		}
		/* A at every abscissa */
		if (st == QAWS_STATUS_OK)
		{
			s->acum[0] = 0;
			for (k = 0; k < nslab; k++)
			{
				double a = 0;
				unsigned int j;
				for (j = s->first[k]; j < s->first[k + 1]; j++)
				{
					rs_piece const* pc = &s->pieces[s->active[j]];
					a -= pc->sigma * rs_int(s, pc, s->xs[k], s->xs[k + 1], NULL);
				}
				s->acum[k + 1] = s->acum[k] + a;
			}
			s->area = s->acum[nslab];
			if (!(s->area > 0))
				st = QAWS_STATUS_DEGENERATE_CURVE;
			if (st == QAWS_STATUS_OK)
			{
				s->rho = desc->density;
				s->rho_user = desc->density_user;
				st = rs_build_density(s);
			}
		}
	}
	free(raw);
	if (st != QAWS_STATUS_OK)
	{
		qaws_region_sampler_destroy(s);
		return st;
	}
	*out = s;
	return QAWS_STATUS_OK;
}

qaws_scalar qaws_region_sampler_get_area(qaws_region_sampler const* s)
{
	return s ? (qaws_scalar)s->area : 0;
}

/* the slab holding x (xs[0] <= x <= xs[n - 1]) */
static unsigned int rs_slab(qaws_region_sampler const* s, double x)
{
	unsigned int lo = 0, hi = s->nxs - 2;
	while (lo < hi)
	{
		unsigned int m = (lo + hi + 1) / 2;
		if (s->xs[m] <= x) lo = m; else hi = m - 1;
	}
	return lo;
}

/* A(x) within slab k, and the chord length there */
static double rs_area_in(qaws_region_sampler const* s, unsigned int k, double x, double* chord)
{
	double a = s->acum[k], l = 0;
	unsigned int j;
	for (j = s->first[k]; j < s->first[k + 1]; j++)
	{
		rs_piece const* pc = &s->pieces[s->active[j]];
		double y;
		a -= pc->sigma * rs_int(s, pc, s->xs[k], x, &y);
		l -= pc->sigma * y;
	}
	if (chord) *chord = l;
	return a;
}

qaws_status qaws_region_sampler_area_left(qaws_scalar* out, qaws_region_sampler const* s, qaws_scalar x)
{
	double xd = (double)x;
	if (!s || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (xd <= s->xs[0]) { *out = 0; return QAWS_STATUS_OK; }
	if (xd >= s->xs[s->nxs - 1]) { *out = (qaws_scalar)s->area; return QAWS_STATUS_OK; }
	*out = (qaws_scalar)rs_area_in(s, rs_slab(s, xd), xd, NULL);
	return QAWS_STATUS_OK;
}

typedef struct rs_cross { double y; int sigma; } rs_cross;

static int rs_cmp_cross(void const* a, void const* b)
{
	double x = ((rs_cross const*)a)->y, y = ((rs_cross const*)b)->y;
	return x < y ? -1 : (x > y ? 1 : 0);
}

/* ------------------------------------------------------------------------ */
/*  Density: the mass left of x and along the chords                        */
/* ------------------------------------------------------------------------ */

/* Gauss-Legendre, 8 points on [-1, 1] */
static double const rs_gx[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
static double const rs_gw[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
	0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };

#define RS_PANELS 256            /* density panels over the region's width */
#define RS_YSEGS 16              /* Gauss segments over the region's height */

/* the sorted crossings of the slab's pieces with the vertical line at x;
   stack holds 64, more are allocated (the caller frees what differs) */
static rs_cross* rs_crossings(qaws_region_sampler const* s, unsigned int k, double x, rs_cross* stack, unsigned int* out_n)
{
	unsigned int nc = s->first[k + 1] - s->first[k], j;
	rs_cross* cr = stack;
	if (nc > 64)
	{
		cr = (rs_cross*)malloc(sizeof(rs_cross) * nc);
		if (!cr) { *out_n = 0; return stack; }
	}
	for (j = 0; j < nc; j++)
	{
		rs_piece const* pc = &s->pieces[s->active[s->first[k] + j]];
		double p[2];
		rs_eval(pc->c, rs_param(pc, x), p, NULL);
		cr[j].y = p[1];
		cr[j].sigma = pc->sigma;
	}
	qsort(cr, nc, sizeof(rs_cross), rs_cmp_cross);
	*out_n = nc;
	return cr;
}

/* Legendre coefficients (times r, the interval's half width) of the degree 7
   polynomial through the values f at the 8 Gauss nodes: its integral over
   the interval is the Gauss sum */
static void rs_legendre(double const* f, double r, double* c)
{
	unsigned int n, g;
	for (n = 0; n < 8; n++) c[n] = 0;
	for (g = 0; g < 8; g++)
	{
		double p0 = 1, p1 = rs_gx[g], wf = rs_gw[g] * f[g] * r;
		c[0] += wf;
		c[1] += wf * p1;
		for (n = 2; n < 8; n++)
		{
			double p2 = ((2 * n - 1) * rs_gx[g] * p1 - (n - 1) * p0) / n;
			c[n] += wf * p2;
			p0 = p1;
			p1 = p2;
		}
	}
	for (n = 0; n < 8; n++) c[n] *= 0.5 * (2 * n + 1);
}

/* the polynomial's integral from -1 to xi (the integral of P_n is
   (P_n+1 - P_n-1) / (2 n + 1)), and its value at xi */
static double rs_legendre_cum(double const* c, double xi, double* d)
{
	double p[9], sum = c[0] * (xi + 1), dv = c[0];
	unsigned int n;
	p[0] = 1;
	p[1] = xi;
	for (n = 1; n < 8; n++) p[n + 1] = ((2 * n + 1) * xi * p[n] - n * p[n - 1]) / (n + 1);
	for (n = 1; n < 8; n++)
	{
		sum += c[n] * (p[n + 1] - p[n - 1]) / (2 * n + 1);
		dv += c[n] * p[n];
	}
	if (d) *d = dv;
	return sum;
}

/* xi in [-1, 1] where the integral from -1 reaches target: Newton, bisection
   when it leaves the bracket */
static double rs_legendre_solve(double const* c, double target)
{
	double lo = -1, hi = 1, total = rs_legendre_cum(c, 1, NULL), xi;
	unsigned int it;
	xi = total > 0 ? -1 + 2 * target / total : 0;
	if (xi < -1) xi = -1;
	if (xi > 1) xi = 1;
	for (it = 0; it < 80; it++)
	{
		double d, f = rs_legendre_cum(c, xi, &d) - target, nx;
		if (f == 0)
			break;
		if (f > 0) hi = xi; else lo = xi;
		nx = d > 0 ? xi - f / d : 0.5 * (lo + hi);
		if (!(nx > lo && nx < hi))
			nx = 0.5 * (lo + hi);
		if (fabs(nx - xi) <= 1e-15 || hi - lo <= 1e-15)
		{
			xi = nx;
			break;
		}
		xi = nx;
	}
	return xi;
}

/* the density's integral over [c0, c1] at x (8-point Gauss), its values at
   the nodes in f */
static double rs_rho_seg(qaws_region_sampler const* s, double x, double c0, double c1, double* rmax, double* f)
{
	double m = 0.5 * (c0 + c1), r = 0.5 * (c1 - c0), sum = 0;
	unsigned int g;
	for (g = 0; g < 8; g++)
	{
		double v = (double)s->rho(s->rho_user, (qaws_scalar)x, (qaws_scalar)(m + r * rs_gx[g]));
		if (rmax && v > *rmax) *rmax = v;
		if (f) f[g] = v;
		sum += rs_gw[g] * r * v;
	}
	return sum;
}

/* the density's integral over [a, b] at x, composite Gauss by the region's height */
static double rs_rho_line(qaws_region_sampler const* s, double x, double a, double b, double* rmax)
{
	double h = (s->yhi - s->ylo) / RS_YSEGS, sum = 0;
	unsigned int n = (unsigned int)ceil((b - a) / h), i;
	if (!(b > a))
		return 0;
	if (n < 1) n = 1;
	for (i = 0; i < n; i++)
		sum += rs_rho_seg(s, x, a + (b - a) * i / n, a + (b - a) * (i + 1) / n, rmax, NULL);
	return sum;
}

/* y in [a, b] with the density's integral from a to y equal to target: the
   segment of rs_rho_line holding it, then its interpolant inverted */
static double rs_rho_solve(qaws_region_sampler const* s, double x, double a, double b, double target)
{
	double h = (s->yhi - s->ylo) / RS_YSEGS;
	unsigned int n = (unsigned int)ceil((b - a) / h), i;
	if (!(b > a))
		return a;
	if (n < 1) n = 1;
	for (i = 0; i < n; i++)
	{
		double c0 = a + (b - a) * i / n, c1 = a + (b - a) * (i + 1) / n, f[8], c[8], m;
		m = rs_rho_seg(s, x, c0, c1, NULL, f);
		if (target > m && i + 1 < n)
		{
			target -= m;
			continue;
		}
		if (!(m > 0))
			return 0.5 * (c0 + c1);
		rs_legendre(f, 0.5 * (c1 - c0), c);
		return 0.5 * (c0 + c1) + 0.5 * (c1 - c0) * rs_legendre_solve(c, target < m ? target : m);
	}
	return b;
}

/* g(x): the density along the chords at x (slab k); without one their length */
static double rs_g(qaws_region_sampler const* s, unsigned int k, double x, double* rmax)
{
	rs_cross stack[64], *cr;
	unsigned int n, j;
	double w = 0, sum = 0;
	cr = rs_crossings(s, k, x, stack, &n);
	for (j = 0; j + 1 < n; j++)
	{
		w += cr[j].sigma;
		if (w > 0) sum += s->rho ? rs_rho_line(s, x, cr[j].y, cr[j + 1].y, rmax) : cr[j + 1].y - cr[j].y;
	}
	if (cr != stack) free(cr);
	return sum;
}

/* x of the slab parameter s in [0, 1], and its derivative */
static double rs_xs(qaws_region_sampler const* s, unsigned int k, double t, double* dx)
{
	double w = s->xs[k + 1] - s->xs[k];
	if (dx) *dx = w * 6 * t * (1 - t);
	return s->xs[k] + w * t * t * (3 - 2 * t);
}

/* the slab parameter of x (the smoothstep inverted, by bisection then Newton) */
static double rs_sx(qaws_region_sampler const* s, unsigned int k, double x)
{
	double w = s->xs[k + 1] - s->xs[k], f = w > 0 ? (x - s->xs[k]) / w : 0, lo = 0, hi = 1, t = f;
	unsigned int it;
	if (f <= 0) return 0;
	if (f >= 1) return 1;
	for (it = 0; it < 80; it++)
	{
		double v = t * t * (3 - 2 * t) - f, d = 6 * t * (1 - t), nt;
		if (v > 0) hi = t; else lo = t;
		nt = d > 0 ? t - v / d : 0.5 * (lo + hi);
		if (!(nt > lo && nt < hi)) nt = 0.5 * (lo + hi);
		if (fabs(nt - t) < 1e-16) return nt;
		t = nt;
	}
	return t;
}

/* the mass of slab k between the parameters ta < tb: Gauss on the smooth
   integrand, its values at the nodes in f */
static double rs_mass_between(qaws_region_sampler const* s, unsigned int k, double ta, double tb, double* rmax, double* f)
{
	double m = 0.5 * (ta + tb), r = 0.5 * (tb - ta), sum = 0;
	unsigned int g;
	for (g = 0; g < 8; g++)
	{
		double dx, x = rs_xs(s, k, m + r * rs_gx[g], &dx), v = rs_g(s, k, x, rmax) * dx;
		if (f) f[g] = v;
		sum += rs_gw[g] * r * v;
	}
	return sum;
}

/* the slab of panel p */
static unsigned int rs_panel_slab(qaws_region_sampler const* s, unsigned int p)
{
	unsigned int c = 0, d = s->nxs - 2;
	while (c < d)
	{
		unsigned int m = (c + d + 1) / 2;
		if (s->pfirst[m] <= p) c = m; else d = m - 1;
	}
	return c;
}

typedef struct rs_build_ctx
{
	qaws_region_sampler* s;
	double rmax[256];           /* per chunk */
} rs_build_ctx;

/* panels [b, e): their masses in pmass[p + 1], their interpolants */
static qaws_status rs_build_job(void* ctx, unsigned int chunk, unsigned int b, unsigned int e)
{
	rs_build_ctx* c = (rs_build_ctx*)ctx;
	qaws_region_sampler* s = c->s;
	unsigned int p, k = b < e ? rs_panel_slab(s, b) : 0;
	for (p = b; p < e; p++)
	{
		unsigned int np, i;
		double f[8];
		while (s->pfirst[k + 1] <= p) k++;
		np = s->pfirst[k + 1] - s->pfirst[k];
		i = p - s->pfirst[k];
		s->pmass[p + 1] = rs_mass_between(s, k, (double)i / np, (double)(i + 1) / np, &c->rmax[chunk], f);
		rs_legendre(f, 0.5 / np, s->pcoef + 8 * p);
	}
	return QAWS_STATUS_OK;
}

static qaws_status rs_build_density(qaws_region_sampler* s)
{
	unsigned int nslab = s->nxs - 1, k, total = 0, p;
	double width = s->xs[s->nxs - 1] - s->xs[0], rmax = 0;
	rs_build_ctx* c;
	qaws_status st;
	s->pfirst = (unsigned int*)malloc(sizeof(unsigned int) * (nslab + 1));
	if (!s->pfirst)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < nslab; k++)
	{
		unsigned int np = (unsigned int)ceil(RS_PANELS * (s->xs[k + 1] - s->xs[k]) / width);
		s->pfirst[k] = total;
		total += np < 2 ? 2 : np;
	}
	s->pfirst[nslab] = total;
	s->pmass = (double*)malloc(sizeof(double) * (total + 1));
	s->pcoef = (double*)malloc(sizeof(double) * 8 * total);
	c = (rs_build_ctx*)calloc(1, sizeof(rs_build_ctx));
	if (!s->pmass || !s->pcoef || !c)
	{
		free(c);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the panels in parallel, then their masses summed in order */
	c->s = s;
	st = qaws_internal_parallel(s->executor, total, 2, rs_build_job, c);
	for (k = 0; k < 256; k++)
		if (c->rmax[k] > rmax) rmax = c->rmax[k];
	free(c);
	if (st != QAWS_STATUS_OK)
		return st;
	s->pmass[0] = 0;
	for (p = 0; p < total; p++)
		s->pmass[p + 1] += s->pmass[p];
	s->mass = s->pmass[total];
	/* the density's bound for rejection: the largest value met, and a grid
	   over the box, with a margin */
	if (s->rho)
	{
		unsigned int i, j;
		for (j = 0; j < 64; j++)
			for (i = 0; i < 64; i++)
			{
				double x = s->xs[0] + width * (i + 0.5) / 64, y = s->ylo + (s->yhi - s->ylo) * (j + 0.5) / 64;
				double v = (double)s->rho(s->rho_user, (qaws_scalar)x, (qaws_scalar)y);
				if (v > rmax) rmax = v;
			}
	}
	s->rhomax = 1.25 * rmax;
	return s->mass > 0 ? QAWS_STATUS_OK : QAWS_STATUS_DEGENERATE_CURVE;
}

/* the mass left of the slab parameter t in slab k */
static double rs_mass_at(qaws_region_sampler const* s, unsigned int k, double t)
{
	unsigned int np = s->pfirst[k + 1] - s->pfirst[k], p = (unsigned int)(t * np);
	if (p >= np) p = np - 1;
	return s->pmass[s->pfirst[k] + p] + rs_mass_between(s, k, (double)p / np, t, NULL, NULL);
}

qaws_scalar qaws_region_sampler_get_mass(qaws_region_sampler const* s)
{
	return s ? (qaws_scalar)(s->rho ? s->mass : s->area) : 0;
}

qaws_status qaws_region_sampler_mass_left(qaws_scalar* out, qaws_region_sampler const* s, qaws_scalar x)
{
	double xd = (double)x;
	unsigned int k;
	if (!s || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!s->rho)
		return qaws_region_sampler_area_left(out, s, x);
	if (xd <= s->xs[0]) { *out = 0; return QAWS_STATUS_OK; }
	if (xd >= s->xs[s->nxs - 1]) { *out = (qaws_scalar)s->mass; return QAWS_STATUS_OK; }
	k = rs_slab(s, xd);
	*out = (qaws_scalar)rs_mass_at(s, k, rs_sx(s, k, xd));
	return QAWS_STATUS_OK;
}

/* x with the mass (the area without a density) left of it u times the
   whole: the panel holding it, then its interpolant inverted; k its slab */
static double rs_map_x(qaws_region_sampler const* s, double u, unsigned int* out_k)
{
	double target = u * s->mass;
	unsigned int a, b, p, k, np;
	if (target < 0) target = 0;
	if (target > s->mass) target = s->mass;
	/* the panel holding the target, then its slab */
	a = 0; b = s->pfirst[s->nxs - 1] - 1;
	while (a < b)
	{
		unsigned int m = (a + b + 1) / 2;
		if (s->pmass[m] <= target) a = m; else b = m - 1;
	}
	p = a;
	k = rs_panel_slab(s, p);
	np = s->pfirst[k + 1] - s->pfirst[k];
	*out_k = k;
	return rs_xs(s, k, (p - s->pfirst[k] + 0.5 * (1 + rs_legendre_solve(s->pcoef + 8 * p, target - s->pmass[p]))) / np, NULL);
}

/* the mass-keeping map: x by the panels' interpolants, y by the chords' */
static qaws_status rs_map_mass(qaws_region_sampler const* s, double u, double v, qaws_vec2* out)
{
	unsigned int k;
	double x = rs_map_x(s, u, &k);
	/* y: the mass along the chords at x */
	{
		rs_cross stack[64], *cr;
		double mstack[64], *mm = mstack;
		unsigned int n, j;
		double w = 0, total = 0, yt;
		cr = rs_crossings(s, k, x, stack, &n);
		if (n > 64)
		{
			mm = (double*)malloc(sizeof(double) * n);
			if (!mm)
			{
				if (cr != stack) free(cr);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
		}
		for (j = 0; j + 1 < n; j++)
		{
			w += cr[j].sigma;
			mm[j] = w > 0 ? rs_rho_line(s, x, cr[j].y, cr[j + 1].y, NULL) : 0;
			total += mm[j];
		}
		yt = v * total;
		out->x = (qaws_scalar)x;
		out->y = n ? (qaws_scalar)cr[0].y : 0;
		w = 0;
		for (j = 0; j + 1 < n; j++)
		{
			w += cr[j].sigma;
			if (!(w > 0))
				continue;
			if (yt > mm[j] && j + 2 < n)
			{
				yt -= mm[j];
				out->y = (qaws_scalar)cr[j + 1].y;
				continue;
			}
			out->y = (qaws_scalar)rs_rho_solve(s, x, cr[j].y, cr[j + 1].y, yt);
			break;
		}
		if (mm != mstack) free(mm);
		if (cr != stack) free(cr);
	}
	return QAWS_STATUS_OK;
}

/* the area-keeping map */
static qaws_status rs_map_area(qaws_region_sampler const* s, qaws_scalar u, qaws_scalar v, qaws_vec2* out)
{
	double x, ytot = 0, yt, w = 0;
	unsigned int k, j, nc;
	rs_cross stackc[64], *cr = stackc;
	/* x: the table of the chord length, its interpolants inverted */
	x = rs_map_x(s, (double)u, &k);
	/* y: the chords at x */
	nc = s->first[k + 1] - s->first[k];
	if (nc > 64)
	{
		cr = (rs_cross*)malloc(sizeof(rs_cross) * nc);
		if (!cr)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (j = 0; j < nc; j++)
	{
		rs_piece const* pc = &s->pieces[s->active[s->first[k] + j]];
		double p[2];
		rs_eval(pc->c, rs_param(pc, x), p, NULL);
		cr[j].y = p[1];
		cr[j].sigma = pc->sigma;
	}
	qsort(cr, nc, sizeof(rs_cross), rs_cmp_cross);
	/* inside between two crossings when the winding below them is positive:
	   a piece running to +x has the region above it */
	for (j = 0; j + 1 < nc; j++)
	{
		w += cr[j].sigma;
		if (w > 0) ytot += cr[j + 1].y - cr[j].y;
	}
	yt = (double)v * ytot;
	out->x = (qaws_scalar)x;
	out->y = nc ? (qaws_scalar)cr[0].y : 0;
	w = 0;
	for (j = 0; j + 1 < nc; j++)
	{
		double len;
		w += cr[j].sigma;
		if (!(w > 0))
			continue;
		len = cr[j + 1].y - cr[j].y;
		out->y = (qaws_scalar)(cr[j].y + (yt < len ? yt : len));
		if (yt <= len)
			break;
		yt -= len;
	}
	if (cr != stackc)
		free(cr);
	return QAWS_STATUS_OK;
}

qaws_status qaws_region_sampler_map(qaws_vec2* out, qaws_region_sampler const* s, qaws_scalar u, qaws_scalar v)
{
	if (!s || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return s->rho ? rs_map_mass(s, (double)u, (double)v, out) : rs_map_area(s, u, v, out);
}

/* splitmix64: the value after z, and the next of a state */
static unsigned long long rs_mix(unsigned long long z)
{
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

static unsigned long long rs_next(unsigned long long* st)
{
	return rs_mix(*st += 0x9E3779B97F4A7C15ull);
}

static double rs_unit(unsigned long long r)
{
	return (double)(r >> 11) * (1.0 / 9007199254740992.0);
}

static double rs_frac(double v)
{
	return v - floor(v);
}

/* Kensler's hashed permutation of [0, l), and hashed float in [0, 1) */
static unsigned int rs_permute(unsigned int i, unsigned int l, unsigned int p)
{
	unsigned int w = l - 1;
	w |= w >> 1; w |= w >> 2; w |= w >> 4; w |= w >> 8; w |= w >> 16;
	do
	{
		i ^= p; i *= 0xe170893du; i ^= p >> 16; i ^= (i & w) >> 4;
		i ^= p >> 8; i *= 0x0929eb3fu; i ^= p >> 23; i ^= (i & w) >> 1;
		i *= 1u | p >> 27; i *= 0x6935fa69u; i ^= (i & w) >> 11; i *= 0x74dcb303u;
		i ^= (i & w) >> 2; i *= 0x9e501cc3u; i ^= (i & w) >> 2; i *= 0xc860a3dfu;
		i &= w; i ^= i >> 5;
	} while (i >= l);
	return (i + p) % l;
}

static double rs_hashfloat(unsigned int i, unsigned int p)
{
	i ^= p; i ^= i >> 17; i ^= i >> 10; i *= 0xb36534e5u; i ^= i >> 12;
	i ^= i >> 21; i *= 0x93fc4795u; i ^= 0xdf6e307fu; i ^= i >> 17; i *= 1u | p >> 18;
	return i * (1.0 / 4294967808.0);
}

static unsigned int rs_reverse_bits(unsigned int x)
{
	x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
	x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
	x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
	x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
	return (x >> 16) | (x << 16);
}

/* Owen scrambling by hashing (Burley, Practical Hash-based Owen Scrambling,
   2020, after Laine and Karras): every digit flipped by a hash of the digits
   above it, so the Sobol strata are kept and the lattice is gone */
static unsigned int rs_owen(unsigned int x, unsigned int seed)
{
	x = rs_reverse_bits(x);
	x += seed;
	x ^= x * 0x6c50b47cu;
	x ^= x * 0xb82f1e52u;
	x ^= x * 0xc7afe638u;
	x ^= x * 0x8d22f6e6u;
	return rs_reverse_bits(x);
}

static double rs_radical(unsigned int i, unsigned int base)
{
	double f = 1.0 / base, r = 0;
	while (i)
	{
		r += f * (i % base);
		i /= base;
		f /= base;
	}
	return r;
}

/* the (u, v) of point i of a pattern of the unit square */
typedef struct rs_pattern
{
	unsigned int kind, count, m, n, p;
	unsigned int shift0, shift1, dir[32];
	double off0, off1;
	unsigned long long state;
} rs_pattern;

static void rs_pattern_init(rs_pattern* q, unsigned int kind, unsigned int count, unsigned long long seed)
{
	unsigned int b;
	memset(q, 0, sizeof(*q));
	q->kind = kind;
	q->count = count;
	q->state = seed;
	q->shift0 = (unsigned int)(rs_next(&q->state) >> 32);
	q->shift1 = (unsigned int)(rs_next(&q->state) >> 32);
	q->off0 = rs_unit(rs_next(&q->state));
	q->off1 = rs_unit(rs_next(&q->state));
	q->p = (unsigned int)(rs_next(&q->state) >> 32);
	if (seed == 0)
	{
		q->shift0 = q->shift1 = 0;
		q->off0 = q->off1 = 0;
	}
	q->dir[0] = 0x80000000u;
	for (b = 1; b < 32; b++) q->dir[b] = q->dir[b - 1] ^ (q->dir[b - 1] >> 1);
	/* the multi-jitter grid: m columns, n rows. A factor pair of count with
	   sides within a factor of two fills every stratum; else the rounded-up
	   square, whose spare strata the permutation leaves out */
	{
		unsigned int c = count ? count : 1, a, lo = (unsigned int)ceil(sqrt(c / 2.0));
		q->m = 0;
		for (a = (unsigned int)floor(sqrt((double)c)); a >= lo && a > 0; a--)
			if (c % a == 0) { q->m = a; q->n = c / a; break; }
		if (!q->m)
		{
			q->m = (unsigned int)ceil(sqrt((double)c));
			q->n = (c + q->m - 1) / q->m;
		}
	}
}

static void rs_pattern_point(rs_pattern const* q, unsigned int i, double* u, double* v)
{
	switch (q->kind)
	{
	case QAWS_REGION_SAMPLE_SOBOL:
	{
		unsigned int x0 = 0, x1 = 0, k = i, b;
		for (b = 0; k; b++, k >>= 1)
			if (k & 1u) { x0 ^= 0x80000000u >> b; x1 ^= q->dir[b]; }
		if (q->shift0 || q->shift1) { x0 = rs_owen(x0, q->shift0); x1 = rs_owen(x1, q->shift1); }
		*u = (x0 + 0.5) / 4294967296.0;
		*v = (x1 + 0.5) / 4294967296.0;
		break;
	}
	case QAWS_REGION_SAMPLE_STRATIFIED:
	{
		/* Kensler, Correlated Multi-Jittered Sampling (2013), for any count */
		unsigned int s = rs_permute(i, q->count, q->p * 0x51633e2du), m = q->m, n = q->n;
		unsigned int sx = rs_permute(s % m, m, q->p * 0xa511e9b3u), sy = rs_permute(s / m, n, q->p * 0x63d83595u);
		double jx = rs_hashfloat(s, q->p * 0xa399d265u), jy = rs_hashfloat(s, q->p * 0x711ad6a5u);
		*u = ((s % m) + (sy + jx) / n) / m;
		*v = ((s / m) + (sx + jy) / m) / n;
		break;
	}
	case QAWS_REGION_SAMPLE_HALTON:
		*u = rs_frac(rs_radical(i + 1, 2) + q->off0);
		*v = rs_frac(rs_radical(i + 1, 3) + q->off1);
		break;
	case QAWS_REGION_SAMPLE_R2:
	{
		/* the plastic number g, g^3 = g + 1 */
		double const g = 1.32471795724474602596;
		*u = rs_frac(0.5 + (i + 1) / g + q->off0);
		*v = rs_frac(0.5 + (i + 1) / (g * g) + q->off1);
		break;
	}
	default:
		/* the stream from the seed, point i its draws 2 i + 1 and 2 i + 2 */
		*u = rs_unit(rs_mix(q->state + (2ull * i + 1) * 0x9E3779B97F4A7C15ull));
		*v = rs_unit(rs_mix(q->state + (2ull * i + 2) * 0x9E3779B97F4A7C15ull));
		break;
	}
	if (*u >= 1) *u = 1 - 1e-16;
	if (*v >= 1) *v = 1 - 1e-16;
}

/* ------------------------------------------------------------------------ */
/*  Batches: chunks of points through the executor, each point a function   */
/*  of its index alone, so the result does not depend on the executor        */
/* ------------------------------------------------------------------------ */

#define RS_GRAIN 64              /* points per chunk */

typedef struct rs_batch
{
	qaws_region_sampler const* s;
	rs_pattern const* q;        /* the points' (u, v): the pattern's from index first, */
	qaws_vec2 const* uv;        /* or given */
	unsigned int first;
	double box[4];              /* x0, w, y0, h: the pattern over the box */
	qaws_vec2* out;             /* mapped points */
	unsigned char* keep;        /* over the box: kept or not */
	double* rad;                /* blue noise: the radii at out, rmax given in rad0 */
	double rad0;
} rs_batch;

static int rs_keep(qaws_region_sampler const* s, double x, double y, unsigned int i, unsigned int p);

static qaws_status rs_batch_job(void* ctx, unsigned int chunk, unsigned int b, unsigned int e)
{
	rs_batch const* c = (rs_batch const*)ctx;
	unsigned int i;
	qaws_status st = QAWS_STATUS_OK;
	(void)chunk;
	for (i = b; i < e && st == QAWS_STATUS_OK; i++)
	{
		double u, v;
		if (c->rad)
		{
			/* radii as 1 / sqrt(rho), rmax at the mean, within a factor of 4 */
			double r = (double)c->s->rho(c->s->rho_user, c->out[i].x, c->out[i].y), f;
			f = r > 0 ? sqrt(c->s->mass / c->s->area / r) : 4;
			if (f > 4) f = 4;
			if (f < 0.25) f = 0.25;
			c->rad[i] = c->rad0 * f;
			continue;
		}
		if (c->uv) { u = (double)c->uv[i].x; v = (double)c->uv[i].y; }
		else rs_pattern_point(c->q, c->first + i, &u, &v);
		if (c->keep)
			c->keep[i] = (unsigned char)rs_keep(c->s, c->box[0] + c->box[1] * u, c->box[2] + c->box[3] * v, c->first + i, c->q->p);
		else
			st = qaws_region_sampler_map(&c->out[i], c->s, (qaws_scalar)u, (qaws_scalar)v);
	}
	return st;
}

static qaws_status rs_batch_run(rs_batch const* c, unsigned int count)
{
	return qaws_internal_parallel(c->s->executor, count, RS_GRAIN, rs_batch_job, (void*)c);
}

/* count pattern points (or given uv) through the map */
static qaws_status rs_map_batch(qaws_region_sampler const* s, rs_pattern const* q, qaws_vec2 const* uv, unsigned int count, qaws_vec2* out)
{
	rs_batch c;
	memset(&c, 0, sizeof(c));
	c.s = s;
	c.q = q;
	c.uv = uv;
	c.out = out;
	return rs_batch_run(&c, count);
}

qaws_status qaws_region_sampler_map_many(qaws_vec2* out, qaws_region_sampler const* s, qaws_vec2 const* uv, unsigned int count)
{
	if (!s || (count && (!uv || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	return rs_map_batch(s, NULL, uv, count, out);
}

/* ------------------------------------------------------------------------ */
/*  Blue noise: weighted sample elimination (Yuksel, Sample Elimination for */
/*  Generating Poisson Disk Sample Sets, 2015)                               */
/* ------------------------------------------------------------------------ */

typedef struct rs_elim
{
	qaws_vec2 const* pts;
	unsigned int n;
	double* w;                  /* weights */
	unsigned int* heap;         /* max-heap of point indices by weight */
	unsigned int* pos;          /* index -> heap position, ~0 once removed */
	unsigned int hn;
	unsigned int* cell_first;   /* grid: cell c holds cell_items[cell_first[c] .. cell_first[c + 1]) */
	unsigned int* cell_items;
	unsigned int gx, gy;
	double lo[2], cs, rmax, rmin;
	double* rad;                /* with a density: each point's radius (rmax at mean density) */
} rs_elim;

static void rs_heap_swap(rs_elim* e, unsigned int a, unsigned int b)
{
	unsigned int t = e->heap[a];
	e->heap[a] = e->heap[b];
	e->heap[b] = t;
	e->pos[e->heap[a]] = a;
	e->pos[e->heap[b]] = b;
}

static void rs_heap_down(rs_elim* e, unsigned int i)
{
	for (;;)
	{
		unsigned int l = 2 * i + 1, r = l + 1, m = i;
		if (l < e->hn && e->w[e->heap[l]] > e->w[e->heap[m]]) m = l;
		if (r < e->hn && e->w[e->heap[r]] > e->w[e->heap[m]]) m = r;
		if (m == i) return;
		rs_heap_swap(e, i, m);
		i = m;
	}
}

static double rs_elim_weight(rs_elim const* e, double d, unsigned int i, unsigned int j)
{
	double t, r = e->rad ? 0.5 * (e->rad[i] + e->rad[j]) : e->rmax, rmin = e->rad ? e->rmin * r / e->rmax : e->rmin;
	if (d < rmin) d = rmin;
	if (d >= 2 * r) return 0;
	t = 1 - d / (2 * r);
	t *= t; t *= t;
	return t * t;                /* alpha = 8 */
}

static unsigned int rs_elim_cell(rs_elim const* e, qaws_vec2 const* p)
{
	unsigned int cx = (unsigned int)((p->x - e->lo[0]) / e->cs), cy = (unsigned int)((p->y - e->lo[1]) / e->cs);
	if (cx >= e->gx) cx = e->gx - 1;
	if (cy >= e->gy) cy = e->gy - 1;
	return cy * e->gx + cx;
}

/* the live neighbours of i within 2 rmax: add their weights to i, or (sub)
   take i's weight off each of them */
static void rs_elim_neighbours(rs_elim* e, unsigned int i, int sub)
{
	qaws_vec2 const* p = &e->pts[i];
	unsigned int c = rs_elim_cell(e, p);
	int cx = (int)(c % e->gx), cy = (int)(c / e->gx), dx, dy;
	for (dy = -1; dy <= 1; dy++)
		for (dx = -1; dx <= 1; dx++)
		{
			int x = cx + dx, y = cy + dy;
			unsigned int cc, k;
			if (x < 0 || y < 0 || x >= (int)e->gx || y >= (int)e->gy)
				continue;
			cc = (unsigned int)y * e->gx + (unsigned int)x;
			for (k = e->cell_first[cc]; k < e->cell_first[cc + 1]; k++)
			{
				unsigned int j = e->cell_items[k];
				double w;
				if (j == i || e->pos[j] == ~0u)
					continue;
				w = rs_elim_weight(e, hypot(e->pts[j].x - p->x, e->pts[j].y - p->y), i, j);
				if (w == 0)
					continue;
				if (sub)
				{
					/* j loses i's push: its weight falls, so it sinks */
					e->w[j] -= w;
					rs_heap_down(e, e->pos[j]);
				}
				else
					e->w[i] += w;
			}
		}
}

/* count of the ncand candidates kept, the most crowded taken out first */
static qaws_status rs_eliminate(qaws_region_sampler const* s, qaws_vec2 const* cand, unsigned int ncand, unsigned int count, qaws_vec2* out)
{
	unsigned int i, k;
	rs_elim e;
	double hi[2];
	qaws_status st = QAWS_STATUS_OK;
	memset(&e, 0, sizeof(e));
	/* the radius of n points packed hexagonally over the area; the weights'
	   lower bound with beta = 0.65, gamma = 1.5 (the paper's values) */
	e.pts = cand;
	e.n = ncand;
	e.rmax = sqrt(s->area / (2 * sqrt(3.0) * count));
	e.rmin = e.rmax * (1 - pow((double)count / ncand, 1.5)) * 0.65;
	e.cs = 2 * e.rmax;
	/* with a density: radii as 1 / sqrt(rho), rmax at the mean, kept within
	   a factor of 4 of it; the grid sized by the largest */
	if (s->rho)
	{
		rs_batch c;
		double big = 0;
		e.rad = (double*)malloc(sizeof(double) * (ncand + 1));
		if (!e.rad)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		memset(&c, 0, sizeof(c));
		c.s = s;
		c.out = (qaws_vec2*)cand;
		c.rad = e.rad;
		c.rad0 = e.rmax;
		st = rs_batch_run(&c, ncand);
		if (st != QAWS_STATUS_OK)
		{
			free(e.rad);
			return st;
		}
		for (i = 0; i < ncand; i++)
			if (e.rad[i] > big) big = e.rad[i];
		e.cs = 2 * big;
	}
	e.lo[0] = e.lo[1] = 1e300; hi[0] = hi[1] = -1e300;
	for (i = 0; i < ncand; i++)
	{
		if (cand[i].x < e.lo[0]) e.lo[0] = cand[i].x;
		if (cand[i].y < e.lo[1]) e.lo[1] = cand[i].y;
		if (cand[i].x > hi[0]) hi[0] = cand[i].x;
		if (cand[i].y > hi[1]) hi[1] = cand[i].y;
	}
	e.gx = (unsigned int)((hi[0] - e.lo[0]) / e.cs) + 1;
	e.gy = (unsigned int)((hi[1] - e.lo[1]) / e.cs) + 1;
	e.cell_first = (unsigned int*)calloc((size_t)e.gx * e.gy + 2, sizeof(unsigned int));
	e.cell_items = (unsigned int*)malloc(sizeof(unsigned int) * (ncand + 1));
	e.w = (double*)calloc(ncand + 1, sizeof(double));
	e.heap = (unsigned int*)malloc(sizeof(unsigned int) * (ncand + 1));
	e.pos = (unsigned int*)malloc(sizeof(unsigned int) * (ncand + 1));
	if (!e.cell_first || !e.cell_items || !e.w || !e.heap || !e.pos)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
	{
		unsigned int* fill = (unsigned int*)calloc((size_t)e.gx * e.gy + 1, sizeof(unsigned int));
		if (!fill)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			for (i = 0; i < ncand; i++) e.cell_first[rs_elim_cell(&e, &cand[i]) + 1]++;
			for (k = 0; k < e.gx * e.gy; k++) e.cell_first[k + 1] += e.cell_first[k];
			for (i = 0; i < ncand; i++)
			{
				unsigned int c = rs_elim_cell(&e, &cand[i]);
				e.cell_items[e.cell_first[c] + fill[c]++] = i;
			}
			free(fill);
		}
	}
	if (st == QAWS_STATUS_OK)
	{
		/* the weights, the heap, then the heaviest point out until count remain */
		for (i = 0; i < ncand; i++) { e.pos[i] = i; e.heap[i] = i; }
		for (i = 0; i < ncand; i++) rs_elim_neighbours(&e, i, 0);
		e.hn = ncand;
		for (i = ncand / 2; i-- > 0;) rs_heap_down(&e, i);
		while (e.hn > count)
		{
			unsigned int top = e.heap[0];
			rs_heap_swap(&e, 0, e.hn - 1);
			e.hn--;
			e.pos[top] = ~0u;
			rs_heap_down(&e, 0);
			rs_elim_neighbours(&e, top, 1);
		}
		for (i = 0, k = 0; i < ncand; i++)
			if (e.pos[i] != ~0u)
				out[k++] = cand[i];
	}
	free(e.cell_first); free(e.cell_items); free(e.w); free(e.heap); free(e.pos); free(e.rad);
	return st;
}

/* blue noise: 5 n candidates uniform over the region, eliminated to n */
static qaws_status rs_blue_noise(qaws_region_sampler const* s, unsigned int count, unsigned long long seed, qaws_vec2* out)
{
	unsigned int ncand = 5 * count;
	qaws_vec2* cand = (qaws_vec2*)malloc(sizeof(qaws_vec2) * (ncand + 1));
	rs_pattern q;
	qaws_status st = QAWS_STATUS_OK;
	if (!cand)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	rs_pattern_init(&q, QAWS_REGION_SAMPLE_RANDOM, ncand, seed);
	st = rs_map_batch(s, &q, NULL, ncand, cand);
	if (st == QAWS_STATUS_OK)
		st = rs_eliminate(s, cand, ncand, count, out);
	free(cand);
	return st;
}

/* inside: the winding below (x, y) of the pieces crossing the slab of x */
static int rs_inside(qaws_region_sampler const* s, double x, double y)
{
	unsigned int k, j;
	int w = 0;
	if (!(x > s->xs[0] && x < s->xs[s->nxs - 1]) || !(y > s->ylo && y < s->yhi))
		return 0;
	k = rs_slab(s, x);
	for (j = s->first[k]; j < s->first[k + 1]; j++)
	{
		rs_piece const* pc = &s->pieces[s->active[j]];
		double p[2];
		rs_eval(pc->c, rs_param(pc, x), p, NULL);
		if (p[1] < y)
			w += pc->sigma;
	}
	return w > 0;
}

qaws_status qaws_region_sampler_contains(int* out_inside, qaws_region_sampler const* s, qaws_vec2 point)
{
	if (!s || !out_inside)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_inside = rs_inside(s, (double)point.x, (double)point.y);
	return QAWS_STATUS_OK;
}

/* inside, and with a density kept at the rate rho / rhomax (a hash of the
   point's index, so a pattern's point is kept or not whatever the count) */
static int rs_keep(qaws_region_sampler const* s, double x, double y, unsigned int i, unsigned int p)
{
	if (!rs_inside(s, x, y))
		return 0;
	if (!s->rho)
		return 1;
	return rs_hashfloat(i, p ^ 0x5bd1e995u) * s->rhomax < (double)s->rho(s->rho_user, (qaws_scalar)x, (qaws_scalar)y);
}

/* the pattern of total points (of kind and seed) over the box: which are
   kept, in keep (grown to total), and how many */
static qaws_status rs_box_count(qaws_region_sampler const* s, unsigned int kind, unsigned int total, unsigned long long seed,
	double const* box, rs_pattern* q, unsigned char** keep, unsigned int* cap, unsigned int* got)
{
	rs_batch c;
	unsigned int i;
	qaws_status st;
	if (total > *cap)
	{
		unsigned char* k = (unsigned char*)realloc(*keep, total);
		if (!k)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		*keep = k;
		*cap = total;
	}
	rs_pattern_init(q, kind, total, seed);
	memset(&c, 0, sizeof(c));
	c.s = s;
	c.q = q;
	memcpy(c.box, box, sizeof(c.box));
	c.keep = *keep;
	st = rs_batch_run(&c, total);
	*got = 0;
	for (i = 0; i < total; i++)
		*got += (*keep)[i];
	return st;
}

/* the pattern over the region's box, the points outside dropped: the
   pattern's own spacing, the same in every direction. A sequence runs
   until count points fall inside, in blocks; the multi-jitter grid is sized
   by the box's share of the area, and grown when it falls short. */
static qaws_status rs_box(qaws_region_sampler const* s, unsigned int count, unsigned int kind, unsigned long long seed, qaws_vec2* out)
{
	double box[4], ratio;
	unsigned int got = 0, i, cap = 0;
	unsigned char* keep = NULL;
	rs_pattern q;
	qaws_status st = QAWS_STATUS_OK;
	box[0] = s->xs[0];
	box[1] = s->xs[s->nxs - 1] - s->xs[0];
	box[2] = s->ylo;
	box[3] = s->yhi - s->ylo;
	ratio = s->rho ? box[1] * box[3] * s->rhomax / s->mass : box[1] * box[3] / s->area;
	if (!(ratio >= 1) || ratio > 1e6)
		return QAWS_STATUS_DEGENERATE_CURVE;
	if (kind == QAWS_REGION_SAMPLE_STRATIFIED)
	{
		/* every stratum inside kept: dropping strata leaves holes, and
		   eliminating the crowded ones favours the boundary. The grid's size
		   is searched until the strata inside number count, or a few more;
		   those few go at random. */
		unsigned int total = (unsigned int)ceil(count * ratio), best_total = 0, best_got = 0, tries;
		memset(&q, 0, sizeof(q));
		for (tries = 0; tries < 24 && st == QAWS_STATUS_OK; tries++)
		{
			long step;
			st = rs_box_count(s, kind, total, seed, box, &q, &keep, &cap, &got);
			if (got >= count && (best_got == 0 || got < best_got))
			{
				best_got = got;
				best_total = total;
			}
			if (got == count)
				break;
			step = (long)floor(((double)count - got) * ratio + 0.5);
			if (step == 0) step = got < count ? 1 : -1;
			if ((long)total + step < (long)count) step = (long)count - (long)total;
			total = (unsigned int)((long)total + step);
		}
		while (st == QAWS_STATUS_OK && best_got == 0)
		{
			/* never enough: grow until there are */
			total += total / 8 + 1;
			st = rs_box_count(s, kind, total, seed, box, &q, &keep, &cap, &got);
			if (got >= count) { best_got = got; best_total = total; }
		}
		if (st == QAWS_STATUS_OK && total != best_total)
			st = rs_box_count(s, kind, best_total, seed, box, &q, &keep, &cap, &got);
		if (st == QAWS_STATUS_OK)
		{
			/* the surplus at random: a hashed permutation's first count of
			   the kept, in order */
			unsigned int* in = (unsigned int*)malloc(sizeof(unsigned int) * (got + 1));
			if (!in)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			else
			{
				unsigned int k = 0;
				for (i = 0; i < best_total; i++)
					if (keep[i]) in[k++] = i;
				for (i = 0; i < count; i++)
				{
					double u, v;
					rs_pattern_point(&q, in[got == count ? i : rs_permute(i, got, q.p * 0x2545f491u)], &u, &v);
					out[i].x = (qaws_scalar)(box[0] + box[1] * u);
					out[i].y = (qaws_scalar)(box[2] + box[3] * v);
				}
				free(in);
			}
		}
		free(keep);
		return st;
	}
	rs_pattern_init(&q, kind, count, seed);
	for (i = 0; got < count && st == QAWS_STATUS_OK;)
	{
		/* a block of the sequence, its flags in parallel, the kept in order */
		rs_batch c;
		double want = (count - got) * ratio * 1.1 + 256;
		unsigned int n = want > 16777216 ? 16777216u : (unsigned int)want, j;
		if (n > 4000000000u - i)
		{
			st = QAWS_STATUS_DEGENERATE_CURVE;
			break;
		}
		if (n > cap)
		{
			unsigned char* k = (unsigned char*)realloc(keep, n);
			if (!k)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			keep = k;
			cap = n;
		}
		memset(&c, 0, sizeof(c));
		c.s = s;
		c.q = &q;
		c.first = i;
		memcpy(c.box, box, sizeof(c.box));
		c.keep = keep;
		st = rs_batch_run(&c, n);
		for (j = 0; j < n && got < count && st == QAWS_STATUS_OK; j++)
			if (keep[j])
			{
				double u, v;
				rs_pattern_point(&q, i + j, &u, &v);
				out[got].x = (qaws_scalar)(box[0] + box[1] * u);
				out[got].y = (qaws_scalar)(box[2] + box[3] * v);
				got++;
			}
		i += n;
	}
	free(keep);
	return st;
}

qaws_status qaws_region_sampler_generate(qaws_vec2* out, qaws_region_sampler const* s, unsigned int count,
	unsigned int kind, unsigned long long seed)
{
	rs_pattern q;
	unsigned int base = kind & ~QAWS_REGION_SAMPLE_BOX;
	if (!s || (count && !out) || base > QAWS_REGION_SAMPLE_BLUE_NOISE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (count == 0)
		return QAWS_STATUS_OK;
	if (base == QAWS_REGION_SAMPLE_BLUE_NOISE)
		return rs_blue_noise(s, count, seed, out);
	if (kind & QAWS_REGION_SAMPLE_BOX)
		return rs_box(s, count, base, seed, out);
	rs_pattern_init(&q, base, count, seed);
	return rs_map_batch(s, &q, NULL, count, out);
}
