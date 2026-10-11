/*
 * Uniform sampling of curved regions (see qaws_region_sample.h).
 *
 *   1. region: the paths united under the fill rule (qaws_clip), so the
 *      boundary is simple: outer loops counter-clockwise, holes clockwise
 *   2. pieces: every curve cut at its span ends and where dx/dt changes
 *      sign; each piece is x-monotone, sigma = +1 when x grows along it
 *   3. slabs: the pieces' end abscissas sorted and merged; each piece spans
 *      a run of slabs, and A at every slab's left side is summed up
 *   4. map: u * area located among the slabs, x by safeguarded Newton on
 *      A(x) (its derivative is the chord length at x); y in the chords at x,
 *      the pieces' y sorted and paired by the winding
 */

#include "qaws_region_sample.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_path.h"
#include "internal/qaws_internal_flatten.h"
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
};

static qaws_status rs_eval(qaws_curve const* c, double t, double* p, double* d)
{
	qaws_eval_result_2d r;
	qaws_status s = qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | (d ? QAWS_EVAL_FLAG_D1 : 0), &r);
	if (s != QAWS_STATUS_OK)
		return s;
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
	free(s);
}

qaws_status qaws_region_sampler_create(qaws_region_sampler_desc const* desc, qaws_region_sampler** out)
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
	st = qaws_clip_boolean(QAWS_CLIP_UNION, desc->fill_rule, desc->paths, desc->path_count, NULL, 0, &s->region);
	np = st == QAWS_STATUS_OK ? qaws_clip_result_get_path_count(s->region) : 0;
	/* 2. x-monotone pieces */
	for (i = 0; st == QAWS_STATUS_OK && i < np; i++)
	{
		qaws_path_2d p;
		qaws_clip_result_get_path(s->region, i, &p);
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
				double h = 1;
				for (i = 0; i < np; i++)
				{
					qaws_path_2d p;
					qaws_clip_result_get_path(s->region, i, &p);
					if (qaws_path_compute_bounds_2d(&p, &blo, &bhi) == QAWS_STATUS_OK && (double)(bhi.y - blo.y) > h)
						h = (double)(bhi.y - blo.y);
				}
				s->tol = desc->tolerance > 0 ? (double)desc->tolerance : RS_AREA_REL * (s->ext > 0 ? s->ext : 1) * h;
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

qaws_status qaws_region_sampler_area_left(qaws_region_sampler const* s, qaws_scalar x, qaws_scalar* out)
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

qaws_status qaws_region_sampler_map(qaws_region_sampler const* s, qaws_scalar u, qaws_scalar v, qaws_vec2* out)
{
	double target, lo, hi, x, fx, chord, ytot = 0, yt, w = 0;
	unsigned int k, it, j, nc;
	rs_cross stackc[64], *cr = stackc;
	if (!s || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	target = (double)u * s->area;
	if (target < 0) target = 0;
	if (target > s->area) target = s->area;
	/* the slab: the last abscissa with A below the target */
	{
		unsigned int a = 0, b = s->nxs - 2;
		while (a < b)
		{
			unsigned int m = (a + b + 1) / 2;
			if (s->acum[m] <= target) a = m; else b = m - 1;
		}
		k = a;
	}
	/* x: Newton on A(x) - target, its derivative the chord, in the slab */
	lo = s->xs[k]; hi = s->xs[k + 1];
	x = s->acum[k + 1] > s->acum[k] ? lo + (hi - lo) * (target - s->acum[k]) / (s->acum[k + 1] - s->acum[k]) : 0.5 * (lo + hi);
	for (it = 0; it < 60; it++)
	{
		double nx;
		fx = rs_area_in(s, k, x, &chord) - target;
		if (fabs(fx) <= s->tol)
			break;
		if (fx > 0) hi = x; else lo = x;
		nx = chord > 0 ? x - fx / chord : 0.5 * (lo + hi);
		if (!(nx > lo && nx < hi))
			nx = 0.5 * (lo + hi);
		if (hi - lo <= RS_MERGE_REL * 1e-3 * s->ext)
			break;
		x = nx;
	}
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

qaws_status qaws_region_sampler_map_many(qaws_region_sampler const* s, qaws_vec2 const* uv, unsigned int count, qaws_vec2* out)
{
	unsigned int i;
	qaws_status st = QAWS_STATUS_OK;
	if (!s || (count && (!uv || !out)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count && st == QAWS_STATUS_OK; i++)
		st = qaws_region_sampler_map(s, uv[i].x, uv[i].y, &out[i]);
	return st;
}

/* splitmix64 */
static unsigned long long rs_next(unsigned long long* st)
{
	unsigned long long z = (*st += 0x9E3779B97F4A7C15ull);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

static double rs_unit(unsigned long long r)
{
	return (double)(r >> 11) * (1.0 / 9007199254740992.0);
}

qaws_status qaws_region_sampler_generate(qaws_region_sampler const* s, unsigned int count, unsigned int kind,
	unsigned long long seed, qaws_vec2* out)
{
	unsigned long long st = seed;
	unsigned int i;
	qaws_status r = QAWS_STATUS_OK;
	if (!s || (count && !out) || kind > QAWS_REGION_SAMPLE_SOBOL)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (kind == QAWS_REGION_SAMPLE_RANDOM)
	{
		for (i = 0; i < count && r == QAWS_STATUS_OK; i++)
		{
			double u = rs_unit(rs_next(&st)), v = rs_unit(rs_next(&st));
			r = qaws_region_sampler_map(s, (qaws_scalar)u, (qaws_scalar)v, &out[i]);
		}
		return r;
	}
	/* Sobol in 2D: van der Corput in base 2, then the direction numbers
	   v_j = v_(j-1) ^ (v_(j-1) >> 1); both shifted digitally by the seed */
	{
		unsigned int shift0 = (unsigned int)(rs_next(&st) >> 32), shift1 = (unsigned int)(rs_next(&st) >> 32);
		unsigned int dir[32], b;
		dir[0] = 0x80000000u;
		for (b = 1; b < 32; b++) dir[b] = dir[b - 1] ^ (dir[b - 1] >> 1);
		if (seed == 0) shift0 = shift1 = 0;
		for (i = 0; i < count && r == QAWS_STATUS_OK; i++)
		{
			unsigned int x0 = 0, x1 = 0, n = i;
			for (b = 0; n; b++, n >>= 1)
				if (n & 1u) { x0 ^= 0x80000000u >> b; x1 ^= dir[b]; }
			x0 ^= shift0; x1 ^= shift1;
			r = qaws_region_sampler_map(s, (qaws_scalar)((x0 + 0.5) / 4294967296.0), (qaws_scalar)((x1 + 0.5) / 4294967296.0), &out[i]);
		}
	}
	return r;
}
