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

static void rs_pattern_point(rs_pattern* q, unsigned int i, double* u, double* v)
{
	switch (q->kind)
	{
	case QAWS_REGION_SAMPLE_SOBOL:
	{
		unsigned int x0 = 0, x1 = 0, k = i, b;
		for (b = 0; k; b++, k >>= 1)
			if (k & 1u) { x0 ^= 0x80000000u >> b; x1 ^= q->dir[b]; }
		x0 ^= q->shift0; x1 ^= q->shift1;
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
		*u = rs_unit(rs_next(&q->state));
		*v = rs_unit(rs_next(&q->state));
		break;
	}
	if (*u >= 1) *u = 1 - 1e-16;
	if (*v >= 1) *v = 1 - 1e-16;
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

static double rs_elim_weight(rs_elim const* e, double d)
{
	double t;
	if (d < e->rmin) d = e->rmin;
	if (d >= 2 * e->rmax) return 0;
	t = 1 - d / (2 * e->rmax);
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
				w = rs_elim_weight(e, hypot(e->pts[j].x - p->x, e->pts[j].y - p->y));
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

static qaws_status rs_blue_noise(qaws_region_sampler const* s, unsigned int count, unsigned long long seed, qaws_vec2* out)
{
	unsigned int ncand = 5 * count, i, k;
	qaws_vec2* cand = (qaws_vec2*)malloc(sizeof(qaws_vec2) * (ncand + 1));
	rs_elim e;
	rs_pattern q;
	double hi[2];
	qaws_status st = QAWS_STATUS_OK;
	memset(&e, 0, sizeof(e));
	if (!cand)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* 5 n candidates, uniform over the region */
	rs_pattern_init(&q, QAWS_REGION_SAMPLE_RANDOM, ncand, seed);
	for (i = 0; i < ncand && st == QAWS_STATUS_OK; i++)
	{
		double u, v;
		rs_pattern_point(&q, i, &u, &v);
		st = qaws_region_sampler_map(s, (qaws_scalar)u, (qaws_scalar)v, &cand[i]);
	}
	if (st != QAWS_STATUS_OK) { free(cand); return st; }
	/* the radius of n points packed hexagonally over the area; the weights'
	   lower bound with beta = 0.65, gamma = 1.5 (the paper's values) */
	e.pts = cand;
	e.n = ncand;
	e.rmax = sqrt(s->area / (2 * sqrt(3.0) * count));
	e.rmin = e.rmax * (1 - pow((double)count / ncand, 1.5)) * 0.65;
	e.cs = 2 * e.rmax;
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
	free(e.cell_first); free(e.cell_items); free(e.w); free(e.heap); free(e.pos);
	free(cand);
	return st;
}

qaws_status qaws_region_sampler_generate(qaws_region_sampler const* s, unsigned int count, unsigned int kind,
	unsigned long long seed, qaws_vec2* out)
{
	rs_pattern q;
	unsigned int i;
	qaws_status r = QAWS_STATUS_OK;
	if (!s || (count && !out) || kind > QAWS_REGION_SAMPLE_BLUE_NOISE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (count == 0)
		return QAWS_STATUS_OK;
	if (kind == QAWS_REGION_SAMPLE_BLUE_NOISE)
		return rs_blue_noise(s, count, seed, out);
	rs_pattern_init(&q, kind, count, seed);
	for (i = 0; i < count && r == QAWS_STATUS_OK; i++)
	{
		double u, v;
		rs_pattern_point(&q, i, &u, &v);
		r = qaws_region_sampler_map(s, (qaws_scalar)u, (qaws_scalar)v, &out[i]);
	}
	return r;
}
