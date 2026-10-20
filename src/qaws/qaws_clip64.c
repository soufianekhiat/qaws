/*
 * Exact Boolean operations on int64 polygons: the arrangement of qaws_clip
 * with exact integer predicates.
 *
 *   1. segments of every path (closed ones closed); an x sweep pairs them
 *   2. hits: a crossing is the rational point p + r tn / den (den, tn: exact
 *      128-bit cross products); collinear segments cut each other at the
 *      end points lying inside the other
 *   3. vertices: points sorted exactly (rational comparisons), equal ones
 *      merged; no tolerance anywhere
 *   4. edges: each segment cut at its points in parameter order; straight
 *      edges between the same two vertices are the same edge, merged with
 *      their winding steps summed; those that cancel dropped
 *   5. graph: half-edges sorted round each vertex by direction (half plane,
 *      then the exact sign of the cross product of integer directions)
 *   6. faces: a cycle turning left at its lowest vertex bounds a face; every
 *      other cycle is a part's outside and belongs to the face of another
 *      part that winds round it (exact crossing count)
 *   7. windings, the clip type on the fill rule, loops cut at touch points,
 *      nesting from face components (as qaws_clip)
 *   8. vertices rounded to the nearest integer; repeated and (unless
 *      preserved) collinear points of the rounded path dropped
 */

#include "qaws_clip64.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_wide.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_clip.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define C64_NONE 0xFFFFFFFFu

/* ======================================================================== */
/*  Data                                                                    */
/* ======================================================================== */

typedef struct c64_seg
{
	int64_t ax, ay, bx, by;
	int64_t lox, hix, loy, hiy;
	unsigned int operand, path;   /* operand 0 subject, 1 clip, 2 open subject */
} c64_seg;

typedef struct c64_pt                /* the rational point (x / d, y / d), d > 0 */
{
	qaws_wide x, y, d;
	double fx, fy;                   /* nearest doubles, for fast comparisons */
} c64_pt;

typedef struct c64_cut               /* a point on a segment at t = tn / td */
{
	qaws_wide tn, td;
	double ft;
	unsigned int seg, pt;
} c64_cut;

typedef struct c64_edge
{
	unsigned int va, vb;
	int64_t dx, dy;                  /* direction va -> vb (the segment's) */
	int ds, dc;
	unsigned int open_path;          /* open subjects: the path; else C64_NONE */
	int alive;
} c64_edge;

typedef struct c64_half
{
	unsigned int next, face, pos;
	int out;
} c64_half;

struct qaws_clip64_result
{
	struct c64_out
	{
		unsigned int first, count;   /* into pts (rounded) */
		unsigned int efirst, ecount; /* into ex (exact, before rounding) */
		unsigned int parent, depth;
		int hole;
	} *paths, *open_paths;
	unsigned int path_count, path_cap, open_count, open_cap;
	int64_t* pts;
	double* exact;
	double* ex;
	unsigned int nex, capex;
	unsigned int npts, cappts;
};

typedef struct c64_ctx
{
	qaws_clip64_desc const* d;
	c64_seg* seg;
	unsigned int nseg, capseg;
	c64_pt* pt;
	unsigned int npt, cappt;
	c64_cut* cut;
	unsigned int ncut, capcut;
	unsigned int* vid;               /* point -> vertex */
	unsigned int nv;
	unsigned int* vrep;              /* vertex -> a representative point */
	double* vf;                      /* vertex -> its point as doubles (x, y), compact */
	c64_edge* e;
	unsigned int ne, cape;
	c64_edge* oe;                    /* open-subject edges */
	unsigned int noe, capoe;
	c64_half* h;
	unsigned int* vfirst;            /* vertex -> first outgoing in vout */
	unsigned int* vcount;
	unsigned int* vout;
	int64_t* vround;                 /* vertex rounded */
	qaws_clip64_result* r;
} c64_ctx;

#define C64_GROW(ptr, n, cap, type) \
	do { if ((n) >= 0x7FFFFFFFu) { return QAWS_STATUS_ALLOCATION_FAILURE; } /* 2 n half-edges count in unsigned int */ \
		if ((n) == (cap)) { unsigned int nc_ = (cap) ? 2 * (cap) : 64; type* np_ = (type*)qaws_internal_realloc((ptr), qaws_internal_mul2(sizeof(type), nc_)); \
		if (!np_) { return QAWS_STATUS_ALLOCATION_FAILURE; } (ptr) = np_; (cap) = nc_; } } while (0)

/* ======================================================================== */
/*  Rational helpers                                                        */
/* ======================================================================== */

static void c64_pt_float(c64_pt* p)
{
	double d = qaws_wide_to_double(&p->d);
	p->fx = qaws_wide_to_double(&p->x) / d;
	p->fy = qaws_wide_to_double(&p->y) / d;
}

/* sign of a.x - b.x (or y): a.x b.d - b.x a.d */
static int c64_cmp_coord(qaws_wide const* an, qaws_wide const* ad, qaws_wide const* bn, qaws_wide const* bd)
{
	qaws_wide u = qaws_wide_mul(an, bd), v = qaws_wide_mul(bn, ad);
	return qaws_wide_cmp(&u, &v);
}

static int c64_cmp_pt(c64_pt const* a, c64_pt const* b)
{
	double tx = 1e-9 * (fabs(a->fx) + fabs(b->fx) + 1), ty = 1e-9 * (fabs(a->fy) + fabs(b->fy) + 1);
	int c;
	if (a->fy < b->fy - ty) return -1;
	if (a->fy > b->fy + ty) return 1;
	c = c64_cmp_coord(&a->y, &a->d, &b->y, &b->d);
	if (c) return c;
	if (a->fx < b->fx - tx) return -1;
	if (a->fx > b->fx + tx) return 1;
	return c64_cmp_coord(&a->x, &a->d, &b->x, &b->d);
}

static int c64_cmp_t(c64_cut const* a, c64_cut const* b)
{
	if (a->ft < b->ft - 1e-9) return -1;
	if (a->ft > b->ft + 1e-9) return 1;
	return c64_cmp_coord(&a->tn, &a->td, &b->tn, &b->td);
}

/* sign of ax by - ay bx for int64 components: in double when the result is
   farther from zero than its rounding error, exactly otherwise */
static int c64_cross_sign(int64_t ax, int64_t ay, int64_t bx, int64_t by)
{
	double p = (double)ax * (double)by, q = (double)ay * (double)bx, d = p - q;
	qaws_wide c;
	if (fabs(d) > 8e-16 * (fabs(p) + fabs(q)) + 1e-300)
		return d > 0 ? 1 : -1;
	c = qaws_wide_cross(ax, ay, bx, by);
	return qaws_wide_sign(&c);
}

/* compact sort keys: two doubles (compared within a tolerance first) and the item */
typedef struct c64_skey
{
	double a, b;
	unsigned int i;
} c64_skey;

typedef int (*c64_skey_fn)(void const* ctx, c64_skey const* p, c64_skey const* q);

/* merge sort of keys (insertion below 16): contiguous, so the comparisons stay in cache */
static void c64_sort_keys(c64_skey* k, c64_skey* tmp, unsigned int n, c64_skey_fn cmp, void const* ctx)
{
	unsigned int w, i;
	for (i = 0; i < n; i += 16)
	{
		unsigned int e = i + 16 < n ? i + 16 : n, a, b;
		for (a = i + 1; a < e; a++)
		{
			c64_skey key = k[a];
			for (b = a; b > i && cmp(ctx, &k[b - 1], &key) > 0; b--)
				k[b] = k[b - 1];
			k[b] = key;
		}
	}
	/* (the merge passes bound runs by differences and cap w at n, so their indices do not wrap) */
	for (w = 16; w < n; w = w <= (n - 1) / 2 ? 2 * w : n)
	{
		for (i = 0; i < n;)
		{
			unsigned int a = i, am = w < n - i ? i + w : n, b = am, bm = w < n - am ? am + w : n, j = i;
			if (am >= bm || cmp(ctx, &k[am - 1], &k[am]) <= 0)
				memcpy(&tmp[i], &k[i], sizeof(c64_skey) * (bm - i));
			else
			{
				while (a < am && b < bm)
					tmp[j++] = cmp(ctx, &k[a], &k[b]) <= 0 ? k[a++] : k[b++];
				while (a < am) tmp[j++] = k[a++];
				while (b < bm) tmp[j++] = k[b++];
			}
			i = bm;
		}
		memcpy(k, tmp, sizeof(c64_skey) * n);
	}
}

static int c64_skey_pt(void const* ctx, c64_skey const* p, c64_skey const* q)
{
	c64_ctx const* x = (c64_ctx const*)ctx;
	double ty = 1e-9 * (fabs(p->a) + fabs(q->a) + 1);
	if (p->a < q->a - ty) return -1;
	if (p->a > q->a + ty) return 1;
	return c64_cmp_pt(&x->pt[p->i], &x->pt[q->i]);
}

static int c64_skey_t(void const* ctx, c64_skey const* p, c64_skey const* q)
{
	c64_ctx const* x = (c64_ctx const*)ctx;
	if (p->a < q->a - 1e-9) return -1;
	if (p->a > q->a + 1e-9) return 1;
	return c64_cmp_t(&x->cut[p->i], &x->cut[q->i]);
}

/* the points sorted (y, then x), as indices */
static qaws_status c64_sort_points(c64_ctx* x, unsigned int* idx)
{
	c64_skey* k = (c64_skey*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(c64_skey), 2, (x->npt + 1)));
	unsigned int i;
	if (!k)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->npt; i++)
	{
		k[i].a = x->pt[i].fy;
		k[i].b = x->pt[i].fx;
		k[i].i = i;
	}
	c64_sort_keys(k, k + x->npt, x->npt, c64_skey_pt, x);
	for (i = 0; i < x->npt; i++)
		idx[i] = k[i].i;
	qaws_internal_free(k);
	return QAWS_STATUS_OK;
}

/* the cuts by segment (counting sort), then along each segment by t */
static qaws_status c64_sort_cuts(c64_ctx* x, unsigned int* idx)
{
	unsigned int* start = (unsigned int*)qaws_internal_scratch_zalloc(x->nseg + 2, sizeof(unsigned int));
	c64_skey* k = (c64_skey*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(c64_skey), 2, (x->ncut + 1)));
	unsigned int i, s;
	if (!start || !k)
	{
		qaws_internal_free(start);
		qaws_internal_free(k);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < x->ncut; i++)
		start[x->cut[i].seg + 1]++;
	for (s = 0; s < x->nseg; s++)
		start[s + 1] += start[s];
	for (i = 0; i < x->ncut; i++)
	{
		unsigned int p = start[x->cut[i].seg]++;
		k[p].a = x->cut[i].ft;
		k[p].b = 0;
		k[p].i = i;
	}
	/* start[s] is now the end of segment s's run */
	for (s = 0, i = 0; s < x->nseg; s++)
	{
		unsigned int e = start[s];
		if (e - i > 1)
			c64_sort_keys(k + i, k + x->ncut + i, e - i, c64_skey_t, x);
		i = e;
	}
	for (i = 0; i < x->ncut; i++)
		idx[i] = k[i].i;
	qaws_internal_free(start);
	qaws_internal_free(k);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  1-2. Segments and hits                                                  */
/* ======================================================================== */

static qaws_status c64_add_point(c64_ctx* x, int64_t px, int64_t py, unsigned int* out)
{
	c64_pt* p;
	C64_GROW(x->pt, x->npt, x->cappt, c64_pt);
	p = &x->pt[x->npt];
	p->x = qaws_wide_from_i64(px);
	p->y = qaws_wide_from_i64(py);
	p->d = qaws_wide_from_i64(1);
	p->fx = (double)px;
	p->fy = (double)py;
	*out = x->npt++;
	return QAWS_STATUS_OK;
}

static qaws_status c64_add_cut(c64_ctx* x, unsigned int seg, qaws_wide const* tn, qaws_wide const* td, unsigned int pt)
{
	c64_cut* c;
	C64_GROW(x->cut, x->ncut, x->capcut, c64_cut);
	c = &x->cut[x->ncut++];
	c->tn = *tn;
	c->td = *td;
	c->ft = qaws_wide_to_double(tn) / qaws_wide_to_double(td);
	c->seg = seg;
	c->pt = pt;
	return QAWS_STATUS_OK;
}

static int c64_ok_coord(int64_t v)
{
	return v <= QAWS_CLIP64_MAX_COORD && v >= -QAWS_CLIP64_MAX_COORD;
}

static qaws_status c64_add_paths(c64_ctx* x, qaws_path64 const* paths, unsigned int count, unsigned int operand)
{
	unsigned int i, k;
	for (i = 0; i < count; i++)
	{
		qaws_path64 const* p = &paths[i];
		unsigned int n = p->point_count, first = x->nseg;
		if (n && !p->points)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (k = 0; k < 2 * n; k++)
			if (!c64_ok_coord(p->points[k]))
				return QAWS_STATUS_OUT_OF_RANGE;
		for (k = 0; k + (operand == 2 ? 1u : 0u) < n; k++)
		{
			unsigned int j = (k + 1) % n;
			c64_seg* s;
			int64_t ax = p->points[2 * k], ay = p->points[2 * k + 1], bx = p->points[2 * j], by = p->points[2 * j + 1];
			if (ax == bx && ay == by)
				continue;
			C64_GROW(x->seg, x->nseg, x->capseg, c64_seg);
			s = &x->seg[x->nseg++];
			s->ax = ax; s->ay = ay; s->bx = bx; s->by = by;
			s->lox = ax < bx ? ax : bx; s->hix = ax < bx ? bx : ax;
			s->loy = ay < by ? ay : by; s->hiy = ay < by ? by : ay;
			s->operand = operand;
			s->path = i;
		}
		/* a region path needs three corners to enclose anything */
		if (operand < 2 && x->nseg - first < 2)
			x->nseg = first;
	}
	return QAWS_STATUS_OK;
}

/* the hits of segments a and b (a < b) */
static qaws_status c64_pair(c64_ctx* x, unsigned int ia, unsigned int ib)
{
	c64_seg const* A = &x->seg[ia];
	c64_seg const* B = &x->seg[ib];
	int64_t rx = A->bx - A->ax, ry = A->by - A->ay, sx = B->bx - B->ax, sy = B->by - B->ay;
	int64_t wx = B->ax - A->ax, wy = B->ay - A->ay;
	qaws_wide den, zero, one;
	qaws_status st;
	if (A->operand == 2 && B->operand == 2)
		return QAWS_STATUS_OK;
	/* filter in double: the products are exact to a relative 2^-52 each, so
	   a cross product off zero by more than its error bound has its sign;
	   pairs that clearly miss are dropped without the exact arithmetic */
	{
		double drx = (double)rx, dry = (double)ry, dsx = (double)sx, dsy = (double)sy, dwx = (double)wx, dwy = (double)wy;
		double fden = drx * dsy - dry * dsx, eden = 4e-16 * (fabs(drx * dsy) + fabs(dry * dsx));
		double ftn = dwx * dsy - dwy * dsx, etn = 4e-16 * (fabs(dwx * dsy) + fabs(dwy * dsx));
		double fun = dwx * dry - dwy * drx, eun = 4e-16 * (fabs(dwx * dry) + fabs(dwy * drx));
		if (fabs(fden) > eden)
		{
			/* t = tn / den and u = un / den must both lie in [0, 1] */
			double sden = fden > 0 ? 1.0 : -1.0, tn = ftn * sden, un = fun * sden, ad = fabs(fden);
			if (tn < -etn || un < -eun || tn - ad > etn + eden || un - ad > eun + eden)
				return QAWS_STATUS_OK;
		}
		else if (fabs(fun) > eun)
			return QAWS_STATUS_OK;   /* parallel, not collinear */
	}
	den = qaws_wide_cross(rx, ry, sx, sy);
	zero = qaws_wide_from_i64(0);
	one = qaws_wide_from_i64(1);
	if (qaws_wide_sign(&den) != 0)
	{
		qaws_wide tn = qaws_wide_cross(wx, wy, sx, sy), un = qaws_wide_cross(wx, wy, rx, ry);
		int at0, at1, bt0, bt1;
		unsigned int p;
		if (qaws_wide_sign(&den) < 0)
		{
			den = qaws_wide_neg(den);
			tn = qaws_wide_neg(tn);
			un = qaws_wide_neg(un);
		}
		if (qaws_wide_sign(&tn) < 0 || qaws_wide_cmp(&tn, &den) > 0 || qaws_wide_sign(&un) < 0 || qaws_wide_cmp(&un, &den) > 0)
			return QAWS_STATUS_OK;
		at0 = qaws_wide_sign(&tn) == 0; at1 = qaws_wide_cmp(&tn, &den) == 0;
		bt0 = qaws_wide_sign(&un) == 0; bt1 = qaws_wide_cmp(&un, &den) == 0;
		/* both at their ends: the end points are points already */
		if ((at0 || at1) && (bt0 || bt1))
			return QAWS_STATUS_OK;
		if (at0 || at1 || bt0 || bt1)
		{
			/* an end point on the other segment's inside: an integer point */
			int64_t px = at0 ? A->ax : at1 ? A->bx : bt0 ? B->ax : B->bx;
			int64_t py = at0 ? A->ay : at1 ? A->by : bt0 ? B->ay : B->by;
			st = c64_add_point(x, px, py, &p);
		}
		else
		{
			/* (a den + r tn) / den */
			c64_pt* q;
			qaws_wide ax = qaws_wide_from_i64(A->ax), ay = qaws_wide_from_i64(A->ay);
			qaws_wide wrx = qaws_wide_from_i64(rx), wry = qaws_wide_from_i64(ry), t1, t2;
			C64_GROW(x->pt, x->npt, x->cappt, c64_pt);
			q = &x->pt[x->npt];
			t1 = qaws_wide_mul(&ax, &den); t2 = qaws_wide_mul(&wrx, &tn);
			q->x = qaws_wide_add(t1, &t2);
			t1 = qaws_wide_mul(&ay, &den); t2 = qaws_wide_mul(&wry, &tn);
			q->y = qaws_wide_add(t1, &t2);
			q->d = den;
			c64_pt_float(q);
			p = x->npt++;
			st = QAWS_STATUS_OK;
		}
		if (st == QAWS_STATUS_OK) st = c64_add_cut(x, ia, &tn, &den, p);
		if (st == QAWS_STATUS_OK) st = c64_add_cut(x, ib, &un, &den, p);
		return st;
	}
	/* parallel: collinear segments cut each other at the end points inside */
	{
		qaws_wide c = qaws_wide_cross(wx, wy, rx, ry);
		qaws_wide la, lb;
		int k;
		if (qaws_wide_sign(&c) != 0)
			return QAWS_STATUS_OK;
		la = qaws_wide_dot(rx, ry, rx, ry);
		lb = qaws_wide_dot(sx, sy, sx, sy);
		for (k = 0; k < 4; k++)
		{
			/* k 0, 1: B's ends on A; k 2, 3: A's ends on B */
			c64_seg const* S = k < 2 ? B : A;
			c64_seg const* O = k < 2 ? A : B;
			int64_t px = (k & 1) ? S->bx : S->ax, py = (k & 1) ? S->by : S->ay;
			int64_t ox = O->bx - O->ax, oy = O->by - O->ay;
			qaws_wide dt = qaws_wide_dot(px - O->ax, py - O->ay, ox, oy);
			qaws_wide const* L = k < 2 ? &la : &lb;
			unsigned int p;
			qaws_wide send = (k & 1) ? one : zero;
			if (qaws_wide_sign(&dt) <= 0 || qaws_wide_cmp(&dt, L) >= 0)
				continue;
			st = c64_add_point(x, px, py, &p);
			if (st == QAWS_STATUS_OK) st = c64_add_cut(x, k < 2 ? ia : ib, &dt, L, p);
			if (st == QAWS_STATUS_OK) st = c64_add_cut(x, k < 2 ? ib : ia, &send, &one, p);
			if (st != QAWS_STATUS_OK)
				return st;
		}
	}
	return QAWS_STATUS_OK;
}

static double c64_down(int64_t v) { double d = (double)v; return d - fabs(d) * 4e-16 - 1e-300; }
static double c64_up(int64_t v) { double d = (double)v; return d + fabs(d) * 4e-16 + 1e-300; }

/* might segments ia and ib meet? The double filter of c64_pair alone (0: proven apart) */
static int c64_maybe(c64_ctx const* x, unsigned int ia, unsigned int ib)
{
	c64_seg const* A = &x->seg[ia];
	c64_seg const* B = &x->seg[ib];
	int64_t rx = A->bx - A->ax, ry = A->by - A->ay, sx = B->bx - B->ax, sy = B->by - B->ay;
	int64_t wx = B->ax - A->ax, wy = B->ay - A->ay;
	double drx = (double)rx, dry = (double)ry, dsx = (double)sx, dsy = (double)sy, dwx = (double)wx, dwy = (double)wy;
	double fden = drx * dsy - dry * dsx, eden = 4e-16 * (fabs(drx * dsy) + fabs(dry * dsx));
	double ftn = dwx * dsy - dwy * dsx, etn = 4e-16 * (fabs(dwx * dsy) + fabs(dwy * dsx));
	double fun = dwx * dry - dwy * drx, eun = 4e-16 * (fabs(dwx * dry) + fabs(dwy * drx));
	if (A->operand == 2 && B->operand == 2)
		return 0;
	if (B->lox > A->hix || B->hix < A->lox || B->loy > A->hiy || B->hiy < A->loy)
		return 0;
	if (fabs(fden) > eden)
	{
		double sden = fden > 0 ? 1.0 : -1.0, tn = ftn * sden, un = fun * sden, ad = fabs(fden);
		return !(tn < -etn || un < -eun || tn - ad > etn + eden || un - ad > eun + eden);
	}
	return !(fabs(fun) > eun);
}

typedef struct c64_pairs
{
	c64_ctx* x;
	unsigned int* p;
	unsigned int n, cap;
} c64_pairs;

static qaws_status c64_visit(void* user, unsigned int i, unsigned int j)
{
	c64_pairs* P = (c64_pairs*)user;
	if (!c64_maybe(P->x, i, j))
		return QAWS_STATUS_OK;
	if (P->n + 2 > P->cap)
	{
		unsigned int nc = P->cap ? 2 * P->cap : 4096;
		unsigned int* g = (unsigned int*)qaws_internal_realloc(P->p, qaws_internal_mul2(sizeof(unsigned int), nc));
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		P->p = g;
		P->cap = nc;
	}
	P->p[P->n++] = i;
	P->p[P->n++] = j;
	return QAWS_STATUS_OK;
}

/* room for n more cuts and m more points, once */
static qaws_status c64_reserve(c64_ctx* x, unsigned int n, unsigned int m)
{
	if (x->ncut + n > x->capcut)
	{
		c64_cut* g = (c64_cut*)qaws_internal_realloc(x->cut, qaws_internal_mul2(sizeof(c64_cut), (x->ncut + n)));
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->cut = g;
		x->capcut = x->ncut + n;
	}
	if (x->npt + m > x->cappt)
	{
		c64_pt* g = (c64_pt*)qaws_internal_realloc(x->pt, qaws_internal_mul2(sizeof(c64_pt), (x->npt + m)));
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->pt = g;
		x->cappt = x->npt + m;
	}
	return QAWS_STATUS_OK;
}

static qaws_status c64_hits(c64_ctx* x)
{
	qaws_bp_box* boxes = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (x->nseg + 1)));
	unsigned int i;
	c64_pairs P;
	qaws_status s;
	qaws_wide zero = qaws_wide_from_i64(0), one = qaws_wide_from_i64(1);
	if (!boxes)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* every segment's end points */
	s = c64_reserve(x, 2 * x->nseg, 2 * x->nseg);
	for (i = 0; i < x->nseg && s == QAWS_STATUS_OK; i++)
	{
		unsigned int p;
		s = c64_add_point(x, x->seg[i].ax, x->seg[i].ay, &p);
		if (s == QAWS_STATUS_OK) s = c64_add_cut(x, i, &zero, &one, p);
		if (s == QAWS_STATUS_OK) s = c64_add_point(x, x->seg[i].bx, x->seg[i].by, &p);
		if (s == QAWS_STATUS_OK) s = c64_add_cut(x, i, &one, &one, p);
		boxes[i].lo[0] = c64_down(x->seg[i].lox); boxes[i].hi[0] = c64_up(x->seg[i].hix);
		boxes[i].lo[1] = c64_down(x->seg[i].loy); boxes[i].hi[1] = c64_up(x->seg[i].hiy);
		boxes[i].lo[2] = boxes[i].hi[2] = 0;
	}
	/* the pairs that may meet (boxes, then the double filter), then exactly with the room reserved */
	memset(&P, 0, sizeof(P));
	P.x = x;
	if (s == QAWS_STATUS_OK && x->nseg > 1)
		s = qaws_internal_broadphase(NULL, boxes, x->nseg, 2, NULL, c64_visit, &P);
	qaws_internal_free(boxes);
	if (s == QAWS_STATUS_OK)
		s = c64_reserve(x, P.n, P.n / 2);
	for (i = 0; i < P.n && s == QAWS_STATUS_OK; i += 2)
		s = c64_pair(x, P.p[i], P.p[i + 1]);
	qaws_internal_free(P.p);
	return s;
}

/* ======================================================================== */
/*  3. Vertices                                                             */
/* ======================================================================== */

static qaws_status c64_vertices(c64_ctx* x)
{
	unsigned int* idx = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (x->npt + 1)));
	unsigned int i;
	if (!idx)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	x->vid = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->npt + 1)));
	x->vrep = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->npt + 1)));
	if (!x->vid || !x->vrep)
	{
		qaws_internal_free(idx);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (c64_sort_points(x, idx) != QAWS_STATUS_OK)
	{
		qaws_internal_free(idx);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	x->nv = 0;
	for (i = 0; i < x->npt; i++)
	{
		if (i == 0 || c64_cmp_pt(&x->pt[idx[i - 1]], &x->pt[idx[i]]) != 0)
			x->vrep[x->nv++] = idx[i];
		x->vid[idx[i]] = x->nv - 1;
	}
	qaws_internal_free(idx);
	x->vf = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (x->nv + 1)));
	if (!x->vf)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->nv; i++)
	{
		x->vf[2 * i] = x->pt[x->vrep[i]].fx;
		x->vf[2 * i + 1] = x->pt[x->vrep[i]].fy;
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  4. Edges                                                                */
/* ======================================================================== */

typedef struct c64_key { unsigned int lo, hi, e; } c64_key;

static int c64_cmp_key(void const* a, void const* b)
{
	c64_key const* p = (c64_key const*)a;
	c64_key const* q = (c64_key const*)b;
	if (p->lo != q->lo) return p->lo < q->lo ? -1 : 1;
	if (p->hi != q->hi) return p->hi < q->hi ? -1 : 1;
	return p->e < q->e ? -1 : (p->e > q->e ? 1 : 0);
}

/* room for one more edge (C64_GROW as a call, for callers that free their own blocks on failure) */
static qaws_status c64_edge_room(c64_edge** e, unsigned int n, unsigned int* cap)
{
	C64_GROW(*e, n, *cap, c64_edge);
	return QAWS_STATUS_OK;
}

static qaws_status c64_edges(c64_ctx* x)
{
	unsigned int* idx = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (x->ncut + 1)));
	unsigned int i, j, n;
	c64_key* k;
	if (!idx)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (c64_sort_cuts(x, idx) != QAWS_STATUS_OK)
	{
		qaws_internal_free(idx);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < x->ncut; i = j)
	{
		unsigned int sg = x->cut[idx[i]].seg, prev = C64_NONE;
		c64_seg const* S = &x->seg[sg];
		for (j = i; j < x->ncut && x->cut[idx[j]].seg == sg; j++)
		{
			unsigned int v = x->vid[x->cut[idx[j]].pt];
			if (prev != C64_NONE && v != prev)
			{
				c64_edge e;
				e.va = prev; e.vb = v;
				e.dx = S->bx - S->ax; e.dy = S->by - S->ay;
				e.ds = S->operand == 0 ? 1 : 0;
				e.dc = S->operand == 1 ? 1 : 0;
				e.open_path = S->operand == 2 ? S->path : C64_NONE;
				e.alive = 1;
				if (S->operand == 2)
				{
					if (c64_edge_room(&x->oe, x->noe, &x->capoe) != QAWS_STATUS_OK)
					{
						qaws_internal_free(idx);
						return QAWS_STATUS_ALLOCATION_FAILURE;
					}
					x->oe[x->noe++] = e;
				}
				else
				{
					if (c64_edge_room(&x->e, x->ne, &x->cape) != QAWS_STATUS_OK)
					{
						qaws_internal_free(idx);
						return QAWS_STATUS_ALLOCATION_FAILURE;
					}
					x->e[x->ne++] = e;
				}
			}
			prev = v;
		}
	}
	qaws_internal_free(idx);
	/* straight edges between the same vertices are the same edge */
	k = (c64_key*)qaws_internal_malloc(qaws_internal_mul2(sizeof(c64_key), (x->ne + 1)));
	if (!k)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->ne; i++)
	{
		k[i].lo = x->e[i].va < x->e[i].vb ? x->e[i].va : x->e[i].vb;
		k[i].hi = x->e[i].va < x->e[i].vb ? x->e[i].vb : x->e[i].va;
		k[i].e = i;
	}
	qsort(k, x->ne, sizeof(c64_key), c64_cmp_key);
	for (i = 0; i < x->ne; i = j)
	{
		c64_edge* a = &x->e[k[i].e];
		for (j = i + 1; j < x->ne && k[j].lo == k[i].lo && k[j].hi == k[i].hi; j++)
		{
			c64_edge* b = &x->e[k[j].e];
			int same = b->va == a->va;
			a->ds += same ? b->ds : -b->ds;
			a->dc += same ? b->dc : -b->dc;
			b->alive = 0;
		}
	}
	qaws_internal_free(k);
	n = 0;
	for (i = 0; i < x->ne; i++)
		if (x->e[i].alive && (x->e[i].ds || x->e[i].dc))
			x->e[n++] = x->e[i];
	x->ne = n;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  5. Graph                                                                */
/* ======================================================================== */

static unsigned int c64_origin(c64_ctx const* x, unsigned int h)
{
	return (h & 1) ? x->e[h >> 1].vb : x->e[h >> 1].va;
}

static void c64_dir(c64_ctx const* x, unsigned int h, int64_t* dx, int64_t* dy)
{
	*dx = (h & 1) ? -x->e[h >> 1].dx : x->e[h >> 1].dx;
	*dy = (h & 1) ? -x->e[h >> 1].dy : x->e[h >> 1].dy;
}

/* CCW from the +x axis: the upper half plane (y > 0, or y = 0 and x > 0) first */
static int c64_before(c64_ctx const* x, unsigned int a, unsigned int b)
{
	int64_t ax, ay, bx, by;
	int ha, hb;
	c64_dir(x, a, &ax, &ay);
	c64_dir(x, b, &bx, &by);
	ha = (ay > 0 || (ay == 0 && ax > 0)) ? 0 : 1;
	hb = (by > 0 || (by == 0 && bx > 0)) ? 0 : 1;
	if (ha != hb)
		return ha < hb;
	return c64_cross_sign(ax, ay, bx, by) > 0;
}

static qaws_status c64_graph(c64_ctx* x)
{
	unsigned int nh = 2 * x->ne, i;
	x->h = (c64_half*)qaws_internal_malloc(qaws_internal_mul2(sizeof(c64_half), (nh + 1)));
	x->vout = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	x->vfirst = (unsigned int*)qaws_internal_zalloc(x->nv + 1, sizeof(unsigned int));
	x->vcount = (unsigned int*)qaws_internal_zalloc(x->nv + 1, sizeof(unsigned int));
	if (!x->h || !x->vout || !x->vfirst || !x->vcount)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nh; i++)
		x->vcount[c64_origin(x, i)]++;
	{
		unsigned int acc = 0;
		for (i = 0; i < x->nv; i++)
		{
			x->vfirst[i] = acc;
			acc += x->vcount[i];
			x->vcount[i] = 0;
		}
	}
	for (i = 0; i < nh; i++)
	{
		unsigned int v = c64_origin(x, i);
		x->vout[x->vfirst[v] + x->vcount[v]++] = i;
		x->h[i].face = x->h[i].next = C64_NONE;
		x->h[i].out = 0;
	}
	for (i = 0; i < x->nv; i++)
	{
		unsigned int* f = &x->vout[x->vfirst[i]];
		unsigned int n = x->vcount[i], a, b;
		for (a = 1; a < n; a++)
		{
			unsigned int key = f[a];
			for (b = a; b > 0 && c64_before(x, key, f[b - 1]); b--)
				f[b] = f[b - 1];
			f[b] = key;
		}
		for (a = 0; a < n; a++)
			x->h[f[a]].pos = a;
	}
	for (i = 0; i < nh; i++)
	{
		unsigned int tw = i ^ 1, v = c64_origin(x, tw), p = x->h[tw].pos;
		x->h[i].next = x->vout[x->vfirst[v] + (p + x->vcount[v] - 1) % x->vcount[v]];
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  6. Faces                                                                */
/* ======================================================================== */

typedef struct c64_cycle
{
	unsigned int h0, low_h, comp;
	int ccw;
	double lo[2], hi[2];
} c64_cycle;

static int c64_winding(c64_ctx const* x, unsigned int h0, c64_pt const* q);

/* the outside of a part: candidates for the cycle holding its lowest vertex */
typedef struct c64_nest_ring
{
	c64_ctx const* x;
	c64_cycle const* cyc;
	qaws_bp_box const* boxes;
	unsigned int const* of;      /* box -> cycle */
	unsigned int i;              /* the outer cycle */
	c64_pt const* q;             /* its lowest vertex */
	double p[3];
	unsigned int best;
	double best_size;
} c64_nest_ring;

static void c64_nest_visit(void* user, unsigned int item)
{
	c64_nest_ring* r = (c64_nest_ring*)user;
	unsigned int j = r->of[item];
	c64_cycle const* C = &r->cyc[j];
	double size;
	if (C->comp == r->cyc[r->i].comp)
		return;
	if (r->p[0] < C->lo[0] || r->p[0] > C->hi[0] || r->p[1] < C->lo[1] || r->p[1] > C->hi[1])
		return;
	size = (C->hi[0] - C->lo[0]) + (C->hi[1] - C->lo[1]);
	/* the smallest wins, the lower cycle on a tie (the order of a plain scan) */
	if (size > r->best_size || (size == r->best_size && j > r->best))
		return;
	if (c64_winding(r->x, C->h0, r->q) != 0)
	{
		r->best = j;
		r->best_size = size;
	}
}

static unsigned int c64_find(unsigned int* p, unsigned int i)
{
	while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
	return i;
}

/* does the cycle from h0 wind round the point q (exact crossing count)? */
static int c64_winding(c64_ctx const* x, unsigned int h0, c64_pt const* q)
{
	int w = 0;
	unsigned int h = h0;
	do
	{
		unsigned int va = c64_origin(x, h), vb = c64_origin(x, h ^ 1);
		c64_pt const* a = &x->pt[x->vrep[va]];
		c64_pt const* b = &x->pt[x->vrep[vb]];
		int a_up = c64_cmp_coord(&a->y, &a->d, &q->y, &q->d) <= 0;
		int b_up = c64_cmp_coord(&b->y, &b->d, &q->y, &q->d) <= 0;
		if (a_up != b_up)
		{
			/* side of q from the edge: d x (q - a), over q.d a.d > 0 */
			int64_t dx, dy;
			qaws_wide t1, t2, u, v, wdx, wdy, side;
			c64_dir(x, h, &dx, &dy);
			t1 = qaws_wide_mul(&q->y, &a->d); t2 = qaws_wide_mul(&a->y, &q->d);
			u = qaws_wide_sub(t1, &t2);
			t1 = qaws_wide_mul(&q->x, &a->d); t2 = qaws_wide_mul(&a->x, &q->d);
			v = qaws_wide_sub(t1, &t2);
			wdx = qaws_wide_from_i64(dx); wdy = qaws_wide_from_i64(dy);
			t1 = qaws_wide_mul(&wdx, &u); t2 = qaws_wide_mul(&wdy, &v);
			side = qaws_wide_sub(t1, &t2);
			if (a_up && qaws_wide_sign(&side) > 0) w++;
			else if (!a_up && qaws_wide_sign(&side) < 0) w--;
		}
		h = x->h[h].next;
	} while (h != h0);
	return w;
}

static qaws_status c64_faces(c64_ctx* x, unsigned int* out_faces)
{
	unsigned int nh = 2 * x->ne, i, ncyc = 0, capc = 0, nfaces = 1;
	c64_cycle* cyc = NULL;
	unsigned int* hcyc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	unsigned int* vcomp = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nv + 1)));
	unsigned int* cface;
	if (!hcyc || !vcomp)
	{
		qaws_internal_free(hcyc); qaws_internal_free(vcomp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < x->nv; i++) vcomp[i] = i;
	for (i = 0; i < x->ne; i++)
	{
		unsigned int a = c64_find(vcomp, x->e[i].va), b = c64_find(vcomp, x->e[i].vb);
		if (a != b) vcomp[a] = b;
	}
	for (i = 0; i < nh; i++) hcyc[i] = C64_NONE;
	for (i = 0; i < nh; i++)
	{
		unsigned int h = i, low = i, low_prev = C64_NONE, last = C64_NONE;
		c64_cycle c;
		if (hcyc[i] != C64_NONE) continue;
		c.lo[0] = c.lo[1] = 1e300; c.hi[0] = c.hi[1] = -1e300;
		do
		{
			unsigned int v = c64_origin(x, h);
			double const* p = &x->vf[2 * v];
			hcyc[h] = ncyc;
			if (p[0] < c.lo[0]) c.lo[0] = p[0];
			if (p[1] < c.lo[1]) c.lo[1] = p[1];
			if (p[0] > c.hi[0]) c.hi[0] = p[0];
			if (p[1] > c.hi[1]) c.hi[1] = p[1];
			if (h != i)
			{
				/* the lowest vertex (y, then x): doubles first, exactly when they tie */
				unsigned int lv = c64_origin(x, low);
				double const* q = &x->vf[2 * lv];
				double ty = 1e-9 * (fabs(p[1]) + fabs(q[1]) + 1);
				if (p[1] < q[1] - ty || (p[1] <= q[1] + ty && c64_cmp_pt(&x->pt[x->vrep[v]], &x->pt[x->vrep[lv]]) < 0))
				{
					low = h;
					low_prev = last;
				}
			}
			last = h;
			h = x->h[h].next;
		} while (h != i);
		if (low_prev == C64_NONE)
			low_prev = last;
		c.h0 = i;
		c.low_h = low;
		c.comp = c64_find(vcomp, c64_origin(x, i));
		/* at the lowest vertex the cycle turns left iff it runs counter-clockwise */
		{
			int64_t ix, iy, ox, oy;
			c64_dir(x, low_prev, &ix, &iy);
			c64_dir(x, low, &ox, &oy);
			c.ccw = c64_cross_sign(ix, iy, ox, oy) > 0;
		}
		if (ncyc == capc)
		{
			unsigned int nc = qaws_internal_grow(capc, (size_t)ncyc + 1, 64);
			c64_cycle* g = nc ? (c64_cycle*)qaws_internal_realloc(cyc, qaws_internal_mul2(sizeof(c64_cycle), nc)) : NULL;
			if (!g) { qaws_internal_free(cyc); qaws_internal_free(hcyc); qaws_internal_free(vcomp); return QAWS_STATUS_ALLOCATION_FAILURE; }
			cyc = g; capc = nc;
		}
		cyc[ncyc++] = c;
	}
	cface = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (ncyc + 1)));
	if (!cface)
	{
		qaws_internal_free(cyc); qaws_internal_free(hcyc); qaws_internal_free(vcomp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < ncyc; i++)
		cface[i] = cyc[i].ccw ? nfaces++ : C64_NONE;
	/* a part's outside: the smallest counter-clockwise cycle of another part
	   winding round its lowest vertex (smallest: nested in all the others,
	   so the one with the smallest box); the candidates from a grid over the
	   counter-clockwise cycles' boxes, in the cell holding the vertex */
	{
		unsigned int nccw = 0, *of = NULL;
		qaws_bp_box* boxes = NULL;
		qaws_bp_grid* grid = NULL;
		c64_nest_ring rg;
		for (i = 0; i < ncyc; i++) nccw += cyc[i].ccw;
		if (nccw)
		{
			boxes = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), nccw));
			of = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), nccw));
			if (!boxes || !of)
			{
				qaws_internal_free(boxes); qaws_internal_free(of); qaws_internal_free(cyc); qaws_internal_free(hcyc); qaws_internal_free(vcomp); qaws_internal_free(cface);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			nccw = 0;
			for (i = 0; i < ncyc; i++)
				if (cyc[i].ccw)
				{
					boxes[nccw].lo[0] = cyc[i].lo[0]; boxes[nccw].lo[1] = cyc[i].lo[1];
					boxes[nccw].hi[0] = cyc[i].hi[0]; boxes[nccw].hi[1] = cyc[i].hi[1];
					boxes[nccw].lo[2] = boxes[nccw].hi[2] = 0;
					of[nccw++] = i;
				}
			if (qaws_internal_grid_create(&grid, boxes, nccw, 2) != QAWS_STATUS_OK)
				grid = NULL;
		}
		rg.x = x;
		rg.cyc = cyc;
		rg.boxes = boxes;
		rg.of = of;
		for (i = 0; i < ncyc; i++)
		{
			if (cface[i] != C64_NONE) continue;
			rg.i = i;
			rg.q = &x->pt[x->vrep[c64_origin(x, cyc[i].low_h)]];
			rg.p[0] = rg.q->fx;
			rg.p[1] = rg.q->fy;
			rg.p[2] = 0;
			rg.best = C64_NONE;
			rg.best_size = 1e300;
			if (grid)
				qaws_internal_grid_ring(grid, rg.p, 0, c64_nest_visit, &rg);
			else
			{
				unsigned int j;
				for (j = 0; j < nccw; j++)
					c64_nest_visit(&rg, j);
			}
			cface[i] = rg.best == C64_NONE ? 0 : cface[rg.best];
		}
		qaws_internal_grid_destroy(grid);
		qaws_internal_free(boxes);
		qaws_internal_free(of);
	}
	for (i = 0; i < nh; i++)
		x->h[i].face = cface[hcyc[i]];
	qaws_internal_free(cyc); qaws_internal_free(hcyc); qaws_internal_free(vcomp); qaws_internal_free(cface);
	*out_faces = nfaces;
	return QAWS_STATUS_OK;
}

static int c64_inside(qaws_clip_type ct, qaws_fill_rule fr, int ws, int wc)
{
	int s = qaws_fill_rule_inside(fr, ws), c = qaws_fill_rule_inside(fr, wc);
	switch (ct)
	{
	case QAWS_CLIP_INTERSECTION: return s && c;
	case QAWS_CLIP_UNION: return s || c;
	case QAWS_CLIP_DIFFERENCE: return s && !c;
	case QAWS_CLIP_XOR: return s != c;
	default: return 0;
	}
}

static qaws_status c64_windings(c64_ctx* x, unsigned int nfaces, int* ws, int* wc, int* known)
{
	unsigned int nh = 2 * x->ne, i, *queue, *first, *list, *fill, qh = 0, qt = 0;
	queue = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	first = (unsigned int*)qaws_internal_zalloc(nfaces + 2, sizeof(unsigned int));
	list = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	fill = (unsigned int*)qaws_internal_zalloc(nfaces + 1, sizeof(unsigned int));
	if (!queue || !first || !list || !fill)
	{
		qaws_internal_free(queue); qaws_internal_free(first); qaws_internal_free(list); qaws_internal_free(fill);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nh; i++) first[x->h[i].face + 1]++;
	for (i = 0; i < nfaces; i++) first[i + 1] += first[i];
	for (i = 0; i < nh; i++)
		list[first[x->h[i].face] + fill[x->h[i].face]++] = i;
	for (i = 0; i < nfaces; i++) known[i] = 0;
	ws[0] = wc[0] = 0;
	known[0] = 1;
	queue[qt++] = 0;
	while (qh < qt)
	{
		unsigned int f = queue[qh++], k;
		for (k = first[f]; k < first[f + 1]; k++)
		{
			unsigned int h = list[k], g = x->h[h ^ 1].face;
			int sg = (h & 1) ? -1 : 1;
			if (known[g]) continue;
			ws[g] = ws[f] - sg * x->e[h >> 1].ds;
			wc[g] = wc[f] - sg * x->e[h >> 1].dc;
			known[g] = 1;
			queue[qt++] = g;
		}
	}
	qaws_internal_free(queue); qaws_internal_free(first); qaws_internal_free(list); qaws_internal_free(fill);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  8. Output                                                               */
/* ======================================================================== */

/* round(n / d), d > 0, to the nearest integer (halves away from zero) */
static int64_t c64_round(qaws_wide const* n, qaws_wide const* d)
{
	double approx = qaws_wide_to_double(n) / qaws_wide_to_double(d);
	int64_t c = (int64_t)floor(approx + 0.5);
	int k;
	/* the double is good to a few ulps: far from a half and well inside
	   2^52, its rounding is the exact one */
	if (fabs(approx) < 4503599627370496.0 / 1024 && fabs(approx - floor(approx) - 0.5) > 1e-6)
		return c;
	/* correct the double guess: |2 (n - c d)| <= d */
	for (k = 0; k < 4; k++)
	{
		qaws_wide wc = qaws_wide_from_i64(c), cd = qaws_wide_mul(&wc, d), r = qaws_wide_sub(*n, &cd), r2, md;
		r2 = qaws_wide_add(r, &r);
		md = qaws_wide_neg(*d);
		if (qaws_wide_cmp(&r2, d) > 0) c++;
		else if (qaws_wide_cmp(&r2, &md) < 0) c--;
		else break;
	}
	return c;
}

static qaws_status c64_out_point(qaws_clip64_result* r, int64_t px, int64_t py, double ex, double ey)
{
	if (r->npts == r->cappts)
	{
		unsigned int nc = qaws_internal_grow(r->cappts, (size_t)r->npts + 1, 256);
		int64_t* p = nc ? (int64_t*)qaws_internal_realloc(r->pts, qaws_internal_mul3(sizeof(int64_t), 2, nc)) : NULL;
		double* e = p ? (double*)qaws_internal_realloc(r->exact, qaws_internal_mul3(sizeof(double), 2, nc)) : NULL;
		if (p) r->pts = p;
		if (e) r->exact = e;
		if (!p || !e)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		r->cappts = nc;
	}
	r->pts[2 * r->npts] = px;
	r->pts[2 * r->npts + 1] = py;
	r->exact[2 * r->npts] = ex;
	r->exact[2 * r->npts + 1] = ey;
	r->npts++;
	return QAWS_STATUS_OK;
}

/* a path of vertices as rounded points: repeats dropped, and collinear
   corners unless preserved; 0 when nothing is left to enclose */
static int c64_emit(c64_ctx* x, unsigned int const* verts, unsigned int n, int closed, struct c64_out* op)
{
	qaws_clip64_result* r = x->r;
	int preserve = (x->d->flags & QAWS_CLIP_PRESERVE_COLLINEAR) != 0;
	unsigned int start = r->npts, i, m, changed = 1;
	op->first = start;
	op->efirst = r->nex;
	for (i = 0; i < n; i++)
	{
		c64_pt const* p = &x->pt[x->vrep[verts[i]]];
		if (i && verts[i] == verts[i - 1])
			continue;
		if (r->nex == r->capex)
		{
			unsigned int nc = qaws_internal_grow(r->capex, (size_t)r->nex + 1, 256);
			double* g = nc ? (double*)qaws_internal_realloc(r->ex, qaws_internal_mul3(sizeof(double), 2, nc)) : NULL;
			if (!g)
				return -1;
			r->ex = g;
			r->capex = nc;
		}
		r->ex[2 * r->nex] = p->fx;
		r->ex[2 * r->nex + 1] = p->fy;
		r->nex++;
	}
	op->ecount = r->nex - op->efirst;
	for (i = 0; i < n; i++)
	{
		unsigned int v = verts[i];
		int64_t px = x->vround[2 * v], py = x->vround[2 * v + 1];
		c64_pt const* p = &x->pt[x->vrep[v]];
		if (r->npts > start && r->pts[2 * r->npts - 2] == px && r->pts[2 * r->npts - 1] == py)
			continue;
		if (c64_out_point(r, px, py, p->fx, p->fy) != QAWS_STATUS_OK)
			return -1;
	}
	m = r->npts - start;
	if (closed)
		while (m > 1 && r->pts[2 * (start + m - 1)] == r->pts[2 * start] && r->pts[2 * (start + m - 1) + 1] == r->pts[2 * start + 1])
			m--;
	/* drop corners collinear with their neighbours (spikes too, always):
	   one corner per pass against its current neighbours */
	while (changed && m >= 3)
	{
		changed = 0;
		for (i = 0; i < m && m >= 3; i++)
		{
			unsigned int pi, ni, q;
			int64_t const *a, *b, *c;
			qaws_wide cr, dt;
			if (!closed && (i == 0 || i == m - 1))
				continue;
			pi = closed ? (i + m - 1) % m : i - 1;
			ni = closed ? (i + 1) % m : i + 1;
			a = &r->pts[2 * (start + pi)];
			b = &r->pts[2 * (start + i)];
			c = &r->pts[2 * (start + ni)];
			cr = qaws_wide_cross(b[0] - a[0], b[1] - a[1], c[0] - b[0], c[1] - b[1]);
			dt = qaws_wide_dot(b[0] - a[0], b[1] - a[1], c[0] - b[0], c[1] - b[1]);
			if (qaws_wide_sign(&cr) != 0 || (preserve && qaws_wide_sign(&dt) >= 0))
				continue;
			for (q = i; q + 1 < m; q++)
			{
				r->pts[2 * (start + q)] = r->pts[2 * (start + q + 1)];
				r->pts[2 * (start + q) + 1] = r->pts[2 * (start + q + 1) + 1];
				r->exact[2 * (start + q)] = r->exact[2 * (start + q + 1)];
				r->exact[2 * (start + q) + 1] = r->exact[2 * (start + q + 1) + 1];
			}
			m--;
			i = i ? i - 1 : 0;
			changed = 1;
		}
	}
	r->npts = start + m;
	op->count = m;
	if (closed && m < 3)
	{
		r->npts = start;
		r->nex = op->efirst;
		return 0;
	}
	if (!closed && m < 2)
	{
		r->npts = start;
		r->nex = op->efirst;
		return 0;
	}
	if (closed && qaws_path64_area2(&r->pts[2 * start], m) == 0)
	{
		r->npts = start;
		r->nex = op->efirst;
		return 0;
	}
	return 1;
}

static qaws_status c64_push(qaws_clip64_result* r, int open, struct c64_out const* op)
{
	if (open)
	{
		C64_GROW(r->open_paths, r->open_count, r->open_cap, struct c64_out);
		r->open_paths[r->open_count++] = *op;
	}
	else
	{
		C64_GROW(r->paths, r->path_count, r->path_cap, struct c64_out);
		r->paths[r->path_count++] = *op;
	}
	return QAWS_STATUS_OK;
}

static qaws_status c64_output(c64_ctx* x, unsigned int nfaces, int const* inside)
{
	unsigned int nh = 2 * x->ne, i, nloop = 0, cap = 0, *seq = NULL, *st = NULL, *vpos = NULL;
	unsigned int *fcomp, *comp_outer, *comp_hole, *loop_path = NULL, *loop_left = NULL, *loop_right = NULL;
	int* loop_ccw = NULL;
	unsigned char* seen;
	unsigned int first_path = x->r->path_count;
	qaws_status s = QAWS_STATUS_OK;
	fcomp = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	comp_outer = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	comp_hole = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	seen = (unsigned char*)qaws_internal_scratch_zalloc(nh + 1, 1);
	vpos = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nv + 1)));
	if (!fcomp || !comp_outer || !comp_hole || !seen || !vpos)
	{
		s = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < nh; i++)
		x->h[i].out = inside[x->h[i].face] && !inside[x->h[i ^ 1].face];
	for (i = 0; i < nfaces; i++) { fcomp[i] = i; comp_outer[i] = comp_hole[i] = C64_NONE; }
	for (i = 0; i < nh; i += 2)
		if (inside[x->h[i].face] == inside[x->h[i + 1].face])
		{
			unsigned int a = c64_find(fcomp, x->h[i].face), b = c64_find(fcomp, x->h[i + 1].face);
			if (a != b) fcomp[a] = b;
		}
	for (i = 0; i < x->nv; i++) vpos[i] = C64_NONE;
	for (i = 0; i < nh && s == QAWS_STATUS_OK; i++)
	{
		unsigned int h = i, n = 0, top = 0, k;
		if (!x->h[i].out || seen[i])
			continue;
		do
		{
			unsigned int tw, v, p, g = C64_NONE;
			seen[h] = 1;
			if ((size_t)n + 1 > cap)
			{
				unsigned int nc = qaws_internal_grow(cap, (size_t)n + 1, 64);
				unsigned int* a1 = nc ? (unsigned int*)qaws_internal_realloc(seq, qaws_internal_mul2(sizeof(unsigned int), nc)) : NULL;
				unsigned int* a2 = a1 ? (unsigned int*)qaws_internal_realloc(st, qaws_internal_mul2(sizeof(unsigned int), nc)) : NULL;
				if (a1) seq = a1;
				if (a2) st = a2;
				if (!a2) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
				cap = nc;
			}
			seq[n++] = h;
			tw = h ^ 1;
			v = c64_origin(x, tw);
			p = x->h[tw].pos;
			for (k = 1; k <= x->vcount[v]; k++)
			{
				unsigned int cand = x->vout[x->vfirst[v] + (p + x->vcount[v] - k) % x->vcount[v]];
				if (x->h[cand].out) { g = cand; break; }
			}
			if (g == C64_NONE) { s = QAWS_STATUS_INTERNAL_ERROR; goto done; }
			h = g;
		} while (h != i && !seen[h]);
		/* cut at repeated vertices into simple loops */
		for (k = 0; k < n && s == QAWS_STATUS_OK; k++)
		{
			unsigned int va = c64_origin(x, seq[k]), vb = c64_origin(x, seq[k] ^ 1), from;
			if (vpos[va] == C64_NONE)
				vpos[va] = top;
			st[top++] = seq[k];
			from = vpos[vb];
			if (from == C64_NONE)
				continue;
			{
				unsigned int len = top - from, a, low = from, *verts;
				int ccw;
				struct c64_out op;
				/* orientation: the turn at the lowest vertex */
				for (a = from; a < top; a++)
				{
					vpos[c64_origin(x, st[a])] = C64_NONE;
					if (c64_cmp_pt(&x->pt[x->vrep[c64_origin(x, st[a])]], &x->pt[x->vrep[c64_origin(x, st[low])]]) < 0)
						low = a;
				}
				{
					unsigned int prev = low == from ? top - 1 : low - 1;
					int64_t ix, iy, ox, oy;
					c64_dir(x, st[prev], &ix, &iy);
					c64_dir(x, st[low], &ox, &oy);
					ccw = c64_cross_sign(ix, iy, ox, oy) > 0;
				}
				verts = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (len + 1)));
				if (!verts) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
				for (a = 0; a < len; a++)
				{
					unsigned int idx = (x->d->flags & QAWS_CLIP_REVERSE_SOLUTION) ? from + len - 1 - a : from + a;
					verts[a] = (x->d->flags & QAWS_CLIP_REVERSE_SOLUTION) ? c64_origin(x, st[idx] ^ 1) : c64_origin(x, st[idx]);
				}
				memset(&op, 0, sizeof(op));
				if (c64_emit(x, verts, len, 1, &op) > 0)
				{
					unsigned int* g1 = (unsigned int*)qaws_internal_realloc(loop_path, qaws_internal_mul2(sizeof(unsigned int), (nloop + 1)));
					unsigned int* g2 = g1 ? (unsigned int*)qaws_internal_realloc(loop_left, qaws_internal_mul2(sizeof(unsigned int), (nloop + 1))) : NULL;
					unsigned int* g3 = g2 ? (unsigned int*)qaws_internal_realloc(loop_right, qaws_internal_mul2(sizeof(unsigned int), (nloop + 1))) : NULL;
					int* g4 = g3 ? (int*)qaws_internal_realloc(loop_ccw, qaws_internal_mul2(sizeof(int), (nloop + 1))) : NULL;
					if (g1) loop_path = g1;
					if (g2) loop_left = g2;
					if (g3) loop_right = g3;
					if (g4) loop_ccw = g4;
					if (!g4) { qaws_internal_free(verts); s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
					loop_path[nloop] = x->r->path_count;
					loop_left[nloop] = c64_find(fcomp, x->h[st[from]].face);
					loop_right[nloop] = c64_find(fcomp, x->h[st[from] ^ 1].face);
					loop_ccw[nloop] = ccw;
					if (ccw) comp_outer[loop_left[nloop]] = x->r->path_count;
					else comp_hole[loop_right[nloop]] = x->r->path_count;
					s = c64_push(x->r, 0, &op);
					nloop++;
				}
				qaws_internal_free(verts);
				top = from;
			}
		}
	}
	/* nesting */
	for (i = 0; s == QAWS_STATUS_OK && i < nloop; i++)
	{
		struct c64_out* op = &x->r->paths[loop_path[i]];
		unsigned int par = loop_ccw[i] ? comp_hole[loop_right[i]] : comp_outer[loop_left[i]];
		op->parent = par;
	}
	for (i = first_path; s == QAWS_STATUS_OK && i < x->r->path_count; i++)
	{
		struct c64_out* op = &x->r->paths[i];
		unsigned int d = 0, p = op->parent;
		while (p != C64_NONE && d <= x->r->path_count)
		{
			d++;
			p = x->r->paths[p].parent;
		}
		op->depth = d;
		op->hole = (d & 1) != 0;
	}
done:
	qaws_internal_free(fcomp); qaws_internal_free(comp_outer); qaws_internal_free(comp_hole); qaws_internal_free(seen); qaws_internal_free(vpos);
	qaws_internal_free(seq); qaws_internal_free(st); qaws_internal_free(loop_path); qaws_internal_free(loop_left); qaws_internal_free(loop_right); qaws_internal_free(loop_ccw);
	return s;
}

/* ======================================================================== */
/*  Open subjects                                                           */
/* ======================================================================== */

/* winding of one operand's closed segments round the point q; *on when q
   lies on one of them */
static int c64_operand_winding(c64_ctx const* x, unsigned int operand, c64_pt const* q, int* on)
{
	int w = 0;
	unsigned int i;
	for (i = 0; i < x->nseg; i++)
	{
		c64_seg const* s = &x->seg[i];
		qaws_wide ay = qaws_wide_from_i64(s->ay), by = qaws_wide_from_i64(s->by), one = qaws_wide_from_i64(1);
		qaws_wide ax = qaws_wide_from_i64(s->ax);
		qaws_wide t1, t2, u, v, wdx, wdy, side;
		int a_up, b_up;
		if (s->operand != operand)
			continue;
		a_up = c64_cmp_coord(&ay, &one, &q->y, &q->d) <= 0;
		b_up = c64_cmp_coord(&by, &one, &q->y, &q->d) <= 0;
		/* side of q: d x (q - a) times q.d */
		t1 = q->y; t2 = qaws_wide_mul(&ay, &q->d); u = qaws_wide_sub(t1, &t2);
		t1 = q->x; t2 = qaws_wide_mul(&ax, &q->d); v = qaws_wide_sub(t1, &t2);
		wdx = qaws_wide_from_i64(s->bx - s->ax); wdy = qaws_wide_from_i64(s->by - s->ay);
		t1 = qaws_wide_mul(&wdx, &u); t2 = qaws_wide_mul(&wdy, &v);
		side = qaws_wide_sub(t1, &t2);
		if (qaws_wide_sign(&side) == 0)
		{
			/* on the segment's line: on the segment when inside its box */
			qaws_wide lx = qaws_wide_from_i64(s->lox), hx = qaws_wide_from_i64(s->hix);
			qaws_wide ly = qaws_wide_from_i64(s->loy), hy = qaws_wide_from_i64(s->hiy);
			if (c64_cmp_coord(&q->x, &q->d, &lx, &one) >= 0 && c64_cmp_coord(&q->x, &q->d, &hx, &one) <= 0 &&
				c64_cmp_coord(&q->y, &q->d, &ly, &one) >= 0 && c64_cmp_coord(&q->y, &q->d, &hy, &one) <= 0)
				*on = 1;
		}
		if (a_up != b_up)
		{
			if (a_up && qaws_wide_sign(&side) > 0) w++;
			else if (!a_up && qaws_wide_sign(&side) < 0) w--;
		}
	}
	return w;
}

static qaws_status c64_open(c64_ctx* x)
{
	qaws_clip64_desc const* d = x->d;
	unsigned int i, n = 0, cap = 0, *verts = NULL, last_path = C64_NONE, last_v = C64_NONE;
	qaws_status s = QAWS_STATUS_OK;
	for (i = 0; i <= x->noe && s == QAWS_STATUS_OK; i++)
	{
		c64_edge const* e = i < x->noe ? &x->oe[i] : NULL;
		int keep = 0;
		if (e)
		{
			/* the middle (a + b) / 2 */
			c64_pt const* a = &x->pt[x->vrep[e->va]];
			c64_pt const* b = &x->pt[x->vrep[e->vb]];
			c64_pt m;
			qaws_wide t1, t2, two = qaws_wide_from_i64(2);
			int on = 0, ws, wc, in_s, in_c;
			t1 = qaws_wide_mul(&a->x, &b->d); t2 = qaws_wide_mul(&b->x, &a->d); m.x = qaws_wide_add(t1, &t2);
			t1 = qaws_wide_mul(&a->y, &b->d); t2 = qaws_wide_mul(&b->y, &a->d); m.y = qaws_wide_add(t1, &t2);
			t1 = qaws_wide_mul(&a->d, &b->d); m.d = qaws_wide_mul(&t1, &two);
			ws = c64_operand_winding(x, 0, &m, &on);
			wc = c64_operand_winding(x, 1, &m, &on);
			in_s = !on && qaws_fill_rule_inside(d->fill_rule, ws);
			in_c = !on && qaws_fill_rule_inside(d->fill_rule, wc);
			switch (d->clip_type)
			{
			case QAWS_CLIP_INTERSECTION: keep = in_c; break;
			case QAWS_CLIP_UNION: keep = !in_s && !in_c; break;
			case QAWS_CLIP_NONE: keep = 0; break;
			default: keep = !in_c; break;
			}
		}
		if (n && (!keep || e->open_path != last_path || e->va != last_v))
		{
			struct c64_out op;
			memset(&op, 0, sizeof(op));
			op.parent = C64_NONE;
			if (c64_emit(x, verts, n, 0, &op) > 0)
				s = c64_push(x->r, 1, &op);
			n = 0;
		}
		if (keep)
		{
			if ((size_t)n + 2 > cap)
			{
				unsigned int nc = qaws_internal_grow(cap, (size_t)n + 2, 64);
				unsigned int* g = nc ? (unsigned int*)qaws_internal_realloc(verts, qaws_internal_mul2(sizeof(unsigned int), nc)) : NULL;
				if (!g) { qaws_internal_free(verts); return QAWS_STATUS_ALLOCATION_FAILURE; }
				verts = g; cap = nc;
			}
			if (n == 0)
				verts[n++] = e->va;
			verts[n++] = e->vb;
			last_path = e->open_path;
			last_v = e->vb;
		}
	}
	qaws_internal_free(verts);
	return s;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

double qaws_path64_area2(int64_t const* p, unsigned int n)
{
	qaws_wide sum = qaws_wide_from_i64(0);
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + 1) % n;
		qaws_wide c = qaws_wide_cross(p[2 * i], p[2 * i + 1], p[2 * j], p[2 * j + 1]);
		sum = qaws_wide_add(sum, &c);
	}
	return qaws_wide_to_double(&sum);
}

static void c64_free(c64_ctx* x)
{
	qaws_internal_free(x->seg); qaws_internal_free(x->pt); qaws_internal_free(x->cut); qaws_internal_free(x->vid); qaws_internal_free(x->vrep); qaws_internal_free(x->vf); qaws_internal_free(x->e); qaws_internal_free(x->oe);
	qaws_internal_free(x->h); qaws_internal_free(x->vfirst); qaws_internal_free(x->vcount); qaws_internal_free(x->vout); qaws_internal_free(x->vround);
}

qaws_status qaws_clip64_execute(qaws_clip64_result** out_result, qaws_clip64_desc const* desc)
{
	c64_ctx x;
	qaws_status s;
	unsigned int nfaces = 0, i;
	int *ws = NULL, *wc = NULL, *known = NULL, *inside = NULL;
	if (!desc || !out_result || (desc->subject_count && !desc->subjects) || (desc->clip_count && !desc->clips) ||
		(desc->open_subject_count && !desc->open_subjects))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_result = NULL;
	memset(&x, 0, sizeof(x));
	x.d = desc;
	x.r = (qaws_clip64_result*)qaws_internal_zalloc(1, sizeof(qaws_clip64_result));
	if (!x.r)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s = c64_add_paths(&x, desc->subjects, desc->subject_count, 0);
	if (s == QAWS_STATUS_OK) s = c64_add_paths(&x, desc->clips, desc->clip_count, 1);
	if (s == QAWS_STATUS_OK) s = c64_add_paths(&x, desc->open_subjects, desc->open_subject_count, 2);
	if (s == QAWS_STATUS_OK && desc->clip_type != QAWS_CLIP_NONE && x.nseg)
	{
		s = c64_hits(&x);
		if (s == QAWS_STATUS_OK) s = c64_vertices(&x);
		if (s == QAWS_STATUS_OK)
		{
			x.vround = (int64_t*)qaws_internal_malloc(qaws_internal_mul3(sizeof(int64_t), 2, (x.nv + 1)));
			if (!x.vround)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; s == QAWS_STATUS_OK && i < x.nv; i++)
			{
				c64_pt const* p = &x.pt[x.vrep[i]];
				x.vround[2 * i] = c64_round(&p->x, &p->d);
				x.vround[2 * i + 1] = c64_round(&p->y, &p->d);
			}
		}
		if (s == QAWS_STATUS_OK) s = c64_edges(&x);
		if (s == QAWS_STATUS_OK) s = c64_graph(&x);
		if (s == QAWS_STATUS_OK) s = c64_faces(&x, &nfaces);
		if (s == QAWS_STATUS_OK)
		{
			ws = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (nfaces + 1)));
			wc = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (nfaces + 1)));
			known = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (nfaces + 1)));
			inside = (int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(int), (nfaces + 1)));
			if (!ws || !wc || !known || !inside)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
		}
		if (s == QAWS_STATUS_OK) s = c64_windings(&x, nfaces, ws, wc, known);
		if (s == QAWS_STATUS_OK)
		{
			for (i = 0; i < nfaces; i++)
				inside[i] = known[i] && c64_inside(desc->clip_type, desc->fill_rule, ws[i], wc[i]);
			s = c64_output(&x, nfaces, inside);
		}
		if (s == QAWS_STATUS_OK && x.noe)
			s = c64_open(&x);
	}
	qaws_internal_free(ws); qaws_internal_free(wc); qaws_internal_free(known); qaws_internal_free(inside);
	c64_free(&x);
	if (s != QAWS_STATUS_OK)
	{
		qaws_clip64_result_destroy(x.r);
		return s;
	}
	*out_result = x.r;
	return QAWS_STATUS_OK;
}

void qaws_clip64_result_destroy(qaws_clip64_result* r)
{
	if (!r) return;
	qaws_internal_free(r->paths); qaws_internal_free(r->open_paths); qaws_internal_free(r->pts); qaws_internal_free(r->exact); qaws_internal_free(r->ex);
	qaws_internal_free(r);
}

unsigned int qaws_clip64_result_get_path_count(qaws_clip64_result const* r) { return r ? r->path_count : 0; }

qaws_status qaws_clip64_result_get_path(int64_t const** out, unsigned int* n, qaws_clip64_result const* r, unsigned int i)
{
	if (!r || !out || !n || i >= r->path_count) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = r->pts + 2 * r->paths[i].first;
	*n = r->paths[i].count;
	return QAWS_STATUS_OK;
}

qaws_status qaws_clip64_result_get_path_exact(double const** out, unsigned int* n, qaws_clip64_result const* r, unsigned int i)
{
	if (!r || !out || !n || i >= r->path_count) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = r->ex + 2 * r->paths[i].efirst;
	*n = r->paths[i].ecount;
	return QAWS_STATUS_OK;
}

unsigned int qaws_clip64_result_get_parent(qaws_clip64_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].parent : C64_NONE;
}

int qaws_clip64_result_is_hole(qaws_clip64_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].hole : 0;
}

unsigned int qaws_clip64_result_get_depth(qaws_clip64_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].depth : 0;
}

unsigned int qaws_clip64_result_get_open_path_count(qaws_clip64_result const* r) { return r ? r->open_count : 0; }

qaws_status qaws_clip64_result_get_open_path(int64_t const** out, unsigned int* n, qaws_clip64_result const* r, unsigned int i)
{
	if (!r || !out || !n || i >= r->open_count) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = r->pts + 2 * r->open_paths[i].first;
	*n = r->open_paths[i].count;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Internal: results put together                                          */
/* ======================================================================== */

qaws_status qaws_internal_clip64_result_empty(qaws_clip64_result** out)
{
	*out = (qaws_clip64_result*)qaws_internal_zalloc(1, sizeof(qaws_clip64_result));
	return *out ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}

qaws_status qaws_internal_clip64_result_append(qaws_clip64_result* dst, qaws_clip64_result* src)
{
	unsigned int ptoff = dst->npts, exoff = dst->nex, poff = dst->path_count, i;
	for (i = 0; i < src->npts; i++)
		if (c64_out_point(dst, src->pts[2 * i], src->pts[2 * i + 1], src->exact[2 * i], src->exact[2 * i + 1]) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < src->nex; i++)
	{
		if (dst->nex == dst->capex)
		{
			unsigned int nc = qaws_internal_grow(dst->capex, (size_t)dst->nex + 1, 256);
			double* g = nc ? (double*)qaws_internal_realloc(dst->ex, qaws_internal_mul3(sizeof(double), 2, nc)) : NULL;
			if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
			dst->ex = g;
			dst->capex = nc;
		}
		dst->ex[2 * dst->nex] = src->ex[2 * i];
		dst->ex[2 * dst->nex + 1] = src->ex[2 * i + 1];
		dst->nex++;
	}
	for (i = 0; i < src->path_count; i++)
	{
		struct c64_out p = src->paths[i];
		p.first += ptoff;
		p.efirst += exoff;
		if (p.parent != C64_NONE) p.parent += poff;
		if (c64_push(dst, 0, &p) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < src->open_count; i++)
	{
		struct c64_out p = src->open_paths[i];
		p.first += ptoff;
		p.efirst += exoff;
		if (c64_push(dst, 1, &p) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	qaws_clip64_result_destroy(src);
	return QAWS_STATUS_OK;
}
