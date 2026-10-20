/*
 * Exact Boolean operations on regions of curves: the planar arrangement of
 * qaws_clip with every topological decision taken on the exact kernel.
 *
 *   1. sources: the curves of every path (polylines and composites piece by
 *      piece), prepared onto one exact lattice; a path whose pieces do not
 *      meet exactly gets the lines between them
 *   2. hits: the certified batch (self-intersections included)
 *   3. events: the hits' parameters on each source, the sources' ends and
 *      the overlaps' ends; exact points where the parameter is rational
 *   4. vertices: events joined when they are proven to be one point (equal
 *      exact points, an exact point on the partner curve, a common root of
 *      the two defining polynomials)
 *   5. edges: sources cut at their vertices; edges along an overlap merged
 *      with their winding steps summed, cancelled ones dropped
 *   6. graph: half-edges sorted round each vertex from the exact sub-curves
 *      leaving it (directions, then curvature; cones at irrational vertices)
 *   7. faces: the cycles; each part's outer cycle lies in the face of the
 *      others hit first by a ray from an exact point of it
 *   8. windings across edges from the unbounded face, the truth tables, and
 *      qaws_clip's output stage
 */

#include "qaws_clip_exact.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_bezier.h"
#include "qaws_bspline.h"
#include "qaws_composite.h"
#include "qaws_path.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_clip.h"
#include "internal/qaws_internal_path.h"
#include "internal/qaws_internal_broadphase.h"
#include "exact/qaws_exact_curve.h"
#include "exact/qaws_exact_poly.h"
#include "exact/qaws_exact_roots.h"
#include "exact/qaws_exact_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CX_NONE QAWS_CLIP_NONE_INDEX
#define TRY(x) do { st = (x); if (st != QAWS_STATUS_OK) return st; } while (0)

typedef qaws_exact_int xint;

/* ======================================================================== */
/*  Data                                                                    */
/* ======================================================================== */

typedef struct cx_src
{
	qaws_curve const* c;
	qaws_exact_curve* e;
	unsigned int operand, path, index;
	double base, scale;
	double t0, t1;
	unsigned int first_cut, cut_count;  /* its cuts, sorted along it */
	unsigned int first_edge;            /* its edges: first_edge + k between cut k and k + 1 */
	double lo[2], hi[2];                /* box, world units */
	int is_line;                        /* one straight span: a, b its ends (world units) */
	double la[2], lb[2];
} cx_src;

/* A homogeneous point (x, y, w), w > 0, lattice units. */
typedef struct cx_pt
{
	xint x, y, w;
} cx_pt;

/* A rational parameter: span k of the source, local s = num / den (0 <= s <= 1, den > 0). */
typedef struct cx_rat
{
	unsigned int k;
	xint num, den;
} cx_rat;

typedef struct cx_event
{
	unsigned int src;
	double lo, hi;                      /* parameter enclosure (lo == hi: exactly that double) */
	unsigned int pt;                    /* exact point (CX_NONE: not known) */
	unsigned int rat;                   /* exact rational parameter (CX_NONE: not known) */
	unsigned int partner;               /* the hit's other event */
	unsigned int parent;                /* union-find: one vertex */
	unsigned int vertex;
} cx_event;

typedef struct cx_cut
{
	unsigned int ev;                    /* representative event (the one with a rational parameter if any) */
	unsigned int vertex;
	double t;                           /* the double for output */
} cx_cut;

typedef struct cx_edge
{
	unsigned int src;
	unsigned int c0, c1;                /* cuts (global index) */
	unsigned int v0, v1;
	int ds, dc;
	unsigned int parent;                /* union-find over overlaps */
	int flip;                           /* runs against the parent */
	unsigned int live;                  /* index among the live edges (CX_NONE: merged or dropped) */
	unsigned int comp;
} cx_edge;

typedef struct cx_ctx
{
	qaws_clip_exact_desc const* d;
	qaws_exact_desc ed;
	cx_src* src;
	unsigned int nsrc, capsrc;
	qaws_curve** temp;
	unsigned int ntemp, captemp;
	cx_pt* pts;
	unsigned int npts, cappts;
	double* ptf;                        /* the points as doubles (x, y), lattice units */
	unsigned int capptf;
	cx_rat* rats;
	unsigned int nrats, caprats;
	cx_event* ev;
	unsigned int nev, capev;
	unsigned int* sev_first;            /* source s's events: sev[sev_first[s] .. sev_first[s + 1]) */
	unsigned int* sev;
	qaws_exact_batch_hit* hits;
	unsigned int nhits;
	cx_cut* cuts;
	unsigned int ncuts;
	cx_edge* e;
	unsigned int ne;
	unsigned int nv;
	double ext, scale;                  /* world extent; lattice unit in world units */
	qaws_clip_exact_stats stats;
} cx_ctx;

#define CX_GROW(ptr, n, cap, type) \
	do { if ((n) >= 0x7FFFFFFFu) { return QAWS_STATUS_ALLOCATION_FAILURE; } /* doubling the capacity, and 2 n half-edges, stay in unsigned int */ \
		if ((n) == (cap)) { unsigned int nc_ = (cap) ? 2 * (cap) : 64; type* np_ = (type*)qaws_internal_realloc((ptr), qaws_internal_mul2(sizeof(type), nc_)); \
		if (!np_) { return QAWS_STATUS_ALLOCATION_FAILURE; } (ptr) = np_; (cap) = nc_; } } while (0)

static unsigned int cx_find(cx_event* ev, unsigned int i)
{
	while (ev[i].parent != i)
	{
		ev[i].parent = ev[ev[i].parent].parent;
		i = ev[i].parent;
	}
	return i;
}

static void cx_union(cx_event* ev, unsigned int a, unsigned int b)
{
	a = cx_find(ev, a);
	b = cx_find(ev, b);
	if (a != b)
		ev[a < b ? b : a].parent = a < b ? a : b;
}

/* ======================================================================== */
/*  Exact helpers                                                           */
/* ======================================================================== */

static qaws_status cx_pt_norm(cx_pt* p)
{
	if (qaws_exact_int_sign(&p->w) < 0)
	{
		qaws_exact_int_neg(&p->x, &p->x);
		qaws_exact_int_neg(&p->y, &p->y);
		qaws_exact_int_neg(&p->w, &p->w);
	}
	return qaws_exact_int_is_zero(&p->w) ? QAWS_STATUS_CERTIFICATION_FAILED : QAWS_STATUS_OK;
}

/* sign of a.c / a.w - b.c / b.w for coordinate c (0 x, 1 y) */
static qaws_status cx_pt_cmp1(cx_pt const* a, cx_pt const* b, int c, int* out)
{
	xint l, r;
	qaws_status st;
	TRY(qaws_exact_int_mul(&l, c ? &a->y : &a->x, &b->w));
	TRY(qaws_exact_int_mul(&r, c ? &b->y : &b->x, &a->w));
	*out = qaws_exact_int_cmp(&l, &r);
	return QAWS_STATUS_OK;
}

static qaws_status cx_pt_cmp(cx_pt const* a, cx_pt const* b, int* out)
{
	qaws_status st;
	TRY(cx_pt_cmp1(a, b, 0, out));
	if (*out == 0)
		TRY(cx_pt_cmp1(a, b, 1, out));
	return QAWS_STATUS_OK;
}

static void cx_pt_world(cx_ctx const* x, cx_pt const* p, double* out)
{
	out[0] = ldexp(qaws_exact_ratio_to_double(&p->x, &p->w), x->ed.space_exp2);
	out[1] = ldexp(qaws_exact_ratio_to_double(&p->y, &p->w), x->ed.space_exp2);
}

/* the point exactly as doubles (world units), when it is one */
static int cx_pt_double(cx_ctx const* x, cx_pt const* p, double* out)
{
	double lo, hi;
	int ex0, ex1;
	qaws_exact_ratio_enclose(&p->x, &p->w, &lo, &hi, &ex0);
	out[0] = ldexp(lo, x->ed.space_exp2);
	qaws_exact_ratio_enclose(&p->y, &p->w, &lo, &hi, &ex1);
	out[1] = ldexp(lo, x->ed.space_exp2);
	return ex0 && ex1 && fabs(out[0]) < 1e300 && fabs(out[1]) < 1e300;
}

/*
 * De Casteljau at s = num / den in integers: L over [0, s], R over [s, 1],
 * every control point at the common scale den^n (so the rational curves are
 * the same). Either output may be NULL.
 */
static qaws_status cx_split(xint const* h, unsigned int n, xint const* num, xint const* den, xint* L, xint* R)
{
	/* low degrees on the stack */
	xint wbuf[3 * 8], pwbuf[8];
	xint* w = n < 8 ? wbuf : (xint*)qaws_internal_malloc(qaws_internal_mul3(sizeof(xint), 3, (n + 1)));
	xint* pw = n < 8 ? pwbuf : (xint*)qaws_internal_malloc(qaws_internal_mul2(sizeof(xint), (n + 1)));
	xint om, t1, t2;
	unsigned int i, k, c;
	qaws_status st = QAWS_STATUS_OK;
	if (!w || !pw)
	{
		if (w != wbuf) qaws_internal_free(w);
		if (pw != pwbuf) qaws_internal_free(pw);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(w, h, sizeof(xint) * 3 * (n + 1));
	st = qaws_exact_int_sub(&om, den, num);
	qaws_exact_int_from_i64(&pw[0], 1);
	for (k = 1; k <= n && st == QAWS_STATUS_OK; k++)
		st = qaws_exact_int_mul(&pw[k], &pw[k - 1], den);
	for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
	{
		if (L) st = qaws_exact_int_mul(&L[c], &w[c], &pw[n]);
		if (R && st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&R[3 * n + c], &w[3 * n + c], &pw[n]);
	}
	for (k = 1; k <= n && st == QAWS_STATUS_OK; k++)
	{
		for (i = 0; i + k <= n && st == QAWS_STATUS_OK; i++)
			for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
			{
				st = qaws_exact_int_mul(&t1, &w[3 * i + c], &om);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t2, &w[3 * (i + 1) + c], num);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_add(&w[3 * i + c], &t1, &t2);
			}
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
		{
			if (L) st = qaws_exact_int_mul(&L[3 * k + c], &w[c], &pw[n - k]);
			if (R && st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&R[3 * (n - k) + c], &w[3 * (n - k) + c], &pw[n - k]);
		}
	}
	if (w != wbuf) qaws_internal_free(w);
	if (pw != pwbuf) qaws_internal_free(pw);
	return st;
}

/* divide every value by their common gcd (keeps the rational curve) */
static void cx_reduce(xint* v, unsigned int count)
{
	xint g;
	unsigned int i;
	qaws_exact_int_zero(&g);
	for (i = 0; i < count; i++)
		qaws_exact_int_gcd(&g, &g, &v[i]);
	if (qaws_exact_int_sign(&g) < 0)
		qaws_exact_int_neg(&g, &g);
	if (qaws_exact_int_is_zero(&g) || qaws_exact_int_bits(&g) <= 1)
		return;
	for (i = 0; i < count; i++)
		qaws_exact_int_divmod(&v[i], NULL, &v[i], &g);
}

/* the span's point at local s = num / den */
static qaws_status cx_span_point(qaws_exact_span const* sp, xint const* num, xint const* den, cx_pt* out)
{
	unsigned int n = sp->degree;
	xint rbuf[3 * 8];
	xint* R = n < 8 ? rbuf : (xint*)qaws_internal_malloc(qaws_internal_mul3(sizeof(xint), 3, (n + 1)));
	qaws_status st;
	if (!R)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = cx_split(sp->h, n, num, den, NULL, R);
	if (st == QAWS_STATUS_OK)
	{
		out->x = R[0];
		out->y = R[1];
		out->w = R[2];
		cx_reduce(&out->x, 3);
		st = cx_pt_norm(out);
	}
	if (R != rbuf) qaws_internal_free(R);
	return st;
}

/*
 * A double parameter t of curve e as span k and local s = num / den: the
 * span that starts at t when side > 0 (and t is a knot), the one that ends
 * there when side < 0.
 */
static qaws_status cx_local(qaws_exact_curve const* e, double t, int side, unsigned int* out_k, xint* num, xint* den)
{
	unsigned int k, n = e->span_count;
	int64_t m;
	int ex, q;
	qaws_exact_span const* sp;
	qaws_status st;
	xint T, a, L;
	for (k = 0; k + 1 < n; k++)
	{
		double b = ldexp((double)e->spans[k].b, -e->param_shift);
		if (t < b || (t == b && side < 0))
			break;
	}
	sp = &e->spans[k];
	*out_k = k;
	if (!qaws_exact_split_double(t, &m, &ex))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* T = t 2^shift = m 2^q */
	q = ex + e->param_shift;
	qaws_exact_int_from_i64(&T, m);
	qaws_exact_int_from_i64(&a, sp->a);
	qaws_exact_int_from_i64(&L, sp->b - sp->a);
	if (q >= 0)
	{
		TRY(qaws_exact_int_shl(&T, &T, (unsigned int)q));
		TRY(qaws_exact_int_sub(num, &T, &a));
		*den = L;
	}
	else
	{
		TRY(qaws_exact_int_shl(&a, &a, (unsigned int)-q));
		TRY(qaws_exact_int_sub(num, &T, &a));
		TRY(qaws_exact_int_shl(den, &L, (unsigned int)-q));
	}
	/* clamp to the span (an enclosure end a few ulps past a knot) */
	if (qaws_exact_int_sign(num) < 0)
		qaws_exact_int_zero(num);
	if (qaws_exact_int_cmp(num, den) > 0)
		*num = *den;
	return QAWS_STATUS_OK;
}

/* sign of a / b - c / d (b, d > 0) */
static qaws_status cx_ratio_cmp(xint const* a, xint const* b, xint const* c, xint const* d, int* out)
{
	xint l, r;
	qaws_status st;
	TRY(qaws_exact_int_mul(&l, a, d));
	TRY(qaws_exact_int_mul(&r, c, b));
	*out = qaws_exact_int_cmp(&l, &r);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  1. Sources                                                              */
/* ======================================================================== */

static double cx_tmin(qaws_curve const* c) { return c->parameter_range.min_value; }
static double cx_tmax(qaws_curve const* c) { return c->parameter_range.max_value; }

static qaws_status cx_line(cx_ctx* x, double const* a, double const* b, qaws_curve** out)
{
	qaws_scalar cp[4];
	qaws_bezier_desc d;
	qaws_status s;
	cp[0] = (qaws_scalar)a[0]; cp[1] = (qaws_scalar)a[1]; cp[2] = (qaws_scalar)b[0]; cp[3] = (qaws_scalar)b[1];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cp;
	d.control_point_count = 2;
	CX_GROW(x->temp, x->ntemp, x->captemp, qaws_curve*);
	s = qaws_curve_create_bezier(out, &d);
	if (s == QAWS_STATUS_OK)
		x->temp[x->ntemp++] = *out;
	return s;
}

static qaws_status cx_add_curve(cx_ctx* x, qaws_curve const* c, unsigned int operand, unsigned int path, unsigned int index, double base, double scale)
{
	cx_src* s;
	qaws_status st;
	if (c->kind == QAWS_CURVE_KIND_BSPLINE && c->degree == 1)
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)c->impl;
		unsigned int k;
		for (k = 0; k + 1 < impl->control_point_count; k++)
		{
			double a[2], b[2], t0 = impl->knots[k + 1], t1 = impl->knots[k + 2];
			qaws_curve* line = NULL;
			a[0] = impl->control_points[2 * k]; a[1] = impl->control_points[2 * k + 1];
			b[0] = impl->control_points[2 * k + 2]; b[1] = impl->control_points[2 * k + 3];
			if ((a[0] == b[0] && a[1] == b[1]) || !(t1 > t0))
				continue;
			TRY(cx_line(x, a, b, &line));
			TRY(cx_add_curve(x, line, operand, path, index, base + scale * t0, scale * (t1 - t0)));
		}
		return QAWS_STATUS_OK;
	}
	if (c->kind == QAWS_CURVE_KIND_COMPOSITE)
	{
		qaws_composite_impl const* impl = (qaws_composite_impl const*)c->impl;
		unsigned int k;
		for (k = 0; k < impl->segment_count; k++)
		{
			qaws_curve const* g = impl->segments[k];
			double a = cx_tmin(g), b = cx_tmax(g);
			TRY(cx_add_curve(x, g, operand, path, index, base + scale * (k - a / (b - a)), scale / (b - a)));
		}
		return QAWS_STATUS_OK;
	}
	CX_GROW(x->src, x->nsrc, x->capsrc, cx_src);
	s = &x->src[x->nsrc];
	memset(s, 0, sizeof(*s));
	s->c = c;
	s->operand = operand;
	s->path = path;
	s->index = index;
	s->base = base;
	s->scale = scale;
	s->t0 = cx_tmin(c);
	s->t1 = cx_tmax(c);
	TRY(qaws_exact_curve_prepare(&s->e, NULL, c, &x->ed));
	if (s->e->span_count == 1 && s->e->spans[0].degree == 1)
	{
		xint const* h = s->e->spans[0].h;
		s->is_line = 1;
		s->la[0] = ldexp(qaws_exact_ratio_to_double(&h[0], &h[2]), x->ed.space_exp2);
		s->la[1] = ldexp(qaws_exact_ratio_to_double(&h[1], &h[2]), x->ed.space_exp2);
		s->lb[0] = ldexp(qaws_exact_ratio_to_double(&h[3], &h[5]), x->ed.space_exp2);
		s->lb[1] = ldexp(qaws_exact_ratio_to_double(&h[4], &h[5]), x->ed.space_exp2);
	}
	x->nsrc++;
	return QAWS_STATUS_OK;
}

static qaws_status cx_end_point(qaws_exact_curve const* e, int end, cx_pt* out)
{
	qaws_exact_span const* sp = &e->spans[end ? e->span_count - 1 : 0];
	unsigned int o = end ? 3 * sp->degree : 0;
	out->x = sp->h[o];
	out->y = sp->h[o + 1];
	out->w = sp->h[o + 2];
	cx_reduce(&out->x, 3);
	return cx_pt_norm(out);
}

/* a line from the end of source a to the start of source b when they do not meet exactly */
static qaws_status cx_bridge(cx_ctx* x, unsigned int a, unsigned int b, unsigned int operand, unsigned int path, unsigned int index)
{
	cx_pt pa, pb;
	int c;
	double qa[2], qb[2];
	qaws_curve* line = NULL;
	qaws_status st;
	TRY(cx_end_point(x->src[a].e, 1, &pa));
	TRY(cx_end_point(x->src[b].e, 0, &pb));
	TRY(cx_pt_cmp(&pa, &pb, &c));
	if (c == 0)
		return QAWS_STATUS_OK;
	cx_pt_world(x, &pa, qa);
	cx_pt_world(x, &pb, qb);
	TRY(cx_line(x, qa, qb, &line));
	return cx_add_curve(x, line, operand, path, index, 0.0, 1.0);
}

static qaws_status cx_add_paths(cx_ctx* x, qaws_path_2d const* paths, unsigned int count, unsigned int operand)
{
	unsigned int i, j;
	qaws_status st;
	for (i = 0; i < count; i++)
	{
		qaws_path_2d const* p = &paths[i];
		unsigned int first = x->nsrc, last_end = CX_NONE;
		if (p->curve_count && !p->curves)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (j = 0; j < p->curve_count; j++)
		{
			unsigned int before = x->nsrc;
			if (!p->curves[j])
				return QAWS_STATUS_INVALID_ARGUMENT;
			if (p->curves[j]->dimension != QAWS_DIMENSION_2D)
				return QAWS_STATUS_INVALID_DIMENSION;
			TRY(cx_add_curve(x, p->curves[j], operand, i, j, 0.0, 1.0));
			/* pieces that do not meet exactly: the line between them */
			if (last_end != CX_NONE && x->nsrc > before)
			{
				unsigned int after = x->nsrc;
				TRY(cx_bridge(x, last_end, before, operand, i, j));
				if (x->nsrc > after)
				{
					/* keep the path order: the bridge goes before the new curve's pieces */
					cx_src br = x->src[after];
					memmove(&x->src[before + 1], &x->src[before], sizeof(cx_src) * (after - before));
					x->src[before] = br;
				}
			}
			if (x->nsrc > first)
				last_end = x->nsrc - 1;
		}
		/* closed: the line back to the start */
		if (x->nsrc > first)
			TRY(cx_bridge(x, x->nsrc - 1, first, operand, i, p->curve_count));
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  2-3. Hits and events                                                    */
/* ======================================================================== */

static qaws_status cx_event_add(cx_ctx* x, unsigned int src, double lo, double hi, unsigned int* out)
{
	cx_event* v;
	CX_GROW(x->ev, x->nev, x->capev, cx_event);
	v = &x->ev[x->nev];
	v->src = src;
	v->lo = lo;
	v->hi = hi;
	v->pt = v->rat = v->partner = v->vertex = CX_NONE;
	v->parent = x->nev;
	*out = x->nev++;
	return QAWS_STATUS_OK;
}

static qaws_status cx_pt_add(cx_ctx* x, cx_pt const* p, unsigned int* out)
{
	if (qaws_internal_add2(qaws_internal_mul2(2, x->npts), 2) > x->capptf)
	{
		unsigned int nc = qaws_internal_grow(x->capptf, qaws_internal_add2(qaws_internal_mul2(2, x->npts), 2), 256);
		double* g = nc ? (double*)qaws_internal_realloc(x->ptf, qaws_internal_mul2(sizeof(double), nc)) : NULL;
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->ptf = g;
		x->capptf = nc;
	}
	x->ptf[2 * x->npts] = qaws_exact_ratio_to_double(&p->x, &p->w);
	x->ptf[2 * x->npts + 1] = qaws_exact_ratio_to_double(&p->y, &p->w);
	CX_GROW(x->pts, x->npts, x->cappts, cx_pt);
	x->pts[x->npts] = *p;
	*out = x->npts++;
	return QAWS_STATUS_OK;
}

static qaws_status cx_rat_add(cx_ctx* x, unsigned int k, xint const* num, xint const* den, unsigned int* out)
{
	cx_rat* r;
	CX_GROW(x->rats, x->nrats, x->caprats, cx_rat);
	r = &x->rats[x->nrats];
	r->k = k;
	r->num = *num;
	r->den = *den;
	*out = x->nrats++;
	return QAWS_STATUS_OK;
}

/* an event at an exact double parameter: its rational parameter and exact point */
static qaws_status cx_exact_event(cx_ctx* x, unsigned int ei, int side)
{
	cx_event* v = &x->ev[ei];
	qaws_exact_curve const* e = x->src[v->src].e;
	unsigned int k;
	xint num, den;
	cx_pt p;
	qaws_status st;
	TRY(cx_local(e, v->lo, side, &k, &num, &den));
	TRY(cx_rat_add(x, k, &num, &den, &x->ev[ei].rat));
	/* at a span's end on the parameter lattice: that end's double, so the
	   order along the curve agrees (a double a few ulps past a curve's start
	   is its start on the lattice, not a vertex after it) */
	if (qaws_exact_int_is_zero(&num) || qaws_exact_int_cmp(&num, &den) == 0)
	{
		v = &x->ev[ei];
		v->lo = v->hi = ldexp((double)(qaws_exact_int_is_zero(&num) ? e->spans[k].a : e->spans[k].b), -e->param_shift);
	}
	TRY(cx_span_point(&e->spans[k], &num, &den, &p));
	return cx_pt_add(x, &p, &x->ev[ei].pt);
}

/* an event at a known exact point: its rational parameter on the source (inverted on the span holding the enclosure) */
static qaws_status cx_point_event(cx_ctx* x, unsigned int ei, cx_pt const* p)
{
	cx_event* v = &x->ev[ei];
	qaws_exact_curve const* e = x->src[v->src].e;
	unsigned int k, kk;
	xint num, den, P[3];
	int on = 0;
	qaws_status st;
	TRY(cx_pt_add(x, p, &x->ev[ei].pt));
	if (v->lo == v->hi)
		return cx_exact_event(x, ei, 0);
	TRY(cx_local(e, 0.5 * (v->lo + v->hi), 0, &k, &num, &den));
	P[0] = p->x;
	P[1] = p->y;
	P[2] = p->w;
	/* the span holding the middle, then its neighbours (the point at their common end) */
	for (kk = 0; kk < 3; kk++)
	{
		unsigned int use = kk == 0 ? k : (kk == 1 ? k - 1 : k + 1);
		if ((kk == 1 && k == 0) || use >= e->span_count)
			continue;
		st = qaws_exact_span_invert(&e->spans[use], 2, P, &num, &den, &on);
		if (st == QAWS_STATUS_OK && on && !qaws_exact_int_is_zero(&den) && qaws_exact_int_sign(&num) >= 0 && qaws_exact_int_cmp(&num, &den) <= 0)
			return cx_rat_add(x, use, &num, &den, &x->ev[ei].rat);
	}
	return QAWS_STATUS_OK;
}

/* the exact crossing point of two line sources */
static qaws_status cx_line_cross(cx_ctx* x, unsigned int a, unsigned int b, cx_pt* out, int* ok)
{
	qaws_exact_span const* A = &x->src[a].e->spans[0];
	qaws_exact_span const* B = &x->src[b].e->spans[0];
	xint la[3], lb[3], t1, t2;
	unsigned int k;
	qaws_status st;
	*ok = 0;
	if (x->src[a].e->span_count != 1 || x->src[b].e->span_count != 1 || A->degree != 1 || B->degree != 1)
		return QAWS_STATUS_OK;
	/* lines through the points (homogeneous cross products), the point where they meet */
	for (k = 0; k < 3; k++)
	{
		unsigned int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
		TRY(qaws_exact_int_mul(&t1, &A->h[k1], &A->h[3 + k2]));
		TRY(qaws_exact_int_mul(&t2, &A->h[k2], &A->h[3 + k1]));
		TRY(qaws_exact_int_sub(&la[k], &t1, &t2));
		TRY(qaws_exact_int_mul(&t1, &B->h[k1], &B->h[3 + k2]));
		TRY(qaws_exact_int_mul(&t2, &B->h[k2], &B->h[3 + k1]));
		TRY(qaws_exact_int_sub(&lb[k], &t1, &t2));
	}
	{
		xint* o[3];
		o[0] = &out->x; o[1] = &out->y; o[2] = &out->w;
		for (k = 0; k < 3; k++)
		{
			unsigned int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
			TRY(qaws_exact_int_mul(&t1, &la[k1], &lb[k2]));
			TRY(qaws_exact_int_mul(&t2, &la[k2], &lb[k1]));
			TRY(qaws_exact_int_sub(o[k], &t1, &t2));
		}
	}
	if (qaws_exact_int_is_zero(&out->w))
		return QAWS_STATUS_OK;   /* parallel */
	cx_reduce(&out->x, 3);
	TRY(cx_pt_norm(out));
	*ok = 1;
	return QAWS_STATUS_OK;
}

static int cx_is_knot(qaws_exact_curve const* e, double t, unsigned int* k, int* end)
{
	unsigned int i;
	for (i = 0; i < e->span_count; i++)
	{
		if (ldexp((double)e->spans[i].a, -e->param_shift) == t) { *k = i; *end = 0; return 1; }
		if (ldexp((double)e->spans[i].b, -e->param_shift) == t) { *k = i; *end = 1; return 1; }
	}
	return 0;
}

static qaws_status cx_span_end_point(qaws_exact_span const* sp, int end, cx_pt* out)
{
	unsigned int o = end ? 3 * sp->degree : 0;
	out->x = sp->h[o];
	out->y = sp->h[o + 1];
	out->w = sp->h[o + 2];
	cx_reduce(&out->x, 3);
	return cx_pt_norm(out);
}

/* the exact point at an overlap end: a knot of one of the two curves */
static qaws_status cx_overlap_end(cx_ctx* x, unsigned int sa, double ta, unsigned int sb, double tb, cx_pt* out, int* ok)
{
	unsigned int k;
	int end;
	*ok = 0;
	if (cx_is_knot(x->src[sa].e, ta, &k, &end))
	{
		*ok = 1;
		return cx_span_end_point(&x->src[sa].e->spans[k], end, out);
	}
	if (cx_is_knot(x->src[sb].e, tb, &k, &end))
	{
		*ok = 1;
		return cx_span_end_point(&x->src[sb].e->spans[k], end, out);
	}
	return QAWS_STATUS_OK;
}

/* room for the events, points and parameters the hits and the sources' ends will make, once */
static qaws_status cx_reserve(cx_ctx* x, unsigned int n)
{
	if (n > x->capev)
	{
		cx_event* g = (cx_event*)qaws_internal_realloc(x->ev, qaws_internal_mul2(sizeof(cx_event), n));
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->ev = g;
		x->capev = n;
	}
	if (n > x->cappts)
	{
		cx_pt* g = (cx_pt*)qaws_internal_realloc(x->pts, qaws_internal_mul2(sizeof(cx_pt), n));
		double* gf;
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->pts = g;
		x->cappts = n;
		gf = (double*)qaws_internal_realloc(x->ptf, qaws_internal_mul3(sizeof(double), 2, n));
		if (!gf) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->ptf = gf;
		x->capptf = 2 * n;
	}
	if (n > x->caprats)
	{
		cx_rat* g = (cx_rat*)qaws_internal_realloc(x->rats, qaws_internal_mul2(sizeof(cx_rat), n));
		if (!g) return QAWS_STATUS_ALLOCATION_FAILURE;
		x->rats = g;
		x->caprats = n;
	}
	return QAWS_STATUS_OK;
}

static qaws_status cx_hits(cx_ctx* x)
{
	qaws_exact_batch_desc bd;
	qaws_exact_curve const** cs;
	qaws_exact_batch_stats bs;
	unsigned int cap = 1024, n = 0, i;
	qaws_status st;
	cs = (qaws_exact_curve const**)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_exact_curve*), (x->nsrc + 1)));
	if (!cs)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->nsrc; i++)
		cs[i] = x->src[i].e;
	memset(&bd, 0, sizeof(bd));
	bd.curves = cs;
	bd.curve_count = x->nsrc;
	bd.flags = QAWS_EXACT_BATCH_SELF;
	bd.executor = x->d->executor;
	for (;;)
	{
		qaws_exact_batch_hit* g = (qaws_exact_batch_hit*)qaws_internal_realloc(x->hits, qaws_internal_mul2(sizeof(qaws_exact_batch_hit), cap));
		if (!g) { qaws_internal_free(cs); return QAWS_STATUS_ALLOCATION_FAILURE; }
		x->hits = g;
		st = qaws_exact_curve_batch_hits(x->hits, &n, &bs, cap, &bd);
		if (st != QAWS_STATUS_BUFFER_TOO_SMALL)
			break;
		cap = n + 16;
	}
	qaws_internal_free(cs);
	if (st != QAWS_STATUS_OK)
		return st;
	x->nhits = n;
	x->stats.hit_count = n;
	TRY(cx_reserve(x, 4 * n + 2 * x->nsrc + 16));
	/* events */
	for (i = 0; i < n; i++)
	{
		qaws_exact_batch_hit const* h = &x->hits[i];
		qaws_exact_pair const* p = &h->pair;
		if (p->kind == QAWS_EXACT_HIT_OVERLAP)
		{
			/* both ends, on both curves */
			unsigned int e;
			x->stats.overlap_count++;
			for (e = 0; e < 2; e++)
			{
				double ta = e ? p->a_hi : p->a_lo;
				double tb = (e != 0) != (p->opposite != 0) ? p->b_hi : p->b_lo;
				unsigned int ea, eb;
				cx_pt q;
				int ok = 0;
				TRY(cx_event_add(x, h->curve_a, ta, ta, &ea));
				TRY(cx_event_add(x, h->curve_b, tb, tb, &eb));
				TRY(cx_overlap_end(x, h->curve_a, ta, h->curve_b, tb, &q, &ok));
				if (!ok)
					return QAWS_STATUS_CERTIFICATION_FAILED;
				/* the ends rounded outward: an enclosure of a few ulps round the double */
				x->ev[ea].lo = nextafter(nextafter(ta, -HUGE_VAL), -HUGE_VAL);
				x->ev[ea].hi = nextafter(nextafter(ta, HUGE_VAL), HUGE_VAL);
				x->ev[eb].lo = nextafter(nextafter(tb, -HUGE_VAL), -HUGE_VAL);
				x->ev[eb].hi = nextafter(nextafter(tb, HUGE_VAL), HUGE_VAL);
				TRY(cx_point_event(x, ea, &q));
				TRY(cx_point_event(x, eb, &q));
				x->ev[ea].partner = eb;
				x->ev[eb].partner = ea;
				cx_union(x->ev, ea, eb);
			}
			continue;
		}
		{
			unsigned int ea, eb;
			cx_pt q;
			int ok = 0;
			if (p->multiplicity > 1 && p->multiplicity % 2 == 0)
				x->stats.tangency_count++;
			TRY(cx_event_add(x, h->curve_a, p->a_lo, p->a_hi, &ea));
			TRY(cx_event_add(x, h->curve_b, p->b_lo, p->b_hi, &eb));
			x->ev[ea].partner = eb;
			x->ev[eb].partner = ea;
			cx_union(x->ev, ea, eb);
			if (p->a_lo == p->a_hi)
			{
				TRY(cx_exact_event(x, ea, 0));
				q = x->pts[x->ev[ea].pt];
				TRY(cx_point_event(x, eb, &q));
			}
			else if (p->b_lo == p->b_hi)
			{
				TRY(cx_exact_event(x, eb, 0));
				q = x->pts[x->ev[eb].pt];
				TRY(cx_point_event(x, ea, &q));
			}
			else
			{
				TRY(cx_line_cross(x, h->curve_a, h->curve_b, &q, &ok));
				if (ok)
				{
					TRY(cx_point_event(x, ea, &q));
					TRY(cx_point_event(x, eb, &q));
				}
			}
		}
	}
	/* the ends of every source */
	for (i = 0; i < x->nsrc; i++)
	{
		unsigned int a, b;
		TRY(cx_event_add(x, i, x->src[i].t0, x->src[i].t0, &a));
		TRY(cx_exact_event(x, a, 1));
		TRY(cx_event_add(x, i, x->src[i].t1, x->src[i].t1, &b));
		TRY(cx_exact_event(x, b, -1));
	}
	return QAWS_STATUS_OK;
}

/* every source's events, by a counting sort */
static qaws_status cx_bucket(cx_ctx* x)
{
	unsigned int i, s;
	qaws_internal_free(x->sev_first);
	qaws_internal_free(x->sev);
	x->sev_first = (unsigned int*)qaws_internal_zalloc(x->nsrc + 2, sizeof(unsigned int));
	x->sev = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nev + 1)));
	if (!x->sev_first || !x->sev)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->nev; i++)
		x->sev_first[x->ev[i].src + 1]++;
	for (s = 0; s < x->nsrc; s++)
		x->sev_first[s + 1] += x->sev_first[s];
	for (i = 0; i < x->nev; i++)
		x->sev[x->sev_first[x->ev[i].src]++] = i;
	for (s = x->nsrc; s > 0; s--)
		x->sev_first[s] = x->sev_first[s - 1];
	x->sev_first[0] = 0;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  4. Vertices                                                             */
/* ======================================================================== */

/* points a and b (indices) ordered by x, then y: their doubles (correctly
   rounded) when they are clearly apart, exactly otherwise */
static qaws_status cx_pt_cmp_idx(cx_ctx const* x, unsigned int a, unsigned int b, int* out)
{
	double const* p = &x->ptf[2 * a];
	double const* q = &x->ptf[2 * b];
	double tx = 1e-12 * (fabs(p[0]) + fabs(q[0])) + 1e-300, ty = 1e-12 * (fabs(p[1]) + fabs(q[1])) + 1e-300;
	if (p[0] < q[0] - tx) { *out = -1; return QAWS_STATUS_OK; }
	if (p[0] > q[0] + tx) { *out = 1; return QAWS_STATUS_OK; }
	if (fabs(p[0] - q[0]) <= tx && 0)
		return QAWS_STATUS_OK;
	{
		qaws_status st;
		TRY(cx_pt_cmp1(&x->pts[a], &x->pts[b], 0, out));
		if (*out)
			return QAWS_STATUS_OK;
	}
	if (p[1] < q[1] - ty) { *out = -1; return QAWS_STATUS_OK; }
	if (p[1] > q[1] + ty) { *out = 1; return QAWS_STATUS_OK; }
	return cx_pt_cmp1(&x->pts[a], &x->pts[b], 1, out);
}

/* sort events by their exact points (merge sort; insertion below 32) */
static qaws_status cx_sort_by_point(cx_ctx* x, unsigned int* idx, unsigned int n)
{
	unsigned int i;
	int j, c;
	qaws_status st;
	if (n > 32)
	{
		unsigned int h = n / 2, a = 0, b = h, k = 0;
		unsigned int* tmp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), n));
		if (!tmp)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		st = cx_sort_by_point(x, idx, h);
		if (st == QAWS_STATUS_OK) st = cx_sort_by_point(x, idx + h, n - h);
		while (st == QAWS_STATUS_OK && (a < h || b < n))
		{
			if (a < h && b < n)
			{
				st = cx_pt_cmp_idx(x, x->ev[idx[a]].pt, x->ev[idx[b]].pt, &c);
				tmp[k++] = c <= 0 ? idx[a++] : idx[b++];
			}
			else
				tmp[k++] = a < h ? idx[a++] : idx[b++];
		}
		if (st == QAWS_STATUS_OK)
			memcpy(idx, tmp, sizeof(unsigned int) * n);
		qaws_internal_free(tmp);
		return st;
	}
	for (i = 1; i < n; i++)
	{
		unsigned int key = idx[i];
		for (j = (int)i - 1; j >= 0; j--)
		{
			TRY(cx_pt_cmp_idx(x, x->ev[idx[j]].pt, x->ev[key].pt, &c));
			if (c <= 0)
				break;
			idx[j + 1] = idx[j];
		}
		idx[j + 1] = key;
	}
	return QAWS_STATUS_OK;
}

/* the span of e holding the enclosure [lo, hi] (CX_NONE when it straddles a knot) */
static unsigned int cx_span_of(qaws_exact_curve const* e, double lo, double hi)
{
	unsigned int k;
	for (k = 0; k < e->span_count; k++)
	{
		double a = ldexp((double)e->spans[k].a, -e->param_shift), b = ldexp((double)e->spans[k].b, -e->param_shift);
		if (lo >= a && hi <= b)
			return k;
	}
	return CX_NONE;
}

/* Bernstein coefficients of the power-basis polynomial p on [s0, s1] (rational ends) */
static qaws_status cx_poly_restrict(qaws_exact_poly const* p, xint const* n0, xint const* d0, xint const* n1, xint const* d1, xint* out)
{
	unsigned int n = p->deg, i;
	xint* b = (xint*)qaws_internal_malloc(qaws_internal_mul3(sizeof(xint), 3, (n + 1)));
	xint* L = (xint*)qaws_internal_malloc(qaws_internal_mul3(sizeof(xint), 3, (n + 1)));
	xint* R = (xint*)qaws_internal_malloc(qaws_internal_mul3(sizeof(xint), 3, (n + 1)));
	xint sn, sd;
	qaws_status st = QAWS_STATUS_OK;
	if (!b || !L || !R)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
	{
		xint* bb = (xint*)qaws_internal_malloc(qaws_internal_mul2(sizeof(xint), (n + 1)));
		if (!bb)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			st = qaws_exact_poly_to_bernstein(p, bb);
			/* as homogeneous points (b, 0, 0) to reuse the split */
			for (i = 0; i <= n && st == QAWS_STATUS_OK; i++)
			{
				b[3 * i] = bb[i];
				qaws_exact_int_zero(&b[3 * i + 1]);
				qaws_exact_int_from_i64(&b[3 * i + 2], 1);
			}
			qaws_internal_free(bb);
		}
	}
	/* [0, s1], then its [s0 / s1, 1] */
	if (st == QAWS_STATUS_OK && n > 0)
		st = cx_split(b, n, n1, d1, L, NULL);
	if (st == QAWS_STATUS_OK && n > 0)
	{
		if (qaws_exact_int_is_zero(n1))
		{
			for (i = 0; i <= n; i++) out[i] = L[3 * i];
		}
		else
		{
			st = qaws_exact_int_mul(&sn, n0, d1);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&sd, d0, n1);
			if (st == QAWS_STATUS_OK) st = cx_split(L, n, &sn, &sd, NULL, R);
			for (i = 0; i <= n && st == QAWS_STATUS_OK; i++) out[i] = R[3 * i];
		}
	}
	else if (st == QAWS_STATUS_OK)
		out[0] = b[0];
	qaws_internal_free(b); qaws_internal_free(L); qaws_internal_free(R);
	return st;
}

static unsigned int cx_variations(xint const* b, unsigned int n)
{
	unsigned int i, v = 0;
	int last = 0;
	for (i = 0; i <= n; i++)
	{
		int s = qaws_exact_int_sign(&b[i]);
		if (!s) continue;
		if (last && s != last) v++;
		last = s;
	}
	return v;
}

/* the square-free part of p */
static qaws_status cx_squarefree(qaws_exact_poly const* p, qaws_exact_poly* out)
{
	qaws_exact_poly* w = (qaws_exact_poly*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_exact_poly), 2));
	qaws_status st;
	if (!w)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = qaws_exact_poly_derivative(p, &w[0]);
	if (st == QAWS_STATUS_OK)
	{
		qaws_exact_poly_trim(&w[0]);
		if (qaws_exact_poly_is_zero(&w[0]))
			*out = *p;
		else
		{
			st = qaws_exact_poly_gcd(p, &w[0], &w[1]);
			if (st == QAWS_STATUS_OK) st = qaws_exact_poly_div_exact(p, &w[1], out);
		}
	}
	qaws_internal_free(w);
	return st;
}

/* roots of the square-free q in the local enclosure (n0/d0, n1/d1): 0, 1, or 2 for "more" */
static qaws_status cx_roots_in(qaws_exact_poly const* q, xint const* n0, xint const* d0, xint const* n1, xint const* d1, unsigned int* out)
{
	xint* b = (xint*)qaws_internal_malloc(qaws_internal_mul2(sizeof(xint), (q->deg + 1)));
	qaws_status st;
	if (!b)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = cx_poly_restrict(q, n0, d0, n1, d1, b);
	if (st == QAWS_STATUS_OK)
	{
		unsigned int v = cx_variations(b, q->deg);
		/* an exact zero at an end counts as inside too */
		if (qaws_exact_int_is_zero(&b[0]) || qaws_exact_int_is_zero(&b[q->deg]))
			v = 2;
		*out = v > 2 ? 2 : v;
	}
	qaws_internal_free(b);
	return st;
}

/* does the partner curve of event ei pass through the exact point p (within the partner's enclosure)? */
static qaws_status cx_point_on_partner(cx_ctx* x, unsigned int ei, cx_pt const* p, int* out)
{
	cx_event const* pe = &x->ev[x->ev[ei].partner];
	qaws_exact_curve const* e = x->src[pe->src].e;
	unsigned int k = cx_span_of(e, pe->lo, pe->hi);
	xint P[3], num, den, n0, d0, n1, d1;
	unsigned int k0, k1;
	int on = 0, c0, c1;
	qaws_status st;
	*out = 0;
	if (k == CX_NONE)
		return QAWS_STATUS_CERTIFICATION_FAILED;
	P[0] = p->x; P[1] = p->y; P[2] = p->w;
	TRY(qaws_exact_span_invert(&e->spans[k], 2, P, &num, &den, &on));
	if (!on || qaws_exact_int_is_zero(&den))
		return QAWS_STATUS_OK;
	/* its parameter inside the partner's enclosure */
	TRY(cx_local(e, pe->lo, 1, &k0, &n0, &d0));
	TRY(cx_local(e, pe->hi, -1, &k1, &n1, &d1));
	if (k0 != k || k1 != k)
		return QAWS_STATUS_CERTIFICATION_FAILED;
	TRY(cx_ratio_cmp(&num, &den, &n0, &d0, &c0));
	TRY(cx_ratio_cmp(&num, &den, &n1, &d1, &c1));
	*out = c0 >= 0 && c1 <= 0;
	return QAWS_STATUS_OK;
}

/*
 * Two events on one source with overlapping enclosures: the same point?
 * *out 1 same, 0 apart; QAWS_STATUS_CERTIFICATION_FAILED when it cannot be
 * proven (or they are apart but cannot be ordered).
 */
static qaws_status cx_same(cx_ctx* x, unsigned int a, unsigned int b, int* out)
{
	cx_event const* ea = &x->ev[a];
	cx_event const* eb = &x->ev[b];
	qaws_status st;
	*out = 0;
	x->stats.identity_tests++;
	if (ea->pt != CX_NONE && eb->pt != CX_NONE)
	{
		/* different exact points: apart; ordered when both parameters are rational */
		if (ea->rat != CX_NONE && eb->rat != CX_NONE)
			return QAWS_STATUS_OK;
		return QAWS_STATUS_CERTIFICATION_FAILED;
	}
	if (ea->pt != CX_NONE || eb->pt != CX_NONE)
	{
		unsigned int alg = ea->pt != CX_NONE ? b : a, ex = ea->pt != CX_NONE ? a : b;
		if (x->ev[alg].partner == CX_NONE)
			return QAWS_STATUS_CERTIFICATION_FAILED;
		TRY(cx_point_on_partner(x, alg, &x->pts[x->ev[ex].pt], out));
		return *out ? QAWS_STATUS_OK : QAWS_STATUS_CERTIFICATION_FAILED;
	}
	/* two algebraic points: a common root of the partners' implicit equations on this span */
	{
		qaws_exact_curve const* e = x->src[ea->src].e;
		cx_event const* pa = &x->ev[ea->partner];
		cx_event const* pb = &x->ev[eb->partner];
		unsigned int k = cx_span_of(e, ea->lo < eb->lo ? ea->lo : eb->lo, ea->hi > eb->hi ? ea->hi : eb->hi);
		unsigned int ka = cx_span_of(x->src[pa->src].e, pa->lo, pa->hi), kb = cx_span_of(x->src[pb->src].e, pb->lo, pb->hi);
		qaws_exact_poly* w;
		xint n0, d0, n1, d1, m0, e0, m1, e1;
		unsigned int kk, ra = 0, rb = 0, rg = 0;
		double lo, hi;
		if (k == CX_NONE || ka == CX_NONE || kb == CX_NONE)
			return QAWS_STATUS_CERTIFICATION_FAILED;
		w = (qaws_exact_poly*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_exact_poly), 6));
		if (!w)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		st = qaws_exact_span_implicit_on(&x->src[pa->src].e->spans[ka], &e->spans[k], &w[0]);
		if (st == QAWS_STATUS_OK) st = qaws_exact_span_implicit_on(&x->src[pb->src].e->spans[kb], &e->spans[k], &w[1]);
		if (st == QAWS_STATUS_OK && (qaws_exact_poly_is_zero(&w[0]) || qaws_exact_poly_is_zero(&w[1])))
			st = QAWS_STATUS_CERTIFICATION_FAILED;
		if (st == QAWS_STATUS_OK) st = cx_squarefree(&w[0], &w[2]);
		if (st == QAWS_STATUS_OK) st = cx_squarefree(&w[1], &w[3]);
		/* each isolates its root in its own enclosure */
		if (st == QAWS_STATUS_OK) st = cx_local(e, ea->lo, 1, &kk, &n0, &d0);
		if (st == QAWS_STATUS_OK) st = cx_local(e, ea->hi, -1, &kk, &n1, &d1);
		if (st == QAWS_STATUS_OK) st = cx_roots_in(&w[2], &n0, &d0, &n1, &d1, &ra);
		if (st == QAWS_STATUS_OK) st = cx_local(e, eb->lo, 1, &kk, &m0, &e0);
		if (st == QAWS_STATUS_OK) st = cx_local(e, eb->hi, -1, &kk, &m1, &e1);
		if (st == QAWS_STATUS_OK) st = cx_roots_in(&w[3], &m0, &e0, &m1, &e1, &rb);
		if (st == QAWS_STATUS_OK && (ra != 1 || rb != 1))
			st = QAWS_STATUS_CERTIFICATION_FAILED;
		/* the common part of the enclosures, and the gcd's roots there */
		lo = ea->lo > eb->lo ? ea->lo : eb->lo;
		hi = ea->hi < eb->hi ? ea->hi : eb->hi;
		if (st == QAWS_STATUS_OK) st = qaws_exact_poly_gcd(&w[2], &w[3], &w[4]);
		if (st == QAWS_STATUS_OK)
		{
			qaws_exact_poly_trim(&w[4]);
			if (w[4].deg == 0)
				rg = 0;
			else
			{
				st = cx_local(e, lo, 1, &kk, &n0, &d0);
				if (st == QAWS_STATUS_OK) st = cx_local(e, hi, -1, &kk, &n1, &d1);
				if (st == QAWS_STATUS_OK) st = cx_roots_in(&w[4], &n0, &d0, &n1, &d1, &rg);
			}
		}
		qaws_internal_free(w);
		if (st != QAWS_STATUS_OK)
			return st;
		if (rg == 1)
		{
			*out = 1;
			return QAWS_STATUS_OK;
		}
		/* apart, but the enclosures overlap: no order without refining them */
		return QAWS_STATUS_CERTIFICATION_FAILED;
	}
}

static int cx_cmp_ev_lo(cx_ctx const* x, unsigned int a, unsigned int b)
{
	double ma = 0.5 * (x->ev[a].lo + x->ev[a].hi), mb = 0.5 * (x->ev[b].lo + x->ev[b].hi);
	return ma < mb ? -1 : (ma > mb ? 1 : 0);
}

/* rational order of two events on one source with rational parameters */
static qaws_status cx_rat_cmp(cx_ctx const* x, unsigned int a, unsigned int b, int* out)
{
	cx_rat const* ra = &x->rats[x->ev[a].rat];
	cx_rat const* rb = &x->rats[x->ev[b].rat];
	if (ra->k != rb->k)
	{
		/* different spans: the knot between, unless both sit on it */
		int za = qaws_exact_int_cmp(&ra->num, &ra->den) == 0, zb = qaws_exact_int_is_zero(&rb->num);
		if (ra->k + 1 == rb->k && za && zb) { *out = 0; return QAWS_STATUS_OK; }
		za = qaws_exact_int_cmp(&rb->num, &rb->den) == 0; zb = qaws_exact_int_is_zero(&ra->num);
		if (rb->k + 1 == ra->k && za && zb) { *out = 0; return QAWS_STATUS_OK; }
		*out = ra->k < rb->k ? -1 : 1;
		return QAWS_STATUS_OK;
	}
	return cx_ratio_cmp(&ra->num, &ra->den, &rb->num, &rb->den, out);
}

static qaws_status cx_vertices(cx_ctx* x)
{
	unsigned int* idx = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nev + 1)));
	unsigned int n = 0, i, s;
	qaws_status st = QAWS_STATUS_OK;
	if (!idx)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* equal exact points */
	for (i = 0; i < x->nev; i++)
		if (x->ev[i].pt != CX_NONE)
			idx[n++] = i;
	st = cx_sort_by_point(x, idx, n);
	for (i = 1; i < n && st == QAWS_STATUS_OK; i++)
	{
		int c = 0;
		st = cx_pt_cmp(&x->pts[x->ev[idx[i - 1]].pt], &x->pts[x->ev[idx[i]].pt], &c);
		if (st == QAWS_STATUS_OK && c == 0)
			cx_union(x->ev, idx[i - 1], idx[i]);
	}
	/* per source: events whose enclosures overlap */
	for (s = 0; s < x->nsrc && st == QAWS_STATUS_OK; s++)
	{
		unsigned int m = 0, a, b;
		for (i = x->sev_first[s]; i < x->sev_first[s + 1]; i++)
			idx[m++] = x->sev[i];
		/* by the middle of the enclosure */
		for (a = 1; a < m; a++)
		{
			unsigned int key = idx[a];
			for (b = a; b > 0 && cx_cmp_ev_lo(x, idx[b - 1], key) > 0; b--)
				idx[b] = idx[b - 1];
			idx[b] = key;
		}
		for (a = 0; a < m && st == QAWS_STATUS_OK; a++)
			for (b = a + 1; b < m && st == QAWS_STATUS_OK && x->ev[idx[b]].lo <= x->ev[idx[a]].hi; b++)
			{
				int same = 0;
				if (cx_find(x->ev, idx[a]) == cx_find(x->ev, idx[b]))
					continue;
				st = cx_same(x, idx[a], idx[b], &same);
				if (st == QAWS_STATUS_OK && same)
					cx_union(x->ev, idx[a], idx[b]);
			}
	}
	qaws_internal_free(idx);
	if (st != QAWS_STATUS_OK)
		return st;
	/* number the vertices */
	x->nv = 0;
	for (i = 0; i < x->nev; i++)
		x->ev[i].vertex = CX_NONE;
	for (i = 0; i < x->nev; i++)
	{
		unsigned int r = cx_find(x->ev, i);
		if (x->ev[r].vertex == CX_NONE)
			x->ev[r].vertex = x->nv++;
		x->ev[i].vertex = x->ev[r].vertex;
	}
	x->stats.vertex_count = x->nv;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  5. Edges                                                                */
/* ======================================================================== */

/* is event a strictly before event b along their source (proven)? */
static qaws_status cx_before(cx_ctx const* x, unsigned int a, unsigned int b, int* out)
{
	cx_event const* ea = &x->ev[a];
	cx_event const* eb = &x->ev[b];
	if (ea->hi < eb->lo) { *out = 1; return QAWS_STATUS_OK; }
	if (eb->hi < ea->lo) { *out = -1; return QAWS_STATUS_OK; }
	if (ea->rat != CX_NONE && eb->rat != CX_NONE)
	{
		int c = 0;
		qaws_status st = cx_rat_cmp(x, a, b, &c);
		*out = c < 0 ? 1 : (c > 0 ? -1 : 0);
		return st;
	}
	*out = 0;   /* not separated */
	return QAWS_STATUS_OK;
}

static qaws_status cx_edges(cx_ctx* x)
{
	unsigned int* idx = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nev + 1)));
	unsigned int s, i, capc = 0, cape = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (!idx)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	x->ncuts = 0;
	x->ne = 0;
	for (s = 0; s < x->nsrc && st == QAWS_STATUS_OK; s++)
	{
		unsigned int m = 0, a, b, first = x->ncuts;
		for (i = x->sev_first[s]; i < x->sev_first[s + 1]; i++)
			idx[m++] = x->sev[i];
		for (a = 1; a < m; a++)
		{
			unsigned int key = idx[a];
			for (b = a; b > 0 && cx_cmp_ev_lo(x, idx[b - 1], key) > 0; b--)
				idx[b] = idx[b - 1];
			idx[b] = key;
		}
		/* exact order where the middles do not decide it */
		for (a = 1; a < m && st == QAWS_STATUS_OK; a++)
		{
			int o = 0;
			st = cx_before(x, idx[a - 1], idx[a], &o);
			if (st == QAWS_STATUS_OK && o < 0)
			{
				unsigned int t = idx[a - 1];
				idx[a - 1] = idx[a];
				idx[a] = t;
				a = a > 1 ? a - 2 : 0;
			}
		}
		/* cuts: runs of one vertex at one parameter */
		for (a = 0; a < m && st == QAWS_STATUS_OK; a++)
		{
			cx_event const* v = &x->ev[idx[a]];
			int o = 1;
			if (x->ncuts > first && x->cuts[x->ncuts - 1].vertex == v->vertex)
			{
				st = cx_before(x, x->cuts[x->ncuts - 1].ev, idx[a], &o);
				if (st != QAWS_STATUS_OK)
					break;
				if (o == 0)
				{
					/* the same cut: keep the event with a rational parameter */
					if (x->ev[x->cuts[x->ncuts - 1].ev].rat == CX_NONE && v->rat != CX_NONE)
						x->cuts[x->ncuts - 1].ev = idx[a];
					continue;
				}
			}
			else if (x->ncuts > first)
			{
				st = cx_before(x, x->cuts[x->ncuts - 1].ev, idx[a], &o);
				if (st == QAWS_STATUS_OK && o != 1)
					st = QAWS_STATUS_CERTIFICATION_FAILED;   /* two vertices not apart along the curve */
				if (st != QAWS_STATUS_OK)
					break;
			}
			if (x->ncuts == capc)
			{
				unsigned int nc = qaws_internal_grow(capc, (size_t)x->ncuts + 1, 256);
				cx_cut* g = nc ? (cx_cut*)qaws_internal_realloc(x->cuts, qaws_internal_mul2(sizeof(cx_cut), nc)) : NULL;
				if (!g) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
				x->cuts = g;
				capc = nc;
			}
			x->cuts[x->ncuts].ev = idx[a];
			x->cuts[x->ncuts].vertex = v->vertex;
			x->cuts[x->ncuts].t = v->lo == v->hi ? v->lo : 0.5 * (v->lo + v->hi);
			x->ncuts++;
		}
		x->src[s].first_cut = first;
		x->src[s].cut_count = x->ncuts - first;
		x->src[s].first_edge = x->ne;
		for (a = first; a + 1 < x->ncuts && st == QAWS_STATUS_OK; a++)
		{
			cx_edge* e;
			if (x->ne == cape || x->ne >= 0x7FFFFFFFu)   /* 2 ne half-edges count in unsigned int */
			{
				unsigned int nc = x->ne >= 0x7FFFFFFFu ? 0 : qaws_internal_grow(cape, (size_t)x->ne + 1, 256);
				cx_edge* g = nc ? (cx_edge*)qaws_internal_realloc(x->e, qaws_internal_mul2(sizeof(cx_edge), nc)) : NULL;
				if (!g) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
				x->e = g;
				cape = nc;
			}
			e = &x->e[x->ne];
			e->src = s;
			e->c0 = a;
			e->c1 = a + 1;
			e->v0 = x->cuts[a].vertex;
			e->v1 = x->cuts[a + 1].vertex;
			e->ds = x->src[s].operand == 0 ? 1 : 0;
			e->dc = x->src[s].operand == 1 ? 1 : 0;
			e->parent = x->ne;
			e->flip = 0;
			e->live = CX_NONE;
			e->comp = CX_NONE;
			x->ne++;
		}
	}
	qaws_internal_free(idx);
	return st;
}

/* union-find over edges with the parity of their directions */
static unsigned int cx_efind(cx_edge* e, unsigned int i, int* flip)
{
	int f = 0;
	while (e[i].parent != i)
	{
		f ^= e[i].flip;
		i = e[i].parent;
	}
	*flip = f;
	return i;
}

static qaws_status cx_eunion(cx_edge* e, unsigned int a, unsigned int b, int opposite)
{
	int fa, fb;
	unsigned int ra = cx_efind(e, a, &fa), rb = cx_efind(e, b, &fb);
	int rel = fa ^ fb ^ (opposite ? 1 : 0);
	if (ra == rb)
		return rel ? QAWS_STATUS_CERTIFICATION_FAILED : QAWS_STATUS_OK;
	if (ra < rb)
	{
		e[rb].parent = ra;
		e[rb].flip = rel;
	}
	else
	{
		e[ra].parent = rb;
		e[ra].flip = rel;
	}
	return QAWS_STATUS_OK;
}

/* edges of source s whose middle lies within [lo, hi] */
static unsigned int cx_edges_in(cx_ctx const* x, unsigned int s, double lo, double hi, unsigned int* out, unsigned int cap)
{
	unsigned int k, n = 0;
	cx_src const* sr = &x->src[s];
	for (k = 0; k + 1 < sr->cut_count && n < cap; k++)
	{
		double m = 0.5 * (x->cuts[sr->first_cut + k].t + x->cuts[sr->first_cut + k + 1].t);
		if (m > lo && m < hi)
			out[n++] = sr->first_edge + k;
	}
	return n;
}

static qaws_status cx_merge(cx_ctx* x)
{
	unsigned int i, cap = 1024, nlive = 0;
	unsigned int* la = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), cap));
	unsigned int* lb = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), cap));
	qaws_status st = QAWS_STATUS_OK;
	if (!la || !lb)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->nhits && st == QAWS_STATUS_OK; i++)
	{
		qaws_exact_batch_hit const* h = &x->hits[i];
		unsigned int na, nb, a, b;
		if (h->pair.kind != QAWS_EXACT_HIT_OVERLAP)
			continue;
		na = cx_edges_in(x, h->curve_a, h->pair.a_lo, h->pair.a_hi, la, cap);
		nb = cx_edges_in(x, h->curve_b, h->pair.b_lo, h->pair.b_hi, lb, cap);
		if (na != nb)
			st = QAWS_STATUS_CERTIFICATION_FAILED;   /* the vertices along the overlap do not match */
		for (a = 0; a < na && st == QAWS_STATUS_OK; a++)
		{
			cx_edge const* ea = &x->e[la[a]];
			int found = 0;
			for (b = 0; b < nb && !found; b++)
			{
				cx_edge const* eb = &x->e[lb[b]];
				if (!h->pair.opposite && eb->v0 == ea->v0 && eb->v1 == ea->v1)
					found = 1;
				else if (h->pair.opposite && eb->v0 == ea->v1 && eb->v1 == ea->v0)
					found = 1;
				if (found)
					st = cx_eunion(x->e, la[a], lb[b], h->pair.opposite);
			}
			if (!found)
				st = QAWS_STATUS_CERTIFICATION_FAILED;
		}
	}
	qaws_internal_free(la);
	qaws_internal_free(lb);
	if (st != QAWS_STATUS_OK)
		return st;
	/* steps summed on the representative; cancelled ones dropped */
	for (i = 0; i < x->ne; i++)
	{
		int f;
		unsigned int r = cx_efind(x->e, i, &f);
		if (r != i)
		{
			x->e[r].ds += f ? -x->e[i].ds : x->e[i].ds;
			x->e[r].dc += f ? -x->e[i].dc : x->e[i].dc;
		}
	}
	for (i = 0; i < x->ne; i++)
		if (x->e[i].parent == i && (x->e[i].ds != 0 || x->e[i].dc != 0))
			x->e[i].live = nlive++;
	x->stats.edge_count = nlive;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  6. Order round a vertex                                                 */
/* ======================================================================== */

typedef struct cx_dir
{
	int exact;                  /* a single direction (and the curvature) */
	xint lo[2], hi[2];          /* the cone of directions, counter-clockwise lo -> hi */
	int has_kappa;
	int kappa_sign;             /* the signed curvature as sign and log2 of its size */
	double kappa_log;
} cx_dir;

static qaws_status cx_cross(xint const* u, xint const* v, int* out)
{
	xint l, r;
	qaws_status st;
	TRY(qaws_exact_int_mul(&l, &u[0], &v[1]));
	TRY(qaws_exact_int_mul(&r, &u[1], &v[0]));
	*out = qaws_exact_int_cmp(&l, &r);
	return QAWS_STATUS_OK;
}

static qaws_status cx_dot(xint const* u, xint const* v, int* out)
{
	xint l, r;
	qaws_status st;
	TRY(qaws_exact_int_mul(&l, &u[0], &v[0]));
	TRY(qaws_exact_int_mul(&r, &u[1], &v[1]));
	TRY(qaws_exact_int_add(&l, &l, &r));
	*out = qaws_exact_int_sign(&l);
	return QAWS_STATUS_OK;
}

/* V_ij = Q_j.xy Q_i.w - Q_i.xy Q_j.w (the direction from point i to point j, times w_i w_j) */
static qaws_status cx_diff(xint const* Q, unsigned int i, unsigned int j, xint* out)
{
	xint t1, t2;
	unsigned int c;
	qaws_status st;
	for (c = 0; c < 2; c++)
	{
		TRY(qaws_exact_int_mul(&t1, &Q[3 * j + c], &Q[3 * i + 2]));
		TRY(qaws_exact_int_mul(&t2, &Q[3 * i + c], &Q[3 * j + 2]));
		TRY(qaws_exact_int_sub(&out[c], &t1, &t2));
	}
	return QAWS_STATUS_OK;
}

static double cx_log2_abs(xint const* v)
{
	unsigned int b = qaws_exact_int_bits(v);
	xint t;
	if (qaws_exact_int_is_zero(v))
		return -HUGE_VAL;
	if (b <= 60)
		return log2(fabs(qaws_exact_int_to_double(v)));
	qaws_exact_int_shr(&t, v, b - 60);
	return log2(fabs(qaws_exact_int_to_double(&t))) + (double)(b - 60);
}

/* a vector's two components as doubles at a common scale 2^-shift */
static void cx_vec_double(xint const* v, double* out, int* shift)
{
	unsigned int b0 = qaws_exact_int_bits(&v[0]), b1 = qaws_exact_int_bits(&v[1]), b = b0 > b1 ? b0 : b1;
	xint t;
	*shift = b > 60 ? (int)(b - 60) : 0;
	qaws_exact_int_shr(&t, &v[0], (unsigned int)*shift);
	out[0] = qaws_exact_int_to_double(&t);
	qaws_exact_int_shr(&t, &v[1], (unsigned int)*shift);
	out[1] = qaws_exact_int_to_double(&t);
}

/*
 * The sub-curve leaving the vertex along the half-edge: Q[0..n], Q[0] at the
 * vertex. An exact parameter gives the exact sub-curve (forward: [s, 1],
 * backward: [0, s] reversed); an enclosure gives the sub-curve over it.
 */
static qaws_status cx_direction(cx_ctx* x, unsigned int ev, int forward, cx_dir* out)
{
	cx_event const* v = &x->ev[ev];
	qaws_exact_curve const* e = x->src[v->src].e;
	unsigned int k, n, i, j, nq;
	xint* Q;
	xint num, den;
	qaws_status st = QAWS_STATUS_OK;
	memset(out, 0, sizeof(*out));
	if (v->rat != CX_NONE)
	{
		cx_rat const* r = &x->rats[v->rat];
		k = r->k;
		/* the span on the half-edge's side of a knot */
		if (forward && qaws_exact_int_cmp(&r->num, &r->den) == 0 && k + 1 < e->span_count)
		{
			k++;
			qaws_exact_int_zero(&num);
			qaws_exact_int_from_i64(&den, 1);
		}
		else if (!forward && qaws_exact_int_is_zero(&r->num) && k > 0)
		{
			k--;
			qaws_exact_int_from_i64(&num, 1);
			qaws_exact_int_from_i64(&den, 1);
		}
		else
		{
			num = r->num;
			den = r->den;
		}
	}
	else
	{
		k = cx_span_of(e, v->lo, v->hi);
		if (k == CX_NONE)
			return QAWS_STATUS_CERTIFICATION_FAILED;
	}
	n = e->spans[k].degree;
	if (n == 1)
	{
		/* a straight span: one direction everywhere, no bending */
		st = cx_diff(e->spans[k].h, 0, 1, out->lo);
		if (st != QAWS_STATUS_OK)
			return st;
		if (qaws_exact_int_is_zero(&out->lo[0]) && qaws_exact_int_is_zero(&out->lo[1]))
			return QAWS_STATUS_CERTIFICATION_FAILED;
		if (!forward)
		{
			qaws_exact_int_neg(&out->lo[0], &out->lo[0]);
			qaws_exact_int_neg(&out->lo[1], &out->lo[1]);
		}
		out->hi[0] = out->lo[0];
		out->hi[1] = out->lo[1];
		out->exact = 1;
		out->has_kappa = 1;
		out->kappa_sign = 0;
		return QAWS_STATUS_OK;
	}
	nq = n + 1;
	Q = (xint*)qaws_internal_malloc(qaws_internal_mul4(sizeof(xint), 3, nq, 2));
	if (!Q)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (v->rat != CX_NONE)
	{
		if (forward)
			st = cx_split(e->spans[k].h, n, &num, &den, NULL, Q);
		else
		{
			xint* L = Q + 3 * nq;
			st = cx_split(e->spans[k].h, n, &num, &den, L, NULL);
			for (i = 0; i < nq && st == QAWS_STATUS_OK; i++)
				memcpy(&Q[3 * i], &L[3 * (n - i)], sizeof(xint) * 3);
		}
		if (st == QAWS_STATUS_OK)
			cx_reduce(Q, 3 * nq);
		/* the first point apart from Q0 gives the direction */
		for (j = 1; j <= n && st == QAWS_STATUS_OK; j++)
		{
			st = cx_diff(Q, 0, j, out->lo);
			if (st == QAWS_STATUS_OK && !(qaws_exact_int_is_zero(&out->lo[0]) && qaws_exact_int_is_zero(&out->lo[1])))
				break;
		}
		if (st == QAWS_STATUS_OK && j > n)
			st = QAWS_STATUS_CERTIFICATION_FAILED;   /* a degenerate span */
		if (st == QAWS_STATUS_OK)
		{
			out->exact = 1;
			out->hi[0] = out->lo[0];
			out->hi[1] = out->lo[1];
		}
		/* the curvature at Q0: A1 = n V01, A2 = n (n - 1) (V02 - 2 V01), kappa = cross(A1, A2) w0^2 / |A1|^3 */
		if (st == QAWS_STATUS_OK && j == 1)
		{
			xint a1[2], a2[2], v02[2], t, c;
			double d1[2];
			int sh;
			out->has_kappa = 1;
			a1[0] = out->lo[0];
			a1[1] = out->lo[1];
			if (n >= 2)
			{
				st = cx_diff(Q, 0, 2, v02);
				for (i = 0; i < 2 && st == QAWS_STATUS_OK; i++)
				{
					st = qaws_exact_int_add(&t, &a1[i], &a1[i]);
					if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&a2[i], &v02[i], &t);
				}
				/* cross(A1, A2) = n * n (n - 1) cross(V01, V02 - 2 V01) */
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&c, &a1[0], &a2[1]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&t, &a1[1], &a2[0]);
				if (st == QAWS_STATUS_OK) st = qaws_exact_int_sub(&c, &c, &t);
			}
			else
				qaws_exact_int_zero(&c);
			if (st == QAWS_STATUS_OK)
			{
				out->kappa_sign = qaws_exact_int_sign(&c);
				if (out->kappa_sign)
				{
					/* log2 |kappa| = log2 |c| + log2 (n^2 (n - 1)) + 2 log2 w0 - 3 (log2 n + log2 |V01|) */
					cx_vec_double(a1, d1, &sh);
					out->kappa_log = cx_log2_abs(&c) + log2((double)n * n * (n - 1)) + 2 * cx_log2_abs(&Q[2]) -
						3 * (log2((double)n) + log2(hypot(d1[0], d1[1])) + sh);
				}
			}
		}
	}
	else
	{
		/* the sub-curve over the enclosure: its tangent cone */
		xint n0, d0, n1, d1, sn, sd;
		xint* L = Q + 3 * nq;
		unsigned int kk;
		st = cx_local(e, v->lo, 1, &kk, &n0, &d0);
		if (st == QAWS_STATUS_OK) st = cx_local(e, v->hi, -1, &kk, &n1, &d1);
		if (st == QAWS_STATUS_OK) st = cx_split(e->spans[k].h, n, &n1, &d1, L, NULL);
		if (st == QAWS_STATUS_OK && !qaws_exact_int_is_zero(&n1))
		{
			st = qaws_exact_int_mul(&sn, &n0, &d1);
			if (st == QAWS_STATUS_OK) st = qaws_exact_int_mul(&sd, &d0, &n1);
			if (st == QAWS_STATUS_OK) st = cx_split(L, n, &sn, &sd, NULL, Q);
		}
		else if (st == QAWS_STATUS_OK)
			memcpy(Q, L, sizeof(xint) * 3 * nq);
		if (st == QAWS_STATUS_OK)
			cx_reduce(Q, 3 * nq);
		/* every V_ij, i < j, within an open half-plane: the extreme ones */
		{
			int have = 0;
			for (i = 0; i < nq && st == QAWS_STATUS_OK; i++)
				for (j = i + 1; j < nq && st == QAWS_STATUS_OK; j++)
				{
					xint vv[2];
					int c = 0, d = 0;
					st = cx_diff(Q, i, j, vv);
					if (st != QAWS_STATUS_OK || (qaws_exact_int_is_zero(&vv[0]) && qaws_exact_int_is_zero(&vv[1])))
						continue;
					if (!have)
					{
						out->lo[0] = out->hi[0] = vv[0];
						out->lo[1] = out->hi[1] = vv[1];
						have = 1;
						continue;
					}
					st = cx_cross(out->lo, vv, &c);
					if (st == QAWS_STATUS_OK && c < 0) { out->lo[0] = vv[0]; out->lo[1] = vv[1]; }
					if (st == QAWS_STATUS_OK) st = cx_cross(out->hi, vv, &d);
					if (st == QAWS_STATUS_OK && d > 0) { out->hi[0] = vv[0]; out->hi[1] = vv[1]; }
				}
			if (st == QAWS_STATUS_OK && !have)
				st = QAWS_STATUS_CERTIFICATION_FAILED;
			if (st == QAWS_STATUS_OK)
			{
				/* within a half-plane: lo -> hi turns left by less than a half turn */
				int c = 0, d = 0;
				st = cx_cross(out->lo, out->hi, &c);
				if (st == QAWS_STATUS_OK) st = cx_dot(out->lo, out->hi, &d);
				if (st == QAWS_STATUS_OK && (c < 0 || (c == 0 && d <= 0)))
					st = QAWS_STATUS_CERTIFICATION_FAILED;
			}
			if (st == QAWS_STATUS_OK && !forward)
			{
				/* backward: the cone turned half round */
				for (i = 0; i < 2; i++)
				{
					qaws_exact_int_neg(&out->lo[i], &out->lo[i]);
					qaws_exact_int_neg(&out->hi[i], &out->hi[i]);
				}
			}
		}
		/* a line's cone is its direction */
		if (st == QAWS_STATUS_OK && n == 1)
		{
			out->exact = 1;
			out->has_kappa = 1;
			out->kappa_sign = 0;
		}
	}
	qaws_internal_free(Q);
	return st;
}

/* in which half-turn from r is v: 0 for [0, pi), 1 for [pi, 2 pi) */
static qaws_status cx_half(xint const* r, xint const* v, int* out)
{
	int c, d;
	qaws_status st;
	TRY(cx_cross(r, v, &c));
	TRY(cx_dot(r, v, &d));
	*out = (c > 0 || (c == 0 && d > 0)) ? 0 : 1;
	return QAWS_STATUS_OK;
}

/* angle from r: u before v (-1), same (0), after (1) */
static qaws_status cx_angle_cmp(xint const* r, xint const* u, xint const* v, int* out)
{
	int hu, hv, c;
	qaws_status st;
	TRY(cx_half(r, u, &hu));
	TRY(cx_half(r, v, &hv));
	if (hu != hv)
	{
		*out = hu < hv ? -1 : 1;
		return QAWS_STATUS_OK;
	}
	TRY(cx_cross(u, v, &c));
	*out = c > 0 ? -1 : (c < 0 ? 1 : 0);
	return QAWS_STATUS_OK;
}

/* does the cone strictly contain direction r (so angles from r would wrap inside it)? */
static qaws_status cx_cone_holds(cx_dir const* d, xint const* r, int* out)
{
	int a, b;
	qaws_status st;
	*out = 0;
	if (d->exact)
		return QAWS_STATUS_OK;
	TRY(cx_cross(d->lo, r, &a));
	TRY(cx_cross(r, d->hi, &b));
	*out = a > 0 && b > 0;
	return QAWS_STATUS_OK;
}

/* sort a vertex's fan counter-clockwise; every neighbour proven apart */
static qaws_status cx_sort_fan(cx_ctx* x, unsigned int* fan, unsigned int n, unsigned int const* h_ev, int const* h_fwd)
{
	cx_dir* d;
	xint r[2];
	unsigned int i, j, tries;
	qaws_status st = QAWS_STATUS_OK;
	static int const refs[8][2] = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 3, 1 }, { -1, 3 }, { -3, -1 }, { 1, -3 } };
	if (n < 3)
		return QAWS_STATUS_OK;   /* two half-edges: either order is the same cycle */
	d = (cx_dir*)qaws_internal_malloc(qaws_internal_mul2(sizeof(cx_dir), n));
	if (!d)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
		st = cx_direction(x, h_ev[fan[i]], h_fwd[fan[i]], &d[i]);
	/* a reference direction inside no cone */
	for (tries = 0; tries < 8 && st == QAWS_STATUS_OK; tries++)
	{
		int held = 0;
		qaws_exact_int_from_i64(&r[0], refs[tries][0]);
		qaws_exact_int_from_i64(&r[1], refs[tries][1]);
		for (i = 0; i < n && !held && st == QAWS_STATUS_OK; i++)
			st = cx_cone_holds(&d[i], r, &held);
		if (!held)
			break;
	}
	if (st == QAWS_STATUS_OK && tries == 8)
		st = QAWS_STATUS_CERTIFICATION_FAILED;
	/* insertion sort by the cones' first directions, ties by curvature */
	for (i = 1; i < n && st == QAWS_STATUS_OK; i++)
	{
		unsigned int key = fan[i];
		cx_dir kd = d[i];
		for (j = i; j > 0 && st == QAWS_STATUS_OK; j--)
		{
			int c = 0;
			st = cx_angle_cmp(r, d[j - 1].lo, kd.lo, &c);
			if (st == QAWS_STATUS_OK && c == 0)
			{
				if (!d[j - 1].exact || !kd.exact || !d[j - 1].has_kappa || !kd.has_kappa)
				{
					st = QAWS_STATUS_CERTIFICATION_FAILED;   /* leaving tangent where only cones are known */
					break;
				}
				/* the one bending more to the left comes later */
				if (d[j - 1].kappa_sign != kd.kappa_sign)
					c = d[j - 1].kappa_sign < kd.kappa_sign ? -1 : 1;
				else if (d[j - 1].kappa_sign == 0 || fabs(d[j - 1].kappa_log - kd.kappa_log) < 1e-9)
				{
					st = QAWS_STATUS_CERTIFICATION_FAILED;   /* the same curvature: higher contact */
					break;
				}
				else
					c = ((d[j - 1].kappa_log < kd.kappa_log) == (kd.kappa_sign > 0)) ? -1 : 1;
			}
			if (st != QAWS_STATUS_OK || c <= 0)
				break;
			fan[j] = fan[j - 1];
			d[j] = d[j - 1];
		}
		fan[j] = key;
		d[j] = kd;
	}
	/* neighbours apart: the end of one cone before the start of the next */
	for (i = 0; i + 1 < n && st == QAWS_STATUS_OK; i++)
	{
		int c = 0;
		if (d[i].exact && d[i + 1].exact)
			continue;   /* single directions: the sort decided (ties by curvature) */
		st = cx_angle_cmp(r, d[i].hi, d[i + 1].lo, &c);
		if (st == QAWS_STATUS_OK && c >= 0)
			st = QAWS_STATUS_CERTIFICATION_FAILED;
	}
	qaws_internal_free(d);
	return st;
}

/* ======================================================================== */
/*  7-8. Faces, windings, output                                            */
/* ======================================================================== */

/* the sources a ray may meet, from the grid over their boxes (each once) */
typedef struct cx_raycol
{
	unsigned int* list;
	unsigned int n, cap;
	unsigned int* stamp;
	unsigned int gen;
	int failed;
} cx_raycol;

static int cx_raycol_visit(void* user, unsigned int const* items, unsigned int count, double t_exit)
{
	cx_raycol* r = (cx_raycol*)user;
	unsigned int k;
	(void)t_exit;
	for (k = 0; k < count; k++)
	{
		unsigned int s = items[k];
		if (r->stamp[s] == r->gen)
			continue;
		r->stamp[s] = r->gen;
		if (r->n == r->cap)
		{
			unsigned int nc = qaws_internal_grow(r->cap, (size_t)r->n + 1, 256);
			unsigned int* g = nc ? (unsigned int*)qaws_internal_realloc(r->list, qaws_internal_mul2(sizeof(unsigned int), nc)) : NULL;
			if (!g) { r->failed = 1; return 1; }
			r->list = g;
			r->cap = nc;
		}
		r->list[r->n++] = s;
	}
	return 0;
}

typedef struct cx_graph
{
	qaws_internal_arr_source* src;
	qaws_internal_arr_vertex* v;
	qaws_internal_arr_edge* e;
	unsigned int* emap;                 /* live edge -> cx edge */
	unsigned int* out;
	unsigned int* pos;
	unsigned int* next;
	unsigned int* cyc;
	unsigned int* face;
	unsigned int* h_ev;                 /* the event at a half-edge's origin */
	int* h_fwd;
	double* area;                       /* per cycle */
	unsigned int* cyc_comp;
	int* vexact;                        /* the vertex is exactly the double point vxy */
	double* vxy;
	unsigned int* vcomp;                /* each vertex's part */
	unsigned int* cvfirst;              /* part (root vertex) r's vertices: cvlist[cvfirst[r] .. cvfirst[r + 1]) */
	unsigned int* cvlist;
	qaws_bp_grid* sgrid;                /* the sources' boxes, for rays */
	cx_raycol ray;
	unsigned int ne, nv, ncyc;
} cx_graph;

static void cx_graph_free(cx_graph* g)
{
	qaws_internal_free(g->src); qaws_internal_free(g->v); qaws_internal_free(g->e); qaws_internal_free(g->emap); qaws_internal_free(g->out); qaws_internal_free(g->pos); qaws_internal_free(g->next); qaws_internal_free(g->cyc);
	qaws_internal_free(g->face); qaws_internal_free(g->h_ev); qaws_internal_free(g->h_fwd); qaws_internal_free(g->area); qaws_internal_free(g->cyc_comp); qaws_internal_free(g->vexact); qaws_internal_free(g->vxy); qaws_internal_free(g->cvfirst); qaws_internal_free(g->cvlist); qaws_internal_grid_destroy(g->sgrid); qaws_internal_free(g->ray.list); qaws_internal_free(g->ray.stamp);
}

static unsigned int cx_vfind(unsigned int* p, unsigned int i)
{
	while (p[i] != i)
	{
		p[i] = p[p[i]];
		i = p[i];
	}
	return i;
}

static double cx_half_area(cx_ctx const* x, cx_graph const* g, unsigned int h)
{
	qaws_internal_arr_edge const* e = &g->e[h >> 1];
	qaws_curve const* c = x->src[e->src].c;
	double a = 0.0;
	if (c->kind == QAWS_CURVE_KIND_BEZIER && c->degree == 1)
	{
		double const* p = g->v[(h & 1) ? e->v1 : e->v0].p;
		double const* q = g->v[(h & 1) ? e->v0 : e->v1].p;
		return 0.5 * (p[0] + q[0]) * (q[1] - p[1]);
	}
	qaws_internal_curve_area_2d(&a, c, (h & 1) ? e->t1 : e->t0, (h & 1) ? e->t0 : e->t1, x->ext * x->ext * 1e-12);
	return a;
}

/* are all of source s's edges that remain in part comp (so a ray from comp may pass it by)? */
static int cx_src_in(cx_ctx* x, unsigned int s, unsigned int comp)
{
	unsigned int k;
	cx_src const* sr = &x->src[s];
	for (k = 0; k + 1 < sr->cut_count; k++)
	{
		int f;
		unsigned int r = cx_efind(x->e, sr->first_edge + k, &f);
		if (x->e[r].live != CX_NONE && x->e[r].comp != comp)
			return 0;
	}
	return 1;
}

/* can the ray p + t dir (t >= 0) meet the box grown by margin? */
static int cx_ray_box(double const* p, double const* dir, double const* lo, double const* hi, double margin)
{
	double t0 = 0, t1 = HUGE_VAL;
	unsigned int c;
	for (c = 0; c < 2; c++)
	{
		double l = lo[c] - margin, h = hi[c] + margin;
		if (fabs(dir[c]) < 1e-300)
		{
			if (p[c] < l || p[c] > h)
				return 0;
			continue;
		}
		{
			double a = (l - p[c]) / dir[c], b = (h - p[c]) / dir[c];
			if (a > b) { double t = a; a = b; b = t; }
			if (a > t0) t0 = a;
			if (b < t1) t1 = b;
			if (t0 > t1)
				return 0;
		}
	}
	return 1;
}

/*
 * Where part `comp` lies among the others: the innermost bounded face of
 * another part holding it, found along a ray from an exact point of it.
 * The crossings are taken nearest first; one with the outside of another
 * part facing the point shows that part does not hold it (all of that part
 * is skipped from then on); the first one with a bounded face facing the
 * point gives the face; none: the unbounded face (0). A ray passing within
 * a margin of a vertex, touching a curve or with two crossings too close to
 * order is given up for the next one.
 */
typedef struct cx_xing
{
	double dist;
	unsigned int half;    /* the half-edge whose left side faces the ray's start */
} cx_xing;

static qaws_status cx_locate(cx_ctx* x, cx_graph* g, unsigned int comp, unsigned int const* cface, unsigned int* out_face)
{
	static double const dirs[4][2] = { { 1, 0.29296875 }, { -0.4140625, 1 }, { -1, -0.33984375 }, { 0.44921875, -1 } };
	unsigned int kv, capk = 64, capx = 64;
	qaws_exact_hit* hk = (qaws_exact_hit*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_exact_hit), capk));
	cx_xing* xs = (cx_xing*)qaws_internal_malloc(qaws_internal_mul2(sizeof(cx_xing), capx));
	unsigned char* skip = (unsigned char*)qaws_internal_zalloc(g->nv + 1, 1);
	double margin = 1e-9 * x->ext;
	qaws_status st = QAWS_STATUS_OK;
	if (!hk || !xs || !skip)
	{
		qaws_internal_free(hk); qaws_internal_free(xs); qaws_internal_free(skip);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	*out_face = 0;
	for (kv = g->cvfirst[comp]; kv < g->cvfirst[comp + 1] && st == QAWS_STATUS_OK; kv++)
	{
		unsigned int di, s, si, vi = g->cvlist[kv];
		double const* p = &g->vxy[2 * vi];
		if (g->v[vi].count == 0 || !g->vexact[vi])
			continue;
		for (di = 0; di < 4 && st == QAWS_STATUS_OK; di++)
		{
			double dir[2], p1[2];
			unsigned int nx = 0, a;
			int bad = 0;
			p1[0] = p[0] + dirs[di][0];
			p1[1] = p[1] + dirs[di][1];
			/* the direction the doubles give */
			dir[0] = p1[0] - p[0];
			dir[1] = p1[1] - p[1];
			/* the sources along the ray */
			g->ray.n = 0;
			g->ray.gen++;
			if (g->sgrid)
				qaws_internal_grid_ray(g->sgrid, p, dir, HUGE_VAL, cx_raycol_visit, &g->ray);
			else
				for (s = 0; s < x->nsrc; s++)
					cx_raycol_visit(&g->ray, &s, 1, 0);
			if (g->ray.failed)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (si = 0; si < g->ray.n && !bad && st == QAWS_STATUS_OK; si++)
			{
				cx_src const* sr = &x->src[s = g->ray.list[si]];
				unsigned int nh = 0, j;
				qaws_status ls;
				/* sources of the part itself, and those the ray cannot meet */
				if (cx_src_in(x, s, comp) || !cx_ray_box(p, dir, sr->lo, sr->hi, margin))
					continue;
				if (sr->is_line)
				{
					/* a straight source: the crossing in double (the ray is located with a margin anyway) */
					double ddx = sr->lb[0] - sr->la[0], ddy = sr->lb[1] - sr->la[1], ax = sr->la[0] - p[0], ay = sr->la[1] - p[1];
					double den = dir[0] * ddy - dir[1] * ddx, u;
					ls = QAWS_STATUS_OK;
					if (fabs(den) <= 1e-12 * hypot(dir[0], dir[1]) * hypot(ddx, ddy))
					{
						/* parallel: along the ray spoils it */
						if (fabs(ax * dir[1] - ay * dir[0]) <= margin * hypot(dir[0], dir[1]))
							ls = QAWS_STATUS_CERTIFICATION_FAILED;
					}
					else
					{
						u = (ax * dir[1] - ay * dir[0]) / den;
						if (u >= -1e-12 && u <= 1 + 1e-12)
						{
							double t = sr->t0 + u * (sr->t1 - sr->t0), w = 1e-12 * (sr->t1 - sr->t0);
							hk[0].kind = QAWS_EXACT_HIT_CROSSING;
							hk[0].t_lo = t - w;
							hk[0].t_hi = t + w;
							nh = 1;
						}
					}
				}
				else
				for (;;)
				{
					/* a coarse enclosure is enough to find the edge and the distance */
					ls = qaws_exact_curve_line_hits(hk, &nh, sr->e, p, p1, 1e-9 * (sr->t1 - sr->t0), capk);
					if (ls != QAWS_STATUS_BUFFER_TOO_SMALL)
						break;
					capk *= 4;
					{
						qaws_exact_hit* gk = (qaws_exact_hit*)qaws_internal_realloc(hk, qaws_internal_mul2(sizeof(qaws_exact_hit), capk));
						if (!gk) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
						hk = gk;
					}
				}
				if (st != QAWS_STATUS_OK)
					break;
				if (ls != QAWS_STATUS_OK)
				{
					bad = 1;   /* the ray touches a curve */
					break;
				}
				for (j = 0; j < nh && !bad; j++)
				{
					double tm = 0.5 * (hk[j].t_lo + hk[j].t_hi), pd[4], dist, left;
					unsigned int c, ei, r, live;
					int f, own = 0;
					if (hk[j].kind == QAWS_EXACT_HIT_OVERLAP) { bad = 1; break; }
					if (sr->is_line)
					{
						double u = (tm - sr->t0) / (sr->t1 - sr->t0);
						pd[0] = sr->la[0] + u * (sr->lb[0] - sr->la[0]);
						pd[1] = sr->la[1] + u * (sr->lb[1] - sr->la[1]);
						pd[2] = sr->lb[0] - sr->la[0];
						pd[3] = sr->lb[1] - sr->la[1];
					}
					else if (qaws_exact_curve_evaluate(pd, NULL, sr->e, tm, 1) != QAWS_STATUS_OK) { bad = 1; break; }
					dist = ((pd[0] - p[0]) * dir[0] + (pd[1] - p[1]) * dir[1]) / hypot(dir[0], dir[1]);
					if (dist < -margin)
						continue;
					/* the edge of s holding the hit, apart from its vertices */
					for (c = 0; c + 1 < sr->cut_count; c++)
					{
						cx_event const* ea = &x->ev[x->cuts[sr->first_cut + c].ev];
						cx_event const* eb = &x->ev[x->cuts[sr->first_cut + c + 1].ev];
						if (hk[j].t_hi >= ea->lo && hk[j].t_lo <= ea->hi)
						{
							/* at a vertex: the part's own one is passed, another one spoils the ray */
							if (cx_vfind(g->vcomp, x->cuts[sr->first_cut + c].vertex) != comp) bad = 1;
							else own = 1;
							break;
						}
						if (hk[j].t_lo > ea->hi && hk[j].t_hi < eb->lo)
							break;
					}
					if (bad)
						break;
					if (own)
						continue;
					if (c + 1 >= sr->cut_count)
					{
						bad = 1;   /* at the source's last vertex */
						break;
					}
					ei = sr->first_edge + c;
					r = cx_efind(x->e, ei, &f);
					live = x->e[r].live;
					if (live == CX_NONE || x->e[r].comp == comp)
						continue;   /* a cancelled edge, or the part itself */
					if (dist <= margin) { bad = 1; break; }
					left = pd[3] * dir[0] - pd[2] * dir[1];
					if (fabs(left) <= 1e-9 * hypot(pd[2], pd[3])) { bad = 1; break; }
					if (f) left = -left;
					if (nx == capx)
					{
						unsigned int nc = qaws_internal_grow(capx, (size_t)nx + 1, 64);
						cx_xing* gx = nc ? (cx_xing*)qaws_internal_realloc(xs, qaws_internal_mul2(sizeof(cx_xing), nc)) : NULL;
						if (!gx) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
						xs = gx;
						capx = nc;
					}
					xs[nx].dist = dist;
					xs[nx].half = left > 0 ? 2 * live : 2 * live + 1;
					nx++;
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
			if (!bad)
			{
				/* nearest first, every two apart */
				for (a = 1; a < nx; a++)
				{
					cx_xing key = xs[a];
					unsigned int b = a;
					while (b > 0 && xs[b - 1].dist > key.dist)
					{
						xs[b] = xs[b - 1];
						b--;
					}
					xs[b] = key;
				}
				/* sources merged into one edge (an overlap) cross the ray there once */
				{
					unsigned int m = 0;
					for (a = 0; a < nx; a++)
						if (m == 0 || (xs[a].half >> 1) != (xs[m - 1].half >> 1) || xs[a].dist - xs[m - 1].dist > margin)
							xs[m++] = xs[a];
					nx = m;
				}
				for (a = 1; a < nx && !bad; a++)
					if (xs[a].dist - xs[a - 1].dist <= margin)
						bad = 1;
			}
			if (!bad)
			{
				unsigned int face = 0;
				memset(skip, 0, g->nv + 1);
				for (a = 0; a < nx; a++)
				{
					unsigned int cy = g->cyc[xs[a].half], c2 = g->cyc_comp[cy];
					if (skip[c2])
						continue;
					if (cface[cy] != CX_NONE)
					{
						face = cface[cy];
						break;
					}
					skip[c2] = 1;   /* the outside of that part: it does not hold the point */
				}
				*out_face = face;
				qaws_internal_free(hk); qaws_internal_free(xs); qaws_internal_free(skip);
				return QAWS_STATUS_OK;
			}
			x->stats.ray_retries++;
		}
	}
	qaws_internal_free(hk); qaws_internal_free(xs); qaws_internal_free(skip);
	return st != QAWS_STATUS_OK ? st : QAWS_STATUS_CERTIFICATION_FAILED;
}

static int cx_inside(qaws_clip_type ct, qaws_fill_rule fr, int ws, int wc)
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

static qaws_status cx_build(cx_ctx* x, qaws_clip_result* result)
{
	cx_graph g;
	unsigned int i, k, nh, *cnt = NULL, *vcomp = NULL, nfaces = 1, *comp_outer = NULL, ncomp = 0;
	int *ws = NULL, *wc = NULL, *known = NULL, *inside = NULL;
	qaws_internal_arrangement a;
	qaws_status st = QAWS_STATUS_OK;
	memset(&g, 0, sizeof(g));
	g.ne = x->stats.edge_count;
	g.nv = x->nv;
	nh = 2 * g.ne;
	g.src = (qaws_internal_arr_source*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_internal_arr_source), (x->nsrc + 1)));
	g.v = (qaws_internal_arr_vertex*)qaws_internal_zalloc(g.nv + 1, sizeof(qaws_internal_arr_vertex));
	g.e = (qaws_internal_arr_edge*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_internal_arr_edge), (g.ne + 1)));
	g.emap = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.ne + 1)));
	g.out = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.pos = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.next = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.cyc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.face = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.h_ev = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	g.h_fwd = (int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(int), (nh + 1)));
	g.vexact = (int*)qaws_internal_zalloc(g.nv + 1, sizeof(int));
	g.vxy = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (g.nv + 1)));
	cnt = (unsigned int*)qaws_internal_zalloc(g.nv + 1, sizeof(unsigned int));
	vcomp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.nv + 1)));
	if (!g.src || !g.v || !g.e || !g.emap || !g.out || !g.pos || !g.next || !g.cyc || !g.face || !g.h_ev || !g.h_fwd || !g.vexact || !g.vxy || !cnt || !vcomp)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < x->nsrc; i++)
	{
		g.src[i].curve = x->src[i].c;
		g.src[i].operand = x->src[i].operand;
		g.src[i].path = x->src[i].path;
		g.src[i].index = x->src[i].index;
		g.src[i].base = x->src[i].base;
		g.src[i].scale = x->src[i].scale;
	}
	/* vertices: the doubles nearest their exact points (else the curve there) */
	for (i = 0; i < g.nv; i++)
	{
		g.v[i].sa = g.v[i].sb = CX_NONE;
		g.v[i].p[0] = g.v[i].p[1] = HUGE_VAL;
	}
	for (i = 0; i < x->nev; i++)
	{
		cx_event const* ev = &x->ev[i];
		qaws_internal_arr_vertex* v = &g.v[ev->vertex];
		double t = ev->lo == ev->hi ? ev->lo : 0.5 * (ev->lo + ev->hi);
		if (ev->pt != CX_NONE && v->p[0] == HUGE_VAL)
			cx_pt_world(x, &x->pts[ev->pt], v->p);
		if (ev->pt != CX_NONE && !g.vexact[ev->vertex])
			g.vexact[ev->vertex] = cx_pt_double(x, &x->pts[ev->pt], &g.vxy[2 * ev->vertex]);
		if (v->sa == CX_NONE)
		{
			v->sa = ev->src;
			v->ta = t;
		}
		else if (v->sb == CX_NONE && ev->src != v->sa)
		{
			v->sb = ev->src;
			v->tb = t;
		}
	}
	for (i = 0; i < g.nv; i++)
		if (g.v[i].p[0] == HUGE_VAL)
		{
			double pd[2] = { 0, 0 };
			if (g.v[i].sa != CX_NONE)
				qaws_exact_curve_evaluate(pd, NULL, x->src[g.v[i].sa].e, g.v[i].ta, 0);
			g.v[i].p[0] = pd[0];
			g.v[i].p[1] = pd[1];
		}
	/* live edges */
	for (i = 0; i < x->ne; i++)
	{
		cx_edge const* e = &x->e[i];
		qaws_internal_arr_edge* o;
		if (e->live == CX_NONE)
			continue;
		o = &g.e[e->live];
		o->src = e->src;
		o->t0 = x->cuts[e->c0].t;
		o->t1 = x->cuts[e->c1].t;
		o->v0 = e->v0;
		o->v1 = e->v1;
		o->ds = e->ds;
		o->dc = e->dc;
		g.emap[e->live] = i;
		g.h_ev[2 * e->live] = x->cuts[e->c0].ev;
		g.h_fwd[2 * e->live] = 1;
		g.h_ev[2 * e->live + 1] = x->cuts[e->c1].ev;
		g.h_fwd[2 * e->live + 1] = 0;
	}
	/* fans */
	for (i = 0; i < nh; i++)
		cnt[(i & 1) ? g.e[i >> 1].v1 : g.e[i >> 1].v0]++;
	for (i = 0, k = 0; i < g.nv; i++)
	{
		g.v[i].first = k;
		k += cnt[i];
		g.v[i].count = 0;
	}
	for (i = 0; i < nh; i++)
	{
		qaws_internal_arr_vertex* v = &g.v[(i & 1) ? g.e[i >> 1].v1 : g.e[i >> 1].v0];
		g.out[v->first + v->count++] = i;
	}
	for (i = 0; i < g.nv && st == QAWS_STATUS_OK; i++)
	{
		st = cx_sort_fan(x, &g.out[g.v[i].first], g.v[i].count, g.h_ev, g.h_fwd);
		for (k = 0; k < g.v[i].count; k++)
			g.pos[g.out[g.v[i].first + k]] = k;
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	/* next(h): clockwise of twin(h) round h's target */
	for (i = 0; i < nh; i++)
	{
		unsigned int tw = i ^ 1;
		qaws_internal_arr_vertex const* v = &g.v[(tw & 1) ? g.e[tw >> 1].v1 : g.e[tw >> 1].v0];
		g.next[i] = g.out[v->first + (g.pos[tw] + v->count - 1) % v->count];
	}
	/* parts (connected components) */
	for (i = 0; i < g.nv; i++) vcomp[i] = i;
	for (i = 0; i < g.ne; i++)
	{
		unsigned int p = cx_vfind(vcomp, g.e[i].v0), q = cx_vfind(vcomp, g.e[i].v1);
		if (p != q) vcomp[p] = q;
	}
	for (i = 0; i < g.ne; i++)
		x->e[g.emap[i]].comp = cx_vfind(vcomp, g.e[i].v0);
	g.vcomp = vcomp;
	/* each part's vertices */
	g.cvfirst = (unsigned int*)qaws_internal_zalloc(g.nv + 2, sizeof(unsigned int));
	g.cvlist = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.nv + 1)));
	if (!g.cvfirst || !g.cvlist) { st = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
	for (i = 0; i < g.nv; i++)
		g.cvfirst[cx_vfind(vcomp, i) + 1]++;
	for (i = 0; i < g.nv; i++)
		g.cvfirst[i + 1] += g.cvfirst[i];
	for (i = 0; i < g.nv; i++)
		g.cvlist[g.cvfirst[cx_vfind(vcomp, i)]++] = i;
	for (i = g.nv; i > 0; i--)
		g.cvfirst[i] = g.cvfirst[i - 1];
	g.cvfirst[0] = 0;
	/* the comp of merged / dropped edges: their representative's */
	for (i = 0; i < x->ne; i++)
	{
		int f;
		unsigned int r = cx_efind(x->e, i, &f);
		x->e[i].comp = x->e[r].live != CX_NONE ? x->e[r].comp : CX_NONE;
	}
	/* cycles and their areas */
	for (i = 0; i < nh; i++) g.cyc[i] = CX_NONE;
	g.ncyc = 0;
	for (i = 0; i < nh; i++)
	{
		unsigned int h = i;
		double area = 0;
		if (g.cyc[i] != CX_NONE) continue;
		do
		{
			g.cyc[h] = g.ncyc;
			area += cx_half_area(x, &g, h);
			h = g.next[h];
		} while (h != i && g.cyc[h] == CX_NONE);
		{
			double* ga = (double*)qaws_internal_realloc(g.area, qaws_internal_mul2(sizeof(double), (g.ncyc + 1)));
			unsigned int* gc = (unsigned int*)qaws_internal_realloc(g.cyc_comp, qaws_internal_mul2(sizeof(unsigned int), (g.ncyc + 1)));
			if (ga) g.area = ga;
			if (gc) g.cyc_comp = gc;
			if (!ga || !gc) { st = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
		}
		g.area[g.ncyc] = area;
		g.cyc_comp[g.ncyc] = x->e[g.emap[i >> 1]].comp;
		g.ncyc++;
	}
	/* each part's outer cycle: its smallest area (negative) */
	comp_outer = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.nv + 1)));
	if (!comp_outer) { st = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
	for (i = 0; i < g.nv; i++) comp_outer[i] = CX_NONE;
	for (i = 0; i < g.ncyc; i++)
	{
		unsigned int c = g.cyc_comp[i];
		if (comp_outer[c] == CX_NONE)
		{
			comp_outer[c] = i;
			ncomp++;
		}
		else if (g.area[i] < g.area[comp_outer[c]])
			comp_outer[c] = i;
	}
	x->stats.component_count = ncomp;
	{
		unsigned int* cface = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.ncyc + 1)));
		if (!cface) { st = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
		for (i = 0; i < g.ncyc; i++)
			cface[i] = comp_outer[g.cyc_comp[i]] == i ? CX_NONE : nfaces++;
		for (i = 0; i < g.ncyc; i++)
			if (cface[i] == CX_NONE && !(g.area[i] < 0))
				st = QAWS_STATUS_CERTIFICATION_FAILED;   /* no cycle turns clockwise round the part */
		for (i = 0; i < nh; i++)
			g.face[i] = cface[g.cyc[i]];
		/* the sources' boxes in a grid, for the rays */
		{
			qaws_bp_box* sb = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (x->nsrc + 1)));
			g.ray.stamp = (unsigned int*)qaws_internal_zalloc(x->nsrc + 1, sizeof(unsigned int));
			if (!sb || !g.ray.stamp)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; st == QAWS_STATUS_OK && i < x->nsrc; i++)
			{
				sb[i].lo[0] = x->src[i].lo[0]; sb[i].lo[1] = x->src[i].lo[1];
				sb[i].hi[0] = x->src[i].hi[0]; sb[i].hi[1] = x->src[i].hi[1];
				sb[i].lo[2] = sb[i].hi[2] = 0;
			}
			if (st == QAWS_STATUS_OK && qaws_internal_grid_create(&g.sgrid, sb, x->nsrc, 2) != QAWS_STATUS_OK)
				g.sgrid = NULL;
			qaws_internal_free(sb);
		}
		/* where each part lies among the others */
		{
			unsigned int* pface = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (g.nv + 1)));
			if (!pface)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (i = 0; i < g.ncyc && st == QAWS_STATUS_OK; i++)
				if (cface[i] == CX_NONE)
					st = cx_locate(x, &g, g.cyc_comp[i], cface, &pface[g.cyc_comp[i]]);
			for (i = 0; i < nh && st == QAWS_STATUS_OK; i++)
				if (cface[g.cyc[i]] == CX_NONE)
					g.face[i] = pface[g.cyc_comp[g.cyc[i]]];
			qaws_internal_free(pface);
		}
		qaws_internal_free(cface);
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	x->stats.face_count = nfaces;
	/* windings from the unbounded face across edges */
	ws = (int*)qaws_internal_zalloc(nfaces + 1, sizeof(int));
	wc = (int*)qaws_internal_zalloc(nfaces + 1, sizeof(int));
	known = (int*)qaws_internal_zalloc(nfaces + 1, sizeof(int));
	inside = (int*)qaws_internal_zalloc(nfaces + 1, sizeof(int));
	if (!ws || !wc || !known || !inside) { st = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
	{
		unsigned int changed = 1;
		known[0] = 1;
		while (changed)
		{
			changed = 0;
			for (i = 0; i < nh; i++)
			{
				unsigned int f = g.face[i], h2 = g.face[i ^ 1];
				qaws_internal_arr_edge const* e = &g.e[i >> 1];
				int sg = (i & 1) ? -1 : 1;
				if (known[f] && !known[h2])
				{
					ws[h2] = ws[f] - sg * e->ds;
					wc[h2] = wc[f] - sg * e->dc;
					known[h2] = 1;
					changed = 1;
				}
			}
		}
	}
	for (i = 0; i < nfaces; i++)
		inside[i] = known[i] && cx_inside(x->d->clip_type, x->d->fill_rule, ws[i], wc[i]);
	/* output */
	{
		qaws_clip_desc cd;
		memset(&cd, 0, sizeof(cd));
		cd.clip_type = x->d->clip_type;
		cd.fill_rule = x->d->fill_rule;
		cd.flags = x->d->flags;
		cd.z_fn = x->d->z_fn;
		cd.z_user = x->d->z_user;
		memset(&a, 0, sizeof(a));
		a.desc = &cd;
		a.sources = g.src;
		a.source_count = x->nsrc;
		a.vertices = g.v;
		a.vertex_count = g.nv;
		a.edges = g.e;
		a.edge_count = g.ne;
		a.out = g.out;
		a.half_face = g.face;
		a.face_count = nfaces;
		a.inside = inside;
		a.tolerance = x->ext * 1e-12;
		a.extent = x->ext;
		st = qaws_internal_clip_output_arrangement(&a, result);
	}
done:
	cx_graph_free(&g);
	qaws_internal_free(cnt); qaws_internal_free(vcomp); qaws_internal_free(comp_outer); qaws_internal_free(ws); qaws_internal_free(wc); qaws_internal_free(known); qaws_internal_free(inside);
	return st;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

static void cx_free(cx_ctx* x)
{
	unsigned int i;
	for (i = 0; i < x->nsrc; i++)
		qaws_exact_curve_destroy(x->src[i].e);
	for (i = 0; i < x->ntemp; i++)
		qaws_curve_destroy(x->temp[i]);
	qaws_internal_free(x->src); qaws_internal_free(x->temp); qaws_internal_free(x->pts); qaws_internal_free(x->ptf); qaws_internal_free(x->sev_first); qaws_internal_free(x->sev); qaws_internal_free(x->rats); qaws_internal_free(x->ev); qaws_internal_free(x->hits); qaws_internal_free(x->cuts); qaws_internal_free(x->e);
}

qaws_status qaws_clip_exact_execute(qaws_clip_result** out_result, qaws_clip_exact_stats* out_stats, qaws_clip_exact_desc const* desc)
{
	cx_ctx x;
	qaws_clip_result* r = NULL;
	unsigned int i;
	qaws_status st;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!desc || !out_result || (desc->subject_count && !desc->subjects) || (desc->clip_count && !desc->clips))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_result = NULL;
	memset(&x, 0, sizeof(x));
	x.d = desc;
	if (desc->exact)
		x.ed = *desc->exact;
	else
		qaws_exact_desc_default(&x.ed);
	st = cx_add_paths(&x, desc->subjects, desc->subject_count, 0);
	if (st == QAWS_STATUS_OK) st = cx_add_paths(&x, desc->clips, desc->clip_count, 1);
	x.stats.source_count = x.nsrc;
	/* boxes and extent */
	if (st == QAWS_STATUS_OK)
	{
		double lo[2] = { HUGE_VAL, HUGE_VAL }, hi[2] = { -HUGE_VAL, -HUGE_VAL };
		for (i = 0; i < x.nsrc; i++)
		{
			unsigned int k, c;
			cx_src* s = &x.src[i];
			s->lo[0] = s->lo[1] = HUGE_VAL;
			s->hi[0] = s->hi[1] = -HUGE_VAL;
			for (k = 0; k < s->e->span_count; k++)
			{
				double bl[3], bh[3];
				qaws_exact_span_box(&s->e->spans[k], 2, bl, bh);
				for (c = 0; c < 2; c++)
				{
					/* lattice units to world, and an unbounded box (a weight not positive) as everything */
					double l = bl[c] <= bh[c] ? ldexp(bl[c], x.ed.space_exp2) : -HUGE_VAL;
					double h = bl[c] <= bh[c] ? ldexp(bh[c], x.ed.space_exp2) : HUGE_VAL;
					if (l < s->lo[c]) s->lo[c] = l;
					if (h > s->hi[c]) s->hi[c] = h;
				}
			}
			for (k = 0; k < 2; k++)
			{
				if (s->lo[k] < lo[k]) lo[k] = s->lo[k];
				if (s->hi[k] > hi[k]) hi[k] = s->hi[k];
			}
		}
		x.ext = x.nsrc ? fmax(hi[0] - lo[0], hi[1] - lo[1]) : 1.0;
		if (!(x.ext > 0) || !isfinite(x.ext))
			x.ext = 1.0;
	}
	if (st == QAWS_STATUS_OK) st = qaws_internal_clip_result_empty(&r);
	if (st == QAWS_STATUS_OK && x.nsrc)
	{
		st = cx_hits(&x);
		if (st == QAWS_STATUS_OK) st = cx_bucket(&x);
		if (st == QAWS_STATUS_OK) st = cx_vertices(&x);
		if (st == QAWS_STATUS_OK) st = cx_edges(&x);
		if (st == QAWS_STATUS_OK) st = cx_merge(&x);
		if (st == QAWS_STATUS_OK) st = cx_build(&x, r);
	}
	if (out_stats)
		*out_stats = x.stats;
	cx_free(&x);
	if (st != QAWS_STATUS_OK)
	{
		qaws_clip_result_destroy(r);
		return st;
	}
	*out_result = r;
	return QAWS_STATUS_OK;
}
