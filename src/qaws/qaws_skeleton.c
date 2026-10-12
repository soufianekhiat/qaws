/*
 * Straight and curved skeletons (see qaws_skeleton.h).
 *
 *   1. region: the paths united under the fill rule (qaws_clip): outer loops
 *      counter-clockwise, holes clockwise, the interior on the left
 *   2. supports: every edge's moving carrier. A line is n . p = d + w t (n
 *      its inward normal); a circular arc is |p - m| = r0 - s t, s = +1
 *      when the region is inside the circle (it shrinks), -1 outside (it
 *      grows). Over (x, y, t) a line is a plane and a circle a cone.
 *   3. corners: the wavefront's vertices in circular lists. Between two
 *      lines, or two supports that stay tangent (a smooth joint), a corner
 *      moves at a constant velocity; otherwise along a conic, one of the two
 *      crossings of its supports (its branch).
 *   4. events, in a heap by time:
 *        - a corner meets a support (its neighbour's, or any, when reflex):
 *          three supports at one point. The differences of two cones are
 *          planes, so this is a linear system and one quadratic in t.
 *        - two supports touch (a circle with a line or a circle): squeeze
 *        - a reflex corner's angle reaches sigma, or its supports turn
 *          tangent: it splices, opening a circle that grows from its point
 *          (sigma = pi everywhere: the medial axis)
 *        - a circle shrinks to its centre
 *   5. the generalised intersection event at each point (Kelly and Wonka):
 *      the corners there and the edges through it make chains; a chain's
 *      inner edges end (an arc whose ends meet round the circle closes into a
 *      loop of its own: a bubble); chains sorted by angle are joined; corners
 *      between coincident opposite supports run along them (spikes)
 *   6. faces: each support's arcs walked from its edge's end back to its start
 */

#include "qaws_skeleton.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_arc.h"
#include "qaws_bezier.h"
#include "qaws_nurbs.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_kinds.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SK_NONE 0xFFFFFFFFu
#if QAWS_SCALAR_IS_FLOAT
#define SK_TOL_REL 1e-6          /* default tolerance, of the extent: the input's float rounding */
#else
#define SK_TOL_REL 1e-9          /* default tolerance, of the extent */
#endif
#define SK_PI 3.14159265358979323846

#define SK_LINE 0
#define SK_CIRCLE 1

typedef struct sk_sup
{
	int kind;
	double n[2], d;             /* line: n . p = d + w t */
	double m[2], r0;            /* circle: |p - m| = r0 - s t */
	int s;
	double w;                   /* speed: 1 for circles */
	unsigned int na, nb;        /* the input edge's nodes; a splice circle: its centre's, twice */
	double start[2];            /* the input edge's start, in full precision */
	double a0, sweep;           /* an input arc: start angle and signed sweep */
	unsigned int path, curve, span;
	int input;
} sk_sup;

#define SK_LINEAR 0              /* p + v (t - t0) */
#define SK_CONIC 1               /* a crossing of its supports */

/* a splitting corner's earliest meetings with the supports ahead of it,
   queued one at a time (all of them would fill the queue with stale ones) */
#define SK_CANDS 8
typedef struct sk_cand
{
	double t, p[2];
} sk_cand;

/* a support's earliest touchings with the others, kept and queued one at a
   time like a splitting corner's meetings (an event carries 1 + its index) */
typedef struct sk_touch
{
	double t, p[2];
	unsigned int b;
} sk_touch;

typedef struct sk_touches
{
	sk_touch v[SK_CANDS];
	unsigned int n, i;          /* how many, the one queued */
	unsigned int from, upto;    /* the supports looked at */
	int more;                   /* later ones not kept */
} sk_touches;

typedef struct sk_corner
{
	int kind;
	double p[2], t0, v[2];
	int branch;
	double tend;                /* conic: when its supports turn tangent */
	unsigned int in, out;       /* the supports of the edges before and after */
	unsigned int prev, next;
	unsigned int node;          /* where it started */
	int alive, smooth, fake;    /* fake: a loop of one circle, the corner a marker */
	int initial;                /* one of the input's corners */
	sk_cand cand[SK_CANDS];     /* the next meetings with supports ahead, by time */
	unsigned int ncand, icand;  /* how many, the one queued */
	int more;                   /* later ones not kept: found again when these run out */
	unsigned int greg, gedge;   /* its own and its out edge's grid entries now */
	double gexp;                /* until when they hold */
} sk_corner;

#define SK_EV_CORNER 0           /* the corner meets a support */
#define SK_EV_SPLICE 1           /* the reflex corner opens a circle */
#define SK_EV_SQUEEZE 2          /* supports a and b touch */
#define SK_EV_TANGENT 3          /* a convex corner's supports turn tangent */

typedef struct sk_event
{
	double t, p[2];
	unsigned int kind, corner, a, b;
	unsigned int seq;
} sk_event;

typedef struct sk_ilist
{
	unsigned int* v;
	unsigned int n, cap;
} sk_ilist;

typedef struct sk_gexp
{
	double t;
	unsigned int i;
} sk_gexp;

/* an arc's path */
#define SK_GEO_LINE 0
#define SK_GEO_CONIC 1
#define SK_GEO_CIRCLE 2

typedef struct sk_geo
{
	int kind;
	unsigned int in, out;       /* conic: the supports, the branch, the times */
	int branch;
	double t0, t1;
	double m[2], rho, a0, sweep;   /* circle: along it */
} sk_geo;

struct qaws_skeleton
{
	sk_sup* sup;                /* the input's supports first, then the splice circles */
	unsigned int nsup, ninput;
	unsigned int* loop_first;   /* loop l: supports loop_first[l] .. loop_first[l + 1], in order */
	unsigned int nloops;
	qaws_skeleton_node* nodes;
	unsigned int nnodes, capnodes;
	qaws_skeleton_arc* arcs;
	sk_geo* geo;
	unsigned int narcs, caparcs, capgeo;
	unsigned int* face_first;
	unsigned int* face_nodes;
	qaws_curve** edges;         /* per input support, its edge as a curve */
	double tol, tmax, sigma, ext;
};

typedef struct sk_run
{
	qaws_skeleton* s;
	int record;
	double stop, now, eps, eps_t;
	double eps_r, far;
	double rad;                 /* near a tangency: one point (a square root of eps_r) */
	double tbound;              /* no point is farther from the boundary than half the extent */
	sk_sup* sup;
	unsigned int nsup, capsup;
	sk_corner* c;
	unsigned int nc, capc;
	unsigned int* refs;         /* per support: alive corners whose out edge is on it */
	unsigned int caprefs;
	sk_touches* sq;             /* per support: its touchings ahead */
	unsigned int capsq;
	sk_event* heap;
	unsigned int nh, caph, seq;
	unsigned int* tag;
	unsigned int captag, stamp;
	unsigned int node_from;
	int merge;                  /* the event joins corners born at its point too */
	unsigned int touch[2];      /* a squeeze's supports (SK_NONE otherwise): both must have edges there */
	unsigned int depth;         /* events within events */
	int collect;                /* events gathered into coll, not queued */
	sk_event* coll;
	unsigned int ncoll, capcoll;
	unsigned int* found;        /* a query's corners */
	unsigned int nfound, capfound;
	/* the corners and edges between lines in a grid, by where they pass
	   until they have moved about a cell (then entered again); the others
	   always looked at */
	int gon;
	double glo[2], gh, grc, gre;
	unsigned int gn;            /* cells per side */
	sk_ilist* gc;               /* per cell: corners, as (index, entry) pairs */
	sk_ilist* ge;               /* per cell: edges by their first corner, as pairs */
	sk_ilist gac, gae;          /* the corners and edges on circles, as pairs */
	sk_gexp* gq;                /* when corners' entries run out */
	unsigned int ngq, capgq;
	qaws_status st;
} sk_run;

static int sk_grow(void** p, unsigned int* cap, unsigned int need, size_t size)
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

/* ------------------------------------------------------------------------ */
/*  Supports                                                                */
/* ------------------------------------------------------------------------ */

static double sk_rho(sk_sup const* e, double t)
{
	return e->r0 - e->s * t;
}

/* zero on the support at time t; > 0 on the unswept side */
static double sk_resid(sk_sup const* e, double const* p, double t)
{
	if (e->kind == SK_LINE)
		return e->n[0] * p[0] + e->n[1] * p[1] - e->d - e->w * t;
	return e->s * (e->r0 - hypot(p[0] - e->m[0], p[1] - e->m[1])) - t;
}

static int sk_cmp_uint(void const* a, void const* b)
{
	unsigned int x = *(unsigned int const*)a, y = *(unsigned int const*)b;
	return x < y ? -1 : x > y;
}

/* the edge's direction at p (the interior on its left) */
static void sk_dir_at(sk_sup const* e, double const* p, double* d)
{
	if (e->kind == SK_LINE)
	{
		d[0] = e->n[1];
		d[1] = -e->n[0];
	}
	else
	{
		double x = p[0] - e->m[0], y = p[1] - e->m[1], l = hypot(x, y);
		if (!(l > 0)) { d[0] = 1; d[1] = 0; return; }
		d[0] = -e->s * y / l;
		d[1] = e->s * x / l;
	}
}

/* the gradient of the support's equation, and its time derivative (for the
   corners' velocities) */
static void sk_grad(sk_sup const* e, double const* p, double t, double* g, double* gt)
{
	if (e->kind == SK_LINE)
	{
		g[0] = e->n[0]; g[1] = e->n[1];
		*gt = -e->w;
	}
	else
	{
		g[0] = 2 * (p[0] - e->m[0]); g[1] = 2 * (p[1] - e->m[1]);
		*gt = 2 * sk_rho(e, t) * e->s;
	}
}

/* the crossings of supports a and b at time t, in branch order: a line and
   a circle along the line; two circles by the side of the line from a's
   centre to b's. Returns 2 (1 for two lines); 0 when they do not meet,
   with pts[0] the nearest approach. */
static int sk_cross(sk_sup const* a, sk_sup const* b, double t, double eps, double* pts)
{
	if (a->kind == SK_LINE && b->kind == SK_LINE)
	{
		double det = a->n[0] * b->n[1] - a->n[1] * b->n[0], da = a->d + a->w * t, db = b->d + b->w * t;
		if (fabs(det) <= 1e-14)
			return 0;
		pts[0] = pts[2] = (da * b->n[1] - a->n[1] * db) / det;
		pts[1] = pts[3] = (a->n[0] * db - da * b->n[0]) / det;
		return 1;
	}
	if (a->kind != b->kind)
	{
		sk_sup const* l = a->kind == SK_LINE ? a : b;
		sk_sup const* c = a->kind == SK_LINE ? b : a;
		double rho = sk_rho(c, t), h = l->n[0] * c->m[0] + l->n[1] * c->m[1] - l->d - l->w * t;
		double k2 = rho * rho - h * h, k, f[2], dir[2];
		int ok = 1;
		f[0] = c->m[0] - h * l->n[0];
		f[1] = c->m[1] - h * l->n[1];
		dir[0] = l->n[1]; dir[1] = -l->n[0];
		if (k2 < 0)
		{
			ok = k2 > -2 * (fabs(rho) + eps) * eps;
			k2 = 0;
		}
		k = sqrt(k2);
		pts[0] = f[0] - k * dir[0]; pts[1] = f[1] - k * dir[1];
		pts[2] = f[0] + k * dir[0]; pts[3] = f[1] + k * dir[1];
		return ok ? 2 : 0;
	}
	{
		double dx = b->m[0] - a->m[0], dy = b->m[1] - a->m[1], dd = hypot(dx, dy);
		double ra = sk_rho(a, t), rb = sk_rho(b, t), x, k2, k, base[2], perp[2];
		int ok = 1;
		if (!(dd > 1e-300))
			return 0;
		x = (ra * ra - rb * rb + dd * dd) / (2 * dd);
		k2 = ra * ra - x * x;
		if (k2 < 0)
		{
			ok = k2 > -2 * (fabs(ra) + eps) * eps;
			k2 = 0;
		}
		k = sqrt(k2);
		base[0] = a->m[0] + x * dx / dd; base[1] = a->m[1] + x * dy / dd;
		perp[0] = -dy / dd; perp[1] = dx / dd;
		pts[0] = base[0] - k * perp[0]; pts[1] = base[1] - k * perp[1];
		pts[2] = base[0] + k * perp[0]; pts[3] = base[1] + k * perp[1];
		return ok ? 2 : 0;
	}
}

/* the times after t at which supports a and b turn tangent (up to 3) */
static unsigned int sk_tangencies(sk_sup const* a, sk_sup const* b, double* out)
{
	unsigned int n = 0;
	if (a->kind == SK_LINE && b->kind == SK_LINE)
		return 0;
	if (a->kind != b->kind)
	{
		sk_sup const* l = a->kind == SK_LINE ? a : b;
		sk_sup const* c = a->kind == SK_LINE ? b : a;
		double h0 = l->n[0] * c->m[0] + l->n[1] * c->m[1] - l->d;
		/* rho = r0 - s t against h = h0 - w t, and against -h */
		if (fabs(c->s - l->w) > 1e-15) out[n++] = (c->r0 - h0) / (c->s - l->w);
		if (fabs(c->s + l->w) > 1e-15) out[n++] = (c->r0 + h0) / (c->s + l->w);
		return n;
	}
	{
		double dd = hypot(b->m[0] - a->m[0], b->m[1] - a->m[1]);
		if (a->s + b->s != 0) out[n++] = (a->r0 + b->r0 - dd) / (a->s + b->s);
		if (a->s - b->s != 0)
		{
			out[n++] = (a->r0 - b->r0 - dd) / (a->s - b->s);
			out[n++] = (b->r0 - a->r0 - dd) / (b->s - a->s);
		}
		return n;
	}
}

/* supports tangent at every time (a smooth joint that stays smooth) */
static int sk_persistent(sk_sup const* a, sk_sup const* b, double eps)
{
	if (a->kind == SK_LINE && b->kind == SK_LINE)
		return 0;
	if (a->kind != b->kind)
	{
		sk_sup const* l = a->kind == SK_LINE ? a : b;
		sk_sup const* c = a->kind == SK_LINE ? b : a;
		double h0 = l->n[0] * c->m[0] + l->n[1] * c->m[1] - l->d;
		return (fabs(c->s - l->w) <= 1e-12 && fabs(c->r0 - h0) <= eps) || (fabs(c->s + l->w) <= 1e-12 && fabs(c->r0 + h0) <= eps);
	}
	{
		double dd = hypot(b->m[0] - a->m[0], b->m[1] - a->m[1]);
		if (dd <= eps)
			return 0;
		if (a->s + b->s == 0 && fabs(a->r0 + b->r0 - dd) <= eps) return 1;
		if (a->s == b->s && (fabs(a->r0 - b->r0 - dd) <= eps || fabs(b->r0 - a->r0 - dd) <= eps)) return 1;
		return 0;
	}
}

/* coincident supports with opposite sides at time t: a spike's two edges */
static int sk_opposite(sk_sup const* a, sk_sup const* b, double t, double eps)
{
	if (a->kind != b->kind)
		return 0;
	if (a->kind == SK_LINE)
		return fabs(a->n[0] * b->n[1] - a->n[1] * b->n[0]) <= 1e-12 && a->n[0] * b->n[0] + a->n[1] * b->n[1] < 0;
	return a->s == -b->s && hypot(a->m[0] - b->m[0], a->m[1] - b->m[1]) <= eps && fabs(sk_rho(a, t) - sk_rho(b, t)) <= eps;
}

static int sk_solve_quad(double a, double b, double c, double* x)
{
	double disc, q;
	if (fabs(a) <= 1e-14 * (fabs(b) + fabs(c)))
	{
		if (fabs(b) <= 1e-300)
			return 0;
		x[0] = -c / b;
		return 1;
	}
	disc = b * b - 4 * a * c;
	if (disc < 0)
	{
		if (disc > -1e-12 * b * b)
			disc = 0;
		else
			return 0;
	}
	q = -0.5 * (b + (b >= 0 ? sqrt(disc) : -sqrt(disc)));
	if (q == 0)
	{
		x[0] = 0;
		return 1;
	}
	x[0] = q / a;
	x[1] = c / q;
	return 2;
}

/* the points (x, y, t) where three supports meet: up to 2 */
static unsigned int sk_meet3(sk_sup const* const* e, double* p, double* t)
{
	unsigned int i, nc = 0, k = 0, n;
	double row[2][4], kq[5];
	sk_sup const* o[2];
	for (i = 0; i < 3; i++) if (e[i]->kind == SK_CIRCLE) nc++;
	if (nc == 0)
	{
		double m[3][4], det, dx, dy, dt;
		for (i = 0; i < 3; i++)
		{
			m[i][0] = e[i]->n[0]; m[i][1] = e[i]->n[1]; m[i][2] = -e[i]->w; m[i][3] = e[i]->d;
		}
		det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
			+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
		if (!(fabs(det) > 1e-14))
			return 0;
		dx = m[0][3] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][3] * m[2][2] - m[1][2] * m[2][3])
			+ m[0][2] * (m[1][3] * m[2][1] - m[1][1] * m[2][3]);
		dy = m[0][0] * (m[1][3] * m[2][2] - m[1][2] * m[2][3]) - m[0][3] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
			+ m[0][2] * (m[1][0] * m[2][3] - m[1][3] * m[2][0]);
		dt = m[0][0] * (m[1][1] * m[2][3] - m[1][3] * m[2][1]) - m[0][1] * (m[1][0] * m[2][3] - m[1][3] * m[2][0])
			+ m[0][3] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
		p[0] = dx / det; p[1] = dy / det; t[0] = dt / det;
		return 1;
	}
	/* K, a cone: x^2 + y^2 + B x + C y - t^2 + E t + F = 0 */
	for (i = 0; i < 3; i++) if (e[i]->kind == SK_CIRCLE) { k = i; break; }
	kq[0] = -2 * e[k]->m[0]; kq[1] = -2 * e[k]->m[1]; kq[2] = 2 * e[k]->r0 * e[k]->s;
	kq[3] = e[k]->m[0] * e[k]->m[0] + e[k]->m[1] * e[k]->m[1] - e[k]->r0 * e[k]->r0;
	for (i = 0, n = 0; i < 3; i++)
		if (i != k) o[n++] = e[i];
	/* the others as planes a x + b y + c t + f = 0 (a cone minus K) */
	for (i = 0; i < 2; i++)
	{
		if (o[i]->kind == SK_LINE)
		{
			row[i][0] = o[i]->n[0]; row[i][1] = o[i]->n[1]; row[i][2] = -o[i]->w; row[i][3] = -o[i]->d;
		}
		else
		{
			double ob = -2 * o[i]->m[0], oc = -2 * o[i]->m[1], oe = 2 * o[i]->r0 * o[i]->s;
			double of = o[i]->m[0] * o[i]->m[0] + o[i]->m[1] * o[i]->m[1] - o[i]->r0 * o[i]->r0;
			row[i][0] = kq[0] - ob; row[i][1] = kq[1] - oc; row[i][2] = kq[2] - oe; row[i][3] = kq[3] - of;
		}
	}
	{
		double det = row[0][0] * row[1][1] - row[0][1] * row[1][0];
		double scale = fabs(row[0][0]) + fabs(row[0][1]) + fabs(row[1][0]) + fabs(row[1][1]);
		if (fabs(det) > 1e-12 * scale * scale)
		{
			/* x = x0 + x1 t, y = y0 + y1 t */
			double x0 = (-row[0][3] * row[1][1] + row[0][1] * row[1][3]) / det, x1 = (-row[0][2] * row[1][1] + row[0][1] * row[1][2]) / det;
			double y0 = (-row[0][0] * row[1][3] + row[0][3] * row[1][0]) / det, y1 = (-row[0][0] * row[1][2] + row[0][2] * row[1][0]) / det;
			double qa = x1 * x1 + y1 * y1 - 1, qb = 2 * (x0 * x1 + y0 * y1) + kq[0] * x1 + kq[1] * y1 + kq[2];
			double qc = x0 * x0 + y0 * y0 + kq[0] * x0 + kq[1] * y0 + kq[3], r[2];
			unsigned int nr = (unsigned int)sk_solve_quad(qa, qb, qc, r), j;
			for (j = 0; j < nr; j++)
			{
				t[j] = r[j];
				p[2 * j] = x0 + x1 * r[j];
				p[2 * j + 1] = y0 + y1 * r[j];
			}
			return nr;
		}
		/* parallel rows: their combination fixes t, then a line meets K */
		{
			unsigned int a = fabs(row[0][0]) + fabs(row[0][1]) >= fabs(row[1][0]) + fabs(row[1][1]) ? 0 : 1, bb = 1 - a, j, nr = 0;
			double mu, ct, cf, tt, la, lb, lc, nn, h, rr, k2, kk, f[2];
			nn = row[a][0] * row[a][0] + row[a][1] * row[a][1];
			if (!(nn > 0))
				return 0;
			mu = (row[bb][0] * row[a][0] + row[bb][1] * row[a][1]) / nn;
			ct = row[bb][2] - mu * row[a][2];
			cf = row[bb][3] - mu * row[a][3];
			if (!(fabs(ct) > 1e-14))
				return 0;
			tt = -cf / ct;
			/* the line la x + lb y + lc = 0 against the circle at tt */
			la = row[a][0]; lb = row[a][1]; lc = row[a][2] * tt + row[a][3];
			nn = sqrt(nn);
			la /= nn; lb /= nn; lc /= nn;
			rr = e[k]->r0 - e[k]->s * tt;
			h = la * e[k]->m[0] + lb * e[k]->m[1] + lc;
			k2 = rr * rr - h * h;
			if (k2 < 0)
				return 0;
			kk = sqrt(k2);
			f[0] = e[k]->m[0] - h * la; f[1] = e[k]->m[1] - h * lb;
			for (j = 0; j < 2; j++)
			{
				double sg = j ? 1 : -1;
				t[nr] = tt;
				p[2 * nr] = f[0] + sg * kk * lb;
				p[2 * nr + 1] = f[1] - sg * kk * la;
				nr++;
			}
			return nr;
		}
	}
}

/* the angle along a circle from a to b in its travel direction, in [0, 2 pi) */
static double sk_angle_along(sk_sup const* e, double const* a, double const* b)
{
	double aa = atan2(a[1] - e->m[1], a[0] - e->m[0]), ab = atan2(b[1] - e->m[1], b[0] - e->m[0]), d = (ab - aa) * e->s;
	while (d < 0) d += 2 * SK_PI;
	while (d >= 2 * SK_PI) d -= 2 * SK_PI;
	return d;
}

/* the distance along support e from a to b in its travel direction */
static double sk_along(sk_sup const* e, double const* a, double const* b, double t)
{
	if (e->kind == SK_LINE)
		return (b[0] - a[0]) * e->n[1] - (b[1] - a[1]) * e->n[0];
	return sk_angle_along(e, a, b) * fabs(sk_rho(e, t));
}

/* ------------------------------------------------------------------------ */
/*  Corners                                                                 */
/* ------------------------------------------------------------------------ */

static void sk_pos(sk_run const* r, sk_corner const* c, double t, double* p)
{
	double pts[4];
	if (c->kind == SK_LINEAR)
	{
		p[0] = c->p[0] + c->v[0] * (t - c->t0);
		p[1] = c->p[1] + c->v[1] * (t - c->t0);
		return;
	}
	sk_cross(&r->sup[c->in], &r->sup[c->out], t, r->eps, pts);
	if (t >= c->tend - 1e3 * r->eps_t)
	{
		/* at its supports' tangency the crossing moves like a square root of
		   the time: the tangent point */
		p[0] = 0.5 * (pts[0] + pts[2]);
		p[1] = 0.5 * (pts[1] + pts[3]);
		return;
	}
	p[0] = pts[2 * c->branch];
	p[1] = pts[2 * c->branch + 1];
}

/* the corner's velocity (implicitly, on its two supports) */
static void sk_vel(sk_run const* r, sk_corner const* c, double t, double* v)
{
	double p[2], ga[2], gb[2], ta, tb, det;
	if (c->kind == SK_LINEAR)
	{
		v[0] = c->v[0]; v[1] = c->v[1];
		return;
	}
	sk_pos(r, c, t, p);
	sk_grad(&r->sup[c->in], p, t, ga, &ta);
	sk_grad(&r->sup[c->out], p, t, gb, &tb);
	det = ga[0] * gb[1] - ga[1] * gb[0];
	if (fabs(det) <= 1e-12 * hypot(ga[0], ga[1]) * hypot(gb[0], gb[1]))
	{
		/* tangent: along the common tangent, fast */
		double l = hypot(ga[0], ga[1]);
		v[0] = -ga[1] / (l > 0 ? l : 1) * 1e6;
		v[1] = ga[0] / (l > 0 ? l : 1) * 1e6;
		return;
	}
	v[0] = (-ta * gb[1] + ga[1] * tb) / det;
	v[1] = (-ga[0] * tb + ta * gb[0]) / det;
}

/* the interior angle at the corner, in (0, 2 pi) */
static double sk_angle(sk_run const* r, sk_corner const* c, double t)
{
	double p[2], din[2], dout[2];
	sk_pos(r, c, t, p);
	sk_dir_at(&r->sup[c->in], p, din);
	sk_dir_at(&r->sup[c->out], p, dout);
	return SK_PI - atan2(din[0] * dout[1] - din[1] * dout[0], din[0] * dout[0] + din[1] * dout[1]);
}

/* ------------------------------------------------------------------------ */
/*  Nodes, arcs, the heap                                                   */
/* ------------------------------------------------------------------------ */

static unsigned int sk_node(sk_run* r, double const* p, double t)
{
	qaws_skeleton* s = r->s;
	unsigned int i;
	if (!r->record)
		return 0;
	for (i = s->nnodes; i-- > r->node_from;)
	{
		qaws_skeleton_node const* n = &s->nodes[i];
		if (fabs((double)n->time - t) <= r->eps_t && hypot((double)n->position.x - p[0], (double)n->position.y - p[1]) <= r->eps)
			return i;
	}
	if (!sk_grow((void**)&s->nodes, &s->capnodes, s->nnodes + 1, sizeof(qaws_skeleton_node)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return 0;
	}
	s->nodes[s->nnodes].position.x = (qaws_scalar)p[0];
	s->nodes[s->nnodes].position.y = (qaws_scalar)p[1];
	s->nodes[s->nnodes].time = (qaws_scalar)t;
	if (t > s->tmax) s->tmax = t;
	return s->nnodes++;
}

static void sk_arc(sk_run* r, unsigned int a, unsigned int b, unsigned int fa, unsigned int fb, int smooth, sk_geo const* g)
{
	qaws_skeleton* s = r->s;
	if (!r->record || a == b)
		return;
	if (!sk_grow((void**)&s->arcs, &s->caparcs, s->narcs + 1, sizeof(qaws_skeleton_arc)) ||
		!sk_grow((void**)&s->geo, &s->capgeo, s->narcs + 1, sizeof(sk_geo)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	s->arcs[s->narcs].node_a = a;
	s->arcs[s->narcs].node_b = b;
	s->arcs[s->narcs].face_a = fa;
	s->arcs[s->narcs].face_b = fb;
	s->arcs[s->narcs].smooth = smooth;
	s->arcs[s->narcs].boundary = s->nodes[a].time == 0 || s->nodes[b].time == 0;
	s->geo[s->narcs] = *g;
	s->narcs++;
}

static int sk_before(sk_event const* a, sk_event const* b)
{
	return a->t < b->t || (a->t == b->t && a->seq < b->seq);
}

static void sk_push(sk_run* r, unsigned int kind, double t, double const* p, unsigned int corner, unsigned int a, unsigned int b)
{
	unsigned int i;
	sk_event e;
	if (!(t <= 1e300) || t > r->stop + r->eps_t)
		return;
	if (r->collect)
	{
		if (!sk_grow((void**)&r->coll, &r->capcoll, r->ncoll + 1, sizeof(sk_event)))
		{
			r->st = QAWS_STATUS_ALLOCATION_FAILURE;
			return;
		}
		e.t = t; e.p[0] = p[0]; e.p[1] = p[1]; e.kind = kind; e.corner = corner; e.a = a; e.b = b; e.seq = 0;
		r->coll[r->ncoll++] = e;
		return;
	}
	if (!sk_grow((void**)&r->heap, &r->caph, r->nh + 1, sizeof(sk_event)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	e.t = t; e.p[0] = p[0]; e.p[1] = p[1]; e.kind = kind; e.corner = corner; e.a = a; e.b = b; e.seq = r->seq++;
	i = r->nh++;
	while (i > 0 && sk_before(&e, &r->heap[(i - 1) / 2]))
	{
		r->heap[i] = r->heap[(i - 1) / 2];
		i = (i - 1) / 2;
	}
	r->heap[i] = e;
}

static sk_event sk_pop(sk_run* r)
{
	sk_event top = r->heap[0], last = r->heap[--r->nh];
	unsigned int i = 0;
	for (;;)
	{
		unsigned int a = 2 * i + 1, b = a + 1, m = i;
		sk_event const* mv = &last;
		if (a < r->nh && sk_before(&r->heap[a], mv)) { m = a; mv = &r->heap[a]; }
		if (b < r->nh && sk_before(&r->heap[b], mv)) { m = b; mv = &r->heap[b]; }
		if (m == i)
			break;
		r->heap[i] = r->heap[m];
		i = m;
	}
	if (r->nh) r->heap[i] = last;
	return top;
}

/* ------------------------------------------------------------------------ */
/*  Making and ending corners                                               */
/* ------------------------------------------------------------------------ */

/* the earliest time after t0 at which the corner's supports turn tangent */
static double sk_tend(sk_run const* r, sk_corner const* c)
{
	double ts[3], best = 1e300;
	unsigned int n = sk_tangencies(&r->sup[c->in], &r->sup[c->out], ts), i;
	for (i = 0; i < n; i++)
		if (ts[i] > c->t0 + r->eps_t && ts[i] < best)
		{
			sk_sup const* a = &r->sup[c->in];
			sk_sup const* b = &r->sup[c->out];
			double q[4];
			if ((a->kind == SK_CIRCLE && sk_rho(a, ts[i]) < -r->eps) || (b->kind == SK_CIRCLE && sk_rho(b, ts[i]) < -r->eps))
				continue;
			/* where they start to cross (a corner born just before, its event
			   within the tolerances of the touching) the crossing goes on */
			if (sk_cross(a, b, ts[i] + 1e3 * r->eps_t, r->eps, q) == 2 && hypot(q[2] - q[0], q[3] - q[1]) > r->eps)
				continue;
			best = ts[i];
		}
	return best;
}

/* ------------------------------------------------------------------------ */
/*  A grid of where corners and edges pass                                  */
/* ------------------------------------------------------------------------ */

static void sk_pair_add(sk_run* r, sk_ilist* l, unsigned int i, unsigned int reg)
{
	if (!sk_grow((void**)&l->v, &l->cap, l->n + 2, sizeof(unsigned int)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	l->v[l->n++] = i;
	l->v[l->n++] = reg;
}

static unsigned int sk_cell(sk_run const* r, double x, int axis)
{
	double u = floor((x - r->glo[axis]) / r->gh);
	return u < 0 ? 0 : u >= r->gn ? r->gn - 1 : (unsigned int)u;
}

/* entry (i, reg) into the cells of the box lo..hi widened by w; a box over
   many cells (a fast corner, a long edge) into the always list instead */
static void sk_grid_box(sk_run* r, sk_ilist* cells, sk_ilist* always, double const* lo, double const* hi, double w, unsigned int i, unsigned int reg, unsigned int most)
{
	unsigned int x0 = sk_cell(r, lo[0] - w, 0), x1 = sk_cell(r, hi[0] + w, 0);
	unsigned int y0 = sk_cell(r, lo[1] - w, 1), y1 = sk_cell(r, hi[1] + w, 1), x, y;
	if ((x1 - x0 + 1) * (y1 - y0 + 1) > most)
	{
		sk_pair_add(r, always, i, reg);
		return;
	}
	for (y = y0; y <= y1; y++)
		for (x = x0; x <= x1; x++)
			sk_pair_add(r, &cells[y * r->gn + x], i, reg);
}

static void sk_gq_push(sk_run* r, double t, unsigned int i)
{
	unsigned int k;
	if (!sk_grow((void**)&r->gq, &r->capgq, r->ngq + 1, sizeof(sk_gexp)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	k = r->ngq++;
	while (k > 0 && r->gq[(k - 1) / 2].t > t)
	{
		r->gq[k] = r->gq[(k - 1) / 2];
		k = (k - 1) / 2;
	}
	r->gq[k].t = t;
	r->gq[k].i = i;
}

static sk_gexp sk_gq_pop(sk_run* r)
{
	sk_gexp top = r->gq[0], last = r->gq[--r->ngq];
	unsigned int k = 0;
	for (;;)
	{
		unsigned int a = 2 * k + 1, b = a + 1, m = k;
		double mt = last.t;
		if (a < r->ngq && r->gq[a].t < mt) { m = a; mt = r->gq[a].t; }
		if (b < r->ngq && r->gq[b].t < mt) m = b;
		if (m == k)
			break;
		r->gq[k] = r->gq[m];
		k = m;
	}
	if (r->ngq) r->gq[k] = last;
	return top;
}

/* a box round the points q[0..2n) */
static void sk_box(double const* q, unsigned int n, double* lo, double* hi)
{
	unsigned int k;
	lo[0] = hi[0] = q[0]; lo[1] = hi[1] = q[1];
	for (k = 1; k < n; k++)
	{
		if (q[2 * k] < lo[0]) lo[0] = q[2 * k];
		if (q[2 * k] > hi[0]) hi[0] = q[2 * k];
		if (q[2 * k + 1] < lo[1]) lo[1] = q[2 * k + 1];
		if (q[2 * k + 1] > hi[1]) hi[1] = q[2 * k + 1];
	}
}

/* the box grown to hold circle e's arc from a to b (round it when whole) */
static void sk_arc_box(sk_sup const* e, double const* a, double const* b, double rho, int whole, double* lo, double* hi)
{
	double a0 = atan2(a[1] - e->m[1], a[0] - e->m[0]), sw = whole ? 2 * SK_PI : sk_angle_along(e, a, b);
	unsigned int k;
	for (k = 0; k < 4; k++)
	{
		/* the axis directions the arc passes */
		double ang = 0.5 * SK_PI * k, d = (ang - a0) * e->s, x, y;
		while (d < 0) d += 2 * SK_PI;
		while (d >= 2 * SK_PI) d -= 2 * SK_PI;
		if (d > sw)
			continue;
		x = e->m[0] + rho * cos(ang);
		y = e->m[1] + rho * sin(ang);
		if (x < lo[0]) lo[0] = x;
		if (x > hi[0]) hi[0] = x;
		if (y < lo[1]) lo[1] = y;
		if (y > hi[1]) hi[1] = y;
	}
}

#define SK_GS 4                  /* samples along a curved path's window */

/* corner i from time t on: its path while it moves about a cell, then
   entered again (a straight one by its ends, a conic by samples with half a
   step to spare); near its supports' tangency, or turning too fast, always
   looked at */
static void sk_grid_corner(sk_run* r, unsigned int i, double t)
{
	sk_corner* c = &r->c[i];
	double q[2 * (SK_GS + 1)], lo[2], hi[2], ta, v, w = 0, pad = r->grc;
	unsigned int k, tries;
	if (!r->gon || !c->alive)
		return;
	c->greg++;
	ta = c->t0 > t ? c->t0 : t;
	if (c->kind == SK_LINEAR)
	{
		v = hypot(c->v[0], c->v[1]);
		c->gexp = ta + (v * 1e300 > r->gh ? r->gh / v : 1e300);
		sk_pos(r, c, ta, q);
		sk_pos(r, c, c->gexp < 1e299 ? c->gexp : ta, q + 2);
		sk_box(q, 2, lo, hi);
	}
	else
	{
		double vv[2];
		int ok = 0;
		sk_vel(r, c, ta, vv);
		v = hypot(vv[0], vv[1]);
		if (v > 0 && v < 1e300)
		{
			w = r->gh / v;
			if (c->tend < 1e299 && ta + w > ta + 0.5 * (c->tend - ta))
				w = 0.5 * (c->tend - ta);
		}
		for (tries = 0; tries < 6 && w > 1e3 * r->eps_t && !ok; tries++, w *= 0.25)
		{
			double step = 0;
			for (k = 0; k <= SK_GS; k++)
			{
				sk_pos(r, c, ta + w * k / SK_GS, q + 2 * k);
				if (k)
				{
					double d = hypot(q[2 * k] - q[2 * k - 2], q[2 * k + 1] - q[2 * k - 1]);
					if (d > step) step = d;
				}
			}
			ok = step <= 0.5 * r->gh;
			pad = r->grc + step;
		}
		if (!ok)
		{
			c->gexp = 1e300;
			sk_pair_add(r, &r->gac, i, c->greg);
			return;
		}
		w *= 4;                  /* the loop quartered it once more */
		c->gexp = ta + w;
		sk_box(q, SK_GS + 1, lo, hi);
	}
	sk_grid_box(r, r->gc, &r->gac, lo, hi, pad, i, c->greg, 64);
	if (c->gexp < 1e299)
		sk_gq_push(r, c->gexp, i);
}

/* the edge out of corner i from time t on, until a corner's entry runs out:
   a line between straight paths is within their four ends' hull; else the
   edge at samples (an arc with the axis points it passes), with its ends'
   largest step between samples and the circle's growth to spare */
static void sk_grid_edge(sk_run* r, unsigned int i, double t)
{
	sk_corner* c = &r->c[i];
	sk_corner const* n;
	sk_sup const* e;
	double q[8], lo[2], hi[2], ta, tb, pad = r->gre;
	unsigned int k;
	if (!r->gon || !c->alive)
		return;
	c->gedge++;
	n = &r->c[c->next];
	e = &r->sup[c->out];
	if (c->gexp >= 1e299 || n->gexp >= 1e299)
	{
		sk_pair_add(r, &r->gae, i, c->gedge);
		return;
	}
	ta = c->t0 > t ? c->t0 : t;
	if (n->t0 > ta) ta = n->t0;
	tb = c->gexp < n->gexp ? c->gexp : n->gexp;
	if (tb < ta) tb = ta;
	if (e->kind == SK_LINE && c->kind == SK_LINEAR && n->kind == SK_LINEAR)
	{
		sk_pos(r, c, ta, q);
		sk_pos(r, c, tb, q + 2);
		sk_pos(r, n, ta, q + 4);
		sk_pos(r, n, tb, q + 6);
		sk_box(q, 4, lo, hi);
	}
	else
	{
		double last[4], step = 0;
		for (k = 0; k <= SK_GS; k++)
		{
			double tk = ta + (tb - ta) * k / SK_GS, ab[4], blo[2], bhi[2];
			sk_pos(r, c, tk, ab);
			sk_pos(r, n, tk, ab + 2);
			if (e->kind == SK_CIRCLE && c->next != i && hypot(ab[2] - ab[0], ab[3] - ab[1]) <= 1e3 * r->eps)
			{
				/* ends together: a whole circle or nothing, by the rounding */
				sk_pair_add(r, &r->gae, i, c->gedge);
				return;
			}
			if (k)
			{
				double d0 = hypot(ab[0] - last[0], ab[1] - last[1]), d1 = hypot(ab[2] - last[2], ab[3] - last[3]);
				if (d0 > step) step = d0;
				if (d1 > step) step = d1;
			}
			memcpy(last, ab, sizeof(last));
			sk_box(ab, 2, blo, bhi);
			if (e->kind == SK_CIRCLE)
				sk_arc_box(e, ab, ab + 2, fabs(sk_rho(e, tk)), c->next == i, blo, bhi);
			if (!k)
			{
				lo[0] = blo[0]; lo[1] = blo[1]; hi[0] = bhi[0]; hi[1] = bhi[1];
			}
			else
			{
				if (blo[0] < lo[0]) lo[0] = blo[0];
				if (blo[1] < lo[1]) lo[1] = blo[1];
				if (bhi[0] > hi[0]) hi[0] = bhi[0];
				if (bhi[1] > hi[1]) hi[1] = bhi[1];
			}
		}
		pad += step + (tb - ta) / SK_GS;
	}
	sk_grid_box(r, r->ge, &r->gae, lo, hi, pad, i, c->gedge, 512);
}

/* the entries run out by t, again from t (with the edges at the corner) */
static void sk_grid_advance(sk_run* r, double t)
{
	while (r->ngq && r->gq[0].t <= t && r->st == QAWS_STATUS_OK)
	{
		sk_gexp g = sk_gq_pop(r);
		sk_corner const* c = &r->c[g.i];
		if (!c->alive || c->gexp != g.t)
			continue;           /* entered again since */
		sk_grid_corner(r, g.i, t);
		sk_grid_edge(r, g.i, t);
		if (r->c[g.i].prev != g.i)
			sk_grid_edge(r, r->c[g.i].prev, t);
	}
}

/* the alive corners below n0 that can be at p (edges: whose out edge can
   pass through p) at time t, in increasing order, into r->found */
static void sk_grid_query(sk_run* r, int edges, double const* p, double t, unsigned int n0)
{
	unsigned int k, m, pass, cell;
	sk_grid_advance(r, t);
	cell = sk_cell(r, p[1], 1) * r->gn + sk_cell(r, p[0], 0);
	r->nfound = 0;
	for (pass = 0; pass < 2; pass++)
	{
		sk_ilist* s = pass ? (edges ? &r->gae : &r->gac) : (edges ? &r->ge[cell] : &r->gc[cell]);
		m = 0;
		for (k = 0; k + 1 < s->n; k += 2)
		{
			unsigned int i = s->v[k], reg = s->v[k + 1];
			sk_corner const* c = &r->c[i];
			if (!c->alive || (edges ? c->gedge : c->greg) != reg)
				continue;       /* ended, or entered again: the entry goes */
			s->v[m++] = i;
			s->v[m++] = reg;
			if (i >= n0)
				continue;
			if (pass)
			{
				/* on a circle: only where p is on it (an edge through p to the
				   events' tolerance, a corner at p to its position's) */
				sk_sup const* a = &r->sup[c->in];
				sk_sup const* b = &r->sup[c->out];
				if (edges ? b->kind == SK_CIRCLE && fabs(sk_resid(b, p, t)) > r->eps_r
					: (a->kind == SK_CIRCLE && fabs(sk_resid(a, p, t)) > r->rad + 10 * r->eps_r) ||
					  (b->kind == SK_CIRCLE && fabs(sk_resid(b, p, t)) > r->rad + 10 * r->eps_r))
					continue;
			}
			if (!sk_grow((void**)&r->found, &r->capfound, r->nfound + 1, sizeof(unsigned int)))
			{
				r->st = QAWS_STATUS_ALLOCATION_FAILURE;
				return;
			}
			r->found[r->nfound++] = i;
		}
		s->n = m;
	}
	if (r->nfound > 1)
	{
		qsort(r->found, r->nfound, sizeof(unsigned int), sk_cmp_uint);
		for (k = 1, m = 1; k < r->nfound; k++)
			if (r->found[k] != r->found[m - 1])
				r->found[m++] = r->found[k];
		r->nfound = m;
	}
}

static int sk_same_circle(sk_run const* r, sk_sup const* a, sk_sup const* b);

/* a new corner at p, time t, between supports in and out, linked between
   prev and next (SK_NONE: itself); *spike when its supports are coincident
   and opposite */
static unsigned int sk_corner_new(sk_run* r, double const* p, double t, unsigned int in, unsigned int out,
	unsigned int prev, unsigned int next, unsigned int node, int* spike)
{
	sk_corner* c;
	sk_sup const *a, *b;
	unsigned int i;
	*spike = 0;
	if (!sk_grow((void**)&r->c, &r->capc, r->nc + 1, sizeof(sk_corner)) ||
		!sk_grow((void**)&r->tag, &r->captag, r->nc + 1, sizeof(unsigned int)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return SK_NONE;
	}
	i = r->nc++;
	r->tag[i] = 0;
	c = &r->c[i];
	memset(c, 0, sizeof(*c));
	c->p[0] = p[0]; c->p[1] = p[1]; c->t0 = t;
	c->in = in; c->out = out;
	c->prev = prev == SK_NONE ? i : prev;
	c->next = next == SK_NONE ? i : next;
	c->node = node;
	c->alive = 1;
	c->tend = 1e300;
	c->kind = SK_LINEAR;
	a = &r->sup[in];
	b = &r->sup[out];
	if (in == out)
	{
		/* a loop of one support: a marker on it, moving with it */
		c->fake = c->smooth = 1;
		if (a->kind == SK_CIRCLE)
		{
			double x = p[0] - a->m[0], y = p[1] - a->m[1], l = hypot(x, y), rho = sk_rho(a, t);
			if (l > 0)
			{
				c->v[0] = -a->s * x / l; c->v[1] = -a->s * y / l;
				/* on its circle (an event's point is good to its tolerance) */
				c->p[0] = a->m[0] + rho * x / l;
				c->p[1] = a->m[1] + rho * y / l;
			}
		}
		else
		{
			c->v[0] = a->w * a->n[0]; c->v[1] = a->w * a->n[1];
		}
	}
	else if (a->kind == SK_LINE && b->kind == SK_LINE)
	{
		double det = a->n[0] * b->n[1] - a->n[1] * b->n[0], dot = a->n[0] * b->n[0] + a->n[1] * b->n[1];
		if (fabs(det) <= 1e-12)
		{
			if (dot < 0)
				*spike = 1;
			else
			{
				/* collinear: the average speed (exact for equal weights) */
				c->v[0] = 0.5 * (a->w + b->w) * a->n[0];
				c->v[1] = 0.5 * (a->w + b->w) * a->n[1];
				c->smooth = a->w == b->w;
			}
		}
		else
		{
			/* on both lines: an event's point is good to the tolerances, and its
			   error would ride along the whole path */
			double da = a->d + a->w * t, db = b->d + b->w * t;
			double x = (da * b->n[1] - a->n[1] * db) / det, y = (a->n[0] * db - da * b->n[0]) / det;
			c->v[0] = (a->w * b->n[1] - a->n[1] * b->w) / det;
			c->v[1] = (a->n[0] * b->w - a->w * b->n[0]) / det;
			if (hypot(x - p[0], y - p[1]) <= r->rad)
			{
				c->p[0] = x;
				c->p[1] = y;
			}
		}
	}
	else if (a->kind == SK_CIRCLE && b->kind == SK_CIRCLE && a->s == b->s && sk_same_circle(r, a, b))
	{
		/* one circle on both sides: a marker on it */
		double x = p[0] - a->m[0], y = p[1] - a->m[1], l = hypot(x, y);
		c->smooth = 1;
		if (l > 0) { c->v[0] = -a->s * x / l; c->v[1] = -a->s * y / l; }
	}
	else if (sk_opposite(a, b, t, r->eps))
		*spike = 1;
	else if (sk_persistent(a, b, r->eps))
	{
		/* a smooth joint: the tangent point moves along the common normal */
		c->smooth = 1;
		if (a->kind != b->kind)
		{
			sk_sup const* l = a->kind == SK_LINE ? a : b;
			c->v[0] = l->w * l->n[0]; c->v[1] = l->w * l->n[1];
		}
		else
		{
			/* along the line of centres, from the circle whose centre is not
			   here (a circle just opened at p has none yet) */
			sk_sup const* e = hypot(p[0] - a->m[0], p[1] - a->m[1]) >= hypot(p[0] - b->m[0], p[1] - b->m[1]) ? a : b;
			double x = p[0] - e->m[0], y = p[1] - e->m[1], l = hypot(x, y);
			if (l > 0) { c->v[0] = -e->s * x / l; c->v[1] = -e->s * y / l; }
		}
	}
	else
	{
		double pts[4], d0, d1;
		int touching, met;
		sk_cross(a, b, t, r->eps, pts);
		c->kind = SK_CONIC;
		d0 = hypot(pts[0] - p[0], pts[1] - p[1]);
		d1 = hypot(pts[2] - p[0], pts[3] - p[1]);
		c->branch = d1 < d0;
		/* tangent at birth (or just past it: the crossings parting like a
		   square root, both about as near the event's point), or an input's
		   corner between supports nearly tangent there (their other crossing
		   just before it on the edge coming in: the sliver between them
		   closes at once) */
		touching = hypot(pts[2] - pts[0], pts[3] - pts[1]) <= (t == 0 ? r->rad : 1e3 * r->eps);
		met = 0;
#if !QAWS_SCALAR_IS_FLOAT
		/* (float builds keep the first branch along the edge: past their
		   failures the convex branch reaches runs whose faces come up short
		   by the vertices' rounding) */
		if (t > 0)
#else
		if (0)
#endif
		{
			/* two fronts just met (on the double's scale, whatever the
			   tolerance): their crossings one point, or parting since a
			   touching a moment ago */
			double apart = hypot(pts[2] - pts[0], pts[3] - pts[1]), fine = 1e-9 * r->s->ext;
			met = apart <= 1e3 * fine;
			if (!met && apart <= r->rad)
			{
				double ts[3];
				unsigned int nt = sk_tangencies(a, b, ts), k;
				for (k = 0; k < nt; k++)
					if (ts[k] <= t + 1e2 * fine && ts[k] >= t - 1e2 * fine)
						met = 1;
			}
		}
		if (met)
		{
			/* two fronts met: where they cross the region between them opens
			   from nothing, a convex corner; the branch whose angle is the
			   smaller a little later */
			double dt = 1e3 * r->eps, a0, a1;
			c->branch = 0;
			a0 = sk_angle(r, c, t + dt);
			c->branch = 1;
			a1 = sk_angle(r, c, t + dt);
			c->branch = a1 < a0;
		}
		else if (touching)
		{
			/* the branch met first along the edge coming in, a little later */
			double dt = 1e3 * r->eps, q[4], from[2];
			if ((prev != SK_NONE || t == 0) && sk_cross(a, b, t + dt, r->eps, q) == 2)
			{
				double u0, u1;
				if (t == 0)
				{
					from[0] = a->start[0]; from[1] = a->start[1];
				}
				else
					sk_pos(r, &r->c[prev], t + dt, from);
				u0 = sk_along(a, from, q, t + dt);
				u1 = sk_along(a, from, q + 2, t + dt);
				if (a->kind == SK_LINE)
					c->branch = (u1 >= 0 && (u0 < 0 || u1 < u0));
				else
					c->branch = u1 < u0;
			}
		}
		c->tend = sk_tend(r, c);
	}
	if (prev != SK_NONE) r->c[prev].next = i;
	if (next != SK_NONE) r->c[next].prev = i;
	if (!sk_grow((void**)&r->refs, &r->caprefs, r->nsup + 1, sizeof(unsigned int)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return i;
	}
	r->refs[out]++;
	sk_grid_corner(r, i, t);
	sk_grid_edge(r, i, t);
	if (prev != SK_NONE)
		sk_grid_edge(r, prev, t);
	return i;
}

/* an edge on e is nowhere behind its tangent: a line, or a growing circle's
   arc (the region outside it) */
static int sk_ahead(sk_sup const* e)
{
	return e->kind == SK_LINE || e->s < 0;
}

/* the corner ends at node, its path an arc */
static void sk_kill(sk_run* r, unsigned int i, unsigned int node, double t)
{
	sk_corner* c = &r->c[i];
	sk_geo g;
	if (!c->alive)
		return;
	c->alive = 0;
	r->refs[c->out]--;
	memset(&g, 0, sizeof(g));
	g.kind = c->kind == SK_CONIC ? SK_GEO_CONIC : SK_GEO_LINE;
	g.in = c->in; g.out = c->out; g.branch = c->branch; g.t0 = c->t0; g.t1 = t;
	sk_arc(r, c->node, node, c->in, c->out, c->smooth, &g);
}

/* ------------------------------------------------------------------------ */
/*  Events of a corner                                                      */
/* ------------------------------------------------------------------------ */

/* where corner c meets support e after now */
static void sk_hits(sk_run* r, unsigned int ci, unsigned int ei)
{
	sk_corner const* c = &r->c[ci];
	sk_sup const* e = &r->sup[ei];
	double tt[2], p[4], lo = r->now - r->eps_t;
	unsigned int n = 0, j;
	if (c->kind == SK_LINEAR)
	{
		double a[2];
		a[0] = c->p[0] - c->v[0] * c->t0; a[1] = c->p[1] - c->v[1] * c->t0;
		if (e->kind == SK_LINE)
		{
			double den = e->n[0] * c->v[0] + e->n[1] * c->v[1] - e->w;
			if (fabs(den) > 1e-14)
			{
				tt[0] = (e->d - e->n[0] * a[0] - e->n[1] * a[1]) / den;
				n = 1;
			}
		}
		else
		{
			double ax = a[0] - e->m[0], ay = a[1] - e->m[1];
			n = (unsigned int)sk_solve_quad(c->v[0] * c->v[0] + c->v[1] * c->v[1] - 1, 2 * (ax * c->v[0] + ay * c->v[1]) + 2 * e->r0 * e->s,
				ax * ax + ay * ay - e->r0 * e->r0, tt);
		}
		for (j = 0; j < n; j++)
		{
			double q[2];
			if (tt[j] < lo || (e->kind == SK_CIRCLE && sk_rho(e, tt[j]) < -r->eps))
				continue;
			q[0] = a[0] + c->v[0] * tt[j];
			q[1] = a[1] + c->v[1] * tt[j];
			if (c->smooth && tt[j] <= c->t0 + r->eps_t && hypot(q[0] - c->p[0], q[1] - c->p[1]) <= r->eps)
				continue;           /* where it was born */
			sk_push(r, SK_EV_CORNER, tt[j] < r->now ? r->now : tt[j], q, ci, 0, 0);
		}
		return;
	}
	{
		sk_sup const* three[3];
		three[0] = &r->sup[c->in]; three[1] = &r->sup[c->out]; three[2] = e;
		n = sk_meet3(three, p, tt);
		for (j = 0; j < n; j++)
		{
			double q[2];
			if (tt[j] < lo || tt[j] > c->tend + r->eps_t)
				continue;
			if (e->kind == SK_CIRCLE && sk_rho(e, tt[j]) < -r->eps)
				continue;
			{
				/* this corner's crossing, not its twin's: the nearer of the two (near
				   their tangency both are close and the positions ill-conditioned) */
				double two[4], d0, d1;
				sk_cross(&r->sup[c->in], &r->sup[c->out], tt[j], r->eps, two);
				d0 = hypot(two[0] - p[2 * j], two[1] - p[2 * j + 1]);
				d1 = hypot(two[2] - p[2 * j], two[3] - p[2 * j + 1]);
				if ((c->branch ? d1 > d0 : d0 > d1) && (c->branch ? d1 : d0) > r->eps)
					continue;
				q[0] = p[2 * j]; q[1] = p[2 * j + 1];
			}
			if (tt[j] <= c->t0 + r->eps_r && hypot(q[0] - c->p[0], q[1] - c->p[1]) <= r->rad)
				continue;
			sk_push(r, SK_EV_CORNER, tt[j] < r->now ? r->now : tt[j], q, ci, 0, 0);
		}
	}
}

/* when a conic corner's angle reaches sigma (it grows towards 2 pi as its
   supports turn tangent) */
static double sk_splice_time(sk_run const* r, sk_corner const* c, double sigma)
{
	double lo = c->t0, hi = c->tend, flo, fhi;
	unsigned int i, n = 48;
	if (c->kind != SK_CONIC || !(hi < 1e299))
		return 1e300;
	if (sigma >= 2 * SK_PI - 1e-12)
		return sk_angle(r, c, hi - 1e4 * r->eps_t - 1e-6 * (hi - lo)) > SK_PI ? hi : 1e300;
	flo = sk_angle(r, c, lo) - sigma;
	for (i = 1; i <= n; i++)
	{
		double t = c->t0 + (c->tend - c->t0) * i / n;
		fhi = sk_angle(r, c, i == n ? t - 1e4 * r->eps_t - 1e-6 * (c->tend - c->t0) : t) - sigma;
		if (flo < 0 && fhi >= 0)
		{
			unsigned int it;
			hi = t;
			for (it = 0; it < 80 && hi - lo > r->eps_t; it++)
			{
				double m = 0.5 * (lo + hi);
				if (sk_angle(r, c, m) - sigma >= 0) hi = m; else lo = m;
			}
			return hi;
		}
		lo = t;
		flo = fhi;
	}
	return 1e300;
}

/* corner i's queued meeting with the supports ahead: its next kept one (a
   carries 1 + its index, so a stale copy is told apart) */
static void sk_push_cand(sk_run* r, unsigned int i)
{
	sk_corner const* c = &r->c[i];
	if (c->icand < c->ncand)
		sk_push(r, SK_EV_CORNER, c->cand[c->icand].t, c->cand[c->icand].p, i, 1 + c->icand, 0);
}

/* the splitting corner i against every support ahead of it, from time
   floor on: the earliest kept, the first queued */
static void sk_split_cands(sk_run* r, unsigned int i, double floor)
{
	sk_corner* c = &r->c[i];
	double here[2];
	unsigned int e, k, n = 0;
	sk_pos(r, c, r->now, here);
	r->collect = 1;
	r->ncoll = 0;
	for (e = 0; e < r->nsup; e++)
	{
		if (!r->refs[e] || e == c->in || e == c->out)
			continue;
		if (sk_resid(&r->sup[e], here, r->now) < -r->eps)
			continue;
		if (c->kind == SK_LINEAR && r->sup[e].kind == SK_LINE)
		{
			sk_sup const* q = &r->sup[e];
			if (q->n[0] * c->v[0] + q->n[1] * c->v[1] - q->w >= 0)
				continue;
		}
		sk_hits(r, i, e);
	}
	r->collect = 0;
	c = &r->c[i];
	c->more = 0;
	for (k = 0; k < r->ncoll; k++)
	{
		sk_event const* f = &r->coll[k];
		unsigned int j;
		if (!(f->t > floor))
			continue;
		if (n == SK_CANDS && f->t >= c->cand[n - 1].t)
		{
			c->more = 1;
			continue;
		}
		if (n == SK_CANDS)
		{
			c->more = 1;
			n--;
		}
		for (j = n; j > 0 && c->cand[j - 1].t > f->t; j--)
			c->cand[j] = c->cand[j - 1];
		c->cand[j].t = f->t;
		c->cand[j].p[0] = f->p[0];
		c->cand[j].p[1] = f->p[1];
		n++;
	}
	c->ncand = n;
	c->icand = 0;
	sk_push_cand(r, i);
}

/* corner i's queued meeting came and went: the next */
static void sk_next_cand(sk_run* r, unsigned int i, unsigned int a)
{
	sk_corner* c = &r->c[i];
	if (!c->alive || a != 1 + c->icand)
		return;
	c->icand++;
	if (c->icand < c->ncand)
		sk_push_cand(r, i);
	else if (c->more)
		sk_split_cands(r, i, c->cand[c->ncand - 1].t);
}

static void sk_candidates(sk_run* r, unsigned int i)
{
	sk_corner const* c = &r->c[i];
	unsigned int j = c->next, k = c->prev, e;
	double theta;
	if (c->fake)
	{
		/* a marker on a circle: it ends at the centre, or where the skeleton
		   meets its ray (a neighbour's support) */
		sk_sup const* a = &r->sup[c->in];
		if (a->kind == SK_CIRCLE && a->s > 0 && a->r0 >= r->now)
			sk_push(r, SK_EV_CORNER, a->r0, a->m, i, 0, 0);
		if (j != i && r->c[j].out != c->in)
			sk_hits(r, i, r->c[j].out);
		if (k != i && r->c[k].in != c->in)
			sk_hits(r, i, r->c[k].in);
		return;
	}
	if (j != i && r->c[j].out != c->in && r->c[j].out != c->out)
		sk_hits(r, i, r->c[j].out);
	if (k != i && r->c[k].in != c->out && r->c[k].in != c->in)
		sk_hits(r, i, r->c[k].in);
	/* a smooth neighbour meets this one's far support too: from this side
	   the meeting is at the neighbour's tangency, a double root */
	if (j != i && r->c[j].smooth && !r->c[j].fake && r->c[j].kind == SK_LINEAR && (r->sup[r->c[j].out].kind == SK_LINE || r->sup[r->c[j].in].kind == SK_LINE) && c->in != r->c[j].out && c->in != r->c[j].in)
		sk_hits(r, j, c->in);
	if (k != i && k != j && r->c[k].smooth && !r->c[k].fake && r->c[k].kind == SK_LINEAR && (r->sup[r->c[k].in].kind == SK_LINE || r->sup[r->c[k].out].kind == SK_LINE) && c->out != r->c[k].in && c->out != r->c[k].out)
		sk_hits(r, k, c->out);
	/* its circles shrinking to their centres */
	for (e = 0; e < 2; e++)
	{
		sk_sup const* a = &r->sup[e ? c->out : c->in];
		if (a->kind == SK_CIRCLE && a->s > 0 && a->r0 >= r->now && a->r0 <= c->tend + r->eps_t)
			sk_push(r, SK_EV_CORNER, a->r0, a->m, i, 0, 0);
	}
	/* a conic corner born at its supports' tangency: its angle just after */
	theta = c->smooth ? SK_PI : sk_angle(r, c, c->kind == SK_CONIC ? c->t0 + 1e3 * r->eps_t : c->t0);
	if (c->kind == SK_CONIC && c->tend < 1e299)
	{
		double ts = sk_splice_time(r, c, r->s->sigma), q[2];
		if (ts < 1e299)
		{
			sk_pos(r, c, ts, q);
			sk_push(r, SK_EV_SPLICE, ts, q, i, 0, 0);
		}
		else
		{
			/* a convex corner whose supports part: it meets its twin there */
			sk_pos(r, c, c->tend, q);
			sk_push(r, SK_EV_TANGENT, c->tend, q, i, 0, 0);
		}
	}
	/* every support ahead of it: a reflex corner splits edges; between
	   lines a convex one never does first (its edges meet the other front
	   no later), but between arcs bending in it can */
	if (!(theta > SK_PI + 1e-9) && (c->smooth || (sk_ahead(&r->sup[c->in]) && sk_ahead(&r->sup[c->out]))))
		return;
	sk_split_cands(r, i, -1e300);
}

/* support a's touchings with supports from .. upto after floor: the earliest
   kept, the first queued */
static void sk_touch_find(sk_run* r, unsigned int a, double floor)
{
	sk_touches* lst = &r->sq[a];
	unsigned int b, from = lst->from, upto = lst->upto;
	lst->n = lst->i = 0;
	lst->more = 0;
	for (b = from; b < upto; b++)
	{
		sk_sup const* x = &r->sup[a];
		sk_sup const* y = &r->sup[b];
		double ts[3], p[2];
		unsigned int n, i;
		if (b == a || (x->kind == SK_LINE && y->kind == SK_LINE) || !r->refs[b])
			continue;           /* (a support without edges has none again) */
		n = sk_tangencies(x, y, ts);
		for (i = 0; i < n; i++)
		{
			double t = ts[i];
			if (!(t > r->now + r->eps_t))
				continue;
			if ((x->kind == SK_CIRCLE && sk_rho(x, t) < r->eps) || (y->kind == SK_CIRCLE && sk_rho(y, t) < r->eps))
				continue;
			{
				/* the supports come into contact (apart before, crossing after);
				   at a parting tangency their crossings' corners just meet */
				double h = 1e-6 * r->s->ext + 1e3 * r->eps_t, q[4];
				if (sk_cross(x, y, t - h, 0, q) != 0 || sk_cross(x, y, t + h, 0, q) == 0)
					continue;
			}
			if (x->kind != y->kind)
			{
				sk_sup const* l = x->kind == SK_LINE ? x : y;
				sk_sup const* c = x->kind == SK_LINE ? y : x;
				double h = l->n[0] * c->m[0] + l->n[1] * c->m[1] - l->d - l->w * t;
				p[0] = c->m[0] - h * l->n[0];
				p[1] = c->m[1] - h * l->n[1];
			}
			else
			{
				double dx = y->m[0] - x->m[0], dy = y->m[1] - x->m[1], dd = hypot(dx, dy), rx = sk_rho(x, t), ry = sk_rho(y, t), u;
				if (!(dd > r->eps))
					continue;
				/* external: rx along the line of centres; x inside y: -rx */
				u = fabs(dd - (ry - rx)) <= fabs(dd - (rx + ry)) && fabs(dd - (ry - rx)) <= fabs(dd - (rx - ry)) ? -rx : rx;
				p[0] = x->m[0] + u * dx / dd;
				p[1] = x->m[1] + u * dy / dd;
			}
			if (!(t > floor))
				continue;
			{
				unsigned int j;
				if (lst->n == SK_CANDS && t >= lst->v[lst->n - 1].t)
				{
					lst->more = 1;
					continue;
				}
				if (lst->n == SK_CANDS)
				{
					lst->more = 1;
					lst->n--;
				}
				for (j = lst->n; j > 0 && lst->v[j - 1].t > t; j--)
					lst->v[j] = lst->v[j - 1];
				lst->v[j].t = t;
				lst->v[j].p[0] = p[0];
				lst->v[j].p[1] = p[1];
				lst->v[j].b = b;
				lst->n++;
			}
		}
	}
	if (lst->n)
		sk_push(r, SK_EV_SQUEEZE, lst->v[0].t, lst->v[0].p, 1, a, lst->v[0].b);
}

/* support a against supports from on (as they are now) */
static void sk_squeezes(sk_run* r, unsigned int a, unsigned int from)
{
	unsigned int n = r->capsq;
	if (!sk_grow((void**)&r->sq, &r->capsq, r->nsup, sizeof(sk_touches)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return;
	}
	if (r->capsq > n)
		memset(r->sq + n, 0, sizeof(sk_touches) * (r->capsq - n));
	r->sq[a].from = from;
	r->sq[a].upto = r->nsup;
	sk_touch_find(r, a, -1e300);
}

/* support a's queued touching came and went (tag: the event's 1 + index): the next */
static void sk_next_touch(sk_run* r, unsigned int a, unsigned int tag)
{
	sk_touches* l;
	if (a >= r->capsq || !r->refs[a])
		return;
	l = &r->sq[a];
	if (tag != 1 + l->i)
		return;
	l->i++;
	if (l->i < l->n)
		sk_push(r, SK_EV_SQUEEZE, l->v[l->i].t, l->v[l->i].p, 1 + l->i, a, l->v[l->i].b);
	else if (l->more)
		sk_touch_find(r, a, l->v[l->n - 1].t);
}

static void sk_event_at(sk_run* r, double const* p, double t);
static double sk_loop_area(sk_run const* r, unsigned int i, double t);

/* a new support: a circle growing from p at time t */
static unsigned int sk_support_circle(sk_run* r, double const* p, double t, unsigned int node)
{
	sk_sup* g;
	if (!sk_grow((void**)&r->sup, &r->capsup, r->nsup + 1, sizeof(sk_sup)) ||
		!sk_grow((void**)&r->refs, &r->caprefs, r->nsup + 1, sizeof(unsigned int)))
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		return SK_NONE;
	}
	g = &r->sup[r->nsup];
	memset(g, 0, sizeof(*g));
	g->kind = SK_CIRCLE;
	g->m[0] = p[0]; g->m[1] = p[1];
	g->s = -1;
	g->r0 = -t;
	g->w = 1;
	g->na = g->nb = node;
	g->path = g->curve = g->span = SK_NONE;
	r->refs[r->nsup] = 0;
	return r->nsup++;
}

/* a corner with a circle among its supports */
static int sk_curved(sk_run const* r, sk_corner const* c)
{
	return r->sup[c->in].kind == SK_CIRCLE || r->sup[c->out].kind == SK_CIRCLE;
}

/* two supports on one circle (to the rounding of a fit: circles this close
   only meet tangentially, and their crossings mean nothing) */
static int sk_same_circle(sk_run const* r, sk_sup const* a, sk_sup const* b)
{
	double tol = r->eps_r + 1e-6 * (fabs(a->r0) + fabs(b->r0));
	return a->kind == SK_CIRCLE && b->kind == SK_CIRCLE && a->s == b->s &&
		hypot(a->m[0] - b->m[0], a->m[1] - b->m[1]) <= tol && fabs(a->r0 - b->r0) <= tol;
}

/* corner c is at p at time t: by its position, or, with a circle among its
   supports, by both supports passing through p (near a tangency the
   position is ill-conditioned: it moves like a square root of the time,
   and the supports' residuals are second order along the tangent) */
static int sk_at(sk_run const* r, sk_corner const* c, double const* p, double t)
{
	double q[2], d;
	sk_pos(r, c, t, q);
	d = hypot(q[0] - p[0], q[1] - p[1]);
	if (d <= r->eps)
		return 1;
	if (c->in == c->out || sk_same_circle(r, &r->sup[c->in], &r->sup[c->out]))
		return d <= r->rad;  /* on one circle: the residual says nothing, the position is well
		                        conditioned (moving radially); to a fit's rounding */
	if (d > r->far || (r->sup[c->in].kind == SK_LINE && r->sup[c->out].kind == SK_LINE))
		return 0;
	{
		/* the residuals to the tolerance; a little past it when the corner is
		   right there (a point solved at a near double root is good to half
		   the digits) */
		double ra = fabs(sk_resid(&r->sup[c->in], p, t)), rb = fabs(sk_resid(&r->sup[c->out], p, t));
		double tol = d <= 10 * r->eps_r ? 10 * r->eps_r : r->eps_r;
		if (!(ra <= tol && rb <= tol))
			return 0;
	}
	if (c->kind == SK_CONIC)
	{
		/* both crossings of the supports pass through p: the corner is there
		   only when its own is the nearer, unless the supports touch now
		   (their crossings one point, apart only by the time's rounding) or
		   are about to (closing in, the twins about to meet); just past a
		   touching they part and each is where it is */
		double x[4], ts[3];
		unsigned int nt, k;
		if (sk_cross(&r->sup[c->in], &r->sup[c->out], t, r->eps, x) == 2)
		{
			double dn = hypot(x[2 * c->branch] - p[0], x[2 * c->branch + 1] - p[1]);
			double df = hypot(x[2 - 2 * c->branch] - p[0], x[3 - 2 * c->branch] - p[1]);
			if (df < dn)
			{
				nt = sk_tangencies(&r->sup[c->in], &r->sup[c->out], ts);
				for (k = 0; k < nt; k++)
					if ((ts[k] >= t && ts[k] <= t + r->rad) ||
						(ts[k] >= t - r->eps_r && hypot(x[2] - x[0], x[3] - x[1]) <= 1e3 * r->eps))
						return 1;
				return 0;
			}
		}
	}
	return 1;
}

static void sk_settle(sk_run* r, unsigned int i, int spike, double t);

/* corners a and b made neighbours: the events of their meeting */
static void sk_link(sk_run* r, unsigned int a, unsigned int b)
{
	r->c[a].next = b;
	r->c[b].prev = a;
	sk_grid_edge(r, a, r->now);
	if (a == b || !r->c[a].alive || !r->c[b].alive)
		return;
	if (r->c[b].out != r->c[a].in && r->c[b].out != r->c[a].out)
		sk_hits(r, a, r->c[b].out);
	if (r->c[a].in != r->c[b].out && r->c[a].in != r->c[b].in)
		sk_hits(r, b, r->c[a].in);
}

/* the reflex corner opens a circle at its point */
static void sk_splice(sk_run* r, unsigned int i, double t)
{
	sk_corner* c = &r->c[i];
	double p[2];
	unsigned int node, g, prev = c->prev, next = c->next, in = c->in, out = c->out, c1, c2;
	int sp1 = 0, sp2 = 0;
	if (t <= c->t0 + r->eps_t)
	{
		p[0] = c->p[0]; p[1] = c->p[1];     /* where it was born */
	}
	else
		sk_pos(r, c, t, p);
	node = sk_node(r, p, t);
	sk_kill(r, i, node, t);
	g = sk_support_circle(r, p, t, node);
	if (g == SK_NONE)
		return;
	sk_squeezes(r, g, 0);
	{
		/* the corners that split edges meet the opened circle too */
		unsigned int j;
		for (j = 0; j < r->nc; j++)
		{
			sk_corner const* d = &r->c[j];
			double here[2], theta;
			if (!d->alive || d->fake)
				continue;
			theta = d->smooth ? SK_PI : sk_angle(r, d, t);
			if (!(theta > SK_PI + 1e-9) && (d->smooth || (sk_ahead(&r->sup[d->in]) && sk_ahead(&r->sup[d->out]))))
				continue;
			sk_pos(r, d, t, here);
			if (sk_resid(&r->sup[g], here, t) < -r->eps)
				continue;
			sk_hits(r, j, g);
		}
	}
	c1 = sk_corner_new(r, p, t, in, g, prev, next, node, &sp1);
	if (c1 == SK_NONE)
		return;
	c2 = sk_corner_new(r, p, t, g, out, c1, next, node, &sp2);
	if (c2 == SK_NONE)
		return;
	{
		/* the opened circle is tangent to both neighbours by construction:
		   smooth corners at the tangent points, whatever the rounding says */
		unsigned int k;
		for (k = 0; k < 2; k++)
		{
			sk_corner* d = &r->c[k ? c2 : c1];
			sk_sup const* e = &r->sup[k ? out : in];
			if (d->smooth)
				continue;
			d->kind = SK_LINEAR;
			d->smooth = 1;
			d->tend = 1e300;
			if (e->kind == SK_LINE)
			{
				d->v[0] = e->w * e->n[0]; d->v[1] = e->w * e->n[1];
			}
			else
			{
				double x = p[0] - e->m[0], y = p[1] - e->m[1], l = hypot(x, y);
				d->v[0] = l > 0 ? -e->s * x / l : 0;
				d->v[1] = l > 0 ? -e->s * y / l : 0;
			}
			if (k) sp2 = 0; else sp1 = 0;
		}
		/* their paths changed: where they pass again */
		sk_grid_corner(r, c1, t);
		sk_grid_corner(r, c2, t);
		sk_grid_edge(r, c1, t);
		sk_grid_edge(r, c2, t);
		sk_grid_edge(r, r->c[c1].prev, t);
	}
	sk_settle(r, c1, sp1, t);
	sk_settle(r, c2, sp2, t);
}

/* a corner between coincident opposite supports at time t: the edges lie on
   one curve, and the corner runs along it to the nearer end; both ends at
   once when they are as near */
static void sk_spike(sk_run* r, unsigned int i, double t)
{
	unsigned int guard = 0;
	while (r->st == QAWS_STATUS_OK && guard++ < 1000000)
	{
		sk_corner* c = &r->c[i];
		unsigned int a = c->prev, b = c->next, m, node, in = c->in, out = c->out;
		double pa[2], pb[2], here[2], sa, sb, q[2];
		sk_sup const* e = &r->sup[out];
		sk_geo g;
		int sp = 0;
		sk_pos(r, c, t, here);
		if (a == i)
		{
			sk_kill(r, i, sk_node(r, here, t), t);
			return;
		}
		sk_pos(r, &r->c[a], t, pa);
		sk_pos(r, &r->c[b], t, pb);
		sa = sk_along(e, here, pa, t);
		sb = sk_along(e, here, pb, t);
		if (a == b || fabs(sa - sb) <= r->eps)
		{
			q[0] = 0.5 * (pa[0] + pb[0]); q[1] = 0.5 * (pa[1] + pb[1]);
			sb = sa;
		}
		else if (sb < sa) { q[0] = pb[0]; q[1] = pb[1]; }
		else { q[0] = pa[0]; q[1] = pa[1]; }
		node = sk_node(r, q, t);
		/* the run along the support */
		memset(&g, 0, sizeof(g));
		if (e->kind == SK_LINE)
			g.kind = SK_GEO_LINE;
		else
		{
			g.kind = SK_GEO_CIRCLE;
			g.m[0] = e->m[0]; g.m[1] = e->m[1];
			g.rho = fabs(sk_rho(e, t));
			g.a0 = atan2(here[1] - e->m[1], here[0] - e->m[0]);
			g.sweep = e->s * sk_angle_along(e, here, q);
		}
		c->alive = 0;
		r->refs[c->out]--;
		sk_arc(r, c->node, node, in, out, 0, &g);
		if (a == b)
		{
			sk_kill(r, a, node, t);
			return;
		}
		if (fabs(sa - sb) <= r->eps)
		{
			unsigned int aa = r->c[a].prev, bb = r->c[b].next, ain = r->c[a].in, bout = r->c[b].out;
			sk_kill(r, a, node, t);
			sk_kill(r, b, node, t);
			if (aa == b)
				return;                    /* the loop was these three */
			if (ain == bout)
			{
				sk_link(r, aa, bb);
				return;
			}
			m = sk_corner_new(r, q, t, ain, bout, aa, bb, node, &sp);
		}
		else if (sb < sa)
		{
			unsigned int bb = r->c[b].next, bout = r->c[b].out;
			sk_kill(r, b, node, t);
			if (bout == in)
			{
				sk_link(r, a, bb);
				return;
			}
			m = sk_corner_new(r, q, t, in, bout, a, bb, node, &sp);
		}
		else
		{
			unsigned int aa = r->c[a].prev, ain = r->c[a].in;
			sk_kill(r, a, node, t);
			if (ain == out)
			{
				sk_link(r, aa, b);
				return;
			}
			m = sk_corner_new(r, q, t, ain, out, aa, b, node, &sp);
		}
		if (m == SK_NONE)
			return;
		if (!sp)
		{
			sk_settle(r, m, 0, t);
			return;
		}
		i = m;
	}
}

static void sk_settle(sk_run* r, unsigned int i, int spike, double t)
{
	sk_corner* c;
	double q[2];
	if (i == SK_NONE || !r->c[i].alive)
		return;
	c = &r->c[i];
	sk_pos(r, c, t, q);
	if (c->next == i && !c->fake)
	{
		sk_kill(r, i, sk_node(r, q, t), t);
		return;
	}
	if (spike)
	{
		sk_spike(r, i, t);
		return;
	}
	if (!c->fake && c->next == c->prev && r->c[c->next].in == c->out && r->c[c->next].out == c->in)
	{
		/* two corners left on the same two supports */
		double pts[4];
		unsigned int node, o = c->next;
		if (r->sup[c->in].kind == SK_LINE && r->sup[c->out].kind == SK_LINE)
		{
			node = sk_node(r, q, t);
			sk_kill(r, i, node, t);
			sk_kill(r, o, node, t);
			return;
		}
		(void)pts;
	}
	/* other corners already here (a spike run to an event's point): one
	   event with them; corners born with this one are not counted */
	{
		unsigned int j, k;
		sk_grid_query(r, 0, q, t, r->nc);
		for (k = 0; k < r->nfound; k++)
		{
			sk_corner const* o;
			double w[2];
			j = r->found[k];
			o = &r->c[j];
			if (j == i || !o->alive || o->fake || (fabs(o->t0 - c->t0) <= r->eps_t && (o->node == c->node || hypot(o->p[0] - c->p[0], o->p[1] - c->p[1]) <= r->eps)))
				continue;
			/* the other crossing of the same two supports (just past their
			   touching both pass through here) */
			if (c->kind == SK_CONIC && o->kind == SK_CONIC && o->in == c->in && o->out == c->out && o->branch != c->branch)
				continue;
			(void)w;
			if (sk_at(r, o, q, t))
			{
				int m = r->merge;
				r->merge = 1;
				sk_event_at(r, q, t);
				r->merge = m;
				return;
			}
		}
	}
	/* a reflex corner past the splice angle opens at once */
	if (!c->smooth && r->s->sigma < 2 * SK_PI - 1e-12 && sk_angle(r, c, c->kind == SK_CONIC ? t + 1e3 * r->eps_t : t) > r->s->sigma + 1e-9)
	{
		sk_splice(r, i, t);
		return;
	}
	sk_candidates(r, i);
}

/* ------------------------------------------------------------------------ */
/*  The generalised intersection event                                      */
/* ------------------------------------------------------------------------ */

typedef struct sk_chain
{
	unsigned int in, out;       /* supports: the first edge, the last */
	unsigned int cin, cout;     /* the corners before the first edge and after the last */
	double ang;                 /* of the ray back along the first edge */
} sk_chain;

static int sk_cmp_chain(void const* a, void const* b)
{
	double x = ((sk_chain const*)a)->ang, y = ((sk_chain const*)b)->ang;
	return x < y ? -1 : (x > y ? 1 : 0);
}

/* a loop of one circle, from p */
static void sk_bubble(sk_run* r, unsigned int e, double const* p, double t, unsigned int node)
{
	int sp;
	unsigned int f;
	sk_sup const* c = &r->sup[e];
	double q[2], l = hypot(p[0] - c->m[0], p[1] - c->m[1]);
	q[0] = p[0]; q[1] = p[1];
	if (l > 0)
	{
		/* its node on the circle */
		q[0] = c->m[0] + sk_rho(c, t) * (p[0] - c->m[0]) / l;
		q[1] = c->m[1] + sk_rho(c, t) * (p[1] - c->m[1]) / l;
		node = sk_node(r, q, t);
	}
	f = sk_corner_new(r, q, t, e, e, SK_NONE, SK_NONE, node, &sp);
	if (f != SK_NONE)
		sk_candidates(r, f);
}

static void sk_event_at(sk_run* r, double const* p_in, double t)
{
	double p[2];
	unsigned int i, ncons = 0, nsplit = 0, nch = 0, node = 0, stamp = ++r->stamp, k, n0 = r->nc, seen;
	unsigned int *cons = NULL, *split = NULL;
	sk_chain* ch = NULL;
	unsigned int* fresh = NULL;
	int* spikes = NULL;
	unsigned int touch[2];
	int merge = r->merge, skipped = 0;
	touch[0] = r->touch[0]; touch[1] = r->touch[1];
	r->touch[0] = r->touch[1] = SK_NONE;
	p[0] = p_in[0];
	if (++r->depth > 64)
	{
		/* events re-triggering one another at one point */
		r->depth--;
		r->st = QAWS_STATUS_NUMERICAL_FAILURE;
		return;
	}
	p[1] = p_in[1];
	cons = (unsigned int*)malloc(sizeof(unsigned int) * (n0 + 1));
	split = (unsigned int*)malloc(sizeof(unsigned int) * (n0 + 1));
	if (!cons || !split)
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* the corners at p; a circle's marker only when the circle ends here */
collect:
	ncons = nsplit = 0;
	stamp = ++r->stamp;
	sk_grid_query(r, 0, p, t, n0);
	for (k = 0; k < r->nfound; k++)
	{
		sk_corner const* c;
		i = r->found[k];
		c = &r->c[i];
		if (!c->alive)
			continue;
		if (c->fake && c->next == i && !(r->sup[c->in].kind == SK_CIRCLE && fabs(sk_rho(&r->sup[c->in], t)) <= r->eps))
			continue;
		/* corners this point's earlier event made at this time are its result */
		if (!merge && !c->fake && !c->initial && (sk_curved(r, c)
			? c->t0 >= t - r->eps_r && hypot(c->p[0] - p[0], c->p[1] - p[1]) <= r->rad
			: c->t0 >= t - r->eps_t && hypot(c->p[0] - p[0], c->p[1] - p[1]) <= r->eps))
		{
			skipped |= sk_at(r, c, p, t);
			continue;
		}
		if (sk_at(r, c, p, t))
		{
			r->tag[i] = stamp;
			cons[ncons++] = i;
		}
	}
	/* the edges through p (between two corners not at p) */
	sk_grid_query(r, 1, p, t, n0);
	for (k = 0; k < r->nfound; k++)
	{
		sk_corner const* c;
		sk_sup const* e;
		double a[2], b[2], s, l;
		i = r->found[k];
		c = &r->c[i];
		if (!c->alive || r->tag[i] == stamp || r->tag[c->next] == stamp)
			continue;
		e = &r->sup[c->out];
		if (fabs(sk_resid(e, p, t)) > (e->kind == SK_LINE ? r->eps : r->eps_r))
			continue;
		if (c->next == i)
		{
			/* a whole circle */
			split[nsplit++] = i;
			continue;
		}
		sk_pos(r, c, t, a);
		sk_pos(r, &r->c[c->next], t, b);
		s = sk_along(e, a, p, t);
		l = sk_along(e, a, b, t);
		if (s > r->eps && s < l - r->eps)
			split[nsplit++] = i;
	}
	if (!(ncons >= 2 || (ncons == 1 && (nsplit > 0 || r->c[cons[0]].next == cons[0])) || nsplit >= 2))
	{
		/* nothing here but corners born a moment ago nearby: a meeting with
		   them, not this point's earlier event's result */
		if (!merge && skipped && ncons + nsplit > 0)
		{
			merge = 1;
			goto collect;
		}
		goto done;
	}
	if (touch[0] != SK_NONE)
	{
		/* supports touching where one has no edge: nothing meets (the corners
		   of a circle opened a moment ago are within the tolerances of each
		   other, not at an event) */
		unsigned int z, y, has[2] = { 0, 0 };
		for (z = 0; z < ncons + nsplit; z++)
		{
			sk_corner const* c = &r->c[z < ncons ? cons[z] : split[z - ncons]];
			for (y = 0; y < 2; y++)
				if (c->out == touch[y] || (z < ncons && c->in == touch[y]))
					has[y] = 1;
		}
		if (!has[0] || !has[1])
			goto done;
	}
	if (ncons == 0 && nsplit == 2)
	{
		/* two edges touching: a squeeze only when they face each other across
		   the unswept side (their inward normals opposite), not one behind
		   the other */
		double na[2], nb[2];
		unsigned int z;
		for (z = 0; z < 2; z++)
		{
			sk_sup const* e = &r->sup[r->c[split[z]].out];
			double* nn = z ? nb : na;
			if (e->kind == SK_LINE) { nn[0] = e->n[0]; nn[1] = e->n[1]; }
			else
			{
				double x = e->m[0] - p[0], y = e->m[1] - p[1], l = hypot(x, y);
				nn[0] = l > 0 ? e->s * x / l : 0;
				nn[1] = l > 0 ? e->s * y / l : 0;
			}
		}
		if (na[0] * nb[0] + na[1] * nb[1] > -0.5)
			goto done;
	}
	{
		/* among circles the point is good to the clustering's tolerance: the
		   exact meeting of three of the supports there, when near */
		unsigned int sup[3], ns = 0, z, curved = 0;
		for (z = 0; z < ncons + nsplit && ns < 3; z++)
		{
			sk_corner const* c = &r->c[z < ncons ? cons[z] : split[z - ncons]];
			unsigned int e2[2], y, w;
			e2[0] = c->in; e2[1] = c->out;
			for (y = 0; y < 2 && ns < 3; y++)
			{
				if (z >= ncons && y == 0) continue;     /* an edge through p: its own support */
				for (w = 0; w < ns && sup[w] != e2[y]; w++) {}
				if (w == ns) sup[ns++] = e2[y];
			}
		}
		for (z = 0; z < ns; z++) if (r->sup[sup[z]].kind == SK_CIRCLE) curved = 1;
		if (ns == 2 && curved)
		{
			/* two supports: their crossings meet at a tangency, the exact one near */
			double ts[3], q[4];
			unsigned int nt = sk_tangencies(&r->sup[sup[0]], &r->sup[sup[1]], ts);
			for (z = 0; z < nt; z++)
				if (fabs(ts[z] - t) <= r->eps_r)
				{
					sk_cross(&r->sup[sup[0]], &r->sup[sup[1]], ts[z], r->eps, q);
					if (hypot(0.5 * (q[0] + q[2]) - p[0], 0.5 * (q[1] + q[3]) - p[1]) <= r->rad)
					{
						p[0] = 0.5 * (q[0] + q[2]);
						p[1] = 0.5 * (q[1] + q[3]);
						if (ts[z] > t) t = ts[z];
					}
					break;
				}
		}
		if (ns == 3 && curved)
		{
			sk_sup const* three[3];
			double q[4], tq[2];
			unsigned int nq, best = 2;
			double bd = r->rad;
			three[0] = &r->sup[sup[0]]; three[1] = &r->sup[sup[1]]; three[2] = &r->sup[sup[2]];
			nq = sk_meet3(three, q, tq);
			for (z = 0; z < nq; z++)
			{
				double d = hypot(q[2 * z] - p[0], q[2 * z + 1] - p[1]);
				if (d <= bd && fabs(tq[z] - t) <= r->eps_r) { bd = d; best = z; }
			}
			if (best < 2)
			{
				/* kept only when it fits the three supports better (near a double
				   root the solve itself is good to half the digits) */
				double before = 0, after = 0, tt = tq[best] > t ? tq[best] : t, w;
				for (z = 0; z < 3; z++)
				{
					w = fabs(sk_resid(three[z], p, t)); if (w > before) before = w;
					w = fabs(sk_resid(three[z], q + 2 * best, tt)); if (w > after) after = w;
				}
				if (after < before)
				{
					p[0] = q[2 * best]; p[1] = q[2 * best + 1];
					t = tt;
				}
			}
		}
	}
	node = sk_node(r, p, t);
	ch = (sk_chain*)malloc(sizeof(sk_chain) * (ncons + nsplit + 1));
	fresh = (unsigned int*)malloc(sizeof(unsigned int) * (ncons + nsplit + 1));
	spikes = (int*)malloc(sizeof(int) * (ncons + nsplit + 1));
	if (!ch || !fresh || !spikes)
	{
		r->st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* runs of consecutive corners at p; a loop all at p closes, but an arc
	   between two of them that runs round its circle stays as a loop (a
	   bubble) */
	for (k = 0; k < ncons; k++)
	{
		unsigned int a = cons[k], b, w;
		if (r->tag[r->c[a].prev] == stamp && r->c[a].prev != a)
		{
			unsigned int n = 0;
			w = r->c[a].prev;
			while (w != a && r->tag[w] == stamp && n++ < r->nc) w = r->c[w].prev;
			if (w != a)
				continue;
			for (w = r->c[a].next; w != a; w = r->c[w].next)
				if (w < a) break;
			if (w != a)
				continue;
			ch[nch].in = SK_NONE;
			ch[nch].cin = a;
			nch++;
			continue;
		}
		b = a;
		while (r->tag[r->c[b].next] == stamp && r->c[b].next != a)
			b = r->c[b].next;
		ch[nch].in = r->c[a].in;
		ch[nch].out = r->c[b].out;
		ch[nch].cin = r->c[a].prev;
		ch[nch].cout = r->c[b].next;
		ch[nch].ang = (double)a;          /* the run's first corner, for the bubbles */
		nch++;
	}
	/* bubbles: arcs between two corners of a run going round their circle */
	for (k = 0; k < nch; k++)
	{
		unsigned int a = ch[k].in == SK_NONE ? ch[k].cin : (unsigned int)ch[k].ang, w = a, n = 0;
		do
		{
			unsigned int nx = r->c[w].next;
			sk_sup const* e = &r->sup[r->c[w].out];
			if (r->tag[nx] != stamp || (ch[k].in != SK_NONE && w == r->c[ch[k].cout].prev && nx == ch[k].cout))
				break;
			if (e->kind == SK_CIRCLE && fabs(sk_rho(e, t)) > r->eps)
			{
				/* the arc's extent shrinking to 0 (it vanishes) or growing to a whole
				   turn (a bubble): its trend before the event, away from where the
				   two corners are too close to order */
				double pa[2], pb[2], q[2], born = r->c[w].t0 > r->c[nx].t0 ? r->c[w].t0 : r->c[nx].t0, t1, t2, e1, e2;
				t1 = t - 0.2 * (t - born);
				t2 = t - 0.1 * (t - born);
				sk_pos(r, &r->c[w], t1, pa);
				sk_pos(r, &r->c[nx], t1, pb);
				e1 = sk_angle_along(e, pa, pb);
				sk_pos(r, &r->c[w], t2, pa);
				sk_pos(r, &r->c[nx], t2, pb);
				e2 = sk_angle_along(e, pa, pb);
				if (e2 > e1 && e2 > SK_PI)
				{
					q[0] = p[0]; q[1] = p[1];
					sk_bubble(r, r->c[w].out, q, t, node);
				}
			}
			w = nx;
		} while (w != a && n++ < r->nc);
	}
	for (k = 0; k < nsplit; k++)
	{
		unsigned int d = split[k];
		ch[nch].in = ch[nch].out = r->c[d].out;
		ch[nch].cin = d;
		ch[nch].cout = r->c[d].next;
		nch++;
	}
	for (k = 0; k < ncons; k++)
		sk_kill(r, cons[k], node, t);
	/* drop the closed loops, then the chains by angle round p */
	for (i = k = 0; i < nch; i++)
		if (ch[i].in != SK_NONE)
		{
			double dir[2];
			sk_sup const* e;
			ch[k] = ch[i];
			e = &r->sup[ch[k].in];
			sk_dir_at(e, p, dir);
			ch[k].ang = atan2(-dir[1], -dir[0]);
			if (e->kind == SK_CIRCLE)
			{
				/* back along the edge by a short chord, not its tangent: edges
				   tangent at p are ordered by their curvature */
				double rho = fabs(sk_rho(e, t)), a0 = atan2(p[1] - e->m[1], p[0] - e->m[0]), da, q[2];
				da = rho > 0 ? 1e-3 * r->s->ext / rho : 0;
				if (da > 0.2) da = 0.2;
				q[0] = e->m[0] + rho * cos(a0 - e->s * da);
				q[1] = e->m[1] + rho * sin(a0 - e->s * da);
				if (rho > 0)
					ch[k].ang = atan2(q[1] - p[1], q[0] - p[0]);
			}
			k++;
		}
	nch = k;
	qsort(ch, nch, sizeof(sk_chain), sk_cmp_chain);
	/* the last edge of each chain joined to the first of the next */
	for (i = 0; i < nch; i++)
	{
		sk_chain const* a = &ch[i];
		sk_chain const* b = &ch[(i + 1) % nch];
		fresh[i] = SK_NONE;
		spikes[i] = 0;
		if (b->in == a->out)
		{
			/* one support on both sides: no corner */
			sk_link(r, b->cin, a->cout);
			continue;
		}
		fresh[i] = sk_corner_new(r, p, t, b->in, a->out, b->cin, a->cout, node, &spikes[i]);
	}
	for (i = 0; i < nch && r->st == QAWS_STATUS_OK; i++)
		if (fresh[i] != SK_NONE && r->c[fresh[i]].alive)
			sk_settle(r, fresh[i], spikes[i], t);
	/* a loop left enclosing nothing that turns inside out at once (a sliver
	   between two supports just past their tangency, cut off at its far
	   crossing) closes here */
	seen = ++r->stamp;
	for (i = 0; i < 3 * nch && r->st == QAWS_STATUS_OK; i++)
	{
		/* the loops of the new corners and of their neighbours (a new
		   corner spliced at once leaves its loop to them) */
		unsigned int f = i % 3 == 1 ? ch[i / 3].cout : i % 3 == 2 || fresh[i / 3] == SK_NONE ? ch[i / 3].cin : fresh[i / 3], w, n = 0;
		double a0, a1, tiny = 1e-8 * r->s->ext * r->s->ext;
		if (f == SK_NONE || !r->c[f].alive || r->tag[f] == seen)
			continue;
		r->tag[f] = seen;
		for (w = r->c[f].next; w != f && n < 16; w = r->c[w].next)
		{
			r->tag[w] = seen;
			n++;
		}
		if (w != f)
			continue;   /* not a small loop */
		n = 0;
		a0 = sk_loop_area(r, f, t);
		if (fabs(a0) > tiny)
			continue;
		/* inside out already a little later, and more so further on (not a
		   loop closing in between) */
		a1 = sk_loop_area(r, f, t + 1e3 * r->eps_t);
		if (!(a1 < 0 && a1 <= a0))
			continue;
		a1 = sk_loop_area(r, f, t + 1e4 * r->eps_t);
		if (!(a1 < -tiny))
			continue;
		{
			/* its corners where they coincide one point: at one point the
			   loop closes there; at two its faces meet along the sliver, a
			   run along one support at this time */
			unsigned int cs[16], at[16], nn = 0, y, z, ends[2] = { 0, 0 }, sup[2] = { 0, 0 };
			double ps[16][2];
			w = f;
			do
			{
				if (n == 16) break;
				cs[n] = w;
				sk_pos(r, &r->c[w], t, ps[n]);
				for (z = 0; z < n && hypot(ps[z][0] - ps[n][0], ps[z][1] - ps[n][1]) > 1e3 * r->eps; z++) {}
				at[n] = z < n ? at[z] : nn++;
				n++;
				w = r->c[w].next;
			} while (w != f);
			if (w != f || nn > 2)
				continue;   /* a small loop about to close: its own events */
			{
				unsigned int nodes[2];
				for (z = 0; z < nn; z++)
					for (y = 0; y < n; y++)
						if (at[y] == z) { nodes[z] = sk_node(r, ps[y], t); break; }
				/* the edges from one point to the other: their supports face
				   each other across the sliver */
				for (y = 0; y < n && nn == 2; y++)
				{
					unsigned int a = at[y], b = at[(y + 1) % n];
					if (a == b)
						continue;
					ends[a] = y;
					sup[a] = r->c[cs[y]].out;
				}
				for (y = 0; y < n; y++)
					sk_kill(r, cs[y], nodes[at[y]], t);
				if (nn == 2)
				{
					/* along a line where there is one, exact on it */
					unsigned int from = r->sup[sup[0]].kind == SK_LINE || r->sup[sup[1]].kind != SK_LINE ? 0 : 1;
					sk_sup const* e = &r->sup[sup[from]];
					double const* pa = ps[ends[from]];
					double const* pb = ps[ends[1 - from]];
					sk_geo g;
					memset(&g, 0, sizeof(g));
					if (e->kind == SK_LINE)
						g.kind = SK_GEO_LINE;
					else
					{
						g.kind = SK_GEO_CIRCLE;
						g.m[0] = e->m[0]; g.m[1] = e->m[1];
						g.rho = fabs(sk_rho(e, t));
						g.a0 = atan2(pa[1] - e->m[1], pa[0] - e->m[0]);
						g.sweep = e->s * sk_angle_along(e, pa, pb);
					}
					g.t0 = g.t1 = t;
					sk_arc(r, nodes[from], nodes[1 - from], sup[from], sup[1 - from], 0, &g);
				}
			}
		}
	}
	/* markers on circles with new neighbours: a marker outside the arc left
	   between them (a whole circle joined to another loop keeps only part of
	   itself) goes; the others get their new events */
	for (i = 0; i < nch && r->st == QAWS_STATUS_OK; i++)
	{
		unsigned int m[2], j;
		m[0] = ch[i].cin; m[1] = ch[i].cout;
		for (j = 0; j < 2; j++)
		{
			unsigned int f = m[j], a, b;
			sk_corner* c = &r->c[f];
			if (!c->alive || !c->fake || c->next == f)
				continue;
			a = c->prev;
			b = c->next;
			if (r->sup[c->in].kind == SK_CIRCLE && a != b)
			{
				sk_sup const* e = &r->sup[c->in];
				double pa[2], pf[2], pb[2], t2 = t + 1e3 * r->eps_t;
				sk_pos(r, &r->c[a], t2, pa);
				sk_pos(r, c, t2, pf);
				sk_pos(r, &r->c[b], t2, pb);
				if (sk_angle_along(e, pa, pf) + sk_angle_along(e, pf, pb) > sk_angle_along(e, pa, pb) + 1e-9)
				{
					double q[2];
					sk_pos(r, c, t, q);
					/* swept by this event's corners (it is right by its point):
					   its path ends at the event's node, the faces' corner */
					sk_kill(r, f, hypot(q[0] - p[0], q[1] - p[1]) <= r->far ? node : sk_node(r, q, t), t);
					sk_link(r, a, b);
					continue;
				}
			}
			sk_candidates(r, f);
		}
	}
done:
	free(cons); free(split); free(ch); free(fresh); free(spikes);
	r->depth--;
}

/* ------------------------------------------------------------------------ */
/*  The run                                                                 */
/* ------------------------------------------------------------------------ */

static void sk_run_free(sk_run* r)
{
	unsigned int e;
	for (e = 0; r->gc && e < r->gn * r->gn; e++)
	{
		free(r->gc[e].v);
		free(r->ge[e].v);
	}
	free(r->gc); free(r->ge); free(r->gac.v); free(r->gae.v); free(r->gq);
	free(r->sq);
	free(r->found);
	free(r->c); free(r->refs); free(r->heap); free(r->tag); free(r->sup); free(r->coll);
}

/* an event among circles (its time and point good to half the digits near
   a tangency); polygons' events are exact to the rounding */
static int sk_curved_event(sk_run const* r, sk_event const* e)
{
	if (e->kind == SK_EV_SQUEEZE)
		return 1;
	return r->sup[r->c[e->corner].in].kind == SK_CIRCLE || r->sup[r->c[e->corner].out].kind == SK_CIRCLE;
}

/* the signed area of corner i's loop at time t (exact for its arcs) */
static double sk_loop_area(sk_run const* r, unsigned int i, double t)
{
	unsigned int j = i, n = 0;
	double a = 0;
	do
	{
		sk_corner const* c = &r->c[j];
		sk_sup const* e = &r->sup[c->out];
		double p[2], q[2];
		sk_pos(r, c, t, p);
		sk_pos(r, &r->c[c->next], t, q);
		if (e->kind == SK_LINE)
			a += p[0] * q[1] - q[0] * p[1];
		else
		{
			/* to the arc's start, along it (twice its integral of x dy - y dx
			   is rho^2 sweep + rho (mx d sin - my d cos)), on to the next corner */
			double rho = fabs(sk_rho(e, t)), a0 = atan2(p[1] - e->m[1], p[0] - e->m[0]);
			double sw = e->s * (c->next == j ? 2 * SK_PI : sk_angle_along(e, p, q)), a1 = a0 + sw;
			double s0[2], s1[2];
			s0[0] = e->m[0] + rho * cos(a0); s0[1] = e->m[1] + rho * sin(a0);
			s1[0] = e->m[0] + rho * cos(a1); s1[1] = e->m[1] + rho * sin(a1);
			a += p[0] * s0[1] - s0[0] * p[1];
			a += rho * rho * sw + rho * (e->m[0] * (sin(a1) - sin(a0)) - e->m[1] * (cos(a1) - cos(a0)));
			a += s1[0] * q[1] - q[0] * s1[1];
		}
		j = c->next;
	} while (j != i && n++ <= r->nc);
	return 0.5 * a;
}

static qaws_status sk_simulate(qaws_skeleton* s, int record, double stop, sk_run* r)
{
	unsigned int l, i, nheld = 0;
	sk_event held[32];
	memset(r, 0, sizeof(*r));
	r->s = s;
	r->record = record;
	r->stop = stop;
	r->eps = s->tol;
	r->touch[0] = r->touch[1] = SK_NONE;
	r->eps_r = s->tol > 1e-7 * s->ext ? s->tol : 1e-7 * s->ext;
	r->far = 10 * sqrt(r->eps_r * (s->ext > 0 ? s->ext : 1));
	r->rad = sqrt(r->eps_r * (s->ext > 0 ? s->ext : 1));
	r->eps_t = s->tol;
	r->st = QAWS_STATUS_OK;
	r->sup = (sk_sup*)malloc(sizeof(sk_sup) * (s->ninput + 1));
	r->refs = (unsigned int*)calloc(s->ninput + 1, sizeof(unsigned int));
	if (!r->sup || !r->refs)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memcpy(r->sup, s->sup, sizeof(sk_sup) * s->ninput);
	r->nsup = s->ninput;
	r->capsup = s->ninput + 1;
	r->caprefs = s->ninput + 1;
	{
		double wmax = 1;
		for (i = 0; i < s->ninput; i++) if (s->sup[i].w > wmax) wmax = s->sup[i].w;
		r->eps_t = s->tol / wmax;
	}
	{
		double wmin = 1;
		for (i = 0; i < s->ninput; i++) if (s->sup[i].w < wmin) wmin = s->sup[i].w;
		r->tbound = (0.5 * s->ext + 10 * s->tol) / wmin;
	}
	{
		/* the grid over the input, about one corner per cell; a corner is
		   there within the tolerance */
		double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, w;
		unsigned int n;
		for (i = 0; i < s->ninput; i++)
			for (l = 0; l < 2; l++)
			{
				if (s->sup[i].start[l] < lo[l]) lo[l] = s->sup[i].start[l];
				if (s->sup[i].start[l] > hi[l]) hi[l] = s->sup[i].start[l];
			}
		w = hi[0] - lo[0] > hi[1] - lo[1] ? hi[0] - lo[0] : hi[1] - lo[1];
		n = (unsigned int)sqrt((double)s->ninput);
		r->gn = n < 4 ? 4 : n > 512 ? 512 : n;
		r->gh = (w > 0 ? w : 1) / r->gn;
		r->glo[0] = lo[0]; r->glo[1] = lo[1];
		r->grc = r->gre = 4 * r->eps;
		r->gc = (sk_ilist*)calloc((size_t)r->gn * r->gn, sizeof(sk_ilist));
		r->ge = (sk_ilist*)calloc((size_t)r->gn * r->gn, sizeof(sk_ilist));
		if (!r->gc || !r->ge)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		r->gon = 1;
	}
	/* the input's corners, at their nodes */
	for (l = 0; l < s->nloops && r->st == QAWS_STATUS_OK; l++)
	{
		unsigned int f = s->loop_first[l], g = s->loop_first[l + 1], first = r->nc;
		for (i = f; i < g; i++)
		{
			unsigned int in = i == f ? g - 1 : i - 1, ci;
			double p[2];
			int sp;
			p[0] = s->sup[i].start[0];
			p[1] = s->sup[i].start[1];
			ci = sk_corner_new(r, p, 0, in, i, i == f ? SK_NONE : r->nc - 1, SK_NONE, s->sup[i].na, &sp);
			if (ci == SK_NONE)
				break;
			r->c[ci].initial = 1;
		}
		if (r->st == QAWS_STATUS_OK && r->nc > first)
		{
			r->c[r->nc - 1].next = first;
			r->c[first].prev = r->nc - 1;
			sk_grid_edge(r, r->nc - 1, 0);
		}
	}
	for (i = 0; i < r->nsup; i++)
		sk_squeezes(r, i, i + 1);
	{
		unsigned int n = r->nc;
		for (i = 0; i < n && r->st == QAWS_STATUS_OK; i++)
			if (r->c[i].alive)
				sk_settle(r, i, 0, 0);
	}
	while (r->nh && r->st == QAWS_STATUS_OK)
	{
		sk_event e = sk_pop(r);
		double q[2];
		if (e.t > stop)
			break;
		if (e.t > r->tbound)
			break;
		if (r->nc > 100 * (s->ninput + 10))
		{
			/* events re-making one another without end */
			r->st = QAWS_STATUS_NUMERICAL_FAILURE;
			break;
		}
		if (e.kind == SK_EV_SQUEEZE)
		{
			if (!r->refs[e.a] || !r->refs[e.b])
			{
				sk_next_touch(r, e.a, e.corner);
				continue;
			}
		}
		else if (!r->c[e.corner].alive)
			continue;
		if (e.t > r->now)
		{
			r->now = e.t;
			if (record) r->node_from = s->nnodes;
		}
		if (e.kind == SK_EV_SPLICE)
		{
			/* a corner alone at its point splices; with others there it meets them */
			unsigned int j;
			int met = 0;
			for (j = 0; j < r->nc && !met; j++)
				met = j != e.corner && r->c[j].alive && !r->c[j].fake && sk_at(r, &r->c[j], e.p, r->now);
			if (met)
				sk_event_at(r, e.p, r->now);
			else
				sk_splice(r, e.corner, r->now);
			continue;
		}
		{
			/* the queued events within the tolerances of this one are this one:
			   near tangencies their times are near-double roots, good to half
			   the digits, and their points move like square roots of them */
			sk_event keep[32];
			unsigned int nkeep = 0, k;
			double rad = r->rad;
			nheld = 0;
			while (r->nh && r->heap[0].t <= e.t + r->eps_r && nkeep < 32)
			{
				sk_event f = sk_pop(r);
				if (f.kind != SK_EV_SPLICE && sk_curved_event(r, &e) && sk_curved_event(r, &f) && hypot(f.p[0] - e.p[0], f.p[1] - e.p[1]) <= rad)
				{
					/* it comes back after this one, in case this one found
					   nothing of it (a corner gone by then if it did, an edge
					   through its point split there already) */
					if (nheld < 32)
						held[nheld++] = f;
					continue;
				}
				keep[nkeep++] = f;
			}
			for (k = 0; k < nkeep; k++)
				sk_push(r, keep[k].kind, keep[k].t, keep[k].p, keep[k].corner, keep[k].a, keep[k].b);
		}
		if (e.kind == SK_EV_SQUEEZE)
		{
			r->touch[0] = e.a;
			r->touch[1] = e.b;
			sk_event_at(r, e.p, r->now);
			sk_next_touch(r, e.a, e.corner);
		}
		else
		{
			if (r->c[e.corner].fake || e.t >= r->now)
			{
				/* the solved meeting point (a conic corner's own position is
				   ill-conditioned there) */
				q[0] = e.p[0]; q[1] = e.p[1];
			}
			else
				sk_pos(r, &r->c[e.corner], r->now, q);
			sk_event_at(r, q, r->now);
			if (e.kind == SK_EV_CORNER && e.a)
				sk_next_cand(r, e.corner, e.a);
		}
		for (i = 0; i < nheld; i++)
			if (held[i].kind == SK_EV_SQUEEZE ? r->refs[held[i].a] && r->refs[held[i].b] : r->c[held[i].corner].alive)
				sk_push(r, held[i].kind, held[i].t, held[i].p, held[i].corner, held[i].a, held[i].b);
			else if (held[i].kind == SK_EV_SQUEEZE)
				sk_next_touch(r, held[i].a, held[i].corner);
	}
	if (r->st == QAWS_STATUS_OK && stop >= r->tbound)
	{
		/* loops left enclosing nothing (slivers closed at tangencies) go;
		   any other leftover is a failure */
		unsigned int* seen = (unsigned int*)calloc(r->nc + 1, sizeof(unsigned int));
		int positive = 0;
		if (!seen)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < r->nc && !positive; i++)
			if (r->c[i].alive && !seen[i])
			{
				unsigned int j = i, n = 0;
				do { seen[j] = 1; j = r->c[j].next; } while (j != i && n++ <= r->nc);
				positive = fabs(sk_loop_area(r, i, r->now)) > 10 * r->eps * r->s->ext;
			}
		free(seen);
		if (!positive)
			for (i = 0; i < r->nc; i++)
				if (r->c[i].alive)
				{
					r->c[i].alive = 0;
					r->refs[r->c[i].out]--;
				}
	}
	return r->st;
}

/* ------------------------------------------------------------------------ */
/*  The input                                                               */
/* ------------------------------------------------------------------------ */

/* the knot span of a B-spline holding t (0 for a Bezier) */
static unsigned int sk_span(qaws_curve const* c, double t)
{
	qaws_scalar buf[256], *k = buf;
	unsigned int n = 0, i, span = 0, deg;
	if (c->kind != QAWS_CURVE_KIND_BSPLINE && c->kind != QAWS_CURVE_KIND_NURBS)
		return 0;
	deg = c->degree;
	if (qaws_bspline_get_knots(c, buf, 256, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		k = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
		if (!k) return 0;
		qaws_bspline_get_knots(c, k, n, &n);
	}
	for (i = deg + 1; i + 1 + deg < n; i++)
		if ((double)k[i] > (double)k[i - 1] && (double)k[i] <= t)
			span++;
	if (k != buf) free(k);
	return span;
}

static int sk_is_line(qaws_curve const* c)
{
	return c->degree == 1 && (c->kind == QAWS_CURVE_KIND_BSPLINE || c->kind == QAWS_CURVE_KIND_BEZIER);
}

void qaws_skeleton_destroy(qaws_skeleton* s)
{
	unsigned int i;
	if (!s)
		return;
	if (s->edges)
		for (i = 0; i < s->ninput; i++)
			if (s->edges[i]) qaws_curve_destroy(s->edges[i]);
	free(s->edges);
	free(s->sup);
	free(s->loop_first);
	free(s->nodes);
	free(s->arcs);
	free(s->geo);
	free(s->face_first);
	free(s->face_nodes);
	free(s);
}

/* the input curves of output curve i of a clip path: the source both its
   vertices lie on */
static void sk_source(qaws_clip_vertex const* v, unsigned int n, unsigned int i, unsigned int* path, unsigned int* curve)
{
	qaws_clip_vertex const* a = &v[i];
	qaws_clip_vertex const* b = &v[(i + 1) % n];
	*path = a->path_a;
	*curve = a->curve_a;
	if ((a->path_a == b->path_a && a->curve_a == b->curve_a) || (a->path_a == b->path_b && a->curve_a == b->curve_b))
		return;
	if (a->curve_b != QAWS_CLIP_NONE_INDEX &&
		((a->path_b == b->path_a && a->curve_b == b->curve_a) || (a->path_b == b->path_b && a->curve_b == b->curve_b)))
	{
		*path = a->path_b;
		*curve = a->curve_b;
	}
}

/* the faces: each support's arcs walked from its edge's end back to its start */
static void sk_arc_eval(qaws_skeleton const* s, unsigned int i, double u, double* p, double* t);

static qaws_status sk_faces(qaws_skeleton* s)
{
	unsigned int f, i, total = 0, *cnt, *first, *list, *used, *back, *walks;
	int open = 0;
	cnt = (unsigned int*)calloc(s->nsup + 1, sizeof(unsigned int));
	first = (unsigned int*)calloc(s->nsup + 1, sizeof(unsigned int));
	list = (unsigned int*)malloc(sizeof(unsigned int) * (2 * s->narcs + 1));
	used = (unsigned int*)calloc(s->narcs + 1, sizeof(unsigned int));
	back = (unsigned int*)calloc(s->narcs + 1, sizeof(unsigned int));
	walks = (unsigned int*)calloc(2 * s->narcs + 2, sizeof(unsigned int));
	s->face_first = (unsigned int*)malloc(sizeof(unsigned int) * (s->nsup + 1));
	s->face_nodes = (unsigned int*)malloc(sizeof(unsigned int) * (4 * s->narcs + 2 * s->nsup + 1));
	if (!cnt || !first || !list || !used || !back || !walks || !s->face_first || !s->face_nodes)
	{
		free(cnt); free(first); free(list); free(used); free(back); free(walks);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < s->narcs; i++)
	{
		cnt[s->arcs[i].face_a]++;
		if (s->arcs[i].face_b != s->arcs[i].face_a) cnt[s->arcs[i].face_b]++;
	}
	for (f = 0; f < s->nsup; f++) { first[f + 1] = first[f] + cnt[f]; cnt[f] = 0; }
	for (i = 0; i < s->narcs; i++)
	{
		unsigned int a = s->arcs[i].face_a, b = s->arcs[i].face_b;
		list[first[a] + cnt[a]++] = i;
		if (b != a) list[first[b] + cnt[b]++] = i;
	}
	for (f = 0; f < s->nsup; f++)
	{
		unsigned int cur = s->sup[f].nb, guard = 0, na = s->sup[f].na;
		double din[2] = { 0, 0 };
		unsigned int last = SK_NONE;
		int have_dir = 0;
		s->face_first[f] = total;
		s->face_nodes[total++] = na;
		if (cur != na)
			s->face_nodes[total++] = cur;
		if (s->sup[f].input)
		{
			/* arriving along the edge */
			double q[2];
			q[0] = (double)s->nodes[cur].position.x; q[1] = (double)s->nodes[cur].position.y;
			sk_dir_at(&s->sup[f], q, din);
			have_dir = 1;
		}
		while ((cur != na || guard == 0) && guard++ <= 2 * (first[f + 1] - first[f]) + 1)
		{
			unsigned int j, best = SK_NONE;
			double bestd = 1e300, q[2];
			q[0] = (double)s->nodes[cur].position.x; q[1] = (double)s->nodes[cur].position.y;
			for (j = first[f]; j < first[f + 1]; j++)
			{
				qaws_skeleton_arc const* a = &s->arcs[list[j]];
				double x[2], tt, ang, dd;
				int fwd;
				if (used[list[j]] == f + 1 || (a->node_a != cur && a->node_b != cur))
					continue;
				fwd = a->node_a == cur;
				if (!have_dir)
				{
					/* a splice's face starts out along the arc of the corner
					   leaving its circle */
					if (fwd && a->face_a == f) { best = list[j]; break; }
					if (best == SK_NONE) best = list[j];
					continue;
				}
				/* the arc first clockwise from the way back keeps the face on the left */
				sk_arc_eval(s, list[j], fwd ? 1e-4 : 1 - 1e-4, x, &tt);
				ang = atan2(-din[1], -din[0]) - atan2(x[1] - q[1], x[0] - q[0]);
				while (ang <= 1e-12) ang += 2 * SK_PI;
				while (ang > 2 * SK_PI) ang -= 2 * SK_PI;
				dd = ang;
				if (dd < bestd) { bestd = dd; best = list[j]; }
			}
			if (best == SK_NONE && have_dir)
			{
				/* a dead end: back along a path into the face (a slit), first the
				   one entered by earlier, by the same turn; the arc just taken last */
				for (j = first[f]; j < first[f + 1]; j++)
				{
					qaws_skeleton_arc const* a = &s->arcs[list[j]];
					double x[2], tt, ang;
					int fwd;
					if (used[list[j]] != f + 1 || back[list[j]] == f + 1 || list[j] == last || (a->node_a != cur && a->node_b != cur))
						continue;
					fwd = a->node_a == cur;
					sk_arc_eval(s, list[j], fwd ? 1e-4 : 1 - 1e-4, x, &tt);
					ang = atan2(-din[1], -din[0]) - atan2(x[1] - q[1], x[0] - q[0]);
					while (ang <= 1e-12) ang += 2 * SK_PI;
					while (ang > 2 * SK_PI) ang -= 2 * SK_PI;
					if (ang < bestd) { bestd = ang; best = list[j]; }
				}
				if (best == SK_NONE && last != SK_NONE && back[last] != f + 1)
					best = last;
				if (best != SK_NONE)
					back[best] = f + 1;
			}
			if (best == SK_NONE)
				break;
			used[best] = f + 1;
			if (s->arcs[best].face_a == f) walks[2 * best]++;
			if (s->arcs[best].face_b == f) walks[2 * best + 1]++;
			last = best;
			{
				int fwd = s->arcs[best].node_a == cur;
				double x[2], tt, y[2], l;
				cur = fwd ? s->arcs[best].node_b : s->arcs[best].node_a;
				sk_arc_eval(s, best, fwd ? 1 - 1e-4 : 1e-4, x, &tt);
				y[0] = (double)s->nodes[cur].position.x; y[1] = (double)s->nodes[cur].position.y;
				l = hypot(y[0] - x[0], y[1] - x[1]);
				if (l > 0) { din[0] = (y[0] - x[0]) / l; din[1] = (y[1] - x[1]) / l; have_dir = 1; }
			}
			if (cur != na)
				s->face_nodes[total++] = cur;
		}
		if (cur != na)
			open = 1;
	}
	s->face_first[s->nsup] = total;
	/* each arc between two faces walked by both */
	for (i = 0; i < s->narcs && !open; i++)
		if (s->arcs[i].face_a != s->arcs[i].face_b && !(walks[2 * i] && walks[2 * i + 1]))
			open = 1;
	free(cnt); free(first); free(list); free(used); free(back); free(walks);
	/* a face whose walk never came back, or an arc left out: arcs that do
	   not meet as a skeleton */
	return open ? QAWS_STATUS_NUMERICAL_FAILURE : QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Free-form curves: G1 arc splines (biarcs)                               */
/* ------------------------------------------------------------------------ */

typedef struct sk_prim
{
	int arc;                    /* else a line from a to b */
	double a[2], b[2];
	double m[2], r, a0, sweep;
} sk_prim;

/* the arc (or line) leaving a along the unit tangent ta and ending at b */
static void sk_arc_from(double const* a, double const* ta, double const* b, sk_prim* o)
{
	double c[2], cr, cc, dt, k;
	c[0] = b[0] - a[0]; c[1] = b[1] - a[1];
	cr = ta[0] * c[1] - ta[1] * c[0];
	dt = ta[0] * c[0] + ta[1] * c[1];
	cc = c[0] * c[0] + c[1] * c[1];
	o->a[0] = a[0]; o->a[1] = a[1]; o->b[0] = b[0]; o->b[1] = b[1];
	if (!(cc > 0) || fabs(cr) <= 1e-12 * cc)
	{
		o->arc = 0;
		return;
	}
	k = 2 * cr / cc;
	o->arc = 1;
	o->m[0] = a[0] - ta[1] / k;
	o->m[1] = a[1] + ta[0] / k;
	o->r = 1 / fabs(k);
	o->a0 = atan2(a[1] - o->m[1], a[0] - o->m[0]);
	o->sweep = 2 * atan2(cr, dt);
}

/* the distance from p to a line or arc */
static double sk_prim_dist(sk_prim const* g, double const* p)
{
	if (!g->arc)
	{
		double d[2], l2, s;
		d[0] = g->b[0] - g->a[0]; d[1] = g->b[1] - g->a[1];
		l2 = d[0] * d[0] + d[1] * d[1];
		s = l2 > 0 ? ((p[0] - g->a[0]) * d[0] + (p[1] - g->a[1]) * d[1]) / l2 : 0;
		if (s < 0) s = 0;
		if (s > 1) s = 1;
		return hypot(g->a[0] + s * d[0] - p[0], g->a[1] + s * d[1] - p[1]);
	}
	{
		double ang = atan2(p[1] - g->m[1], p[0] - g->m[0]) - g->a0, two = 2 * SK_PI;
		/* the angle from the arc's start in its sweep's direction */
		if (g->sweep < 0) ang = -ang;
		while (ang < 0) ang += two;
		while (ang >= two) ang -= two;
		if (ang <= fabs(g->sweep))
			return fabs(hypot(p[0] - g->m[0], p[1] - g->m[1]) - g->r);
		return fmin(hypot(p[0] - g->a[0], p[1] - g->a[1]), hypot(p[0] - g->b[0], p[1] - g->b[1]));
	}
}

/* the point and unit tangent of the curve at t (the chord's direction where
   the derivative vanishes) */
static int sk_point_tangent(qaws_curve const* c, double t, double tother, double* p, double* tg)
{
	qaws_eval_result_2d e, f;
	double l;
	if (qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &e) != QAWS_STATUS_OK)
		return 0;
	p[0] = (double)e.position.x; p[1] = (double)e.position.y;
	tg[0] = (double)e.d1.x; tg[1] = (double)e.d1.y;
	l = hypot(tg[0], tg[1]);
	if (!(l > 1e-300))
	{
		if (qaws_curve_evaluate_2d(c, (qaws_scalar)(t + 1e-6 * (tother - t)), QAWS_EVAL_FLAG_POSITION, &f) != QAWS_STATUS_OK)
			return 0;
		tg[0] = (double)f.position.x - p[0]; tg[1] = (double)f.position.y - p[1];
		if (tother < t) { tg[0] = -tg[0]; tg[1] = -tg[1]; }
		l = hypot(tg[0], tg[1]);
		if (!(l > 0))
			return 0;
	}
	tg[0] /= l; tg[1] /= l;
	return 1;
}

/* the equal-tangent biarc from (p0, t0) to (p1, t1); 0 when it does not exist */
static int sk_biarc(double const* p0, double const* t0, double const* p1, double const* t1, sk_prim* o)
{
	double v[2], vt, tt, d, q1[2], q2[2], j[2], jt[2], l;
	v[0] = p1[0] - p0[0]; v[1] = p1[1] - p0[1];
	vt = v[0] * (t0[0] + t1[0]) + v[1] * (t0[1] + t1[1]);
	tt = 1 - (t0[0] * t1[0] + t0[1] * t1[1]);
	if (tt < 1e-12)
	{
		double den = 4 * (v[0] * t1[0] + v[1] * t1[1]);
		if (!(den > 0))
			return 0;
		d = (v[0] * v[0] + v[1] * v[1]) / den;
	}
	else
		d = (-vt + sqrt(vt * vt + 2 * tt * (v[0] * v[0] + v[1] * v[1]))) / (2 * tt);
	if (!(d > 0))
		return 0;
	q1[0] = p0[0] + d * t0[0]; q1[1] = p0[1] + d * t0[1];
	q2[0] = p1[0] - d * t1[0]; q2[1] = p1[1] - d * t1[1];
	j[0] = 0.5 * (q1[0] + q2[0]); j[1] = 0.5 * (q1[1] + q2[1]);
	jt[0] = q2[0] - q1[0]; jt[1] = q2[1] - q1[1];
	l = hypot(jt[0], jt[1]);
	if (!(l > 0))
	{
		jt[0] = t0[0]; jt[1] = t0[1];
	}
	else
	{
		jt[0] /= l; jt[1] /= l;
	}
	sk_arc_from(p0, t0, j, &o[0]);
	sk_arc_from(j, jt, p1, &o[1]);
	return 1;
}

typedef struct sk_prims
{
	sk_prim* g;
	unsigned int n, cap;
} sk_prims;

/* the piece [ta, tb] of c as biarcs within tol, halved until it fits */
static qaws_status sk_fit(qaws_curve const* c, double ta, double tb, double tol, unsigned int depth, sk_prims* out)
{
	double p0[2], p1[2], t0[2], t1[2], worst = 0;
	sk_prim g[2];
	unsigned int k, i;
	int ok = sk_point_tangent(c, ta, tb, p0, t0) && sk_point_tangent(c, tb, ta, p1, t1) && sk_biarc(p0, t0, p1, t1, g);
	if (ok)
	{
		for (k = 1; k < 16 && worst <= tol; k++)
		{
			qaws_eval_result_2d e;
			double q[2], d0, d1;
			if (qaws_curve_evaluate_2d(c, (qaws_scalar)(ta + (tb - ta) * k / 16), QAWS_EVAL_FLAG_POSITION, &e) != QAWS_STATUS_OK)
				return QAWS_STATUS_NUMERICAL_FAILURE;
			q[0] = (double)e.position.x; q[1] = (double)e.position.y;
			d0 = sk_prim_dist(&g[0], q);
			d1 = sk_prim_dist(&g[1], q);
			if ((d0 < d1 ? d0 : d1) > worst) worst = d0 < d1 ? d0 : d1;
		}
	}
	if ((!ok || worst > tol) && depth < 24 && tb - ta > 1e-12 * (fabs(ta) + fabs(tb) + 1))
	{
		qaws_status st = sk_fit(c, ta, 0.5 * (ta + tb), tol, depth + 1, out);
		if (st == QAWS_STATUS_OK)
			st = sk_fit(c, 0.5 * (ta + tb), tb, tol, depth + 1, out);
		return st;
	}
	if (!ok)
	{
		/* a chord where no biarc exists (a vanishing piece) */
		g[0].arc = 0;
		g[0].a[0] = p0[0]; g[0].a[1] = p0[1]; g[0].b[0] = p1[0]; g[0].b[1] = p1[1];
		g[1] = g[0];
		g[1].a[0] = p1[0]; g[1].a[1] = p1[1];
	}
	for (i = 0; i < 2; i++)
	{
		if (hypot(g[i].b[0] - g[i].a[0], g[i].b[1] - g[i].a[1]) <= 1e-12 * tol + 1e-300)
			continue;
		if (out->n && g[i].arc && out->g[out->n - 1].arc)
		{
			/* on the circle of the arc before (a biarc of a circular piece): one arc */
			sk_prim* h = &out->g[out->n - 1];
			if (hypot(h->m[0] - g[i].m[0], h->m[1] - g[i].m[1]) <= 1e-9 * (h->r + g[i].r) && fabs(h->r - g[i].r) <= 1e-9 * (h->r + g[i].r) &&
				(h->sweep > 0) == (g[i].sweep > 0))
			{
				h->b[0] = g[i].b[0]; h->b[1] = g[i].b[1];
				h->sweep += g[i].sweep;
				continue;
			}
		}
		if (!sk_grow((void**)&out->g, &out->cap, out->n + 1, sizeof(sk_prim)))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		out->g[out->n++] = g[i];
	}
	return QAWS_STATUS_OK;
}

/* a whole curve as biarcs: span by span (a span's end can be a corner) */
static qaws_status sk_fit_curve(qaws_curve const* c, double tol, sk_prims* out)
{
	qaws_range rg = qaws_curve_get_parameter_range(c);
	unsigned int k, n = c->span_count ? c->span_count : 1;
	qaws_status st = QAWS_STATUS_OK;
	for (k = 0; k < n && st == QAWS_STATUS_OK; k++)
	{
		double a = c->span_count ? (double)c->span_boundaries[k] : (double)rg.min_value;
		double b = c->span_count ? (double)c->span_boundaries[k + 1] : (double)rg.max_value;
		if (b > a)
			st = sk_fit(c, a, b, tol, 0, out);
	}
	return st;
}

static qaws_status sk_add_support(qaws_skeleton* s, unsigned int* cap, sk_sup const* e, double const* start)
{
	if (!sk_grow((void**)&s->sup, cap, s->nsup + 1, sizeof(sk_sup)) ||
		!sk_grow((void**)&s->nodes, &s->capnodes, s->nnodes + 1, sizeof(qaws_skeleton_node)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s->sup[s->nsup] = *e;
	s->sup[s->nsup].start[0] = start[0];
	s->sup[s->nsup].start[1] = start[1];
	s->nodes[s->nnodes].position.x = (qaws_scalar)start[0];
	s->nodes[s->nnodes].position.y = (qaws_scalar)start[1];
	s->nodes[s->nnodes].time = 0;
	s->sup[s->nsup].na = s->nnodes++;
	s->nsup++;
	return QAWS_STATUS_OK;
}

static qaws_status sk_edge_curve(sk_sup const* e, qaws_skeleton_node const* a, qaws_skeleton_node const* b, qaws_curve** out)
{
	if (e->kind == SK_LINE)
	{
		qaws_scalar cp[4];
		qaws_bezier_desc d;
		cp[0] = a->position.x; cp[1] = a->position.y; cp[2] = b->position.x; cp[3] = b->position.y;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
		return qaws_curve_create_bezier(&d, out);
	}
	{
		qaws_arc_segment g;
		qaws_arc_desc d;
		memset(&g, 0, sizeof(g));
		g.center[0] = (qaws_scalar)e->m[0]; g.center[1] = (qaws_scalar)e->m[1];
		g.radius = (qaws_scalar)e->r0;
		g.angle_start = (qaws_scalar)e->a0;
		g.angle_end = (qaws_scalar)(e->a0 + e->sweep);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
		return qaws_curve_create_arc(&d, out);
	}
}

qaws_status qaws_skeleton_create(qaws_skeleton_desc const* desc, qaws_skeleton** out)
{
	qaws_skeleton* s;
	qaws_clip_result* region = NULL;
	qaws_status st;
	unsigned int np, pi, capsup = 0;
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, fit_tol;
	int fitted = 0;
	sk_run r;
	if (!desc || !out || (desc->path_count && !desc->paths))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->splice_angle != 0 && !((double)desc->splice_angle >= SK_PI - 1e-12 && (double)desc->splice_angle <= 2 * SK_PI + 1e-12))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = qaws_clip_boolean(QAWS_CLIP_UNION, desc->fill_rule, desc->paths, desc->path_count, NULL, 0, &region);
	if (st != QAWS_STATUS_OK)
		return st;
	s = (qaws_skeleton*)calloc(1, sizeof(qaws_skeleton));
	if (!s)
	{
		qaws_clip_result_destroy(region);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	s->sigma = desc->splice_angle != 0 ? (double)desc->splice_angle : 2 * SK_PI;
	np = qaws_clip_result_get_path_count(region);
	{
		/* the free-form curves' fit tolerance, of the region's extent */
		double ext0 = 0;
		for (pi = 0; pi < np; pi++)
		{
			qaws_path_2d path;
			qaws_vec2 a, b;
			qaws_clip_result_get_path(region, pi, &path);
			if (qaws_path_compute_bounds_2d(&path, &a, &b) == QAWS_STATUS_OK && hypot((double)b.x - a.x, (double)b.y - a.y) > ext0)
				ext0 = hypot((double)b.x - a.x, (double)b.y - a.y);
		}
		fit_tol = desc->curve_tolerance > 0 ? (double)desc->curve_tolerance : 1e-4 * (ext0 > 0 ? ext0 : 1);
	}
	s->loop_first = (unsigned int*)malloc(sizeof(unsigned int) * (np + 1));
	if (!s->loop_first)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* the supports, loop by loop; a node at every input corner */
	for (pi = 0; pi < np && st == QAWS_STATUS_OK; pi++)
	{
		qaws_path_2d path;
		qaws_clip_vertex const* vx = NULL;
		unsigned int nv = 0, ci, first = s->nsup, first_node = s->nnodes, k;
		qaws_clip_result_get_path(region, pi, &path);
		qaws_clip_result_get_vertices(region, 0, pi, &vx, &nv);
		s->loop_first[s->nloops] = first;
		for (ci = 0; ci < path.curve_count && st == QAWS_STATUS_OK; ci++)
		{
			qaws_curve const* c = path.curves[ci];
			unsigned int sp = 0, sc = 0;
			if (vx && nv == path.curve_count)
				sk_source(vx, nv, ci, &sp, &sc);
			if (sk_is_line(c))
			{
				qaws_scalar buf[2 * 64], *cp = buf;
				unsigned int n = 0, j;
				if (qaws_curve_get_control_points(c, buf, 64, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
				{
					cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
					if (!cp) { st = QAWS_STATUS_ALLOCATION_FAILURE; break; }
					qaws_curve_get_control_points(c, cp, n, &n);
				}
				for (j = 0; j + 1 < n && st == QAWS_STATUS_OK; j++)
				{
					double a[2], b[2], len;
					sk_sup e;
					a[0] = (double)cp[2 * j]; a[1] = (double)cp[2 * j + 1];
					b[0] = (double)cp[2 * j + 2]; b[1] = (double)cp[2 * j + 3];
					len = hypot(b[0] - a[0], b[1] - a[1]);
					if (!(len > 0))
						continue;
					memset(&e, 0, sizeof(e));
					e.kind = SK_LINE;
					e.n[0] = -(b[1] - a[1]) / len;
					e.n[1] = (b[0] - a[0]) / len;
					e.d = e.n[0] * a[0] + e.n[1] * a[1];
					e.w = 1;
					e.input = 1;
					e.path = sp; e.curve = sc; e.span = 0;
					if (sp < desc->path_count && sc < desc->paths[sp].curve_count)
					{
						qaws_curve const* src = desc->paths[sp].curves[sc];
						qaws_vec2 mid;
						qaws_scalar t = 0;
						mid.x = (qaws_scalar)(0.5 * (a[0] + b[0]));
						mid.y = (qaws_scalar)(0.5 * (a[1] + b[1]));
						if (qaws_curve_find_closest_parameter_2d(src, mid, &t) == QAWS_STATUS_OK)
							e.span = sk_span(src, (double)t);
					}
					if (desc->weight_fn)
					{
						e.w = (double)desc->weight_fn(desc->weight_user, e.path, e.curve, e.span);
						if (!(e.w > 0))
							st = QAWS_STATUS_INVALID_ARGUMENT;
					}
					if (st == QAWS_STATUS_OK)
						st = sk_add_support(s, &capsup, &e, a);
					if (a[0] < lo[0]) lo[0] = a[0];
					if (a[1] < lo[1]) lo[1] = a[1];
					if (a[0] > hi[0]) hi[0] = a[0];
					if (a[1] > hi[1]) hi[1] = a[1];
				}
				if (cp != buf) free(cp);
			}
			else if (c->kind == QAWS_CURVE_KIND_ARC && c->dimension == QAWS_DIMENSION_2D)
			{
				qaws_arc_impl const* impl = (qaws_arc_impl const*)c->impl;
				unsigned int j;
				for (j = 0; j < impl->segment_count && st == QAWS_STATUS_OK; j++)
				{
					qaws_arc_segment const* g = &impl->segments[j];
					double sw = (double)g->angle_end - (double)g->angle_start, a[2];
					sk_sup e;
					if (!(fabs(sw) > 0) || !((double)g->radius > 0))
						continue;
					memset(&e, 0, sizeof(e));
					e.kind = SK_CIRCLE;
					e.m[0] = (double)g->center[0]; e.m[1] = (double)g->center[1];
					e.r0 = (double)g->radius;
					e.s = sw > 0 ? 1 : -1;
					e.w = 1;
					e.a0 = (double)g->angle_start;
					e.sweep = sw;
					e.input = 1;
					e.path = sp; e.curve = sc; e.span = j;
					a[0] = e.m[0] + e.r0 * cos(e.a0);
					a[1] = e.m[1] + e.r0 * sin(e.a0);
					st = sk_add_support(s, &capsup, &e, a);
					{
						/* the extent from points on the arc (a flat arc's circle is huge) */
						unsigned int q;
						for (q = 0; q <= 8; q++)
						{
							double ang = e.a0 + e.sweep * q / 8, x = e.m[0] + e.r0 * cos(ang), y = e.m[1] + e.r0 * sin(ang);
							if (x < lo[0]) lo[0] = x;
							if (y < lo[1]) lo[1] = y;
							if (x > hi[0]) hi[0] = x;
							if (y > hi[1]) hi[1] = y;
						}
					}
				}
			}
			else
			{
				/* any other curve: G1 arc splines within the curve tolerance */
				sk_prims pr;
				unsigned int j;
				memset(&pr, 0, sizeof(pr));
				fitted = 1;
				st = sk_fit_curve(c, fit_tol, &pr);
				for (j = 0; j < pr.n && st == QAWS_STATUS_OK; j++)
				{
					sk_prim const* g = &pr.g[j];
					sk_sup e;
					memset(&e, 0, sizeof(e));
					if (!g->arc)
					{
						double len = hypot(g->b[0] - g->a[0], g->b[1] - g->a[1]);
						e.kind = SK_LINE;
						e.n[0] = -(g->b[1] - g->a[1]) / len;
						e.n[1] = (g->b[0] - g->a[0]) / len;
						e.d = e.n[0] * g->a[0] + e.n[1] * g->a[1];
					}
					else
					{
						e.kind = SK_CIRCLE;
						e.m[0] = g->m[0]; e.m[1] = g->m[1];
						e.r0 = g->r;
						e.s = g->sweep > 0 ? 1 : -1;
						e.a0 = g->a0;
						e.sweep = g->sweep;
					}
					e.w = 1;
					e.input = 1;
					e.path = sp; e.curve = sc; e.span = j;
					st = sk_add_support(s, &capsup, &e, g->a);
					if (g->a[0] < lo[0]) lo[0] = g->a[0];
					if (g->a[1] < lo[1]) lo[1] = g->a[1];
					if (g->a[0] > hi[0]) hi[0] = g->a[0];
					if (g->a[1] > hi[1]) hi[1] = g->a[1];
				}
				free(pr.g);
			}
		}
		/* a full circle is one support; a polygon needs three */
		for (k = first; k < s->nsup && s->sup[k].kind != SK_CIRCLE; k++) {}
		if (st == QAWS_STATUS_OK && s->nsup - first < 3 && k == s->nsup)
		{
			s->nsup = first;
			s->nnodes = first_node;
		}
		if (st == QAWS_STATUS_OK && s->nsup > first)
		{
			for (k = first; k < s->nsup; k++)
				s->sup[k].nb = s->sup[k + 1 < s->nsup ? k + 1 : first].na;
			s->nloops++;
		}
	}
	qaws_clip_result_destroy(region);
	s->ninput = s->nsup;
	if (st == QAWS_STATUS_OK)
	{
		double ext = hypot(hi[0] - lo[0], hi[1] - lo[1]);
		s->loop_first[s->nloops] = s->nsup;
		/* fitted arcs are exact doubles: events among them resolve to the
		   double default even where float vertices need the float one */
		s->tol = desc->tolerance > 0 ? (double)desc->tolerance : (fitted ? 1e-9 : SK_TOL_REL) * (ext > 0 ? ext : 1);
		s->ext = ext > 0 ? ext : 1;
		/* the input's edges as curves, for the faces */
		s->edges = (qaws_curve**)calloc(s->ninput + 1, sizeof(qaws_curve*));
		if (!s->edges)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (pi = 0; pi < s->ninput && st == QAWS_STATUS_OK; pi++)
			st = sk_edge_curve(&s->sup[pi], &s->nodes[s->sup[pi].na], &s->nodes[s->sup[pi].nb], &s->edges[pi]);
	}
	if (st == QAWS_STATUS_OK)
	{
		/* a loop that never closed is a near-coincidence the tolerance cut
		   the wrong way: the default tolerance is tried coarser and finer */
#if QAWS_SCALAR_IS_FLOAT
		/* not below the vertices' rounding: their meetings come apart */
		static double const scales[2] = { 1, 8 };
#else
		static double const scales[4] = { 1, 0.125, 8, 1.0 / 64 };
#endif
		unsigned int attempt, nodes0 = s->nnodes, attempts = desc->tolerance > 0 ? 1 : (unsigned int)(sizeof(scales) / sizeof(scales[0]));
		double tol0 = s->tol;
		sk_sup* inputs = (sk_sup*)malloc(sizeof(sk_sup) * (s->ninput + 1));
		if (!inputs)
			attempts = 0, st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
			memcpy(inputs, s->sup, sizeof(sk_sup) * s->ninput);
		for (attempt = 0; attempt < attempts; attempt++)
		{
			if (attempt)
			{
				/* back to the input */
				sk_sup* fresh = (sk_sup*)malloc(sizeof(sk_sup) * (s->ninput + 1));
				if (!fresh)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				memcpy(fresh, inputs, sizeof(sk_sup) * s->ninput);
				free(s->sup);
				s->sup = fresh;
				s->nsup = s->ninput;
				free(s->face_first); s->face_first = NULL;
				free(s->face_nodes); s->face_nodes = NULL;
			}
			s->tol = tol0 * scales[attempt];
			s->nnodes = nodes0;
			s->narcs = 0;
			s->tmax = 0;
			st = sk_simulate(s, 1, 1e300, &r);
			if (st == QAWS_STATUS_OK)
			{
				unsigned int i;
				for (i = 0; i < r.nc; i++)
					if (r.c[i].alive)
					{
						st = QAWS_STATUS_NUMERICAL_FAILURE;   /* a loop that never closed */
						break;
					}
			}
			if (st == QAWS_STATUS_OK)
			{
				/* keep the splice circles */
				free(s->sup);
				s->sup = r.sup;
				s->nsup = r.nsup;
				r.sup = NULL;
			}
			sk_run_free(&r);
			if (st == QAWS_STATUS_OK)
				st = sk_faces(s);
			if (st != QAWS_STATUS_NUMERICAL_FAILURE)
				break;
		}
		free(inputs);
	}
	if (st != QAWS_STATUS_OK)
	{
		qaws_skeleton_destroy(s);
		return st;
	}
	*out = s;
	return QAWS_STATUS_OK;
}

unsigned int qaws_skeleton_get_node_count(qaws_skeleton const* s) { return s ? s->nnodes : 0; }
unsigned int qaws_skeleton_get_arc_count(qaws_skeleton const* s) { return s ? s->narcs : 0; }
unsigned int qaws_skeleton_get_face_count(qaws_skeleton const* s) { return s ? s->nsup : 0; }
qaws_scalar qaws_skeleton_get_max_time(qaws_skeleton const* s) { return s ? (qaws_scalar)s->tmax : 0; }

qaws_status qaws_skeleton_get_node(qaws_skeleton const* s, unsigned int i, qaws_skeleton_node* out)
{
	if (!s || !out || i >= s->nnodes)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = s->nodes[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_get_arc(qaws_skeleton const* s, unsigned int i, qaws_skeleton_arc* out)
{
	if (!s || !out || i >= s->narcs)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = s->arcs[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_skeleton_get_face(qaws_skeleton const* s, unsigned int i, qaws_skeleton_face* out)
{
	if (!s || !out || i >= s->nsup)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out->nodes = s->face_nodes + s->face_first[i];
	out->node_count = s->face_first[i + 1] - s->face_first[i];
	out->path = s->sup[i].path;
	out->curve = s->sup[i].curve;
	out->span = s->sup[i].span;
	out->weight = (qaws_scalar)s->sup[i].w;
	out->edge = i < s->ninput ? s->edges[i] : NULL;
	out->circle = s->sup[i].kind == SK_CIRCLE;
	out->center.x = (qaws_scalar)s->sup[i].m[0];
	out->center.y = (qaws_scalar)s->sup[i].m[1];
	out->radius = (qaws_scalar)s->sup[i].r0;
	out->side = s->sup[i].s;
	out->normal.x = (qaws_scalar)s->sup[i].n[0];
	out->normal.y = (qaws_scalar)s->sup[i].n[1];
	out->offset = (qaws_scalar)s->sup[i].d;
	return QAWS_STATUS_OK;
}

/* the point of arc i at u in [0, 1], and its time */
static void sk_arc_eval(qaws_skeleton const* s, unsigned int i, double u, double* p, double* t)
{
	sk_geo const* g = &s->geo[i];
	qaws_skeleton_node const* a = &s->nodes[s->arcs[i].node_a];
	qaws_skeleton_node const* b = &s->nodes[s->arcs[i].node_b];
	if (g->kind == SK_GEO_CIRCLE)
	{
		double ang = g->a0 + g->sweep * u;
		p[0] = g->m[0] + g->rho * cos(ang);
		p[1] = g->m[1] + g->rho * sin(ang);
		*t = (double)a->time;
		return;
	}
	if (g->kind == SK_GEO_CONIC)
	{
		double pts[4];
		*t = g->t0 + (g->t1 - g->t0) * u;
		sk_cross(&s->sup[g->in], &s->sup[g->out], *t, s->tol, pts);
		p[0] = pts[2 * g->branch];
		p[1] = pts[2 * g->branch + 1];
		if (u == 0) { p[0] = (double)a->position.x; p[1] = (double)a->position.y; }
		if (u == 1) { p[0] = (double)b->position.x; p[1] = (double)b->position.y; }
		return;
	}
	p[0] = (double)a->position.x + ((double)b->position.x - a->position.x) * u;
	p[1] = (double)a->position.y + ((double)b->position.y - a->position.y) * u;
	*t = (double)a->time + ((double)b->time - a->time) * u;
}

qaws_status qaws_skeleton_evaluate_arc(qaws_skeleton const* s, unsigned int i, qaws_scalar u, qaws_vec2* point, qaws_scalar* time)
{
	double p[2], t;
	if (!s || i >= s->narcs || !point)
		return QAWS_STATUS_INVALID_ARGUMENT;
	sk_arc_eval(s, i, (double)u, p, &t);
	point->x = (qaws_scalar)p[0];
	point->y = (qaws_scalar)p[1];
	if (time) *time = (qaws_scalar)t;
	return QAWS_STATUS_OK;
}

/* a conic piece of arc i over [u0, u1] as a rational quadratic: the end
   tangents meet at the middle control point, whose weight puts the curve
   through the conic's point on the line from the chord's middle */
static int sk_conic_piece(qaws_skeleton const* s, unsigned int i, double u0, double u1, double* cp, double* w)
{
	sk_geo const* g = &s->geo[i];
	double p0[2], p2[2], t0, t2, d0[2], d2[2], det, a, q[2], lo = u0, hi = u1, f0, m[2];
	unsigned int it;
	sk_arc_eval(s, i, u0, p0, &t0);
	sk_arc_eval(s, i, u1, p2, &t2);
	{
		/* tangents from close points on either side */
		double h = 1e-7 * (u1 - u0) + 1e-12, x[2], tt;
		sk_arc_eval(s, i, u0 + h, x, &tt);
		d0[0] = x[0] - p0[0]; d0[1] = x[1] - p0[1];
		sk_arc_eval(s, i, u1 - h, x, &tt);
		d2[0] = p2[0] - x[0]; d2[1] = p2[1] - x[1];
	}
	(void)g;
	det = d0[0] * d2[1] - d0[1] * d2[0];
	if (!(fabs(det) > 1e-14 * hypot(d0[0], d0[1]) * hypot(d2[0], d2[1])))
		return 0;
	a = ((p2[0] - p0[0]) * d2[1] - (p2[1] - p0[1]) * d2[0]) / det;
	if (!(a > 0))
		return 0;
	cp[0] = p0[0]; cp[1] = p0[1];
	cp[2] = p0[0] + a * d0[0]; cp[3] = p0[1] + a * d0[1];
	cp[4] = p2[0]; cp[5] = p2[1];
	m[0] = 0.5 * (p0[0] + p2[0]); m[1] = 0.5 * (p0[1] + p2[1]);
	/* the shoulder: the conic's point on the segment m -> P1 */
	{
		double x[2], tt;
		sk_arc_eval(s, i, lo, x, &tt);
		f0 = (cp[2] - m[0]) * (x[1] - m[1]) - (cp[3] - m[1]) * (x[0] - m[0]);
		for (it = 0; it < 100; it++)
		{
			double mid = 0.5 * (lo + hi), fm;
			sk_arc_eval(s, i, mid, x, &tt);
			fm = (cp[2] - m[0]) * (x[1] - m[1]) - (cp[3] - m[1]) * (x[0] - m[0]);
			if ((fm > 0) == (f0 > 0)) { lo = mid; f0 = fm; } else hi = mid;
		}
		sk_arc_eval(s, i, 0.5 * (lo + hi), q, &tt);
	}
	{
		double ex = q[0] - m[0], ey = q[1] - m[1], fx = cp[2] - q[0], fy = cp[3] - q[1], nn = fx * fx + fy * fy;
		if (!(nn > 0))
			return 0;
		*w = (ex * fx + ey * fy) / nn;
	}
	return *w > 0;
}

qaws_status qaws_skeleton_get_arc_curve(qaws_skeleton const* s, unsigned int i, qaws_curve** out)
{
	sk_geo const* g;
	if (!s || i >= s->narcs || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	g = &s->geo[i];
	if (g->kind == SK_GEO_CIRCLE)
	{
		qaws_arc_segment a;
		qaws_arc_desc d;
		memset(&a, 0, sizeof(a));
		a.center[0] = (qaws_scalar)g->m[0]; a.center[1] = (qaws_scalar)g->m[1];
		a.radius = (qaws_scalar)g->rho;
		a.angle_start = (qaws_scalar)g->a0;
		a.angle_end = (qaws_scalar)(g->a0 + g->sweep);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = &a; d.segment_count = 1;
		return qaws_curve_create_arc(&d, out);
	}
	if (g->kind == SK_GEO_CONIC)
	{
		/* rational quadratic pieces, halved until each is one */
		double cps[3 * 2 * 16], ws[16], u[17];
		unsigned int n = 1, k;
		for (;;)
		{
			int ok = 1;
			for (k = 0; k <= n; k++) u[k] = (double)k / n;
			for (k = 0; k < n && ok; k++)
				ok = sk_conic_piece(s, i, u[k], u[k + 1], cps + 6 * k, &ws[k]);
			if (ok)
				break;
			if (n == 16)
			{
				n = 0;
				break;
			}
			n *= 2;
		}
		if (n)
		{
			qaws_scalar cp[2 * 33], kw[33], kn[36];
			qaws_nurbs_desc d;
			unsigned int j = 0;
			for (k = 0; k < n; k++)
			{
				if (k == 0) { cp[2 * j] = (qaws_scalar)cps[0]; cp[2 * j + 1] = (qaws_scalar)cps[1]; kw[j++] = 1; }
				cp[2 * j] = (qaws_scalar)cps[6 * k + 2]; cp[2 * j + 1] = (qaws_scalar)cps[6 * k + 3]; kw[j++] = (qaws_scalar)ws[k];
				cp[2 * j] = (qaws_scalar)cps[6 * k + 4]; cp[2 * j + 1] = (qaws_scalar)cps[6 * k + 5]; kw[j++] = 1;
			}
			/* clamped knots, each joint doubled */
			kn[0] = kn[1] = kn[2] = 0;
			for (k = 1; k < n; k++) { kn[1 + 2 * k] = kn[2 + 2 * k] = (qaws_scalar)((double)k / n); }
			kn[2 * n + 1] = kn[2 * n + 2] = kn[2 * n + 3] = 1;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.degree = 2;
			d.control_points = cp; d.control_point_count = j;
			d.knots = kn; d.knot_count = j + 3;
			d.weights = kw; d.weight_count = j;
			return qaws_curve_create_nurbs(&d, out);
		}
	}
	{
		qaws_scalar cp[4];
		qaws_bezier_desc d;
		qaws_skeleton_node const* a = &s->nodes[s->arcs[i].node_a];
		qaws_skeleton_node const* b = &s->nodes[s->arcs[i].node_b];
		cp[0] = a->position.x; cp[1] = a->position.y; cp[2] = b->position.x; cp[3] = b->position.y;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
		return qaws_curve_create_bezier(&d, out);
	}
}

/* ------------------------------------------------------------------------ */
/*  Offsets                                                                 */
/* ------------------------------------------------------------------------ */

qaws_status qaws_skeleton_offset(qaws_skeleton const* s, qaws_scalar time, qaws_clip_result** out)
{
	sk_run r;
	qaws_status st;
	unsigned int i, nloops = 0, ncur = 0;
	qaws_curve** curves = NULL;
	qaws_curve const** cv = NULL;
	unsigned int* lp = NULL;
	qaws_path_2d* paths = NULL;
	double T = (double)time;
	if (!s || !out || !(time >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = sk_simulate((qaws_skeleton*)s, 0, T, &r);
	if (st == QAWS_STATUS_OK)
	{
		unsigned int* seen = (unsigned int*)calloc(r.nc + 1, sizeof(unsigned int));
		curves = (qaws_curve**)calloc(r.nc + 1, sizeof(qaws_curve*));
		lp = (unsigned int*)malloc(sizeof(unsigned int) * (r.nc + 2));
		if (!seen || !curves || !lp)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < r.nc && st == QAWS_STATUS_OK; i++)
		{
			unsigned int j = i, n = 0, start = ncur;
			double area = 0;
			if (!r.c[i].alive || seen[i])
				continue;
			do
			{
				sk_corner const* c = &r.c[j];
				sk_sup const* e = &r.sup[c->out];
				double a[2], b[2];
				seen[j] = 1;
				sk_pos(&r, c, T, a);
				sk_pos(&r, &r.c[c->next], T, b);
				if (e->kind == SK_LINE)
				{
					if (hypot(b[0] - a[0], b[1] - a[1]) > 0)
					{
						qaws_scalar cp[4];
						qaws_bezier_desc d;
						cp[0] = (qaws_scalar)a[0]; cp[1] = (qaws_scalar)a[1]; cp[2] = (qaws_scalar)b[0]; cp[3] = (qaws_scalar)b[1];
						memset(&d, 0, sizeof(d));
						d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
						st = qaws_curve_create_bezier(&d, &curves[ncur]);
						if (st == QAWS_STATUS_OK) ncur++;
					}
					area += a[0] * b[1] - b[0] * a[1];
				}
				else
				{
					double rho = fabs(sk_rho(e, T)), sw = c->next == j ? 2 * SK_PI : sk_angle_along(e, a, b);
					if (rho > 0 && sw > 0)
					{
						qaws_arc_segment g;
						qaws_arc_desc d;
						memset(&g, 0, sizeof(g));
						g.center[0] = (qaws_scalar)e->m[0]; g.center[1] = (qaws_scalar)e->m[1];
						g.radius = (qaws_scalar)rho;
						g.angle_start = (qaws_scalar)atan2(a[1] - e->m[1], a[0] - e->m[0]);
						g.angle_end = (qaws_scalar)((double)g.angle_start + e->s * sw);
						memset(&d, 0, sizeof(d));
						d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
						st = qaws_curve_create_arc(&d, &curves[ncur]);
						if (st == QAWS_STATUS_OK) ncur++;
						area += e->s * sw * rho * rho + (a[0] * b[1] - b[0] * a[1]);
					}
				}
				j = c->next;
				n++;
			} while (j != i && n <= r.nc && st == QAWS_STATUS_OK);
			if (ncur > start && fabs(area) > 0)
				lp[nloops++] = start;
			else
			{
				while (ncur > start) qaws_curve_destroy(curves[--ncur]);
			}
		}
		lp[nloops] = ncur;
		free(seen);
	}
	sk_run_free(&r);
	if (st == QAWS_STATUS_OK)
	{
		cv = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (ncur + 1));
		paths = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (nloops + 1));
		if (!cv || !paths)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < ncur && st == QAWS_STATUS_OK; i++) cv[i] = curves[i];
		for (i = 0; i < nloops && st == QAWS_STATUS_OK; i++)
		{
			paths[i].curves = cv + lp[i];
			paths[i].curve_count = lp[i + 1] - lp[i];
			paths[i].closed = 1;
		}
		if (st == QAWS_STATUS_OK)
			st = qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_POSITIVE, paths, nloops, NULL, 0, out);
	}
	for (i = 0; curves && i < ncur; i++)
		if (curves[i]) qaws_curve_destroy(curves[i]);
	free(curves); free(lp); free(cv); free(paths);
	return st;
}
