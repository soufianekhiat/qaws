/*
 * The surface of a stack of contours as a 3D triangle mesh: sections
 * sampled on their curves, contours paired between heights and stitched,
 * fans at tips, caps where the topology changes and at the ends.
 */

#include "qaws_stack.h"
#include "qaws_triangulate.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_parallel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SM_NONE 0xFFFFFFFFu
#define SM_DENSE 32             /* samples per span to measure arc length */
#define STM_CORNER 0.9            /* a turn over about 25 degrees between dense points: a corner */

struct qaws_mesh_3d
{
	double* v;                  /* 3 per vertex */
	unsigned int nv, capv;
	unsigned int* t;            /* 3 per triangle */
	unsigned int nt, capt;
	int failed;
};

static unsigned int sm_vert(qaws_mesh_3d* m, double x, double y, double z)
{
	if (m->nv == m->capv)
	{
		unsigned int nc = m->capv ? 2 * m->capv : 1024;
		double* g = (double*)realloc(m->v, sizeof(double) * 3 * nc);
		if (!g) { m->failed = 1; return SM_NONE; }
		m->v = g;
		m->capv = nc;
	}
	m->v[3 * m->nv] = x; m->v[3 * m->nv + 1] = y; m->v[3 * m->nv + 2] = z;
	return m->nv++;
}

static void sm_tri(qaws_mesh_3d* m, unsigned int a, unsigned int b, unsigned int c)
{
	if (a == b || b == c || a == c) return;
	if (m->nt == m->capt)
	{
		unsigned int nc = m->capt ? 2 * m->capt : 2048;
		unsigned int* g = (unsigned int*)realloc(m->t, sizeof(unsigned int) * 3 * nc);
		if (!g) { m->failed = 1; return; }
		m->t = g;
		m->capt = nc;
	}
	m->t[3 * m->nt] = a; m->t[3 * m->nt + 1] = b; m->t[3 * m->nt + 2] = c;
	m->nt++;
}

/* ======================================================================== */
/*  Sections and their contours                                             */
/* ======================================================================== */

typedef struct sm_contour
{
	qaws_path_2d path;
	double length, area;        /* signed area: counter-clockwise outer, clockwise hole */
	double c[2], lo[2], hi[2];
	unsigned int n;             /* samples */
	unsigned int first;         /* its first vertex in the mesh */
	unsigned int up, down;      /* banded to the next / previous height */
	unsigned int* corner;      /* its corners (sample indices) */
	unsigned int nk;
	double area_up, area_dn;    /* the area banded to above / below (its share at a branch) */
} sm_contour;

typedef struct sm_section
{
	double z;
	qaws_clip_result* r;
	sm_contour* c;
	unsigned int nc;
} sm_section;

typedef struct sm_job
{
	qaws_stack_2d const* st;
	sm_section* s;
} sm_job;

static qaws_status sm_section_job(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sm_job* j = (sm_job*)ctx;
	unsigned int i;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		qaws_status s = qaws_stack_section_2d(&j->s[i].r, j->st, (qaws_scalar)j->s[i].z);
		if (s != QAWS_STATUS_OK) return s;
	}
	return QAWS_STATUS_OK;
}

static void sm_eval(qaws_curve const* c, double t, double* p)
{
	qaws_eval_result_2d r;
	qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION);
	p[0] = r.position.x;
	p[1] = r.position.y;
}

/* a contour's dense polyline: (curve, parameter, cumulative length) per point */
typedef struct sm_dense
{
	unsigned int* curve;
	double* t;
	double* s;
	double* p;
	unsigned int n;
} sm_dense;

static int sm_densify(qaws_path_2d const* p, sm_dense* d)
{
	unsigned int k, total = 1, n = 0;
	for (k = 0; k < p->curve_count; k++) total += SM_DENSE * p->curves[k]->span_count;
	d->curve = (unsigned int*)malloc(sizeof(unsigned int) * total);
	d->t = (double*)malloc(sizeof(double) * total);
	d->s = (double*)malloc(sizeof(double) * total);
	d->p = (double*)malloc(sizeof(double) * 2 * total);
	if (!d->curve || !d->t || !d->s || !d->p) return 0;
	for (k = 0; k < p->curve_count; k++)
	{
		qaws_curve const* c = p->curves[k];
		unsigned int m = SM_DENSE * c->span_count, j;
		double t0 = c->parameter_range.min_value, t1 = c->parameter_range.max_value;
		for (j = (k ? 1u : 0u); j <= m; j++)
		{
			double t = t0 + (t1 - t0) * j / m;
			d->curve[n] = k;
			d->t[n] = t;
			sm_eval(c, t, d->p + 2 * n);
			d->s[n] = n ? d->s[n - 1] + hypot(d->p[2 * n] - d->p[2 * n - 2], d->p[2 * n + 1] - d->p[2 * n - 1]) : 0;
			n++;
		}
	}
	d->n = n;
	return 1;
}

static void sm_dense_free(sm_dense* d)
{
	free(d->curve); free(d->t); free(d->s); free(d->p);
	memset(d, 0, sizeof(*d));
}

static int sm_measure(sm_contour* c)
{
	sm_dense d;
	unsigned int i;
	double a = 0, cx = 0, cy = 0;
	memset(&d, 0, sizeof(d));
	if (!sm_densify(&c->path, &d)) { sm_dense_free(&d); return 0; }
	c->length = d.n ? d.s[d.n - 1] : 0;
	c->lo[0] = c->lo[1] = 1e300;
	c->hi[0] = c->hi[1] = -1e300;
	for (i = 0; i < d.n; i++)
	{
		double const* p = d.p + 2 * i;
		double const* q = d.p + 2 * ((i + 1) % d.n);
		double w = p[0] * q[1] - q[0] * p[1];
		a += w;
		cx += (p[0] + q[0]) * w;
		cy += (p[1] + q[1]) * w;
		if (p[0] < c->lo[0]) c->lo[0] = p[0];
		if (p[1] < c->lo[1]) c->lo[1] = p[1];
		if (p[0] > c->hi[0]) c->hi[0] = p[0];
		if (p[1] > c->hi[1]) c->hi[1] = p[1];
	}
	c->area = 0.5 * a;
	if (fabs(a) > 0) { c->c[0] = cx / (3 * a); c->c[1] = cy / (3 * a); }
	else { c->c[0] = 0.5 * (c->lo[0] + c->hi[0]); c->c[1] = 0.5 * (c->lo[1] + c->hi[1]); }
	sm_dense_free(&d);
	return 1;
}

/* the point at arc length s along the dense polyline, on its curve */
static void sm_at(sm_contour const* c, sm_dense const* d, double s, double* p)
{
	unsigned int lo = 0, hi = d->n - 1;
	double f;
	while (hi - lo > 1)
	{
		unsigned int mid = (lo + hi) / 2;
		if (d->s[mid] <= s) lo = mid; else hi = mid;
	}
	f = d->s[hi] > d->s[lo] ? (s - d->s[lo]) / (d->s[hi] - d->s[lo]) : 0;
	if (d->curve[lo] == d->curve[hi])
		sm_eval(c->path.curves[d->curve[lo]], d->t[lo] + f * (d->t[hi] - d->t[lo]), p);
	else
	{
		p[0] = d->p[2 * lo] + f * (d->p[2 * hi] - d->p[2 * lo]);
		p[1] = d->p[2 * lo + 1] + f * (d->p[2 * hi + 1] - d->p[2 * lo + 1]);
	}
}

/* the contour's samples into the mesh: on its corners (turns over STM_CORNER
   between dense points, where pieces meet or polylines kink), and by arc
   length between them; about n in all (c->n is updated) */
static int sm_emit(qaws_mesh_3d* m, sm_contour* c, double z)
{
	sm_dense d;
	double* corner;
	unsigned int i, k = 0, total = 0;
	memset(&d, 0, sizeof(d));
	if (!sm_densify(&c->path, &d)) { sm_dense_free(&d); return 0; }
	corner = (double*)malloc(sizeof(double) * d.n);
	if (!corner) { sm_dense_free(&d); return 0; }
	/* d's last point closes onto its first */
	for (i = 0; d.n >= 3 && i + 1 < d.n; i++)
	{
		unsigned int a = i ? i - 1 : d.n - 2, b = i + 1;
		double ux = d.p[2 * i] - d.p[2 * a], uy = d.p[2 * i + 1] - d.p[2 * a + 1];
		double vx = d.p[2 * b] - d.p[2 * i], vy = d.p[2 * b + 1] - d.p[2 * i + 1];
		double lu = hypot(ux, uy), lv = hypot(vx, vy);
		if (lu <= 1e-12 * c->length || lv <= 1e-12 * c->length) continue;
		if ((ux * vx + uy * vy) / (lu * lv) < STM_CORNER) corner[k++] = d.s[i];
	}
	c->first = m->nv;
	if (k == 0)
	{
		for (i = 0; i < c->n; i++)
		{
			double p[2];
			sm_at(c, &d, c->length * i / c->n, p);
			if (sm_vert(m, p[0], p[1], z) == SM_NONE) { free(corner); sm_dense_free(&d); return 0; }
		}
	}
	else
	{
		unsigned int g;
		c->corner = (unsigned int*)malloc(sizeof(unsigned int) * k);
		if (!c->corner) { free(corner); sm_dense_free(&d); return 0; }
		c->nk = k;
		for (g = 0; g < k; g++)
		{
			double s0 = corner[g], len = g + 1 < k ? corner[g + 1] - s0 : c->length - s0 + corner[0];
			unsigned int segs = (unsigned int)floor(c->n * len / c->length + 0.5), j;
			if (segs < 1) segs = 1;
			for (j = 0; j < segs; j++)
			{
				double s = s0 + len * j / segs, p[2];
				if (s >= c->length) s -= c->length;
				sm_at(c, &d, s, p);
				if (sm_vert(m, p[0], p[1], z) == SM_NONE) { free(corner); sm_dense_free(&d); return 0; }
			}
			c->corner[g] = total;
			total += segs;
		}
		c->n = total;
	}
	free(corner);
	sm_dense_free(&d);
	return 1;
}

/* ======================================================================== */
/*  Pairing, bands, fans, caps                                              */
/* ======================================================================== */

static double sm_d2(qaws_mesh_3d const* m, unsigned int a, unsigned int b)
{
	double dx = m->v[3 * a] - m->v[3 * b], dy = m->v[3 * a + 1] - m->v[3 * b + 1], dz = m->v[3 * a + 2] - m->v[3 * b + 2];
	return dx * dx + dy * dy + dz * dz;
}

static double sm_d2xy(qaws_mesh_3d const* m, unsigned int a, unsigned int b)
{
	double dx = m->v[3 * a] - m->v[3 * b], dy = m->v[3 * a + 1] - m->v[3 * b + 1];
	return dx * dx + dy * dy;
}

/* a band between the closed loops P (below, n vertices) and Q (above, k),
   both in the contours' travel: Q's start nearest P's, then the greedy
   stitch (the shorter diagonal first), n + k triangles */
static void sm_stitch(qaws_mesh_3d* m, unsigned int const* P, unsigned int n, unsigned int const* Q, unsigned int k)
{
	unsigned int s = 0, i = 0, j = 0, t;
	double bd = 1e300;
	for (t = 0; t < k; t++)
	{
		double d = sm_d2xy(m, P[0], Q[t]);
		if (d < bd) { bd = d; s = t; }
	}
	while (i < n || j < k)
	{
		unsigned int p0 = P[i % n], p1 = P[(i + 1) % n], q0 = Q[(s + j) % k], q1 = Q[(s + j + 1) % k];
		int along_p = i == n ? 0 : (j == k ? 1 : sm_d2(m, p1, q0) <= sm_d2(m, p0, q1));
		if (along_p) { sm_tri(m, p0, p1, q0); i++; }
		else { sm_tri(m, p0, q1, q0); j++; }
	}
}

static unsigned int* sm_loop(sm_contour const* c)
{
	unsigned int* l = (unsigned int*)malloc(sizeof(unsigned int) * (c->n + 1)), i;
	if (l) for (i = 0; i < c->n; i++) l[i] = c->first + i;
	return l;
}

/* a strip between the open chains P (below, np vertices) and Q (above, nq):
   the greedy stitch, np + nq - 2 triangles */
static void sm_strip(qaws_mesh_3d* m, unsigned int const* P, unsigned int np, unsigned int const* Q, unsigned int nq)
{
	unsigned int i = 0, j = 0;
	while (i + 1 < np || j + 1 < nq)
	{
		int along_p = i + 1 == np ? 0 : (j + 1 == nq ? 1 : sm_d2(m, P[i + 1], Q[j]) <= sm_d2(m, P[i], Q[j + 1]));
		if (along_p) { sm_tri(m, P[i], P[i + 1], Q[j]); i++; }
		else { sm_tri(m, P[i], Q[j + 1], Q[j]); j++; }
	}
}

/* the chain of c's samples from corner g to corner g + 1 (both included) */
static unsigned int sm_chain(sm_contour const* c, unsigned int g, unsigned int* out)
{
	unsigned int a = c->corner[g], b = c->corner[(g + 1) % c->nk], len = 0;
	if (b <= a) b += c->n;
	for (; a <= b; a++) out[len++] = c->first + a % c->n;
	return len;
}

/* the band between p and q: corner to corner when they have as many corners
   (the creases run straight up), else one stitch round */
static void sm_band(qaws_mesh_3d* m, sm_contour const* p, sm_contour const* q)
{
	if (p->nk && p->nk == q->nk)
	{
		unsigned int *a = (unsigned int*)malloc(sizeof(unsigned int) * (p->n + 2));
		unsigned int *b = (unsigned int*)malloc(sizeof(unsigned int) * (q->n + 2));
		unsigned int r = 0, g, t;
		double bd = 1e300;
		if (!a || !b) { free(a); free(b); m->failed = 1; return; }
		for (t = 0; t < q->nk; t++)
		{
			/* over every pair: corners come in close twins at creases */
			double d = 0;
			for (g = 0; g < p->nk; g++) d += sm_d2xy(m, p->first + p->corner[g], q->first + q->corner[(g + t) % q->nk]);
			if (d < bd) { bd = d; r = t; }
		}
		for (g = 0; g < p->nk; g++)
		{
			unsigned int na = sm_chain(p, g, a), nb = sm_chain(q, (g + r) % q->nk, b);
			sm_strip(m, a, na, b, nb);
		}
		free(a);
		free(b);
		return;
	}
	{
		unsigned int *a = sm_loop(p), *b = sm_loop(q);
		if (a && b) sm_stitch(m, a, p->n, b, q->n);
		else m->failed = 1;
		free(a);
		free(b);
	}
}

/* is (x, y) inside contour c (crossings over its samples)? */
static int sm_inside(qaws_mesh_3d const* m, sm_contour const* c, double x, double y)
{
	unsigned int i;
	int in = 0;
	for (i = 0; i < c->n; i++)
	{
		double const* a = m->v + 3 * (c->first + (i + c->n - 1) % c->n);
		double const* b = m->v + 3 * (c->first + i);
		if ((a[1] > y) != (b[1] > y) && x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])) in = !in;
	}
	return in;
}

/* do contours a and b (neighbouring heights) overlap in plan: the same
   orientation, and a sample of one inside the other */
static int sm_overlap(qaws_mesh_3d const* m, sm_contour const* a, sm_contour const* b)
{
	unsigned int i, step;
	if ((a->area > 0) != (b->area > 0)) return 0;
	if (a->hi[0] < b->lo[0] || b->hi[0] < a->lo[0] || a->hi[1] < b->lo[1] || b->hi[1] < a->lo[1]) return 0;
	if (sm_inside(m, a, b->c[0], b->c[1]) || sm_inside(m, b, a->c[0], a->c[1])) return 1;
	step = b->n > 64 ? b->n / 64 : 1;
	for (i = 0; i < b->n; i += step)
		if (sm_inside(m, a, m->v[3 * (b->first + i)], m->v[3 * (b->first + i) + 1])) return 1;
	step = a->n > 64 ? a->n / 64 : 1;
	for (i = 0; i < a->n; i += step)
		if (sm_inside(m, b, m->v[3 * (a->first + i)], m->v[3 * (a->first + i) + 1])) return 1;
	return 0;
}

/* contour x splitting into (x below) or merging from (x above) the
   contours ys: x's samples go to their nearest partner, each run becomes a
   loop closed across the neck by a bridge that the next loop shares, each
   loop is banded to its partner; more than two partners also close the
   polygon of the cut points. Returns 0 when the runs do not give one per
   partner (the caller falls back). */
static int sm_branch(qaws_mesh_3d* m, sm_contour const* x, sm_contour const* const* ys, unsigned int k, int x_below)
{
	unsigned int* own = (unsigned int*)malloc(sizeof(unsigned int) * (x->n + 1));
	unsigned int* start = (unsigned int*)malloc(sizeof(unsigned int) * (x->n + 1));
	unsigned int* who = (unsigned int*)malloc(sizeof(unsigned int) * (x->n + 1));
	unsigned int i, j, runs = 0, ok = 1;
	if (!own || !start || !who) { free(own); free(start); free(who); m->failed = 1; return 0; }
	for (i = 0; i < x->n; i++)
	{
		double bd = 1e300;
		own[i] = 0;
		for (j = 0; j < k; j++)
		{
			unsigned int s, step = ys[j]->n > 256 ? ys[j]->n / 256 : 1;
			for (s = 0; s < ys[j]->n; s += step)
			{
				double d = sm_d2xy(m, x->first + i, ys[j]->first + s);
				if (d < bd) { bd = d; own[i] = j; }
			}
		}
	}
	/* the runs; flickers at the neck (more runs than partners) go to the run
	   before them, the shortest first */
	for (;;)
	{
		unsigned int shortest = SM_NONE, len_min = SM_NONE;
		runs = 0;
		for (i = 0; i < x->n; i++)
			if (own[i] != own[(i + x->n - 1) % x->n]) { start[runs] = i; who[runs] = own[i]; runs++; }
		if (runs <= k) break;
		for (i = 0; i < runs; i++)
		{
			unsigned int len = (start[(i + 1) % runs] + x->n - start[i]) % x->n;
			if (len < len_min) { len_min = len; shortest = i; }
		}
		for (i = 0; i < len_min; i++) own[(start[shortest] + i) % x->n] = who[(shortest + runs - 1) % runs];
	}
	/* one run per partner */
	if (runs != k) ok = 0;
	for (i = 0; ok && i < runs; i++)
		for (j = i + 1; j < runs; j++) if (who[i] == who[j]) ok = 0;
	if (ok)
	{
		unsigned int* loop = (unsigned int*)malloc(sizeof(unsigned int) * (x->n + 2));
		if (!loop) { m->failed = 1; ok = 0; }
		for (i = 0; ok && i < runs; i++)
		{
			unsigned int a = start[i], b = start[(i + 1) % runs], len = 0, v = a;
			unsigned int* yl = sm_loop(ys[who[i]]);
			/* x[a .. b], b the next run's first: its last edge b -> a is the bridge */
			for (;;)
			{
				loop[len++] = x->first + v;
				if (v == b) break;
				v = (v + 1) % x->n;
			}
			if (!yl) { m->failed = 1; break; }
			if (x_below) sm_stitch(m, loop, len, yl, ys[who[i]]->n);
			else sm_stitch(m, yl, ys[who[i]]->n, loop, len);
			free(yl);
		}
		if (ok && runs >= 3)
			for (i = 1; i + 1 < runs; i++)
			{
				unsigned int c0 = x->first + start[0], c1 = x->first + start[i], c2 = x->first + start[i + 1];
				if (x_below) sm_tri(m, c0, c1, c2);
				else sm_tri(m, c0, c2, c1);
			}
		free(loop);
	}
	free(own); free(start); free(who);
	return ok;
}

/* a branch's areas: the single contour sees the others' total, each of the
   others its share of the single one */
static void sm_shares(sm_section* a, unsigned int const* qa, unsigned int nqa, sm_section* b, unsigned int const* qb, unsigned int nqb)
{
	double ta = 0, tb = 0;
	unsigned int g;
	for (g = 0; g < nqa; g++) ta += a->c[qa[g]].area;
	for (g = 0; g < nqb; g++) tb += b->c[qb[g]].area;
	for (g = 0; g < nqa; g++)
	{
		a->c[qa[g]].up = 1;
		a->c[qa[g]].area_up = ta != 0 ? tb * a->c[qa[g]].area / ta : tb;
	}
	for (g = 0; g < nqb; g++)
	{
		b->c[qb[g]].down = 1;
		b->c[qb[g]].area_dn = tb != 0 ? ta * b->c[qb[g]].area / tb : ta;
	}
}

/* the surface between sections a (below) and b: overlapping contours
   grouped; one to one banded, one to several (a split) or several to one
   (a merge) branched, anything else paired by nearest centres. Linked
   contours are marked (up on a's, down on b's). */
static int sm_link(qaws_mesh_3d* m, sm_section* a, sm_section* b)
{
	unsigned int na = a->nc, nb = b->nc, i, j;
	unsigned char* ov;
	unsigned int *qa, *qb, *seen_a, *seen_b;
	if (!na || !nb) return 1;
	ov = (unsigned char*)calloc((size_t)na * nb, 1);
	qa = (unsigned int*)malloc(sizeof(unsigned int) * (na + 1));
	qb = (unsigned int*)malloc(sizeof(unsigned int) * (nb + 1));
	seen_a = (unsigned int*)calloc(na + 1, sizeof(unsigned int));
	seen_b = (unsigned int*)calloc(nb + 1, sizeof(unsigned int));
	if (!ov || !qa || !qb || !seen_a || !seen_b) { free(ov); free(qa); free(qb); free(seen_a); free(seen_b); return 0; }
	for (i = 0; i < na; i++)
		for (j = 0; j < nb; j++)
			ov[(size_t)i * nb + j] = (unsigned char)sm_overlap(m, &a->c[i], &b->c[j]);
	for (i = 0; i < na; i++)
	{
		unsigned int nqa = 0, nqb = 0, h;
		if (seen_a[i]) continue;
		/* the component of a's contour i */
		seen_a[i] = 1;
		qa[nqa++] = i;
		for (h = 0; h < nqa || h < nqb;)
		{
			unsigned int progress = 0;
			for (; h < nqa; h++)
				for (j = 0; j < nb; j++)
					if (ov[(size_t)qa[h] * nb + j] && !seen_b[j]) { seen_b[j] = 1; qb[nqb++] = j; progress = 1; }
			{
				unsigned int g;
				for (g = 0; g < nqb; g++)
					for (j = 0; j < na; j++)
						if (ov[(size_t)j * nb + qb[g]] && !seen_a[j]) { seen_a[j] = 1; qa[nqa++] = j; progress = 1; }
			}
			if (!progress) break;
		}
		if (nqb == 0) continue;
		if (nqa == 1 && nqb == 1)
		{
			sm_band(m, &a->c[qa[0]], &b->c[qb[0]]);
			a->c[qa[0]].up = 1; b->c[qb[0]].down = 1;
			a->c[qa[0]].area_up = b->c[qb[0]].area; b->c[qb[0]].area_dn = a->c[qa[0]].area;
			continue;
		}
		if (nqa == 1 || nqb == 1)
		{
			sm_contour const* ys[64];
			unsigned int k = nqa == 1 ? nqb : nqa, g;
			if (k <= 64)
			{
				for (g = 0; g < k; g++) ys[g] = nqa == 1 ? &b->c[qb[g]] : &a->c[qa[g]];
				if (sm_branch(m, nqa == 1 ? &a->c[qa[0]] : &b->c[qb[0]], ys, k, nqa == 1))
				{
					sm_shares(a, qa, nqa, b, qb, nqb);
					continue;
				}
			}
		}
		/* fallback: nearest centres, one to one; the rest are capped */
		{
			unsigned int g;
			for (g = 0; g < nqa; g++)
			{
				sm_contour* p = &a->c[qa[g]];
				unsigned int best = SM_NONE, t;
				double bd = 1e300;
				for (t = 0; t < nqb; t++)
				{
					sm_contour const* q = &b->c[qb[t]];
					double d = hypot(p->c[0] - q->c[0], p->c[1] - q->c[1]);
					if (q->down || !ov[(size_t)qa[g] * nb + qb[t]]) continue;
					if (d < bd) { bd = d; best = qb[t]; }
				}
				if (best == SM_NONE) continue;
				sm_band(m, p, &b->c[best]);
				p->up = 1; b->c[best].down = 1;
				p->area_up = b->c[best].area; b->c[best].area_dn = p->area;
			}
		}
	}
	free(ov); free(qa); free(qb); free(seen_a); free(seen_b);
	return !m->failed;
}

/* where contour c, left without a partner towards height z_next (dz away),
   closes: shrinking towards a point at the pace set by its partner on the
   other side (gap away) gives a tip (returned: its distance, at most dz);
   otherwise 0, a flat cap. Without a partner, a tip at z_next when nothing
   is there. */
static double sm_tip(sm_contour const* c, int linked_back, double area_back, double gap, double dz, int next_empty)
{
	double r, rate, reach;
	if (!linked_back) return next_empty ? dz : 0;
	r = sqrt(fabs(c->area));
	rate = (sqrt(fabs(area_back)) - r) / gap;
	if (rate <= 0) return 0;
	reach = r / rate;
	if (reach > 4 * dz) return 0;
	return reach < dz ? reach : dz;
}

/* a fan from contour p to its centre at height z (up: the tip above p;
   the triangles keep the contour's travel, so holes need no turn) */
static void sm_fan(qaws_mesh_3d* m, sm_contour const* p, double z, int up)
{
	unsigned int apex = sm_vert(m, p->c[0], p->c[1], z), i;
	if (apex == SM_NONE) return;
	for (i = 0; i < p->n; i++)
	{
		unsigned int a = p->first + i, b = p->first + (i + 1) % p->n;
		if (up) sm_tri(m, a, b, apex);
		else sm_tri(m, b, a, apex);
	}
}

/* the samples strictly between u and w along their contour (in either
   travel) when they all lie on the segment u w (a straight run the
   triangulation merged): their count, the first into *from, the step */
static unsigned int sm_between(qaws_mesh_3d const* m, sm_contour const* const* cs, unsigned int n, unsigned int u, unsigned int w,
	unsigned int* from, int* step)
{
	unsigned int c, dir;
	for (c = 0; c < n; c++)
	{
		unsigned int f = cs[c]->first, cn = cs[c]->n;
		if (u < f || u >= f + cn || w < f || w >= f + cn) continue;
		for (dir = 0; dir < 2; dir++)
		{
			unsigned int iu = u - f, iw = w - f, gap = dir ? (iu + cn - iw) % cn : (iw + cn - iu) % cn, j;
			double const* a = m->v + 3 * u;
			double const* b = m->v + 3 * w;
			double ex = b[0] - a[0], ey = b[1] - a[1], el = ex * ex + ey * ey;
			int on = gap > 1 && el > 0;
			for (j = 1; on && j < gap; j++)
			{
				double const* p = m->v + 3 * (f + (dir ? iu + cn - j : iu + j) % cn);
				double t = ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / el, cr = (p[0] - a[0]) * ey - (p[1] - a[1]) * ex;
				if (t <= 0 || t >= 1 || cr * cr > 1e-16 * el * el) on = 0;
			}
			if (on)
			{
				*from = f + (dir ? iu + cn - 1 : iu + 1) % cn;
				*step = dir ? -1 : 1;
				return gap - 1;
			}
		}
		return 0;
	}
	return 0;
}

/* a cap triangle, its edges split again through the samples on them (a
   fan from the opposite corner), so the cap meets the band edge for edge */
static void sm_cap_tri(qaws_mesh_3d* m, sm_contour const* const* cs, unsigned int n, unsigned int a, unsigned int b, unsigned int c)
{
	unsigned int t[3], e;
	t[0] = a; t[1] = b; t[2] = c;
	for (e = 0; e < 3; e++)
	{
		unsigned int p = t[e], q = t[(e + 1) % 3], r = t[(e + 2) % 3], from = 0, cnt, j, prev = p;
		int step = 1;
		cnt = sm_between(m, cs, n, p, q, &from, &step);
		if (!cnt) continue;
		for (j = 0; j < cnt; j++)
		{
			unsigned int c0 = 0, s;
			int local, cn;
			while (from < cs[c0]->first || from >= cs[c0]->first + cs[c0]->n) c0++;
			cn = (int)cs[c0]->n;
			local = ((int)(from - cs[c0]->first) + (int)j * step) % cn;
			s = cs[c0]->first + (unsigned int)(local < 0 ? local + cn : local);
			sm_cap_tri(m, cs, n, prev, s, r);
			prev = s;
		}
		sm_cap_tri(m, cs, n, prev, q, r);
		return;
	}
	sm_tri(m, a, b, c);
}

/* flat faces over the contours cs (their samples' rings) at their height:
   facing up when up, down otherwise; the triangulation uses the samples as
   its only vertices */
static int sm_cap(qaws_mesh_3d* m, sm_contour const* const* cs, unsigned int n, int up)
{
	qaws_curve** rings = (qaws_curve**)calloc(n + 1, sizeof(qaws_curve*));
	qaws_curve const** views = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (n + 1));
	qaws_path_2d* paths = (qaws_path_2d*)malloc(sizeof(qaws_path_2d) * (n + 1));
	qaws_triangulate_desc d;
	qaws_mesh_2d* mesh = NULL;
	unsigned int i, k, total = 0;
	int ok = 1;
	if (!rings || !views || !paths) ok = 0;
	for (i = 0; ok && i < n; i++)
	{
		qaws_scalar* xy = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * cs[i]->n);
		if (!xy) { ok = 0; break; }
		for (k = 0; k < cs[i]->n; k++)
		{
			xy[2 * k] = (qaws_scalar)m->v[3 * (cs[i]->first + k)];
			xy[2 * k + 1] = (qaws_scalar)m->v[3 * (cs[i]->first + k) + 1];
		}
		if (qaws_curve_create_polyline_2d(&rings[i], xy, cs[i]->n, 1) != QAWS_STATUS_OK) ok = 0;
		free(xy);
		views[i] = rings[i];
		paths[i].curves = &views[i];
		paths[i].curve_count = 1;
		paths[i].closed = 1;
		total += cs[i]->n;
	}
	if (ok)
	{
		memset(&d, 0, sizeof(d));
		d.paths = paths;
		d.path_count = n;
		d.fill_rule = QAWS_FILL_EVEN_ODD;
		d.tolerance = (qaws_scalar)1e300;
		if (qaws_triangulate_region_2d(&mesh, &d) != QAWS_STATUS_OK) ok = 0;
	}
	if (ok && mesh)
	{
		unsigned int nv = qaws_mesh_2d_get_vertex_count(mesh), nt = qaws_mesh_2d_get_triangle_count(mesh);
		unsigned int* map = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
		double z = m->v[3 * cs[0]->first + 2];
		if (!map) ok = 0;
		/* each vertex back to its sample (the same point up to the lattice) */
		for (i = 0; ok && i < nv; i++)
		{
			qaws_vec2 q = qaws_mesh_2d_get_vertex(mesh, i);
			unsigned int best = SM_NONE, c, s;
			double bd = 1e300;
			for (c = 0; c < n; c++)
				for (s = 0; s < cs[c]->n; s++)
				{
					unsigned int v = cs[c]->first + s;
					double dx = m->v[3 * v] - (double)q.x, dy = m->v[3 * v + 1] - (double)q.y, dd = dx * dx + dy * dy;
					if (dd < bd) { bd = dd; best = v; }
				}
			/* a point the triangulation made (where rings cross): its own vertex */
			if (best == SM_NONE || bd > 1e-18 * (1 + fabs((double)q.x) + fabs((double)q.y)))
				best = sm_vert(m, q.x, q.y, z);
			map[i] = best;
		}
		for (i = 0; ok && i < nt; i++)
		{
			unsigned int t[3];
			qaws_mesh_2d_get_triangle(t, mesh, i);
			if (up) sm_cap_tri(m, cs, n, map[t[0]], map[t[1]], map[t[2]]);
			else sm_cap_tri(m, cs, n, map[t[0]], map[t[2]], map[t[1]]);
		}
		free(map);
	}
	(void)total;
	qaws_mesh_2d_destroy(mesh);
	for (i = 0; rings && i < n; i++) qaws_curve_destroy(rings[i]);
	free(rings); free((void*)views); free(paths);
	return ok && !m->failed;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

qaws_status qaws_stack_mesh_3d(qaws_mesh_3d** out, qaws_stack_mesh_desc const* desc)
{
	qaws_stack_2d const* st;
	qaws_mesh_3d* m;
	sm_section* s = NULL;
	unsigned int ns = 0, i, j, per;
	double spacing, ext = 0;
	sm_job job;
	qaws_status status = QAWS_STATUS_OK;
	if (!desc || !out || !desc->stack) return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = desc->stack;
	m = (qaws_mesh_3d*)calloc(1, sizeof(qaws_mesh_3d));
	if (!m) return QAWS_STATUS_ALLOCATION_FAILURE;
	if (st->level_count == 0) { *out = m; return QAWS_STATUS_OK; }
	/* the heights: every level, and the sections between */
	per = desc->sections + 1;
	s = (sm_section*)calloc((st->level_count - 1) * per + 1, sizeof(sm_section));
	if (!s) { free(m); return QAWS_STATUS_ALLOCATION_FAILURE; }
	for (i = 0; i < st->level_count; i++)
	{
		s[ns++].z = (double)st->levels[i].z;
		if (i + 1 < st->level_count)
			for (j = 1; j < per; j++)
				s[ns++].z = (double)st->levels[i].z + ((double)st->levels[i + 1].z - (double)st->levels[i].z) * j / per;
	}
	job.st = st;
	job.s = s;
	status = qaws_internal_parallel(ns, 1, sm_section_job, desc->executor, &job);
	/* contours, measured */
	for (i = 0; status == QAWS_STATUS_OK && i < ns; i++)
	{
		unsigned int nc = s[i].r ? qaws_clip_result_get_path_count(s[i].r) : 0;
		s[i].c = (sm_contour*)calloc(nc + 1, sizeof(sm_contour));
		if (!s[i].c) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
		s[i].nc = nc;
		for (j = 0; j < nc; j++)
		{
			sm_contour* c = &s[i].c[j];
			qaws_clip_result_get_path(&c->path, s[i].r, j);
			c->up = c->down = 0;
			if (!sm_measure(c)) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
			if (c->hi[0] - c->lo[0] > ext) ext = c->hi[0] - c->lo[0];
			if (c->hi[1] - c->lo[1] > ext) ext = c->hi[1] - c->lo[1];
		}
	}
	spacing = desc->spacing > 0 ? (double)desc->spacing : (ext > 0 ? ext / 64 : 1);
	/* every contour sampled by its own length */
	for (i = 0; i < ns && status == QAWS_STATUS_OK; i++)
		for (j = 0; j < s[i].nc; j++)
		{
			unsigned int w = (unsigned int)ceil(s[i].c[j].length / spacing);
			s[i].c[j].n = w < 8 ? 8 : w;
			if (!sm_emit(m, &s[i].c[j], s[i].z)) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
		}
	/* bands between neighbouring heights (with splits and merges) */
	for (i = 0; status == QAWS_STATUS_OK && i + 1 < ns; i++)
		if (!sm_link(m, &s[i], &s[i + 1])) status = QAWS_STATUS_ALLOCATION_FAILURE;
	/* tips and caps for the contours left: a tip where the contour closes to
	   a point before the next height, a flat cap otherwise */
	for (i = 0; status == QAWS_STATUS_OK && i < ns; i++)
	{
		for (j = 0; j < s[i].nc; j++)
		{
			sm_contour const* c = &s[i].c[j];
			sm_contour const* one = c;
			int hole = c->area < 0;
			if (i + 1 < ns && !c->up)
			{
				double tip = sm_tip(c, i > 0 && c->down, c->area_dn, i > 0 ? s[i].z - s[i - 1].z : 1, s[i + 1].z - s[i].z, s[i + 1].nc == 0);
				if (tip > 0) sm_fan(m, c, s[i].z + tip, 1);
				else if (!sm_cap(m, &one, 1, !hole)) status = QAWS_STATUS_ALLOCATION_FAILURE;
			}
			if (i > 0 && !c->down && status == QAWS_STATUS_OK)
			{
				double tip = sm_tip(c, i + 1 < ns && c->up, c->area_up, i + 1 < ns ? s[i + 1].z - s[i].z : 1, s[i].z - s[i - 1].z, s[i - 1].nc == 0);
				if (tip > 0) sm_fan(m, c, s[i].z - tip, 0);
				else if (!sm_cap(m, &one, 1, hole)) status = QAWS_STATUS_ALLOCATION_FAILURE;
			}
		}
	}
	if (status == QAWS_STATUS_OK && desc->caps)
	{
		unsigned int e;
		for (e = 0; e < 2 && status == QAWS_STATUS_OK; e++)
		{
			sm_section const* sec = &s[e ? ns - 1 : 0];
			sm_contour const** cs;
			if (!sec->nc || (e == 1 && ns == 1)) continue;
			cs = (sm_contour const**)malloc(sizeof(sm_contour*) * sec->nc);
			if (!cs) { status = QAWS_STATUS_ALLOCATION_FAILURE; break; }
			for (j = 0; j < sec->nc; j++) cs[j] = &sec->c[j];
			if (!sm_cap(m, cs, sec->nc, e == 1)) status = QAWS_STATUS_ALLOCATION_FAILURE;
			free((void*)cs);
		}
	}
	if (m->failed && status == QAWS_STATUS_OK) status = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; s && i < ns; i++)
	{
		qaws_clip_result_destroy(s[i].r);
		for (j = 0; j < s[i].nc; j++) free(s[i].c[j].corner);
		free(s[i].c);
	}
	free(s);
	if (status != QAWS_STATUS_OK) { qaws_mesh_3d_destroy(m); return status; }
	*out = m;
	return QAWS_STATUS_OK;
}

void qaws_mesh_3d_destroy(qaws_mesh_3d* m)
{
	if (!m) return;
	free(m->v);
	free(m->t);
	free(m);
}

unsigned int qaws_mesh_3d_get_vertex_count(qaws_mesh_3d const* m) { return m ? m->nv : 0; }

qaws_vec3 qaws_mesh_3d_get_vertex(qaws_mesh_3d const* m, unsigned int i)
{
	qaws_vec3 r;
	r.x = (qaws_scalar)m->v[3 * i];
	r.y = (qaws_scalar)m->v[3 * i + 1];
	r.z = (qaws_scalar)m->v[3 * i + 2];
	return r;
}

unsigned int qaws_mesh_3d_get_triangle_count(qaws_mesh_3d const* m) { return m ? m->nt : 0; }

void qaws_mesh_3d_get_triangle(unsigned int out[3], qaws_mesh_3d const* m, unsigned int t)
{
	out[0] = m->t[3 * t];
	out[1] = m->t[3 * t + 1];
	out[2] = m->t[3 * t + 2];
}
