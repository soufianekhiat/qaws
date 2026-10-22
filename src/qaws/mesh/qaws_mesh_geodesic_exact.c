/*
 * Exact polyhedral geodesic distances (qaws_mesh_geodesic.h): Chen and
 * Han's windows with Xin and Wang's filtering (ICH, 2009). A window is an
 * interval of an edge seen straight from a (pseudo) source unfolded into
 * the plane of the triangle the window enters; windows are propagated
 * cheapest first, each triangle's far edges taking the part of the window's
 * rays that crosses them; a saddle vertex (its angles summing to more than
 * 2 pi) or a border vertex reached is a new pseudo-source, expanded once
 * when the queue reaches its distance. A window is dropped when one of the
 * three vertices round it already reaches both its ends shorter; a window
 * splitting at a vertex no sooner than an earlier one through the same edge
 * keeps only the child away from that one (one angle, one split).
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_geodesic.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define GX_NONE 0xFFFFFFFFu
#define GX_EPS 1e-12

typedef struct gx_win
{
	unsigned int f, k;      /* the triangle it enters and its edge (corner k to k + 1), t measured from corner k's vertex */
	double t0, t1;          /* the interval */
	double sx, sy;          /* the source in the edge's frame (x along it, y into the triangle: sy <= 0) */
	double sigma;           /* the pseudo-source's own distance */
	double key;             /* the least distance over the interval */
} gx_win;

typedef struct gx_ctx
{
	unsigned int nv, nt;
	double* x;              /* 3 per vertex */
	unsigned int* t;        /* 3 per triangle */
	unsigned int* adj;      /* per triangle edge: the neighbour's edge 3 g + j, GX_NONE on the border */
	unsigned int *v0, *vc;  /* per vertex its corners (3 t + k) */
	unsigned char* saddle;  /* per vertex: angles summing to more than 2 pi, or on the border */
	double* d;              /* per vertex its distance */
	double *split_d, *split_t; /* per triangle edge: the shortest distance a window entering there gave the third vertex, where its ray crossed */
	gx_win* w;
	unsigned int nw, capw;
	unsigned int *heap, nh, caph;
} gx_ctx;

static double gx_len(double const* a, double const* b)
{
	return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

/* triangle f's edge k as a frame: its start U, its unit direction e, the unit n toward the third vertex in the plane, its
   length; the third vertex's coordinates (wx, wy) */
static void gx_frame(gx_ctx const* c, unsigned int f, unsigned int k, double* U, double* e, double* n, double* L, double* wx, double* wy)
{
	double const *u = c->x + 3 * (size_t)c->t[3 * f + k], *v = c->x + 3 * (size_t)c->t[3 * f + (k + 1) % 3], *w = c->x + 3 * (size_t)c->t[3 * f + (k + 2) % 3];
	double a[3], l, px, h;
	unsigned int i;
	for (i = 0; i < 3; i++) { U[i] = u[i]; e[i] = v[i] - u[i]; a[i] = w[i] - u[i]; }
	l = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	*L = l;
	for (i = 0; i < 3; i++) e[i] = l > 0 ? e[i] / l : 0;
	px = a[0] * e[0] + a[1] * e[1] + a[2] * e[2];
	for (i = 0; i < 3; i++) n[i] = a[i] - px * e[i];
	h = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (i = 0; i < 3; i++) n[i] = h > 0 ? n[i] / h : 0;
	*wx = px;
	*wy = h;
}

static double gx_dist(gx_win const* w, double t)
{
	return w->sigma + sqrt((t - w->sx) * (t - w->sx) + w->sy * w->sy);
}

static int gx_heap_less(gx_ctx const* c, unsigned int a, unsigned int b)
{
	return c->w[a].key < c->w[b].key;
}

static int gx_push(gx_ctx* c, gx_win const* w)
{
	unsigned int i;
	if (c->nw == c->capw)
	{
		unsigned int cap = c->capw ? 2 * c->capw : 1024;
		gx_win* g = (gx_win*)qaws_internal_realloc(c->w, sizeof(gx_win) * (size_t)cap);
		if (!g)
			return 0;
		c->w = g;
		c->capw = cap;
	}
	if (c->nh == c->caph)
	{
		unsigned int cap = c->caph ? 2 * c->caph : 1024;
		unsigned int* g = (unsigned int*)qaws_internal_realloc(c->heap, sizeof(unsigned int) * (size_t)cap);
		if (!g)
			return 0;
		c->heap = g;
		c->caph = cap;
	}
	c->w[c->nw] = *w;
	i = c->nh++;
	c->heap[i] = c->nw++;
	while (i > 0 && gx_heap_less(c, c->heap[i], c->heap[(i - 1) / 2]))
	{
		unsigned int p = (i - 1) / 2, s = c->heap[i];
		c->heap[i] = c->heap[p];
		c->heap[p] = s;
		i = p;
	}
	return 1;
}

static unsigned int gx_pop(gx_ctx* c)
{
	unsigned int top = c->heap[0], i = 0;
	c->heap[0] = c->heap[--c->nh];
	for (;;)
	{
		unsigned int l = 2 * i + 1, r = l + 1, m = i, s;
		if (l < c->nh && gx_heap_less(c, c->heap[l], c->heap[m])) m = l;
		if (r < c->nh && gx_heap_less(c, c->heap[r], c->heap[m])) m = r;
		if (m == i)
			break;
		s = c->heap[i];
		c->heap[i] = c->heap[m];
		c->heap[m] = s;
		i = m;
	}
	return top;
}

/* A window from a source at (sx, sy) in triangle f's edge k's frame (sy <= 0: behind the edge), of distance sigma, entering
   f over [t0, t1]; dropped when one of f's vertices reaches both its ends shorter (Xin and Wang's filter). 0 on allocation
   failure. */
static int gx_emit(gx_ctx* c, unsigned int f, unsigned int k, double sx, double sy, double sigma, double t0, double t1)
{
	double U[3], e[3], n[3], L, wx, wy, d0, d1;
	gx_win w;
	unsigned int a = c->t[3 * f + k], b = c->t[3 * f + (k + 1) % 3], z = c->t[3 * f + (k + 2) % 3];
	gx_frame(c, f, k, U, e, n, &L, &wx, &wy);
	/* (a source on the edge's line sends no ray into the triangle) */
	if (!(L > 0) || !(t1 - t0 > GX_EPS * L) || !(sy < -GX_EPS * L))
		return 1;
	w.f = f;
	w.k = k;
	w.t0 = t0 < 0 ? 0 : t0;
	w.t1 = t1 > L ? L : t1;
	w.sx = sx;
	w.sy = sy;
	w.sigma = sigma;
	d0 = gx_dist(&w, w.t0);
	d1 = gx_dist(&w, w.t1);
	if ((c->d[a] + w.t0 < d0 && c->d[a] + w.t1 < d1) || (c->d[b] + (L - w.t1) < d1 && c->d[b] + (L - w.t0) < d0))
		return 1;
	{
		double z0 = sqrt((wx - w.t0) * (wx - w.t0) + wy * wy), z1 = sqrt((wx - w.t1) * (wx - w.t1) + wy * wy);
		if (c->d[z] + z0 < d0 && c->d[z] + z1 < d1)
			return 1;
	}
	w.key = w.sx >= w.t0 && w.sx <= w.t1 ? sigma - w.sy : fmin(d0, d1);
	return gx_push(c, &w);
}

/* the window across triangle f's edge k (entering the neighbour) over [tf0, tf1] along it, from S (a point of f's plane):
   unfolded about the edge, S's place along it measured from the other end and its height off it turned to the other side */
static int gx_emit_across(gx_ctx* c, unsigned int f, unsigned int k, double const* S, double sigma, double tf0, double tf1)
{
	unsigned int o = c->adj[3 * f + k], i;
	double U[3], e[3], n[3], L, wx, wy, u = 0, v = 0;
	if (o == GX_NONE)
		return 1;
	gx_frame(c, f, k, U, e, n, &L, &wx, &wy);
	for (i = 0; i < 3; i++)
	{
		u += (S[i] - U[i]) * e[i];
		v += (S[i] - U[i]) * n[i];
	}
	return gx_emit(c, o / 3, o % 3, L - u, -v, sigma, L - tf1, L - tf0);
}

/* vertex v reached at distance dist: kept if shorter; a saddle or border vertex then a pseudo-source, its expansion queued
   at its distance (f GX_NONE, k the vertex) so that a vertex reached shorter before then is expanded once. 0 on allocation
   failure. */
static int gx_reach(gx_ctx* c, unsigned int v, double dist)
{
	gx_win w;
	if (!(dist < c->d[v] - GX_EPS * (1 + dist)))
		return 1;
	c->d[v] = dist;
	if (!c->saddle[v])
		return 1;
	memset(&w, 0, sizeof w);
	w.f = GX_NONE;
	w.k = v;
	w.sigma = w.key = dist;
	return gx_push(c, &w);
}

/* a pseudo-source: its windows across the edges facing it in its triangles */
static int gx_expand(gx_ctx* c, unsigned int v, double dist)
{
	unsigned int j;
	for (j = c->v0[v]; j < c->v0[v + 1]; j++)
	{
		/* (the facing edge's ends are reached when its window is, in order: reaching them here would cascade) */
		unsigned int f = c->vc[j] / 3, k = c->vc[j] % 3, a = c->t[3 * f + (k + 1) % 3], b = c->t[3 * f + (k + 2) % 3];
		if (!gx_emit_across(c, f, (k + 1) % 3, c->x + 3 * (size_t)v, dist, 0, gx_len(c->x + 3 * (size_t)a, c->x + 3 * (size_t)b)))
			return 0;
	}
	return 1;
}

/* a window across its triangle: the third vertex if its ray crosses the interval, the far edges' parts the rays reach */
static int gx_propagate(gx_ctx* c, gx_win const* w)
{
	double U[3], e[3], n[3], L, wx, wy, S[3], sy = w->sy < -GX_EPS * 1e-3 ? w->sy : -GX_EPS * 1e-3;
	unsigned int f = w->f, k = w->k, i, ce, drop = 0;
	gx_frame(c, f, k, U, e, n, &L, &wx, &wy);
	/* the source back in space (unfolded into this triangle's plane) */
	for (i = 0; i < 3; i++) S[i] = U[i] + w->sx * e[i] + w->sy * n[i];
	/* the interval's ends reached */
	if (w->t0 <= GX_EPS * L && !gx_reach(c, c->t[3 * f + k], gx_dist(w, 0)))
		return 0;
	if (w->t1 >= L * (1 - GX_EPS) && !gx_reach(c, c->t[3 * f + (k + 1) % 3], gx_dist(w, L)))
		return 0;
	{
		/* the third vertex: its ray from the source crosses y = 0 at tw */
		double tw = w->sx + (wx - w->sx) * (-sy) / (wy - sy);
		if (tw >= w->t0 - GX_EPS * L && tw <= w->t1 + GX_EPS * L)
		{
			/* one angle, one split (Xin and Wang): a window reaching the third vertex no shorter than an earlier one through
			   this edge keeps only its child on the side away from that one's ray (the other's paths would cross it) */
			double dw = w->sigma + sqrt((wx - w->sx) * (wx - w->sx) + (wy - w->sy) * (wy - w->sy));
			size_t a = 3 * (size_t)f + k;
			if (dw < c->split_d[a])
			{
				c->split_d[a] = dw;
				c->split_t[a] = tw;
			}
			else if (dw > c->split_d[a] + GX_EPS * (1 + dw))
				drop = tw < c->split_t[a] ? 1u : tw > c->split_t[a] ? 2u : 0u;
			if (!gx_reach(c, c->t[3 * f + (k + 2) % 3], dw))
				return 0;
		}
	}
	/* the far edges: V -> W (edge k + 1) and W -> U (edge k + 2), each point X(l) = P + l (Q - P) seen through the interval */
	for (ce = 1; ce <= 2; ce++)
	{
		double P[2], Q[2], dx, dy, lam[2], lo, hi;
		unsigned int q;
		if (ce == drop)
			continue;
		if (ce == 1) { P[0] = L; P[1] = 0; Q[0] = wx; Q[1] = wy; }
		else { P[0] = wx; P[1] = wy; Q[0] = 0; Q[1] = 0; }
		dx = Q[0] - P[0];
		dy = Q[1] - P[1];
		{
			/* the edge's ends seen on y = 0 (the map from the edge to the line is monotonic: the edge lies above the source), the
			   window clipped to them, then mapped back (no pole there) */
			double tP = w->sx + (P[0] - w->sx) * (-sy) / (P[1] - sy), tQ = w->sx + (Q[0] - w->sx) * (-sy) / (Q[1] - sy);
			double tmin = fmin(tP, tQ), tmax = fmax(tP, tQ), ta = fmax(w->t0, tmin), tb = fmin(w->t1, tmax);
			if (!(tb > ta))
				continue;
			for (q = 0; q < 2; q++)
			{
				double tau = q ? tb : ta, den = -dx * sy - (tau - w->sx) * dy;
				lam[q] = fabs(den) > 1e-300 ? ((tau - w->sx) * (P[1] - sy) + (P[0] - w->sx) * sy) / den : (tau == tP ? 0 : 1);
			}
		}
		lo = fmax(0, fmin(lam[0], lam[1]));
		hi = fmin(1, fmax(lam[0], lam[1]));
		if (!(hi > lo))
			continue;
		{
			unsigned int ek = (k + ce) % 3;
			double Le = sqrt(dx * dx + dy * dy);
			/* (along this triangle's edge ek, from its start P: l Le) */
			if (!gx_emit_across(c, f, ek, S, w->sigma, lo * Le, hi * Le))
				return 0;
		}
	}
	return 1;
}

static int gx_cmp_key(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y ? 1 : 0;
}

qaws_status qaws_mesh_3d_geodesic_exact(qaws_scalar* out_distances, qaws_mesh_3d const* mesh, qaws_mesh_location const* sources, unsigned int source_count,
	qaws_scalar max_distance)
{
	gx_ctx c;
	qaws_scalar* p = NULL;
	double* ang = NULL;
	unsigned int i, k, s;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	double R = (double)max_distance;
	if (!out_distances || !mesh || (source_count && !sources) || !(R >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	c.nv = qaws_mesh_3d_get_vertex_count(mesh);
	c.nt = qaws_mesh_3d_get_triangle_count(mesh);
	for (s = 0; s < source_count; s++)
		if (sources[s].patch >= c.nt || !isfinite((double)sources[s].a) || !isfinite((double)sources[s].b))
			return QAWS_STATUS_INVALID_ARGUMENT;
	c.x = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nv + 1));
	c.t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.nt + 1));
	c.adj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.nt + 1));
	c.v0 = (unsigned int*)qaws_internal_zalloc((size_t)c.nv + 2, sizeof(unsigned int));
	c.vc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.nt + 1));
	c.saddle = (unsigned char*)qaws_internal_zalloc((size_t)c.nv + 1, 1);
	c.d = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)c.nv + 1));
	c.split_d = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nt + 1));
	c.split_t = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nt + 1));
	p = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	ang = (double*)qaws_internal_scratch_zalloc((size_t)c.nv + 1, sizeof(double));
	if (!c.x || !c.t || !c.adj || !c.v0 || !c.vc || !c.saddle || !c.d || !c.split_d || !c.split_t || !p || !ang)
		goto done;
	qaws_mesh_3d_get_positions(p, mesh);
	for (i = 0; i < 3 * c.nv; i++) c.x[i] = (double)p[i];
	qaws_mesh_3d_get_triangles(c.t, mesh);
	/* the neighbours: each directed edge's (smaller, larger) vertex pair and corner sorted, pairs of two matched */
	for (i = 0; i < 3 * c.nt; i++)
	{
		c.adj[i] = GX_NONE;
		c.split_d[i] = HUGE_VAL;
		c.split_t[i] = 0;
	}
	{
		unsigned long long* kc = (unsigned long long*)qaws_internal_scratch_malloc(sizeof(unsigned long long) * 2 * (3 * (size_t)c.nt + 1));
		if (!kc)
			goto done;
		for (i = 0; i < 3 * c.nt; i++)
		{
			unsigned long long a = c.t[i], b = c.t[3 * (i / 3) + (i % 3 + 1) % 3];
			kc[2 * (size_t)i] = ((a < b ? a : b) << 32) | (a < b ? b : a);
			kc[2 * (size_t)i + 1] = i;
		}
		qsort(kc, 3 * (size_t)c.nt, sizeof(unsigned long long) * 2, gx_cmp_key);
		for (i = 0; i + 1 < 3 * c.nt; i++)
			if (kc[2 * (size_t)i] == kc[2 * (size_t)i + 2] && (i + 2 >= 3 * c.nt || kc[2 * (size_t)i + 4] != kc[2 * (size_t)i]))
			{
				unsigned int a = (unsigned int)kc[2 * (size_t)i + 1], b = (unsigned int)kc[2 * (size_t)i + 3];
				c.adj[a] = b;
				c.adj[b] = a;
				i++;
			}
		qaws_internal_free(kc);
	}
	/* each vertex's corners, its angle sum; saddles and border vertices */
	for (i = 0; i < 3 * c.nt; i++) c.v0[c.t[i] + 1]++;
	for (i = 0; i < c.nv; i++) c.v0[i + 1] += c.v0[i];
	{
		unsigned int* fill = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)c.nv + 1));
		if (!fill)
			goto done;
		memcpy(fill, c.v0, sizeof(unsigned int) * c.nv);
		for (i = 0; i < 3 * c.nt; i++) c.vc[fill[c.t[i]]++] = i;
		qaws_internal_free(fill);
	}
	for (i = 0; i < 3 * c.nt; i++)
	{
		double const *v = c.x + 3 * (size_t)c.t[i], *a = c.x + 3 * (size_t)c.t[3 * (i / 3) + (i % 3 + 1) % 3], *b = c.x + 3 * (size_t)c.t[3 * (i / 3) + (i % 3 + 2) % 3];
		double u[3], w2[3], la, lb, cs;
		for (k = 0; k < 3; k++) { u[k] = a[k] - v[k]; w2[k] = b[k] - v[k]; }
		la = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
		lb = sqrt(w2[0] * w2[0] + w2[1] * w2[1] + w2[2] * w2[2]);
		cs = la > 0 && lb > 0 ? (u[0] * w2[0] + u[1] * w2[1] + u[2] * w2[2]) / (la * lb) : 1;
		ang[c.t[i]] += acos(cs > 1 ? 1 : cs < -1 ? -1 : cs);
		if (c.adj[i] == GX_NONE)
		{
			c.saddle[c.t[i]] = 1;
			c.saddle[c.t[3 * (i / 3) + (i % 3 + 1) % 3]] = 1;
		}
	}
	for (i = 0; i < c.nv; i++)
	{
		c.saddle[i] |= ang[i] > 2 * 3.14159265358979323846 + 1e-9;
		c.d[i] = HUGE_VAL;
	}
	/* the sources: their triangle's vertices reached straight, a window across each of its edges (a source on a vertex: that
	   vertex a pseudo-source) */
	for (s = 0; s < source_count; s++)
	{
		unsigned int f = sources[s].patch;
		double a = (double)sources[s].a, b = (double)sources[s].b, S[3];
		for (k = 0; k < 3; k++)
			S[k] = (1 - a - b) * c.x[3 * (size_t)c.t[3 * f] + k] + a * c.x[3 * (size_t)c.t[3 * f + 1] + k] + b * c.x[3 * (size_t)c.t[3 * f + 2] + k];
		for (k = 0; k < 3; k++)
		{
			unsigned int v = c.t[3 * f + k];
			double dv = gx_len(S, c.x + 3 * (size_t)v);
			if (dv <= GX_EPS)
			{
				/* (on the vertex: every triangle round it at once) */
				unsigned char was = c.saddle[v];
				c.saddle[v] = 1;
				if (!gx_reach(&c, v, 0))
					goto done;
				c.saddle[v] = was;
			}
			else if (!gx_reach(&c, v, dv))
				goto done;
			if (!gx_emit_across(&c, f, k, S, 0, 0, gx_len(c.x + 3 * (size_t)c.t[3 * f + k], c.x + 3 * (size_t)c.t[3 * f + (k + 1) % 3])))
				goto done;
		}
	}
	while (c.nh)
	{
		gx_win w = c.w[gx_pop(&c)];
		if (w.key > R)
			break;   /* (every window left is farther: the distances within R are settled) */
		if (w.f == GX_NONE ? w.sigma == c.d[w.k] && !gx_expand(&c, w.k, w.sigma) : !gx_propagate(&c, &w))
			goto done;
	}
	for (i = 0; i < c.nv; i++) out_distances[i] = (qaws_scalar)(c.d[i] <= R ? c.d[i] : HUGE_VAL);
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(c.x); qaws_internal_free(c.t); qaws_internal_free(c.adj); qaws_internal_free(c.v0); qaws_internal_free(c.vc);
	qaws_internal_free(c.saddle); qaws_internal_free(c.d); qaws_internal_free(c.split_d); qaws_internal_free(c.split_t); qaws_internal_free(c.w); qaws_internal_free(c.heap);
	qaws_internal_free(p); qaws_internal_free(ang);
	return st;
}
