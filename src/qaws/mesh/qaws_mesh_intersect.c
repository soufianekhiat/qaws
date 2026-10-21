/*
 * Self-intersections resolved numerically (qaws_mesh_intersect.h): pairs
 * from the broad phase cut each other into keyed points and carriers (an
 * input edge, or the line where two faces cross: each keeps its points
 * sorted along it), every face settles the points lying on or crossing its
 * carriers, then each cut face is triangulated again in its plane.
 */

#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_broadphase.h"
#include "../qaws_mesh_intersect.h"
#include "../qaws_mesh_repair.h"
#include "../qaws_exact.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SX_NONE 0xFFFFFFFFu

/* key kinds: points made by an edge through a face's plane, three faces' planes, two edges; carriers of an edge, of two faces */
enum { SX_LPI = 1, SX_TPI, SX_EE, SX_EDGE, SX_CUT };

typedef struct sx_node
{
	unsigned int a, b, c, next;
} sx_node;

typedef struct sx_carrier
{
	unsigned int kind;                 /* SX_EDGE: points a < b; SX_CUT: faces a < b */
	unsigned int a, b;
	double dir[3];                     /* its points sorted along it */
	unsigned int points;               /* list: node.a the point */
	unsigned int faces;                /* faces triangulated through it: node.a */
} sx_carrier;

typedef struct sx_face
{
	unsigned int c[3];                 /* corner points */
	unsigned int e[3];                 /* carriers of the edges c0 c1, c1 c2, c2 c0 */
	unsigned int src;
	unsigned int cuts;                 /* list: node.a carrier, node.b and node.c its ends in the face */
	unsigned int touch;                /* isolated points: node.a */
	double n[3], d;                    /* unit normal and offset; n zero: degenerate */
	int queued;
} sx_face;

typedef struct sx_key
{
	unsigned int k[5];
	unsigned int id;
} sx_key;

typedef struct sx_ctx
{
	double* X;                         /* points */
	unsigned int* chain;               /* next point in the same spatial cell */
	unsigned int np, capx, capc;
	unsigned int* cell;                /* spatial hash heads */
	unsigned int ncell;                /* power of two */
	sx_key* keys;
	unsigned int nkey, capk;           /* capk power of two */
	sx_node* nodes;
	unsigned int nn, capn;
	sx_carrier* car;
	unsigned int ncar, capcar;
	sx_face* F;
	unsigned int nf;
	unsigned int* work;
	unsigned int nw, capw;
	double eps, inv;
	qaws_intersect_report* r;
	int oom;
} sx_ctx;

static int sx_reserve(void** p, unsigned int* cap, size_t need, size_t elem)
{
	void* q;
	unsigned int nc;
	if (need <= *cap)
		return 1;
	nc = qaws_internal_grow(*cap, need, 64);
	q = nc ? qaws_internal_realloc(*p, qaws_internal_mul2(nc, elem)) : NULL;
	if (!q)
		return 0;
	*p = q;
	*cap = nc;
	return 1;
}

static double sx_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void sx_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double sx_dist2(double const* a, double const* b)
{
	double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
	return x * x + y * y + z * z;
}

static unsigned int sx_hash(unsigned int const* k, unsigned int n)
{
	unsigned int h = 2166136261u, i;
	for (i = 0; i < n; i++)
		h = (h ^ k[i]) * 16777619u;
	return h ^ (h >> 15);
}

static unsigned int sx_cell_of(sx_ctx const* c, double const* p, int dx, int dy, int dz)
{
	long long q[3];
	unsigned int k[6];
	q[0] = (long long)floor(p[0] * c->inv) + dx;
	q[1] = (long long)floor(p[1] * c->inv) + dy;
	q[2] = (long long)floor(p[2] * c->inv) + dz;
	k[0] = (unsigned int)q[0];
	k[1] = (unsigned int)(q[0] >> 32);
	k[2] = (unsigned int)q[1];
	k[3] = (unsigned int)(q[1] >> 32);
	k[4] = (unsigned int)q[2];
	k[5] = (unsigned int)(q[2] >> 32);
	return sx_hash(k, 6) & (c->ncell - 1);
}

static void sx_rehash(sx_ctx* c)
{
	unsigned int i;
	for (i = 0; i < c->ncell; i++)
		c->cell[i] = SX_NONE;
	for (i = 0; i < c->np; i++)
	{
		unsigned int h = sx_cell_of(c, c->X + 3 * (size_t)i, 0, 0, 0);
		c->chain[i] = c->cell[h];
		c->cell[h] = i;
	}
}

/* the first point within the tolerance of p, else a new one */
static unsigned int sx_point(sx_ctx* c, double const* p)
{
	int dx, dy, dz;
	unsigned int id, h;
	double e2 = c->eps * c->eps;
	for (dx = -1; dx <= 1; dx++)
		for (dy = -1; dy <= 1; dy++)
			for (dz = -1; dz <= 1; dz++)
				for (id = c->cell[sx_cell_of(c, p, dx, dy, dz)]; id != SX_NONE; id = c->chain[id])
					if (sx_dist2(c->X + 3 * (size_t)id, p) <= e2)
						return id;
	if (!sx_reserve((void**)&c->X, &c->capx, 3 * ((size_t)c->np + 1), sizeof(double))
		|| !sx_reserve((void**)&c->chain, &c->capc, (size_t)c->np + 1, sizeof(unsigned int)))
	{
		c->oom = 1;
		return 0;
	}
	memcpy(c->X + 3 * (size_t)c->np, p, 3 * sizeof(double));
	h = sx_cell_of(c, p, 0, 0, 0);
	c->chain[c->np] = c->cell[h];
	c->cell[h] = c->np;
	c->np++;
	if (c->np > c->ncell && c->ncell < (1u << 30))
	{
		unsigned int* nc = (unsigned int*)qaws_internal_realloc(c->cell, sizeof(unsigned int) * 2 * (size_t)c->ncell);
		if (!nc)
			c->oom = 1;
		else
		{
			c->cell = nc;
			c->ncell *= 2;
			sx_rehash(c);
		}
	}
	return c->np - 1;
}

/* the slot of key k (its id SX_NONE when absent) */
static sx_key* sx_find(sx_ctx* c, unsigned int const* k)
{
	unsigned int h, i;
	if (2 * ((size_t)c->nkey + 1) > c->capk)
	{
		unsigned int cap = c->capk ? 2 * c->capk : 1024;
		sx_key* nk = (sx_key*)qaws_internal_malloc(qaws_internal_mul2(cap, sizeof(sx_key)));
		if (!nk)
		{
			c->oom = 1;
			return NULL;
		}
		for (i = 0; i < cap; i++)
			nk[i].id = SX_NONE;
		for (i = 0; i < c->capk; i++)
			if (c->keys[i].id != SX_NONE)
			{
				unsigned int j = sx_hash(c->keys[i].k, 5) & (cap - 1);
				while (nk[j].id != SX_NONE)
					j = (j + 1) & (cap - 1);
				nk[j] = c->keys[i];
			}
		qaws_internal_free(c->keys);
		c->keys = nk;
		c->capk = cap;
	}
	h = sx_hash(k, 5) & (c->capk - 1);
	while (c->keys[h].id != SX_NONE && memcmp(c->keys[h].k, k, sizeof c->keys[h].k))
		h = (h + 1) & (c->capk - 1);
	return c->keys + h;
}

static unsigned int sx_node_new(sx_ctx* c, unsigned int a, unsigned int b, unsigned int d, unsigned int next)
{
	if (!sx_reserve((void**)&c->nodes, &c->capn, (size_t)c->nn + 1, sizeof(sx_node)))
	{
		c->oom = 1;
		return next;
	}
	c->nodes[c->nn].a = a;
	c->nodes[c->nn].b = b;
	c->nodes[c->nn].c = d;
	c->nodes[c->nn].next = next;
	return c->nn++;
}

static void sx_queue(sx_ctx* c, unsigned int f)
{
	if (c->F[f].queued)
		return;
	if (!sx_reserve((void**)&c->work, &c->capw, (size_t)c->nw + 1, sizeof(unsigned int)))
	{
		c->oom = 1;
		return;
	}
	c->F[f].queued = 1;
	c->work[c->nw++] = f;
}

/* face f triangulated through carrier k */
static void sx_use(sx_ctx* c, unsigned int k, unsigned int f)
{
	unsigned int i;
	for (i = c->car[k].faces; i != SX_NONE; i = c->nodes[i].next)
		if (c->nodes[i].a == f)
			return;
	i = sx_node_new(c, f, 0, 0, c->car[k].faces);
	c->car[k].faces = i;
}

/* point p on carrier k, sorted along it; its faces queued when new */
static int sx_on(sx_ctx* c, unsigned int k, unsigned int p)
{
	sx_carrier* K;
	unsigned int *at, i, nd;
	double s;
	if (c->oom)
		return 0;
	K = c->car + k;
	for (i = K->points; i != SX_NONE; i = c->nodes[i].next)
		if (c->nodes[i].a == p)
			return 0;
	s = sx_dot(c->X + 3 * (size_t)p, K->dir);
	at = &K->points;
	while (*at != SX_NONE && sx_dot(c->X + 3 * (size_t)c->nodes[*at].a, K->dir) < s)
		at = &c->nodes[*at].next;
	nd = sx_node_new(c, p, 0, 0, *at);
	if (c->oom)
		return 0;
	K = c->car + k;
	at = &K->points;
	while (*at != SX_NONE && *at != c->nodes[nd].next)
		at = &c->nodes[*at].next;
	*at = nd;
	for (i = K->faces; i != SX_NONE; i = c->nodes[i].next)
		sx_queue(c, c->nodes[i].a);
	return 1;
}

static unsigned int sx_carrier_new(sx_ctx* c, unsigned int kind, unsigned int a, unsigned int b, double const* dir)
{
	sx_carrier* K;
	if (!sx_reserve((void**)&c->car, &c->capcar, (size_t)c->ncar + 1, sizeof(sx_carrier)))
	{
		c->oom = 1;
		return 0;
	}
	K = c->car + c->ncar;
	K->kind = kind;
	K->a = a;
	K->b = b;
	memcpy(K->dir, dir, sizeof K->dir);
	K->points = SX_NONE;
	K->faces = SX_NONE;
	return c->ncar++;
}

static unsigned int sx_edge(sx_ctx* c, unsigned int p, unsigned int q)
{
	unsigned int k[5] = { SX_EDGE, p < q ? p : q, p < q ? q : p, 0, 0 }, id;
	sx_key* s = sx_find(c, k);
	double dir[3];
	if (!s)
		return 0;
	if (s->id != SX_NONE)
		return s->id;
	dir[0] = c->X[3 * (size_t)k[2]] - c->X[3 * (size_t)k[1]];
	dir[1] = c->X[3 * (size_t)k[2] + 1] - c->X[3 * (size_t)k[1] + 1];
	dir[2] = c->X[3 * (size_t)k[2] + 2] - c->X[3 * (size_t)k[1] + 2];
	id = sx_carrier_new(c, SX_EDGE, k[1], k[2], dir);
	if (c->oom)
		return 0;
	memcpy(s->k, k, sizeof k);
	s->id = id;
	c->nkey++;
	sx_on(c, id, k[1]);
	sx_on(c, id, k[2]);
	return id;
}

static unsigned int sx_cut(sx_ctx* c, unsigned int f, unsigned int g)
{
	unsigned int k[5] = { SX_CUT, f < g ? f : g, f < g ? g : f, 0, 0 }, id;
	sx_key* s = sx_find(c, k);
	double dir[3], l;
	if (!s)
		return 0;
	if (s->id != SX_NONE)
		return s->id;
	sx_cross(c->F[k[1]].n, c->F[k[2]].n, dir);
	l = sqrt(sx_dot(dir, dir));
	dir[0] /= l;
	dir[1] /= l;
	dir[2] /= l;
	id = sx_carrier_new(c, SX_CUT, k[1], k[2], dir);
	if (c->oom)
		return 0;
	memcpy(s->k, k, sizeof k);
	s->id = id;
	c->nkey++;
	sx_use(c, id, f);
	sx_use(c, id, g);
	return id;
}

/* the point keyed by k: computed once (canonically from the key), else the fallback position */
static unsigned int sx_keyed(sx_ctx* c, unsigned int const* k, double const* fallback)
{
	sx_key* s = sx_find(c, k);
	double p[3];
	unsigned int id;
	if (!s)
		return 0;
	if (s->id != SX_NONE)
		return s->id;
	memcpy(p, fallback, sizeof p);
	if (k[0] == SX_LPI)
	{
		/* edge k1 k2 through the plane of face k3 */
		double const *a = c->X + 3 * (size_t)k[1], *b = c->X + 3 * (size_t)k[2];
		sx_face const* F = c->F + k[3];
		double da = sx_dot(F->n, a) - F->d, db = sx_dot(F->n, b) - F->d;
		if (fabs(da - db) > 0)
		{
			double t = fmin(fmax(da / (da - db), 0.0), 1.0);
			p[0] = a[0] + t * (b[0] - a[0]);
			p[1] = a[1] + t * (b[1] - a[1]);
			p[2] = a[2] + t * (b[2] - a[2]);
		}
	}
	else if (k[0] == SX_TPI)
	{
		sx_face const *A = c->F + k[1], *B = c->F + k[2], *C = c->F + k[3];
		double bc[3], ca[3], ab[3], det;
		sx_cross(B->n, C->n, bc);
		sx_cross(C->n, A->n, ca);
		sx_cross(A->n, B->n, ab);
		det = sx_dot(A->n, bc);
		if (fabs(det) > 1e-9)
		{
			unsigned int i;
			for (i = 0; i < 3; i++)
				p[i] = (A->d * bc[i] + B->d * ca[i] + C->d * ab[i]) / det;
			if (sx_dist2(p, fallback) > 1e4 * c->eps * c->eps)
				memcpy(p, fallback, sizeof p);   /* nearly dependent planes: the crossing as seen */
		}
	}
	else if (k[0] == SX_EE)
	{
		/* the closest points of edges k1 k2 and k3 k4, their midpoint */
		double const *a = c->X + 3 * (size_t)k[1], *b = c->X + 3 * (size_t)k[2], *q = c->X + 3 * (size_t)k[3], *r = c->X + 3 * (size_t)k[4];
		double u[3], v[3], w[3], uu, uv, vv, uw, vw, den;
		unsigned int i;
		for (i = 0; i < 3; i++)
		{
			u[i] = b[i] - a[i];
			v[i] = r[i] - q[i];
			w[i] = a[i] - q[i];
		}
		uu = sx_dot(u, u);
		uv = sx_dot(u, v);
		vv = sx_dot(v, v);
		uw = sx_dot(u, w);
		vw = sx_dot(v, w);
		den = uu * vv - uv * uv;
		if (den > 1e-12 * uu * vv)
		{
			double sa = fmin(fmax((uv * vw - vv * uw) / den, 0.0), 1.0), t = fmin(fmax((uu * vw - uv * uw) / den, 0.0), 1.0);
			for (i = 0; i < 3; i++)
				p[i] = 0.5 * (a[i] + sa * u[i] + q[i] + t * v[i]);
		}
	}
	id = sx_point(c, p);
	if (c->oom)
		return 0;
	s = sx_find(c, k);   /* (the table may have grown) */
	if (!s)
		return 0;
	memcpy(s->k, k, sizeof s->k);
	s->id = id;
	c->nkey++;
	return id;
}

static unsigned int sx_lpi(sx_ctx* c, unsigned int edge, unsigned int f)
{
	unsigned int k[5] = { SX_LPI, c->car[edge].a, c->car[edge].b, f, 0 }, id;
	double mid[3];
	unsigned int i;
	for (i = 0; i < 3; i++)
		mid[i] = 0.5 * (c->X[3 * (size_t)k[1] + i] + c->X[3 * (size_t)k[2] + i]);
	id = sx_keyed(c, k, mid);
	sx_on(c, edge, id);
	return id;
}

static unsigned int sx_ee(sx_ctx* c, unsigned int e0, unsigned int e1, double const* fallback)
{
	unsigned int lo = e0 < e1 ? e0 : e1, hi = e0 < e1 ? e1 : e0;
	unsigned int k[5] = { SX_EE, c->car[lo].a, c->car[lo].b, c->car[hi].a, c->car[hi].b }, id = sx_keyed(c, k, fallback);
	sx_on(c, e0, id);
	sx_on(c, e1, id);
	return id;
}

static void sx_touch(sx_ctx* c, unsigned int f, unsigned int p)
{
	sx_face* F = c->F + f;
	unsigned int i;
	if (p == F->c[0] || p == F->c[1] || p == F->c[2])
		return;
	for (i = F->touch; i != SX_NONE; i = c->nodes[i].next)
		if (c->nodes[i].a == p)
			return;
	i = sx_node_new(c, p, 0, 0, c->F[f].touch);
	c->F[f].touch = i;
	sx_queue(c, f);
}

static void sx_add_cut(sx_ctx* c, unsigned int f, unsigned int k, unsigned int p, unsigned int q)
{
	unsigned int i;
	for (i = c->F[f].cuts; i != SX_NONE; i = c->nodes[i].next)
		if (c->nodes[i].a == k && ((c->nodes[i].b == p && c->nodes[i].c == q) || (c->nodes[i].b == q && c->nodes[i].c == p)))
			return;
	i = sx_node_new(c, k, p, q, c->F[f].cuts);
	c->F[f].cuts = i;
	sx_use(c, k, f);
	sx_on(c, k, p);
	sx_on(c, k, q);
	sx_queue(c, f);
}

/* where face g's boundary meets the plane of face f (signs of g's corners s): up to two points */
static unsigned int sx_section(sx_ctx* c, unsigned int g, unsigned int f, int const* s, unsigned int* out)
{
	sx_face const* G = c->F + g;
	unsigned int n = 0, k;
	for (k = 0; k < 3; k++)
		if (!s[k] && n < 2)
			out[n++] = G->c[k];
	for (k = 0; k < 3; k++)
		if (s[k] * s[(k + 1) % 3] < 0 && n < 2)
			out[n++] = sx_lpi(c, G->e[k], f);
	return n;
}

/* the edges of face g clipped to face f, in their common plane: cuts of f */
static int sx_coplanar(sx_ctx* c, unsigned int f, unsigned int g)
{
	sx_face const *A = c->F + f, *B = c->F + g;
	double const* n = A->n;
	unsigned int ax = fabs(n[0]) > fabs(n[1]) ? (fabs(n[0]) > fabs(n[2]) ? 0 : 2) : (fabs(n[1]) > fabs(n[2]) ? 1 : 2);
	unsigned int x = ax == 0 ? 1 : 0, y = ax == 2 ? 1 : 2, k, j, any = 0;
	double a2[3][2], o;
	for (k = 0; k < 3; k++)
	{
		a2[k][0] = c->X[3 * (size_t)A->c[k] + x];
		a2[k][1] = c->X[3 * (size_t)A->c[k] + y];
	}
	o = (a2[1][0] - a2[0][0]) * (a2[2][1] - a2[0][1]) - (a2[1][1] - a2[0][1]) * (a2[2][0] - a2[0][0]) > 0 ? 1 : -1;
	for (k = 0; k < 3; k++)
	{
		unsigned int pi = B->c[k], qi = B->c[(k + 1) % 3], j0 = SX_NONE, j1 = SX_NONE, P, Q;
		double const *p3 = c->X + 3 * (size_t)pi, *q3 = c->X + 3 * (size_t)qi;
		double p[2] = { p3[x], p3[y] }, q[2] = { q3[x], q3[y] }, t0 = 0, t1 = 1, len = sqrt(sx_dist2(p3, q3));
		int empty = 0;
		if (len <= c->eps)
			continue;
		for (j = 0; j < 3 && !empty; j++)
		{
			double const *u = a2[j], *v = a2[(j + 1) % 3];
			double ex = v[0] - u[0], ey = v[1] - u[1], el = sqrt(ex * ex + ey * ey);
			double hp, hq;
			if (el <= 0)
				continue;
			hp = o * (ex * (p[1] - u[1]) - ey * (p[0] - u[0])) / el;
			hq = o * (ex * (q[1] - u[1]) - ey * (q[0] - u[0])) / el;
			if (hp < -c->eps && hq < -c->eps)
				empty = 1;
			else if (hp < -c->eps || hq < -c->eps)
			{
				double t = hp / (hp - hq);
				if (hp < -c->eps && t > t0)
				{
					t0 = t;
					j0 = j;
				}
				else if (hq < -c->eps && t < t1)
				{
					t1 = t;
					j1 = j;
				}
			}
		}
		if (empty || (t1 - t0) * len < -c->eps)
			continue;
		{
			double m0[3], m1[3];
			unsigned int i;
			for (i = 0; i < 3; i++)
			{
				m0[i] = p3[i] + t0 * (q3[i] - p3[i]);
				m1[i] = p3[i] + t1 * (q3[i] - p3[i]);
			}
			P = j0 == SX_NONE ? pi : sx_ee(c, B->e[k], A->e[j0], m0);
			Q = j1 == SX_NONE ? qi : sx_ee(c, B->e[k], A->e[j1], m1);
		}
		if (c->oom)
			return 0;
		if (P == Q || (t1 - t0) * len <= c->eps)
			sx_touch(c, f, P);
		else
			sx_add_cut(c, f, B->e[k], P, Q);
		any = 1;
	}
	return any;
}

static qaws_status sx_pair(void* user, unsigned int f, unsigned int g)
{
	sx_ctx* c = (sx_ctx*)user;
	sx_face const *A = c->F + f, *B = c->F + g;
	int sa[3], sb[3], za = 0, zb = 0, pa = 0, pb = 0, na = 0, nb = 0;
	unsigned int shared = 0, i, j;
	if (c->oom)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (A->n[0] == 0 && A->n[1] == 0 && A->n[2] == 0)
		return QAWS_STATUS_OK;
	if (B->n[0] == 0 && B->n[1] == 0 && B->n[2] == 0)
		return QAWS_STATUS_OK;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			shared += A->c[i] == B->c[j];
	if (shared >= 3)
		return QAWS_STATUS_OK;
	for (i = 0; i < 3; i++)
	{
		double da = sx_dot(B->n, c->X + 3 * (size_t)A->c[i]) - B->d, db = sx_dot(A->n, c->X + 3 * (size_t)B->c[i]) - A->d;
		sa[i] = da > c->eps ? 1 : da < -c->eps ? -1 : 0;
		sb[i] = db > c->eps ? 1 : db < -c->eps ? -1 : 0;
		/* a shared corner is on both planes */
		for (j = 0; j < 3; j++)
		{
			if (A->c[i] == B->c[j])
				sa[i] = 0;
			if (B->c[i] == A->c[j])
				sb[i] = 0;
		}
		za += !sa[i];
		zb += !sb[i];
		pa += sa[i] > 0;
		na += sa[i] < 0;
		pb += sb[i] > 0;
		nb += sb[i] < 0;
	}
	if (pa == 3 || na == 3 || pb == 3 || nb == 3)
		return QAWS_STATUS_OK;
	if (za == 3 || zb == 3)
	{
		int hit;
		if (shared == 2)
		{
			/* coplanar on an edge: they overlap only when folded (both on one side of it) */
			double s0 = 0, s1 = 0;
			unsigned int oa = 0, ob = 0, a = SX_NONE, b = SX_NONE;
			double e[3], w0[3], w1[3], x0[3], x1[3];
			for (i = 0; i < 3; i++)
			{
				int in = A->c[i] == B->c[0] || A->c[i] == B->c[1] || A->c[i] == B->c[2];
				if (!in)
					oa = A->c[i];
				else if (a == SX_NONE)
					a = A->c[i];
				else
					b = A->c[i];
				if (B->c[i] != A->c[0] && B->c[i] != A->c[1] && B->c[i] != A->c[2])
					ob = B->c[i];
			}
			for (i = 0; i < 3; i++)
			{
				e[i] = c->X[3 * (size_t)b + i] - c->X[3 * (size_t)a + i];
				w0[i] = c->X[3 * (size_t)oa + i] - c->X[3 * (size_t)a + i];
				w1[i] = c->X[3 * (size_t)ob + i] - c->X[3 * (size_t)a + i];
			}
			sx_cross(e, w0, x0);
			sx_cross(e, w1, x1);
			s0 = sx_dot(x0, A->n);
			s1 = sx_dot(x1, A->n);
			if (s0 * s1 <= 0)
				return QAWS_STATUS_OK;
		}
		hit = sx_coplanar(c, f, g);
		hit |= sx_coplanar(c, g, f);
		if (hit)
		{
			c->r->pairs++;
			c->r->coplanar_pairs++;
		}
		return c->oom ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	}
	if (shared == 2)
		return QAWS_STATUS_OK;   /* two planes through a shared edge meet only along it */
	{
		unsigned int ia[2], ib[2], n0, n1, P, Q;
		double dir[3], l, s[4], sP, sQ;
		n0 = sx_section(c, f, g, sa, ia);
		n1 = sx_section(c, g, f, sb, ib);
		if (c->oom)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		if (!n0 || !n1)
			return QAWS_STATUS_OK;
		if (n0 == 1)
			ia[1] = ia[0];
		if (n1 == 1)
			ib[1] = ib[0];
		sx_cross(A->n, B->n, dir);
		l = sqrt(sx_dot(dir, dir));
		if (!(l > 0))
			return QAWS_STATUS_OK;
		s[0] = sx_dot(c->X + 3 * (size_t)ia[0], dir) / l;
		s[1] = sx_dot(c->X + 3 * (size_t)ia[1], dir) / l;
		s[2] = sx_dot(c->X + 3 * (size_t)ib[0], dir) / l;
		s[3] = sx_dot(c->X + 3 * (size_t)ib[1], dir) / l;
		if (s[0] > s[1])
		{
			unsigned int t = ia[0];
			double u = s[0];
			ia[0] = ia[1];
			ia[1] = t;
			s[0] = s[1];
			s[1] = u;
		}
		if (s[2] > s[3])
		{
			unsigned int t = ib[0];
			double u = s[2];
			ib[0] = ib[1];
			ib[1] = t;
			s[2] = s[3];
			s[3] = u;
		}
		P = s[0] >= s[2] ? ia[0] : ib[0];
		sP = fmax(s[0], s[2]);
		Q = s[1] <= s[3] ? ia[1] : ib[1];
		sQ = fmin(s[1], s[3]);
		if (sP > sQ + c->eps)
			return QAWS_STATUS_OK;
		if (P == Q || sQ - sP <= c->eps)
		{
			/* touching at a point (a shared corner touching is no meeting) */
			int corner_a = P == A->c[0] || P == A->c[1] || P == A->c[2], corner_b = P == B->c[0] || P == B->c[1] || P == B->c[2];
			if (corner_a && corner_b)
				return QAWS_STATUS_OK;
			sx_touch(c, f, P);
			sx_touch(c, g, P);
			c->r->pairs++;
		}
		else
		{
			unsigned int k = sx_cut(c, f, g);
			if (!c->oom)
			{
				sx_add_cut(c, f, k, P, Q);
				sx_add_cut(c, g, k, P, Q);
			}
			c->r->pairs++;
		}
	}
	return c->oom ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
}

/* ---- settling a face: points on its pieces registered, crossing pieces cut at keyed points ---- */

typedef struct sx_piece
{
	unsigned int k, p, q;
} sx_piece;

typedef struct sx_buf
{
	sx_piece* s;
	unsigned int ns, caps;
	unsigned int* p;
	unsigned int npt, capp;
} sx_buf;

static int sx_has(unsigned int const* a, unsigned int n, unsigned int v)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		if (a[i] == v)
			return 1;
	return 0;
}

static void sx_add_point(sx_ctx* c, sx_buf* b, unsigned int v)
{
	if (sx_has(b->p, b->npt, v))
		return;
	if (!sx_reserve((void**)&b->p, &b->capp, (size_t)b->npt + 1, sizeof(unsigned int)))
	{
		c->oom = 1;
		return;
	}
	b->p[b->npt++] = v;
}

/* carrier k's points from p to q, as pieces */
static void sx_expand(sx_ctx* c, sx_buf* b, unsigned int k, unsigned int p, unsigned int q)
{
	unsigned int i, in = 0, last = SX_NONE, end = SX_NONE;
	for (i = c->car[k].points; i != SX_NONE; i = c->nodes[i].next)
	{
		unsigned int v = c->nodes[i].a;
		if (!in && (v == p || v == q))
		{
			in = 1;
			end = v == p ? q : p;
			last = v;
			sx_add_point(c, b, v);
			if (p == q)
				return;
			continue;
		}
		if (!in)
			continue;
		sx_add_point(c, b, v);
		if (!sx_reserve((void**)&b->s, &b->caps, (size_t)b->ns + 1, sizeof(sx_piece)))
		{
			c->oom = 1;
			return;
		}
		b->s[b->ns].k = k;
		b->s[b->ns].p = last;
		b->s[b->ns].q = v;
		b->ns++;
		last = v;
		if (v == end)
			return;
	}
}

static void sx_gather(sx_ctx* c, sx_buf* b, unsigned int f)
{
	sx_face const* F = c->F + f;
	unsigned int k, i;
	b->ns = 0;
	b->npt = 0;
	for (k = 0; k < 3; k++)
		sx_add_point(c, b, F->c[k]);
	for (k = 0; k < 3; k++)
		sx_expand(c, b, F->e[k], F->c[k], F->c[(k + 1) % 3]);
	for (i = F->cuts; i != SX_NONE; i = c->nodes[i].next)
		sx_expand(c, b, c->nodes[i].a, c->nodes[i].b, c->nodes[i].c);
	for (i = F->touch; i != SX_NONE; i = c->nodes[i].next)
		sx_add_point(c, b, c->nodes[i].a);
}

/* the parameter of p's foot on segment ab and its squared distance */
static double sx_foot(double const* p, double const* a, double const* b, double* d2)
{
	double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, w[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
	double uu = sx_dot(u, u), t = uu > 0 ? sx_dot(u, w) / uu : 0, q[3];
	q[0] = a[0] + t * u[0] - p[0];
	q[1] = a[1] + t * u[1] - p[1];
	q[2] = a[2] + t * u[2] - p[2];
	*d2 = sx_dot(q, q);
	return t;
}

/* one pass over face f; returns whether a carrier gained a point */
static int sx_settle(sx_ctx* c, sx_buf* b, unsigned int f)
{
	unsigned int i, j;
	double e2 = c->eps * c->eps;
	int changed = 0;
	sx_gather(c, b, f);
	if (c->oom)
		return 0;
	for (i = 0; i < b->npt; i++)
	{
		unsigned int v = b->p[i];
		double const* x = c->X + 3 * (size_t)v;
		for (j = 0; j < b->ns; j++)
		{
			sx_piece const* s = b->s + j;
			double const *a = c->X + 3 * (size_t)s->p, *q = c->X + 3 * (size_t)s->q;
			double d2, t;
			if (v == s->p || v == s->q)
				continue;
			t = sx_foot(x, a, q, &d2);
			if (t > 0 && t < 1 && d2 <= e2 && sx_dist2(x, a) > e2 && sx_dist2(x, q) > e2)
				changed |= sx_on(c, s->k, v);
		}
	}
	for (i = 0; i < b->ns; i++)
		for (j = i + 1; j < b->ns; j++)
		{
			sx_piece const *s = b->s + i, *t = b->s + j;
			sx_carrier const *K = c->car + s->k, *L = c->car + t->k;
			double const *a = c->X + 3 * (size_t)s->p, *a1 = c->X + 3 * (size_t)s->q, *q = c->X + 3 * (size_t)t->p, *q1 = c->X + 3 * (size_t)t->q;
			double u[3], v[3], w[3], uu, uv, vv, uw, vw, den, su, tv, m[3], n[3], x[3];
			unsigned int key[5], id, h;
			if (s->k == t->k || s->p == t->p || s->p == t->q || s->q == t->p || s->q == t->q)
				continue;
			for (h = 0; h < 3; h++)
			{
				u[h] = a1[h] - a[h];
				v[h] = q1[h] - q[h];
				w[h] = a[h] - q[h];
			}
			uu = sx_dot(u, u);
			uv = sx_dot(u, v);
			vv = sx_dot(v, v);
			uw = sx_dot(u, w);
			vw = sx_dot(v, w);
			den = uu * vv - uv * uv;
			if (!(den > 1e-14 * uu * vv))
				continue;   /* parallel: overlaps show as points on pieces */
			su = (uv * vw - vv * uw) / den;
			tv = (uu * vw - uv * uw) / den;
			if (!(su > 0 && su < 1 && tv > 0 && tv < 1))
				continue;
			for (h = 0; h < 3; h++)
			{
				m[h] = a[h] + su * u[h];
				n[h] = q[h] + tv * v[h];
				x[h] = 0.5 * (m[h] + n[h]);
			}
			if (sx_dist2(m, n) > e2 || sx_dist2(x, a) <= e2 || sx_dist2(x, a1) <= e2 || sx_dist2(x, q) <= e2 || sx_dist2(x, q1) <= e2)
				continue;
			memset(key, 0, sizeof key);
			if (K->kind == SX_EDGE && L->kind == SX_EDGE)
			{
				sx_ee(c, s->k, t->k, x);
				changed = 1;
				continue;
			}
			if (K->kind == SX_CUT && L->kind == SX_CUT)
			{
				unsigned int fs[4] = { K->a, K->b, L->a, L->b }, g3[3], ng = 0, z;
				for (z = 0; z < 4; z++)
					if (!sx_has(g3, ng, fs[z]))
					{
						if (ng == 3)
							break;
						g3[ng++] = fs[z];
					}
				if (ng != 3)
					continue;
				/* sorted */
				if (g3[0] > g3[1]) { z = g3[0]; g3[0] = g3[1]; g3[1] = z; }
				if (g3[1] > g3[2]) { z = g3[1]; g3[1] = g3[2]; g3[2] = z; }
				if (g3[0] > g3[1]) { z = g3[0]; g3[0] = g3[1]; g3[1] = z; }
				key[0] = SX_TPI;
				key[1] = g3[0];
				key[2] = g3[1];
				key[3] = g3[2];
			}
			else
			{
				/* an edge through the other face of a cut */
				sx_carrier const* E = K->kind == SX_EDGE ? K : L;
				sx_carrier const* C = K->kind == SX_EDGE ? L : K;
				key[0] = SX_LPI;
				key[1] = E->a;
				key[2] = E->b;
				key[3] = C->a == f ? C->b : C->a;
			}
			id = sx_keyed(c, key, x);
			if (c->oom)
				return 0;
			changed |= sx_on(c, s->k, id);
			changed |= sx_on(c, t->k, id);
		}
	return changed;
}

/* ---- constrained Delaunay triangulation of a face in its reference triangle ---- */

typedef struct sx_tri
{
	unsigned int v[3], n[3];           /* n[i]: across edge v[i] v[i+1] */
	unsigned char fix[3];              /* that edge is a cut */
} sx_tri;

typedef struct sx_cdt
{
	double* P;                         /* 2 per local point: reference coordinates (exact on the edges), for orientations */
	double* Q;                         /* the same points in the face's own metric, for circles */
	sx_tri* t;
	unsigned int nt, capt, capp, capq;
} sx_cdt;

static int sx_o2(double const* P, unsigned int a, unsigned int b, unsigned int c)
{
	qaws_exact_sign s = QAWS_EXACT_ZERO;
	qaws_exact_orient2d(&s, NULL, P + 2 * (size_t)a, P + 2 * (size_t)b, P + 2 * (size_t)c);
	return (int)s;
}

static double sx_incircle(double const* P, unsigned int ia, unsigned int ib, unsigned int ic, unsigned int id)
{
	double const *a = P + 2 * (size_t)ia, *b = P + 2 * (size_t)ib, *c = P + 2 * (size_t)ic, *d = P + 2 * (size_t)id;
	double adx = a[0] - d[0], ady = a[1] - d[1], bdx = b[0] - d[0], bdy = b[1] - d[1], cdx = c[0] - d[0], cdy = c[1] - d[1];
	double ad = adx * adx + ady * ady, bd = bdx * bdx + bdy * bdy, cd = cdx * cdx + cdy * cdy;
	return adx * (bdy * cd - bd * cdy) - ady * (bdx * cd - bd * cdx) + ad * (bdx * cdy - bdy * cdx);
}

static void sx_back(sx_cdt* T, unsigned int nb, unsigned int from, unsigned int to)
{
	unsigned int k;
	if (nb == SX_NONE)
		return;
	for (k = 0; k < 3; k++)
		if (T->t[nb].n[k] == from)
		{
			T->t[nb].n[k] = to;
			return;
		}
}

static unsigned int sx_tri_new(sx_cdt* T)
{
	return T->nt++;   /* (capacity reserved for the face beforehand) */
}

static void sx_set(sx_tri* t, unsigned int a, unsigned int b, unsigned int c, unsigned int na, unsigned int nb, unsigned int nc,
	int fa, int fb, int fc)
{
	t->v[0] = a;
	t->v[1] = b;
	t->v[2] = c;
	t->n[0] = na;
	t->n[1] = nb;
	t->n[2] = nc;
	t->fix[0] = (unsigned char)fa;
	t->fix[1] = (unsigned char)fb;
	t->fix[2] = (unsigned char)fc;
}

/* flip the edge i of triangle u (t = a b c, its neighbour b a d): a d c and d b c */
static void sx_flip(sx_cdt* T, unsigned int u, unsigned int i)
{
	sx_tri* t = T->t + u;
	unsigned int w = t->n[i], j, a, b, cc, d, nbc, nca, nad, ndb;
	int fbc, fca, fad, fdb;
	sx_tri* s = T->t + w;
	for (j = 0; j < 3 && s->n[j] != u; j++)
		;
	a = t->v[i];
	b = t->v[(i + 1) % 3];
	cc = t->v[(i + 2) % 3];
	d = s->v[(j + 2) % 3];
	nbc = t->n[(i + 1) % 3];
	fbc = t->fix[(i + 1) % 3];
	nca = t->n[(i + 2) % 3];
	fca = t->fix[(i + 2) % 3];
	nad = s->n[(j + 1) % 3];
	fad = s->fix[(j + 1) % 3];
	ndb = s->n[(j + 2) % 3];
	fdb = s->fix[(j + 2) % 3];
	sx_set(t, a, d, cc, nad, w, nca, fad, 0, fca);
	sx_set(s, d, b, cc, ndb, nbc, u, fdb, fbc, 0);
	sx_back(T, nad, w, u);
	sx_back(T, nbc, u, w);
}

/* the opposite corner across edge i of triangle u: its local index in the neighbour, and the vertex */
static unsigned int sx_opp(sx_cdt const* T, unsigned int u, unsigned int i)
{
	unsigned int w = T->t[u].n[i], j;
	for (j = 0; j < 3 && T->t[w].n[j] != u; j++)
		;
	return T->t[w].v[(j + 2) % 3];
}

static int sx_should_flip(sx_cdt const* T, unsigned int u, unsigned int i)
{
	sx_tri const* t = T->t + u;
	unsigned int a, b, cc, d;
	if (t->n[i] == SX_NONE || t->fix[i])
		return 0;
	a = t->v[i];
	b = t->v[(i + 1) % 3];
	cc = t->v[(i + 2) % 3];
	d = sx_opp(T, u, i);
	if (sx_incircle(T->Q, a, b, cc, d) <= 0)
		return 0;
	/* the quad a d b c must be convex */
	return sx_o2(T->P, cc, d, a) * sx_o2(T->P, cc, d, b) < 0;
}

static void sx_legalize(sx_cdt* T, unsigned int* stack, unsigned int top, unsigned int guard)
{
	/* stack holds pairs (triangle, edge), the edge opposite the new point (at v[edge + 2]) */
	while (top && guard--)
	{
		unsigned int i = stack[--top], u = stack[--top];
		unsigned int w;
		if (!sx_should_flip(T, u, i))
			continue;
		w = T->t[u].n[i];
		sx_flip(T, u, i);
		/* u = a d c, w = d b c with c the new point: their edges a d (u, 0) and d b (w, 0) */
		stack[top++] = u;
		stack[top++] = 0;
		stack[top++] = w;
		stack[top++] = 0;
	}
}

/* insert local point p; triangles reserved by the caller; stack of 4 n entries */
static void sx_insert(sx_cdt* T, unsigned int p, unsigned int* stack, unsigned int guard)
{
	unsigned int u, k, top = 0;
	for (u = 0; u < T->nt; u++)
	{
		sx_tri* t = T->t + u;
		int o[3];
		for (k = 0; k < 3; k++)
			o[k] = sx_o2(T->P, t->v[k], t->v[(k + 1) % 3], p);
		if (o[0] < 0 || o[1] < 0 || o[2] < 0)
			continue;
		for (k = 0; k < 3 && o[k]; k++)
			;
		if (k == 3)
		{
			/* inside: a b c into a b p, b c p, c a p */
			unsigned int a = t->v[0], b = t->v[1], cc = t->v[2], nab = t->n[0], nbc = t->n[1], nca = t->n[2];
			int fab = t->fix[0], fbc = t->fix[1], fca = t->fix[2];
			unsigned int u1 = sx_tri_new(T), u2 = sx_tri_new(T);
			sx_set(T->t + u, a, b, p, nab, u1, u2, fab, 0, 0);
			sx_set(T->t + u1, b, cc, p, nbc, u2, u, fbc, 0, 0);
			sx_set(T->t + u2, cc, a, p, nca, u, u1, fca, 0, 0);
			sx_back(T, nbc, u, u1);
			sx_back(T, nca, u, u2);
			stack[top++] = u;
			stack[top++] = 0;
			stack[top++] = u1;
			stack[top++] = 0;
			stack[top++] = u2;
			stack[top++] = 0;
		}
		else
		{
			/* on edge k (a b, c opposite; neighbour b a d): a p c, p b c, b p d, p a d */
			unsigned int a = t->v[k], b = t->v[(k + 1) % 3], cc = t->v[(k + 2) % 3], w = t->n[k];
			unsigned int nbc = t->n[(k + 1) % 3], nca = t->n[(k + 2) % 3], ub = sx_tri_new(T), wa = SX_NONE;
			int fix = t->fix[k], fbc = t->fix[(k + 1) % 3], fca = t->fix[(k + 2) % 3];
			if (w != SX_NONE)
			{
				sx_tri* s = T->t + w;
				unsigned int j, d, nad, ndb;
				int fad, fdb;
				for (j = 0; j < 3 && s->n[j] != u; j++)
					;
				d = s->v[(j + 2) % 3];
				nad = s->n[(j + 1) % 3];
				fad = s->fix[(j + 1) % 3];
				ndb = s->n[(j + 2) % 3];
				fdb = s->fix[(j + 2) % 3];
				wa = sx_tri_new(T);
				sx_set(T->t + w, b, p, d, ub, wa, ndb, fix, 0, fdb);
				sx_set(T->t + wa, p, a, d, u, nad, w, fix, fad, 0);
				sx_back(T, nad, w, wa);
				stack[top++] = w;
				stack[top++] = 2;
				stack[top++] = wa;
				stack[top++] = 1;
			}
			sx_set(T->t + u, a, p, cc, wa, ub, nca, fix, 0, fca);
			sx_set(T->t + ub, p, b, cc, w, nbc, u, fix, fbc, 0);
			sx_back(T, nbc, u, ub);
			stack[top++] = u;
			stack[top++] = 2;
			stack[top++] = ub;
			stack[top++] = 1;
		}
		break;
	}
	sx_legalize(T, stack, top, guard);   /* (every entry has p at v[edge + 2]) */
}

static int sx_has_edge(sx_cdt* T, unsigned int a, unsigned int b)
{
	unsigned int u, k;
	for (u = 0; u < T->nt; u++)
		for (k = 0; k < 3; k++)
		{
			unsigned int p = T->t[u].v[k], q = T->t[u].v[(k + 1) % 3];
			if ((p == a && q == b) || (p == b && q == a))
			{
				unsigned int w = T->t[u].n[k], j;
				T->t[u].fix[k] = 1;
				if (w != SX_NONE)
					for (j = 0; j < 3; j++)
						if (T->t[w].n[j] == u)
							T->t[w].fix[j] = 1;
				return 1;
			}
		}
	return 0;
}

/* the cut a b recovered by flips (Sloan); 0 when it cannot be */
static int sx_recover(sx_cdt* T, unsigned int a, unsigned int b, unsigned int depth)
{
	unsigned int guard = 4 * T->nt + 16;
	while (guard--)
	{
		unsigned int u, k, found = 0;
		if (sx_has_edge(T, a, b))
			return 1;
		for (u = 0; u < T->nt && !found; u++)
			for (k = 0; k < 3 && !found; k++)
			{
				sx_tri const* t = T->t + u;
				unsigned int p = t->v[k], q = t->v[(k + 1) % 3];
				int o1, o2;
				if (p != a && p != b)
				{
					/* a point exactly on the open cut: split the cut there */
					if (!sx_o2(T->P, a, b, p) && depth < 8)
					{
						double const *A = T->P + 2 * (size_t)a, *B = T->P + 2 * (size_t)b, *X = T->P + 2 * (size_t)p;
						double dt = (X[0] - A[0]) * (B[0] - A[0]) + (X[1] - A[1]) * (B[1] - A[1]);
						double ll = (B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]);
						if (dt > 0 && dt < ll)
							return sx_recover(T, a, p, depth + 1) & sx_recover(T, p, b, depth + 1);
					}
				}
				if (t->n[k] == SX_NONE || p == a || p == b || q == a || q == b || u > t->n[k])
					continue;
				o1 = sx_o2(T->P, a, b, p);
				o2 = sx_o2(T->P, a, b, q);
				if (o1 * o2 >= 0)
					continue;
				if (sx_o2(T->P, p, q, a) * sx_o2(T->P, p, q, b) >= 0)
					continue;
				if (t->fix[k])
					return 0;   /* crosses another cut */
				{
					unsigned int cc = t->v[(k + 2) % 3], d = sx_opp(T, u, k);
					if (sx_o2(T->P, cc, d, p) * sx_o2(T->P, cc, d, q) < 0)
					{
						sx_flip(T, u, k);
						found = 1;
					}
				}
			}
		if (!found)
			return 0;
	}
	return 0;
}

/* Delaunay again away from the cuts */
static void sx_relax(sx_cdt* T)
{
	unsigned int pass, u, k;
	for (pass = 0; pass < 32; pass++)
	{
		unsigned int flips = 0;
		for (u = 0; u < T->nt; u++)
			for (k = 0; k < 3; k++)
				if (T->t[u].n[k] != SX_NONE && T->t[u].n[k] > u && sx_should_flip(T, u, k))
				{
					sx_flip(T, u, k);
					flips++;
				}
		if (!flips)
			return;
	}
}

/* ---- the whole ---- */

static qaws_status sx_bp_visit(void* user, unsigned int i, unsigned int j)
{
	return sx_pair(user, i, j);
}

static double sx_snap52(double t)
{
	t = fmin(fmax(t, 0.0), 1.0);
	return floor(t * 4503599627370496.0 + 0.5) / 4503599627370496.0;
}

qaws_status qaws_mesh_resolve_intersections(qaws_mesh_3d** out, qaws_intersect_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_intersect_desc const* desc)
{
	qaws_intersect_report local;
	sx_ctx c;
	sx_buf b;
	sx_cdt T;
	qaws_status st = QAWS_STATUS_OK;
	unsigned int *alias = NULL, *loc = NULL, *otri = NULL, *stack = NULL, *used = NULL;
	qaws_scalar *opos = NULL, *osrc = NULL;
	qaws_bp_box* box = NULL;
	unsigned int nt = 0, off = 0, f, k, i, ntri = 0, capt = 0, caps = 0, caps2 = 0, onv = 0;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, diag;
	memset(&c, 0, sizeof c);
	memset(&b, 0, sizeof b);
	memset(&T, 0, sizeof T);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!face_count || !face_vertices || (vertex_count && !positions) || (desc && !(desc->tolerance >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * vertex_count; i++)
		if (!isfinite((double)positions[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (f = 0; f < face_count; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3;
		if (s < 3)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (k = 0; k < s; k++)
		{
			unsigned int v = face_vertices[off + k];
			if (v >= vertex_count)
				return QAWS_STATUS_INVALID_ARGUMENT;
			for (i = 0; i < 3; i++)
			{
				lo[i] = fmin(lo[i], (double)positions[3 * (size_t)v + i]);
				hi[i] = fmax(hi[i], (double)positions[3 * (size_t)v + i]);
			}
		}
		nt += s - 2;
		off += s;
	}
	diag = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	c.eps = desc && desc->tolerance > 0 ? (double)desc->tolerance : 1e-9 * diag;
	if (!(c.eps > 0))
		c.eps = 1e-300;
	report->tolerance = (qaws_scalar)c.eps;
	c.inv = 1 / c.eps;
	c.r = report;
	c.ncell = 1024;
	c.cell = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * c.ncell);
	alias = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)vertex_count + 1));
	c.F = (sx_face*)qaws_internal_zalloc((size_t)nt + 1, sizeof(sx_face));
	box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (size_t)nt + 1));
	if (!c.cell || !alias || !c.F || !box)
		goto oom;
	for (i = 0; i < c.ncell; i++)
		c.cell[i] = SX_NONE;
	/* vertices: those within the tolerance of an earlier one are welded */
	for (i = 0; i < vertex_count; i++)
	{
		double p[3] = { (double)positions[3 * (size_t)i], (double)positions[3 * (size_t)i + 1], (double)positions[3 * (size_t)i + 2] };
		unsigned int before = c.np;
		alias[i] = sx_point(&c, p);
		if (c.oom)
			goto oom;
		report->welded_vertices += c.np == before;
	}
	onv = c.np;
	/* triangles (polygons fanned) with their planes and edge carriers */
	for (f = 0, off = 0; f < face_count; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3;
		for (k = 1; k + 1 < s; k++)
		{
			sx_face* F = c.F + c.nf;
			double u[3], v[3], l;
			F->c[0] = alias[face_vertices[off]];
			F->c[1] = alias[face_vertices[off + k]];
			F->c[2] = alias[face_vertices[off + k + 1]];
			F->src = f;
			F->cuts = SX_NONE;
			F->touch = SX_NONE;
			if (F->c[0] == F->c[1] || F->c[1] == F->c[2] || F->c[2] == F->c[0])
				continue;   /* welded away */
			for (i = 0; i < 3; i++)
			{
				u[i] = c.X[3 * (size_t)F->c[1] + i] - c.X[3 * (size_t)F->c[0] + i];
				v[i] = c.X[3 * (size_t)F->c[2] + i] - c.X[3 * (size_t)F->c[0] + i];
			}
			sx_cross(u, v, F->n);
			l = sqrt(sx_dot(F->n, F->n));
			/* lower than the tolerance over its longest edge: no plane, cut by nothing */
			if (l > c.eps * sqrt(fmax(fmax(sx_dot(u, u), sx_dot(v, v)), sx_dist2(c.X + 3 * (size_t)F->c[1], c.X + 3 * (size_t)F->c[2]))))
			{
				F->n[0] /= l;
				F->n[1] /= l;
				F->n[2] /= l;
			}
			else
				F->n[0] = F->n[1] = F->n[2] = 0;
			F->d = sx_dot(F->n, c.X + 3 * (size_t)F->c[0]);
			for (i = 0; i < 3; i++)
			{
				unsigned int p = F->c[i], q = F->c[(i + 1) % 3];
				F->e[i] = sx_edge(&c, p, q);
				if (c.oom)
					goto oom;
				sx_use(&c, F->e[i], c.nf);
			}
			for (i = 0; i < 3; i++)
			{
				box[c.nf].lo[i] = fmin(fmin(c.X[3 * (size_t)F->c[0] + i], c.X[3 * (size_t)F->c[1] + i]), c.X[3 * (size_t)F->c[2] + i]) - c.eps;
				box[c.nf].hi[i] = fmax(fmax(c.X[3 * (size_t)F->c[0] + i], c.X[3 * (size_t)F->c[1] + i]), c.X[3 * (size_t)F->c[2] + i]) + c.eps;
			}
			c.nf++;
		}
		off += s;
	}
	if (c.oom)
		goto oom;
	/* the pairs cut each other */
	if (c.nf)
		st = qaws_internal_broadphase(NULL, box, c.nf, 3, NULL, sx_bp_visit, &c);
	if (st != QAWS_STATUS_OK || c.oom)
		goto fail;
	/* faces settled until no carrier gains a point */
	{
		unsigned long long guard = 64ull * c.nf + 1024, head = 0;
		while (head < c.nw && guard--)
		{
			unsigned int g = c.work[head++];
			c.F[g].queued = 0;
			sx_settle(&c, &b, g);
			if (c.oom)
				goto oom;
			if (head > 4096 && 2 * head > c.nw)
			{
				memmove(c.work, c.work + head, sizeof(unsigned int) * (size_t)(c.nw - head));
				c.nw -= (unsigned int)head;
				head = 0;
			}
		}
		if (head < c.nw)
			report->unresolved += (unsigned int)(c.nw - head);
	}
	/* each face triangulated again through its pieces */
	loc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)c.np + 1));
	used = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)c.np + 1));
	if (!loc || !used)
		goto oom;
	for (i = 0; i < c.np; i++)
	{
		loc[i] = SX_NONE;
		used[i] = SX_NONE;
	}
	for (f = 0; f < c.nf; f++)
	{
		sx_face const* F = c.F + f;
		unsigned int np, need;
		if (F->c[0] == F->c[1] || F->c[1] == F->c[2] || F->c[2] == F->c[0])
			continue;
		sx_gather(&c, &b, f);
		if (c.oom)
			goto oom;
		np = b.npt;
		need = np > 3 ? 2 * np : 1;
		if (!sx_reserve((void**)&otri, &capt, 3 * ((size_t)ntri + need + 1), sizeof(unsigned int))
			|| !sx_reserve((void**)&osrc, &caps, (size_t)ntri + need + 1, sizeof(qaws_scalar)))
			goto oom;
		if (np == 3)
		{
			memcpy(otri + 3 * (size_t)ntri, F->c, 3 * sizeof(unsigned int));
			osrc[ntri++] = (qaws_scalar)F->src;
			continue;
		}
		report->split_faces++;
		if (!sx_reserve((void**)&T.P, &T.capp, 2 * (size_t)np, sizeof(double))
			|| !sx_reserve((void**)&T.Q, &T.capq, 2 * (size_t)np, sizeof(double))
			|| !sx_reserve((void**)&T.t, &T.capt, 2 * (size_t)np + 4, sizeof(sx_tri))
			|| !sx_reserve((void**)&stack, &caps2, 64 * (size_t)np + 64, sizeof(unsigned int)))
			goto oom;
		/* local coordinates: the corners at (0, 0), (1, 0), (0, 1); points on an edge exactly on it */
		{
			double e1[3], e2[3], g11, g12, g22, det;
			unsigned int z;
			for (z = 0; z < 3; z++)
			{
				e1[z] = c.X[3 * (size_t)F->c[1] + z] - c.X[3 * (size_t)F->c[0] + z];
				e2[z] = c.X[3 * (size_t)F->c[2] + z] - c.X[3 * (size_t)F->c[0] + z];
			}
			g11 = sx_dot(e1, e1);
			g12 = sx_dot(e1, e2);
			g22 = sx_dot(e2, e2);
			det = g11 * g22 - g12 * g12;
			for (i = 0; i < np; i++)
				loc[b.p[i]] = i;
			for (i = 0; i < np; i++)
			{
				double* q = T.P + 2 * (size_t)i;
				q[0] = -1;
			}
			T.P[0] = 0; T.P[1] = 0;
			T.P[2] = 1; T.P[3] = 0;
			T.P[4] = 0; T.P[5] = 1;
			for (k = 0; k < 3; k++)
			{
				unsigned int j;
				for (j = c.car[F->e[k]].points; j != SX_NONE; j = c.nodes[j].next)
				{
					unsigned int v = c.nodes[j].a, li = loc[v];
					double d2, t;
					if (li == SX_NONE || li < 3)
						continue;
					t = sx_snap52(sx_foot(c.X + 3 * (size_t)v, c.X + 3 * (size_t)F->c[k], c.X + 3 * (size_t)F->c[(k + 1) % 3], &d2));
					if (k == 0)
					{
						T.P[2 * (size_t)li] = t;
						T.P[2 * (size_t)li + 1] = 0;
					}
					else if (k == 1)
					{
						T.P[2 * (size_t)li] = 1 - t;
						T.P[2 * (size_t)li + 1] = t;
					}
					else
					{
						T.P[2 * (size_t)li] = 0;
						T.P[2 * (size_t)li + 1] = 1 - t;
					}
				}
			}
			for (i = 3; i < np; i++)
				if (T.P[2 * (size_t)i] < 0)
				{
					double w[3], r1, r2, uu, vv;
					for (z = 0; z < 3; z++)
						w[z] = c.X[3 * (size_t)b.p[i] + z] - c.X[3 * (size_t)F->c[0] + z];
					r1 = sx_dot(e1, w);
					r2 = sx_dot(e2, w);
					uu = det > 0 ? (g22 * r1 - g12 * r2) / det : 1.0 / 3;
					vv = det > 0 ? (g11 * r2 - g12 * r1) / det : 1.0 / 3;
					/* strictly inside: a point on an edge was registered on it */
					uu = fmax(uu, 0x1p-40);
					vv = fmax(vv, 0x1p-40);
					if (uu + vv > 1 - 0x1p-40)
					{
						double sc = (1 - 0x1p-40) / (uu + vv);
						uu *= sc;
						vv *= sc;
					}
					T.P[2 * (size_t)i] = uu;
					T.P[2 * (size_t)i + 1] = vv;
				}
			/* (u, v) to lengths: e1 along x */
			{
				double l1 = sqrt(g11), h2 = det > 0 ? sqrt(det) / l1 : 0;
				for (i = 0; i < np; i++)
				{
					double uu = T.P[2 * (size_t)i], vv = T.P[2 * (size_t)i + 1];
					T.Q[2 * (size_t)i] = uu * l1 + vv * g12 / l1;
					T.Q[2 * (size_t)i + 1] = vv * h2;
				}
			}
		}
		T.nt = 1;
		sx_set(T.t, 0, 1, 2, SX_NONE, SX_NONE, SX_NONE, 0, 0, 0);
		for (i = 3; i < np; i++)
			sx_insert(&T, i, stack, 16 * np + 64);
		for (i = 0; i < b.ns; i++)
		{
			unsigned int pa = loc[b.s[i].p], pb = loc[b.s[i].q];
			if (pa == pb || b.s[i].k == F->e[0] || b.s[i].k == F->e[1] || b.s[i].k == F->e[2])
				continue;   /* the face's own edges are its border */
			if (!sx_recover(&T, pa, pb, 0))
				report->unresolved++;
		}
		sx_relax(&T);
		for (i = 0; i < T.nt; i++)
		{
			otri[3 * (size_t)ntri] = b.p[T.t[i].v[0]];
			otri[3 * (size_t)ntri + 1] = b.p[T.t[i].v[1]];
			otri[3 * (size_t)ntri + 2] = b.p[T.t[i].v[2]];
			osrc[ntri++] = (qaws_scalar)F->src;
		}
		for (i = 0; i < np; i++)
			loc[b.p[i]] = SX_NONE;
	}
	/* the vertices used, compacted */
	{
		unsigned int nv = 0;
		for (i = 0; i < 3 * ntri; i++)
		{
			unsigned int v = otri[i];
			if (used[v] == SX_NONE)
			{
				used[v] = nv++;
				report->points += v >= onv;
			}
			otri[i] = used[v];
		}
		opos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
		if (!opos)
			goto oom;
		for (i = 0; i < c.np; i++)
			if (used[i] != SX_NONE)
				for (k = 0; k < 3; k++)
					opos[3 * (size_t)used[i] + k] = (qaws_scalar)c.X[3 * (size_t)i + k];
		/* checked by exact predicates */
		{
			qaws_mesh_repair_desc rd;
			qaws_mesh_repair_report rr;
			qaws_mesh_3d* chk = NULL;
			memset(&rd, 0, sizeof rd);
			rd.weld_tolerance = -1;
			rd.degenerate_tolerance = -1;
			rd.gap_tolerance = -1;
			rd.hole_perimeter = -1;
			rd.flags = QAWS_MESH_REPAIR_FIND_INTERSECTIONS;
			st = qaws_mesh_repair(&chk, &rr, opos, otri, NULL, nv, ntri, &rd);
			qaws_mesh_3d_destroy(chk);
			if (st != QAWS_STATUS_OK)
				goto fail;
			report->remaining = rr.self_intersections;
		}
		{
			qaws_mesh_3d_desc md;
			qaws_attribute_desc ad;
			unsigned int attr;
			memset(&md, 0, sizeof md);
			md.positions = opos;
			md.vertex_count = nv;
			md.face_vertices = otri;
			md.face_count = ntri;
			st = qaws_mesh_3d_create(out, &md);
			if (st != QAWS_STATUS_OK)
				goto fail;
			memset(&ad, 0, sizeof ad);
			ad.name = "source";
			ad.element = QAWS_MESH_FACE;
			ad.components = 1;
			ad.rule = QAWS_ATTRIBUTE_NEAREST;
			st = qaws_mesh_3d_add_attribute(&attr, *out, osrc, &ad);
			if (st != QAWS_STATUS_OK)
			{
				qaws_mesh_3d_destroy(*out);
				*out = NULL;
			}
		}
	}
	goto fail;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
fail:
	qaws_internal_free(c.X); qaws_internal_free(c.chain); qaws_internal_free(c.cell); qaws_internal_free(c.keys); qaws_internal_free(c.nodes);
	qaws_internal_free(c.car); qaws_internal_free(c.F); qaws_internal_free(c.work); qaws_internal_free(b.s); qaws_internal_free(b.p);
	qaws_internal_free(T.P); qaws_internal_free(T.Q); qaws_internal_free(T.t); qaws_internal_free(alias); qaws_internal_free(loc); qaws_internal_free(used);
	qaws_internal_free(otri); qaws_internal_free(stack); qaws_internal_free(opos); qaws_internal_free(osrc);
	qaws_internal_free(box);
	return st;
}
