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
#include "../qaws_mesh_boolean.h"
#include "../qaws_mesh_bvh.h"
#include "../qaws_mesh_repair.h"
#include "../qaws_exact.h"
#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_parallel.h"
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
	unsigned int e[3];                 /* carriers of the edges c0 c1, c1 c2, c2 c0; SX_NONE until the face meets something */
	unsigned int src;
	unsigned int cuts;                 /* list: node.a carrier, node.b and node.c its ends in the face */
	unsigned int touch;                /* isolated points: node.a */
	double n[3], d;                    /* unit normal and offset; n zero: degenerate */
	int queued;
	int touched;                       /* triangulated through a carrier (its own edges' carriers made on demand) */
} sx_face;

/* a face as the pair filter reads it, in one place (128 bytes) */
typedef struct sx_geo
{
	double x[3][3];                    /* corners */
	double n[3], d;                    /* the face's plane (n zero: degenerate) */
	unsigned int c[3], pad;
} sx_geo;

typedef struct sx_key
{
	unsigned int k[5];
	unsigned int id;
} sx_key;

typedef struct sx_ctx
{
	double* X;                         /* points */
	unsigned int* chain;               /* next point in the same spatial cell */
	unsigned int* hash;                /* each point's cell hash */
	unsigned int np, capx, capc, caph;
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
	unsigned int *vf_start, *vf;       /* faces round each point */
	sx_geo* geo;
	unsigned int* work;
	unsigned int nw, capw;
	double eps, inv;
	qaws_intersect_report* r;
	unsigned char* copl;               /* per input face: in a coplanar overlap */
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

static double sx_max3(double a, double b, double c)
{
	double m = a > b ? a : b;
	return m > c ? m : c;
}

static unsigned int sx_hash(unsigned int const* k, unsigned int n)
{
	unsigned int h = 2166136261u, i;
	for (i = 0; i < n; i++)
		h = (h ^ k[i]) * 16777619u;
	return h ^ (h >> 15);
}

/* the hash of a point's cell (masked to the table at use) */
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
	return sx_hash(k, 6);
}

static void sx_rehash(sx_ctx* c)
{
	unsigned int i;
	for (i = 0; i < c->ncell; i++)
		c->cell[i] = SX_NONE;
	for (i = 0; i < c->np; i++)
	{
		unsigned int h = c->hash[i] & (c->ncell - 1);
		c->chain[i] = c->cell[h];
		c->cell[h] = i;
	}
}

/* the hash of p's cell, and the sides of the cell p is within the tolerance of (bit 2 k: low on axis k, 2 k + 1: high) */
static unsigned int sx_cell_sides(sx_ctx const* c, double const* p, unsigned int* sides)
{
	unsigned int k;
	*sides = 0;
	for (k = 0; k < 3; k++)
	{
		double f = p[k] * c->inv - floor(p[k] * c->inv);
		*sides |= (unsigned int)(f < 1.0 / 64) << 2 * k | (unsigned int)(f > 1 - 1.0 / 64) << (2 * k + 1);
	}
	return sx_cell_of(c, p, 0, 0, 0);
}

/* the first point within the tolerance of p (own its cell's hash, sides as sx_cell_sides), else a new one */
static unsigned int sx_point_at(sx_ctx* c, double const* p, unsigned int own, unsigned int sides)
{
	int dx, dy, dz, lo[3], hi[3];
	unsigned int id, h, k;
	double e2 = c->eps * c->eps;
	/* the neighbouring cells only on the sides the point is within the tolerance of */
	for (k = 0; k < 3; k++)
	{
		lo[k] = (sides >> 2 * k) & 1 ? -1 : 0;
		hi[k] = (sides >> (2 * k + 1)) & 1 ? 1 : 0;
	}
	for (dx = lo[0]; dx <= hi[0]; dx++)
		for (dy = lo[1]; dy <= hi[1]; dy++)
			for (dz = lo[2]; dz <= hi[2]; dz++)
				for (id = c->cell[((dx | dy | dz) ? sx_cell_of(c, p, dx, dy, dz) : own) & (c->ncell - 1)]; id != SX_NONE; id = c->chain[id])
					if (sx_dist2(c->X + 3 * (size_t)id, p) <= e2)
						return id;
	if (!sx_reserve((void**)&c->X, &c->capx, 3 * ((size_t)c->np + 1), sizeof(double))
		|| !sx_reserve((void**)&c->chain, &c->capc, (size_t)c->np + 1, sizeof(unsigned int))
		|| !sx_reserve((void**)&c->hash, &c->caph, (size_t)c->np + 1, sizeof(unsigned int)))
	{
		c->oom = 1;
		return 0;
	}
	memcpy(c->X + 3 * (size_t)c->np, p, 3 * sizeof(double));
	c->hash[c->np] = own;   /* (kept for the rehash as the table grows) */
	h = own & (c->ncell - 1);
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

static unsigned int sx_point(sx_ctx* c, double const* p)
{
	unsigned int sides, own = sx_cell_sides(c, p, &sides);
	return sx_point_at(c, p, own, sides);
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
	c->F[f].touched = 1;
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

/* face f's edge carriers, made when it first meets something; every face on each edge is triangulated through it */
static void sx_face_edges(sx_ctx* c, unsigned int f)
{
	unsigned int i, j;
	if (c->F[f].e[0] != SX_NONE)
		return;
	for (i = 0; i < 3 && !c->oom; i++)
	{
		unsigned int p = c->F[f].c[i], q = c->F[f].c[(i + 1) % 3], e = sx_edge(c, p, q);
		if (c->oom)
			return;
		c->F[f].e[i] = e;
		for (j = c->vf_start[p]; j < c->vf_start[p + 1]; j++)
		{
			sx_face const* G = c->F + c->vf[j];
			if (G->c[0] == q || G->c[1] == q || G->c[2] == q)
				sx_use(c, e, c->vf[j]);
		}
	}
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

/* edge a b through the plane of face f (p kept when the edge is parallel to it) */
static void sx_lpi_pos(sx_ctx const* c, unsigned int ia, unsigned int ib, unsigned int f, double* p)
{
	double const *a = c->X + 3 * (size_t)ia, *b = c->X + 3 * (size_t)ib;
	sx_face const* F = c->F + f;
	double da = sx_dot(F->n, a) - F->d, db = sx_dot(F->n, b) - F->d;
	if (fabs(da - db) > 0)
	{
		double t = fmin(fmax(da / (da - db), 0.0), 1.0);
		p[0] = a[0] + t * (b[0] - a[0]);
		p[1] = a[1] + t * (b[1] - a[1]);
		p[2] = a[2] + t * (b[2] - a[2]);
	}
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
		sx_lpi_pos(c, k[1], k[2], k[3], p);
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

/* a point where a face's boundary meets another's plane, not made yet: a corner, or an edge through the plane */
typedef struct sx_cand
{
	unsigned int vertex, edge, face;   /* vertex SX_NONE: edge (a carrier) through face's plane */
	double x[3];
	double s;                          /* along the two planes' line */
} sx_cand;

/* where face g's boundary meets the plane of face f (signs of g's corners s): up to two candidates; nothing is
   made until the pair is known to meet (an edge through the plane beside the other face is no point) */
static unsigned int sx_section(sx_ctx const* c, unsigned int g, unsigned int f, int const* s, sx_cand* out)
{
	sx_face const* G = c->F + g;
	unsigned int n = 0, k;
	for (k = 0; k < 3; k++)
		if (!s[k] && n < 2)
		{
			out[n].vertex = G->c[k];
			out[n].edge = SX_NONE;
			memcpy(out[n].x, c->X + 3 * (size_t)G->c[k], sizeof out[n].x);
			n++;
		}
	for (k = 0; k < 3; k++)
		if (s[k] * s[(k + 1) % 3] < 0 && n < 2)
		{
			sx_carrier const* E = c->car + G->e[k];
			out[n].vertex = SX_NONE;
			out[n].edge = G->e[k];
			out[n].face = f;
			sx_lpi_pos(c, E->a, E->b, f, out[n].x);
			n++;
		}
	return n;
}

static unsigned int sx_make(sx_ctx* c, sx_cand const* a)
{
	return a->vertex != SX_NONE ? a->vertex : sx_lpi(c, a->edge, a->face);
}

/* coplanar faces f and g kept apart by a line through an edge of either, within the tolerance (their interiors overlap
   by less than it: neighbours touching at a corner or along an edge are not cut) */
static int sx_separated(sx_ctx const* c, unsigned int f, unsigned int g)
{
	sx_face const* T[2] = { c->F + f, c->F + g };
	double const* n = T[0]->n;
	unsigned int ax = fabs(n[0]) > fabs(n[1]) ? (fabs(n[0]) > fabs(n[2]) ? 0 : 2) : (fabs(n[1]) > fabs(n[2]) ? 1 : 2);
	unsigned int x = ax == 0 ? 1 : 0, y = ax == 2 ? 1 : 2, s, j, k;
	for (s = 0; s < 2; s++)
	{
		sx_face const *A = T[s], *B = T[1 - s];
		double a2[3][2], o;
		for (k = 0; k < 3; k++)
		{
			a2[k][0] = c->X[3 * (size_t)A->c[k] + x];
			a2[k][1] = c->X[3 * (size_t)A->c[k] + y];
		}
		o = (a2[1][0] - a2[0][0]) * (a2[2][1] - a2[0][1]) - (a2[1][1] - a2[0][1]) * (a2[2][0] - a2[0][0]) > 0 ? 1 : -1;
		for (j = 0; j < 3; j++)
		{
			double const *u = a2[j], *v = a2[(j + 1) % 3];
			double ex = v[0] - u[0], ey = v[1] - u[1], el = sqrt(ex * ex + ey * ey);
			int out = el > 0;
			for (k = 0; k < 3 && out; k++)
			{
				double const* q = c->X + 3 * (size_t)B->c[k];
				out = o * (ex * (q[y] - u[1]) - ey * (q[x] - u[0])) / el <= c->eps;
			}
			if (out)
				return 1;
		}
	}
	return 0;
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
	sx_face_edges(c, f);
	sx_face_edges(c, g);
	if (c->oom)
		return QAWS_STATUS_ALLOCATION_FAILURE;
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
		if (sx_separated(c, f, g))
			return QAWS_STATUS_OK;
		hit = sx_coplanar(c, f, g);
		hit |= sx_coplanar(c, g, f);
		if (hit)
		{
			c->r->pairs++;
			c->r->coplanar_pairs++;
			c->copl[c->F[f].src] = c->copl[c->F[g].src] = 1;
		}
		return c->oom ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	}
	if (shared == 2)
		return QAWS_STATUS_OK;   /* two planes through a shared edge meet only along it */
	{
		sx_cand ia[2], ib[2], t;
		unsigned int n0, n1, P, Q;
		double dir[3], l, sP, sQ;
		n0 = sx_section(c, f, g, sa, ia);
		n1 = sx_section(c, g, f, sb, ib);
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
		for (i = 0; i < 2; i++)
		{
			ia[i].s = sx_dot(ia[i].x, dir) / l;
			ib[i].s = sx_dot(ib[i].x, dir) / l;
		}
		if (ia[0].s > ia[1].s)
		{
			t = ia[0];
			ia[0] = ia[1];
			ia[1] = t;
		}
		if (ib[0].s > ib[1].s)
		{
			t = ib[0];
			ib[0] = ib[1];
			ib[1] = t;
		}
		sP = fmax(ia[0].s, ib[0].s);
		sQ = fmin(ia[1].s, ib[1].s);
		if (sP > sQ + c->eps)
			return QAWS_STATUS_OK;
		/* they meet: only now the two ends are made */
		P = sx_make(c, ia[0].s >= ib[0].s ? ia : ib);
		Q = sx_make(c, ia[1].s <= ib[1].s ? ia + 1 : ib + 1);
		if (c->oom)
			return QAWS_STATUS_ALLOCATION_FAILURE;
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

typedef struct sx_span
{
	double lo, hi;
	unsigned int i;
} sx_span;

typedef struct sx_buf
{
	sx_piece* s;
	unsigned int ns, caps;
	unsigned int* p;
	unsigned int npt, capp;
	sx_span* sp;                       /* the pieces' then the points' spans along the face's longest axis */
	unsigned int capsp;
} sx_buf;

static int sx_span_cmp(void const* p, void const* q)
{
	sx_span const *a = (sx_span const*)p, *b = (sx_span const*)q;
	return a->lo < b->lo ? -1 : a->lo > b->lo ? 1 : a->i < b->i ? -1 : a->i > b->i;
}

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

/* the carriers of face f's edges that exist (SX_NONE: an edge nothing reached, its end points only) */
static void sx_edges_found(sx_ctx* c, unsigned int f, unsigned int* e)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		unsigned int p = c->F[f].c[k], q = c->F[f].c[(k + 1) % 3];
		unsigned int key[5] = { SX_EDGE, p < q ? p : q, p < q ? q : p, 0, 0 };
		sx_key const* s = c->F[f].e[k] != SX_NONE ? NULL : sx_find(c, key);
		e[k] = c->F[f].e[k] != SX_NONE ? c->F[f].e[k] : s ? s->id : SX_NONE;
	}
}

static void sx_gather(sx_ctx* c, sx_buf* b, unsigned int f, unsigned int const* e)
{
	sx_face const* F = c->F + f;
	unsigned int k, i;
	b->ns = 0;
	b->npt = 0;
	for (k = 0; k < 3; k++)
		sx_add_point(c, b, F->c[k]);
	for (k = 0; k < 3; k++)
		if (e[k] != SX_NONE)
			sx_expand(c, b, e[k], F->c[k], F->c[(k + 1) % 3]);
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

/* pieces s and t of face f crossing (within eps, e2 its square): cut at a keyed point on both; whether a carrier gained it */
static int sx_settle_cross(sx_ctx* c, unsigned int f, sx_piece const* s, sx_piece const* t, double e2)
{
	sx_carrier const *K = c->car + s->k, *L = c->car + t->k;
	double const *a = c->X + 3 * (size_t)s->p, *a1 = c->X + 3 * (size_t)s->q, *q = c->X + 3 * (size_t)t->p, *q1 = c->X + 3 * (size_t)t->q;
	double u[3], v[3], w[3], uu, uv, vv, uw, vw, den, su, tv, m[3], n[3], x[3];
	unsigned int key[5], id, h;
	if (s->k == t->k || s->p == t->p || s->p == t->q || s->q == t->p || s->q == t->q)
		return 0;
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
		return 0;   /* parallel: overlaps show as points on pieces */
	su = (uv * vw - vv * uw) / den;
	tv = (uu * vw - uv * uw) / den;
	if (!(su > 0 && su < 1 && tv > 0 && tv < 1))
		return 0;
	for (h = 0; h < 3; h++)
	{
		m[h] = a[h] + su * u[h];
		n[h] = q[h] + tv * v[h];
		x[h] = 0.5 * (m[h] + n[h]);
	}
	if (sx_dist2(m, n) > e2 || sx_dist2(x, a) <= e2 || sx_dist2(x, a1) <= e2 || sx_dist2(x, q) <= e2 || sx_dist2(x, q1) <= e2)
		return 0;
	memset(key, 0, sizeof key);
	if (K->kind == SX_EDGE && L->kind == SX_EDGE)
	{
		sx_ee(c, s->k, t->k, x);
		return 1;
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
			return 0;
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
	return sx_on(c, s->k, id) | sx_on(c, t->k, id);
}

/* one pass over face f; returns whether a carrier gained a point. Points and pieces are met along the face's longest axis
   (spans sorted, only those overlapping tested): a face under many overlapping ones holds thousands of pieces */
static int sx_settle(sx_ctx* c, sx_buf* b, unsigned int f)
{
	unsigned int i, j, k, ax = 0;
	double e2 = c->eps * c->eps, ext = -1;
	int changed = 0;
	sx_span *ps, *pt;
	sx_face_edges(c, f);
	sx_gather(c, b, f, c->F[f].e);
	if (c->oom)
		return 0;
	for (k = 0; k < 3; k++)
	{
		double lo = HUGE_VAL, hi = -HUGE_VAL;
		for (i = 0; i < 3; i++)
		{
			lo = fmin(lo, c->X[3 * (size_t)c->F[f].c[i] + k]);
			hi = fmax(hi, c->X[3 * (size_t)c->F[f].c[i] + k]);
		}
		if (hi - lo > ext)
		{
			ext = hi - lo;
			ax = k;
		}
	}
	if (!sx_reserve((void**)&b->sp, &b->capsp, (size_t)b->ns + b->npt + 1, sizeof(sx_span)))
	{
		c->oom = 1;
		return 0;
	}
	ps = b->sp;
	pt = b->sp + b->ns;
	for (j = 0; j < b->ns; j++)
	{
		double a = c->X[3 * (size_t)b->s[j].p + ax], q = c->X[3 * (size_t)b->s[j].q + ax];
		ps[j].lo = fmin(a, q) - c->eps;
		ps[j].hi = fmax(a, q) + c->eps;
		ps[j].i = j;
	}
	for (i = 0; i < b->npt; i++)
	{
		pt[i].lo = pt[i].hi = c->X[3 * (size_t)b->p[i] + ax];
		pt[i].i = i;
	}
	qsort(ps, b->ns, sizeof *ps, sx_span_cmp);
	qsort(pt, b->npt, sizeof *pt, sx_span_cmp);
	/* points on pieces */
	for (j = 0; j < b->ns; j++)
	{
		sx_piece const* s = b->s + ps[j].i;
		double const *a = c->X + 3 * (size_t)s->p, *q = c->X + 3 * (size_t)s->q;
		unsigned int lo = 0, hi = b->npt;
		while (lo < hi)
		{
			unsigned int mid = (lo + hi) / 2;
			if (pt[mid].lo < ps[j].lo)
				lo = mid + 1;
			else
				hi = mid;
		}
		for (i = lo; i < b->npt && pt[i].lo <= ps[j].hi; i++)
		{
			unsigned int v = b->p[pt[i].i];
			double const* x = c->X + 3 * (size_t)v;
			double d2, t;
			if (v == s->p || v == s->q)
				continue;
			t = sx_foot(x, a, q, &d2);
			if (t > 0 && t < 1 && d2 <= e2 && sx_dist2(x, a) > e2 && sx_dist2(x, q) > e2)
				changed |= sx_on(c, s->k, v);
		}
	}
	/* pieces crossing */
	for (i = 0; i < b->ns; i++)
		for (j = i + 1; j < b->ns && ps[j].lo <= ps[i].hi; j++)
		{
			unsigned int u = ps[i].i < ps[j].i ? ps[i].i : ps[j].i, w = ps[i].i < ps[j].i ? ps[j].i : ps[i].i;
			changed |= sx_settle_cross(c, f, b->s + u, b->s + w, e2);
			if (c->oom)
				return 0;
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
	unsigned int* vt;                  /* per local point a triangle holding it (kept through flips once vton) */
	unsigned int* st;                  /* scratch: a point's triangles */
	unsigned int* qe;                  /* scratch: the edges a cut crosses, 2 points each */
	unsigned int capvt, capst, capqe;
	int vton;
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
	if (T->vton)
	{
		T->vt[a] = u;
		T->vt[d] = u;
		T->vt[cc] = u;
		T->vt[b] = w;
	}
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

static int sx_recover2(sx_cdt* T, unsigned int a, unsigned int x, unsigned int b, unsigned int depth);

/* the triangles round local point a (T->vt kept): into T->st, their count */
static unsigned int sx_star(sx_cdt* T, unsigned int a)
{
	unsigned int u0 = T->vt[a], u, n = 0, i, dir;
	for (dir = 0; dir < 2; dir++)
	{
		u = u0;
		while (u != SX_NONE && n < T->nt)
		{
			for (i = 0; i < 3 && T->t[u].v[i] != a; i++)
				;
			if (i == 3)
				break;
			if (dir == 0 || u != u0)
				T->st[n++] = u;
			u = T->t[u].n[dir ? i : (i + 2) % 3];   /* (round a one way, then the other from where it started) */
			if (u == u0)
				return n;
		}
	}
	return n;
}

/* the triangle and edge index of edge a b (SX_NONE when it is not there) */
static unsigned int sx_find_edge(sx_cdt* T, unsigned int a, unsigned int b, unsigned int* k)
{
	unsigned int n = sx_star(T, a), i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < 3; j++)
		{
			unsigned int p = T->t[T->st[i]].v[j], q = T->t[T->st[i]].v[(j + 1) % 3];
			if ((p == a && q == b) || (p == b && q == a))
			{
				*k = j;
				return T->st[i];
			}
		}
	return SX_NONE;
}

/* edge a b there: marked a cut (both sides) */
static int sx_has_edge(sx_cdt* T, unsigned int a, unsigned int b)
{
	unsigned int k = 0, u = sx_find_edge(T, a, b, &k), w, j;
	if (u == SX_NONE)
		return 0;
	T->t[u].fix[k] = 1;
	w = T->t[u].n[k];
	if (w != SX_NONE)
		for (j = 0; j < 3; j++)
			if (T->t[w].n[j] == u)
				T->t[w].fix[j] = 1;
	return 1;
}

/* local point x strictly between a and b (on their line) */
static int sx_between(sx_cdt const* T, unsigned int a, unsigned int b, unsigned int x)
{
	double const *A = T->P + 2 * (size_t)a, *B = T->P + 2 * (size_t)b, *X = T->P + 2 * (size_t)x;
	double dt = (X[0] - A[0]) * (B[0] - A[0]) + (X[1] - A[1]) * (B[1] - A[1]);
	double ll = (B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]);
	return dt > 0 && dt < ll;
}

static int sx_push_edge(sx_cdt* T, unsigned int* nq, unsigned int p, unsigned int q)
{
	if (!sx_reserve((void**)&T->qe, &T->capqe, 2 * (size_t)*nq + 2, sizeof(unsigned int)))
		return 0;
	T->qe[2 * (size_t)*nq] = p;
	T->qe[2 * (size_t)*nq + 1] = q;
	(*nq)++;
	return 1;
}

/* the cut a b recovered by flips (Sloan): the edges it crosses found by walking from a to b, flipped in turn (a flip
   whose new edge still crosses goes back in line, a quad not convex waits); split at a point exactly on it. 0 when it
   cannot be (it crosses another cut), -1 out of memory */
static int sx_recover(sx_cdt* T, unsigned int a, unsigned int b, unsigned int depth)
{
	unsigned int n, i, u = SX_NONE, p = SX_NONE, q = SX_NONE, nq = 0, head = 0, guard;
	if (sx_has_edge(T, a, b))
		return 1;
	/* the triangle of a's star the cut leaves by */
	n = sx_star(T, a);
	for (i = 0; i < n && u == SX_NONE; i++)
	{
		sx_tri const* t = T->t + T->st[i];
		unsigned int ia, v1, v2, z;
		int s1, s2;
		for (ia = 0; ia < 3 && t->v[ia] != a; ia++)
			;
		v1 = t->v[(ia + 1) % 3];
		v2 = t->v[(ia + 2) % 3];
		s1 = sx_o2(T->P, a, b, v1);
		s2 = sx_o2(T->P, a, b, v2);
		for (z = 0; z < 2; z++)
		{
			unsigned int x = z ? v2 : v1;
			if ((z ? s2 : s1) == 0 && sx_between(T, a, b, x))
				return depth < 8 ? sx_recover2(T, a, x, b, depth) : 0;
		}
		if (s1 * s2 < 0 && sx_o2(T->P, v1, v2, a) * sx_o2(T->P, v1, v2, b) < 0)
		{
			u = T->st[i];
			p = v1;
			q = v2;
		}
	}
	if (u == SX_NONE)
		return 0;
	/* the walk: each crossed edge in line, to the triangle holding b */
	for (guard = T->nt + 4; guard--;)
	{
		unsigned int k, w, j, d;
		int sd;
		for (k = 0; k < 3; k++)
		{
			unsigned int x = T->t[u].v[k], y = T->t[u].v[(k + 1) % 3];
			if ((x == p && y == q) || (x == q && y == p))
				break;
		}
		if (k == 3 || T->t[u].fix[k] || (w = T->t[u].n[k]) == SX_NONE)
			return 0;   /* (another cut crossed, or the border) */
		if (!sx_push_edge(T, &nq, p, q))
			return -1;
		for (j = 0; j < 3 && T->t[w].n[j] != u; j++)
			;
		d = T->t[w].v[(j + 2) % 3];
		if (d == b)
			break;
		sd = sx_o2(T->P, a, b, d);
		if (sd == 0)
			return sx_between(T, a, b, d) && depth < 8 ? sx_recover2(T, a, d, b, depth) : 0;
		if (sd == sx_o2(T->P, a, b, p))
			p = d;
		else
			q = d;
		u = w;
	}
	/* the flips */
	for (guard = 8 * T->nt + 64; head < nq && guard--;)
	{
		unsigned int x = T->qe[2 * (size_t)head], y = T->qe[2 * (size_t)head + 1], k = 0, c, d;
		head++;
		u = sx_find_edge(T, x, y, &k);
		if (u == SX_NONE)
			continue;
		if (T->t[u].fix[k] || T->t[u].n[k] == SX_NONE)
			return 0;
		c = T->t[u].v[(k + 2) % 3];
		d = sx_opp(T, u, k);
		if (sx_o2(T->P, c, d, x) * sx_o2(T->P, c, d, y) < 0)
		{
			sx_flip(T, u, k);
			if (sx_o2(T->P, a, b, c) * sx_o2(T->P, a, b, d) < 0 && sx_o2(T->P, c, d, a) * sx_o2(T->P, c, d, b) < 0 && !sx_push_edge(T, &nq, c, d))
				return -1;
		}
		else if (!sx_push_edge(T, &nq, x, y))
			return -1;
		if (head > 4096 && 2 * head > nq)
		{
			memmove(T->qe, T->qe + 2 * (size_t)head, sizeof(unsigned int) * 2 * (size_t)(nq - head));
			nq -= head;
			head = 0;
		}
	}
	return sx_has_edge(T, a, b);
}

/* the cut a b split at x on it: both halves (-1 out of memory) */
static int sx_recover2(sx_cdt* T, unsigned int a, unsigned int x, unsigned int b, unsigned int depth)
{
	int r = sx_recover(T, a, x, depth + 1), s;
	if (r < 0)
		return r;
	s = sx_recover(T, x, b, depth + 1);
	return s < 0 ? s : r && s;
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

typedef struct sx_chk
{
	qaws_scalar const* x;
	unsigned int const *t, *cut;
	qaws_scalar const* src;
	unsigned int* found;               /* per chunk */
	double m;                          /* the margin a floating point separation must pass */
} sx_chk;

/* certainly apart (floating point, by more than m): a plane of one triangle with the other's vertices it does not share
   strictly on one side, or a plane through an edge of one, square to it, with the other's unshared vertices strictly
   beyond (its shared ones on that edge); 0 leaves the pair to the exact test. Coplanar stacks of pieces (overlapping
   sweeps) pair millions of neighbours that only touch */
static int sx_chk_apart(double const (*P)[3], unsigned int const* ip, double const (*Q)[3], unsigned int const* iq, double m)
{
	unsigned int r, k, e;
	for (r = 0; r < 2; r++)
	{
		double const (*A)[3] = r ? Q : P, (*B)[3] = r ? P : Q;
		unsigned int const *ia = r ? iq : ip, *ib = r ? ip : iq;
		double u[3], v[3], n[3], l;
		int side = 0, ok = 1, any = 0;
		for (k = 0; k < 3; k++)
		{
			u[k] = A[1][k] - A[0][k];
			v[k] = A[2][k] - A[0][k];
		}
		sx_cross(u, v, n);
		l = sqrt(sx_dot(n, n));
		if (!(l > 0))
			continue;
		for (k = 0; k < 3 && ok; k++)
		{
			double w[3], d;
			unsigned int z;
			if (ib[k] == ia[0] || ib[k] == ia[1] || ib[k] == ia[2])
				continue;
			for (z = 0; z < 3; z++)
				w[z] = B[k][z] - A[0][z];
			d = sx_dot(n, w) / l;
			if (fabs(d) <= m || (side && (d > 0 ? 1 : -1) != side))
				ok = 0;
			side = d > 0 ? 1 : -1;
			any = 1;
		}
		if (ok && any)
			return 1;
		for (e = 0; e < 3; e++)
		{
			unsigned int a0 = ia[e], a1 = ia[(e + 1) % 3], a2 = ia[(e + 2) % 3], cnt = 0, z;
			double w[3], c[3], s0, lw, q[3];
			int good = 1;
			for (z = 0; z < 3; z++)
			{
				w[z] = A[(e + 1) % 3][z] - A[e][z];
				q[z] = A[(e + 2) % 3][z] - A[e][z];
			}
			sx_cross(n, w, c);
			lw = sqrt(sx_dot(c, c));
			if (!(lw > 0))
				continue;
			s0 = sx_dot(c, q) / lw;
			for (k = 0; k < 3 && good; k++)
			{
				double d;
				if (ib[k] == a0 || ib[k] == a1)
					continue;
				if (ib[k] == a2)
				{
					good = 0;
					break;
				}
				for (z = 0; z < 3; z++)
					q[z] = B[k][z] - A[e][z];
				d = sx_dot(c, q) / lw;
				good = fabs(d) > m && d * s0 < 0;
				cnt++;
			}
			if (good && cnt)
				return 1;
		}
	}
	return 0;
}

static qaws_status sx_chk_visit(void* user, unsigned int chunk, unsigned int i, unsigned int j)
{
	sx_chk const* C = (sx_chk const*)user;
	unsigned int const *ip = C->t + 3 * (size_t)C->cut[i], *iq = C->t + 3 * (size_t)C->cut[j];
	double P[3][3], Q[3][3];
	unsigned int k, l;
	if (C->src[C->cut[i]] == C->src[C->cut[j]])
		return QAWS_STATUS_OK;
	for (k = 0; k < 3; k++)
		for (l = 0; l < 3; l++)
		{
			P[k][l] = (double)C->x[3 * (size_t)ip[k] + l];
			Q[k][l] = (double)C->x[3 * (size_t)iq[k] + l];
		}
	if (sx_chk_apart((double const(*)[3])P, ip, (double const(*)[3])Q, iq, C->m))
		return QAWS_STATUS_OK;
	C->found[chunk] += (unsigned int)rp_tri_tri_at((double const(*)[3])P, (double const(*)[3])Q, ip, iq);
	return QAWS_STATUS_OK;
}

/* pairs of pieces (triangles cut[0 .. n)) crossing, by exact predicates; pieces of one face skipped */
static qaws_status sx_check(unsigned int* out, qaws_scalar const* x, unsigned int const* t, qaws_scalar const* src, unsigned int const* cut,
	unsigned int n, double eps, qaws_batch_executor const* executor)
{
	qaws_bp_box* box = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (size_t)n + 1));
	qaws_bp_grid* grid = NULL;
	sx_chk C;
	unsigned int i, k, nch;
	qaws_status st;
	*out = 0;
	if (!box)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n; i++)
		for (k = 0; k < 3; k++)
		{
			double a = (double)x[3 * (size_t)t[3 * (size_t)cut[i]] + k], b = (double)x[3 * (size_t)t[3 * (size_t)cut[i] + 1] + k],
				c = (double)x[3 * (size_t)t[3 * (size_t)cut[i] + 2] + k];
			box[i].lo[k] = fmin(fmin(a, b), c) - eps;
			box[i].hi[k] = fmax(fmax(a, b), c) + eps;
		}
	st = n > 1 ? qaws_internal_grid_create_fit(&grid, box, n, 3, executor) : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK && grid)
	{
		nch = qaws_internal_grid_pair_chunks(grid);
		C.x = x;
		C.t = t;
		C.cut = cut;
		C.src = src;
		C.m = 1e-3 * eps;
		C.found = (unsigned int*)qaws_internal_zalloc((size_t)nch + 1, sizeof(unsigned int));
		if (!C.found)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		else
		{
			st = qaws_internal_grid_pairs(NULL, grid, box, NULL, NULL, sx_chk_visit, &C, executor);
			for (k = 0; k < nch; k++)
				*out += C.found[k];
		}
		qaws_internal_free(C.found);
	}
	qaws_internal_grid_destroy(grid);
	qaws_internal_free(box);
	return st;
}

static int sx_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* ---- the pairs, filtered in parallel (read only), then cut in order ---- */

/* face g's boundary through face f's plane: up to two points, positions only */
static unsigned int sx_section_at(sx_ctx const* c, unsigned int g, unsigned int f, int const* s, double (*out)[3])
{
	sx_face const* G = c->F + g;
	unsigned int n = 0, k;
	for (k = 0; k < 3; k++)
		if (!s[k] && n < 2)
			memcpy(out[n++], c->X + 3 * (size_t)G->c[k], 3 * sizeof(double));
	for (k = 0; k < 3; k++)
		if (s[k] * s[(k + 1) % 3] < 0 && n < 2)
		{
			unsigned int p = G->c[k], q = G->c[(k + 1) % 3];
			memcpy(out[n], c->X + 3 * (size_t)p, 3 * sizeof(double));
			sx_lpi_pos(c, p < q ? p : q, p < q ? q : p, f, out[n]);
			n++;
		}
	return n;
}

/* faces f and g may meet beyond their shared corners and edges (the tests of sx_pair, nothing made) */
/* the corners of b not shared with a all strictly on one side of a's plane */
static int sx_geo_aside(sx_geo const* a, sx_geo const* b, double eps)
{
	int pos = 0, neg = 0, k;
	for (k = 0; k < 3; k++)
	{
		double d;
		if (b->c[k] == a->c[0] || b->c[k] == a->c[1] || b->c[k] == a->c[2])
			continue;
		d = sx_dot(a->n, b->x[k]) - a->d;
		if (d > eps)
			pos = 1;
		else if (d < -eps)
			neg = 1;
		else
			return 0;
	}
	return pos != neg;
}

static int sx_may_meet(sx_ctx const* c, unsigned int f, unsigned int g)
{
	sx_face const *A = c->F + f, *B = c->F + g;
	int sa[3], sb[3], za = 0, zb = 0, pa = 0, pb = 0, na = 0, nb = 0;
	unsigned int shared = 0, i, j, n0, n1;
	double xa[2][3], xb[2][3], dir[3], l, s[4];
	{
		/* quick: from the compact records, either face's free corners on one side of the other's plane (with a shared
		   corner, they then meet at most there) */
		sx_geo const *ga = c->geo + f, *gb = c->geo + g;
		if ((ga->n[0] == 0 && ga->n[1] == 0 && ga->n[2] == 0) || (gb->n[0] == 0 && gb->n[1] == 0 && gb->n[2] == 0))
			return 0;
		if (sx_geo_aside(ga, gb, c->eps) || sx_geo_aside(gb, ga, c->eps))
			return 0;
	}
	if ((A->n[0] == 0 && A->n[1] == 0 && A->n[2] == 0) || (B->n[0] == 0 && B->n[1] == 0 && B->n[2] == 0))
		return 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			shared += A->c[i] == B->c[j];
	if (shared >= 3)
		return 0;
	for (i = 0; i < 3; i++)
	{
		double da = sx_dot(B->n, c->X + 3 * (size_t)A->c[i]) - B->d, db = sx_dot(A->n, c->X + 3 * (size_t)B->c[i]) - A->d;
		sa[i] = da > c->eps ? 1 : da < -c->eps ? -1 : 0;
		sb[i] = db > c->eps ? 1 : db < -c->eps ? -1 : 0;
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
		return 0;
	if (za == 3 || zb == 3)
		return !sx_separated(c, f, g);   /* (the fold test of sx_pair stays there) */
	if (shared == 2)
		return 0;
	n0 = sx_section_at(c, f, g, sa, xa);
	n1 = sx_section_at(c, g, f, sb, xb);
	if (!n0 || !n1)
		return 0;
	sx_cross(A->n, B->n, dir);
	l = sqrt(sx_dot(dir, dir));
	if (!(l > 0))
		return 0;
	s[0] = sx_dot(xa[0], dir) / l;
	s[1] = sx_dot(xa[n0 - 1], dir) / l;
	s[2] = sx_dot(xb[0], dir) / l;
	s[3] = sx_dot(xb[n1 - 1], dir) / l;
	{
		/* the two intervals along the planes' line overlap (within the tolerance) */
		double lo = fmax(fmin(s[0], s[1]), fmin(s[2], s[3])), hi = fmin(fmax(s[0], s[1]), fmax(s[2], s[3]));
		double const* x;
		if (lo > hi + 2 * c->eps)
			return 0;
		if (hi - lo > c->eps || !shared)
			return 1;
		/* touching at a point: no meeting when it is a shared corner (neighbours round a vertex) */
		x = fmin(s[0], s[1]) >= fmin(s[2], s[3]) ? xa[s[0] <= s[1] ? 0 : n0 - 1] : xb[s[2] <= s[3] ? 0 : n1 - 1];
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				if (A->c[i] == B->c[j] && sx_dist2(x, c->X + 3 * (size_t)A->c[i]) <= 4 * c->eps * c->eps)
					return 0;
		return 1;
	}
}

typedef struct sx_found
{
	unsigned long long* p;            /* f << 32 | g */
	unsigned int n, cap;
} sx_found;

typedef struct sx_filter
{
	sx_ctx const* c;
	sx_found* chunk;
} sx_filter;

static qaws_status sx_filter_visit(void* user, unsigned int chunk, unsigned int i, unsigned int j)
{
	sx_filter* F = (sx_filter*)user;
	sx_found* o = F->chunk + chunk;
	if (!sx_may_meet(F->c, i, j))
		return QAWS_STATUS_OK;
	if (!sx_reserve((void**)&o->p, &o->cap, (size_t)o->n + 1, sizeof(unsigned long long)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	o->p[o->n++] = (unsigned long long)(i < j ? i : j) << 32 | (i < j ? j : i);
	return QAWS_STATUS_OK;
}

/* every pair of boxes overlapping, filtered through the executor; the survivors cut in order (results the same
   whatever the executor) */
static qaws_status sx_pairs(sx_ctx* c, qaws_bp_box const* box, qaws_batch_executor const* executor)
{
	qaws_bp_grid* grid = NULL;
	sx_filter F;
	unsigned int nch, k, i, n = 0;
	unsigned long long* all = NULL;
	qaws_status st = qaws_internal_grid_create_fit(&grid, box, c->nf, 3, executor);
	if (st != QAWS_STATUS_OK)
		return st;
	nch = qaws_internal_grid_pair_chunks(grid);
	F.c = c;
	F.chunk = (sx_found*)qaws_internal_zalloc((size_t)nch + 1, sizeof(sx_found));
	if (!F.chunk)
	{
		qaws_internal_grid_destroy(grid);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_internal_grid_pairs(NULL, grid, box, NULL, NULL, sx_filter_visit, &F, executor);
	for (k = 0; k < nch; k++)
		n += F.chunk[k].n;
	if (st == QAWS_STATUS_OK)
	{
		all = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * ((size_t)n + 1));
		if (!all)
			st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (st == QAWS_STATUS_OK)
	{
		for (k = 0, n = 0; k < nch; k++)
		{
			memcpy(all + n, F.chunk[k].p, sizeof(unsigned long long) * F.chunk[k].n);
			n += F.chunk[k].n;
		}
		qsort(all, n, sizeof *all, sx_cmp_u64);
		for (i = 0; i < n && st == QAWS_STATUS_OK; i++)
			st = sx_pair(c, (unsigned int)(all[i] >> 32), (unsigned int)all[i]);
	}
	for (k = 0; k < nch; k++)
		qaws_internal_free(F.chunk[k].p);
	qaws_internal_free(F.chunk);
	qaws_internal_free(all);
	qaws_internal_grid_destroy(grid);
	return st;
}


static double sx_snap52(double t)
{
	t = fmin(fmax(t, 0.0), 1.0);
	return floor(t * 4503599627370496.0 + 0.5) / 4503599627370496.0;
}

/* ---- the input read by chunks of faces (vertices for the finite check): checked and bounded, then, once welded, its
   triangles laid out at each chunk's first ---- */

#define SX_LOAD_GRAIN 16384

typedef struct sx_load
{
	qaws_scalar const* x;
	unsigned int const *fv, *fs;
	unsigned int nv;
	unsigned int* corner;              /* per chunk: its corners, then its first */
	unsigned int* tri;                 /* per chunk: its triangles (those kept once welded), then its first kept */
	double* lohi;                      /* per chunk: its corners' bounds */
	int* bad;                          /* per chunk */
	unsigned int* cells;               /* per vertex: its cell's hash and sides (sx_cell_sides) */
	sx_ctx* c;
	unsigned int const* alias;
	qaws_bp_box* box;
} sx_load;

static qaws_status sx_load_finite(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	unsigned int i;
	int bad = 0;
	for (i = 3 * begin; i < 3 * end; i++)
		bad |= !isfinite((double)L->x[i]);
	L->bad[chunk] = bad;
	return QAWS_STATUS_OK;
}

static qaws_status sx_load_sizes(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	unsigned int f, n = 0, t = 0;
	int bad = 0;
	for (f = begin; f < end; f++)
	{
		unsigned int s = L->fs ? L->fs[f] : 3;
		bad |= s < 3;
		n += s;
		t += s - 2;
	}
	L->corner[chunk] = n;
	L->tri[chunk] = t;
	L->bad[chunk] = bad;
	return QAWS_STATUS_OK;
}

static qaws_status sx_load_check(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	double* b = L->lohi + 6 * (size_t)chunk;
	unsigned int f, k, i, off = L->corner[chunk];
	for (i = 0; i < 3; i++)
	{
		b[i] = HUGE_VAL;
		b[3 + i] = -HUGE_VAL;
	}
	L->bad[chunk] = 0;
	for (f = begin; f < end; f++)
	{
		unsigned int s = L->fs ? L->fs[f] : 3;
		for (k = 0; k < s; k++)
		{
			unsigned int v = L->fv[off + k];
			if (v >= L->nv)
			{
				L->bad[chunk] = 1;
				return QAWS_STATUS_OK;
			}
			for (i = 0; i < 3; i++)
			{
				double x = (double)L->x[3 * (size_t)v + i];
				b[i] = x < b[i] ? x : b[i];
				b[3 + i] = x > b[3 + i] ? x : b[3 + i];
			}
		}
		off += s;
	}
	return QAWS_STATUS_OK;
}

static qaws_status sx_load_cells(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	unsigned int i;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		double p[3] = { (double)L->x[3 * (size_t)i], (double)L->x[3 * (size_t)i + 1], (double)L->x[3 * (size_t)i + 2] };
		L->cells[2 * (size_t)i] = sx_cell_sides(L->c, p, L->cells + 2 * (size_t)i + 1);
	}
	return QAWS_STATUS_OK;
}

/* the triangles of a chunk's faces left once welded (polygons fanned) */
static qaws_status sx_load_count(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	unsigned int f, k, off = L->corner[chunk], n = 0;
	for (f = begin; f < end; f++)
	{
		unsigned int s = L->fs ? L->fs[f] : 3;
		for (k = 1; k + 1 < s; k++)
		{
			unsigned int a = L->alias[L->fv[off]], b = L->alias[L->fv[off + k]], d = L->alias[L->fv[off + k + 1]];
			n += a != b && b != d && d != a;
		}
		off += s;
	}
	L->tri[chunk] = n;
	return QAWS_STATUS_OK;
}

/* a chunk's triangles with their planes, boxes and the compact records of the pair filter */
static qaws_status sx_load_fill(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_load const* L = (sx_load const*)user;
	sx_ctx* c = L->c;
	unsigned int f, k, i, off = L->corner[chunk], nf = L->tri[chunk];
	for (f = begin; f < end; f++)
	{
		unsigned int s = L->fs ? L->fs[f] : 3;
		for (k = 1; k + 1 < s; k++)
		{
			sx_face* F = c->F + nf;
			sx_geo* G = c->geo + nf;
			double u[3], v[3], l;
			unsigned int p0 = L->alias[L->fv[off]], p1 = L->alias[L->fv[off + k]], p2 = L->alias[L->fv[off + k + 1]];
			if (p0 == p1 || p1 == p2 || p2 == p0)
				continue;   /* welded away */
			F->c[0] = p0;
			F->c[1] = p1;
			F->c[2] = p2;
			F->src = f;
			F->cuts = SX_NONE;
			F->touch = SX_NONE;
			F->queued = F->touched = 0;
			for (i = 0; i < 3; i++)
			{
				u[i] = c->X[3 * (size_t)p1 + i] - c->X[3 * (size_t)p0 + i];
				v[i] = c->X[3 * (size_t)p2 + i] - c->X[3 * (size_t)p0 + i];
			}
			sx_cross(u, v, F->n);
			l = sqrt(sx_dot(F->n, F->n));
			/* lower than the tolerance over its longest edge: no plane, cut by nothing */
			if (l > c->eps * sqrt(sx_max3(sx_dot(u, u), sx_dot(v, v), sx_dist2(c->X + 3 * (size_t)p1, c->X + 3 * (size_t)p2))))
			{
				F->n[0] /= l;
				F->n[1] /= l;
				F->n[2] /= l;
			}
			else
				F->n[0] = F->n[1] = F->n[2] = 0;
			F->d = sx_dot(F->n, c->X + 3 * (size_t)p0);
			F->e[0] = F->e[1] = F->e[2] = SX_NONE;
			for (i = 0; i < 3; i++)
			{
				double x0 = c->X[3 * (size_t)p0 + i], x1 = c->X[3 * (size_t)p1 + i], x2 = c->X[3 * (size_t)p2 + i];
				double l0 = x0 < x1 ? x0 : x1, h0 = x0 < x1 ? x1 : x0;
				L->box[nf].lo[i] = (l0 < x2 ? l0 : x2) - c->eps;
				L->box[nf].hi[i] = (h0 > x2 ? h0 : x2) + c->eps;
			}
			for (i = 0; i < 3; i++)
			{
				memcpy(G->x[i], c->X + 3 * (size_t)F->c[i], 3 * sizeof(double));
				G->c[i] = F->c[i];
			}
			memcpy(G->n, F->n, sizeof G->n);
			G->d = F->d;
			G->pad = 0;
			nf++;
		}
		off += s;
	}
	return QAWS_STATUS_OK;
}

/* ---- the faces' triangles laid out by chunks: a face as given, or its pieces made aside ---- */

#define SX_CUT_PIECES 0x80000000u    /* (a face's pieces cut it) */

typedef struct sx_lay
{
	sx_face const* F;
	unsigned int const* oc;            /* per face touched: its pieces, | SX_CUT_PIECES when cut */
	unsigned int const* piece;         /* the touched faces' pieces, face by face */
	unsigned int* at;                  /* per chunk: triangles, cut pieces, pieces; then the first of each */
	unsigned int *tri, *cut;
	qaws_scalar* src;
} sx_lay;

static qaws_status sx_lay_count(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_lay const* Y = (sx_lay const*)user;
	unsigned int f, n[3] = { 0, 0, 0 };
	for (f = begin; f < end; f++)
	{
		sx_face const* F = Y->F + f;
		unsigned int k;
		if (F->c[0] == F->c[1] || F->c[1] == F->c[2] || F->c[2] == F->c[0])
			continue;
		k = F->touched ? Y->oc[f] & ~SX_CUT_PIECES : 1;
		n[0] += k;
		n[1] += F->touched && (Y->oc[f] & SX_CUT_PIECES) ? k : 0;
		n[2] += F->touched ? k : 0;
	}
	memcpy(Y->at + 3 * (size_t)chunk, n, sizeof n);
	return QAWS_STATUS_OK;
}

static qaws_status sx_lay_fill(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_lay const* Y = (sx_lay const*)user;
	unsigned int f, t = Y->at[3 * (size_t)chunk], q = Y->at[3 * (size_t)chunk + 1], p = Y->at[3 * (size_t)chunk + 2];
	for (f = begin; f < end; f++)
	{
		sx_face const* F = Y->F + f;
		unsigned int k, i;
		if (F->c[0] == F->c[1] || F->c[1] == F->c[2] || F->c[2] == F->c[0])
			continue;
		if (!F->touched)
		{
			/* met nothing, nothing on its edges: as given */
			memcpy(Y->tri + 3 * (size_t)t, F->c, 3 * sizeof(unsigned int));
			Y->src[t++] = (qaws_scalar)F->src;
			continue;
		}
		k = Y->oc[f] & ~SX_CUT_PIECES;
		memcpy(Y->tri + 3 * (size_t)t, Y->piece + 3 * (size_t)p, 3 * sizeof(unsigned int) * k);
		for (i = 0; i < k; i++)
		{
			if (Y->oc[f] & SX_CUT_PIECES)
				Y->cut[q++] = t;
			Y->src[t++] = (qaws_scalar)F->src;
		}
		p += k;
	}
	return QAWS_STATUS_OK;
}

/* the arrangement as arrays (owned by the caller) */
typedef struct sx_out
{
	qaws_scalar* pos;
	unsigned int* tri;
	qaws_scalar* src;                  /* input face per triangle */
	unsigned char* copl;               /* per input face: in a coplanar overlap with another */
	unsigned int nv, nt;
	double eps, diag;
} sx_out;

static qaws_status sx_resolve(sx_out* o, qaws_intersect_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_intersect_desc const* desc)
{
	qaws_intersect_report local;
	sx_ctx c;
	sx_buf b;
	sx_cdt T;
	qaws_status st = QAWS_STATUS_OK;
	unsigned int *alias = NULL, *loc = NULL, *otri = NULL, *stack = NULL, *used = NULL, *ocut = NULL, *oc = NULL, *ptri = NULL, ncut = 0, npiece = 0;
	qaws_scalar *opos = NULL, *osrc = NULL;
	qaws_bp_box* box = NULL;
	qaws_batch_executor const* ex = desc ? desc->executor : NULL;
	sx_load L;
	unsigned int nt = 0, off = 0, f, k, i, ntri = 0, capt = 0, caps2 = 0, onv = 0, nchf, nchv;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, diag;
	memset(&L, 0, sizeof L);
	memset(&c, 0, sizeof c);
	memset(&b, 0, sizeof b);
	memset(&T, 0, sizeof T);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	memset(o, 0, sizeof *o);
	if (!face_count || !face_vertices || (vertex_count && !positions) || (desc && !(desc->tolerance >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* the input checked by chunks: finite, sizes, indices (bounds alongside) */
	nchf = qaws_internal_chunk_count(face_count, SX_LOAD_GRAIN);
	nchv = qaws_internal_chunk_count(vertex_count, SX_LOAD_GRAIN);
	L.x = positions;
	L.fv = face_vertices;
	L.fs = face_sizes;
	L.nv = vertex_count;
	L.corner = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nchf + 1));
	L.tri = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nchf + 1));
	L.lohi = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)nchf + 1));
	L.bad = (int*)qaws_internal_zalloc((size_t)(nchf > nchv ? nchf : nchv) + 1, sizeof(int));
	L.c = &c;
	L.alias = NULL;
	L.box = NULL;
	if (!L.corner || !L.tri || !L.lohi || !L.bad)
		goto oom;
	st = qaws_internal_parallel(vertex_count, SX_LOAD_GRAIN, sx_load_finite, &L, ex);
	for (k = 0; st == QAWS_STATUS_OK && k < nchv; k++)
		if (L.bad[k])
			st = QAWS_STATUS_INVALID_ARGUMENT;
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_parallel(face_count, SX_LOAD_GRAIN, sx_load_sizes, &L, ex);
	for (k = 0; st == QAWS_STATUS_OK && k < nchf; k++)
	{
		unsigned int n = L.corner[k];
		if (L.bad[k])
			st = QAWS_STATUS_INVALID_ARGUMENT;
		L.corner[k] = off;
		off += n;
		nt += L.tri[k];
	}
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_parallel(face_count, SX_LOAD_GRAIN, sx_load_check, &L, ex);
	for (k = 0; st == QAWS_STATUS_OK && k < nchf; k++)
	{
		if (L.bad[k])
			st = QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < 3; i++)
		{
			lo[i] = L.lohi[6 * (size_t)k + i] < lo[i] ? L.lohi[6 * (size_t)k + i] : lo[i];
			hi[i] = L.lohi[6 * (size_t)k + 3 + i] > hi[i] ? L.lohi[6 * (size_t)k + 3 + i] : hi[i];
		}
	}
	if (st != QAWS_STATUS_OK)
		goto fail;
	diag = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	c.eps = desc && desc->tolerance > 0 ? (double)desc->tolerance : 1e-9 * diag;
	if (!(c.eps > 0))
		c.eps = 1e-300;
	report->tolerance = (qaws_scalar)c.eps;
	o->eps = c.eps;
	o->diag = diag;
	c.inv = 1 / (64 * c.eps);   /* cells of 64 tolerances: a point's neighbours within it mostly in its own cell */
	c.r = report;
	c.copl = (unsigned char*)qaws_internal_zalloc((size_t)face_count + 1, 1);
	c.ncell = 1024;
	c.cell = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * c.ncell);
	alias = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)vertex_count + 1));
	box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (size_t)nt + 1));
	if (!c.copl || !c.cell || !alias || !box || !sx_reserve((void**)&c.X, &c.capx, 3 * ((size_t)vertex_count + 1), sizeof(double))
		|| !sx_reserve((void**)&c.chain, &c.capc, (size_t)vertex_count + 1, sizeof(unsigned int))
		|| !sx_reserve((void**)&c.hash, &c.caph, (size_t)vertex_count + 1, sizeof(unsigned int)))
		goto oom;
	for (i = 0; i < c.ncell; i++)
		c.cell[i] = SX_NONE;
	/* vertices: those within the tolerance of an earlier one are welded (their cells found by chunks first) */
	L.cells = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)vertex_count + 1));
	if (!L.cells)
		goto oom;
	st = qaws_internal_parallel(vertex_count, SX_LOAD_GRAIN, sx_load_cells, &L, ex);
	if (st != QAWS_STATUS_OK)
		goto fail;
	for (i = 0; i < vertex_count; i++)
	{
		double p[3] = { (double)positions[3 * (size_t)i], (double)positions[3 * (size_t)i + 1], (double)positions[3 * (size_t)i + 2] };
		unsigned int before = c.np;
		alias[i] = sx_point_at(&c, p, L.cells[2 * (size_t)i], L.cells[2 * (size_t)i + 1]);
		if (c.oom)
			goto oom;
		report->welded_vertices += c.np == before;
	}
	onv = c.np;
	/* triangles (polygons fanned) with their planes, boxes and the compact records of the pair filter, by chunks
	   laid out in order */
	L.alias = alias;
	L.box = box;
	st = qaws_internal_parallel(face_count, SX_LOAD_GRAIN, sx_load_count, &L, ex);
	if (st != QAWS_STATUS_OK)
		goto fail;
	for (k = 0; k < nchf; k++)
	{
		unsigned int n = L.tri[k];
		L.tri[k] = c.nf;
		c.nf += n;
	}
	c.F = (sx_face*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sx_face), (size_t)c.nf + 1));
	c.geo = (sx_geo*)qaws_internal_malloc(qaws_internal_mul2(sizeof(sx_geo), (size_t)c.nf + 1));
	if (!c.F || !c.geo)
		goto oom;
	memset(c.F + c.nf, 0, sizeof(sx_face));
	st = qaws_internal_parallel(face_count, SX_LOAD_GRAIN, sx_load_fill, &L, ex);
	if (st != QAWS_STATUS_OK)
		goto fail;
	/* the faces round each point (counting sort) */
	c.vf_start = (unsigned int*)qaws_internal_zalloc((size_t)c.np + 2, sizeof(unsigned int));
	c.vf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)c.nf + 1));
	if (!c.vf_start || !c.vf)
		goto oom;
	for (f = 0; f < c.nf; f++)
		for (k = 0; k < 3; k++)
			c.vf_start[c.F[f].c[k] + 2]++;
	for (i = 2; i < c.np + 2; i++)
		c.vf_start[i] += c.vf_start[i - 1];
	for (f = 0; f < c.nf; f++)
		for (k = 0; k < 3; k++)
			c.vf[c.vf_start[c.F[f].c[k] + 1]++] = f;
	/* the pairs cut each other */
	if (c.nf > 1)
		st = sx_pairs(&c, box, desc ? desc->executor : NULL);
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
	/* the faces touched first, their pieces aside (the others are as given: laid out with them by chunks after) */
	oc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)c.nf + 1));
	if (!oc)
		goto oom;
	for (f = 0; f < c.nf; f++)
	{
		sx_face const* F = c.F + f;
		unsigned int np, need, le[3];
		if (F->c[0] == F->c[1] || F->c[1] == F->c[2] || F->c[2] == F->c[0] || !F->touched)
			continue;
		sx_edges_found(&c, f, le);   /* (none made here: a face's edges would reach its neighbours in turn) */
		sx_gather(&c, &b, f, le);
		if (c.oom)
			goto oom;
		np = b.npt;
		need = np > 3 ? 2 * np : 1;
		if (!sx_reserve((void**)&ptri, &capt, 3 * ((size_t)npiece + need + 1), sizeof(unsigned int)))
			goto oom;
		if (np == 3)
		{
			memcpy(ptri + 3 * (size_t)npiece++, F->c, 3 * sizeof(unsigned int));
			oc[f] = 1;
			continue;
		}
		report->split_faces++;
		if (!sx_reserve((void**)&T.P, &T.capp, 2 * (size_t)np, sizeof(double))
			|| !sx_reserve((void**)&T.Q, &T.capq, 2 * (size_t)np, sizeof(double))
			|| !sx_reserve((void**)&T.t, &T.capt, 2 * (size_t)np + 4, sizeof(sx_tri))
			|| !sx_reserve((void**)&T.vt, &T.capvt, (size_t)np + 1, sizeof(unsigned int))
			|| !sx_reserve((void**)&T.st, &T.capst, 2 * (size_t)np + 4, sizeof(unsigned int))
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
				for (j = le[k] == SX_NONE ? SX_NONE : c.car[le[k]].points; j != SX_NONE; j = c.nodes[j].next)
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
		T.vton = 0;
		sx_set(T.t, 0, 1, 2, SX_NONE, SX_NONE, SX_NONE, 0, 0, 0);
		for (i = 3; i < np; i++)
			sx_insert(&T, i, stack, 16 * np + 64);
		/* each point's triangle, kept through the cuts' flips */
		for (i = 0; i < T.nt; i++)
		{
			T.vt[T.t[i].v[0]] = i;
			T.vt[T.t[i].v[1]] = i;
			T.vt[T.t[i].v[2]] = i;
		}
		T.vton = 1;
		for (i = 0; i < b.ns; i++)
		{
			unsigned int pa = loc[b.s[i].p], pb = loc[b.s[i].q];
			int r;
			if (pa == pb || b.s[i].k == le[0] || b.s[i].k == le[1] || b.s[i].k == le[2])
				continue;   /* the face's own edges are its border */
			if ((r = sx_recover(&T, pa, pb, 0)) < 0)
				goto oom;
			if (!r)
				report->unresolved++;
		}
		T.vton = 0;
		sx_relax(&T);
		for (i = 0; i < T.nt; i++)
		{
			ptri[3 * (size_t)npiece] = b.p[T.t[i].v[0]];
			ptri[3 * (size_t)npiece + 1] = b.p[T.t[i].v[1]];
			ptri[3 * (size_t)npiece + 2] = b.p[T.t[i].v[2]];
			npiece++;
		}
		oc[f] = T.nt | SX_CUT_PIECES;
		for (i = 0; i < np; i++)
			loc[b.p[i]] = SX_NONE;
	}
	{
		/* every face's triangles in order: counted by chunks, then laid out */
		unsigned int nch = qaws_internal_chunk_count(c.nf, SX_LOAD_GRAIN), sum[3] = { 0, 0, 0 };
		sx_lay Y;
		Y.F = c.F;
		Y.oc = oc;
		Y.piece = ptri;
		Y.at = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nch + 1));
		if (!Y.at)
			goto oom;
		st = qaws_internal_parallel(c.nf, SX_LOAD_GRAIN, sx_lay_count, &Y, ex);
		for (k = 0; k < nch; k++)
			for (i = 0; i < 3; i++)
			{
				unsigned int n = Y.at[3 * (size_t)k + i];
				Y.at[3 * (size_t)k + i] = sum[i];
				sum[i] += n;
			}
		ntri = sum[0];
		ncut = sum[1];
		otri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)ntri + 1));
		osrc = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)ntri + 1));
		ocut = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)ncut + 1));
		Y.tri = otri;
		Y.src = osrc;
		Y.cut = ocut;
		if (st == QAWS_STATUS_OK && (!otri || !osrc || !ocut))
			st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (st == QAWS_STATUS_OK)
			st = qaws_internal_parallel(c.nf, SX_LOAD_GRAIN, sx_lay_fill, &Y, ex);
		qaws_internal_free(Y.at);
		if (st != QAWS_STATUS_OK)
			goto fail;
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
		/* the cut faces' pieces checked by exact predicates (the others are input faces the pairs tested; pieces of
		   one face are one planar triangulation) */
		if (ncut)
		{
			st = sx_check(&report->remaining, opos, otri, osrc, ocut, ncut, c.eps, desc ? desc->executor : NULL);
			if (st != QAWS_STATUS_OK)
				goto fail;
		}
		o->pos = opos;
		o->tri = otri;
		o->src = osrc;
		o->copl = c.copl;
		c.copl = NULL;
		o->nv = nv;
		o->nt = ntri;
		opos = NULL;
		otri = NULL;
		osrc = NULL;
	}
	goto fail;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
fail:
	qaws_internal_free(c.copl);
	qaws_internal_free(c.X); qaws_internal_free(c.chain); qaws_internal_free(c.hash); qaws_internal_free(c.cell); qaws_internal_free(c.keys); qaws_internal_free(c.nodes);
	qaws_internal_free(c.car); qaws_internal_free(c.F); qaws_internal_free(c.vf_start); qaws_internal_free(c.vf); qaws_internal_free(c.geo); qaws_internal_free(c.work); qaws_internal_free(b.s); qaws_internal_free(b.sp); qaws_internal_free(b.p);
	qaws_internal_free(T.P); qaws_internal_free(T.Q); qaws_internal_free(T.t); qaws_internal_free(T.vt); qaws_internal_free(T.st); qaws_internal_free(T.qe); qaws_internal_free(alias); qaws_internal_free(loc); qaws_internal_free(used);
	qaws_internal_free(otri); qaws_internal_free(ocut); qaws_internal_free(stack); qaws_internal_free(opos); qaws_internal_free(osrc);
	qaws_internal_free(oc); qaws_internal_free(ptri);
	qaws_internal_free(box);
	qaws_internal_free(L.corner); qaws_internal_free(L.tri); qaws_internal_free(L.lohi); qaws_internal_free(L.bad); qaws_internal_free(L.cells);
	return st;
}

/* triangles and their source faces as a mesh with the "source" attribute; opposite (may be NULL): every corner's
   opposite, known (the mesh made without looking its edges up) */
static qaws_status sx_mesh(qaws_mesh_3d** out, qaws_scalar const* pos, unsigned int nv, unsigned int const* tri, unsigned int const* opposite,
	qaws_scalar const* src, unsigned int nt)
{
	qaws_mesh_3d_desc md;
	qaws_attribute_desc ad;
	unsigned int attr;
	qaws_status st;
	memset(&md, 0, sizeof md);
	md.positions = pos;
	md.vertex_count = nv;
	md.face_vertices = tri;
	md.face_count = nt;
	if (opposite)
	{
		double* x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
		qaws_mesh_3d* m = (qaws_mesh_3d*)qaws_internal_zalloc(1, sizeof *m);
		unsigned int i;
		*out = NULL;
		if (!x || !m)
		{
			qaws_internal_free(x);
			qaws_internal_free(m);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < 3 * nv; i++)
			x[i] = (double)pos[i];
		st = hm_build_paired(&m->hm, 3, x, nv, tri, opposite, nt);
		qaws_internal_free(x);
		if (st != QAWS_STATUS_OK)
		{
			qaws_internal_free(m);
			return st;
		}
		*out = m;
	}
	else
		st = qaws_mesh_3d_create(out, &md);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&ad, 0, sizeof ad);
	ad.name = "source";
	ad.element = QAWS_MESH_FACE;
	ad.components = 1;
	ad.rule = QAWS_ATTRIBUTE_NEAREST;
	st = qaws_mesh_3d_add_attribute(&attr, *out, src, &ad);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	return st;
}

static void sx_out_free(sx_out* o)
{
	qaws_internal_free(o->pos);
	qaws_internal_free(o->tri);
	qaws_internal_free(o->src);
	qaws_internal_free(o->copl);
	memset(o, 0, sizeof *o);
}

qaws_status qaws_mesh_resolve_intersections(qaws_mesh_3d** out, qaws_intersect_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_intersect_desc const* desc)
{
	sx_out o;
	qaws_status st;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	st = sx_resolve(&o, report, positions, face_vertices, face_sizes, vertex_count, face_count, desc);
	if (st == QAWS_STATUS_OK)
		st = sx_mesh(out, o.pos, o.nv, o.tri, NULL, o.src, o.nt);
	sx_out_free(&o);
	return st;
}

/* ---- Booleans on the arrangement ---- */

/* ---- winding numbers summed directly (few points): exact solid angles over an operand's triangles ---- */

/* squared distance from p to triangle abc (Ericson's regions) */
static double sx_tri_d2(double const* p, double const* a, double const* b, double const* c)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], q[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, den;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = b[k] - a[k];
		ac[k] = c[k] - a[k];
		ap[k] = p[k] - a[k];
		bp[k] = p[k] - b[k];
		cp[k] = p[k] - c[k];
	}
	d1 = sx_dot(ab, ap);
	d2 = sx_dot(ac, ap);
	if (d1 <= 0 && d2 <= 0)
		return sx_dist2(p, a);
	d3 = sx_dot(ab, bp);
	d4 = sx_dot(ac, bp);
	if (d3 >= 0 && d4 <= d3)
		return sx_dist2(p, b);
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0)
	{
		v = d1 / (d1 - d3);
		for (k = 0; k < 3; k++)
			q[k] = a[k] + v * ab[k];
		return sx_dist2(p, q);
	}
	d5 = sx_dot(ab, cp);
	d6 = sx_dot(ac, cp);
	if (d6 >= 0 && d5 <= d6)
		return sx_dist2(p, c);
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0)
	{
		w = d2 / (d2 - d6);
		for (k = 0; k < 3; k++)
			q[k] = a[k] + w * ac[k];
		return sx_dist2(p, q);
	}
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0)
	{
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		for (k = 0; k < 3; k++)
			q[k] = b[k] + w * (c[k] - b[k]);
		return sx_dist2(p, q);
	}
	den = 1 / (va + vb + vc);
	v = vb * den;
	w = vc * den;
	for (k = 0; k < 3; k++)
		q[k] = a[k] + ab[k] * v + ac[k] * w;
	return sx_dist2(p, q);
}

typedef struct sx_direct
{
	qaws_scalar const* x;              /* the operand's positions */
	unsigned int const* t;             /* its triangles */
	double const* p;                   /* points, 3 each */
	unsigned int np;
	double r2;                         /* near: within sqrt(r2) of a triangle (0: not asked) */
	double* w;                         /* per chunk, np solid-angle sums */
	unsigned char* near;               /* per chunk, np flags */
} sx_direct;

static qaws_status sx_direct_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_direct const* D = (sx_direct const*)user;
	double* w = D->w + (size_t)chunk * D->np;
	unsigned char* nr = D->near + (size_t)chunk * D->np;
	unsigned int f, i, k;
	/* point by point, summed locally: the chunks' sums share cache lines */
	for (i = 0; i < D->np; i++)
	{
		double const* q = D->p + 3 * (size_t)i;
		double s = 0;
		unsigned char hit = 0;
		for (f = begin; f < end; f++)
		{
			double A[3], B[3], C[3], a[3], b[3], c[3], bc[3], la, lb, lc, num, den;
			for (k = 0; k < 3; k++)
			{
				A[k] = (double)D->x[3 * (size_t)D->t[3 * (size_t)f] + k];
				B[k] = (double)D->x[3 * (size_t)D->t[3 * (size_t)f + 1] + k];
				C[k] = (double)D->x[3 * (size_t)D->t[3 * (size_t)f + 2] + k];
				a[k] = A[k] - q[k];
				b[k] = B[k] - q[k];
				c[k] = C[k] - q[k];
			}
			/* the triangle's solid angle (Van Oosterom and Strackee) */
			sx_cross(b, c, bc);
			la = sqrt(sx_dot(a, a));
			lb = sqrt(sx_dot(b, b));
			lc = sqrt(sx_dot(c, c));
			num = sx_dot(a, bc);
			den = la * lb * lc + sx_dot(a, b) * lc + sx_dot(b, c) * la + sx_dot(c, a) * lb;
			s += 2 * atan2(num, den);
			if (D->r2 > 0 && !hit && sx_tri_d2(q, A, B, C) <= D->r2)
				hit = 1;
		}
		w[i] = s;
		nr[i] = hit;
	}
	return QAWS_STATUS_OK;
}

#define SX_DIRECT_GRAIN 4096

/* winding numbers of the triangles t (n) at the points p (np; out_w may be NULL), and which points are within r of
   them (out_near may be NULL); the chunks' sums added in order, so the result does not depend on the executor */
static qaws_status sx_direct_winding(qaws_scalar* out_w, unsigned char* out_near, qaws_scalar const* x, unsigned int const* t, unsigned int n,
	double const* p, unsigned int np, double r, qaws_batch_executor const* executor)
{
	sx_direct D;
	unsigned int nch = qaws_internal_chunk_count(n, SX_DIRECT_GRAIN), c, i;
	qaws_status st;
	D.x = x;
	D.t = t;
	D.p = p;
	D.np = np;
	D.r2 = out_near ? r * r : 0;
	D.w = (double*)qaws_internal_zalloc(qaws_internal_mul2((size_t)nch + 1, np), sizeof(double));
	D.near = (unsigned char*)qaws_internal_zalloc(qaws_internal_mul2((size_t)nch + 1, np), 1);
	if (!D.w || !D.near)
	{
		qaws_internal_free(D.w);
		qaws_internal_free(D.near);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	st = qaws_internal_parallel(n, SX_DIRECT_GRAIN, sx_direct_job, &D, executor);
	for (i = 0; i < np; i++)
	{
		double s = 0;
		unsigned char nr = 0;
		for (c = 0; c < nch; c++)
		{
			s += D.w[(size_t)c * np + i];
			nr |= D.near[(size_t)c * np + i];
		}
		if (out_w)
			out_w[i] = (qaws_scalar)(s / (4 * 3.14159265358979323846));
		if (out_near)
			out_near[i] = nr;
	}
	qaws_internal_free(D.w);
	qaws_internal_free(D.near);
	return st;
}

#define SX_RADIX_BITS 11
#define SX_RADIX_SIZE (1u << SX_RADIX_BITS)
#define SX_RADIX_GRAIN 65536

typedef struct sx_sort
{
	unsigned long long *src, *dst;
	unsigned long long* diff;          /* per chunk: the bits where its keys differ from the first */
	unsigned int* count;               /* per chunk: the digit's buckets, then where its records go */
	unsigned int stride, shift;
} sx_sort;

static qaws_status sx_sort_diff(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_sort const* S = (sx_sort const*)user;
	unsigned long long first = S->src[0], d = 0;
	unsigned int i;
	for (i = begin; i < end; i++)
		d |= S->src[(size_t)i * S->stride] ^ first;
	S->diff[chunk] = d;
	return QAWS_STATUS_OK;
}

static qaws_status sx_sort_count(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_sort const* S = (sx_sort const*)user;
	unsigned int *c = S->count + (size_t)chunk * SX_RADIX_SIZE, i;
	memset(c, 0, sizeof(unsigned int) * SX_RADIX_SIZE);
	for (i = begin; i < end; i++)
		c[(S->src[(size_t)i * S->stride] >> S->shift) & (SX_RADIX_SIZE - 1)]++;
	return QAWS_STATUS_OK;
}

static qaws_status sx_sort_scatter(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_sort const* S = (sx_sort const*)user;
	unsigned int *c = S->count + (size_t)chunk * SX_RADIX_SIZE, i;
	unsigned long long const* src = S->src;
	unsigned long long* dst = S->dst;
	if (S->stride == 1)
		for (i = begin; i < end; i++)
			dst[c[(src[i] >> S->shift) & (SX_RADIX_SIZE - 1)]++] = src[i];
	else
		for (i = begin; i < end; i++)
		{
			unsigned int o = c[(src[2 * (size_t)i] >> S->shift) & (SX_RADIX_SIZE - 1)]++;
			dst[2 * (size_t)o] = src[2 * (size_t)i];
			dst[2 * (size_t)o + 1] = src[2 * (size_t)i + 1];
		}
	return QAWS_STATUS_OK;
}

/* LSD radix sort of n records of stride (1 or 2) 64-bit words by their first word (stable, so the same whatever the
   executor): digits only over the bits where keys differ (an edge key's zero bits between its vertices skipped), each
   pass counted and scattered by chunks; 0 when out of memory */
static int sx_radix(unsigned long long* a, unsigned int n, unsigned int stride, qaws_batch_executor const* executor)
{
	sx_sort S;
	unsigned long long *tmp, mask = 0;
	unsigned int nch = qaws_internal_chunk_count(n, SX_RADIX_GRAIN), i, k;
	int ok = 1;
	if (n < 2)
		return 1;
	tmp = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), stride, n));
	S.count = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), nch, SX_RADIX_SIZE));
	S.diff = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * nch);
	S.src = a;
	S.dst = tmp;
	S.stride = stride;
	S.shift = 0;
	ok = tmp && S.count && S.diff && qaws_internal_parallel(n, SX_RADIX_GRAIN, sx_sort_diff, &S, executor) == QAWS_STATUS_OK;
	for (k = 0; ok && k < nch; k++)
		mask |= S.diff[k];
	while (ok && mask)
	{
		unsigned long long* sw;
		unsigned int sum = 0;
		while (!((mask >> S.shift) & 1))
			S.shift++;
		ok = qaws_internal_parallel(n, SX_RADIX_GRAIN, sx_sort_count, &S, executor) == QAWS_STATUS_OK;
		/* bucket by bucket, chunk by chunk: each chunk's records of a digit after the earlier chunks' */
		for (i = 0; ok && i < SX_RADIX_SIZE; i++)
			for (k = 0; k < nch; k++)
			{
				unsigned int t = S.count[(size_t)k * SX_RADIX_SIZE + i];
				S.count[(size_t)k * SX_RADIX_SIZE + i] = sum;
				sum += t;
			}
		ok = ok && qaws_internal_parallel(n, SX_RADIX_GRAIN, sx_sort_scatter, &S, executor) == QAWS_STATUS_OK;
		sw = S.src;
		S.src = S.dst;
		S.dst = sw;
		mask = S.shift + SX_RADIX_BITS < 64 ? mask >> (S.shift + SX_RADIX_BITS) << (S.shift + SX_RADIX_BITS) : 0;
		S.shift += SX_RADIX_BITS;
	}
	if (ok && S.src != a)
		memcpy(a, S.src, qaws_internal_mul3(sizeof(unsigned long long), stride, n));
	qaws_internal_free(tmp);
	qaws_internal_free(S.count);
	qaws_internal_free(S.diff);
	return ok;
}

/* a mesh's faces fanned into triangles, for its BVH */
/* half the distance between patch i / 2's two samples (its offset) */
static double sx_offset_of(qaws_scalar const* pts, unsigned int i)
{
	unsigned int pi = i / 2, k;
	double s = 0;
	for (k = 0; k < 3; k++)
	{
		double d = (double)pts[6 * (size_t)pi + k] - (double)pts[6 * (size_t)pi + 3 + k];
		s += d * d;
	}
	return 0.5 * sqrt(s);
}

static qaws_status sx_bvh(qaws_mesh_bvh** out, qaws_mesh_3d_desc const* m, unsigned int** tri)
{
	qaws_mesh_bvh_desc d;
	unsigned int f, k, off = 0, n = 0;
	for (f = 0; f < m->face_count; f++)
		n += (m->face_sizes ? m->face_sizes[f] : 3) - 2;
	*tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)n + 1));
	if (!*tri)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0, n = 0; f < m->face_count; f++)
	{
		unsigned int s = m->face_sizes ? m->face_sizes[f] : 3;
		for (k = 1; k + 1 < s; k++, n++)
		{
			(*tri)[3 * (size_t)n] = m->face_vertices[off];
			(*tri)[3 * (size_t)n + 1] = m->face_vertices[off + k];
			(*tri)[3 * (size_t)n + 2] = m->face_vertices[off + k + 1];
		}
		off += s;
	}
	memset(&d, 0, sizeof d);
	d.positions = m->positions;
	d.vertex_count = m->vertex_count;
	d.dimension = 3;
	d.simplices = *tri;
	d.simplex_size = 3;
	d.simplex_count = n;
	return qaws_mesh_bvh_create(out, &d);
}

static unsigned int sx_root(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
		i = par[i] = par[par[i]];
	return i;
}

static int sx_region(qaws_boolean_op op, int a, int b, int two)
{
	if (!two)
		return a;
	return op == QAWS_BOOLEAN_UNION ? (a || b) : op == QAWS_BOOLEAN_INTERSECTION ? (a && b) : (a && !b);
}

static int sx_cmp_edge(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

static int sx_cmp_quad(void const* p, void const* q)
{
	unsigned int const *a = (unsigned int const*)p, *b = (unsigned int const*)q;
	unsigned int k;
	for (k = 0; k < 4; k++)
		if (a[k] != b[k])
			return a[k] < b[k] ? -1 : 1;
	return 0;
}

/* dup[i] = 1 for every triangle with the same vertices in the same cyclic order as one from a lower input face (a coincident copy
   facing the same way: an operand lying on itself); 0 on allocation failure */
static int sx_dup_mark(unsigned char* dup, unsigned int const* tri, qaws_scalar const* src, unsigned int nt)
{
	unsigned int* q = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 5, (size_t)nt + 1));
	unsigned int i, r;
	if (!q)
		return 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int const* t = tri + 3 * (size_t)i;
		r = t[0] < t[1] ? (t[0] < t[2] ? 0 : 2) : (t[1] < t[2] ? 1 : 2);   /* (rotated to its smallest index first) */
		q[5 * (size_t)i] = t[r];
		q[5 * (size_t)i + 1] = t[(r + 1) % 3];
		q[5 * (size_t)i + 2] = t[(r + 2) % 3];
		q[5 * (size_t)i + 3] = (unsigned int)src[i];
		q[5 * (size_t)i + 4] = i;
	}
	qsort(q, nt, 5 * sizeof(unsigned int), sx_cmp_quad);
	memset(dup, 0, nt);
	for (i = 1; i < nt; i++)
		if (q[5 * (size_t)i] == q[5 * (size_t)i - 5] && q[5 * (size_t)i + 1] == q[5 * (size_t)i - 4] && q[5 * (size_t)i + 2] == q[5 * (size_t)i - 3])
			dup[q[5 * (size_t)i + 4]] = 1;
	qaws_internal_free(q);
	return 1;
}

/* an edge set (open addressing; 0 empty, ~0 removed) */
static unsigned long long* sx_eset_slot(unsigned long long* s, unsigned int cap, unsigned long long k, int insert)
{
	unsigned int h = (unsigned int)((k * 0x9E3779B97F4A7C15ull) >> 32) & (cap - 1);
	unsigned long long* tomb = NULL;
	while (s[h])
	{
		if (s[h] == k)
			return s + h;
		if (s[h] == ~0ull && !tomb)
			tomb = s + h;
		h = (h + 1) & (cap - 1);
	}
	return insert ? (tomb ? tomb : s + h) : NULL;
}

static unsigned long long sx_ekey(unsigned int p, unsigned int q)
{
	return p < q ? ((unsigned long long)p << 32 | q) + 1 : ((unsigned long long)q << 32 | p) + 1;
}

static double sx_angle(double const* o, double const* p, double const* q)
{
	double u[3] = { p[0] - o[0], p[1] - o[1], p[2] - o[2] }, v[3] = { q[0] - o[0], q[1] - o[1], q[2] - o[2] }, w[3];
	sx_cross(u, v, w);
	return atan2(sqrt(sx_dot(w, w)), sx_dot(u, v));
}

/* triangle corner positions as doubles */
static void sx_tri_x(qaws_scalar const* pos, unsigned int const* t, unsigned int f, double (*P)[3])
{
	unsigned int k, j;
	for (k = 0; k < 3; k++)
		for (j = 0; j < 3; j++)
			P[k][j] = (double)pos[3 * (size_t)t[3 * (size_t)f + k] + j];
}

/* edge e = 3 f + i (f's corners a b c, the neighbour's d) worth turning: both faces in one plane, the quad convex, the
   opposite angles over pi, both faces of one operand (the caller checks the edge c d is new) */
static int sx_flat_flip_ok(unsigned int const* t, qaws_scalar const* src, unsigned int const* nb, qaws_scalar const* pos, unsigned int e, double eps,
	unsigned int split, int one)
{
	unsigned int f = e / 3, i = e % 3, g, j;
	double P[3][3], Q[3][3], n1[3], n2[3], u[3], v[3], l, a[3], b[3], c[3], d[3], s1, s2, w[3];
	unsigned int k;
	if (nb[e] == SX_NONE)
		return 0;
	g = nb[e] / 3;
	j = nb[e] % 3;
	if (((unsigned int)src[f] < split) != ((unsigned int)src[g] < split) || (one && src[f] != src[g]))
		return 0;
	sx_tri_x(pos, t, f, P);
	sx_tri_x(pos, t, g, Q);
	for (k = 0; k < 3; k++)
	{
		a[k] = P[i][k];
		b[k] = P[(i + 1) % 3][k];
		c[k] = P[(i + 2) % 3][k];
		d[k] = Q[(j + 2) % 3][k];
		u[k] = P[1][k] - P[0][k];
		v[k] = P[2][k] - P[0][k];
	}
	sx_cross(u, v, n1);
	l = sqrt(sx_dot(n1, n1));
	if (!(l > 0))
		return 0;
	for (k = 0; k < 3; k++)
	{
		n1[k] /= l;
		u[k] = Q[1][k] - Q[0][k];
		v[k] = Q[2][k] - Q[0][k];
		w[k] = d[k] - a[k];
	}
	sx_cross(u, v, n2);
	if (sx_dot(n1, n2) <= 0 || fabs(sx_dot(n1, w)) > eps)
		return 0;
	/* convex: c d separates a and b */
	for (k = 0; k < 3; k++)
	{
		u[k] = d[k] - c[k];
		v[k] = a[k] - c[k];
	}
	sx_cross(u, v, w);
	s1 = sx_dot(n1, w);
	for (k = 0; k < 3; k++)
		v[k] = b[k] - c[k];
	sx_cross(u, v, w);
	s2 = sx_dot(n1, w);
	if (!(s1 * s2 < 0) || fabs(s1) <= eps * sqrt(sx_dot(u, u)) || fabs(s2) <= eps * sqrt(sx_dot(u, u)))
		return 0;
	return sx_angle(c, a, b) + sx_angle(d, a, b) > 3.14159265358979323846 + 1e-9;
}

/* whether vertex c (corner ic of face f, on deg faces) has an edge to d, by the faces round c: 1 or 0, -1 when the
   walk meets an edge without one neighbour (a border, a non-manifold edge) or c has faces past its fan */
static int sx_flat_has_edge(unsigned int const* t, unsigned int const* nb, unsigned int f, unsigned int ic, unsigned int d, unsigned int deg)
{
	unsigned int g = f, j = ic, steps = 0;
	do
	{
		unsigned int e = nb[3 * g + (j + 2) % 3];   /* (across the edge into c) */
		if (t[3 * (size_t)g + (j + 1) % 3] == d || t[3 * (size_t)g + (j + 2) % 3] == d)
			return 1;
		if (e == SX_NONE || ++steps > 4096)
			return -1;
		g = e / 3;
		j = e % 3;   /* (that edge leaves c in g) */
	} while (g != f);
	return steps == deg ? 0 : -1;
}

typedef struct sx_flat_scan
{
	unsigned int const* t;
	qaws_scalar const* src;
	unsigned int const* nb;
	qaws_scalar const* pos;
	double eps;
	unsigned int split;
	int one;
	unsigned char* ok;
} sx_flat_scan;

static qaws_status sx_flat_scan_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_flat_scan const* S = (sx_flat_scan const*)user;
	unsigned int e;
	(void)chunk;
	for (e = begin; e < end; e++)
		S->ok[e] = S->nb[e] != SX_NONE && S->nb[e] > e && sx_flat_flip_ok(S->t, S->src, S->nb, S->pos, e, S->eps, S->split, S->one);
	return QAWS_STATUS_OK;
}

#define SX_FLAT_GRAIN 65536

/* the kept faces' neighbours from the arrangement's edge uses: a use kept is corner 3 k + i of kept face k, i turned
   with the face (its direction too) */
typedef struct sx_flat_pair
{
	unsigned long long const* grp;     /* per use: edge key << 1 | direction, corner; by key */
	unsigned int n;
	unsigned int const* kmap;          /* per piece: 2 k + turned (k its kept face), SX_NONE when dropped */
	unsigned int* nb;
	int* open;                         /* per chunk: an edge whose kept uses do not cancel */
} sx_flat_pair;

static qaws_status sx_flat_pair_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_flat_pair const* P = (sx_flat_pair const*)user;
	unsigned int i = begin, j;
	int open = 0;
	while (i > 0 && i < P->n && P->grp[2 * (size_t)i] >> 1 == P->grp[2 * (size_t)i - 2] >> 1)
		i++;   /* (an edge's uses are the chunk's where they start) */
	for (; i < end; i = j)
	{
		long s = 0;
		unsigned int m = 0, use[2] = { 0, 0 };
		int dir[2] = { 0, 0 };
		for (j = i; j < P->n && P->grp[2 * (size_t)j] >> 1 == P->grp[2 * (size_t)i] >> 1; j++)
		{
			unsigned int c = (unsigned int)P->grp[2 * (size_t)j + 1], f = P->kmap[c / 3];
			int d;
			if (f == SX_NONE)
				continue;
			d = (int)(P->grp[2 * (size_t)j] & 1) ^ (int)(f & 1);
			s += d ? -1 : 1;
			if (m < 2)
			{
				use[m] = 3 * (f >> 1) + ((f & 1) ? 2 - c % 3 : c % 3);
				dir[m] = d;
			}
			m++;
		}
		open |= s != 0;
		if (m == 2 && dir[0] != dir[1])
		{
			P->nb[use[0]] = use[1];
			P->nb[use[1]] = use[0];
		}
	}
	P->open[chunk] = open;
	return QAWS_STATUS_OK;
}

/* triangles' edge uses: the edge's key << 1 | direction, the corner */
typedef struct sx_uses
{
	unsigned int const* tri;
	unsigned long long* grp;
} sx_uses;

static qaws_status sx_uses_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	sx_uses const* U = (sx_uses const*)user;
	unsigned int j;
	(void)chunk;
	for (j = begin; j < end; j++)
	{
		unsigned long long p = U->tri[j], q = U->tri[3 * (j / 3) + (j % 3 + 1) % 3];
		U->grp[2 * (size_t)j] = (p < q ? (p << 32 | q) : (q << 32 | p)) << 1 | (p > q);
		U->grp[2 * (size_t)j + 1] = j;
	}
	return QAWS_STATUS_OK;
}

/* edges inside flat regions turned toward Delaunay: a planar region's pieces were triangulated face by face, its
   input faces' edges kept; the shape does not change. The kept faces t come from the arrangement's pieces through
   kmap, whose edge uses grp (ngrp) give the neighbours nb (3 nt, SX_NONE across an edge not of exactly two faces turned
   opposite ways; as turned at the end). Returns the edges turned, and whether every directed edge is cancelled by one
   turned the other way (closed; turning keeps that). one: only edges between pieces of one source face (attributes
   carried from it are then never extrapolated). */
static unsigned int sx_flatten(unsigned int* t, qaws_scalar* src, unsigned int nt, qaws_scalar const* pos, double eps, unsigned int split,
	int one, unsigned long long const* grp, unsigned int ngrp, unsigned int const* kmap, unsigned int* nb, int* closed, qaws_batch_executor const* executor,
	int* oom)
{
	unsigned int n = 3 * nt, cap = 16, i, top = 0, flips = 0, guard = 10 * nt + 64, nv = 0;
	unsigned int *stack, *deg;
	unsigned long long* set = NULL;
	sx_flat_pair P;
	while (cap < 4 * n)
		cap *= 2;
	for (i = 0; i < n; i++)
		nv = t[i] >= nv ? t[i] + 1 : nv;
	stack = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 4 * (size_t)guard + 1));
	deg = (unsigned int*)qaws_internal_zalloc((size_t)nv + 1, sizeof(unsigned int));   /* (faces on each vertex) */
	P.grp = grp;
	P.n = ngrp;
	P.kmap = kmap;
	P.nb = nb;
	P.open = (int*)qaws_internal_malloc(sizeof(int) * ((size_t)qaws_internal_chunk_count(ngrp, SX_FLAT_GRAIN) + 1));
	if (!stack || !deg || !P.open)
	{
		*oom = 1;
		goto out;
	}
	/* neighbours across edges of exactly two faces turned opposite ways, by the arrangement's edge uses (sorted for
	   the patches) kept */
	for (i = 0; i < n; i++)
	{
		nb[i] = SX_NONE;
		deg[t[i]]++;
	}
	if (qaws_internal_parallel(ngrp, SX_FLAT_GRAIN, sx_flat_pair_job, &P, executor) != QAWS_STATUS_OK)
	{
		*oom = 1;
		goto out;
	}
	*closed = 1;
	for (i = 0; i < qaws_internal_chunk_count(ngrp, SX_FLAT_GRAIN); i++)
		if (P.open[i])
			*closed = 0;
	{
		/* the edges worth turning as given, found in parallel (most meshes have none) */
		sx_flat_scan S;
		S.t = t;
		S.src = src;
		S.nb = nb;
		S.pos = pos;
		S.eps = eps;
		S.split = split;
		S.one = one;
		S.ok = (unsigned char*)qaws_internal_malloc((size_t)n + 1);
		if (!S.ok || qaws_internal_parallel(n, 4096, sx_flat_scan_job, &S, executor) != QAWS_STATUS_OK)
		{
			qaws_internal_free(S.ok);
			*oom = 1;
			goto out;
		}
		for (i = 0; i < n; i++)
			if (S.ok[i])
				stack[top++] = i;
		qaws_internal_free(S.ok);
	}
	while (top && flips < guard)
	{
		unsigned int e = stack[--top], f, g, a, b, c, d, nbc, nca, nad, ndb;
		int k;
		unsigned int ei, ej;
		if (!sx_flat_flip_ok(t, src, nb, pos, e, eps, split, one))
			continue;
		f = e / 3;
		ei = e % 3;
		g = nb[e] / 3;
		ej = nb[e] % 3;
		a = t[3 * (size_t)f + ei];
		b = t[3 * (size_t)f + (ei + 1) % 3];
		c = t[3 * (size_t)f + (ei + 2) % 3];
		d = t[3 * (size_t)g + (ej + 2) % 3];
		k = sx_flat_has_edge(t, nb, f, (ei + 2) % 3, d, deg[c]);
		if (k == 1)
			continue;   /* (the edge c d is there already) */
		if (k < 0 && !set)
		{
			/* the mesh's edges, made when first needed (the walk round c rarely meets a border) */
			set = (unsigned long long*)qaws_internal_zalloc(cap, sizeof(unsigned long long));
			if (!set)
			{
				*oom = 1;
				goto out;
			}
			for (i = 0; i < n; i++)
			{
				unsigned long long ek = sx_ekey(t[i], t[3 * (i / 3) + (i % 3 + 1) % 3]);
				*sx_eset_slot(set, cap, ek, 1) = ek;
			}
		}
		if (k < 0 && sx_eset_slot(set, cap, sx_ekey(c, d), 0))
			continue;
		nbc = nb[3 * f + (ei + 1) % 3];
		nca = nb[3 * f + (ei + 2) % 3];
		nad = nb[3 * g + (ej + 1) % 3];
		ndb = nb[3 * g + (ej + 2) % 3];
		/* f = a d c, g = d b c */
		t[3 * (size_t)f] = a;
		t[3 * (size_t)f + 1] = d;
		t[3 * (size_t)f + 2] = c;
		t[3 * (size_t)g] = d;
		t[3 * (size_t)g + 1] = b;
		t[3 * (size_t)g + 2] = c;
		nb[3 * f] = nad;
		nb[3 * f + 1] = 3 * g + 2;
		nb[3 * f + 2] = nca;
		nb[3 * g] = ndb;
		nb[3 * g + 1] = nbc;
		nb[3 * g + 2] = 3 * f + 1;
		if (nad != SX_NONE) nb[nad] = 3 * f;
		if (nca != SX_NONE) nb[nca] = 3 * f + 2;
		if (ndb != SX_NONE) nb[ndb] = 3 * g;
		if (nbc != SX_NONE) nb[nbc] = 3 * g + 1;
		deg[a]--;
		deg[b]--;
		deg[c]++;
		deg[d]++;
		if (set)
		{
			*sx_eset_slot(set, cap, sx_ekey(a, b), 0) = ~0ull;
			*sx_eset_slot(set, cap, sx_ekey(c, d), 1) = sx_ekey(c, d);
		}
		stack[top++] = 3 * f;
		stack[top++] = 3 * f + 2;
		stack[top++] = 3 * g;
		stack[top++] = 3 * g + 1;
		flips++;
	}
out:
	qaws_internal_free(deg);
	qaws_internal_free(P.open);
	qaws_internal_free(set);
	qaws_internal_free(stack);
	return flips;
}

qaws_status qaws_mesh_boolean(qaws_mesh_3d** out, qaws_boolean_report* report, qaws_mesh_3d_desc const* a, qaws_mesh_3d_desc const* b,
	qaws_boolean_op op, qaws_intersect_desc const* desc)
{
	return sx_boolean(out, report, a, b, op, desc, 0, NULL);
}

/* the pieces round edge uses g (m of them, sorted as sx_uses_job keys them): per use its patch (z), its side toward the
   turn (s: +1 when it runs the edge from its smaller vertex, its front then faces the turn), and their order round the edge
   (ord); 0 when two lie at one angle (coincident: the order is not told) */
static int sx_round(sx_out const* o, unsigned long long const* g, unsigned int m, unsigned int const* pid, unsigned int* ord, unsigned int* z, int* s)
{
	unsigned int v0 = (unsigned int)(g[0] >> 33), v1 = (unsigned int)((g[0] >> 1) & 0xFFFFFFFFu), q, r, k;
	double d[3], e1[3], e2[3], an[16], l = 0;
	for (k = 0; k < 3; k++)
	{
		d[k] = (double)o->pos[3 * (size_t)v1 + k] - (double)o->pos[3 * (size_t)v0 + k];
		l += d[k] * d[k];
	}
	if (!(l > 0))
		return 0;
	for (k = 0; k < 3; k++)
		d[k] /= sqrt(l);
	for (q = 0; q < m; q++)
	{
		unsigned int c = (unsigned int)g[2 * (size_t)q + 1], t = c / 3, w = o->tri[3 * (size_t)t + (c % 3 + 2) % 3];
		double rr[3], dd = 0, rl = 0;
		for (k = 0; k < 3; k++)
		{
			rr[k] = (double)o->pos[3 * (size_t)w + k] - (double)o->pos[3 * (size_t)v0 + k];
			dd += rr[k] * d[k];
		}
		for (k = 0; k < 3; k++)
		{
			rr[k] -= dd * d[k];
			rl += rr[k] * rr[k];
		}
		if (!(rl > 0))
			return 0;
		if (q == 0)
		{
			for (k = 0; k < 3; k++)
				e1[k] = rr[k] / sqrt(rl);
			sx_cross(d, e1, e2);
		}
		an[q] = atan2(sx_dot(rr, e2), sx_dot(rr, e1));
		z[q] = pid[t];
		s[q] = (g[2 * (size_t)q] & 1) ? -1 : 1;
		ord[q] = q;
	}
	for (q = 1; q < m; q++)
		for (r = q; r > 0 && an[ord[r - 1]] > an[ord[r]]; r--)
		{
			unsigned int x = ord[r];
			ord[r] = ord[r - 1];
			ord[r - 1] = x;
		}
	for (q = 0; q + 1 < m; q++)
		if (an[ord[q + 1]] - an[ord[q]] < 1e-12)
			return 0;
	return 1;
}

/* Windings sampled a hair off a patch can land past a surface nearer still (pieces crossing at a shallow angle). Round an
   edge of three or more pieces, in turn by angle, two neighbours share the winding numbers of the gap between them, and
   crossing a piece changes its own operand's by one (its back the inside). A small patch (its fattest piece's inradius
   under 1e-3 of the box) whose neighbours in turn agree with each other across it is distrusted when its own samples
   differ by other than that one, or when it disagrees with both of them (edges by coincident pieces left alone), then each distrusted one takes its windings from a trusted neighbour across such an edge, repeatedly: a
   patch sampled right never changes. wa, wb: per patch its front and back windings (wb NULL for a self-union); nfa: A's
   faces. Returns the patches mended; -1 on allocation failure. */
static int sx_mend(sx_out const* o, unsigned long long const* grp, unsigned int const* pid, double const* fat, unsigned int const* rep,
	unsigned int npatch, qaws_scalar* wa, qaws_scalar* wb, unsigned int nfa)
{
	unsigned int n = 3 * o->nt, i, j, q, pass, mended = 0;
	unsigned char* bad = (unsigned char*)qaws_internal_zalloc((size_t)npatch + 1, 1);
	if (!bad)
		return -1;
	for (pass = 0; pass < 32; pass++)
	{
		int changed = 0;
		for (i = 0; i < n; i = j)
		{
			unsigned int m, ord[16], z[16];
			int s[16];
			for (j = i + 1; j < n && grp[2 * (size_t)j] >> 1 == grp[2 * (size_t)i] >> 1; j++)
				;
			m = j - i;
			if (m < 3 || m > 16 || !sx_round(o, grp + 2 * (size_t)i, m, pid, ord, z, s))
				continue;
			for (q = 0; q < m; q++)
				if (o->copl[(unsigned int)o->src[(unsigned int)grp[2 * ((size_t)i + q) + 1] / 3]])
					break;
			if (q < m)
				continue;   /* (coincident pieces: a winding may jump by more than one there) */
			if (pass == 0)
			{
				/* a patch disagreeing with both its neighbours in turn, which agree with each other across it: sampled wrong */
				for (q = 0; q < m; q++)
				{
					unsigned int pc = ord[q], pn = ord[(q + 1) % m], pp = ord[(q + m - 1) % m], P = z[pc], k2;
					unsigned int iaft = s[pc] > 0 ? 0 : 1, ibef = 1 - iaft, inb = s[pn] > 0 ? 1 : 0, ipb = s[pp] > 0 ? 0 : 1;
					unsigned int src = (unsigned int)o->src[(unsigned int)grp[2 * ((size_t)i + pc) + 1] / 3];
					int ok2 = z[pn] != P && z[pp] != P && fat[rep[P]] < 2e-3 * o->diag, self = 0;
					/* (its own two samples apart by other than its piece's one: sampled wrong for certain) */
					for (k2 = 0; k2 < (wb ? 2u : 1u); k2++)
					{
						qaws_scalar const* w = k2 ? wb : wa;
						double jump = floor((double)w[2 * (size_t)P + 1] + 0.5) - floor((double)w[2 * (size_t)P] + 0.5);
						self |= jump != ((k2 ? src >= nfa : src < nfa) ? 1 : 0);
					}
					for (k2 = 0; k2 < (wb ? 2u : 1u) && ok2; k2++)
					{
						qaws_scalar const* w = k2 ? wb : wa;
						double aft = floor((double)w[2 * (size_t)z[pn] + inb] + 0.5), bef = floor((double)w[2 * (size_t)z[pp] + ipb] + 0.5);
						double own = (k2 ? src >= nfa : src < nfa) ? (iaft == 0 ? 1 : -1) : 0;
						/* (the neighbours' gaps across P: its back one more than its front for its operand) */
						ok2 = bef == aft + own;
					}
					if (!ok2)
						continue;
					if (self)
						bad[P] = 1;
					for (k2 = 0; k2 < (wb ? 2u : 1u); k2++)
					{
						qaws_scalar const* w = k2 ? wb : wa;
						if (floor((double)w[2 * (size_t)P + iaft] + 0.5) != floor((double)w[2 * (size_t)z[pn] + inb] + 0.5)
							&& floor((double)w[2 * (size_t)P + ibef] + 0.5) != floor((double)w[2 * (size_t)z[pp] + ipb] + 0.5))
							bad[P] = 1;
					}
				}
				continue;
			}
			for (q = 0; q < m; q++)
			{
				/* P's side toward the turn and Q's side back from it bound one gap */
				unsigned int pa = ord[q], qb = ord[(q + 1) % m], P = z[pa], Q = z[qb];
				unsigned int ip = s[pa] > 0 ? 0 : 1, iq = s[qb] > 0 ? 1 : 0;
				int same = floor((double)wa[2 * (size_t)P + ip] + 0.5) == floor((double)wa[2 * (size_t)Q + iq] + 0.5)
					&& (!wb || floor((double)wb[2 * (size_t)P + ip] + 0.5) == floor((double)wb[2 * (size_t)Q + iq] + 0.5));
				unsigned int to, from, it, ifr, src;
				double va, vb;
				if (P == Q || same)
					continue;
				if (bad[P] == bad[Q])
					continue;
				to = bad[P] ? P : Q;
				from = bad[P] ? Q : P;
				it = bad[P] ? ip : iq;
				ifr = bad[P] ? iq : ip;
				va = floor((double)wa[2 * (size_t)from + ifr] + 0.5);
				vb = wb ? floor((double)wb[2 * (size_t)from + ifr] + 0.5) : 0;
				/* the other side across the patch's own piece: its operand's one more behind, one less in front */
				src = (unsigned int)o->src[rep[to]];
				wa[2 * (size_t)to + it] = (qaws_scalar)va;
				wa[2 * (size_t)to + 1 - it] = (qaws_scalar)(va + (src < nfa ? (it == 0 ? 1 : -1) : 0));
				if (wb)
				{
					wb[2 * (size_t)to + it] = (qaws_scalar)vb;
					wb[2 * (size_t)to + 1 - it] = (qaws_scalar)(vb + (src >= nfa ? (it == 0 ? 1 : -1) : 0));
				}
				bad[to] = 0;
				mended++;
				changed = 1;
			}
		}
		if (pass > 0 && !changed)
			break;
	}
	qaws_internal_free(bad);
	return (int)mended;
}

/* the opposites of the kept faces' corners on edges of four or more of them (solids touching along an edge): round the
   edge in turn by angle, each face's back (the inside) faces its neighbour's on that side, which bounds the same lobe of the
   result; those two are paired, so each lobe closes on its own. t: the kept faces (outward), opp: per corner its opposite
   (SX_NONE where not known; filled here where the order tells). 0 on allocation failure. */
static int sx_pair_round(unsigned int const* t, qaws_scalar const* pos, unsigned int nt, unsigned int* opp)
{
	unsigned int n = 3 * nt, i, j, q;
	unsigned long long* g = (unsigned long long*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned long long), 2, (size_t)n + 1));
	if (!g)
		return 0;
	for (i = 0; i < n; i++)
	{
		unsigned long long a = t[i], b = t[3 * (i / 3) + (i % 3 + 1) % 3];
		g[2 * (size_t)i] = (a < b ? (a << 32 | b) : (b << 32 | a)) << 1 | (a > b);
		g[2 * (size_t)i + 1] = i;
	}
	qsort(g, n, 2 * sizeof(unsigned long long), sx_cmp_edge);
	for (i = 0; i < n; i = j)
	{
		unsigned int m, ord[16], v0, v1, r, k, missing = 0;
		double d[3], e1[3], e2[3], an[16], l = 0;
		int ok = 1;
		for (j = i + 1; j < n && g[2 * (size_t)j] >> 1 == g[2 * (size_t)i] >> 1; j++)
			;
		m = j - i;
		for (q = i; q < j; q++)
			missing += opp[g[2 * (size_t)q + 1]] == SX_NONE;
		if (m < 4 || m > 16 || (m & 1) || !missing)
			continue;
		v0 = (unsigned int)(g[2 * (size_t)i] >> 33);
		v1 = (unsigned int)((g[2 * (size_t)i] >> 1) & 0xFFFFFFFFu);
		for (k = 0; k < 3; k++)
		{
			d[k] = (double)pos[3 * (size_t)v1 + k] - (double)pos[3 * (size_t)v0 + k];
			l += d[k] * d[k];
		}
		if (!(l > 0))
			continue;
		for (k = 0; k < 3; k++)
			d[k] /= sqrt(l);
		for (q = 0; q < m && ok; q++)
		{
			unsigned int c = (unsigned int)g[2 * ((size_t)i + q) + 1], w = t[3 * (c / 3) + (c % 3 + 2) % 3];
			double rr[3], dd = 0, rl = 0;
			for (k = 0; k < 3; k++)
			{
				rr[k] = (double)pos[3 * (size_t)w + k] - (double)pos[3 * (size_t)v0 + k];
				dd += rr[k] * d[k];
			}
			for (k = 0; k < 3; k++)
			{
				rr[k] -= dd * d[k];
				rl += rr[k] * rr[k];
			}
			if (!(rl > 0))
			{
				ok = 0;
				break;
			}
			if (q == 0)
			{
				for (k = 0; k < 3; k++)
					e1[k] = rr[k] / sqrt(rl);
				sx_cross(d, e1, e2);
			}
			an[q] = atan2(sx_dot(rr, e2), sx_dot(rr, e1));
			ord[q] = q;
		}
		if (!ok)
			continue;
		for (q = 1; q < m; q++)
			for (r = q; r > 0 && an[ord[r - 1]] > an[ord[r]]; r--)
			{
				unsigned int x = ord[r];
				ord[r] = ord[r - 1];
				ord[r - 1] = x;
			}
		for (q = 0; q + 1 < m && ok; q++)
			ok = an[ord[q + 1]] - an[ord[q]] > 1e-12;
		/* a face running the edge from v0 has its front toward the turn, its back toward the one before it */
		for (q = 0; q < m && ok; q++)
		{
			/* (by position in turn: the partner of position q, and the partner's own is q) */
			unsigned int pq = (g[2 * ((size_t)i + ord[q])] & 1) ? (q + 1) % m : (q + m - 1) % m;
			unsigned int pp = (g[2 * ((size_t)i + ord[pq])] & 1) ? (pq + 1) % m : (pq + m - 1) % m;
			ok = pp == q && (g[2 * ((size_t)i + ord[q])] & 1) != (g[2 * ((size_t)i + ord[pq])] & 1);
		}
		for (q = 0; q < m && ok; q++)
		{
			unsigned int a = ord[q], b = (g[2 * ((size_t)i + a)] & 1) ? ord[(q + 1) % m] : ord[(q + m - 1) % m];
			opp[g[2 * ((size_t)i + a) + 1]] = (unsigned int)g[2 * ((size_t)i + b) + 1];
		}
	}
	qaws_internal_free(g);
	return 1;
}

qaws_status sx_boolean(qaws_mesh_3d** out, qaws_boolean_report* report, qaws_mesh_3d_desc const* a, qaws_mesh_3d_desc const* b,
	qaws_boolean_op op, qaws_intersect_desc const* desc, int one, unsigned char** out_turned)
{
	qaws_boolean_report local;
	qaws_scalar *x = NULL, *kpos = NULL, *keep_src = NULL, *pts = NULL, *wa = NULL, *wb = NULL;
	unsigned int *fv = NULL, *fs = NULL, *ta = NULL, *tb = NULL, *keep = NULL, *used = NULL;
	qaws_mesh_bvh *ba = NULL, *bb = NULL;
	qaws_mesh_hit* hits = NULL;
	unsigned int *par = NULL, *rep = NULL, *pid = NULL, *kmap = NULL, *kopp = NULL, npatch = 0;
	unsigned long long* grp = NULL;
	double *fat = NULL, *mid = NULL;
	unsigned char *near = NULL, *dup = NULL;
	sx_out o;
	qaws_status st;
	unsigned int na = 0, nb = 0, nv, nf, i, j, k, nk = 0, two = b != NULL;
	memset(&o, 0, sizeof o);
	if (out_turned)
		*out_turned = NULL;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!a || !a->face_count || !a->face_vertices || !a->positions || (b && (!b->face_count || !b->face_vertices || !b->positions))
		|| (unsigned int)op > (unsigned int)QAWS_BOOLEAN_DIFFERENCE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* the operands as one soup: B's vertices and faces after A's */
	for (i = 0; i < a->face_count; i++)
		na += a->face_sizes ? a->face_sizes[i] : 3;
	for (i = 0; two && i < b->face_count; i++)
		nb += b->face_sizes ? b->face_sizes[i] : 3;
	nv = a->vertex_count + (two ? b->vertex_count : 0);
	nf = a->face_count + (two ? b->face_count : 0);
	x = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)na + nb + 1));
	fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	if (!x || !fv || !fs)
		goto oom;
	memcpy(x, a->positions, sizeof(qaws_scalar) * 3 * (size_t)a->vertex_count);
	for (i = 0; i < na; i++)
		fv[i] = a->face_vertices[i] < a->vertex_count ? a->face_vertices[i] : nv;   /* (out of range: refused by the resolution) */
	for (i = 0; i < a->face_count; i++)
		fs[i] = a->face_sizes ? a->face_sizes[i] : 3;
	if (two)
	{
		memcpy(x + 3 * (size_t)a->vertex_count, b->positions, sizeof(qaws_scalar) * 3 * (size_t)b->vertex_count);
		for (i = 0; i < nb; i++)
			fv[na + i] = b->face_vertices[i] < b->vertex_count ? b->face_vertices[i] + a->vertex_count : nv;
		for (i = 0; i < b->face_count; i++)
			fs[a->face_count + i] = b->face_sizes ? b->face_sizes[i] : 3;
	}
	st = sx_resolve(&o, &report->intersect, x, fv, fs, nv, nf, desc);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* patches: pieces joined across edges of two faces turned opposite ways are on the same sides of everything */
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)o.nt + 1));
	pid = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)o.nt + 1));
	rep = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)o.nt + 1));
	fat = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)o.nt + 1));
	grp = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 6, (size_t)o.nt + 1));
	kmap = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)o.nt + 1));
	keep = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)o.nt + 1));
	keep_src = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)o.nt + 1));
	used = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)o.nv + 1));
	kpos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)o.nv + 1));
	if (!par || !pid || !rep || !fat || !grp || !kmap || !keep || !keep_src || !used || !kpos)
		goto oom;
	{
		/* the edge uses sorted with their corners (kept for the flat regions): key the edge, then the direction bit */
		unsigned int n = 3 * o.nt;
		sx_uses U;
		U.tri = o.tri;
		U.grp = grp;
		if (qaws_internal_parallel(n, SX_FLAT_GRAIN, sx_uses_job, &U, desc ? desc->executor : NULL) != QAWS_STATUS_OK)
			goto oom;
		if (!sx_radix(grp, n, 2, desc ? desc->executor : NULL))
			qsort(grp, n, 2 * sizeof(unsigned long long), sx_cmp_edge);
		for (i = 0; i < o.nt; i++)
			par[i] = i;
		for (i = 0; i < n; i = j)
		{
			for (j = i + 1; j < n && grp[2 * (size_t)j] >> 1 == grp[2 * (size_t)i] >> 1; j++)
				;
			if (j == i + 2 && (grp[2 * (size_t)i] & 1) != (grp[2 * (size_t)i + 2] & 1))
			{
				unsigned int x0 = sx_root(par, (unsigned int)grp[2 * (size_t)i + 1] / 3), x1 = sx_root(par, (unsigned int)grp[2 * (size_t)i + 3] / 3);
				if (x0 != x1)
					par[x0 < x1 ? x1 : x0] = x0 < x1 ? x0 : x1;
			}
		}
	}
	/* each patch's fattest piece stands for it */
	for (i = 0; i < o.nt; i++)
	{
		double P[3][3], u[3], v[3], n[3], per = 0, l;
		unsigned int z = sx_root(par, i);
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				P[k][j] = (double)o.pos[3 * (size_t)o.tri[3 * (size_t)i + k] + j];
		for (k = 0; k < 3; k++)
		{
			u[k] = P[1][k] - P[0][k];
			v[k] = P[2][k] - P[0][k];
			per += sqrt(sx_dist2(P[k], P[(k + 1) % 3]));
		}
		sx_cross(u, v, n);
		l = sqrt(sx_dot(n, n));
		fat[i] = per > 0 ? l / per : 0;   /* twice the inradius */
		if (z == i)
		{
			pid[i] = npatch++;
			rep[pid[i]] = i;
		}
	}
	for (i = 0; i < o.nt; i++)
	{
		unsigned int z = pid[sx_root(par, i)];
		pid[i] = z;
		if (fat[i] > fat[rep[z]])
			rep[z] = i;
	}
	report->patches = npatch;
	/* each patch sampled just in front of and behind its piece */
	pts = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 6, (size_t)npatch + 1));
	wa = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 2 * (size_t)npatch + 1));
	wb = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 2 * (size_t)npatch + 1));
	hits = (qaws_mesh_hit*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_hit), (size_t)npatch + 1));
	if (!pts || !wa || !wb || !hits)
		goto oom;
	for (i = 0; i < npatch; i++)
	{
		unsigned int t = rep[i];
		double P[3][3], u[3], v[3], n[3], g[3], l, dl;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				P[k][j] = (double)o.pos[3 * (size_t)o.tri[3 * (size_t)t + k] + j];
		for (k = 0; k < 3; k++)
		{
			u[k] = P[1][k] - P[0][k];
			v[k] = P[2][k] - P[0][k];
			g[k] = (P[0][k] + P[1][k] + P[2][k]) / 3;
		}
		sx_cross(u, v, n);
		l = sqrt(sx_dot(n, n));
		/* off the piece by a tenth of its inradius, at most 1e-6 D (nearer than any other surface, mostly: see below) */
		dl = l > 0 ? fmin(0.1 * fat[t], 1e-6 * o.diag) / l : 0;
		for (k = 0; k < 3; k++)
		{
			pts[6 * (size_t)i + k] = (qaws_scalar)(g[k] + dl * n[k]);
			pts[6 * (size_t)i + 3 + k] = (qaws_scalar)(g[k] - dl * n[k]);
		}
	}
	/* sheets nearer than that (overlapping sweeps lay walls within 1e-5 D of one another): each sample brought within 0.4 of
	   the first surface its ray meets (from 1e-9 D off the piece, so not the piece itself) */
	{
		unsigned int rop;
		for (rop = 0; rop < (two ? 2u : 1u); rop++)
		{
			qaws_scalar* org = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 6, (size_t)npatch + 1));
			qaws_scalar* dir = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 6, (size_t)npatch + 1));
			qaws_mesh_hit* rh = (qaws_mesh_hit*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_hit), 2 * (size_t)npatch + 1));
			double e0 = 1e-9 * o.diag, far = 1e-6 * o.diag;
			/* (the operand's BVH, kept for its winding numbers below) */
			if (!org || !dir || !rh || sx_bvh(rop ? &bb : &ba, rop ? b : a, rop ? &tb : &ta) != QAWS_STATUS_OK)
			{
				qaws_internal_free(org); qaws_internal_free(dir); qaws_internal_free(rh);
				goto oom;
			}
			for (i = 0; i < 2 * npatch; i++)
			{
				unsigned int pi = i / 2, s = i % 2;
				double g[3], d[3], l2 = 0;
				for (k = 0; k < 3; k++)
				{
					g[k] = 0.5 * ((double)pts[6 * (size_t)pi + k] + (double)pts[6 * (size_t)pi + 3 + k]);
					d[k] = (double)pts[6 * (size_t)pi + 3 * s + k] - g[k];
					l2 += d[k] * d[k];
				}
				l2 = sqrt(l2);
				for (k = 0; k < 3; k++)
				{
					d[k] = l2 > 0 ? d[k] / l2 : 0;
					org[3 * (size_t)i + k] = (qaws_scalar)(g[k] + e0 * d[k]);
					dir[3 * (size_t)i + k] = (qaws_scalar)d[k];
				}
			}
			st = qaws_mesh_bvh_ray(rh, rop ? bb : ba, org, dir, 2 * npatch, (qaws_scalar)far, 0, desc ? desc->executor : NULL);
			for (i = 0; st == QAWS_STATUS_OK && i < 2 * npatch; i++)
				if (rh[i].element != QAWS_MESH_NONE && sx_offset_of(pts, i) > 0)
				{
					unsigned int pi = i / 2, s = i % 2;
					double g[3], d[3], lo = sx_offset_of(pts, i), want = 0.4 * ((double)rh[i].distance + e0);
					if (want >= lo)
						continue;
					for (k = 0; k < 3; k++)
					{
						g[k] = 0.5 * ((double)pts[6 * (size_t)pi + k] + (double)pts[6 * (size_t)pi + 3 + k]);
						d[k] = ((double)pts[6 * (size_t)pi + 3 * s + k] - g[k]) / lo;
					}
					/* (both samples brought in alike, so they stay about the piece's centre) */
					for (k = 0; k < 3; k++)
					{
						pts[6 * (size_t)pi + 3 * s + k] = (qaws_scalar)(g[k] + want * d[k]);
						pts[6 * (size_t)pi + 3 * (1 - s) + k] = (qaws_scalar)(g[k] - want * d[k]);
					}
				}
			qaws_internal_free(org); qaws_internal_free(dir); qaws_internal_free(rh);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
	}
	/* the operands' winding numbers there, and B's patches lying on A (coplanar overlaps: A's pieces stand for both) */
	{
		qaws_batch_executor const* ex = desc ? desc->executor : NULL;
		unsigned int nta = 0, ntb = 0;
		double work;
		for (i = 0; i < a->face_count; i++)
			nta += (a->face_sizes ? a->face_sizes[i] : 3) - 2;
		for (i = 0; two && i < b->face_count; i++)
			ntb += (b->face_sizes ? b->face_sizes[i] : 3) - 2;
		work = 2.0 * npatch * ((double)nta + ntb) + (double)npatch * nta;
		near = (unsigned char*)qaws_internal_zalloc((size_t)npatch + 1, 1);
		mid = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)npatch + 1));
		if (!near || !mid)
			goto oom;
		for (i = 0; i < 6 * npatch; i++)
			mid[i] = (double)pts[i];
		if (work <= 2e6)
		{
			/* few patches: the solid angles summed over the triangles, no BVH built (past about 2e6 solid angles the BVH's
			   queries are cheaper: a rounded cube's 2,330 patches 250 ms to 80) */
			st = QAWS_STATUS_OK;   /* (ta and tb made with the BVHs above) */
			if (st == QAWS_STATUS_OK)
				st = sx_direct_winding(wa, NULL, a->positions, ta, nta, mid, 2 * npatch, 0, ex);
			if (st == QAWS_STATUS_OK && two)
				st = sx_direct_winding(wb, NULL, b->positions, tb, ntb, mid, 2 * npatch, 0, ex);
			if (st == QAWS_STATUS_OK && two)
			{
				for (i = 0; i < npatch; i++)
					for (k = 0; k < 3; k++)
						mid[3 * (size_t)i + k] = 0.5 * (mid[6 * (size_t)i + k] + mid[6 * (size_t)i + 3 + k]);
				st = sx_direct_winding(NULL, near, a->positions, ta, nta, mid, npatch, 4 * o.eps, ex);
			}
		}
		else
		{
			st = QAWS_STATUS_OK;   /* (the BVHs made above) */
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_bvh_winding_number(wa, ba, pts, 2 * npatch, ex);
			if (st == QAWS_STATUS_OK && two)
				st = qaws_mesh_bvh_winding_number(wb, bb, pts, 2 * npatch, ex);
			if (st == QAWS_STATUS_OK && two)
			{
				for (i = 0; i < npatch; i++)
					for (k = 0; k < 3; k++)
						pts[3 * (size_t)i + k] = (qaws_scalar)(0.5 * ((double)pts[6 * (size_t)i + k] + (double)pts[6 * (size_t)i + 3 + k]));
				st = qaws_mesh_bvh_closest(hits, ba, pts, npatch, (qaws_scalar)(4 * o.eps), ex);
				for (i = 0; st == QAWS_STATUS_OK && i < npatch; i++)
					near[i] = hits[i].element != QAWS_MESH_NONE;
			}
		}
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	if (sx_mend(&o, grp, pid, fat, rep, npatch, wa, two ? wb : NULL, a->face_count) < 0)
		goto oom;
	/* coincident pieces (an operand overlapping itself in a plane) kept once: identical pieces from the lowest input face, and
	   in a self-union a patch lying on a lower input face of its plane (either way round) dropped: the lowest face's pieces stand
	   for the overlap, their windings deciding which way it faces */
	dup = (unsigned char*)qaws_internal_malloc((size_t)o.nt + (size_t)npatch + 1);
	if (!dup || !sx_dup_mark(dup, o.tri, o.src, o.nt))
		goto oom;
	memset(dup + o.nt, 0, npatch);
	if (!two && ta)
	{
		unsigned int f, tri = 0, nta = 0, *tf;
		for (f = 0; f < a->face_count; f++)
			nta += (a->face_sizes ? a->face_sizes[f] : 3) - 2;
		tf = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nta + 1));
		if (!tf)
			goto oom;
		for (f = 0; f < a->face_count; f++)
			for (k = 2; k < (a->face_sizes ? a->face_sizes[f] : 3); k++)
				tf[tri++] = f;
		for (i = 0; i < npatch; i++)
		{
			unsigned int t = rep[i], src = (unsigned int)o.src[t], q;
			double P[3][3], u[3], v[3], n[3], g[3];
			if (!o.copl[src])
				continue;
			for (k = 0; k < 3; k++)
				for (j = 0; j < 3; j++)
					P[k][j] = (double)o.pos[3 * (size_t)o.tri[3 * (size_t)t + k] + j];
			for (k = 0; k < 3; k++)
			{
				u[k] = P[1][k] - P[0][k];
				v[k] = P[2][k] - P[0][k];
				g[k] = (P[0][k] + P[1][k] + P[2][k]) / 3;
			}
			sx_cross(u, v, n);
			for (q = 0; q < nta && !dup[o.nt + i]; q++)
			{
				double Q[3][3], e1[3], e2[3], m[3], w[3], mm, d, b1, b2, d11, d12, d22, w1, w2, den;
				if (tf[q] >= src || !o.copl[tf[q]])
					continue;
				for (k = 0; k < 3; k++)
					for (j = 0; j < 3; j++)
						Q[k][j] = (double)a->positions[3 * (size_t)ta[3 * (size_t)q + k] + j];
				for (k = 0; k < 3; k++)
				{
					e1[k] = Q[1][k] - Q[0][k];
					e2[k] = Q[2][k] - Q[0][k];
					w[k] = g[k] - Q[0][k];
				}
				sx_cross(e1, e2, m);
				mm = sqrt(sx_dot(m, m));
				if (!(mm > 0) || fabs(sx_dot(m, n)) < 0.999 * mm * sqrt(sx_dot(n, n)))
					continue;
				d = sx_dot(w, m) / mm;
				if (fabs(d) > 4 * o.eps)
					continue;
				d11 = sx_dot(e1, e1); d12 = sx_dot(e1, e2); d22 = sx_dot(e2, e2);
				w1 = sx_dot(w, e1); w2 = sx_dot(w, e2);
				den = d11 * d22 - d12 * d12;
				b1 = (w1 * d22 - w2 * d12) / den;
				b2 = (w2 * d11 - w1 * d12) / den;
				if (b1 > -1e-9 && b2 > -1e-9 && b1 + b2 < 1 + 1e-9)
					dup[o.nt + i] = 1;
			}
		}
		qaws_internal_free(tf);
	}
	for (i = 0; i < o.nt; i++)
	{
		unsigned int z = pid[i];
		int fa = wa[2 * (size_t)z] >= 0.5, ka = wa[2 * (size_t)z + 1] >= 0.5;
		int fb = two && wb[2 * (size_t)z] >= 0.5, kb = two && wb[2 * (size_t)z + 1] >= 0.5;
		int front = sx_region(op, fa, fb, (int)two), back = sx_region(op, ka, kb, (int)two);
		unsigned int s = (unsigned int)o.src[i];
		if (front == back || dup[i] || dup[o.nt + z] || (two && s >= a->face_count && near[z]))
		{
			kmap[i] = SX_NONE;
			report->dropped++;
			continue;
		}
		/* the region behind: turned when it is in front */
		kmap[i] = nk << 1 | (unsigned int)front;
		keep[3 * (size_t)nk] = o.tri[3 * (size_t)i];
		keep[3 * (size_t)nk + 1] = o.tri[3 * (size_t)i + (front ? 2 : 1)];
		keep[3 * (size_t)nk + 2] = o.tri[3 * (size_t)i + (front ? 1 : 2)];
		keep_src[nk++] = o.src[i];
		report->flipped += (unsigned int)front;
	}
	report->kept = nk;
	if (out_turned)
	{
		*out_turned = (unsigned char*)qaws_internal_malloc((size_t)nk + 1);
		if (!*out_turned)
			goto oom;
		for (i = 0; i < o.nt; i++)
			if (kmap[i] != SX_NONE)
				(*out_turned)[kmap[i] >> 1] = (unsigned char)(kmap[i] & 1);
	}
	/* the vertices kept, compacted */
	{
		unsigned int m = 0;
		for (i = 0; i < o.nv; i++)
			used[i] = SX_NONE;
		for (i = 0; i < 3 * nk; i++)
		{
			if (used[keep[i]] == SX_NONE)
			{
				memcpy(kpos + 3 * (size_t)m, o.pos + 3 * (size_t)keep[i], 3 * sizeof(qaws_scalar));
				used[keep[i]] = m++;
			}
			keep[i] = used[keep[i]];
		}
		{
			int bad = 0;
			kopp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nk + 1));
			if (!kopp)
				goto oom;
			report->flat_flips = sx_flatten(keep, keep_src, nk, kpos, o.eps, a->face_count, one, grp, 3 * o.nt, kmap, kopp, &report->closed,
				desc ? desc->executor : NULL, &bad);
			if (bad)
				goto oom;
		}
		/* every edge of exactly two kept faces turned opposite ways: they are the mesh's pairs; edges of more (solids touching
		   along them) paired lobe by lobe */
		if (!sx_pair_round(keep, kpos, nk, kopp))
			goto oom;
		for (i = 0; i < 3 * nk && kopp[i] != SX_NONE; i++) {}
		st = nk ? sx_mesh(out, kpos, m, keep, i == 3 * nk ? kopp : NULL, keep_src, nk) : QAWS_STATUS_OK;
		if (st == QAWS_STATUS_OK && nk && *out && two && bb)
		{
			/* A's faces lying on B (coplanar overlaps: they stand for both, B's there dropped): those of an input face in a
			   coplanar pair whose centre is within 4 eps of B, the "overlap" face attribute 1 there (none without any) */
			unsigned int nq = 0, f2;
			unsigned int* qi = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nk + 1));
			qaws_scalar* qp = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nk + 1));
			qaws_scalar* ov = (qaws_scalar*)qaws_internal_scratch_zalloc((size_t)nk + 1, sizeof(qaws_scalar));
			qaws_mesh_hit* qh = NULL;
			if (!qi || !qp || !ov)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (f2 = 0; st == QAWS_STATUS_OK && f2 < nk; f2++)
			{
				unsigned int s2 = (unsigned int)keep_src[f2];
				if (s2 >= a->face_count || !o.copl[s2])
					continue;
				for (k = 0; k < 3; k++)
					qp[3 * (size_t)nq + k] = (qaws_scalar)(((double)kpos[3 * (size_t)keep[3 * (size_t)f2] + k] + (double)kpos[3 * (size_t)keep[3 * (size_t)f2 + 1] + k] +
						(double)kpos[3 * (size_t)keep[3 * (size_t)f2 + 2] + k]) / 3);
				qi[nq++] = f2;
			}
			if (st == QAWS_STATUS_OK && nq)
			{
				if (!(qh = (qaws_mesh_hit*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_mesh_hit), (size_t)nq))))
					st = QAWS_STATUS_ALLOCATION_FAILURE;
				else
					st = qaws_mesh_bvh_closest(qh, bb, qp, nq, (qaws_scalar)(4 * o.eps), desc ? desc->executor : NULL);
				for (f2 = 0; st == QAWS_STATUS_OK && f2 < nq; f2++)
					if (qh[f2].element != QAWS_MESH_NONE)
					{
						ov[qi[f2]] = 1;
						report->overlapped++;
					}
			}
			if (st == QAWS_STATUS_OK && report->overlapped)
			{
				qaws_attribute_desc ad;
				unsigned int attr;
				memset(&ad, 0, sizeof ad);
				ad.name = "overlap";
				ad.element = QAWS_MESH_FACE;
				ad.components = 1;
				ad.rule = QAWS_ATTRIBUTE_NEAREST;
				st = qaws_mesh_3d_add_attribute(&attr, *out, ov, &ad);
			}
			qaws_internal_free(qi); qaws_internal_free(qp); qaws_internal_free(ov); qaws_internal_free(qh);
			if (st != QAWS_STATUS_OK)
			{
				qaws_mesh_3d_destroy(*out);
				*out = NULL;
			}
		}
		if (st == QAWS_STATUS_OK && nk && *out)
		{
			/* (closed as the mesh made: an edge of four faces passes the edge count, then splits open) */
			qaws_mesh_topology tp;
			if (qaws_mesh_3d_get_topology(&tp, *out) == QAWS_STATUS_OK && !tp.closed)
				report->closed = 0;
		}
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_mesh_bvh_destroy(ba);
	qaws_mesh_bvh_destroy(bb);
	sx_out_free(&o);
	qaws_internal_free(x); qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(ta); qaws_internal_free(tb);
	qaws_internal_free(pts); qaws_internal_free(wa); qaws_internal_free(wb); qaws_internal_free(hits); qaws_internal_free(keep);
	qaws_internal_free(keep_src); qaws_internal_free(used); qaws_internal_free(kpos);
	qaws_internal_free(par); qaws_internal_free(rep); qaws_internal_free(pid); qaws_internal_free(grp); qaws_internal_free(kmap); qaws_internal_free(kopp); qaws_internal_free(fat); qaws_internal_free(near); qaws_internal_free(mid);
	qaws_internal_free(dup);
	return st;
}
