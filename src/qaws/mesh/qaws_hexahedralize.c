/*
 * Hexahedralization, method GRID (qaws_hexahedralize.h): a grid over the
 * box, the cubes inside it kept and made manifold, a layer of hexahedra
 * padding them to the surface, smoothing guarded by the scaled Jacobian.
 *
 * shortcut: a uniform grid (no octree: one cell size, the padding layer
 * absorbs the surface's detail) and features kept only where the staircase
 * turns near them; Gao et al.'s adaptive octree and feature alignment to
 * come.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_hexahedralize.h"
#include "../qaws_tetrahedralize.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct hx_grid
{
	double o[3], h[3];                  /* the origin, the steps */
	unsigned int n[3];                  /* cells per axis */
	unsigned char *in, *keep;           /* per grid vertex: inside (with a margin); per cell: kept */
} hx_grid;

static size_t hx_v(hx_grid const* g, unsigned int i, unsigned int j, unsigned int k)
{
	return (size_t)i + ((size_t)g->n[0] + 1) * ((size_t)j + ((size_t)g->n[1] + 1) * (size_t)k);
}

/* is cell (i, j, k) kept (outside the grid: no) */
static int hx_kept(hx_grid const* g, long i, long j, long k)
{
	if (i < 0 || j < 0 || k < 0 || i >= (long)g->n[0] || j >= (long)g->n[1] || k >= (long)g->n[2])
		return 0;
	return g->keep[(size_t)i + (size_t)g->n[0] * ((size_t)j + (size_t)g->n[1] * (size_t)k)];
}

static void hx_set(hx_grid* g, long i, long j, long k, unsigned char v)
{
	if (i < 0 || j < 0 || k < 0 || i >= (long)g->n[0] || j >= (long)g->n[1] || k >= (long)g->n[2])
		return;
	g->keep[(size_t)i + (size_t)g->n[0] * ((size_t)j + (size_t)g->n[1] * (size_t)k)] = v;
}

/* the cells of the 2 x 2 x 2 block round grid vertex (i, j, k) kept, as 8 bits (bit dx + 2 dy + 4 dz) */
static unsigned int hx_octants(hx_grid const* g, long i, long j, long k)
{
	unsigned int m = 0, o;
	for (o = 0; o < 8; o++)
		m |= (unsigned int)hx_kept(g, i - 1 + (long)(o & 1), j - 1 + (long)((o >> 1) & 1), k - 1 + (long)(o >> 2)) << o;
	return m;
}

/* the face-connected component of the bits of mask holding octant o */
static unsigned int hx_component(unsigned int mask, unsigned int o)
{
	unsigned int comp = 1u << o, grown = 1;
	while (grown)
	{
		unsigned int q;
		grown = 0;
		for (q = 0; q < 8; q++)
			if ((comp >> q) & 1)
			{
				unsigned int b;
				for (b = 1; b < 8; b <<= 1)
					if (((mask >> (q ^ b)) & 1) && !((comp >> (q ^ b)) & 1))
					{
						comp |= 1u << (q ^ b);
						grown = 1;
					}
			}
	}
	return comp;
}

/* the inside, along each line parallel to x: the parity of the crossings with the triangles (strict tests; a line
   whose count comes out odd, through an edge or a vertex, nudged and counted again) */
static int hx_inside(hx_grid* g, double const* x, unsigned int const* tri, unsigned int nt)
{
	unsigned int nl = (g->n[1] + 1) * (g->n[2] + 1), t, i, *l0 = NULL, *lt = NULL, *fill = NULL, cap = 16;
	double* cross = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), cap));
	int ok = 0;
	l0 = (unsigned int*)qaws_internal_zalloc((size_t)nl + 1, sizeof(unsigned int));
	if (!cross || !l0)
		goto done;
	/* the triangles by the lines their (y, z) box covers */
	for (t = 0; t < nt; t++)
	{
		double lo[2] = { HUGE_VAL, HUGE_VAL }, hi[2] = { -HUGE_VAL, -HUGE_VAL };
		long j, k, j0, j1, k0, k1;
		for (i = 0; i < 3; i++)
		{
			double const* p = x + 3 * (size_t)tri[3 * (size_t)t + i];
			lo[0] = fmin(lo[0], p[1]); hi[0] = fmax(hi[0], p[1]);
			lo[1] = fmin(lo[1], p[2]); hi[1] = fmax(hi[1], p[2]);
		}
		j0 = (long)ceil((lo[0] - g->o[1]) / g->h[1] - 1e-9) - 1; j1 = (long)floor((hi[0] - g->o[1]) / g->h[1] + 1e-9) + 1;
		k0 = (long)ceil((lo[1] - g->o[2]) / g->h[2] - 1e-9) - 1; k1 = (long)floor((hi[1] - g->o[2]) / g->h[2] + 1e-9) + 1;
		for (k = k0 < 0 ? 0 : k0; k <= k1 && k <= (long)g->n[2]; k++)
			for (j = j0 < 0 ? 0 : j0; j <= j1 && j <= (long)g->n[1]; j++)
				l0[(size_t)j + ((size_t)g->n[1] + 1) * (size_t)k + 1]++;
	}
	for (i = 0; i < nl; i++)
		l0[i + 1] += l0[i];
	lt = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)l0[nl] + 1)));
	fill = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nl + 1)));
	if (!lt || !fill)
		goto done;
	memcpy(fill, l0, sizeof(unsigned int) * nl);
	for (t = 0; t < nt; t++)
	{
		double lo[2] = { HUGE_VAL, HUGE_VAL }, hi[2] = { -HUGE_VAL, -HUGE_VAL };
		long j, k, j0, j1, k0, k1;
		for (i = 0; i < 3; i++)
		{
			double const* p = x + 3 * (size_t)tri[3 * (size_t)t + i];
			lo[0] = fmin(lo[0], p[1]); hi[0] = fmax(hi[0], p[1]);
			lo[1] = fmin(lo[1], p[2]); hi[1] = fmax(hi[1], p[2]);
		}
		j0 = (long)ceil((lo[0] - g->o[1]) / g->h[1] - 1e-9) - 1; j1 = (long)floor((hi[0] - g->o[1]) / g->h[1] + 1e-9) + 1;
		k0 = (long)ceil((lo[1] - g->o[2]) / g->h[2] - 1e-9) - 1; k1 = (long)floor((hi[1] - g->o[2]) / g->h[2] + 1e-9) + 1;
		for (k = k0 < 0 ? 0 : k0; k <= k1 && k <= (long)g->n[2]; k++)
			for (j = j0 < 0 ? 0 : j0; j <= j1 && j <= (long)g->n[1]; j++)
				lt[fill[(size_t)j + ((size_t)g->n[1] + 1) * (size_t)k]++] = t;
	}
	for (i = 0; i < nl; i++)
	{
		unsigned int j = i % (g->n[1] + 1), k = i / (g->n[1] + 1), tries, q, nc = 0, a;
		for (tries = 0; tries < 4; tries++)
		{
			double y = g->o[1] + j * g->h[1] + tries * 1.3e-7 * g->h[1], z = g->o[2] + k * g->h[2] + tries * 0.7e-7 * g->h[2];
			nc = 0;
			for (q = l0[i]; q < l0[i + 1]; q++)
			{
				unsigned int const* T = tri + 3 * (size_t)lt[q];
				double const *A = x + 3 * (size_t)T[0], *B = x + 3 * (size_t)T[1], *C = x + 3 * (size_t)T[2];
				double s0 = (B[1] - y) * (C[2] - z) - (B[2] - z) * (C[1] - y), s1 = (C[1] - y) * (A[2] - z) - (C[2] - z) * (A[1] - y);
				double s2 = (A[1] - y) * (B[2] - z) - (A[2] - z) * (B[1] - y), s = s0 + s1 + s2;
				if (!((s0 > 0 && s1 > 0 && s2 > 0) || (s0 < 0 && s1 < 0 && s2 < 0)))
					continue;
				if (nc == cap)
				{
					unsigned int nc2 = qaws_internal_grow(cap, (size_t)nc + 1, 16);
					double* gr = nc2 ? (double*)qaws_internal_realloc(cross, qaws_internal_mul2(sizeof(double), nc2)) : NULL;
					if (!gr)
						goto done;
					cross = gr;
					cap = nc2;
				}
				cross[nc++] = (s0 * A[0] + s1 * B[0] + s2 * C[0]) / s;
			}
			if (!(nc & 1))
				break;
		}
		/* a vertex is inside when an odd number of crossings lie before it */
		for (a = 0; a <= g->n[0]; a++)
		{
			double px = g->o[0] + a * g->h[0];
			unsigned int before = 0;
			for (q = 0; q < nc; q++)
				before += cross[q] < px;
			g->in[hx_v(g, a, j, k)] = (unsigned char)(!(nc & 1) && (before & 1));
		}
	}
	ok = 1;
done:
	qaws_internal_free(cross); qaws_internal_free(l0); qaws_internal_free(lt); qaws_internal_free(fill);
	return ok;
}

static unsigned int hx_bits(unsigned int m)
{
	unsigned int n = 0;
	for (; m; m >>= 1)
		n += m & 1;
	return n;
}

/* the cubes made a manifold set: round every grid vertex the kept octants one face-connected piece and the empty ones
   too (else cubes meet at an edge or a vertex only); the smaller kept pieces dropped, or one cube where the empty
   octants are split; repeated until none */
static void hx_manifold(hx_grid* g)
{
	unsigned int pass, changed = 1;
	for (pass = 0; pass < 64 && changed; pass++)
	{
		long i, j, k;
		changed = 0;
		for (k = 0; k <= (long)g->n[2]; k++)
			for (j = 0; j <= (long)g->n[1]; j++)
				for (i = 0; i <= (long)g->n[0]; i++)
				{
					unsigned int m = hx_octants(g, i, j, k), e = ~m & 255u, o, best = 0, drop;
					if (m == 0 || m == 255)
						continue;
					for (o = 0; !((m >> o) & 1); o++) {}
					if (hx_component(m, o) == m)
					{
						for (o = 0; !((e >> o) & 1); o++) {}
						if (hx_component(e, o) == e)
							continue;
						for (o = 0; !((m >> o) & 1); o++) {}
						drop = 1u << o;
					}
					else
					{
						for (o = 0; o < 8; o++)
							if (((m >> o) & 1) && hx_bits(hx_component(m, o)) > hx_bits(best))
								best = hx_component(m, o);
						drop = m & ~best;
					}
					for (o = 0; o < 8; o++)
						if ((drop >> o) & 1)
							hx_set(g, i - 1 + (long)(o & 1), j - 1 + (long)((o >> 1) & 1), k - 1 + (long)(o >> 2), 0);
					changed = 1;
				}
	}
}

/* ---- smoothing ------------------------------------------------------------ */

typedef struct hx_csr
{
	unsigned int *v0, *cell;            /* per vertex: its cells */
} hx_csr;

static int hx_csr_build(hx_csr* s, unsigned int nv, unsigned int const* cells, unsigned int nc)
{
	unsigned int i, *fill;
	s->v0 = (unsigned int*)qaws_internal_zalloc((size_t)nv + 1, sizeof(unsigned int));
	s->cell = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (8 * (size_t)nc + 1)));
	fill = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nv + 1)));
	if (!s->v0 || !s->cell || !fill)
	{
		qaws_internal_free(fill);
		return 0;
	}
	for (i = 0; i < 8 * nc; i++)
		s->v0[cells[i] + 1]++;
	for (i = 0; i < nv; i++)
		s->v0[i + 1] += s->v0[i];
	memcpy(fill, s->v0, sizeof(unsigned int) * nv);
	for (i = 0; i < 8 * nc; i++)
		s->cell[fill[cells[i]]++] = i / 8;
	qaws_internal_free(fill);
	return 1;
}

static double hx_worst(double const* x, unsigned int const* cells, hx_csr const* s, unsigned int v)
{
	double q = HUGE_VAL;
	unsigned int i;
	for (i = s->v0[v]; i < s->v0[v + 1]; i++)
		q = fmin(q, vm_hex_quality(x, cells + 8 * (size_t)s->cell[i]));
	return q;
}

/* ======================================================================== */
/*  FRAME_FIELD: an integer-grid map of a field without singular edges      */
/* ======================================================================== */

/* a sparse symmetric matrix in compressed rows */
typedef struct hx_mat
{
	unsigned int *r0, *col;
	double* val;
	unsigned int n;
} hx_mat;

typedef struct hx_trip
{
	unsigned long long key;             /* row << 32 | column */
	double v;
} hx_trip;

static int hx_trip_cmp(void const* x, void const* y)
{
	unsigned long long a = ((hx_trip const*)x)->key, b = ((hx_trip const*)y)->key;
	return a < b ? -1 : a > b ? 1 : 0;
}

/* the tetrahedra's barycentric gradients: g[4][3] of cell c, its volume returned */
static double hx_grads(qaws_volume_mesh const* t, unsigned int c, double g[4][3])
{
	double const* p[4];
	double u[3], w[3], d[3], cr[3], V;
	unsigned int i, k, n;
	for (i = 0; i < 4; i++)
		p[i] = t->x + 3 * (size_t)t->cv[t->c0[c] + i];
	for (k = 0; k < 3; k++)
	{
		u[k] = p[1][k] - p[0][k];
		w[k] = p[2][k] - p[0][k];
		d[k] = p[3][k] - p[0][k];
	}
	cr[0] = w[1] * d[2] - w[2] * d[1];
	cr[1] = w[2] * d[0] - w[0] * d[2];
	cr[2] = w[0] * d[1] - w[1] * d[0];
	V = (u[0] * cr[0] + u[1] * cr[1] + u[2] * cr[2]) / 6;
	/* the gradient of the coordinate of vertex i: minus its opposite facet's outward area vector over 3 V */
	for (i = 0; i < 4; i++)
	{
		unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, i, &n);
		double const *a = p[loc[0]], *b = p[loc[1]], *e = p[loc[2]];
		double s[3], q[3], N[3];
		for (k = 0; k < 3; k++)
		{
			s[k] = b[k] - a[k];
			q[k] = e[k] - a[k];
		}
		N[0] = s[1] * q[2] - s[2] * q[1];
		N[1] = s[2] * q[0] - s[0] * q[2];
		N[2] = s[0] * q[1] - s[1] * q[0];
		for (k = 0; k < 3; k++)
			g[i][k] = V > 0 ? -N[k] / (6 * V) : 0;
	}
	return V;
}

/* the stiffness of the tetrahedra (sum of V g_i . g_j) */
static int hx_stiffness(qaws_volume_mesh const* t, hx_mat* K)
{
	hx_trip* tr = (hx_trip*)qaws_internal_malloc(qaws_internal_mul2(sizeof(hx_trip), (16 * (size_t)t->nc + 1)));
	unsigned int c, i, j, n = 0, k;
	memset(K, 0, sizeof *K);
	if (!tr)
		return 0;
	for (c = 0; c < t->nc; c++)
	{
		double g[4][3], V = hx_grads(t, c, g);
		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
			{
				tr[n].key = ((unsigned long long)t->cv[t->c0[c] + i] << 32) | t->cv[t->c0[c] + j];
				tr[n].v = V * (g[i][0] * g[j][0] + g[i][1] * g[j][1] + g[i][2] * g[j][2]);
				n++;
			}
	}
	qsort(tr, n, sizeof(hx_trip), hx_trip_cmp);
	K->n = t->nv;
	K->r0 = (unsigned int*)qaws_internal_zalloc((size_t)t->nv + 1, sizeof(unsigned int));
	K->col = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)n + 1)));
	K->val = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n + 1)));
	if (!K->r0 || !K->col || !K->val)
	{
		qaws_internal_free(tr);
		return 0;
	}
	for (i = 0, k = 0; i < n; i++)
	{
		if (k && tr[i].key == tr[i - 1].key)
		{
			K->val[k - 1] += tr[i].v;
			continue;
		}
		K->col[k] = (unsigned int)(tr[i].key & 0xFFFFFFFFu);
		K->val[k] = tr[i].v;
		K->r0[(tr[i].key >> 32) + 1]++;
		k++;
	}
	for (i = 0; i < t->nv; i++)
		K->r0[i + 1] += K->r0[i];
	qaws_internal_free(tr);
	return 1;
}

/* K u = b where not fixed (u kept where fixed): Jacobi-preconditioned conjugate gradients; 0 on allocation failure */
static int hx_cg(hx_mat const* K, double const* b, unsigned char const* fixed, double* u)
{
	unsigned int n = K->n, i, it, j;
	double *r = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n + 1))), *p = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n + 1)));
	double *Ap = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n + 1))), *dinv = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n + 1))), rz = 0, r0;
	if (!r || !p || !Ap || !dinv)
	{
		qaws_internal_free(r); qaws_internal_free(p); qaws_internal_free(Ap); qaws_internal_free(dinv);
		return 0;
	}
	for (i = 0; i < n; i++)
	{
		double s = b[i], d = 0;
		for (j = K->r0[i]; j < K->r0[i + 1]; j++)
		{
			s -= K->val[j] * u[K->col[j]];
			if (K->col[j] == i)
				d = K->val[j];
		}
		r[i] = fixed[i] ? 0 : s;
		dinv[i] = !fixed[i] && d > 0 ? 1 / d : 0;
		p[i] = dinv[i] * r[i];
		rz += r[i] * p[i];
	}
	r0 = rz;
	for (it = 0; it < 20000 && rz > 1e-28 * r0 && rz > 0; it++)
	{
		double pAp = 0, alpha, rz2 = 0;
		for (i = 0; i < n; i++)
		{
			double s = 0;
			for (j = K->r0[i]; j < K->r0[i + 1]; j++)
				s += K->val[j] * p[K->col[j]];
			Ap[i] = fixed[i] ? 0 : s;
			pAp += p[i] * Ap[i];
		}
		if (!(pAp > 0))
			break;
		alpha = rz / pAp;
		for (i = 0; i < n; i++)
		{
			u[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			rz2 += r[i] * dinv[i] * r[i];
		}
		for (i = 0; i < n; i++)
			p[i] = dinv[i] * r[i] + rz2 / rz * p[i];
		rz = rz2;
	}
	qaws_internal_free(r); qaws_internal_free(p); qaws_internal_free(Ap); qaws_internal_free(dinv);
	return 1;
}

static unsigned int hx_uf(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

/* integer lattice points to indices, open addressing */
typedef struct hx_pmap
{
	unsigned long long* key;
	unsigned int* val;
	unsigned int cap, n;
} hx_pmap;

static unsigned long long hx_pkey(long i, long j, long k)
{
	return ((unsigned long long)(i + 1048576) << 42) | ((unsigned long long)(j + 1048576) << 21) | (unsigned long long)(k + 1048576);
}

static unsigned int hx_pget(hx_pmap const* h, unsigned long long key)
{
	unsigned int s = (unsigned int)(key * 0x9E3779B97F4A7C15ull >> 32) & (h->cap - 1);
	while (h->key[s] != ~0ull)
	{
		if (h->key[s] == key)
			return h->val[s];
		s = (s + 1) & (h->cap - 1);
	}
	return HM_NONE;
}

/* 0 on allocation failure; an existing key keeps its value */
static int hx_pput(hx_pmap* h, unsigned long long key, unsigned int val)
{
	unsigned int s;
	if (qaws_internal_mul2(2, (size_t)h->n + 1) > h->cap)
	{
		hx_pmap g;
		unsigned int i;
		g.cap = qaws_internal_grow(h->cap, qaws_internal_mul2(2, (size_t)h->n + 1), 1024);
		if (!g.cap || (g.cap & (g.cap - 1)))
			return 0;   /* (a power of two: the probe masks with cap - 1) */
		g.n = 0;
		g.key = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned long long), g.cap));
		g.val = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), g.cap));
		if (!g.key || !g.val)
		{
			qaws_internal_free(g.key);
			qaws_internal_free(g.val);
			return 0;
		}
		memset(g.key, 0xFF, sizeof(unsigned long long) * g.cap);
		for (i = 0; i < h->cap; i++)
			if (h->key[i] != ~0ull)
				hx_pput(&g, h->key[i], h->val[i]);
		qaws_internal_free(h->key);
		qaws_internal_free(h->val);
		*h = g;
	}
	s = (unsigned int)(key * 0x9E3779B97F4A7C15ull >> 32) & (h->cap - 1);
	while (h->key[s] != ~0ull)
	{
		if (h->key[s] == key)
			return 1;
		s = (s + 1) & (h->cap - 1);
	}
	h->key[s] = key;
	h->val[s] = val;
	h->n++;
	return 1;
}

/* Hexahedra by an integer-grid map (CubeCover's idea, for fields without singular edges): the inside tetrahedralized,
   its frame field combed into one gauge, three Poisson solves fitting a map u whose gradient follows the frames at the
   cell size (each border patch, facets whose normal follows one axis, held on a whole-number plane of that axis), the
   hexahedra read off the integer lattice in the map's image. */
static qaws_status hx_frame_field(qaws_mesh_3d const* boundary, double h, qaws_volume_mesh** out_mesh, qaws_mesh_report* rep)
{
	qaws_tetrahedralize_desc td;
	qaws_volume_mesh* t = NULL;
	hx_mat K;
	hx_pmap pts, cen;
	double *R = NULL, *U = NULL, *b = NULL, *u = NULL, *px = NULL;
	unsigned char *sing = NULL, *fixed = NULL;
	unsigned int *bf = NULL, *uf = NULL, *lab = NULL, *cells = NULL, nb = 0, c, f, i, k, np = 0, capp = 0, ncell = 0, capc = 0;
	qaws_status st;
	memset(&K, 0, sizeof K);
	memset(&pts, 0, sizeof pts);
	memset(&cen, 0, sizeof cen);
	memset(&td, 0, sizeof td);
	td.target_length = (qaws_scalar)(0.5 * h);
	st = qaws_volume_mesh_tetrahedralize(boundary, &td, &t, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	R = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (9 * (size_t)t->nc + 1)));
	sing = (unsigned char*)qaws_internal_scratch_malloc((size_t)t->nc + 1);
	if (!R || !sing)
		goto done;
	st = ff_frames(t, 60, R, sing);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_UNSUPPORTED_OPERATION;   /* (a field with singular edges, or turning round a hole) */
	for (c = 0; c < t->nc; c++)
		if (sing[c])
			goto done;
	if (!ff_comb(t, R))
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	/* the border's facets labelled by the axis their normal follows (and its side), patches of like neighbours */
	bf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)t->nc + 1)));
	lab = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)t->nc + 1)));
	uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)t->nc + 1)));
	if (!bf || !lab || !uf)
		goto done;
	for (c = 0; c < t->nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int n, best = 0;
			unsigned char const* loc;
			double const *a, *bb, *e;
			double s[3], q[3], N[3], dmax = -1;
			if (t->fo[t->f0[c] + f] != HM_NONE)
				continue;
			loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n);
			a = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[0]];
			bb = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[1]];
			e = t->x + 3 * (size_t)t->cv[t->c0[c] + loc[2]];
			for (k = 0; k < 3; k++)
			{
				s[k] = bb[k] - a[k];
				q[k] = e[k] - a[k];
			}
			N[0] = s[1] * q[2] - s[2] * q[1];
			N[1] = s[2] * q[0] - s[0] * q[2];
			N[2] = s[0] * q[1] - s[1] * q[0];
			for (k = 0; k < 3; k++)
			{
				double d = R[9 * (size_t)c + 3 * k] * N[0] + R[9 * (size_t)c + 3 * k + 1] * N[1] + R[9 * (size_t)c + 3 * k + 2] * N[2];
				if (fabs(d) > dmax)
				{
					dmax = fabs(d);
					best = 2 * k + (d > 0);
				}
			}
			bf[nb] = 4 * c + f;
			lab[nb] = best;
			uf[nb] = nb;
			nb++;
		}
	{
		/* facets sharing an edge with the same label: one patch */
		hx_pmap em;
		memset(&em, 0, sizeof em);
		for (i = 0; i < nb; i++)
		{
			unsigned int n, e2;
			unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, bf[i] % 4, &n);
			for (e2 = 0; e2 < 3; e2++)
			{
				unsigned int va = t->cv[t->c0[bf[i] / 4] + loc[e2]], vb = t->cv[t->c0[bf[i] / 4] + loc[(e2 + 1) % 3]], other;
				unsigned long long key = va < vb ? ((unsigned long long)va << 32) | vb : ((unsigned long long)vb << 32) | va;
				other = em.cap ? hx_pget(&em, key) : HM_NONE;
				if (other == HM_NONE)
				{
					if (!hx_pput(&em, key, i))
					{
						qaws_internal_free(em.key);
						qaws_internal_free(em.val);
						goto done;
					}
				}
				else if (lab[other] == lab[i])
				{
					unsigned int ra = hx_uf(uf, other), rb = hx_uf(uf, i);
					if (ra != rb)
						uf[ra] = rb;
				}
			}
		}
		qaws_internal_free(em.key);
		qaws_internal_free(em.val);
	}
	/* the map: a free solve per axis (one vertex pinned), each patch's plane its mean rounded, then solved again with them */
	if (!hx_stiffness(t, &K))
		goto done;
	U = (double*)qaws_internal_scratch_zalloc(3 * (size_t)t->nv + 1, sizeof(double));
	b = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)t->nv + 1)));
	u = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)t->nv + 1)));
	fixed = (unsigned char*)qaws_internal_scratch_malloc((size_t)t->nv + 1);
	if (!U || !b || !u || !fixed)
		goto done;
	for (k = 0; k < 3; k++)
	{
		double *sum = (double*)qaws_internal_scratch_zalloc((size_t)nb + 1, sizeof(double)), *cnt = (double*)qaws_internal_zalloc((size_t)nb + 1, sizeof(double));
		int pass;
		if (!sum || !cnt)
		{
			qaws_internal_free(sum);
			qaws_internal_free(cnt);
			goto done;
		}
		/* the right-hand side: sum over cells of V g_i . (axis k / h) */
		memset(b, 0, sizeof(double) * t->nv);
		for (c = 0; c < t->nc; c++)
		{
			double g[4][3], V = hx_grads(t, c, g);
			double const* a = R + 9 * (size_t)c + 3 * k;
			for (i = 0; i < 4; i++)
				b[t->cv[t->c0[c] + i]] += V * (g[i][0] * a[0] + g[i][1] * a[1] + g[i][2] * a[2]) / h;
		}
		for (pass = 0; pass < 2; pass++)
		{
			memset(fixed, 0, t->nv);
			if (pass == 0)
			{
				memset(u, 0, sizeof(double) * t->nv);
				fixed[0] = 1;
			}
			else
			{
				/* the patches of this axis: their planes at the nearest whole numbers */
				for (i = 0; i < nb; i++)
					if (lab[i] / 2 == k)
					{
						unsigned int n, q;
						unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, bf[i] % 4, &n);
						for (q = 0; q < 3; q++)
						{
							sum[hx_uf(uf, i)] += u[t->cv[t->c0[bf[i] / 4] + loc[q]]];
							cnt[hx_uf(uf, i)] += 1;
						}
					}
				for (i = 0; i < nb; i++)
					if (lab[i] / 2 == k)
					{
						unsigned int n, q, r = hx_uf(uf, i);
						unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, bf[i] % 4, &n);
						double plane = floor(sum[r] / cnt[r] + 0.5);
						for (q = 0; q < 3; q++)
						{
							unsigned int v = t->cv[t->c0[bf[i] / 4] + loc[q]];
							if (!fixed[v])
							{
								fixed[v] = 1;
								u[v] = plane;
							}
						}
					}
			}
			if (!hx_cg(&K, b, fixed, u))
			{
				qaws_internal_free(sum);
				qaws_internal_free(cnt);
				goto done;
			}
		}
		for (i = 0; i < t->nv; i++)
			U[3 * (size_t)i + k] = u[i];
		qaws_internal_free(sum);
		qaws_internal_free(cnt);
	}
	/* the lattice: each cell's integer points and cube centres in its image, located by barycentric coordinates */
	for (c = 0; c < t->nc; c++)
	{
		double const *U0 = U + 3 * (size_t)t->cv[t->c0[c]], *X[4];
		double M[9], inv[9], det, lo[3], hi[3];
		unsigned int q;
		int half;
		for (q = 0; q < 4; q++)
			X[q] = t->x + 3 * (size_t)t->cv[t->c0[c] + q];
		for (q = 0; q < 3; q++)
			for (k = 0; k < 3; k++)
				M[3 * q + k] = U[3 * (size_t)t->cv[t->c0[c] + q + 1] + k] - U0[k];   /* (column q: u of vertex q + 1 less u of vertex 0) */
		/* the inverse by cofactors */
		inv[0] = M[4] * M[8] - M[5] * M[7]; inv[3] = M[5] * M[6] - M[3] * M[8]; inv[6] = M[3] * M[7] - M[4] * M[6];
		inv[1] = M[7] * M[2] - M[8] * M[1]; inv[4] = M[8] * M[0] - M[6] * M[2]; inv[7] = M[6] * M[1] - M[7] * M[0];
		inv[2] = M[1] * M[5] - M[2] * M[4]; inv[5] = M[2] * M[3] - M[0] * M[5]; inv[8] = M[0] * M[4] - M[1] * M[3];
		det = M[0] * inv[0] + M[3] * inv[1] + M[6] * inv[2];
		if (!(det > 1e-12))
			continue;   /* (turned over or flat in the map) */
		for (k = 0; k < 3; k++)
		{
			lo[k] = U0[k];
			hi[k] = U0[k];
			for (q = 1; q < 4; q++)
			{
				lo[k] = fmin(lo[k], U[3 * (size_t)t->cv[t->c0[c] + q] + k]);
				hi[k] = fmax(hi[k], U[3 * (size_t)t->cv[t->c0[c] + q] + k]);
			}
		}
		for (half = 0; half < 2; half++)
		{
			long i0 = (long)ceil(lo[0] - 0.5 * half - 1e-9), i1 = (long)floor(hi[0] - 0.5 * half + 1e-9);
			long j0 = (long)ceil(lo[1] - 0.5 * half - 1e-9), j1 = (long)floor(hi[1] - 0.5 * half + 1e-9);
			long k0 = (long)ceil(lo[2] - 0.5 * half - 1e-9), k1 = (long)floor(hi[2] - 0.5 * half + 1e-9);
			long a, bb2, e;
			for (a = i0; a <= i1; a++)
				for (bb2 = j0; bb2 <= j1; bb2++)
					for (e = k0; e <= k1; e++)
					{
						double d[3], l[3], l0;
						unsigned long long key = hx_pkey(a, bb2, e);
						d[0] = a + 0.5 * half - U0[0];
						d[1] = bb2 + 0.5 * half - U0[1];
						d[2] = e + 0.5 * half - U0[2];
						for (q = 0; q < 3; q++)
							l[q] = (inv[q] * d[0] + inv[3 + q] * d[1] + inv[6 + q] * d[2]) / det;
						l0 = 1 - l[0] - l[1] - l[2];
						if (l0 < -1e-9 || l[0] < -1e-9 || l[1] < -1e-9 || l[2] < -1e-9)
							continue;
						if (half)
						{
							if (!hx_pput(&cen, key, 0))
								goto done;
							continue;
						}
						if (pts.cap && hx_pget(&pts, key) != HM_NONE)
							continue;
						if (np + 1 > capp)
						{
							unsigned int nc = qaws_internal_grow(capp, (size_t)np + 1, 1024);
							double* g2 = nc ? (double*)qaws_internal_realloc(px, qaws_internal_mul3(sizeof(double), 3, nc)) : NULL;
							if (!g2)
								goto done;
							px = g2;
							capp = nc;
						}
						for (q = 0; q < 3; q++)
							px[3 * (size_t)np + q] = l0 * X[0][q] + l[0] * X[1][q] + l[1] * X[2][q] + l[2] * X[3][q];
						if (!hx_pput(&pts, key, np))
							goto done;
						np++;
					}
		}
	}
	/* a hexahedron for each cube whose centre and eight corners are in the image */
	for (i = 0; i < cen.cap; i++)
	{
		static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
		long a, bb2, e;
		unsigned int q, v[8], all = 1;
		unsigned long long key = cen.key[i];
		if (key == ~0ull)
			continue;
		a = (long)((key >> 42) & 0x1FFFFF) - 1048576;
		bb2 = (long)((key >> 21) & 0x1FFFFF) - 1048576;
		e = (long)(key & 0x1FFFFF) - 1048576;
		for (q = 0; q < 8 && all; q++)
		{
			v[q] = pts.cap ? hx_pget(&pts, hx_pkey(a + CO[q][0], bb2 + CO[q][1], e + CO[q][2])) : HM_NONE;
			all = v[q] != HM_NONE;
		}
		if (!all)
			continue;
		if (ncell + 1 > capc)
		{
			unsigned int nc = qaws_internal_grow(capc, (size_t)ncell + 1, 256);
			unsigned int* g2 = nc ? (unsigned int*)qaws_internal_realloc(cells, qaws_internal_mul3(sizeof(unsigned int), 8, nc)) : NULL;
			if (!g2)
				goto done;
			cells = g2;
			capc = nc;
		}
		memcpy(cells + 8 * (size_t)ncell++, v, sizeof v);
	}
	if (!ncell)
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;   /* (no cube fits: a cell size too large) */
		goto done;
	}
	{
		/* the vertices the hexahedra use, compacted */
		unsigned int* map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)np + 1))), nv2 = 0;
		double* x2 = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)np + 1)));
		unsigned char* types = (unsigned char*)qaws_internal_scratch_malloc((size_t)ncell + 1);
		if (!map || !x2 || !types)
		{
			qaws_internal_free(map); qaws_internal_free(x2); qaws_internal_free(types);
			goto done;
		}
		for (i = 0; i < np; i++)
			map[i] = HM_NONE;
		for (i = 0; i < 8 * ncell; i++)
		{
			if (map[cells[i]] == HM_NONE)
			{
				memcpy(x2 + 3 * (size_t)nv2, px + 3 * (size_t)cells[i], sizeof(double) * 3);
				map[cells[i]] = nv2++;
			}
			cells[i] = map[cells[i]];
		}
		memset(types, QAWS_CELL_HEXAHEDRON, ncell);
		st = vm_create(x2, nv2, cells, types, ncell, out_mesh);
		qaws_internal_free(map); qaws_internal_free(x2); qaws_internal_free(types);
	}
	if (st == QAWS_STATUS_OK)
	{
		qaws_mesh_3d* bd = NULL;
		if (qaws_volume_mesh_get_boundary(*out_mesh, &bd) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(bd, &rep->quality);
		qaws_mesh_3d_destroy(bd);
	}
done:
	qaws_internal_free(R); qaws_internal_free(sing); qaws_internal_free(bf); qaws_internal_free(uf); qaws_internal_free(lab); qaws_internal_free(U); qaws_internal_free(b); qaws_internal_free(u); qaws_internal_free(fixed); qaws_internal_free(px); qaws_internal_free(cells);
	qaws_internal_free(K.r0); qaws_internal_free(K.col); qaws_internal_free(K.val);
	qaws_internal_free(pts.key); qaws_internal_free(pts.val); qaws_internal_free(cen.key); qaws_internal_free(cen.val);
	qaws_volume_mesh_destroy(t);
	return st;
}

/* ---- the operator ----------------------------------------------------------- */

qaws_status qaws_volume_mesh_hexahedralize(qaws_mesh_3d const* boundary, qaws_hexahedralize_desc const* desc, qaws_volume_mesh** out_mesh,
	qaws_mesh_report* out_report)
{
	static unsigned char const F[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static int const FN[6][3] = { { 0, 0, -1 }, { 0, 0, 1 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 } };
	static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
	qaws_hm const* hm;
	qaws_mesh_reference* ref = NULL;
	mr_scratch s;
	hx_grid g;
	hx_csr csr;
	qaws_mesh_report rep;
	double *x = NULL, lo[3], hi[3], L, sum = 0, *cx = NULL;
	unsigned int *tri = NULL, nt, *vid = NULL, *cells = NULL, *outer = NULL, ncell = 0, nvx = 0, nb = 0, i, k, c, ne = 0, capc = 0, *corners = NULL, ncorner = 0;
	unsigned char *fixed = NULL, *kind = NULL;
	size_t ngv;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!boundary || !desc || (unsigned int)desc->method > QAWS_HEXES_FRAME_FIELD || desc->cell_size < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->method == QAWS_HEXES_POLYCUBE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	hm = &boundary->hm;
	for (c = 0; c < hm->nc; c++)
		if (hm->cf[c] != HM_NONE && hm->co[c] == HM_NONE)
			return QAWS_STATUS_INVALID_TOPOLOGY;
	if (!hm->nf)
		return QAWS_STATUS_INVALID_TOPOLOGY;
	memset(&g, 0, sizeof g);
	memset(&s, 0, sizeof s);
	memset(&csr, 0, sizeof csr);
	memset(&rep, 0, sizeof rep);
	/* the surface's triangles, its box, its mean edge */
	nt = qaws_mesh_3d_get_triangle_count(boundary);
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (3 * (size_t)nt + 1)));
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)hm->nv + 1)));
	if (!tri || !x)
		goto done;
	qaws_mesh_3d_get_triangles(boundary, tri);
	memcpy(x, hm->x, sizeof(double) * 3 * hm->nv);
	for (k = 0; k < 3; k++)
	{
		lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	for (i = 0; i < hm->nv; i++)
		if (hm->vc[i] != HM_NONE)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], x[3 * (size_t)i + k]);
				hi[k] = fmax(hi[k], x[3 * (size_t)i + k]);
			}
	for (i = 0; i < nt; i++)
		for (k = 0; k < 3; k++)
		{
			double const *a = x + 3 * (size_t)tri[3 * i + k], *b = x + 3 * (size_t)tri[3 * i + (k + 1) % 3];
			sum += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
			ne++;
		}
	L = desc->cell_size > 0 ? (double)desc->cell_size : sum / (ne ? ne : 1);
	if (!(L > 0) || !(hi[0] >= lo[0]))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto done;
	}
	if (desc->method == QAWS_HEXES_FRAME_FIELD)
	{
		st = hx_frame_field(boundary, L, out_mesh, &rep);
		goto done;
	}
	/* the grid: the box in whole cells per axis (the step adjusted, so a box's faces fall on grid planes) */
	for (k = 0; k < 3; k++)
	{
		double ext = hi[k] - lo[k], n = fmax(1, floor(ext / L + 0.5));
		if (n > 4096)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;   /* (a cell size too small for the box) */
			goto done;
		}
		g.n[k] = (unsigned int)n;
		g.h[k] = ext > 0 ? ext / n : L;
		g.o[k] = lo[k];
	}
	ngv = ((size_t)g.n[0] + 1) * ((size_t)g.n[1] + 1) * ((size_t)g.n[2] + 1);
	if (ngv > ((size_t)1 << 26))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto done;
	}
	g.in = (unsigned char*)qaws_internal_zalloc(ngv, 1);
	g.keep = (unsigned char*)qaws_internal_zalloc(qaws_internal_add2(qaws_internal_mul3(g.n[0], g.n[1], g.n[2]), 1), 1);
	if (!g.in || !g.keep || !hx_inside(&g, x, tri, nt))
		goto done;
	/* the reference: projections, the margin from the surface */
	st = qaws_mesh_reference_create(boundary, NULL, &ref);
	if (st != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!mr_scratch_init(&s, ref->nt))
		goto done;
	{
		double hmin = fmin(g.h[0], fmin(g.h[1], g.h[2])), margin2 = 0.3 * hmin * 0.3 * hmin;
		unsigned int a, b, d;
		for (d = 0; d <= g.n[2]; d++)
			for (b = 0; b <= g.n[1]; b++)
				for (a = 0; a <= g.n[0]; a++)
				{
					size_t v = hx_v(&g, a, b, d);
					double p[3];
					qaws_mesh_location l;
					if (!g.in[v])
						continue;
					p[0] = g.o[0] + a * g.h[0];
					p[1] = g.o[1] + b * g.h[1];
					p[2] = g.o[2] + d * g.h[2];
					if (mr_project(ref, p, HM_NONE, &s, &l) < margin2)
						g.in[v] = 0;
				}
		for (d = 0; d < g.n[2]; d++)
			for (b = 0; b < g.n[1]; b++)
				for (a = 0; a < g.n[0]; a++)
				{
					unsigned int o, all = 1;
					for (o = 0; o < 8; o++)
						all &= g.in[hx_v(&g, a + (unsigned int)CO[o][0], b + (unsigned int)CO[o][1], d + (unsigned int)CO[o][2])];
					g.keep[(size_t)a + (size_t)g.n[0] * ((size_t)b + (size_t)g.n[1] * (size_t)d)] = (unsigned char)all;
				}
	}
	hx_manifold(&g);
	/* the staircase: its cells and vertices */
	vid = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ngv));
	if (!vid)
		goto done;
	for (i = 0; i < ngv; i++)
		vid[i] = HM_NONE;
	{
		unsigned int a, b, d, capv = 1024;
		cx = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, capv));
		if (!cx)
			goto done;
		for (d = 0; d < g.n[2]; d++)
			for (b = 0; b < g.n[1]; b++)
				for (a = 0; a < g.n[0]; a++)
				{
					unsigned int o, *cv;
					if (!hx_kept(&g, a, b, d))
						continue;
					if (ncell + 1 > capc)
					{
						unsigned int nc = qaws_internal_grow(capc, (size_t)ncell + 1, 256);
						unsigned int* gr = nc ? (unsigned int*)qaws_internal_realloc(cells, qaws_internal_mul3(sizeof(unsigned int), 8, nc)) : NULL;
						if (!gr)
							goto done;
						cells = gr;
						capc = nc;
					}
					cv = cells + 8 * (size_t)ncell++;
					for (o = 0; o < 8; o++)
					{
						size_t gv = hx_v(&g, a + (unsigned int)CO[o][0], b + (unsigned int)CO[o][1], d + (unsigned int)CO[o][2]);
						if (vid[gv] == HM_NONE)
						{
							if (nvx == capv)
							{
								unsigned int nc = qaws_internal_grow(capv, (size_t)nvx + 1, 1024);
								double* gr = nc ? (double*)qaws_internal_realloc(cx, qaws_internal_mul3(sizeof(double), 3, nc)) : NULL;
								if (!gr)
									goto done;
								cx = gr;
								capv = nc;
							}
							cx[3 * (size_t)nvx] = g.o[0] + (a + CO[o][0]) * g.h[0];
							cx[3 * (size_t)nvx + 1] = g.o[1] + (b + CO[o][1]) * g.h[1];
							cx[3 * (size_t)nvx + 2] = g.o[2] + (d + CO[o][2]) * g.h[2];
							vid[gv] = nvx++;
						}
						cv[o] = vid[gv];
					}
				}
		if (!ncell)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;   /* (no cube fits inside) */
			goto done;
		}
		/* the staircase's border quads, each padded out to the surface by a hexahedron (at most one projection a vertex) */
		{
			double* gr = (double*)qaws_internal_realloc(cx, qaws_internal_mul3(sizeof(double), 6, ((size_t)nvx + 1)));
			if (!gr)
				goto done;
			cx = gr;
		}
		outer = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (2 * (size_t)nvx + 1)));
		kind = (unsigned char*)qaws_internal_zalloc(2 * (size_t)nvx + 1, 1);
		if (!outer || !kind)
			goto done;
		for (i = 0; i < 2 * nvx; i++)
			outer[i] = HM_NONE;
		for (k = 0; k < ref->mesh->hm.nv; k++)
			if (ref->mesh->hm.vflag[k] & HM_CORNER)
				ncorner++;
		corners = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)ncorner + 1)));
		if (!corners)
			goto done;
		for (k = 0, ncorner = 0; k < ref->mesh->hm.nv; k++)
			if (ref->mesh->hm.vflag[k] & HM_CORNER)
				corners[ncorner++] = k;
		nb = ncell;
		for (c = 0; c < nb; c++)
		{
			unsigned int f;
			/* the cell's grid position back from its first vertex */
			double const* p0 = cx + 3 * (size_t)cells[8 * (size_t)c];
			long ca = (long)floor((p0[0] - g.o[0]) / g.h[0] + 0.5), cb = (long)floor((p0[1] - g.o[1]) / g.h[1] + 0.5), cd = (long)floor((p0[2] - g.o[2]) / g.h[2] + 0.5);
			for (f = 0; f < 6; f++)
			{
				unsigned int q, *cv;
				if (hx_kept(&g, ca + FN[f][0], cb + FN[f][1], cd + FN[f][2]))
					continue;
				/* its corners' projections, made once each */
				for (q = 0; q < 4; q++)
				{
					unsigned int v = cells[8 * (size_t)c + F[f][q]];
					if (outer[v] != HM_NONE)
						continue;
					{
						double const* pv = cx + 3 * (size_t)v;
						long ga = (long)floor((pv[0] - g.o[0]) / g.h[0] + 0.5), gb = (long)floor((pv[1] - g.o[1]) / g.h[1] + 0.5), gd = (long)floor((pv[2] - g.o[2]) / g.h[2] + 0.5);
						unsigned int cnt = hx_bits(hx_octants(&g, ga, gb, gd)), z, kd = 0;
						qaws_mesh_location l;
						double y[3], dc = HUGE_VAL, best = HUGE_VAL;
						l.patch = HM_NONE;
						/* where the staircase turns twice (a convex or concave corner): a corner of the surface nearby;
						   where it turns once: a crease nearby; else the closest point */
						if (cnt == 1 || cnt == 7)
							for (z = 0; z < ncorner; z++)
							{
								double const* pc = ref->mesh->hm.x + 3 * (size_t)corners[z];
								double dd = (pc[0] - pv[0]) * (pc[0] - pv[0]) + (pc[1] - pv[1]) * (pc[1] - pv[1]) + (pc[2] - pv[2]) * (pc[2] - pv[2]);
								if (dd < best)
								{
									best = dd;
									memcpy(y, pc, sizeof y);
								}
							}
						if (best < 4 * L * L)
							kd = 2;
						else
						{
							if (cnt == 1 || cnt == 7 || cnt == 2 || cnt == 6)
								dc = mr_project_crease(ref, pv, HM_NONE, &s, &l);
							if (dc < 2.25 * L * L && l.patch != HM_NONE)
								kd = 1;
							else
								mr_project(ref, pv, HM_NONE, &s, &l);
							mr_eval(ref, &l, y, NULL);
						}
						memcpy(cx + 3 * (size_t)nvx, y, sizeof y);
						kind[nvx] = (unsigned char)(4 + kd);   /* (4 and up: on the surface; 5 a crease's, 6 a corner's) */
						outer[v] = nvx++;
					}
				}
				if (ncell + 1 > capc)
				{
					unsigned int nc = qaws_internal_grow(capc, (size_t)ncell + 1, 256);
					unsigned int* gr = nc ? (unsigned int*)qaws_internal_realloc(cells, qaws_internal_mul3(sizeof(unsigned int), 8, nc)) : NULL;
					if (!gr)
						goto done;
					cells = gr;
					capc = nc;
				}
				cv = cells + 8 * (size_t)ncell++;
				for (q = 0; q < 4; q++)
				{
					cv[q] = cells[8 * (size_t)c + F[f][q]];
					cv[4 + q] = outer[cv[q]];
				}
			}
		}
	}
	/* smoothed: each free vertex to the centre of its neighbours (a surface vertex then back onto the surface), kept when
	   the worst scaled Jacobian round it rises */
	fixed = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nvx + 1, 1);
	if (!fixed || !hx_csr_build(&csr, nvx, cells, ncell))
		goto done;
	for (i = 0; i < nvx; i++)
		fixed[i] = kind[i] >= 5;    /* (the surface's corners and crease vertices stay) */
	for (rep.iterations = 0; rep.iterations < 10; rep.iterations++)
		for (i = 0; i < nvx; i++)
		{
			double keep[3], cen[3] = { 0, 0, 0 }, q0, q1, y[3];
			unsigned int j, n = 0;
			if (fixed[i] || csr.v0[i] == csr.v0[i + 1])
				continue;
			q0 = hx_worst(cx, cells, &csr, i);
			if (q0 > 0.8)
				continue;
			/* the neighbours along the cells' edges: those sharing a cell and an axis with it (one coordinate apart) */
			for (j = csr.v0[i]; j < csr.v0[i + 1]; j++)
			{
				unsigned int const* cv = cells + 8 * (size_t)csr.cell[j];
				unsigned int o;
				static unsigned char const nb3[8][3] = { { 1, 3, 4 }, { 2, 0, 5 }, { 3, 1, 6 }, { 0, 2, 7 }, { 7, 5, 0 }, { 4, 6, 1 }, { 5, 7, 2 }, { 6, 4, 3 } };
				for (o = 0; o < 8 && cv[o] != i; o++) {}
				for (k = 0; k < 3 && o < 8; k++)
				{
					unsigned int w = cv[nb3[o][k]];
					if (kind[i] >= 4 && kind[w] < 4)
						continue;   /* (a surface vertex: its surface neighbours only) */
					cen[0] += cx[3 * (size_t)w];
					cen[1] += cx[3 * (size_t)w + 1];
					cen[2] += cx[3 * (size_t)w + 2];
					n++;
				}
			}
			if (!n)
				continue;
			memcpy(keep, cx + 3 * (size_t)i, sizeof keep);
			for (k = 0; k < 3; k++)
				y[k] = 0.5 * (keep[k] + cen[k] / n);
			if (kind[i] >= 4)
			{
				qaws_mesh_location l;
				mr_project(ref, y, HM_NONE, &s, &l);
				mr_eval(ref, &l, y, NULL);
			}
			memcpy(cx + 3 * (size_t)i, y, sizeof y);
			q1 = hx_worst(cx, cells, &csr, i);
			if (q1 > q0)
				rep.moves++;
			else
				memcpy(cx + 3 * (size_t)i, keep, sizeof keep);
		}
	{
		unsigned char* types = (unsigned char*)qaws_internal_scratch_malloc((size_t)ncell + 1);
		if (!types)
			goto done;
		memset(types, QAWS_CELL_HEXAHEDRON, ncell);
		st = vm_create(cx, nvx, cells, types, ncell, out_mesh);
		qaws_internal_free(types);
	}
	if (st == QAWS_STATUS_OK)
	{
		qaws_mesh_3d* b = NULL;
		if (qaws_volume_mesh_get_boundary(*out_mesh, &b) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(b, &rep.quality);
		qaws_mesh_3d_destroy(b);
	}
done:
	qaws_internal_free(x); qaws_internal_free(tri); qaws_internal_free(vid); qaws_internal_free(cells); qaws_internal_free(outer); qaws_internal_free(cx); qaws_internal_free(fixed); qaws_internal_free(kind); qaws_internal_free(corners);
	qaws_internal_free(g.in); qaws_internal_free(g.keep);
	qaws_internal_free(csr.v0); qaws_internal_free(csr.cell);
	mr_scratch_free(&s);
	qaws_mesh_reference_destroy(ref);
	if (st != QAWS_STATUS_OK)
	{
		qaws_volume_mesh_destroy(*out_mesh);
		*out_mesh = NULL;
		return st;
	}
	if (out_report)
		*out_report = rep;
	return QAWS_STATUS_OK;
}
