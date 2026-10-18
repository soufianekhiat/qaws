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
#include "../qaws_hexahedralize.h"
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
	double* cross = (double*)malloc(sizeof(double) * cap);
	int ok = 0;
	l0 = (unsigned int*)calloc((size_t)nl + 1, sizeof(unsigned int));
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
	lt = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)l0[nl] + 1));
	fill = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)nl + 1));
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
					double* gr = (double*)realloc(cross, sizeof(double) * 2 * cap);
					if (!gr)
						goto done;
					cross = gr;
					cap *= 2;
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
	free(cross); free(l0); free(lt); free(fill);
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
	s->v0 = (unsigned int*)calloc((size_t)nv + 1, sizeof(unsigned int));
	s->cell = (unsigned int*)malloc(sizeof(unsigned int) * (8 * (size_t)nc + 1));
	fill = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)nv + 1));
	if (!s->v0 || !s->cell || !fill)
	{
		free(fill);
		return 0;
	}
	for (i = 0; i < 8 * nc; i++)
		s->v0[cells[i] + 1]++;
	for (i = 0; i < nv; i++)
		s->v0[i + 1] += s->v0[i];
	memcpy(fill, s->v0, sizeof(unsigned int) * nv);
	for (i = 0; i < 8 * nc; i++)
		s->cell[fill[cells[i]]++] = i / 8;
	free(fill);
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
	if (desc->method != QAWS_HEXES_GRID)
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
	tri = (unsigned int*)malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	x = (double*)malloc(sizeof(double) * (3 * (size_t)hm->nv + 1));
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
	g.in = (unsigned char*)calloc(ngv, 1);
	g.keep = (unsigned char*)calloc((size_t)g.n[0] * g.n[1] * g.n[2] + 1, 1);
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
	vid = (unsigned int*)malloc(sizeof(unsigned int) * ngv);
	if (!vid)
		goto done;
	for (i = 0; i < ngv; i++)
		vid[i] = HM_NONE;
	{
		unsigned int a, b, d, capv = 1024;
		cx = (double*)malloc(sizeof(double) * 3 * capv);
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
						unsigned int* gr = (unsigned int*)realloc(cells, sizeof(unsigned int) * 8 * (capc ? 2 * capc : 256));
						if (!gr)
							goto done;
						cells = gr;
						capc = capc ? 2 * capc : 256;
					}
					cv = cells + 8 * (size_t)ncell++;
					for (o = 0; o < 8; o++)
					{
						size_t gv = hx_v(&g, a + (unsigned int)CO[o][0], b + (unsigned int)CO[o][1], d + (unsigned int)CO[o][2]);
						if (vid[gv] == HM_NONE)
						{
							if (nvx == capv)
							{
								double* gr = (double*)realloc(cx, sizeof(double) * 6 * capv);
								if (!gr)
									goto done;
								cx = gr;
								capv *= 2;
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
			double* gr = (double*)realloc(cx, sizeof(double) * 6 * ((size_t)nvx + 1));
			if (!gr)
				goto done;
			cx = gr;
		}
		outer = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)nvx + 1));
		kind = (unsigned char*)calloc(2 * (size_t)nvx + 1, 1);
		if (!outer || !kind)
			goto done;
		for (i = 0; i < 2 * nvx; i++)
			outer[i] = HM_NONE;
		for (k = 0; k < ref->mesh->hm.nv; k++)
			if (ref->mesh->hm.vflag[k] & HM_CORNER)
				ncorner++;
		corners = (unsigned int*)malloc(sizeof(unsigned int) * ((size_t)ncorner + 1));
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
					unsigned int* gr = (unsigned int*)realloc(cells, sizeof(unsigned int) * 8 * 2 * capc);
					if (!gr)
						goto done;
					cells = gr;
					capc *= 2;
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
	fixed = (unsigned char*)calloc((size_t)nvx + 1, 1);
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
		unsigned char* types = (unsigned char*)malloc((size_t)ncell + 1);
		if (!types)
			goto done;
		memset(types, QAWS_CELL_HEXAHEDRON, ncell);
		st = vm_create(cx, nvx, cells, types, ncell, out_mesh);
		free(types);
	}
	if (st == QAWS_STATUS_OK)
	{
		qaws_mesh_3d* b = NULL;
		if (qaws_volume_mesh_get_boundary(*out_mesh, &b) == QAWS_STATUS_OK)
			qaws_mesh_3d_measure_quality(b, &rep.quality);
		qaws_mesh_3d_destroy(b);
	}
done:
	free(x); free(tri); free(vid); free(cells); free(outer); free(cx); free(fixed); free(kind); free(corners);
	free(g.in); free(g.keep);
	free(csr.v0); free(csr.cell);
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
