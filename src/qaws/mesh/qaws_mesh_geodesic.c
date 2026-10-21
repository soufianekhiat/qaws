/*
 * Geodesics on meshes (qaws_mesh_geodesic.h): the heat method on the
 * cotangent Laplacian, its matrices factored once by the sparse LDL^T;
 * paths down the distance's gradient across the triangles.
 */

#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_ldlt.h"
#include "../qaws_mesh_geodesic.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define GD_NONE 0xFFFFFFFFu

struct qaws_mesh_geodesic
{
	unsigned int nv, nt, nparts;
	double* x;                        /* 3 per vertex */
	unsigned int* t;                  /* 3 per triangle */
	unsigned int* adj;                /* per triangle edge (corner k to k + 1): the neighbour's edge 3 g + j, GD_NONE on the border */
	unsigned int* it;                 /* the intrinsic Delaunay triangles (3 each, as many as t) the heat and Poisson steps use */
	double* il;                       /* per intrinsic triangle edge (corner k to k + 1): its length */
	double* cot;                      /* per intrinsic corner: the cotangent of its angle (0 on a triangle without area) */
	double* mass;                     /* per vertex: a third of its intrinsic triangles' areas */
	unsigned int* part;               /* per vertex: its connected part (GD_NONE: on no triangle) */
	unsigned int *vt_start, *vt;      /* the corners (3 t + k) round each vertex */
	unsigned int *hrow, *drow, *prow; /* per vertex its row in the heat, Dirichlet heat and Poisson systems (GD_NONE: none) */
	unsigned int nh, nd, np;
	int dirichlet;                    /* the Dirichlet heat is factored (a border, and asked for) */
	int neumann;                      /* the Neumann heat is used */
	qaws_internal_ldlt heat, heat_d, poisson;
};

static void gd_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static double gd_dist2c(double const* a, double const* b)
{
	double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
	return x * x + y * y + z * z;
}

static double gd_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void gd_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static int gd_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

static unsigned int gd_find(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
	{
		par[i] = par[par[i]];
		i = par[i];
	}
	return i;
}

/* the gradients of triangle f's barycentric coordinates (3 each); 0 when it has no area */
static int gd_bary_grad(qaws_mesh_geodesic const* g, unsigned int f, double (*gl)[3])
{
	double const *x0 = g->x + 3 * (size_t)g->t[3 * (size_t)f], *x1 = g->x + 3 * (size_t)g->t[3 * (size_t)f + 1];
	double const* x2 = g->x + 3 * (size_t)g->t[3 * (size_t)f + 2];
	double u[3], v[3], n[3], e[3], nn;
	unsigned int k;
	gd_sub(x1, x0, u);
	gd_sub(x2, x0, v);
	gd_cross(u, v, n);
	nn = gd_dot(n, n);
	if (!(nn > 0))
		return 0;
	/* grad l_k = n x e_k / |n|^2, e_k the edge opposite corner k turning counter-clockwise */
	for (k = 0; k < 3; k++)
	{
		double const *a = g->x + 3 * (size_t)g->t[3 * (size_t)f + (k + 1) % 3], *b = g->x + 3 * (size_t)g->t[3 * (size_t)f + (k + 2) % 3];
		gd_sub(b, a, e);
		gd_cross(n, e, gl[k]);
		gl[k][0] /= nn;
		gl[k][1] /= nn;
		gl[k][2] /= nn;
	}
	return 1;
}

/* intrinsic triangle f laid out in the plane from its edge lengths, counter-clockwise; its area */
static double gd_layout(qaws_mesh_geodesic const* g, unsigned int f, double (*P)[2])
{
	double l0 = g->il[3 * (size_t)f], l1 = g->il[3 * (size_t)f + 1], l2 = g->il[3 * (size_t)f + 2], x;
	P[0][0] = P[0][1] = 0;
	P[1][0] = l0;
	P[1][1] = 0;
	x = l0 > 0 ? (l0 * l0 + l2 * l2 - l1 * l1) / (2 * l0) : 0;
	P[2][0] = x;
	P[2][1] = l2 * l2 - x * x > 0 ? sqrt(l2 * l2 - x * x) : 0;
	return 0.5 * l0 * P[2][1];
}

/* the cotangent at corner c of intrinsic triangle f */
static double gd_cot_at(qaws_mesh_geodesic const* g, unsigned int f, unsigned int c)
{
	double P[3][2], a = gd_layout(g, f, P), u[2], w[2];
	unsigned int k;
	if (!(a > 0))
		return 0;
	for (k = 0; k < 2; k++)
	{
		u[k] = P[(c + 1) % 3][k] - P[c][k];
		w[k] = P[(c + 2) % 3][k] - P[c][k];
	}
	return (u[0] * w[0] + u[1] * w[1]) / (2 * a);
}

static int gd_cmp_edge_w(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* The intrinsic Delaunay triangulation (Bobenko, Springborn 2007; the robust heat method of Sharp, Crane 2020): the mesh's
   edges flipped by their lengths alone, the surface unchanged, until the angles facing each edge sum to at most pi. Its
   Laplacian's weights are then non-negative, so the heat stays positive however stretched or obtuse the triangles are.
   Then the corners' cotangents, the masses and the edges' weights (ek, ew; ne of them) from it. */
static int gd_intrinsic(qaws_mesh_geodesic* g, unsigned long long** out_ek, double** out_ew, unsigned int* out_ne)
{
	unsigned int n = 3 * g->nt, i, j, k, f, top = 0, flips = 0, guard = 64 * g->nt + 64, ne = 0;
	unsigned int *adj = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1)), *stack;
	unsigned char* queued = (unsigned char*)qaws_internal_zalloc((size_t)n + 1, 1);
	unsigned long long* he;
	g->it = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	g->il = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)n + 1));
	stack = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	if (!adj || !queued || !g->it || !g->il || !stack)
	{
		qaws_internal_free(adj); qaws_internal_free(queued); qaws_internal_free(stack);
		return 0;
	}
	memcpy(g->it, g->t, sizeof(unsigned int) * n);
	memcpy(adj, g->adj, sizeof(unsigned int) * n);
	for (i = 0; i < n; i++)
	{
		g->il[i] = sqrt(gd_dist2c(g->x + 3 * (size_t)g->t[i], g->x + 3 * (size_t)g->t[3 * (i / 3) + (i % 3 + 1) % 3]));
		if (adj[i] != GD_NONE && adj[i] > i)
		{
			stack[top++] = i;
			queued[i] = 1;
		}
	}
	while (top && flips < guard)
	{
		unsigned int e = stack[--top], o = adj[e], ff, gg, kk, jj, A, B, C, D, nbc, nca, nad, ndb, q[4];
		double lAB, lBC, lCA, lAD, lDB, cx, cy, dx, dy, lCD;
		queued[e] = 0;
		if (o == GD_NONE)
			continue;
		ff = e / 3;
		kk = e % 3;
		gg = o / 3;
		jj = o % 3;
		if (gd_cot_at(g, ff, (kk + 2) % 3) + gd_cot_at(g, gg, (jj + 2) % 3) >= -1e-12)
			continue;
		A = g->it[e];
		B = g->it[3 * ff + (kk + 1) % 3];
		C = g->it[3 * ff + (kk + 2) % 3];
		D = g->it[3 * gg + (jj + 2) % 3];
		nbc = adj[3 * ff + (kk + 1) % 3];
		nca = adj[3 * ff + (kk + 2) % 3];
		nad = adj[3 * gg + (jj + 1) % 3];
		ndb = adj[3 * gg + (jj + 2) % 3];
		/* (a quad folded onto itself, round a vertex of two or three triangles: left as it is) */
		if (C == D || (nbc != GD_NONE && (nbc / 3 == ff || nbc / 3 == gg)) || (nca != GD_NONE && (nca / 3 == ff || nca / 3 == gg))
			|| (nad != GD_NONE && (nad / 3 == ff || nad / 3 == gg)) || (ndb != GD_NONE && (ndb / 3 == ff || ndb / 3 == gg)))
			continue;
		lAB = g->il[e];
		lBC = g->il[3 * ff + (kk + 1) % 3];
		lCA = g->il[3 * ff + (kk + 2) % 3];
		lAD = g->il[3 * gg + (jj + 1) % 3];
		lDB = g->il[3 * gg + (jj + 2) % 3];
		/* the quad unfolded: A at 0, B on the x axis, C above, D below */
		cx = (lAB * lAB + lCA * lCA - lBC * lBC) / (2 * lAB);
		cy = lCA * lCA - cx * cx > 0 ? sqrt(lCA * lCA - cx * cx) : 0;
		dx = (lAB * lAB + lAD * lAD - lDB * lDB) / (2 * lAB);
		dy = lAD * lAD - dx * dx > 0 ? -sqrt(lAD * lAD - dx * dx) : 0;
		lCD = sqrt((cx - dx) * (cx - dx) + (cy - dy) * (cy - dy));
		if (!(lCD > 0) || !(lAB > 0))
			continue;
		/* ff = A D C, gg = D B C */
		g->it[3 * ff] = A; g->it[3 * ff + 1] = D; g->it[3 * ff + 2] = C;
		g->it[3 * gg] = D; g->it[3 * gg + 1] = B; g->it[3 * gg + 2] = C;
		g->il[3 * ff] = lAD; g->il[3 * ff + 1] = lCD; g->il[3 * ff + 2] = lCA;
		g->il[3 * gg] = lDB; g->il[3 * gg + 1] = lBC; g->il[3 * gg + 2] = lCD;
		adj[3 * ff] = nad; adj[3 * ff + 1] = 3 * gg + 2; adj[3 * ff + 2] = nca;
		adj[3 * gg] = ndb; adj[3 * gg + 1] = nbc; adj[3 * gg + 2] = 3 * ff + 1;
		if (nad != GD_NONE) adj[nad] = 3 * ff;
		if (nca != GD_NONE) adj[nca] = 3 * ff + 2;
		if (ndb != GD_NONE) adj[ndb] = 3 * gg;
		if (nbc != GD_NONE) adj[nbc] = 3 * gg + 1;
		q[0] = 3 * ff;
		q[1] = 3 * ff + 2;
		q[2] = 3 * gg;
		q[3] = 3 * gg + 1;
		for (k = 0; k < 4; k++)
		{
			/* (an edge queued once, by its lower half) */
			unsigned int x = adj[q[k]] != GD_NONE && adj[q[k]] < q[k] ? adj[q[k]] : q[k];
			if (adj[x] != GD_NONE && !queued[x])
			{
				queued[x] = 1;
				stack[top++] = x;
			}
		}
		flips++;
	}
	qaws_internal_free(adj);
	qaws_internal_free(queued);
	qaws_internal_free(stack);
	/* cotangents, masses, and each pair of vertices' weight (half the cotangents facing its edges) */
	memset(g->mass, 0, sizeof(double) * g->nv);
	for (f = 0; f < g->nt; f++)
	{
		double P[3][2], a = gd_layout(g, f, P);
		for (k = 0; k < 3; k++)
		{
			g->cot[3 * (size_t)f + k] = gd_cot_at(g, f, k);
			g->mass[g->it[3 * (size_t)f + k]] += a / 3;
		}
	}
	he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(2 * sizeof(unsigned long long), 3, (size_t)g->nt + 1));
	*out_ek = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 3, (size_t)g->nt + 1));
	*out_ew = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)g->nt + 1));
	if (!he || !*out_ek || !*out_ew)
	{
		qaws_internal_free(he);
		return 0;
	}
	for (i = 0; i < n; i++)
	{
		unsigned long long p = g->it[i], q = g->it[3 * (i / 3) + (i % 3 + 1) % 3];
		double w = 0.5 * g->cot[3 * (i / 3) + (i % 3 + 2) % 3];
		he[2 * (size_t)i] = p < q ? (p << 32 | q) : (q << 32 | p);
		memcpy(he + 2 * (size_t)i + 1, &w, sizeof w);
	}
	qsort(he, n, 2 * sizeof(unsigned long long), gd_cmp_edge_w);
	for (i = 0; i < n; i = j)
	{
		double w = 0;
		for (j = i; j < n && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
		{
			double wj;
			memcpy(&wj, he + 2 * (size_t)j + 1, sizeof wj);
			w += wj;
		}
		if ((unsigned int)(he[2 * (size_t)i] >> 32) == (unsigned int)he[2 * (size_t)i])
			continue;   /* (an edge from a vertex to itself weighs nothing) */
		(*out_ek)[ne] = he[2 * (size_t)i];
		(*out_ew)[ne++] = w;
	}
	qaws_internal_free(he);
	*out_ne = ne;
	return 1;
}

/* rows of a symmetric matrix over the vertices with row[v] != GD_NONE: diagonal diag[v] + s K_vv, off-diagonal s K_vw;
   CSR with sorted columns, positions for the ordering */
static qaws_status gd_factor(qaws_internal_ldlt* f, qaws_mesh_geodesic const* g, unsigned long long const* ek, double const* ew,
	unsigned int ne, unsigned int const* row, unsigned int n, double const* diag, double s)
{
	unsigned int *rp = NULL, *col = NULL, *fill = NULL, i, v;
	double *val = NULL, *xyz = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(f, 0, sizeof *f);
	if (!n)
		return QAWS_STATUS_OK;
	rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 2, sizeof(unsigned int));
	fill = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	xyz = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)n + 1));
	if (!rp || !fill || !xyz)
		goto done;
	for (i = 0; i < n; i++)
		rp[i + 1] = 1;
	for (i = 0; i < ne; i++)
	{
		unsigned int a = row[(unsigned int)(ek[i] >> 32)], b = row[(unsigned int)ek[i]];
		if (a != GD_NONE && b != GD_NONE)
		{
			rp[a + 1]++;
			rp[b + 1]++;
		}
	}
	for (i = 0; i < n; i++)
		rp[i + 1] += rp[i];
	col = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)rp[n] + 1));
	val = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)rp[n] + 1));
	if (!col || !val)
		goto done;
	for (v = 0; v < g->nv; v++)
		if (row[v] != GD_NONE)
		{
			unsigned int r = row[v];
			fill[r] = rp[r] + 1;
			col[rp[r]] = r;
			val[rp[r]] = diag[v];
			memcpy(xyz + 3 * (size_t)r, g->x + 3 * (size_t)v, 3 * sizeof(double));
		}
	for (i = 0; i < ne; i++)
	{
		unsigned int va = (unsigned int)(ek[i] >> 32), vb = (unsigned int)ek[i], a = row[va], b = row[vb];
		if (a == GD_NONE || b == GD_NONE)
		{
			/* (a row next to a held vertex keeps the held side's weight on its diagonal: K's own) */
			continue;
		}
		col[fill[a]] = b;
		val[fill[a]++] = -s * ew[i];
		col[fill[b]] = a;
		val[fill[b]++] = -s * ew[i];
	}
	/* the diagonals: diag + s times the sum of the weights, held neighbours included */
	for (i = 0; i < ne; i++)
	{
		unsigned int a = row[(unsigned int)(ek[i] >> 32)], b = row[(unsigned int)ek[i]];
		if (a != GD_NONE)
			val[rp[a]] += s * ew[i];
		if (b != GD_NONE)
			val[rp[b]] += s * ew[i];
	}
	/* columns sorted within each row (a few entries each) */
	for (i = 0; i < n; i++)
	{
		unsigned int p, q;
		for (p = rp[i] + 1; p < rp[i + 1]; p++)
		{
			unsigned int c = col[p];
			double w = val[p];
			for (q = p; q > rp[i] && col[q - 1] > c; q--)
			{
				col[q] = col[q - 1];
				val[q] = val[q - 1];
			}
			col[q] = c;
			val[q] = w;
		}
	}
	st = qaws_internal_ldlt_factor(f, rp, col, val, xyz, n, 0);
done:
	qaws_internal_free(rp); qaws_internal_free(col); qaws_internal_free(fill); qaws_internal_free(val); qaws_internal_free(xyz);
	return st;
}

void qaws_mesh_geodesic_destroy(qaws_mesh_geodesic* g)
{
	if (!g)
		return;
	qaws_internal_ldlt_free(&g->heat);
	qaws_internal_ldlt_free(&g->heat_d);
	qaws_internal_ldlt_free(&g->poisson);
	qaws_internal_free(g->it); qaws_internal_free(g->il);
	qaws_internal_free(g->x); qaws_internal_free(g->t); qaws_internal_free(g->adj); qaws_internal_free(g->cot); qaws_internal_free(g->mass);
	qaws_internal_free(g->part); qaws_internal_free(g->vt_start); qaws_internal_free(g->vt); qaws_internal_free(g->hrow);
	qaws_internal_free(g->drow); qaws_internal_free(g->prow);
	qaws_internal_free(g);
}

qaws_status qaws_mesh_geodesic_create(qaws_mesh_geodesic** out, qaws_mesh_3d const* mesh, qaws_mesh_geodesic_desc const* desc)
{
	qaws_mesh_geodesic* g;
	qaws_scalar* pos = NULL;
	unsigned long long *he = NULL, *ek = NULL;
	double *ew = NULL, *zero = NULL, h = 0, area = 0, factor;
	unsigned int *par = NULL, *pin = NULL, nh = 0, ne = 0, i, j, k, f, v;
	unsigned char* onborder = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh || (desc && (!(desc->time >= 0) || (unsigned int)desc->boundary > (unsigned int)QAWS_GEODESIC_BOUNDARY_DIRICHLET)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	factor = desc && desc->time > 0 ? (double)desc->time : 1;
	g = (qaws_mesh_geodesic*)qaws_internal_zalloc(1, sizeof *g);
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	g->nv = qaws_mesh_3d_get_vertex_count(mesh);
	g->nt = qaws_mesh_3d_get_triangle_count(mesh);
	g->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)g->nv + 1));
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)g->nv + 1));
	g->t = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)g->nt + 1));
	g->adj = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)g->nt + 1));
	g->cot = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)g->nt + 1));
	g->mass = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	g->part = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	g->vt_start = (unsigned int*)qaws_internal_zalloc((size_t)g->nv + 2, sizeof(unsigned int));
	g->vt = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)g->nt + 1));
	g->hrow = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	g->drow = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	g->prow = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(2 * sizeof(unsigned long long), 3, (size_t)g->nt + 1));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	onborder = (unsigned char*)qaws_internal_zalloc((size_t)g->nv + 1, 1);
	if (!g->x || !pos || !g->t || !g->adj || !g->cot || !g->mass || !g->part || !g->vt_start || !g->vt || !g->hrow || !g->drow || !g->prow
		|| !he || !par || !onborder)
		goto fail;
	if ((st = qaws_mesh_3d_get_positions(pos, mesh)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_get_triangles(g->t, mesh)) != QAWS_STATUS_OK)
		goto fail;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * g->nv; i++)
		g->x[i] = (double)pos[i];
	/* corners' cotangents, the lumped masses, the mean edge length */
	for (f = 0; f < g->nt; f++)
	{
		double const* X[3];
		double n[3], u[3], w[3], a2;
		for (k = 0; k < 3; k++)
			X[k] = g->x + 3 * (size_t)g->t[3 * (size_t)f + k];
		gd_sub(X[1], X[0], u);
		gd_sub(X[2], X[0], w);
		gd_cross(u, w, n);
		a2 = sqrt(gd_dot(n, n));
		for (k = 0; k < 3; k++)
		{
			double e[3];
			gd_sub(X[(k + 1) % 3], X[k], e);
			h += sqrt(gd_dot(e, e));
			gd_sub(X[(k + 1) % 3], X[k], u);
			gd_sub(X[(k + 2) % 3], X[k], w);
			g->cot[3 * (size_t)f + k] = a2 > 0 ? gd_dot(u, w) / a2 : 0;
			g->mass[g->t[3 * (size_t)f + k]] += a2 / 6;
		}
		area += a2 / 2;
	}
	if (!(area > 0))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto fail;
	}
	h /= 3.0 * g->nt;
	/* half-edges by their edge (key, then the corner): neighbours across edges of exactly two triangles turned opposite ways */
	for (f = 0; f < g->nt; f++)
		for (k = 0; k < 3; k++)
		{
			unsigned long long p = g->t[3 * (size_t)f + k], q = g->t[3 * (size_t)f + (k + 1) % 3];
			he[2 * (3 * (size_t)f + k)] = p < q ? (p << 32 | q) : (q << 32 | p);
			he[2 * (3 * (size_t)f + k) + 1] = 3 * (unsigned long long)f + k;
			g->adj[3 * (size_t)f + k] = GD_NONE;
		}
	qsort(he, 3 * (size_t)g->nt, 2 * sizeof(unsigned long long), gd_cmp_u64);
	ek = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 3, (size_t)g->nt + 1));
	ew = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)g->nt + 1));
	if (!ek || !ew)
		goto fail;
	for (i = 0; i < 3 * g->nt; i = j)
	{
		double w = 0;
		for (j = i; j < 3 * g->nt && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
		{
			unsigned int c = (unsigned int)he[2 * (size_t)j + 1];
			w += 0.5 * g->cot[3 * (c / 3) + (c % 3 + 2) % 3];   /* (the corner opposite the edge) */
		}
		if (j == i + 1)
		{
			onborder[(unsigned int)(he[2 * (size_t)i] >> 32)] = 1;
			onborder[(unsigned int)he[2 * (size_t)i]] = 1;
		}
		if (j == i + 2)
		{
			unsigned int c0 = (unsigned int)he[2 * (size_t)i + 1], c1 = (unsigned int)he[2 * (size_t)i + 3];
			if (g->t[c0] != g->t[c1])
			{
				g->adj[c0] = c1;
				g->adj[c1] = c0;
			}
		}
		ek[ne] = he[2 * (size_t)i];
		ew[ne++] = w;
	}
	qaws_internal_free(he);
	he = NULL;
	/* the corners round each vertex, the parts */
	for (i = 0; i < 3 * g->nt; i++)
		g->vt_start[g->t[i] + 2]++;
	for (v = 0; v < g->nv; v++)
		g->vt_start[v + 2] += g->vt_start[v + 1];
	for (i = 0; i < 3 * g->nt; i++)
		g->vt[g->vt_start[g->t[i] + 1]++] = i;
	for (v = 0; v < g->nv; v++)
		par[v] = v;
	for (i = 0; i < ne; i++)
	{
		unsigned int a = gd_find(par, (unsigned int)(ek[i] >> 32)), b = gd_find(par, (unsigned int)ek[i]);
		if (a != b)
			par[a < b ? b : a] = a < b ? a : b;
	}
	{
		/* (the parts and the border from the mesh's own edges; the systems from the intrinsic Delaunay ones) */
		unsigned long long* iek = NULL;
		double* iew = NULL;
		unsigned int ine = 0;
		if (!gd_intrinsic(g, &iek, &iew, &ine))
		{
			qaws_internal_free(iek);
			qaws_internal_free(iew);
			goto fail;
		}
		qaws_internal_free(ek);
		qaws_internal_free(ew);
		ek = iek;
		ew = iew;
		ne = ine;
	}
	pin = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)g->nv + 1));
	zero = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	if (!pin || !zero)
		goto fail;
	for (v = 0; v < g->nv; v++)
	{
		g->part[v] = GD_NONE;
		pin[v] = GD_NONE;
	}
	for (v = 0; v < g->nv; v++)
		if (g->vt_start[v + 1] > g->vt_start[v])
		{
			unsigned int r = gd_find(par, v);
			if (pin[r] == GD_NONE)
				pin[r] = g->nparts++;
			g->part[v] = pin[r];
		}
	/* rows: the heat over the vertices on triangles; the Dirichlet heat without the border's; the Poisson system without
	   each part's first vertex */
	g->neumann = 1;
	g->dirichlet = 0;
	for (v = 0; v < g->nv; v++)
	{
		g->hrow[v] = g->part[v] != GD_NONE ? nh++ : GD_NONE;
		g->dirichlet |= g->part[v] != GD_NONE && onborder[v];
	}
	if (desc && desc->boundary == QAWS_GEODESIC_BOUNDARY_NEUMANN)
		g->dirichlet = 0;
	if (g->dirichlet && desc && desc->boundary == QAWS_GEODESIC_BOUNDARY_DIRICHLET)
		g->neumann = 0;
	g->nh = nh;
	for (v = 0; v < g->nv; v++)
		g->drow[v] = g->dirichlet && g->part[v] != GD_NONE && !onborder[v] ? g->nd++ : GD_NONE;
	{
		/* (each part's first vertex holds its Poisson solution at 0) */
		unsigned char* held = (unsigned char*)qaws_internal_zalloc((size_t)g->nparts + 1, 1);
		if (!held)
			goto fail;
		for (v = 0; v < g->nv; v++)
		{
			if (g->part[v] == GD_NONE || !held[g->part[v]])
			{
				if (g->part[v] != GD_NONE)
					held[g->part[v]] = 1;
				g->prow[v] = GD_NONE;
				continue;
			}
			g->prow[v] = g->np++;
		}
		qaws_internal_free(held);
	}
	{
		double t = factor * h * h;
		if (g->neumann && (st = gd_factor(&g->heat, g, ek, ew, ne, g->hrow, g->nh, g->mass, t)) != QAWS_STATUS_OK)
			goto fail;
		if (g->dirichlet && (st = gd_factor(&g->heat_d, g, ek, ew, ne, g->drow, g->nd, g->mass, t)) != QAWS_STATUS_OK)
			goto fail;
		if ((st = gd_factor(&g->poisson, g, ek, ew, ne, g->prow, g->np, zero, 1)) != QAWS_STATUS_OK)
			goto fail;
	}
	qaws_internal_free(pos); qaws_internal_free(ek); qaws_internal_free(ew); qaws_internal_free(par); qaws_internal_free(pin);
	qaws_internal_free(zero); qaws_internal_free(onborder);
	*out = g;
	return QAWS_STATUS_OK;
fail:
	qaws_internal_free(pos); qaws_internal_free(he); qaws_internal_free(ek); qaws_internal_free(ew); qaws_internal_free(par);
	qaws_internal_free(pin); qaws_internal_free(zero); qaws_internal_free(onborder);
	qaws_mesh_geodesic_destroy(g);
	return st;
}

static int gd_location_ok(qaws_mesh_geodesic const* g, qaws_mesh_location const* l)
{
	double a = (double)l->a, b = (double)l->b;
	return l->patch < g->nt && isfinite(a) && isfinite(b) && a >= -1e-6 && b >= -1e-6 && a + b <= 1 + 1e-6;
}

/* x = A^-1 b on a factor shared by threads: the solve's workspace its own */
static int gd_solve(qaws_internal_ldlt const* f, double* x)
{
	qaws_internal_ldlt local = *f;
	if (!f->n)
		return 1;
	local.w = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 2 * (size_t)f->n + 1));
	if (!local.w)
		return 0;
	qaws_internal_ldlt_solve(x, &local, x);
	qaws_internal_free(local.w);
	return 1;
}

/* the heat from the sources' unit masses (per vertex) through the system with rows row */
static int gd_heat(qaws_mesh_geodesic const* g, qaws_internal_ldlt const* f, unsigned int const* row, unsigned int n, double const* delta,
	double* out)
{
	double* b = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)n + 1));
	unsigned int v;
	if (!b)
		return 0;
	for (v = 0; v < g->nv; v++)
		if (row[v] != GD_NONE)
			b[row[v]] = delta[v];
	if (!gd_solve(f, b))
	{
		qaws_internal_free(b);
		return 0;
	}
	for (v = 0; v < g->nv; v++)
		out[v] = row[v] != GD_NONE ? b[row[v]] : 0;
	qaws_internal_free(b);
	return 1;
}

qaws_status qaws_mesh_geodesic_distance(qaws_scalar* out, qaws_mesh_geodesic const* g, qaws_mesh_location const* sources, unsigned int ns)
{
	double *delta = NULL, *u = NULL, *ud = NULL, *b = NULL, *sb = NULL, *sm = NULL, *lo = NULL;
	unsigned int i, k, f, v;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out || !g || !sources || !ns)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < ns; i++)
		if (!gd_location_ok(g, sources + i))
			return QAWS_STATUS_INVALID_ARGUMENT;
	delta = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	u = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	ud = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	b = (double*)qaws_internal_zalloc((size_t)g->nv + 1, sizeof(double));
	sb = (double*)qaws_internal_zalloc((size_t)g->nparts + 1, sizeof(double));
	sm = (double*)qaws_internal_zalloc((size_t)g->nparts + 1, sizeof(double));
	lo = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)g->nparts + 1));
	if (!delta || !u || !ud || !b || !sb || !sm || !lo)
		goto done;
	/* each source a unit of heat on its triangle's corners */
	for (i = 0; i < ns; i++)
	{
		unsigned int const* c = g->t + 3 * (size_t)sources[i].patch;
		double a = (double)sources[i].a, bb = (double)sources[i].b;
		delta[c[0]] += 1 - a - bb;
		delta[c[1]] += a;
		delta[c[2]] += bb;
	}
	if (g->neumann && !gd_heat(g, &g->heat, g->hrow, g->nh, delta, u))
		goto done;
	if (g->dirichlet && !gd_heat(g, &g->heat_d, g->drow, g->nd, delta, ud))
		goto done;
	for (v = 0; v < g->nv; v++)
		u[v] = g->neumann && g->dirichlet ? 0.5 * (u[v] + ud[v]) : g->neumann ? u[v] : ud[v];
	/* the heat's unit descent per intrinsic triangle (laid out in its plane), its integrated divergence at the corners */
	for (f = 0; f < g->nt; f++)
	{
		double P[3][2], X[2] = { 0, 0 }, l, a = gd_layout(g, f, P);
		unsigned int const* c = g->it + 3 * (size_t)f;
		if (!(a > 0))
			continue;
		for (k = 0; k < 3; k++)
		{
			/* grad l_k: the edge opposite k turned a quarter left, over twice the area */
			double ex = P[(k + 2) % 3][0] - P[(k + 1) % 3][0], ey = P[(k + 2) % 3][1] - P[(k + 1) % 3][1];
			X[0] += u[c[k]] * ey / (2 * a);
			X[1] -= u[c[k]] * ex / (2 * a);
		}
		l = sqrt(X[0] * X[0] + X[1] * X[1]);
		if (!(l > 0))
			continue;
		X[0] /= l;
		X[1] /= l;
		for (k = 0; k < 3; k++)
		{
			unsigned int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
			double e1 = (P[k1][0] - P[k][0]) * X[0] + (P[k1][1] - P[k][1]) * X[1], e2 = (P[k2][0] - P[k][0]) * X[0] + (P[k2][1] - P[k][1]) * X[1];
			/* cot of the angle opposite e1 (at k2) and opposite e2 (at k1) */
			b[c[k]] += 0.5 * (g->cot[3 * (size_t)f + k2] * e1 + g->cot[3 * (size_t)f + k1] * e2);
		}
	}
	/* Poisson: L phi = div X with L the cotangent Laplacian (-K); the divergence's mean taken out per part (on a closed part it
	   sums to 0, up to rounding; on a border, the flux out) so the held vertex takes no source */
	for (v = 0; v < g->nv; v++)
		if (g->part[v] != GD_NONE)
		{
			sb[g->part[v]] += b[v];
			sm[g->part[v]] += g->mass[v];
		}
	for (v = 0; v < g->nv; v++)
		if (g->part[v] != GD_NONE && sm[g->part[v]] > 0)
			b[v] -= g->mass[v] * sb[g->part[v]] / sm[g->part[v]];
	memset(u, 0, sizeof(double) * g->nv);
	{
		double* r = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)g->np + 1));
		if (!r)
			goto done;
		for (v = 0; v < g->nv; v++)
			if (g->prow[v] != GD_NONE)
				r[g->prow[v]] = -b[v];
		if (!gd_solve(&g->poisson, r))
		{
			qaws_internal_free(r);
			goto done;
		}
		for (v = 0; v < g->nv; v++)
			u[v] = g->prow[v] != GD_NONE ? r[g->prow[v]] : 0;
		qaws_internal_free(r);
	}
	/* each part shifted so its nearest source is at 0 */
	for (i = 0; i < g->nparts; i++)
		lo[i] = HUGE_VAL;
	for (i = 0; i < ns; i++)
	{
		unsigned int const* c = g->t + 3 * (size_t)sources[i].patch;
		double a = (double)sources[i].a, bb = (double)sources[i].b, p = (1 - a - bb) * u[c[0]] + a * u[c[1]] + bb * u[c[2]];
		unsigned int q = g->part[c[0]];
		lo[q] = p < lo[q] ? p : lo[q];
	}
	for (v = 0; v < g->nv; v++)
	{
		double d = g->part[v] != GD_NONE && lo[g->part[v]] < HUGE_VAL ? u[v] - lo[g->part[v]] : HUGE_VAL;
		out[v] = (qaws_scalar)(d > 0 ? d : 0);
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(delta); qaws_internal_free(u); qaws_internal_free(ud); qaws_internal_free(b); qaws_internal_free(sb);
	qaws_internal_free(sm); qaws_internal_free(lo);
	return st;
}

/* ---- paths -------------------------------------------------------------- */

/* the path as traced: fixed points (the start, vertices, the source) and between them the edges crossed */
typedef struct gd_node
{
	unsigned int f;                   /* the triangle the point is in (a crossing: the one it leaves) */
	unsigned int e;                   /* a crossing: f's edge; GD_NONE: a fixed point */
	double l[3];                      /* its barycentric coordinates in f */
} gd_node;

typedef struct gd_path
{
	qaws_mesh_geodesic const* g;
	qaws_scalar const* d;
	qaws_scalar* pts;
	qaws_mesh_location* loc;
	unsigned int n, cap;
	double last[3];
	gd_node* node;
	unsigned int nn, cn;
	int oom;
} gd_path;

static void gd_node_add(gd_path* P, unsigned int f, unsigned int e, double const* l)
{
	if (P->nn == P->cn)
	{
		unsigned int cap = P->cn ? 2 * P->cn : 64;
		gd_node* q = (gd_node*)qaws_internal_realloc(P->node, qaws_internal_mul2(sizeof(gd_node), cap));
		if (!q)
		{
			P->oom = 1;
			return;
		}
		P->node = q;
		P->cn = cap;
	}
	P->node[P->nn].f = f;
	P->node[P->nn].e = e;
	memcpy(P->node[P->nn++].l, l, 3 * sizeof(double));
}

/* a point of the path at barycentric l in triangle f (one equal to the last is not repeated) */
static void gd_emit(gd_path* P, unsigned int f, double const* l)
{
	unsigned int const* c = P->g->t + 3 * (size_t)f;
	double p[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
		p[k] = l[0] * P->g->x[3 * (size_t)c[0] + k] + l[1] * P->g->x[3 * (size_t)c[1] + k] + l[2] * P->g->x[3 * (size_t)c[2] + k];
	if (P->n && (p[0] - P->last[0]) * (p[0] - P->last[0]) + (p[1] - P->last[1]) * (p[1] - P->last[1]) + (p[2] - P->last[2]) * (p[2] - P->last[2])
		<= 0)
		return;
	memcpy(P->last, p, sizeof p);
	if (P->n < P->cap)
	{
		if (P->pts)
			for (k = 0; k < 3; k++)
				P->pts[3 * (size_t)P->n + k] = (qaws_scalar)p[k];
		if (P->loc)
		{
			P->loc[P->n].patch = f;
			P->loc[P->n].a = (qaws_scalar)l[1];
			P->loc[P->n].b = (qaws_scalar)l[2];
		}
	}
	P->n++;
}

/* the descent -grad d in triangle f, and its rates on f's barycentric coordinates; 0 without area */
static int gd_descent(gd_path const* P, unsigned int f, double* dir, double* dl)
{
	double gl[3][3];
	unsigned int k;
	if (!gd_bary_grad(P->g, f, gl))
		return 0;
	dir[0] = dir[1] = dir[2] = 0;
	for (k = 0; k < 3; k++)
	{
		double dk = (double)P->d[P->g->t[3 * (size_t)f + k]];
		dir[0] -= dk * gl[k][0];
		dir[1] -= dk * gl[k][1];
		dir[2] -= dk * gl[k][2];
	}
	for (k = 0; k < 3; k++)
		dl[k] = gd_dot(gl[k], dir);
	return 1;
}

/* the corner of vertex v in triangle f */
static unsigned int gd_corner(qaws_mesh_geodesic const* g, unsigned int f, unsigned int v)
{
	return g->t[3 * (size_t)f] == v ? 0 : g->t[3 * (size_t)f + 1] == v ? 1 : 2;
}

/* the path through vertex v from triangle fa to fb (both round v): the last node (v's, in fa) replaced by the edges round v
   crossed on the side of the smaller angle, so the run goes on through v (the shortening then pulls it off v where that
   side is under pi); 0 when neither side reaches fb (a border between) */
static int gd_fan(gd_path* P, unsigned int fa, unsigned int v, unsigned int fb)
{
	qaws_mesh_geodesic const* g = P->g;
	unsigned int n = g->vt_start[v + 1] - g->vt_start[v], side, best = GD_NONE, s, f;
	double ang[2] = { HUGE_VAL, HUGE_VAL };
	if (fa == fb)
	{
		P->nn--;
		return 1;
	}
	/* side 0 crosses the edge leaving v (corner c to c + 1), side 1 the edge into v (c + 2 to c) */
	for (side = 0; side < 2; side++)
	{
		double sum = 0;
		for (f = fa, s = 0; s <= n; s++)
		{
			unsigned int c = gd_corner(g, f, v), e = side ? (c + 2) % 3 : c, nb = g->adj[3 * (size_t)f + e];
			if (nb == GD_NONE)
				break;
			f = nb / 3;
			if (f == fb)
			{
				ang[side] = sum;
				break;
			}
			{
				double const* X = g->x + 3 * (size_t)v;
				double u[3], w[3], q[3];
				unsigned int cc = gd_corner(g, f, v);
				gd_sub(g->x + 3 * (size_t)g->t[3 * (size_t)f + (cc + 1) % 3], X, u);
				gd_sub(g->x + 3 * (size_t)g->t[3 * (size_t)f + (cc + 2) % 3], X, w);
				gd_cross(u, w, q);
				sum += atan2(sqrt(gd_dot(q, q)), gd_dot(u, w));
			}
		}
	}
	if (ang[0] < HUGE_VAL || ang[1] < HUGE_VAL)
		best = ang[0] <= ang[1] ? 0 : 1;
	if (best == GD_NONE)
		return 0;
	P->nn--;
	for (f = fa; f != fb && !P->oom;)
	{
		unsigned int c = gd_corner(g, f, v), e = best ? (c + 2) % 3 : c;
		double l[3] = { 0, 0, 0 };
		l[c] = 1;
		gd_node_add(P, f, e, l);
		f = g->adj[3 * (size_t)f + e] / 3;
	}
	return 1;
}

/* twice the signed area of a b c, negative when c is left of a -> b (the funnel's convention) */
static double gd_tri2(double const* a, double const* b, double const* c)
{
	return (c[0] - a[0]) * (b[1] - a[1]) - (b[0] - a[0]) * (c[1] - a[1]);
}

/* corner X of a triangle unfolded left of the directed edge L -> R, given its squared distances to L and R */
static void gd_unfold(double const* L, double const* R, double dl2, double dr2, double* X)
{
	double u[2] = { R[0] - L[0], R[1] - L[1] }, len = sqrt(u[0] * u[0] + u[1] * u[1]), a, b;
	if (!(len > 0))
	{
		X[0] = L[0];
		X[1] = L[1];
		return;
	}
	u[0] /= len;
	u[1] /= len;
	a = (dl2 - dr2 + len * len) / (2 * len);
	b = dl2 - a * a > 0 ? sqrt(dl2 - a * a) : 0;
	X[0] = L[0] + a * u[0] - b * u[1];
	X[1] = L[1] + a * u[1] + b * u[0];
}

static double gd_dist2(double const* a, double const* b)
{
	double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
	return x * x + y * y + z * z;
}

/* the run from fixed node a to fixed node z through the crossings between, shortened: its triangles unfolded into the plane
   and the shortest way through their shared edges taken (the funnel algorithm of Lee and Preparata, in Mononen's form); 0 when
   the run's triangles do not chain (its points are then kept as traced) */
static int gd_straighten(gd_path* P, unsigned int a, unsigned int z)
{
	qaws_mesh_geodesic const* g = P->g;
	unsigned int k = z - a - 1, i, j, f = P->node[a].f, *pf = NULL, *pe = NULL, *apex = NULL, na = 0;
	double (*Lp)[2] = NULL, (*Rp)[2] = NULL, c2[3][2], s2[2], e2[2];
	unsigned char* side = NULL;
	int ok = 0;
	if (P->node[a + 1].f != f)
		return 0;
	Lp = (double(*)[2])qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)k + 2));
	Rp = (double(*)[2])qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)k + 2));
	pf = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)k + 2));
	pe = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)k + 2));
	apex = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)k + 3));
	side = (unsigned char*)qaws_internal_malloc((size_t)k + 3);
	if (!Lp || !Rp || !pf || !pe || !apex || !side)
	{
		P->oom = 1;
		goto done;
	}
	{
		/* the first triangle in the plane, counter-clockwise */
		double const *x0 = g->x + 3 * (size_t)g->t[3 * (size_t)f], *x1 = g->x + 3 * (size_t)g->t[3 * (size_t)f + 1];
		double const* x2 = g->x + 3 * (size_t)g->t[3 * (size_t)f + 2];
		c2[0][0] = c2[0][1] = 0;
		c2[1][0] = sqrt(gd_dist2(x0, x1));
		c2[1][1] = 0;
		gd_unfold(c2[0], c2[1], gd_dist2(x0, x2), gd_dist2(x1, x2), c2[2]);
	}
	for (j = 0; j < 2; j++)
		s2[j] = P->node[a].l[0] * c2[0][j] + P->node[a].l[1] * c2[1][j] + P->node[a].l[2] * c2[2][j];
	/* portals: per crossing, the edge's ends left and right of the way across; the next triangle unfolded past it */
	for (i = 0; i < k; i++)
	{
		gd_node const* N = P->node + a + 1 + i;
		unsigned int e = N->e, nb, h, jn, cl, cr, x;
		double nc[3][2];
		if (N->f != f || (nb = g->adj[3 * (size_t)f + e]) == GD_NONE)
			goto done;
		cr = e;
		cl = (e + 1) % 3;
		memcpy(Lp[i + 1], c2[cl], sizeof Lp[i + 1]);
		memcpy(Rp[i + 1], c2[cr], sizeof Rp[i + 1]);
		pf[i + 1] = f;
		pe[i + 1] = e;
		h = nb / 3;
		jn = nb % 3;
		/* in h: corner jn is f's left end, jn + 1 its right one, jn + 2 the new corner, left of L -> R as h turns */
		x = g->t[3 * (size_t)h + (jn + 2) % 3];
		memcpy(nc[jn], c2[cl], sizeof nc[jn]);
		memcpy(nc[(jn + 1) % 3], c2[cr], sizeof nc[jn]);
		gd_unfold(c2[cl], c2[cr], gd_dist2(g->x + 3 * (size_t)x, g->x + 3 * (size_t)g->t[3 * (size_t)f + cl]),
			gd_dist2(g->x + 3 * (size_t)x, g->x + 3 * (size_t)g->t[3 * (size_t)f + cr]), nc[(jn + 2) % 3]);
		memcpy(c2, nc, sizeof c2);
		f = h;
	}
	if (P->node[z].f != f)
		goto done;
	for (j = 0; j < 2; j++)
		e2[j] = P->node[z].l[0] * c2[0][j] + P->node[z].l[1] * c2[1][j] + P->node[z].l[2] * c2[2][j];
	memcpy(Lp[0], s2, sizeof s2);
	memcpy(Rp[0], s2, sizeof s2);
	memcpy(Lp[k + 1], e2, sizeof e2);
	memcpy(Rp[k + 1], e2, sizeof e2);
	{
		/* the funnel: the apexes the shortest way bends at, each a portal's end */
		double ap[2], fl[2], fr[2];
		unsigned int ia = 0, il = 0, ir = 0;
		memcpy(ap, s2, sizeof ap);
		memcpy(fl, s2, sizeof fl);
		memcpy(fr, s2, sizeof fr);
		apex[na] = 0;
		side[na++] = 0;
		for (i = 1; i <= k + 1; i++)
		{
			double const *L = Lp[i], *R = Rp[i];
			if (gd_tri2(ap, fr, R) <= 0)
			{
				if ((ap[0] == fr[0] && ap[1] == fr[1]) || gd_tri2(ap, fl, R) > 0)
				{
					memcpy(fr, R, sizeof fr);
					ir = i;
				}
				else
				{
					/* the right side crosses the left: the left end is a bend */
					apex[na] = il;
					side[na++] = 1;
					memcpy(ap, fl, sizeof ap);
					ia = il;
					memcpy(fr, ap, sizeof fr);
					ir = ia;
					i = ia;
					continue;
				}
			}
			if (gd_tri2(ap, fl, L) >= 0)
			{
				if ((ap[0] == fl[0] && ap[1] == fl[1]) || gd_tri2(ap, fr, L) < 0)
				{
					memcpy(fl, L, sizeof fl);
					il = i;
				}
				else
				{
					apex[na] = ir;
					side[na++] = 2;
					memcpy(ap, fr, sizeof ap);
					ia = ir;
					memcpy(fl, ap, sizeof fl);
					il = ia;
					i = ia;
					continue;
				}
			}
		}
		apex[na] = k + 1;
		side[na++] = 0;
	}
	/* the points: each bend, and between bends the crossing of each edge by the straight way */
	gd_emit(P, P->node[a].f, P->node[a].l);
	for (i = 0; i + 1 < na; i++)
	{
		double A[2], B[2];
		unsigned int ja = apex[i], jb = apex[i + 1];
		memcpy(A, side[i] == 1 ? Lp[ja] : Rp[ja], sizeof A);
		memcpy(B, side[i + 1] == 1 ? Lp[jb] : Rp[jb], sizeof B);
		for (j = ja + 1; j < jb; j++)
		{
			/* t along the portal from its right end to its left one where A -> B meets it */
			double const *L = Lp[j], *R = Rp[j];
			double d1 = gd_tri2(A, B, R), d2 = gd_tri2(A, B, L), t = d1 - d2 != 0 ? d1 / (d1 - d2) : 0.5, l[3] = { 0, 0, 0 };
			t = t < 0 ? 0 : t > 1 ? 1 : t;
			l[pe[j]] = 1 - t;
			l[(pe[j] + 1) % 3] = t;
			gd_emit(P, pf[j], l);
		}
		if (side[i + 1] && jb > 0 && jb <= k)
		{
			double l[3] = { 0, 0, 0 };
			l[side[i + 1] == 1 ? (pe[jb] + 1) % 3 : pe[jb]] = 1;
			gd_emit(P, pf[jb], l);
		}
	}
	gd_emit(P, P->node[z].f, P->node[z].l);
	ok = 1;
done:
	qaws_internal_free(Lp); qaws_internal_free(Rp); qaws_internal_free(pf); qaws_internal_free(pe); qaws_internal_free(apex);
	qaws_internal_free(side);
	return ok;
}

qaws_status qaws_mesh_geodesic_path(qaws_scalar* out_points, qaws_mesh_location* out_locations, unsigned int* out_count,
	qaws_mesh_geodesic const* g, qaws_scalar const* d, qaws_mesh_location const* sources, unsigned int ns, qaws_mesh_location from,
	unsigned int capacity)
{
	gd_path P;
	double l[3];
	unsigned int f = from.patch, v = GD_NONE, i, k;
	unsigned long long steps = 0, guard;
	if (out_count)
		*out_count = 0;
	if (!out_count || !g || !d || !sources || !ns || !gd_location_ok(g, &from) || (capacity && !out_points && !out_locations))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < ns; i++)
		if (!gd_location_ok(g, sources + i))
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&P, 0, sizeof P);
	P.g = g;
	P.d = d;
	P.pts = out_points;
	P.loc = out_locations;
	P.cap = capacity;
	guard = 8ull * g->nt + 64;
	l[1] = (double)from.a;
	l[2] = (double)from.b;
	l[0] = 1 - l[1] - l[2];
	gd_node_add(&P, f, GD_NONE, l);
	while (steps++ < guard && !P.oom)
	{
		double dir[3], dl[3], s = HUGE_VAL;
		unsigned int ke = GD_NONE;
		if (v == GD_NONE)
		{
			/* in triangle f at l: straight to a source in it, else down to its border */
			for (i = 0; i < ns && sources[i].patch != f; i++)
				;
			if (i < ns)
			{
				double ls[3] = { 1 - (double)sources[i].a - (double)sources[i].b, (double)sources[i].a, (double)sources[i].b };
				gd_node_add(&P, f, GD_NONE, ls);
				break;
			}
			if (!gd_descent(&P, f, dir, dl) || !(gd_dot(dir, dir) > 0))
				break;
			for (k = 0; k < 3; k++)
				if (dl[k] < 0 && l[k] / -dl[k] < s)
				{
					s = l[k] / -dl[k];
					ke = k;
				}
			if (ke == GD_NONE)
				break;
			for (k = 0; k < 3; k++)
				l[k] = k == ke ? 0 : l[k] + s * dl[k] > 0 ? l[k] + s * dl[k] : 0;
			s = l[0] + l[1] + l[2];
			for (k = 0; k < 3; k++)
				l[k] /= s;
			for (k = 0; k < 3 && l[k] < 1 - 1e-9; k++)
				;
			if (k < 3)
			{
				/* (at a vertex) */
				l[0] = l[1] = l[2] = 0;
				l[k] = 1;
				gd_node_add(&P, f, GD_NONE, l);
				v = g->t[3 * (size_t)f + k];
				continue;
			}
			{
				/* across the edge opposite ke into the neighbour, or along the edge when both sides lead onto it */
				unsigned int e = (ke + 1) % 3, nb = g->adj[3 * (size_t)f + e], h, j;
				if (nb != GD_NONE)
				{
					double lh[3], dirh[3], dlh[3];
					h = nb / 3;
					j = nb % 3;
					lh[j] = l[(e + 1) % 3];
					lh[(j + 1) % 3] = l[e];
					lh[(j + 2) % 3] = 0;
					for (i = 0; i < ns && sources[i].patch != h; i++)
						;
					if (i < ns || (gd_descent(&P, h, dirh, dlh) && dlh[(j + 2) % 3] > 0))
					{
						/* (into a source's triangle whatever the descent: the path ends straight there) */
						gd_node_add(&P, f, e, l);
						f = h;
						memcpy(l, lh, sizeof l);
						continue;
					}
				}
				{
					unsigned int a = g->t[3 * (size_t)f + e], b = g->t[3 * (size_t)f + (e + 1) % 3];
					double ln[3] = { 0, 0, 0 };
					v = d[a] <= d[b] ? a : b;
					ln[v == a ? e : (e + 1) % 3] = 1;
					gd_node_add(&P, f, GD_NONE, ln);
				}
			}
			continue;
		}
		{
			/* at vertex v: on a source's triangle, straight to the source; else into the triangle round v whose descent leaves v
			   the steepest, or along the steepest edge down; nothing lower: the end */
			double best = 0, best_edge = 0;
			unsigned int bf = GD_NONE, bk = 0, bu = GD_NONE, bu_f = 0, bu_k = 0, c;
			for (i = 0; i < ns; i++)
			{
				unsigned int const* sc = g->t + 3 * (size_t)sources[i].patch;
				if (sc[0] == v || sc[1] == v || sc[2] == v)
					break;
			}
			if (i < ns)
			{
				double ls[3] = { 1 - (double)sources[i].a - (double)sources[i].b, (double)sources[i].a, (double)sources[i].b };
				gd_fan(&P, P.node[P.nn - 1].f, v, sources[i].patch);
				gd_node_add(&P, sources[i].patch, GD_NONE, ls);
				break;
			}
			for (c = g->vt_start[v]; c < g->vt_start[v + 1]; c++)
			{
				unsigned int h = g->vt[c] / 3, kc = g->vt[c] % 3, m;
				if (gd_descent(&P, h, dir, dl) && dl[kc] < 0 && dl[(kc + 1) % 3] >= 0 && dl[(kc + 2) % 3] >= 0)
				{
					double slope = sqrt(gd_dot(dir, dir));
					if (slope > best)
					{
						best = slope;
						bf = h;
						bk = kc;
					}
				}
				for (m = 1; m < 3; m++)
				{
					unsigned int w = g->t[3 * (size_t)h + (kc + m) % 3];
					double e[3], len, slope;
					gd_sub(g->x + 3 * (size_t)w, g->x + 3 * (size_t)v, e);
					len = sqrt(gd_dot(e, e));
					slope = len > 0 ? ((double)d[v] - (double)d[w]) / len : 0;
					if (slope > best_edge)
					{
						best_edge = slope;
						bu = w;
						bu_f = h;
						bu_k = (kc + m) % 3;
					}
				}
			}
			if (bf != GD_NONE && best >= best_edge)
			{
				/* (the run from v starts in bf) */
				f = bf;
				l[0] = l[1] = l[2] = 0;
				l[bk] = 1;
				if (!gd_fan(&P, P.node[P.nn - 1].f, v, bf))
				{
					P.node[P.nn - 1].f = bf;
					memcpy(P.node[P.nn - 1].l, l, sizeof l);
				}
				v = GD_NONE;
				continue;
			}
			if (bu == GD_NONE)
				break;
			{
				double ln[3] = { 0, 0, 0 };
				ln[bu_k] = 1;
				gd_fan(&P, P.node[P.nn - 1].f, v, bu_f);
				gd_node_add(&P, bu_f, GD_NONE, ln);
				v = bu;
			}
		}
	}
	/* each run between fixed points shortened, or its points as traced */
	for (i = 0; i < P.nn && !P.oom;)
	{
		unsigned int z = i + 1;
		while (z < P.nn && P.node[z].e != GD_NONE)
			z++;
		if (z < P.nn && z > i + 1 && gd_straighten(&P, i, z))
		{
			i = z;
			continue;
		}
		for (k = i; k < z && k < P.nn; k++)
			gd_emit(&P, P.node[k].f, P.node[k].l);
		if (z + 1 >= P.nn)
		{
			if (z < P.nn)
				gd_emit(&P, P.node[z].f, P.node[z].l);
			break;
		}
		i = z;
	}
	qaws_internal_free(P.node);
	if (P.oom)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	*out_count = P.n;
	return P.n > capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}
