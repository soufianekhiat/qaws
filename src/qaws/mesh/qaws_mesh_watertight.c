/*
 * Watertight surfaces (qaws_mesh_watertight.h): the winding number's inside
 * on a grid, dual contouring with Hermite data from rays against the input,
 * then qaws_mesh_repair to split non-manifold places.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_watertight.h"
#include "../qaws_mesh_bvh.h"
#include "../qaws_mesh_repair.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct wtg_grid
{
	qaws_mesh_bvh* b;
	double lo[3], h, off;
	unsigned int n[3];
} wtg_grid;

static size_t wtg_id(wtg_grid const* g, unsigned int i, unsigned int j, unsigned int k)
{
	return i + (size_t)g->n[0] * (j + (size_t)g->n[1] * k);
}

static void wtg_point(wtg_grid const* g, size_t id, double* p)
{
	unsigned int i = (unsigned int)(id % g->n[0]), j = (unsigned int)(id / g->n[0] % g->n[1]), k = (unsigned int)(id / ((size_t)g->n[0] * g->n[1]));
	p[0] = g->lo[0] + g->h * i;
	p[1] = g->lo[1] + g->h * j;
	p[2] = g->lo[2] + g->h * k;
}

/* inside: winding number at least 1/2, or within the offset of the input */
static int wtg_inside(wtg_grid const* g, double const* x)
{
	qaws_scalar p[3] = { (qaws_scalar)x[0], (qaws_scalar)x[1], (qaws_scalar)x[2] }, w = 0;
	if (g->off > 0)
	{
		qaws_mesh_hit h;
		if (qaws_mesh_bvh_closest(&h, g->b, p, 1, (qaws_scalar)g->off, NULL) == QAWS_STATUS_OK && h.element != QAWS_MESH_NONE)
			return 1;
	}
	qaws_mesh_bvh_winding_number(&w, g->b, p, 1, NULL);
	return w >= 0.5;
}

/* the symmetric 3 x 3 A's eigenvalues and eigenvectors (columns of V) by Jacobi rotations */
void wt_sym_eigen3(double* A, double* V)
{
	unsigned int it, p, q, k;
	for (k = 0; k < 9; k++)
		V[k] = k % 4 == 0;
	for (it = 0; it < 32; it++)
	{
		double off = A[1] * A[1] + A[2] * A[2] + A[5] * A[5];
		if (off < 1e-30)
			break;
		for (p = 0; p < 2; p++)
			for (q = p + 1; q < 3; q++)
			{
				double apq = A[3 * p + q], app, aqq, th, t, c, s;
				if (fabs(apq) < 1e-300)
					continue;
				app = A[4 * p];
				aqq = A[4 * q];
				th = (aqq - app) / (2 * apq);
				t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
				c = 1 / sqrt(t * t + 1);
				s = t * c;
				for (k = 0; k < 3; k++)
				{
					double akp = A[3 * k + p], akq = A[3 * k + q];
					A[3 * k + p] = c * akp - s * akq;
					A[3 * k + q] = s * akp + c * akq;
				}
				for (k = 0; k < 3; k++)
				{
					double apk = A[3 * p + k], aqk = A[3 * q + k];
					A[3 * p + k] = c * apk - s * aqk;
					A[3 * q + k] = s * apk + c * aqk;
				}
				for (k = 0; k < 3; k++)
				{
					double vkp = V[3 * k + p], vkq = V[3 * k + q];
					V[3 * k + p] = c * vkp - s * vkq;
					V[3 * k + q] = s * vkp + c * vkq;
				}
			}
	}
}

/* Dual contouring of the crossings on a grid (lo, cell h, n points per axis; in per point: 1 inside): a vertex per crossed cell at
   the least squares point of its crossings' planes (Ju et al.), kept within the cell; a quad per crossing round its edge's four
   cells, turned toward the outside, split along its shorter diagonal; then qaws_mesh_repair (identical points only, non-manifold
   places split, slivers dropped). *out_sharp += the vertices on an edge or a corner. */
qaws_status wt_dual_contour(qaws_mesh_3d** out, unsigned int* out_sharp, int* out_closed, double const* lo, double h, unsigned int const* n,
	unsigned char const* in, wt_cross const* X, size_t nc)
{
	unsigned int *cellv = NULL, *ftri = NULL, nfv = 0, nft = 0, i, k;
	qaws_scalar* fx = NULL;
	size_t q, e, ncell;
	qaws_status st = QAWS_STATUS_OK;
	/* a vertex per crossed cell: the least squares point of its crossings' planes (Ju et al.), kept within the cell */
	ncell = (size_t)(n[0] - 1) * (n[1] - 1) * (n[2] - 1);
	cellv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ncell));
	{
		double* acc = (double*)qaws_internal_scratch_zalloc(qaws_internal_mul2(ncell, 13), sizeof(double));   /* A (6), b (3), mass (3), count */
		if (!cellv || !acc)
		{
			qaws_internal_free(acc);
			goto oom;
		}
		for (q = 0; q < ncell; q++)
			cellv[q] = QAWS_MESH_NONE;
		for (e = 0; e < nc; e++)
		{
			wt_cross const* x = X + e;
			unsigned int ci = (unsigned int)(x->p % n[0]), cj = (unsigned int)(x->p / n[0] % n[1]), ck = (unsigned int)(x->p / ((size_t)n[0] * n[1]));
			int s, t;
			for (s = -1; s <= 0; s++)
				for (t = -1; t <= 0; t++)
				{
					int c3[3];
					size_t cell;
					double const* nr = x->nrm;
					double dn;
					c3[0] = (int)ci;
					c3[1] = (int)cj;
					c3[2] = (int)ck;
					c3[(x->axis + 1) % 3] += s;
					c3[(x->axis + 2) % 3] += t;
					if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= (int)n[0] - 1 || c3[1] >= (int)n[1] - 1 || c3[2] >= (int)n[2] - 1)
						continue;
					cell = (size_t)c3[0] + (size_t)(n[0] - 1) * ((size_t)c3[1] + (size_t)(n[1] - 1) * (size_t)c3[2]);
					dn = nr[0] * x->x[0] + nr[1] * x->x[1] + nr[2] * x->x[2];
					acc[13 * cell] += nr[0] * nr[0];
					acc[13 * cell + 1] += nr[0] * nr[1];
					acc[13 * cell + 2] += nr[0] * nr[2];
					acc[13 * cell + 3] += nr[1] * nr[1];
					acc[13 * cell + 4] += nr[1] * nr[2];
					acc[13 * cell + 5] += nr[2] * nr[2];
					for (i = 0; i < 3; i++)
					{
						acc[13 * cell + 6 + i] += nr[i] * dn;
						acc[13 * cell + 9 + i] += x->x[i];
					}
					acc[13 * cell + 12] += 1;
				}
		}
		for (q = 0, e = 0; q < ncell; q++)
			e += acc[13 * q + 12] > 0;
		fx = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, e + 1));
		if (!fx)
		{
			qaws_internal_free(acc);
			goto oom;
		}
		for (q = 0; q < ncell; q++)
		{
			double const* c = acc + 13 * q;
			double A[9], V[9], m[3], r[3], y[3] = { 0, 0, 0 }, lmax = 0, cl[3];
			unsigned int rank = 0;
			if (c[12] == 0)
				continue;
			for (i = 0; i < 3; i++)
				m[i] = c[9 + i] / c[12];
			A[0] = c[0]; A[1] = c[1]; A[2] = c[2];
			A[3] = c[1]; A[4] = c[3]; A[5] = c[4];
			A[6] = c[2]; A[7] = c[4]; A[8] = c[5];
			/* the right-hand side about the mass point: b - A m */
			for (i = 0; i < 3; i++)
				r[i] = c[6 + i] - (A[3 * i] * m[0] + A[3 * i + 1] * m[1] + A[3 * i + 2] * m[2]);
			wt_sym_eigen3(A, V);
			for (i = 0; i < 3; i++)
				lmax = fmax(lmax, fabs(A[4 * i]));
			for (k = 0; k < 3; k++)
			{
				double l = A[4 * k], proj;
				if (!(fabs(l) > 0.1 * lmax))
					continue;
				rank++;
				proj = (V[k] * r[0] + V[3 + k] * r[1] + V[6 + k] * r[2]) / l;
				for (i = 0; i < 3; i++)
					y[i] += proj * V[3 * i + k];
			}
			*out_sharp += rank >= 2;
			{
				unsigned int ci = (unsigned int)(q % (n[0] - 1)), cj = (unsigned int)(q / (n[0] - 1) % (n[1] - 1)),
					ck = (unsigned int)(q / ((size_t)(n[0] - 1) * (n[1] - 1)));
				cl[0] = lo[0] + h * ci;
				cl[1] = lo[1] + h * cj;
				cl[2] = lo[2] + h * ck;
			}
			for (i = 0; i < 3; i++)
			{
				y[i] += m[i];
				if (y[i] < cl[i] - 0.5 * h || y[i] > cl[i] + 1.5 * h)
				{
					memcpy(y, m, sizeof y);   /* past the cell: the mass point */
					break;
				}
			}
			cellv[q] = nfv;
			fx[3 * (size_t)nfv] = (qaws_scalar)y[0];
			fx[3 * (size_t)nfv + 1] = (qaws_scalar)y[1];
			fx[3 * (size_t)nfv + 2] = (qaws_scalar)y[2];
			nfv++;
		}
		qaws_internal_free(acc);
	}
	/* a quad per crossing round its edge's four cells, turned toward the outside, split along its shorter diagonal */
	ftri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, nc + 1));
	if (!ftri)
		goto oom;
	for (e = 0; e < nc; e++)
	{
		wt_cross const* x = X + e;
		unsigned int ci = (unsigned int)(x->p % n[0]), cj = (unsigned int)(x->p / n[0] % n[1]), ck = (unsigned int)(x->p / ((size_t)n[0] * n[1])), v[4];
		static int const ds[4][2] = { { -1, -1 }, { 0, -1 }, { 0, 0 }, { -1, 0 } };
		int ok = 1;
		for (k = 0; k < 4; k++)
		{
			int c3[3];
			c3[0] = (int)ci;
			c3[1] = (int)cj;
			c3[2] = (int)ck;
			c3[(x->axis + 1) % 3] += ds[k][0];
			c3[(x->axis + 2) % 3] += ds[k][1];
			if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= (int)n[0] - 1 || c3[1] >= (int)n[1] - 1 || c3[2] >= (int)n[2] - 1)
			{
				ok = 0;
				break;
			}
			v[k] = cellv[(size_t)c3[0] + (size_t)(n[0] - 1) * ((size_t)c3[1] + (size_t)(n[1] - 1) * (size_t)c3[2])];
		}
		if (!ok)
			continue;
		if (in[x->p] != 1)
		{
			unsigned int t = v[1];
			v[1] = v[3];
			v[3] = t;
		}
		{
			double d02 = 0, d13 = 0;
			for (i = 0; i < 3; i++)
			{
				double u = (double)fx[3 * (size_t)v[0] + i] - (double)fx[3 * (size_t)v[2] + i], w = (double)fx[3 * (size_t)v[1] + i] - (double)fx[3 * (size_t)v[3] + i];
				d02 += u * u;
				d13 += w * w;
			}
			if (d02 <= d13)
			{
				unsigned int t6[6] = { v[0], v[1], v[2], v[0], v[2], v[3] };
				memcpy(ftri + 3 * (size_t)nft, t6, sizeof t6);
			}
			else
			{
				unsigned int t6[6] = { v[0], v[1], v[3], v[1], v[2], v[3] };
				memcpy(ftri + 3 * (size_t)nft, t6, sizeof t6);
			}
			nft += 2;
		}
	}
	/* through the repair: identical points only, non-manifold places split, slivers dropped, nothing filled */
	{
		qaws_mesh_repair_desc rd;
		qaws_mesh_repair_report rr;
		memset(&rd, 0, sizeof rd);
		rd.weld_tolerance = -1;
		rd.gap_tolerance = -1;
		rd.hole_perimeter = -1;
		rd.flags = QAWS_MESH_REPAIR_SPLIT_NON_MANIFOLD;
		st = qaws_mesh_repair(out, &rr, fx, ftri, NULL, nfv, nft, &rd);
		*out_closed = rr.closed;
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(cellv); qaws_internal_free(fx); qaws_internal_free(ftri);
	return st;
}

qaws_status qaws_mesh_make_watertight(qaws_mesh_3d** out, qaws_watertight_report* report, qaws_scalar const* positions,
	unsigned int const* face_vertices, unsigned int const* face_sizes, unsigned int vertex_count, unsigned int face_count,
	qaws_watertight_desc const* desc)
{
	qaws_watertight_report local;
	wtg_grid g;
	qaws_mesh_bvh_desc bd;
	unsigned int *tri = NULL, nt = 0, f, k, res, i, a;
	size_t off = 0, N = 0, q, nc = 0, cc = 0, e;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, ext = 0;
	unsigned char* in = NULL;     /* per grid point: 0 unknown, 1 inside, 2 outside */
	unsigned char* near = NULL;
	qaws_scalar *gp = NULL, *gv = NULL;
	qaws_mesh_hit* hits = NULL;
	size_t* queue = NULL;
	wt_cross* X = NULL;
	qaws_status st = QAWS_STATUS_OK;
	qaws_scalar const* pos = positions;
	qaws_scalar* ppos = NULL;
	unsigned int pnv = vertex_count;
	memset(&g, 0, sizeof g);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!face_count || !face_vertices || (vertex_count && !positions) || (desc && !(desc->offset >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * vertex_count; i++)
		if (!isfinite((double)positions[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	/* the faces fanned into triangles, the box of the vertices used */
	for (f = 0; f < face_count; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3;
		if (s < 3)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (k = 0; k < s; k++)
			if (face_vertices[off + k] >= vertex_count)
				return QAWS_STATUS_INVALID_ARGUMENT;
		nt += s - 2;
		off += s;
	}
	tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, nt));
	if (!tri)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0, off = 0, nt = 0; f < face_count; f++)
	{
		unsigned int s = face_sizes ? face_sizes[f] : 3;
		for (k = 1; k + 1 < s; k++, nt++)
		{
			tri[3 * (size_t)nt] = face_vertices[off];
			tri[3 * (size_t)nt + 1] = face_vertices[off + k];
			tri[3 * (size_t)nt + 2] = face_vertices[off + k + 1];
		}
		off += s;
	}
	/* holes filled first, tangent to the surface round them (not with an offset: a sheet's border is no hole there) */
	if (!(desc && (desc->flags & QAWS_WATERTIGHT_WINDING_CAPS)) && !(desc && desc->offset > 0))
	{
		qaws_mesh_repair_report rr;
		qaws_mesh_3d* pre = NULL;
		st = qaws_mesh_repair(&pre, &rr, positions, face_vertices, face_sizes, vertex_count, face_count, NULL);
		if (st != QAWS_STATUS_OK)
			goto done;
		report->holes_filled = rr.holes_filled;
		if (!rr.holes_filled)
			qaws_mesh_3d_destroy(pre);   /* (no hole: the input as given, its faces exact) */
		else
		{
		pnv = qaws_mesh_3d_get_vertex_count(pre);
		nt = qaws_mesh_3d_get_triangle_count(pre);
		ppos = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)pnv + 1));
		qaws_internal_free(tri);
		tri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
		if (!ppos || !tri)
		{
			qaws_mesh_3d_destroy(pre);
			goto oom;
		}
		qaws_mesh_3d_get_positions(ppos, pre);
		qaws_mesh_3d_get_triangles(tri, pre);
		qaws_mesh_3d_destroy(pre);
		pos = ppos;
		}
	}
	for (f = 0; f < 3 * nt; f++)
		for (i = 0; i < 3; i++)
		{
			double v = (double)pos[3 * (size_t)tri[f] + i];
			lo[i] = fmin(lo[i], v);
			hi[i] = fmax(hi[i], v);
		}
	memset(&bd, 0, sizeof bd);
	bd.positions = pos;
	bd.vertex_count = pnv;
	bd.dimension = 3;
	bd.simplices = tri;
	bd.simplex_size = 3;
	bd.simplex_count = nt;
	st = qaws_mesh_bvh_create(&g.b, &bd);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* the grid: the box, padded by the offset and two cells */
	res = desc && desc->resolution ? desc->resolution : 64;
	g.off = desc ? (double)desc->offset : 0;
	for (i = 0; i < 3; i++)
		ext = fmax(ext, hi[i] - lo[i] + 2 * g.off);
	g.h = ext > 0 ? ext / res : 1;
	for (i = 0; i < 3; i++)
	{
		g.lo[i] = lo[i] - g.off - 2 * g.h;
		g.n[i] = (unsigned int)ceil((hi[i] - lo[i] + 2 * g.off + 4 * g.h) / g.h) + 1;
		report->grid[i] = g.n[i];
	}
	report->cell = (qaws_scalar)g.h;
	N = (size_t)g.n[0] * g.n[1] * g.n[2];
	in = (unsigned char*)qaws_internal_scratch_zalloc(N, 1);
	near = (unsigned char*)qaws_internal_scratch_zalloc(N, 1);
	gp = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, N));
	gv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), N));
	hits = (qaws_mesh_hit*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_mesh_hit), N));
	queue = (size_t*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(size_t), N));
	if (!in || !near || !gp || !gv || !hits || !queue)
		goto oom;
	for (q = 0; q < N; q++)
	{
		double p[3];
		wtg_point(&g, q, p);
		gp[3 * q] = (qaws_scalar)p[0];
		gp[3 * q + 1] = (qaws_scalar)p[1];
		gp[3 * q + 2] = (qaws_scalar)p[2];
	}
	/* the points within a cell diagonal (and the offset) of the input: their side computed each */
	st = qaws_mesh_bvh_closest(hits, g.b, gp, (unsigned int)N, (qaws_scalar)(g.h * 1.7321 + g.off), NULL);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (q = 0, nc = 0; q < N; q++)
		if (hits[q].element != QAWS_MESH_NONE)
		{
			near[q] = 1;
			memcpy(gp + 3 * nc, gp + 3 * q, sizeof(qaws_scalar) * 3);
			queue[nc++] = q;
		}
	st = qaws_mesh_bvh_winding_number(gv, g.b, gp, (unsigned int)nc, NULL);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (q = 0; q < nc; q++)
		in[queue[q]] = gv[q] >= 0.5 || (g.off > 0 && hits[queue[q]].distance <= g.off) ? 1 : 2;
	report->evaluated = (unsigned int)nc;
	/* the far points: winding numbers on every other point (and the last of each row); a block of 2 x 2 x 2 cells with no near
	   point whose corners agree by more than 0.1 takes their side, any other has its points evaluated (the inside can turn far
	   from the input, across a hole) */
	{
		double* wv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), N));
		size_t m = 0, pass;
		unsigned int bi, bj, bk, ii, jj, kk;
		if (!wv)
			goto oom;
		for (pass = 0; pass < 2; pass++)
		{
			m = 0;
			if (!pass)
			{
				for (q = 0; q < N; q++)
				{
					unsigned int ci = (unsigned int)(q % g.n[0]), cj = (unsigned int)(q / g.n[0] % g.n[1]), ck = (unsigned int)(q / ((size_t)g.n[0] * g.n[1]));
					if (!in[q] && (!(ci & 1) || ci + 1 == g.n[0]) && (!(cj & 1) || cj + 1 == g.n[1]) && (!(ck & 1) || ck + 1 == g.n[2]))
						queue[m++] = q;
				}
			}
			else
				for (bk = 0; bk + 1 < g.n[2]; bk += 2)
					for (bj = 0; bj + 1 < g.n[1]; bj += 2)
						for (bi = 0; bi + 1 < g.n[0]; bi += 2)
						{
							unsigned int ei = bi + 2 < g.n[0] ? bi + 2 : g.n[0] - 1, ej = bj + 2 < g.n[1] ? bj + 2 : g.n[1] - 1, ek = bk + 2 < g.n[2] ? bk + 2 : g.n[2] - 1;
							int amb = 0, lo4 = 1, hi6 = 1;
							unsigned char side;
							for (kk = bk; kk <= ek && !amb; kk++)
								for (jj = bj; jj <= ej && !amb; jj++)
									for (ii = bi; ii <= ei && !amb; ii++)
										amb = near[wtg_id(&g, ii, jj, kk)];
							for (kk = 0; kk < 8 && !amb; kk++)
							{
								double v = wv[wtg_id(&g, kk & 1 ? ei : bi, kk & 2 ? ej : bj, kk & 4 ? ek : bk)];
								lo4 &= v < 0.4;
								hi6 &= v > 0.6;
							}
							amb |= !lo4 && !hi6;
							side = hi6 ? 1 : 2;
							for (kk = bk; kk <= ek; kk++)
								for (jj = bj; jj <= ej; jj++)
									for (ii = bi; ii <= ei; ii++)
									{
										size_t id = wtg_id(&g, ii, jj, kk);
										if (in[id])
											continue;
										if (amb)
										{
											in[id] = 3;   /* (pending) */
											queue[m++] = id;
										}
										else
											in[id] = side;
									}
						}
			for (q = 0; q < m; q++)
			{
				double p[3];
				wtg_point(&g, queue[q], p);
				gp[3 * q] = (qaws_scalar)p[0];
				gp[3 * q + 1] = (qaws_scalar)p[1];
				gp[3 * q + 2] = (qaws_scalar)p[2];
			}
			st = qaws_mesh_bvh_winding_number(gv, g.b, gp, (unsigned int)m, NULL);
			if (st != QAWS_STATUS_OK)
			{
				qaws_internal_free(wv);
				goto done;
			}
			for (q = 0; q < m; q++)
			{
				wv[queue[q]] = gv[q];
				in[queue[q]] = gv[q] >= 0.5 ? 1 : 2;
			}
			report->evaluated += (unsigned int)m;
		}
		qaws_internal_free(wv);
	}
	/* the crossings: grid edges whose ends lie on two sides */
	nc = 0;
	for (q = 0; q < N; q++)
	{
		unsigned int ci = (unsigned int)(q % g.n[0]), cj = (unsigned int)(q / g.n[0] % g.n[1]), ck = (unsigned int)(q / ((size_t)g.n[0] * g.n[1]));
		size_t step[3];
		step[0] = 1;
		step[1] = g.n[0];
		step[2] = (size_t)g.n[0] * g.n[1];
		for (a = 0; a < 3; a++)
		{
			unsigned int c3[3];
			c3[0] = ci;
			c3[1] = cj;
			c3[2] = ck;
			if (c3[a] + 1 >= g.n[a] || in[q] == in[q + step[a]])
				continue;
			if (nc == cc)
			{
				size_t cap = cc ? 2 * cc : 1024;
				wt_cross* t = (wt_cross*)qaws_internal_realloc(X, qaws_internal_mul2(sizeof(wt_cross), cap));
				if (!t)
					goto oom;
				X = t;
				cc = cap;
			}
			X[nc].p = q;
			X[nc++].axis = a;
		}
	}
	report->crossings = (unsigned int)nc;
	/* their Hermite data: a ray along the edge from the inside end onto the input, else bisection; normals from inside out */
	for (e = 0; e < nc; e++)
	{
		wt_cross* x = X + e;
		double p0[3], p1[3], pin[3], pout[3], d[3], l;
		size_t q1 = x->p + (x->axis == 0 ? 1 : x->axis == 1 ? g.n[0] : (size_t)g.n[0] * g.n[1]);
		int found = 0;
		wtg_point(&g, x->p, p0);
		wtg_point(&g, q1, p1);
		memcpy(pin, in[x->p] == 1 ? p0 : p1, sizeof pin);
		memcpy(pout, in[x->p] == 1 ? p1 : p0, sizeof pout);
		for (i = 0; i < 3; i++)
			d[i] = pout[i] - pin[i];
		if (g.off == 0)
		{
			qaws_scalar o[3] = { (qaws_scalar)pin[0], (qaws_scalar)pin[1], (qaws_scalar)pin[2] }, dd[3] = { (qaws_scalar)d[0], (qaws_scalar)d[1], (qaws_scalar)d[2] };
			qaws_mesh_hit h;
			if (qaws_mesh_bvh_ray(&h, g.b, o, dd, 1, 1, 0, NULL) == QAWS_STATUS_OK && h.element != QAWS_MESH_NONE)
			{
				double const* A;
				double e1[3], e2[3], P[3][3];
				for (k = 0; k < 3; k++)
					for (i = 0; i < 3; i++)
						P[k][i] = (double)pos[3 * (size_t)h.vertices[k] + i];
				A = P[0];
				for (i = 0; i < 3; i++)
				{
					e1[i] = P[1][i] - A[i];
					e2[i] = P[2][i] - A[i];
					x->x[i] = (double)h.point[i];
				}
				x->nrm[0] = e1[1] * e2[2] - e1[2] * e2[1];
				x->nrm[1] = e1[2] * e2[0] - e1[0] * e2[2];
				x->nrm[2] = e1[0] * e2[1] - e1[1] * e2[0];
				found = 1;
			}
		}
		if (!found)
		{
			/* bisection of the side; the normal toward the closest point's other side */
			double a0[3], a1[3], mid[3];
			unsigned int it;
			qaws_scalar mp[3];
			qaws_mesh_hit h;
			memcpy(a0, pin, sizeof a0);
			memcpy(a1, pout, sizeof a1);
			for (it = 0; it < 12; it++)
			{
				for (i = 0; i < 3; i++)
					mid[i] = 0.5 * (a0[i] + a1[i]);
				if (wtg_inside(&g, mid))
					memcpy(a0, mid, sizeof a0);
				else
					memcpy(a1, mid, sizeof a1);
			}
			for (i = 0; i < 3; i++)
			{
				x->x[i] = 0.5 * (a0[i] + a1[i]);
				mp[i] = (qaws_scalar)x->x[i];
				x->nrm[i] = d[i];
			}
			if (g.off > 0 && qaws_mesh_bvh_closest(&h, g.b, mp, 1, 0, NULL) == QAWS_STATUS_OK && h.element != QAWS_MESH_NONE
				&& h.distance > 0.5 * g.off)
				for (i = 0; i < 3; i++)
					x->nrm[i] = x->x[i] - (double)h.point[i];
			else if (qaws_mesh_bvh_closest(&h, g.b, mp, 1, (qaws_scalar)(1e-3 * g.h), NULL) == QAWS_STATUS_OK && h.element != QAWS_MESH_NONE)
			{
				/* on the input (a grid point on a face's plane, the ray missed it): that triangle's normal */
				double e1[3], e2[3];
				for (i = 0; i < 3; i++)
				{
					e1[i] = (double)pos[3 * (size_t)h.vertices[1] + i] - (double)pos[3 * (size_t)h.vertices[0] + i];
					e2[i] = (double)pos[3 * (size_t)h.vertices[2] + i] - (double)pos[3 * (size_t)h.vertices[0] + i];
				}
				x->nrm[0] = e1[1] * e2[2] - e1[2] * e2[1];
				x->nrm[1] = e1[2] * e2[0] - e1[0] * e2[2];
				x->nrm[2] = e1[0] * e2[1] - e1[1] * e2[0];
			}
			else
			{
				/* across a hole: down the winding number's gradient (central differences at a quarter cell) */
				qaws_scalar sp[18], sw[6];
				for (k = 0; k < 6; k++)
					for (i = 0; i < 3; i++)
						sp[3 * k + i] = (qaws_scalar)(x->x[i] + (i == k / 2 ? (k & 1 ? -0.25 : 0.25) * g.h : 0));
				if (qaws_mesh_bvh_winding_number(sw, g.b, sp, 6, NULL) == QAWS_STATUS_OK)
					for (i = 0; i < 3; i++)
						x->nrm[i] = (double)sw[2 * i + 1] - (double)sw[2 * i];
			}
		}
		l = sqrt(x->nrm[0] * x->nrm[0] + x->nrm[1] * x->nrm[1] + x->nrm[2] * x->nrm[2]);
		if (!(l > 0))
		{
			memcpy(x->nrm, d, sizeof d);
			l = g.h;
		}
		for (i = 0; i < 3; i++)
			x->nrm[i] /= l;
		if (x->nrm[0] * d[0] + x->nrm[1] * d[1] + x->nrm[2] * d[2] < 0)
			for (i = 0; i < 3; i++)
				x->nrm[i] = -x->nrm[i];
	}
	st = wt_dual_contour(out, &report->sharp, &report->closed, g.lo, g.h, g.n, in, X, nc);
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_mesh_bvh_destroy(g.b);
	qaws_internal_free(tri); qaws_internal_free(in); qaws_internal_free(near); qaws_internal_free(gp); qaws_internal_free(gv); qaws_internal_free(hits);
	qaws_internal_free(queue); qaws_internal_free(X); qaws_internal_free(ppos);
	return st;
}

/* ---- point clouds ------------------------------------------------------- */

#define WTP_K 10u          /* neighbours a point's area is measured to */
#define WTP_CHUNK 4096u    /* winding numbers per batch */

typedef struct wtp_field
{
	qaws_mesh_bvh* b;
	double const* nrm;     /* per disk (point) its normal */
} wtp_field;

/* phi = 1/2 - the winding number of the points' disks (negative inside); its gradient's direction the nearest point's
   normal (the winding number's own gradient is all near field there: the 1/2 level hugs the disks; the nearest, not a
   mean, keeps a sampled edge sharp) */
static qaws_status wtp_eval(void* ctx, double const* x, size_t n, double* phi, double* grad)
{
	wtp_field const* f = (wtp_field const*)ctx;
	unsigned int per = 1u, k;
	qaws_scalar *q = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 3 * (size_t)WTP_CHUNK * per), *w = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (size_t)WTP_CHUNK * per);
	size_t c0, i;
	qaws_status st = QAWS_STATUS_OK;
	if (!q || !w)
	{
		qaws_internal_free(q);
		qaws_internal_free(w);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (c0 = 0; c0 < n && st == QAWS_STATUS_OK; c0 += WTP_CHUNK)
	{
		unsigned int m = (unsigned int)(n - c0 < WTP_CHUNK ? n - c0 : WTP_CHUNK), s;
		for (i = 0; i < m; i++)
			for (s = 0; s < per; s++)
				for (k = 0; k < 3; k++)
					q[3 * ((size_t)i * per + s) + k] = (qaws_scalar)x[3 * (c0 + i) + k];
		st = qaws_mesh_bvh_winding_number(w, f->b, q, m * per, NULL);
		for (i = 0; st == QAWS_STATUS_OK && i < m; i++)
			phi[c0 + i] = 0.5 - (double)w[(size_t)i * per];
		if (st == QAWS_STATUS_OK && grad)
		{
			qaws_mesh_hit* hit = (qaws_mesh_hit*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_hit) * m);
			if (!hit)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				break;
			}
			st = qaws_mesh_bvh_closest(hit, f->b, q, m, 0, NULL);
			for (i = 0; st == QAWS_STATUS_OK && i < m; i++)
				for (k = 0; k < 3; k++)
					grad[3 * (c0 + i) + k] = hit[i].element != QAWS_MESH_NONE ? f->nrm[3 * (size_t)hit[i].element + k] : 0;   /* (any length: its direction is read) */
			qaws_internal_free(hit);
		}
	}
	qaws_internal_free(q);
	qaws_internal_free(w);
	return st;
}

static int wtp_dcmp(void const* p, void const* q)
{
	double a = *(double const*)p, b = *(double const*)q;
	return a < b ? -1 : a > b ? 1 : 0;
}

qaws_status qaws_mesh_make_watertight_points(qaws_mesh_3d** out, qaws_watertight_report* report, qaws_scalar const* positions,
	qaws_scalar const* normals, unsigned int count, qaws_watertight_desc const* desc)
{
	qaws_watertight_report local;
	wtp_field f;
	wt_contour_info ci;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, ext = 0, cs, *area = NULL, *d2 = NULL, spacing = 0;
	unsigned int res = desc && desc->resolution ? desc->resolution : 64, i, k, g[3], *cell0 = NULL, *cellp = NULL, *tri = NULL, nt = 0;
	qaws_scalar* tx = NULL;
	double* dn = NULL;
	qaws_mesh_bvh_desc bd;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!positions || !normals || count < 4 || (desc && desc->offset != 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * count; i++)
		if (!isfinite((double)positions[i]) || !isfinite((double)normals[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&f, 0, sizeof f);
	for (i = 0; i < count; i++)
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(lo[k], (double)positions[3 * (size_t)i + k]);
			hi[k] = fmax(hi[k], (double)positions[3 * (size_t)i + k]);
		}
	for (k = 0; k < 3; k++) ext = fmax(ext, hi[k] - lo[k]);
	if (!(ext > 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* each point's area: pi r^2 / K, r its K-th neighbour's distance, found on a grid of about one point per cell */
	cs = ext / fmax(1, cbrt((double)count));
	for (k = 0; k < 3; k++) g[k] = (unsigned int)fmin(1024, floor((hi[k] - lo[k]) / cs) + 1);
	area = (double*)qaws_internal_scratch_malloc(sizeof(double) * ((size_t)count + 1));
	d2 = (double*)qaws_internal_scratch_malloc(sizeof(double) * ((size_t)count + 1));
	cell0 = (unsigned int*)qaws_internal_scratch_zalloc((size_t)g[0] * g[1] * g[2] + 2, sizeof(unsigned int));
	cellp = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)count + 1));
	tri = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)count + 1));
	tx = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (9 * (size_t)count + 1));
	if (!area || !d2 || !cell0 || !cellp || !tri || !tx)
		goto done;
	{
		size_t nc = (size_t)g[0] * g[1] * g[2], cc;
		unsigned int* fill;
		for (i = 0; i < count; i++)
		{
			unsigned int c3[3];
			for (k = 0; k < 3; k++) c3[k] = (unsigned int)fmin(g[k] - 1, floor(((double)positions[3 * (size_t)i + k] - lo[k]) / cs));
			cell0[c3[0] + (size_t)g[0] * (c3[1] + (size_t)g[1] * c3[2]) + 1]++;
		}
		for (cc = 0; cc < nc; cc++) cell0[cc + 1] += cell0[cc];
		fill = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (nc + 1));
		if (!fill)
			goto done;
		memcpy(fill, cell0, sizeof(unsigned int) * nc);
		for (i = 0; i < count; i++)
		{
			unsigned int c3[3];
			for (k = 0; k < 3; k++) c3[k] = (unsigned int)fmin(g[k] - 1, floor(((double)positions[3 * (size_t)i + k] - lo[k]) / cs));
			cellp[fill[c3[0] + (size_t)g[0] * (c3[1] + (size_t)g[1] * c3[2])]++] = i;
		}
		qaws_internal_free(fill);
	}
	for (i = 0; i < count; i++)
	{
		/* rings of cells round the point until K neighbours are found and the ring's distance passes the K-th */
		double best[WTP_K];
		unsigned int nb = 0, ring, c3[3];
		double pp[3];
		for (k = 0; k < 3; k++)
		{
			pp[k] = (double)positions[3 * (size_t)i + k];
			c3[k] = (unsigned int)fmin(g[k] - 1, floor((pp[k] - lo[k]) / cs));
		}
		for (ring = 0; ring < 1024; ring++)
		{
			long a0, b0, e0;
			for (a0 = (long)c3[0] - (long)ring; a0 <= (long)c3[0] + (long)ring; a0++)
				for (b0 = (long)c3[1] - (long)ring; b0 <= (long)c3[1] + (long)ring; b0++)
					for (e0 = (long)c3[2] - (long)ring; e0 <= (long)c3[2] + (long)ring; e0++)
					{
						size_t id;
						unsigned int j;
						if (a0 < 0 || b0 < 0 || e0 < 0 || a0 >= (long)g[0] || b0 >= (long)g[1] || e0 >= (long)g[2])
							continue;
						if (labs(a0 - (long)c3[0]) != (long)ring && labs(b0 - (long)c3[1]) != (long)ring && labs(e0 - (long)c3[2]) != (long)ring)
							continue;   /* (the ring's shell only) */
						id = (size_t)a0 + (size_t)g[0] * ((size_t)b0 + (size_t)g[1] * (size_t)e0);
						for (j = cell0[id]; j < cell0[id + 1]; j++)
						{
							unsigned int o = cellp[j], s;
							double dd = 0;
							if (o == i)
								continue;
							for (k = 0; k < 3; k++) dd += ((double)positions[3 * (size_t)o + k] - pp[k]) * ((double)positions[3 * (size_t)o + k] - pp[k]);
							if (nb < WTP_K)
								best[nb++] = dd;
							else
							{
								double worst = best[0];
								unsigned int wi = 0;
								for (s = 1; s < WTP_K; s++)
									if (best[s] > worst) { worst = best[s]; wi = s; }
								if (dd < worst)
									best[wi] = dd;
							}
						}
					}
			if (nb == WTP_K)
			{
				double worst = 0;
				unsigned int s;
				for (s = 0; s < WTP_K; s++) worst = fmax(worst, best[s]);
				if (worst <= (ring * cs) * (ring * cs))
					break;   /* (nothing past this ring can be nearer) */
			}
		}
		{
			double worst = 0;
			unsigned int s;
			for (s = 0; s < nb; s++) worst = fmax(worst, best[s]);
			d2[i] = worst;
			area[i] = nb ? 3.14159265358979323846 * worst / nb : 0;
		}
	}
	/* the spacing: the median K-th distance */
	{
		double* tmp = (double*)qaws_internal_scratch_malloc(sizeof(double) * ((size_t)count + 1));
		if (!tmp)
			goto done;
		memcpy(tmp, d2, sizeof(double) * count);
		qsort(tmp, count, sizeof(double), wtp_dcmp);
		spacing = sqrt(tmp[count / 2] / WTP_K) * 1.7724538509;   /* (sqrt pi r^2 / K: a point's share of the surface, its side) */
		qaws_internal_free(tmp);
	}
	{
		/* outward normals: sum a n . (p - centre) is three times the volume they bound (the divergence theorem), positive */
		double cen[3] = { 0, 0, 0 }, s = 0;
		for (i = 0; i < count; i++)
			for (k = 0; k < 3; k++) cen[k] += (double)positions[3 * (size_t)i + k] / count;
		for (i = 0; i < count; i++)
		{
			double l = 0, dp = 0;
			for (k = 0; k < 3; k++)
			{
				l += (double)normals[3 * (size_t)i + k] * (double)normals[3 * (size_t)i + k];
				dp += (double)normals[3 * (size_t)i + k] * ((double)positions[3 * (size_t)i + k] - cen[k]);
			}
			s += l > 0 ? area[i] * dp / sqrt(l) : 0;
		}
		if (!(s > 0))
		{
			st = QAWS_STATUS_INVALID_TOPOLOGY;
			goto done;
		}
	}
	/* each point a triangle of its area in its tangent plane, turned by its normal: its dipole as the winding number sees it */
	for (i = 0; i < count; i++)
	{
		double n[3], l = 0, t1[3], t2[3], R;
		unsigned int c;
		for (k = 0; k < 3; k++) { n[k] = (double)normals[3 * (size_t)i + k]; l += n[k] * n[k]; }
		l = sqrt(l);
		if (!(l > 0) || !(area[i] > 0))
			continue;
		for (k = 0; k < 3; k++) n[k] /= l;
		if (fabs(n[0]) < 0.9) { t1[0] = 0; t1[1] = -n[2]; t1[2] = n[1]; }
		else { t1[0] = n[2]; t1[1] = 0; t1[2] = -n[0]; }
		l = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
		for (k = 0; k < 3; k++) t1[k] /= l;
		t2[0] = n[1] * t1[2] - n[2] * t1[1];
		t2[1] = n[2] * t1[0] - n[0] * t1[2];
		t2[2] = n[0] * t1[1] - n[1] * t1[0];
		R = sqrt(4 * area[i] / (3 * 1.7320508075688772));
		for (c = 0; c < 3; c++)
		{
			double an = 2.0943951023931953 * c;
			for (k = 0; k < 3; k++)
				tx[9 * (size_t)nt + 3 * c + k] = (qaws_scalar)((double)positions[3 * (size_t)i + k] + R * (cos(an) * t1[k] + sin(an) * t2[k]));
			tri[3 * nt + c] = 3 * nt + c;
		}
		nt++;
	}
	if (!nt)
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto done;
	}
	memset(&bd, 0, sizeof bd);
	bd.positions = tx;
	bd.vertex_count = 3 * nt;
	bd.dimension = 3;
	bd.simplices = tri;
	bd.simplex_size = 3;
	bd.simplex_count = nt;
	if ((st = qaws_mesh_bvh_create(&f.b, &bd)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	dn = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nt + 1));
	if (!dn)
		goto done;
	for (i = 0; i < nt; i++)
	{
		/* each disk's normal weighted by its area: the triangle's cross product halved */
		double const u[3] = { (double)tx[9 * (size_t)i + 3] - (double)tx[9 * (size_t)i], (double)tx[9 * (size_t)i + 4] - (double)tx[9 * (size_t)i + 1],
			(double)tx[9 * (size_t)i + 5] - (double)tx[9 * (size_t)i + 2] };
		double const w2[3] = { (double)tx[9 * (size_t)i + 6] - (double)tx[9 * (size_t)i], (double)tx[9 * (size_t)i + 7] - (double)tx[9 * (size_t)i + 1],
			(double)tx[9 * (size_t)i + 8] - (double)tx[9 * (size_t)i + 2] };
		dn[3 * (size_t)i] = 0.5 * (u[1] * w2[2] - u[2] * w2[1]);
		dn[3 * (size_t)i + 1] = 0.5 * (u[2] * w2[0] - u[0] * w2[2]);
		dn[3 * (size_t)i + 2] = 0.5 * (u[0] * w2[1] - u[1] * w2[0]);
	}
	f.nrm = dn;
	memset(&ci, 0, sizeof ci);
	st = wt_contour_field(out, &ci, wtp_eval, &f, lo, hi, 2 * spacing, res);
	memcpy(report->grid, ci.grid, sizeof ci.grid);
	report->cell = (qaws_scalar)ci.cell;
	report->crossings = ci.crossings;
	report->sharp = ci.sharp;
	report->closed = ci.closed;
	report->evaluated = report->grid[0] * report->grid[1] * report->grid[2];
	if (st == QAWS_STATUS_OK && !*out)
		st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (nothing inside: normals turned inward, or no surface the points sample) */
done:
	qaws_mesh_bvh_destroy(f.b);
	qaws_internal_free(area); qaws_internal_free(d2); qaws_internal_free(cell0); qaws_internal_free(cellp); qaws_internal_free(tri); qaws_internal_free(tx);
	qaws_internal_free(dn);
	return st;
}
