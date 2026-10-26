/*
 * Offsets of meshes (qaws_mesh_offset.h): the distance's level set dual
 * contoured (the watertight surfaces' contouring, wt_dual_contour), or the
 * vertices moved along mitred directions and the result's self-intersections
 * resolved by the self-union.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_offset.h"
#include "../qaws_mesh_bvh.h"
#include "../qaws_mesh_boolean.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OFS_CHUNK 4096u   /* points per distance query (bounds the hits' memory) */
#define OFS_FALSI 5u      /* regula falsi rounds per crossing */
#define OFS_MITRE 4.0     /* the largest mitre direction's length (needle tips) */
#define OFS_NONE 0xFFFFFFFFu

typedef struct ofs_field
{
	qaws_mesh_bvh* b;
	double const* P;               /* the input's positions */
	int sgn;                       /* signed distance (SURFACE, SHELL) or unsigned (THICKEN) */
	double level;                  /* phi = distance - level */
	qaws_batch_executor const* ex; /* the queries and the Booleans through it; NULL: this thread */
} ofs_field;

/* phi at n points (3 each) and, if grad, its unit gradient (3 each), in chunks */
static qaws_status ofs_eval(ofs_field const* f, double const* x, size_t n, double* phi, double* grad)
{
	qaws_mesh_hit* hits = (qaws_mesh_hit*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_hit) * OFS_CHUNK);
	qaws_scalar* xs = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * OFS_CHUNK);
	qaws_scalar* dist = xs + 3 * OFS_CHUNK;
	size_t c0, i;
	unsigned int k;
	qaws_status st = QAWS_STATUS_OK;
	if (!hits || !xs)
	{
		qaws_internal_free(hits);
		qaws_internal_free(xs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (c0 = 0; c0 < n && st == QAWS_STATUS_OK; c0 += OFS_CHUNK)
	{
		unsigned int m = (unsigned int)(n - c0 < OFS_CHUNK ? n - c0 : OFS_CHUNK);
		for (i = 0; i < 3 * (size_t)m; i++) xs[i] = (qaws_scalar)x[3 * c0 + i];
		st = f->sgn ? qaws_mesh_bvh_signed_distance(dist, hits, f->b, xs, m, f->ex) : qaws_mesh_bvh_closest(hits, f->b, xs, m, 0, f->ex);
		for (i = 0; st == QAWS_STATUS_OK && i < m; i++)
		{
			double d = f->sgn ? (double)dist[i] : (double)hits[i].distance, g[3], l = 0;
			phi[c0 + i] = d - f->level;
			if (!grad)
				continue;
			for (k = 0; k < 3; k++) { g[k] = x[3 * (c0 + i) + k] - (double)hits[i].point[k]; l += g[k] * g[k]; }
			l = sqrt(l);
			if (l > 1e-12)
			{
				/* the distance's gradient: away from the closest point (into the input where the signed distance is negative) */
				for (k = 0; k < 3; k++) g[k] /= f->sgn && d < 0 ? -l : l;
			}
			else
			{
				/* on the input: its triangle's normal */
				double const *a = f->P + 3 * (size_t)hits[i].vertices[0], *b = f->P + 3 * (size_t)hits[i].vertices[1], *c = f->P + 3 * (size_t)hits[i].vertices[2];
				double u[3], w[3];
				for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; w[k] = c[k] - a[k]; }
				g[0] = u[1] * w[2] - u[2] * w[1]; g[1] = u[2] * w[0] - u[0] * w[2]; g[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
				for (k = 0; k < 3; k++) g[k] = l > 0 ? g[k] / l : 0;
			}
			for (k = 0; k < 3; k++) grad[3 * (c0 + i) + k] = g[k];
		}
	}
	qaws_internal_free(hits);
	qaws_internal_free(xs);
	return st;
}

/* the field's evaluation for wt_contour_field */
static qaws_status ofs_field_eval(void* ctx, double const* x, size_t n, double* phi, double* grad)
{
	return ofs_eval((ofs_field const*)ctx, x, n, phi, grad);
}

/* grid point q's position (lowest corner glo, cell h, n points per axis) */
static void ofs_grid_point(double* x, double const* glo, double h, unsigned int const* n, size_t q)
{
	x[0] = glo[0] + h * (double)(q % n[0]);
	x[1] = glo[1] + h * (double)(q / n[0] % n[1]);
	x[2] = glo[2] + h * (double)(q / ((size_t)n[0] * n[1]));
}

/* The level set phi = 0 of a field (f: phi, negative inside, and its gradient at points) contoured on a grid over the box
   lo..hi grown by margin, res cells along its longest side: grid points evaluated, each grid edge whose ends lie on two
   sides a crossing found by regula falsi on phi, its normal phi's gradient there, then dual contouring (wt_dual_contour).
   *out NULL when it is empty. */
qaws_status wt_contour_field(qaws_mesh_3d** out, wt_contour_info* report, wt_field_fn fe, void* ctx, double const* lo, double const* hi,
	double margin, unsigned int res, double lip)
{
	double glo[3], h, ext = 0, *phi = NULL, *A = NULL, *B = NULL, *pa = NULL, *pb = NULL, *xm = NULL, *fm = NULL, *gm = NULL;
	unsigned int n[3], a, i, it;
	size_t N, q, nc = 0, cc = 0, e;
	unsigned char* in = NULL;
	wt_cross* X = NULL;
	qaws_status st;
	*out = NULL;
	for (i = 0; i < 3; i++) ext = fmax(ext, hi[i] - lo[i] + 2 * margin);
	h = ext / res;
	for (i = 0; i < 3; i++)
	{
		glo[i] = lo[i] - margin - 2 * h;
		n[i] = (unsigned int)ceil((hi[i] - lo[i] + 2 * margin + 4 * h) / h) + 1;
		report->grid[i] = n[i];
	}
	report->cell = h;
	N = (size_t)n[0] * n[1] * n[2];
	/* (per grid point its value and side only: the positions made when needed) */
	phi = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), N));
	in = (unsigned char*)qaws_internal_scratch_malloc(N);
	if (!phi || !in)
		goto oom;
	if (lip > 0)
	{
		/* blocks halved from a coarse stride down to one cell, the field evaluated at their corners level by level: a block
		   with a corner farther from the level than lip times its diagonal holds no crossing (its points take that corner's
		   sign), the others split; every crossing edge ends as two corners of a one-cell block, evaluated */
		unsigned char* done2 = (unsigned char*)qaws_internal_scratch_zalloc(N + 1, 1);   /* 1 evaluated, 2 given a sign */
		st = QAWS_STATUS_OK;
		size_t *blk = NULL, *nxt = NULL, nb = 0, nbn = 0, capb = 0, capn = 0, m, r, *sq = NULL;
		double *sx = NULL, *sv = NULL;
		size_t capq = 0;
		unsigned int s = 1;
		int fail = !done2;
		while (s < n[0] || s < n[1] || s < n[2]) s *= 2;
		s /= 2;
		if (s < 1)
			s = 1;
		report->evaluated = 0;
		/* (the blocks: their lowest corner's grid point; the first level tiling the grid) */
		{
			unsigned int b0, b1, b2;
			for (b2 = 0; !fail && b2 + 1 < n[2]; b2 += s)
				for (b1 = 0; !fail && b1 + 1 < n[1]; b1 += s)
					for (b0 = 0; !fail && b0 + 1 < n[0]; b0 += s)
					{
						if (nb == capb)
						{
							size_t cap = capb ? 2 * capb : 1024;
							size_t* g2 = (size_t*)qaws_internal_realloc(blk, sizeof(size_t) * cap);
							if (!g2) { fail = 1; break; }
							blk = g2;
							capb = cap;
						}
						blk[nb++] = (size_t)b0 + n[0] * ((size_t)b1 + (size_t)n[1] * b2);
					}
		}
		while (!fail && nb)
		{
			/* this level's corners not yet evaluated, in one batch */
			m = 0;
			for (r = 0; r < nb && !fail; r++)
			{
				unsigned int c0 = (unsigned int)(blk[r] % n[0]), c1 = (unsigned int)(blk[r] / n[0] % n[1]), c2 = (unsigned int)(blk[r] / ((size_t)n[0] * n[1])), i0;
				for (i0 = 0; i0 < 8; i0++)
				{
					unsigned int x0 = (i0 & 1) ? (c0 + s < n[0] ? c0 + s : n[0] - 1) : c0, x1 = (i0 & 2) ? (c1 + s < n[1] ? c1 + s : n[1] - 1) : c1;
					unsigned int x2 = (i0 & 4) ? (c2 + s < n[2] ? c2 + s : n[2] - 1) : c2;
					size_t cq = (size_t)x0 + n[0] * ((size_t)x1 + (size_t)n[1] * x2);
					if (done2[cq] == 1)
						continue;
					done2[cq] = 1;
					if (m == capq)
					{
						size_t cap = capq ? 2 * capq : 4096;
						size_t* g2 = (size_t*)qaws_internal_realloc(sq, sizeof(size_t) * cap);
						double* g3 = g2 ? (double*)qaws_internal_realloc(sx, sizeof(double) * 3 * cap) : NULL;
						double* g4 = g3 ? (double*)qaws_internal_realloc(sv, sizeof(double) * cap) : NULL;
						if (g2) sq = g2;
						if (g3) sx = g3;
						if (!g2 || !g3 || !g4) { fail = 1; break; }
						sv = g4;
						capq = cap;
					}
					ofs_grid_point(sx + 3 * m, glo, h, n, cq);
					sq[m++] = cq;
				}
			}
			if (fail)
				break;
			if (m && (st = fe(ctx, sx, m, sv, NULL)) != QAWS_STATUS_OK)
				break;
			for (r = 0; r < m; r++) phi[sq[r]] = sv[r];
			report->evaluated += m;
			/* far blocks fill their points by sign; near ones split (at one cell: done, their corners evaluated) */
			nbn = 0;
			for (r = 0; r < nb && !fail; r++)
			{
				unsigned int c0 = (unsigned int)(blk[r] % n[0]), c1 = (unsigned int)(blk[r] / n[0] % n[1]), c2 = (unsigned int)(blk[r] / ((size_t)n[0] * n[1])), i0, i1, i2;
				unsigned int e0 = c0 + s < n[0] ? c0 + s : n[0] - 1, e1 = c1 + s < n[1] ? c1 + s : n[1] - 1, e2 = c2 + s < n[2] ? c2 + s : n[2] - 1;
				double diag = h * sqrt((double)(e0 - c0) * (e0 - c0) + (double)(e1 - c1) * (e1 - c1) + (double)(e2 - c2) * (e2 - c2)), val = 0;
				int far = 0;
				for (i0 = 0; i0 < 8; i0++)
				{
					size_t cq = (size_t)((i0 & 1) ? e0 : c0) + n[0] * ((size_t)((i0 & 2) ? e1 : c1) + (size_t)n[1] * ((i0 & 4) ? e2 : c2));
					if (fabs(phi[cq]) > lip * diag)
					{
						far = 1;
						val = phi[cq];
					}
				}
				if (far)
				{
					for (i2 = c2; i2 <= e2; i2++)
						for (i1 = c1; i1 <= e1; i1++)
							for (i0 = c0; i0 <= e0; i0++)
							{
								size_t fq = (size_t)i0 + n[0] * ((size_t)i1 + (size_t)n[1] * i2);
								if (!done2[fq])
								{
									phi[fq] = val;
									done2[fq] = 2;
								}
							}
					continue;
				}
				if (s == 1)
					continue;
				for (i0 = 0; i0 < 8; i0++)
				{
					unsigned int h0 = s / 2, x0 = c0 + ((i0 & 1) ? h0 : 0), x1 = c1 + ((i0 & 2) ? h0 : 0), x2 = c2 + ((i0 & 4) ? h0 : 0);
					if (x0 >= e0 || x1 >= e1 || x2 >= e2)
						continue;   /* (past the grid's last point: no child there) */
					if (nbn == capn)
					{
						size_t cap = capn ? 2 * capn : 1024;
						size_t* g2 = (size_t*)qaws_internal_realloc(nxt, sizeof(size_t) * cap);
						if (!g2) { fail = 1; break; }
						nxt = g2;
						capn = cap;
					}
					nxt[nbn++] = (size_t)x0 + n[0] * ((size_t)x1 + (size_t)n[1] * x2);
				}
			}
			{
				size_t* sw = blk;
				size_t cs = capb;
				blk = nxt;
				capb = capn;
				nxt = sw;
				capn = cs;
				nb = nbn;
			}
			s /= 2;
			if (s == 0)
				break;
		}
		qaws_internal_free(done2); qaws_internal_free(blk); qaws_internal_free(nxt); qaws_internal_free(sq); qaws_internal_free(sx); qaws_internal_free(sv);
		if (st != QAWS_STATUS_OK)
			goto done;
		if (fail)
			goto oom;
	}
	else
	{
		/* every point, in chunks */
		double* cx = (double*)qaws_internal_scratch_malloc(sizeof(double) * 3 * 65536);
		size_t q0;
		st = QAWS_STATUS_OK;
		if (!cx)
			goto oom;
		for (q0 = 0; q0 < N; q0 += 65536)
		{
			size_t cn = N - q0 < 65536 ? N - q0 : 65536;
			for (q = 0; q < cn; q++) ofs_grid_point(cx + 3 * q, glo, h, n, q0 + q);
			if ((st = fe(ctx, cx, cn, phi + q0, NULL)) != QAWS_STATUS_OK)
				break;
		}
		qaws_internal_free(cx);
		if (st != QAWS_STATUS_OK)
			goto done;
		report->evaluated = N;
	}
	for (q = 0; q < N; q++) in[q] = phi[q] < 0 ? 1 : 2;
	/* the crossings: grid edges whose ends lie on two sides */
	for (q = 0; q < N; q++)
	{
		unsigned int c3[3];
		size_t step[3];
		c3[0] = (unsigned int)(q % n[0]); c3[1] = (unsigned int)(q / n[0] % n[1]); c3[2] = (unsigned int)(q / ((size_t)n[0] * n[1]));
		step[0] = 1; step[1] = n[0]; step[2] = (size_t)n[0] * n[1];
		for (a = 0; a < 3; a++)
		{
			if (c3[a] + 1 >= n[a] || in[q] == in[q + step[a]])
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
	if (!nc)
	{
		st = QAWS_STATUS_OK;
		goto done;
	}
	/* each crossing by regula falsi on phi between its inside end A and outside end B, its normal phi's gradient there */
	A = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, nc));
	B = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, nc));
	pa = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), nc));
	pb = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), nc));
	xm = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, nc));
	fm = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), nc));
	gm = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, nc));
	if (!A || !B || !pa || !pb || !xm || !fm || !gm)
		goto oom;
	for (e = 0; e < nc; e++)
	{
		size_t q1 = X[e].p + (X[e].axis == 0 ? 1 : X[e].axis == 1 ? n[0] : (size_t)n[0] * n[1]), qi = in[X[e].p] == 1 ? X[e].p : q1, qo = qi == q1 ? X[e].p : q1;
		ofs_grid_point(A + 3 * e, glo, h, n, qi);
		ofs_grid_point(B + 3 * e, glo, h, n, qo);
		pa[e] = phi[qi];
		pb[e] = phi[qo];
	}
	for (it = 0; it <= OFS_FALSI; it++)
	{
		for (e = 0; e < nc; e++)
		{
			double t = pa[e] / (pa[e] - pb[e]);
			if (!(t >= 0 && t <= 1)) t = 0.5;
			for (i = 0; i < 3; i++) xm[3 * e + i] = A[3 * e + i] + t * (B[3 * e + i] - A[3 * e + i]);
		}
		if (it == OFS_FALSI)
			break;
		if ((st = fe(ctx, xm, nc, fm, NULL)) != QAWS_STATUS_OK)
			goto done;
		for (e = 0; e < nc; e++)
		{
			if (fm[e] < 0) { memcpy(A + 3 * e, xm + 3 * e, sizeof(double) * 3); pa[e] = fm[e]; }
			else { memcpy(B + 3 * e, xm + 3 * e, sizeof(double) * 3); pb[e] = fm[e]; }
		}
	}
	if ((st = fe(ctx, xm, nc, fm, gm)) != QAWS_STATUS_OK)
		goto done;
	for (e = 0; e < nc; e++)
	{
		double d[3], s = 0;
		for (i = 0; i < 3; i++)
		{
			X[e].x[i] = xm[3 * e + i];
			X[e].nrm[i] = gm[3 * e + i];
			d[i] = B[3 * e + i] - A[3 * e + i];
			s += X[e].nrm[i] * d[i];
		}
		if (!(s > 0))
		{
			/* (no usable gradient: along the edge, inside to outside) */
			double l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			for (i = 0; i < 3; i++) X[e].nrm[i] = d[i] / l;
		}
	}
	st = wt_dual_contour(out, &report->sharp, &report->closed, glo, h, n, in, X, nc);
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(phi); qaws_internal_free(in); qaws_internal_free(X);
	qaws_internal_free(A); qaws_internal_free(B); qaws_internal_free(pa); qaws_internal_free(pb); qaws_internal_free(xm); qaws_internal_free(fm); qaws_internal_free(gm);
	return st;
}

/* per vertex the mitre direction m: least squares n_f . m = 1 over its triangles, weighted by their angles there (3 each) */
static void ofs_mitre(double* m, double const* P, unsigned int nv, unsigned int const* t, unsigned int nt)
{
	double* acc = (double*)qaws_internal_scratch_zalloc(9 * (size_t)nv + 1, sizeof(double));   /* A (6), b (3) */
	unsigned int f, c, i, k;
	if (!acc)
	{
		memset(m, 0, sizeof(double) * 3 * nv);
		return;
	}
	for (f = 0; f < nt; f++)
	{
		double const* p[3];
		double u[3], w[3], nr[3], l;
		for (c = 0; c < 3; c++) p[c] = P + 3 * (size_t)t[3 * f + c];
		for (k = 0; k < 3; k++) { u[k] = p[1][k] - p[0][k]; w[k] = p[2][k] - p[0][k]; }
		nr[0] = u[1] * w[2] - u[2] * w[1]; nr[1] = u[2] * w[0] - u[0] * w[2]; nr[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
		if (!(l > 0))
			continue;
		for (k = 0; k < 3; k++) nr[k] /= l;
		for (c = 0; c < 3; c++)
		{
			double const *o = p[c], *a = p[(c + 1) % 3], *b = p[(c + 2) % 3];
			double e1[3], e2[3], d = 0, l1 = 0, l2 = 0, ang;
			double* s = acc + 9 * (size_t)t[3 * f + c];
			for (k = 0; k < 3; k++) { e1[k] = a[k] - o[k]; e2[k] = b[k] - o[k]; d += e1[k] * e2[k]; l1 += e1[k] * e1[k]; l2 += e2[k] * e2[k]; }
			ang = l1 > 0 && l2 > 0 ? acos(fmax(-1.0, fmin(1.0, d / sqrt(l1 * l2)))) : 0;
			s[0] += ang * nr[0] * nr[0]; s[1] += ang * nr[0] * nr[1]; s[2] += ang * nr[0] * nr[2];
			s[3] += ang * nr[1] * nr[1]; s[4] += ang * nr[1] * nr[2]; s[5] += ang * nr[2] * nr[2];
			for (k = 0; k < 3; k++) s[6 + k] += ang * nr[k];
		}
	}
	for (i = 0; i < nv; i++)
	{
		double const* s = acc + 9 * (size_t)i;
		double M[9], V[9], lmax = 0, l = 0, y[3] = { 0, 0, 0 };
		M[0] = s[0]; M[1] = s[1]; M[2] = s[2]; M[3] = s[1]; M[4] = s[3]; M[5] = s[4]; M[6] = s[2]; M[7] = s[4]; M[8] = s[5];
		wt_sym_eigen3(M, V);
		for (k = 0; k < 3; k++) lmax = fmax(lmax, fabs(M[4 * k]));
		/* the pseudo-inverse: directions the faces do not span (a flat vertex's tangents) not moved */
		for (k = 0; k < 3; k++)
		{
			double lk = M[4 * k], pr;
			if (!(fabs(lk) > 1e-3 * lmax))
				continue;
			pr = (V[k] * s[6] + V[3 + k] * s[7] + V[6 + k] * s[8]) / lk;
			for (c = 0; c < 3; c++) y[c] += pr * V[3 * c + k];
		}
		l = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
		for (k = 0; k < 3; k++) m[3 * (size_t)i + k] = l > OFS_MITRE ? y[k] * OFS_MITRE / l : y[k];
	}
	qaws_internal_free(acc);
}

static int ofs_edge_cmp(void const* a, void const* b)
{
	unsigned int const *x = (unsigned int const*)a, *y = (unsigned int const*)b;
	return x[0] != y[0] ? (x[0] < y[0] ? -1 : 1) : x[1] != y[1] ? (x[1] < y[1] ? -1 : 1) : 0;
}

/* the border's directed edges (a, b) of the triangles (b, a absent): 2 each into *out, their count; ~0u on allocation failure */
static unsigned int ofs_border(unsigned int** out, unsigned int const* t, unsigned int nt)
{
	unsigned int* e = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, (size_t)nt + 1));
	unsigned int* r = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, (size_t)nt + 1));
	unsigned int f, c, nb = 0;
	*out = r;
	if (!e || !r)
	{
		qaws_internal_free(e);
		return ~0u;
	}
	for (f = 0; f < nt; f++)
		for (c = 0; c < 3; c++)
		{
			e[2 * (3 * f + c)] = t[3 * f + c];
			e[2 * (3 * f + c) + 1] = t[3 * f + (c + 1) % 3];
		}
	qsort(e, 3 * (size_t)nt, sizeof(unsigned int) * 2, ofs_edge_cmp);
	for (f = 0; f < 3 * nt; f++)
	{
		unsigned int key[2];
		key[0] = e[2 * f + 1]; key[1] = e[2 * f];
		if (!bsearch(key, e, 3 * (size_t)nt, sizeof(unsigned int) * 2, ofs_edge_cmp))
		{
			r[2 * nb] = e[2 * f];
			r[2 * nb + 1] = e[2 * f + 1];
			nb++;
		}
	}
	qaws_internal_free(e);
	return nb;
}

/* the offset surface's vertices' distances to the input */
static qaws_status ofs_measure(qaws_mesh_offset_report* report, ofs_field const* f, double const* x, size_t n)
{
	ofs_field u = *f;
	double* d = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), n + 1));
	double lo = HUGE_VAL, hi = 0;
	size_t i;
	qaws_status st;
	if (!d)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	u.sgn = 0;
	u.level = 0;
	st = ofs_eval(&u, x, n, d, NULL);
	for (i = 0; st == QAWS_STATUS_OK && i < n; i++) { lo = fmin(lo, d[i]); hi = fmax(hi, d[i]); }
	report->min_distance = (qaws_scalar)(n ? lo : 0);
	report->max_distance = (qaws_scalar)hi;
	qaws_internal_free(d);
	return st;
}

static qaws_status ofs_measure_mesh(qaws_mesh_offset_report* report, ofs_field const* f, qaws_mesh_3d const* m)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* xs = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)n + 1));
	double* x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)n + 1));
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (xs && x)
	{
		qaws_mesh_3d_get_positions(xs, m);
		for (i = 0; i < 3 * n; i++) x[i] = (double)xs[i];
		st = ofs_measure(report, f, x, n);
	}
	qaws_internal_free(xs);
	qaws_internal_free(x);
	return st;
}

/* ROUND's pieces: closed convex solids, their union the points within r of the input's surface */
typedef struct ofs_soup
{
	double* x;
	unsigned int nx, capx;
	unsigned int* t;
	unsigned int nt, capt;
} ofs_soup;

/* a point of the soup, its index; OFS_NONE on allocation failure */
static unsigned int ofs_soup_point(ofs_soup* s, double const* p)
{
	if (s->nx == s->capx)
	{
		unsigned int cap = s->capx ? 2 * s->capx : 1024;
		double* g = (double*)qaws_internal_realloc(s->x, sizeof(double) * 3 * (size_t)cap);
		if (!g)
			return OFS_NONE;
		s->x = g;
		s->capx = cap;
	}
	memcpy(s->x + 3 * (size_t)s->nx, p, sizeof(double) * 3);
	return s->nx++;
}

static int ofs_soup_tri(ofs_soup* s, unsigned int a, unsigned int b, unsigned int c)
{
	if (a == OFS_NONE || b == OFS_NONE || c == OFS_NONE)
		return 0;
	if (s->nt == s->capt)
	{
		unsigned int cap = s->capt ? 2 * s->capt : 1024;
		unsigned int* g = (unsigned int*)qaws_internal_realloc(s->t, sizeof(unsigned int) * 3 * (size_t)cap);
		if (!g)
			return 0;
		s->t = g;
		s->capt = cap;
	}
	s->t[3 * (size_t)s->nt] = a;
	s->t[3 * (size_t)s->nt + 1] = b;
	s->t[3 * (size_t)s->nt + 2] = c;
	s->nt++;
	return 1;
}

/* union-find root of x (paths halved) */
static unsigned int ofs_root(unsigned int* up, unsigned int x)
{
	while (up[x] != x)
	{
		up[x] = up[up[x]];
		x = up[x];
	}
	return x;
}

/* the cylinder of radius r round segment a b, m sides, its first side's direction along u (a unit vector across it), the
   ends closed */
static int ofs_cylinder(ofs_soup* s, double const* a, double const* b, double const* u, double r, unsigned int m)
{
	double e[3], w[3], l;
	unsigned int base = s->nx, i, k;
	for (k = 0; k < 3; k++) e[k] = b[k] - a[k];
	l = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	if (!(l > 0))
		return 1;
	for (k = 0; k < 3; k++) e[k] /= l;
	w[0] = e[1] * u[2] - e[2] * u[1];
	w[1] = e[2] * u[0] - e[0] * u[2];
	w[2] = e[0] * u[1] - e[1] * u[0];
	for (i = 0; i < 2 * m; i++)
	{
		double q[3], th = 6.283185307179586 * (i % m) / m;
		for (k = 0; k < 3; k++) q[k] = (i < m ? a[k] : b[k]) + r * (cos(th) * u[k] + sin(th) * w[k]);
		if (ofs_soup_point(s, q) == OFS_NONE)
			return 0;
	}
	/* (u, w, e right-handed: the rings run counterclockwise about e, the sides out, the ends fanned) */
	for (i = 0; i < m; i++)
	{
		unsigned int j = (i + 1) % m;
		if (!ofs_soup_tri(s, base + i, base + j, base + m + j) || !ofs_soup_tri(s, base + i, base + m + j, base + m + i))
			return 0;
	}
	for (i = 1; i + 1 < m; i++)
		if (!ofs_soup_tri(s, base, base + i + 1, base + i) || !ofs_soup_tri(s, base + m, base + m + i, base + m + i + 1))
			return 0;
	return 1;
}

/* the ball of radius r about c: m meridians, m / 2 bands */
static int ofs_ball(ofs_soup* s, double const* c, double r, unsigned int m)
{
	unsigned int h = m / 2, base = s->nx, i, j;
	double q[3];
	q[0] = c[0];
	q[1] = c[1];
	q[2] = c[2] - r;
	if (ofs_soup_point(s, q) == OFS_NONE)
		return 0;
	for (j = 1; j < h; j++)
		for (i = 0; i < m; i++)
		{
			double ph = 3.141592653589793 * j / h - 1.5707963267948966, th = 6.283185307179586 * i / m;
			q[0] = c[0] + r * cos(ph) * cos(th);
			q[1] = c[1] + r * cos(ph) * sin(th);
			q[2] = c[2] + r * sin(ph);
			if (ofs_soup_point(s, q) == OFS_NONE)
				return 0;
		}
	q[0] = c[0];
	q[1] = c[1];
	q[2] = c[2] + r;
	if (ofs_soup_point(s, q) == OFS_NONE)
		return 0;
	/* (rings 1 .. h - 1 of m from base + 1, the poles first and last, outward) */
	for (i = 0; i < m; i++)
	{
		unsigned int i1 = (i + 1) % m, top = base + 1 + (h - 1) * m;
		if (!ofs_soup_tri(s, base, base + 1 + i1, base + 1 + i) || !ofs_soup_tri(s, top, base + 1 + (h - 2) * m + i, base + 1 + (h - 2) * m + i1))
			return 0;
		for (j = 1; j + 1 < h; j++)
		{
			unsigned int p0 = base + 1 + (j - 1) * m, p1 = base + 1 + j * m;
			if (!ofs_soup_tri(s, p0 + i, p0 + i1, p1 + i1) || !ofs_soup_tri(s, p0 + i, p1 + i1, p1 + i))
				return 0;
		}
	}
	return 1;
}

typedef struct ofs_hkey
{
	unsigned int lo, hi, f;
} ofs_hkey;

static int ofs_hkey_cmp(void const* a, void const* b)
{
	ofs_hkey const *x = (ofs_hkey const*)a, *y = (ofs_hkey const*)b;
	if (x->lo != y->lo)
		return x->lo < y->lo ? -1 : 1;
	if (x->hi != y->hi)
		return x->hi < y->hi ? -1 : 1;
	return x->f < y->f ? -1 : (x->f > y->f ? 1 : 0);
}

/* a triangle of a piece facing away from inside point o; 0 on allocation failure */
static int ofs_soup_tri_out(ofs_soup* s, unsigned int a, unsigned int b, unsigned int c, double const* o)
{
	double const *p = s->x + 3 * (size_t)a, *q = s->x + 3 * (size_t)b, *w = s->x + 3 * (size_t)c;
	double e1[3], e2[3], n[3], m[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		e1[k] = q[k] - p[k];
		e2[k] = w[k] - p[k];
		m[k] = (p[k] + q[k] + w[k]) / 3 - o[k];
	}
	n[0] = e1[1] * e2[2] - e1[2] * e2[1];
	n[1] = e1[2] * e2[0] - e1[0] * e2[2];
	n[2] = e1[0] * e2[1] - e1[1] * e2[0];
	return n[0] * m[0] + n[1] * m[1] + n[2] * m[2] >= 0 ? ofs_soup_tri(s, a, b, c) : ofs_soup_tri(s, a, c, b);
}

/* the wedge of the cylinder of radius r round segment a b over the directions cos(phi) e0 + sin(phi) e1 (e0, e1 unit,
   across it), phi from -widen to th + widen, steps of step at most: the points within r of the segment that lie that way */
static int ofs_wedge(ofs_soup* s, double const* a, double const* b, double const* e0, double const* e1, double th, double widen, double r, double step)
{
	double span = th + 2 * widen, mid[3], o[3];
	unsigned int k = (unsigned int)ceil(span / step), base, i, j;
	if (k < 1)
		k = 1;
	base = s->nx;
	if (ofs_soup_point(s, a) == OFS_NONE || ofs_soup_point(s, b) == OFS_NONE)
		return 0;
	for (i = 0; i <= k; i++)
	{
		double ph = -widen + span * i / k, pa[3], pb[3];
		for (j = 0; j < 3; j++)
		{
			double u = cos(ph) * e0[j] + sin(ph) * e1[j];
			pa[j] = a[j] + r * u;
			pb[j] = b[j] + r * u;
		}
		if (ofs_soup_point(s, pa) == OFS_NONE || ofs_soup_point(s, pb) == OFS_NONE)
			return 0;
	}
	/* (an inside point: on the middle of the wedge's middle direction, half way out; the wedge spans less than a full turn) */
	for (j = 0; j < 3; j++)
	{
		double u = cos(0.5 * th) * e0[j] + sin(0.5 * th) * e1[j];
		mid[j] = 0.5 * (a[j] + b[j]);
		o[j] = mid[j] + 0.5 * r * u * (span < 3.14159 ? 1 : 0.2);
	}
	for (i = 0; i < k; i++)
	{
		unsigned int A0 = base + 2 + 2 * i, B0 = A0 + 1, A1 = A0 + 2, B1 = A0 + 3;
		if (!ofs_soup_tri_out(s, A0, A1, B1, o) || !ofs_soup_tri_out(s, A0, B1, B0, o) || !ofs_soup_tri_out(s, base, A0, A1, o)
			|| !ofs_soup_tri_out(s, base + 1, B0, B1, o))
			return 0;
	}
	return ofs_soup_tri_out(s, base, base + 1, base + 3, o) && ofs_soup_tri_out(s, base, base + 3, base + 2, o)
		&& ofs_soup_tri_out(s, base, base + 1, base + 3 + 2 * k, o) && ofs_soup_tri_out(s, base, base + 3 + 2 * k, base + 2 + 2 * k, o);
}

/* unit d turned away from unit c by angle a (d itself when along c) */
static void ofs_turn_away(double* out, double const* d, double const* c, double a)
{
	double t[3], dc = d[0] * c[0] + d[1] * c[1] + d[2] * c[2], l;
	unsigned int k;
	for (k = 0; k < 3; k++) t[k] = d[k] - dc * c[k];
	l = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
	/* (in the plane of c and d: d = cos(al) c + sin(al) t / l, turned to al + a) */
	a += atan2(l, dc);
	for (k = 0; k < 3; k++) out[k] = l > 1e-12 ? cos(a) * c[k] + sin(a) * t[k] / l : d[k];
}

/* unit direction from a to b, the fraction f of the way along the great arc */
static void ofs_slerp(double* out, double const* a, double const* b, double f)
{
	double cr[3], th, sn;
	unsigned int k;
	cr[0] = a[1] * b[2] - a[2] * b[1];
	cr[1] = a[2] * b[0] - a[0] * b[2];
	cr[2] = a[0] * b[1] - a[1] * b[0];
	th = atan2(sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]), a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	sn = sin(th);
	for (k = 0; k < 3; k++) out[k] = sn > 1e-12 ? (sin((1 - f) * th) * a[k] + sin(f * th) * b[k]) / sn : a[k];
}

static double ofs_angle(double const* a, double const* b)
{
	double cr[3];
	cr[0] = a[1] * b[2] - a[2] * b[1];
	cr[1] = a[2] * b[0] - a[0] * b[2];
	cr[2] = a[0] * b[1] - a[1] * b[0];
	return atan2(sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]), a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* a convex corner's piece: the cone from v over the directions in the spherical convex hull of its faces' normals D (n of
   them), widened by widen, radius r: the points within r of v nearest to v. 2 when the hull is too wide for a cone (the
   caller then puts a ball), 3 when it is a single arc (nothing to add), 0 on allocation failure */
static int ofs_cone(ofs_soup* s, double const* v, double const* D, unsigned int n, double widen, double r, double step)
{
	double c[3] = { 0, 0, 0 }, t1[3], t2[3], l, o[3], *q = NULL;
	unsigned int i, j, k, nh = 0, *ord = NULL, nb = 0, Q, base, *row = NULL;
	int ok = 0;
	double* g = (double*)qaws_internal_scratch_malloc(sizeof(double) * (2 * (size_t)n + 1));
	ord = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (2 * (size_t)n + 2));
	if (!g || !ord)
		goto done;
	for (i = 0; i < n; i++)
		for (k = 0; k < 3; k++) c[k] += D[3 * (size_t)i + k];
	l = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
	if (!(l > 1e-9))
	{
		ok = 2;
		goto done;
	}
	for (k = 0; k < 3; k++) c[k] /= l;
	/* the gnomonic projection on the plane touching at c: great circles lines, the spherical hull the plane's */
	t1[0] = fabs(c[0]) < 0.9 ? 1 : 0;
	t1[1] = 1 - t1[0];
	t1[2] = 0;
	l = t1[0] * c[0] + t1[1] * c[1];
	for (k = 0; k < 3; k++) t1[k] -= l * c[k];
	l = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
	for (k = 0; k < 3; k++) t1[k] /= l;
	t2[0] = c[1] * t1[2] - c[2] * t1[1];
	t2[1] = c[2] * t1[0] - c[0] * t1[2];
	t2[2] = c[0] * t1[1] - c[1] * t1[0];
	for (i = 0; i < n; i++)
	{
		double const* d = D + 3 * (size_t)i;
		double dc = d[0] * c[0] + d[1] * c[1] + d[2] * c[2];
		if (!(dc > 0.2))
		{
			ok = 2;   /* (a normal past about 78 degrees from the middle: no cone) */
			goto done;
		}
		g[2 * i] = (d[0] * t1[0] + d[1] * t1[1] + d[2] * t1[2]) / dc;
		g[2 * i + 1] = (d[0] * t2[0] + d[1] * t2[1] + d[2] * t2[2]) / dc;
		ord[i] = i;
	}
	/* (Andrew's monotone chain, counterclockwise in (t1, t2): the directions' order about c) */
	for (i = 1; i < n; i++)
		for (j = i; j > 0 && (g[2 * ord[j]] < g[2 * ord[j - 1]] || (g[2 * ord[j]] == g[2 * ord[j - 1]] && g[2 * ord[j] + 1] < g[2 * ord[j - 1] + 1])); j--)
		{
			unsigned int x = ord[j];
			ord[j] = ord[j - 1];
			ord[j - 1] = x;
		}
	row = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (2 * (size_t)n + 2));
	if (!row)
		goto done;
	for (j = 0; j < 2 * n - 1; j++)
	{
		unsigned int p = ord[j < n ? j : 2 * n - 2 - j], low = j < n ? 2 : nb + 1;
		if (j == n)
			nb = nh;   /* (the lower chain's length: the upper may not eat it) */
		low = j < n ? 2 : nb + 1;
		while (nh >= low && (g[2 * row[nh - 1]] - g[2 * row[nh - 2]]) * (g[2 * p + 1] - g[2 * row[nh - 2] + 1]) -
			(g[2 * row[nh - 1] + 1] - g[2 * row[nh - 2] + 1]) * (g[2 * p] - g[2 * row[nh - 2]]) <= 1e-12)
			nh--;
		row[nh++] = p;
	}
	nh--;   /* (the chain ends on its first point) */
	nb = 0;
	if (nh < 3)
	{
		ok = 3;
		goto done;
	}
	/* the hull's directions widened, their arcs sampled */
	q = (double*)qaws_internal_scratch_malloc(sizeof(double) * 3 * ((size_t)nh * ((size_t)ceil(3.2 / step) + 2) + 1));
	if (!q)
		goto done;
	for (i = 0; i < nh; i++)
	{
		double a0[3], a1[3];
		unsigned int m;
		ofs_turn_away(a0, D + 3 * (size_t)row[i], c, widen);
		ofs_turn_away(a1, D + 3 * (size_t)row[(i + 1) % nh], c, widen);
		m = (unsigned int)ceil(ofs_angle(a0, a1) / step);
		if (m < 1)
			m = 1;
		for (j = 0; j < m; j++) ofs_slerp(q + 3 * (size_t)nb++, a0, a1, (double)j / m);
	}
	/* the cap: rings from c out to the boundary's points, Q steps */
	for (i = 0, l = 0; i < nb; i++) l = fmax(l, ofs_angle(c, q + 3 * (size_t)i));
	Q = (unsigned int)ceil(l / step);
	if (Q < 1)
		Q = 1;
	base = s->nx;
	{
		double p[3];
		for (k = 0; k < 3; k++) p[k] = v[k];
		if (ofs_soup_point(s, p) == OFS_NONE)
			goto done;
		for (k = 0; k < 3; k++) p[k] = v[k] + r * c[k];
		if (ofs_soup_point(s, p) == OFS_NONE)
			goto done;
	}
	for (j = 1; j <= Q; j++)
		for (i = 0; i < nb; i++)
		{
			double d[3], p[3];
			ofs_slerp(d, c, q + 3 * (size_t)i, (double)j / Q);
			for (k = 0; k < 3; k++) p[k] = v[k] + r * d[k];
			if (ofs_soup_point(s, p) == OFS_NONE)
				goto done;
		}
	for (k = 0; k < 3; k++) o[k] = v[k] + 0.5 * r * c[k];
	for (i = 0; i < nb; i++)
	{
		unsigned int i1 = (i + 1) % nb, last = base + 2 + (Q - 1) * nb;
		if (!ofs_soup_tri_out(s, base + 1, base + 2 + i, base + 2 + i1, v) || !ofs_soup_tri_out(s, base, last + i, last + i1, o))
			goto done;
		for (j = 1; j < Q; j++)
		{
			unsigned int r0 = base + 2 + (j - 1) * nb, r1 = base + 2 + j * nb;
			if (!ofs_soup_tri_out(s, r0 + i, r1 + i, r1 + i1, v) || !ofs_soup_tri_out(s, r0 + i, r1 + i1, r0 + i1, v))
				goto done;
		}
	}
	ok = 1;
done:
	qaws_internal_free(g);
	qaws_internal_free(ord);
	qaws_internal_free(row);
	qaws_internal_free(q);
	return ok;
}

/* ROUND: the points within r of the triangles' surface on the given sides (1 along the normals, 2 against them, 3 both) as
   closed pieces whose union they are, a soup for the Boolean: each triangle's slab +-r; at each edge whose faces turn,
   the cylinder's wedge between their normals on the side it is convex to (the side it is concave to the slabs cover); a
   border or non-manifold edge's whole cylinder; at each vertex the cone over its normals' spherical hull where all its
   turning edges are convex to a side, a ball where they turn both ways or a cylinder ends there. Wedges and cones are
   widened by a quarter step so no two pieces share a plane. 0 on allocation failure */
static int ofs_pieces(ofs_soup* s, double const* P, unsigned int nv, unsigned int const* t, unsigned int nt, double r, double step, unsigned int sides)
{
	ofs_hkey* key = (ofs_hkey*)qaws_internal_scratch_malloc(sizeof(ofs_hkey) * (3 * (size_t)nt + 1));
	double* N = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nt + 1));
	unsigned char* vf = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nv + 1, 1);   /* 1, 2 convex to side 1, 2; 4, 8 concave; 16 ball */
	unsigned char* flat = (unsigned char*)qaws_internal_scratch_zalloc(3 * (size_t)nt + 1, 1);   /* per corner: its edge flat */
	unsigned int* mate = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));   /* the face across it */
	unsigned int *co = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 2, sizeof(unsigned int)), *cf = NULL, *fill = NULL;
	double* D = NULL;
	double widen = 0.25 * step;
	unsigned int f, i, k, m = (unsigned int)ceil(6.283185307179586 / step), side;
	int ok = 0;
	if (m < 8)
		m = 8;
	m += m & 1;
	if (!key || !N || !vf || !co || !flat || !mate)
		goto done;
	for (f = 0; f < nt; f++)
	{
		double const *a = P + 3 * (size_t)t[3 * f], *b = P + 3 * (size_t)t[3 * f + 1], *c = P + 3 * (size_t)t[3 * f + 2];
		double e1[3], e2[3], *n = N + 3 * (size_t)f, l;
		for (k = 0; k < 3; k++)
		{
			e1[k] = b[k] - a[k];
			e2[k] = c[k] - a[k];
		}
		n[0] = e1[1] * e2[2] - e1[2] * e2[1];
		n[1] = e1[2] * e2[0] - e1[0] * e2[2];
		n[2] = e1[0] * e2[1] - e1[1] * e2[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; k < 3; k++) n[k] = l > 0 ? n[k] / l : 0;
		for (k = 0; k < 3; k++)
		{
			unsigned int u = t[3 * f + k], w = t[3 * f + (k + 1) % 3];
			key[3 * f + k].lo = u < w ? u : w;
			key[3 * f + k].hi = u < w ? w : u;
			key[3 * f + k].f = f;
			co[t[3 * f + k] + 1]++;
		}
	}
	qsort(key, 3 * (size_t)nt, sizeof(ofs_hkey), ofs_hkey_cmp);
	for (i = 0; i < 3 * nt;)
	{
		unsigned int j = i, g, lo = key[i].lo, hi = key[i].hi, a, b, cg = 0;
		double const *nf, *ng;
		double dt, e0[3], e1[3], l, side_dot = 0;
		while (j < 3 * nt && key[j].lo == lo && key[j].hi == hi) j++;
		f = key[i].f;
		nf = N + 3 * (size_t)f;
		/* a and b as they run in f */
		for (k = 0; k < 3 && !(t[3 * f + k] == lo && (t[3 * f + (k + 1) % 3] == hi)); k++)
			;
		a = k < 3 ? lo : hi;
		b = k < 3 ? hi : lo;
		if (j - i != 2)
		{
			/* a border (or a non-manifold edge): its whole cylinder, balls at its ends */
			double u[3], e[3], d = 0, el = 0;
			for (k = 0; k < 3; k++) e[k] = P[3 * (size_t)b + k] - P[3 * (size_t)a + k];
			for (k = 0; k < 3; k++) { d += nf[k] * e[k]; el += e[k] * e[k]; }
			for (k = 0; k < 3; k++) u[k] = nf[k] - (el > 0 ? d / el : 0) * e[k];
			l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
			if (l > 0)
			{
				for (k = 0; k < 3; k++) u[k] /= l;
				if (!ofs_cylinder(s, P + 3 * (size_t)a, P + 3 * (size_t)b, u, r, m))
					goto done;
			}
			vf[a] |= 16;
			vf[b] |= 16;
			i = j;
			continue;
		}
		g = key[i + 1].f;
		ng = N + 3 * (size_t)g;
		dt = nf[0] * ng[0] + nf[1] * ng[1] + nf[2] * ng[2];
		if (dt >= 1 - 1e-12)
		{
			/* flat: one slab for both faces (below) covers it */
			unsigned int hf = OFS_NONE, hg = OFS_NONE;
			for (k = 0; k < 3; k++)
			{
				unsigned int u = t[3 * f + k], w = t[3 * f + (k + 1) % 3], x = t[3 * g + k], y = t[3 * g + (k + 1) % 3];
				if ((u == lo && w == hi) || (u == hi && w == lo))
					hf = 3 * f + k;
				if ((x == lo && y == hi) || (x == hi && y == lo))
					hg = 3 * g + k;
			}
			if (hf != OFS_NONE && hg != OFS_NONE)
			{
				flat[hf] = flat[hg] = 1;
				mate[hf] = g;
				mate[hg] = f;
			}
			i = j;
			continue;
		}
		/* convex to side 1 when g's third vertex lies behind f's plane */
		for (k = 0; k < 3; k++)
			if (t[3 * g + k] != lo && t[3 * g + k] != hi)
				cg = t[3 * g + k];
		for (k = 0; k < 3; k++) side_dot += nf[k] * (P[3 * (size_t)cg + k] - P[3 * (size_t)a + k]);
		side = side_dot < 0 ? 1u : 2u;
		vf[a] |= (unsigned char)(side == 1 ? 1 | 8 : 2 | 4);
		vf[b] |= (unsigned char)(side == 1 ? 1 | 8 : 2 | 4);
		if (!(sides & side))
		{
			i = j;
			continue;
		}
		/* the wedge from n_f to n_g on that side (both turned for side 2) */
		for (k = 0; k < 3; k++)
		{
			e0[k] = side == 1 ? nf[k] : -nf[k];
			e1[k] = (side == 1 ? ng[k] : -ng[k]) - dt * e0[k];
		}
		l = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
		if (l > 0)
		{
			for (k = 0; k < 3; k++) e1[k] /= l;
			if (!ofs_wedge(s, P + 3 * (size_t)a, P + 3 * (size_t)b, e0, e1, acos(dt < -1 ? -1 : dt), widen, r, step))
				goto done;
		}
		i = j;
	}
	{
		/* the slabs: the triangles joined across flat edges into planar groups, each group's points moved -r and +r along
		   its normal (a prism, its sides on the group's border only: no two slabs share a side) */
		unsigned int *up = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
		unsigned int *at = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nt + 2, sizeof(unsigned int));
		unsigned int *byg = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
		unsigned int *stamp = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nv + 1));
		unsigned int *pt = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (2 * (size_t)nv + 1)), q;
		int fine = up && at && byg && stamp && pt;
		if (fine)
		{
			for (f = 0; f < nt; f++) up[f] = f;
			for (f = 0; f < 3 * nt; f++)
				if (flat[f])
				{
					unsigned int x = ofs_root(up, f / 3), y = ofs_root(up, mate[f]);
					if (x != y)
						up[x] = y;
				}
			/* (the faces by their group: counted, then placed) */
			for (f = 0; f < nt; f++)
			{
				up[f] = ofs_root(up, f);
				at[up[f] + 1]++;
			}
			for (f = 0; f < nt; f++) at[f + 1] += at[f];
			for (f = 0; f < nt; f++) byg[at[up[f]]++] = f;
			for (i = 0; i < nv; i++) stamp[i] = OFS_NONE;
		}
		for (q = 0; fine && q < nt; q++)
		{
			unsigned int ff = byg[q], gr = up[ff], v3[6];
			double const* n = N + 3 * (size_t)gr;
			if (n[0] == 0 && n[1] == 0 && n[2] == 0)
				continue;
			for (k = 0; k < 3; k++)
			{
				unsigned int v = t[3 * ff + k];
				if (stamp[v] != gr)
				{
					double x[3];
					unsigned int c2;
					stamp[v] = gr;
					for (c2 = 0; c2 < 3; c2++) x[c2] = P[3 * (size_t)v + c2] - r * n[c2];
					pt[2 * v] = ofs_soup_point(s, x);
					for (c2 = 0; c2 < 3; c2++) x[c2] = P[3 * (size_t)v + c2] + r * n[c2];
					pt[2 * v + 1] = ofs_soup_point(s, x);
				}
				v3[k] = pt[2 * v];
				v3[3 + k] = pt[2 * v + 1];
			}
			fine = ofs_soup_tri(s, v3[3], v3[4], v3[5]) && ofs_soup_tri(s, v3[0], v3[2], v3[1]);
			for (k = 0; k < 3 && fine; k++)
				if (!flat[3 * ff + k])
					fine = ofs_soup_tri(s, v3[k], v3[(k + 1) % 3], v3[3 + (k + 1) % 3]) && ofs_soup_tri(s, v3[k], v3[3 + (k + 1) % 3], v3[3 + k]);
		}
		qaws_internal_free(up);
		qaws_internal_free(at);
		qaws_internal_free(byg);
		qaws_internal_free(stamp);
		qaws_internal_free(pt);
		if (!fine)
			goto done;
	}
	/* the vertices: their faces' normals */
	cf = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	fill = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nv + 1, sizeof(unsigned int));
	D = (double*)qaws_internal_scratch_malloc(sizeof(double) * (9 * (size_t)nt + 1));
	if (!cf || !fill || !D)
		goto done;
	for (i = 0; i < nv; i++) co[i + 1] += co[i];
	for (f = 0; f < nt; f++)
		for (k = 0; k < 3; k++) cf[co[t[3 * f + k]] + fill[t[3 * f + k]]++] = f;
	for (i = 0; i < nv; i++)
	{
		int ball = (vf[i] & 16) != 0;
		for (side = 1; side <= 2 && !ball; side++)
		{
			unsigned int q, nd = 0;
			int res;
			if (!(sides & side) || !(vf[i] & side))
				continue;   /* (no edge convex to this side: nothing nearest to the vertex there) */
			if (vf[i] & (side == 1 ? 4 : 8))
			{
				ball = 1;   /* (a saddle: turning both ways) */
				break;
			}
			for (q = co[i]; q < co[i + 1]; q++, nd++)
				for (k = 0; k < 3; k++) D[3 * (size_t)nd + k] = (side == 1 ? 1 : -1) * N[3 * (size_t)cf[q] + k];
			res = ofs_cone(s, P + 3 * (size_t)i, D, nd, widen, r, step);
			if (!res)
				goto done;
			ball = res == 2;   /* (3: the normals on one great circle, a straight crease the wedges cover) */
		}
		if (ball && !ofs_ball(s, P + 3 * (size_t)i, r, m))
			goto done;
	}
	ok = 1;
done:
	qaws_internal_free(key);
	qaws_internal_free(N);
	qaws_internal_free(vf);
	qaws_internal_free(flat);
	qaws_internal_free(mate);
	qaws_internal_free(co);
	qaws_internal_free(cf);
	qaws_internal_free(fill);
	qaws_internal_free(D);
	return ok;
}

/* ROUND: SURFACE outward the input united with the pieces at d, inward the input less the pieces at |d|; SHELL the input
   meeting the pieces at |d|; THICKEN the pieces at |d| / 2 alone. The pieces' circles split every step at most */
static qaws_status ofs_round(qaws_mesh_3d** out, qaws_boolean_report* br, double const* P, qaws_scalar const* ps, unsigned int nv,
	unsigned int const* t, unsigned int nt, qaws_mesh_offset_kind kind, double d, double step, qaws_intersect_desc const* id)
{
	ofs_soup s;
	unsigned int i;
	double r = kind == QAWS_MESH_OFFSET_THICKEN ? 0.5 * fabs(d) : fabs(d);
	qaws_scalar* xs = NULL;
	qaws_mesh_3d_desc a, b;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&s, 0, sizeof s);
	if (!ofs_pieces(&s, P, nv, t, nt, r, step, kind == QAWS_MESH_OFFSET_THICKEN ? 3u : d > 0 ? 1u : 2u))
		goto done;
	xs = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)s.nx + 1));
	if (!xs)
		goto done;
	for (i = 0; i < 3 * s.nx; i++) xs[i] = (qaws_scalar)s.x[i];
	memset(&a, 0, sizeof a);
	a.positions = ps;
	a.vertex_count = nv;
	a.face_vertices = t;
	a.face_count = nt;
	memset(&b, 0, sizeof b);
	b.positions = xs;
	b.vertex_count = s.nx;
	b.face_vertices = s.t;
	b.face_count = s.nt;
	if (kind == QAWS_MESH_OFFSET_THICKEN)
		st = qaws_mesh_boolean(out, br, &b, NULL, QAWS_BOOLEAN_UNION, id);
	else
		st = qaws_mesh_boolean(out, br, &a, &b, kind == QAWS_MESH_OFFSET_SHELL ? QAWS_BOOLEAN_INTERSECTION : d > 0 ? QAWS_BOOLEAN_UNION : QAWS_BOOLEAN_DIFFERENCE, id);
done:
	qaws_internal_free(s.x);
	qaws_internal_free(s.t);
	qaws_internal_free(xs);
	return st;
}

/* SHELL: the input's triangles, then the inner surface c's turned inward, into out (c destroyed) */
static qaws_status ofs_shell(qaws_mesh_3d** out, qaws_scalar const* ps, unsigned int nv, unsigned int const* t, unsigned int nt, qaws_mesh_3d* c)
{
	unsigned int cv = qaws_mesh_3d_get_vertex_count(c), ct = qaws_mesh_3d_get_triangle_count(c), i;
	qaws_mesh_3d_desc md;
	qaws_scalar* Ys = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + cv + 1));
	unsigned int* T = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + ct + 1));
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (Ys && T)
	{
		memcpy(Ys, ps, sizeof(qaws_scalar) * 3 * nv);
		qaws_mesh_3d_get_positions(Ys + 3 * (size_t)nv, c);
		memcpy(T, t, sizeof(unsigned int) * 3 * nt);
		qaws_mesh_3d_get_triangles(T + 3 * (size_t)nt, c);
		for (i = nt; i < nt + ct; i++)
		{
			unsigned int a = T[3 * i + 1];
			T[3 * i] += nv; T[3 * i + 1] = T[3 * i + 2] + nv; T[3 * i + 2] = a + nv;
		}
		memset(&md, 0, sizeof md);
		md.positions = Ys;
		md.vertex_count = nv + cv;
		md.face_vertices = T;
		md.face_count = nt + ct;
		st = qaws_mesh_3d_create(out, &md);
	}
	qaws_mesh_3d_destroy(c);
	qaws_internal_free(Ys);
	qaws_internal_free(T);
	return st;
}

qaws_status qaws_mesh_3d_offset(qaws_mesh_3d** out, qaws_mesh_offset_report* report, qaws_mesh_3d const* mesh, qaws_mesh_offset_desc const* desc)
{
	qaws_mesh_offset_report local;
	ofs_field f;
	qaws_intersect_desc idesc;
	unsigned int nv, nt, i, k;
	double d, lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
	double *P = NULL, *m = NULL, *Y = NULL;
	unsigned int *t = NULL, *T = NULL, *bd = NULL;
	qaws_scalar *ps = NULL, *Ys = NULL;
	qaws_status st;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh || !desc || !(qaws_mesh_3d_get_triangle_count(mesh) > 0) || !isfinite((double)desc->distance) || desc->distance == 0 ||
		(unsigned int)desc->kind > QAWS_MESH_OFFSET_THICKEN || (unsigned int)desc->method > QAWS_MESH_OFFSET_ROUND ||
		!(desc->round_angle == 0 || (desc->round_angle >= 0.01 && desc->round_angle <= 1.6)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	d = (double)desc->distance;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	if (desc->distances)
	{
		/* (MOVE only: a contour of a varying distance needs a level per point) */
		if (desc->method != QAWS_MESH_OFFSET_MOVE)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (i = 0; i < nv; i++)
			if (!isfinite((double)desc->distances[i]))
				return QAWS_STATUS_INVALID_ARGUMENT;
	}
	memset(&f, 0, sizeof f);
	P = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	ps = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	if (!P || !ps || !t)
		goto oom;
	qaws_mesh_3d_get_positions(ps, mesh);
	qaws_mesh_3d_get_triangles(t, mesh);
	for (i = 0; i < 3 * nv; i++) P[i] = (double)ps[i];
	for (i = 0; i < 3 * nt; i++)
		for (k = 0; k < 3; k++) { lo[k] = fmin(lo[k], P[3 * (size_t)t[i] + k]); hi[k] = fmax(hi[k], P[3 * (size_t)t[i] + k]); }
	if ((st = qaws_mesh_bvh_create_mesh_3d(&f.b, mesh, NULL)) != QAWS_STATUS_OK)
		goto done;
	f.P = P;
	f.sgn = desc->kind != QAWS_MESH_OFFSET_THICKEN;
	f.ex = desc->executor;
	memset(&idesc, 0, sizeof idesc);
	idesc.executor = desc->executor;
	if (desc->method == QAWS_MESH_OFFSET_CONTOUR)
	{
		unsigned int res = desc->resolution ? desc->resolution : 64;
		qaws_mesh_3d* c = NULL;
		wt_contour_info ci;
		/* SURFACE: phi = s - d; SHELL: the inner surface s = -|d|; THICKEN: |dist| = |d| / 2 */
		f.level = desc->kind == QAWS_MESH_OFFSET_SURFACE ? d : desc->kind == QAWS_MESH_OFFSET_SHELL ? -fabs(d) : 0.5 * fabs(d);
		memset(&ci, 0, sizeof ci);
		st = wt_contour_field(&c, &ci, ofs_field_eval, &f, lo, hi, f.level > 0 ? f.level : 0, res, 1);
		memcpy(report->grid, ci.grid, sizeof ci.grid);
		report->cell = (qaws_scalar)ci.cell;
		report->crossings = ci.crossings;
		report->sharp = ci.sharp;
		report->closed = ci.closed;
		report->evaluated = (unsigned int)ci.evaluated;
		if (st == QAWS_STATUS_OK && c)
			st = ofs_measure_mesh(report, &f, c);
		if (st != QAWS_STATUS_OK || !c || desc->kind != QAWS_MESH_OFFSET_SHELL)
		{
			if (st != QAWS_STATUS_OK)
				qaws_mesh_3d_destroy(c);
			else
				*out = c;
			goto done;
		}
		st = ofs_shell(out, ps, nv, t, nt, c);
		goto done;
	}
	if (desc->method == QAWS_MESH_OFFSET_ROUND)
	{
		/* SURFACE at d, SHELL's inner surface at -|d| (then joined to the input), THICKEN within |d| / 2, the border rounded */
		qaws_mesh_3d* c = NULL;
		qaws_boolean_report br;
		memset(&br, 0, sizeof br);
		st = ofs_round(&c, &br, P, ps, nv, t, nt, desc->kind == QAWS_MESH_OFFSET_THICKEN ? QAWS_MESH_OFFSET_THICKEN : QAWS_MESH_OFFSET_SURFACE,
			desc->kind == QAWS_MESH_OFFSET_SHELL ? -fabs(d) : d, desc->round_angle > 0 ? (double)desc->round_angle : 3.14159265358979323846 / 8, &idesc);
		report->intersect = br.intersect;
		report->closed = br.closed;
		if (st == QAWS_STATUS_OK && c)
			st = ofs_measure_mesh(report, &f, c);
		if (st != QAWS_STATUS_OK || !c || desc->kind != QAWS_MESH_OFFSET_SHELL)
		{
			if (st != QAWS_STATUS_OK)
				qaws_mesh_3d_destroy(c);
			else
				*out = c;
			goto done;
		}
		st = ofs_shell(out, ps, nv, t, nt, c);
		goto done;
	}
	/* MOVE: the mitred vertices; SURFACE one copy at d, SHELL the input and a copy at -|d| turned, THICKEN copies at +-|d|/2 */
	{
		unsigned int copies = desc->kind == QAWS_MESH_OFFSET_SURFACE ? 1 : 2, nb = 0, nf, c;
		double sgn0 = desc->kind == QAWS_MESH_OFFSET_SURFACE ? d : desc->kind == QAWS_MESH_OFFSET_SHELL ? 0 : 0.5 * fabs(d);
		double sgn1 = desc->kind == QAWS_MESH_OFFSET_SHELL ? -fabs(d) : -0.5 * fabs(d);
		qaws_mesh_3d_desc md;
		qaws_boolean_report br;
		m = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
		if (!m)
			goto oom;
		ofs_mitre(m, P, nv, t, nt);
		if (desc->kind == QAWS_MESH_OFFSET_THICKEN && (nb = ofs_border(&bd, t, nt)) == ~0u)
			goto oom;
		nf = copies * nt + 2 * nb;
		Y = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)copies * nv + 1));
		Ys = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)copies * nv + 1));
		T = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nf + 1));
		if (!Y || !Ys || !T)
			goto oom;
		for (c = 0; c < copies; c++)
			for (i = 0; i < nv; i++)
				for (k = 0; k < 3; k++)
				{
					/* (a distance per vertex: the copies' at its own, as the uniform ones at d) */
					double s = c ? sgn1 : sgn0;
					if (desc->distances)
						s = s / (desc->kind == QAWS_MESH_OFFSET_SURFACE ? d : fabs(d)) * (desc->kind == QAWS_MESH_OFFSET_SURFACE ? (double)desc->distances[i] : fabs((double)desc->distances[i]));
					Y[3 * ((size_t)c * nv + i) + k] = P[3 * (size_t)i + k] + s * m[3 * (size_t)i + k];
				}
		for (i = 0; i < 3 * copies * nv; i++) Ys[i] = (qaws_scalar)Y[i];
		memcpy(T, t, sizeof(unsigned int) * 3 * nt);
		for (i = 0; copies == 2 && i < nt; i++)
		{
			/* the second copy turned */
			T[3 * (nt + i)] = t[3 * i] + nv; T[3 * (nt + i) + 1] = t[3 * i + 2] + nv; T[3 * (nt + i) + 2] = t[3 * i + 1] + nv;
		}
		for (i = 0; i < nb; i++)
		{
			/* the rim: the border edge a -> b on the first copy, a quad b, a, a', b' */
			unsigned int a = bd[2 * i], b = bd[2 * i + 1], *q = T + 3 * ((size_t)2 * nt + 2 * i);
			q[0] = b; q[1] = a; q[2] = a + nv;
			q[3] = b; q[4] = a + nv; q[5] = b + nv;
		}
		/* the offset surface measured as moved (the shell's input copy aside) */
		if (desc->kind == QAWS_MESH_OFFSET_SHELL)
			st = ofs_measure(report, &f, Y + 3 * (size_t)nv, nv);
		else
			st = ofs_measure(report, &f, Y, (size_t)copies * nv);
		if (st != QAWS_STATUS_OK)
			goto done;
		memset(&md, 0, sizeof md);
		md.positions = Ys;
		md.vertex_count = copies * nv;
		md.face_vertices = T;
		md.face_count = nf;
		memset(&br, 0, sizeof br);
		st = qaws_mesh_boolean(out, &br, &md, NULL, QAWS_BOOLEAN_UNION, &idesc);
		report->intersect = br.intersect;
		report->closed = br.closed;
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_mesh_bvh_destroy(f.b);
	qaws_internal_free(P); qaws_internal_free(ps); qaws_internal_free(t); qaws_internal_free(m); qaws_internal_free(Y); qaws_internal_free(Ys);
	qaws_internal_free(T); qaws_internal_free(bd);
	return st;
}
