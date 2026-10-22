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

typedef struct ofs_field
{
	qaws_mesh_bvh* b;
	double const* P;               /* the input's positions */
	int sgn;                       /* signed distance (SURFACE, SHELL) or unsigned (THICKEN) */
	double level;                  /* phi = distance - level */
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
		st = f->sgn ? qaws_mesh_bvh_signed_distance(dist, hits, f->b, xs, m, NULL) : qaws_mesh_bvh_closest(hits, f->b, xs, m, 0, NULL);
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

/* the level set phi = 0 contoured on a grid over the box lo..hi grown by margin; *out NULL when it is empty */
static qaws_status ofs_contour(qaws_mesh_3d** out, qaws_mesh_offset_report* report, ofs_field const* f, double const* lo, double const* hi,
	double margin, unsigned int res)
{
	double glo[3], h, ext = 0, *gx = NULL, *phi = NULL, *A = NULL, *B = NULL, *pa = NULL, *pb = NULL, *xm = NULL, *fm = NULL, *gm = NULL;
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
	report->cell = (qaws_scalar)h;
	N = (size_t)n[0] * n[1] * n[2];
	gx = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, N));
	phi = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), N));
	in = (unsigned char*)qaws_internal_scratch_malloc(N);
	if (!gx || !phi || !in)
		goto oom;
	for (q = 0; q < N; q++)
	{
		gx[3 * q] = glo[0] + h * (double)(q % n[0]);
		gx[3 * q + 1] = glo[1] + h * (double)(q / n[0] % n[1]);
		gx[3 * q + 2] = glo[2] + h * (double)(q / ((size_t)n[0] * n[1]));
	}
	if ((st = ofs_eval(f, gx, N, phi, NULL)) != QAWS_STATUS_OK)
		goto done;
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
		memcpy(A + 3 * e, gx + 3 * qi, sizeof(double) * 3);
		memcpy(B + 3 * e, gx + 3 * qo, sizeof(double) * 3);
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
		if ((st = ofs_eval(f, xm, nc, fm, NULL)) != QAWS_STATUS_OK)
			goto done;
		for (e = 0; e < nc; e++)
		{
			if (fm[e] < 0) { memcpy(A + 3 * e, xm + 3 * e, sizeof(double) * 3); pa[e] = fm[e]; }
			else { memcpy(B + 3 * e, xm + 3 * e, sizeof(double) * 3); pb[e] = fm[e]; }
		}
	}
	if ((st = ofs_eval(f, xm, nc, fm, gm)) != QAWS_STATUS_OK)
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
	qaws_internal_free(gx); qaws_internal_free(phi); qaws_internal_free(in); qaws_internal_free(X);
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

qaws_status qaws_mesh_3d_offset(qaws_mesh_3d** out, qaws_mesh_offset_report* report, qaws_mesh_3d const* mesh, qaws_mesh_offset_desc const* desc)
{
	qaws_mesh_offset_report local;
	ofs_field f;
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
		(unsigned int)desc->kind > QAWS_MESH_OFFSET_THICKEN || (unsigned int)desc->method > QAWS_MESH_OFFSET_MOVE)
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
	if (desc->method == QAWS_MESH_OFFSET_CONTOUR)
	{
		unsigned int res = desc->resolution ? desc->resolution : 64;
		qaws_mesh_3d* c = NULL;
		/* SURFACE: phi = s - d; SHELL: the inner surface s = -|d|; THICKEN: |dist| = |d| / 2 */
		f.level = desc->kind == QAWS_MESH_OFFSET_SURFACE ? d : desc->kind == QAWS_MESH_OFFSET_SHELL ? -fabs(d) : 0.5 * fabs(d);
		st = ofs_contour(&c, report, &f, lo, hi, f.level > 0 ? f.level : 0, res);
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
		{
			/* the shell: the input's triangles, then the inner surface's turned inward */
			unsigned int cv = qaws_mesh_3d_get_vertex_count(c), ct = qaws_mesh_3d_get_triangle_count(c);
			qaws_mesh_3d_desc md;
			Ys = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + cv + 1));
			T = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + ct + 1));
			if (!Ys || !T)
			{
				qaws_mesh_3d_destroy(c);
				goto oom;
			}
			memcpy(Ys, ps, sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_positions(Ys + 3 * (size_t)nv, c);
			memcpy(T, t, sizeof(unsigned int) * 3 * nt);
			qaws_mesh_3d_get_triangles(T + 3 * (size_t)nt, c);
			for (i = nt; i < nt + ct; i++)
			{
				unsigned int a = T[3 * i + 1];
				T[3 * i] += nv; T[3 * i + 1] = T[3 * i + 2] + nv; T[3 * i + 2] = a + nv;
			}
			qaws_mesh_3d_destroy(c);
			memset(&md, 0, sizeof md);
			md.positions = Ys;
			md.vertex_count = nv + cv;
			md.face_vertices = T;
			md.face_count = nt + ct;
			st = qaws_mesh_3d_create(out, &md);
		}
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
		st = qaws_mesh_boolean(out, &br, &md, NULL, QAWS_BOOLEAN_UNION, NULL);
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
