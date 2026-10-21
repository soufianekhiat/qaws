/*
 * Decimation (qaws_mesh_decimate.h): quadrics per vertex, a heap of edges
 * by their collapse's error (stale entries skipped by the vertices'
 * stamps), collapses on the corner store without compacting until the end.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_decimate.h"
#include "../qaws_mesh_reference.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* a quadric: the symmetric 4 x 4 [a b c d]^T [a b c d] summed, its 10 entries a2 ab ac ad b2 bc bd c2 cd d2 */
static void dc_plane(double* Q, double const* n, double d, double w)
{
	double p[4] = { n[0], n[1], n[2], d };
	unsigned int i, j, k = 0;
	for (i = 0; i < 4; i++)
		for (j = i; j < 4; j++)
			Q[k++] += w * p[i] * p[j];
}

static double dc_eval(double const* Q, double const* x)
{
	return Q[0] * x[0] * x[0] + 2 * Q[1] * x[0] * x[1] + 2 * Q[2] * x[0] * x[2] + 2 * Q[3] * x[0] + Q[4] * x[1] * x[1] + 2 * Q[5] * x[1] * x[2]
		+ 2 * Q[6] * x[1] + Q[7] * x[2] * x[2] + 2 * Q[8] * x[2] + Q[9];
}

/* the point where Q is least, by Cramer's rule; 0 when Q's 3 x 3 part is near singular (flat, or along a line) */
static int dc_optimal(double const* Q, double* x)
{
	double a = Q[0], b = Q[1], c = Q[2], e = Q[4], f = Q[5], i = Q[7], r[3] = { -Q[3], -Q[6], -Q[8] };
	double det = a * (e * i - f * f) - b * (b * i - f * c) + c * (b * f - e * c), s = fabs(a) + fabs(e) + fabs(i);
	if (!(fabs(det) > 1e-10 * s * s * s))
		return 0;
	x[0] = (r[0] * (e * i - f * f) - b * (r[1] * i - f * r[2]) + c * (r[1] * f - e * r[2])) / det;
	x[1] = (a * (r[1] * i - f * r[2]) - r[0] * (b * i - f * c) + c * (b * r[2] - r[1] * c)) / det;
	x[2] = (a * (e * r[2] - r[1] * f) - b * (b * r[2] - r[1] * c) + r[0] * (b * f - e * c)) / det;
	return 1;
}

typedef struct dc_entry
{
	double cost;
	unsigned int a, b, sa, sb;
} dc_entry;

typedef struct dc_heap
{
	dc_entry* e;
	unsigned int n, cap;
} dc_heap;

static int dc_push(dc_heap* H, dc_entry const* x)
{
	unsigned int i;
	if (H->n == H->cap)
	{
		unsigned int cap = H->cap ? 2 * H->cap : 1024;
		dc_entry* g = (dc_entry*)qaws_internal_realloc(H->e, qaws_internal_mul2(sizeof(dc_entry), cap));
		if (!g)
			return 0;
		H->e = g;
		H->cap = cap;
	}
	for (i = H->n++; i && H->e[(i - 1) / 2].cost > x->cost; i = (i - 1) / 2)
		H->e[i] = H->e[(i - 1) / 2];
	H->e[i] = *x;
	return 1;
}

static void dc_pop(dc_heap* H, dc_entry* out)
{
	dc_entry last;
	unsigned int i = 0;
	*out = H->e[0];
	last = H->e[--H->n];
	for (;;)
	{
		unsigned int c = 2 * i + 1;
		if (c >= H->n)
			break;
		if (c + 1 < H->n && H->e[c + 1].cost < H->e[c].cost)
			c++;
		if (H->e[c].cost >= last.cost)
			break;
		H->e[i] = H->e[c];
		i = c;
	}
	if (H->n)
		H->e[i] = last;
}

typedef struct dc_work
{
	qaws_hm* m;
	double *Q, *W;
	unsigned int* stamp;
	unsigned char* border;
	int lock;
} dc_work;

/* a merged into b: allowed, and where b goes and the error there; 0 when refused */
static int dc_candidate(dc_work const* w, unsigned int a, unsigned int b, double* x, double* cost)
{
	qaws_hm const* m = w->m;
	unsigned int c = hm_find_edge(m, a, b), k;
	double Q[10], best = HUGE_VAL, tries[4][3], wsum = w->W[a] + w->W[b];
	int fixed, nt = 0, i;
	if (c == HM_NONE || (m->vflag[a] & HM_CORNER) || (hm_feature_degree(m, a) && !(m->cflag[c] & HM_FEATURE)) || (w->lock && w->border[a])
		|| !hm_can_collapse(m, c, a))
		return 0;
	for (k = 0; k < 10; k++)
		Q[k] = w->Q[10 * (size_t)a + k] + w->Q[10 * (size_t)b + k];
	/* b stays where it is when it may not leave its place: a corner, a locked border vertex, a feature vertex the edge leaves */
	fixed = (m->vflag[b] & HM_CORNER) || (w->lock && w->border[b]) || (hm_feature_degree(m, b) && !(m->cflag[c] & HM_FEATURE));
	memcpy(tries[nt++], m->x + 3 * (size_t)b, 3 * sizeof(double));
	if (!fixed)
	{
		memcpy(tries[nt++], m->x + 3 * (size_t)a, 3 * sizeof(double));
		for (k = 0; k < 3; k++)
			tries[nt][k] = 0.5 * (m->x[3 * (size_t)a + k] + m->x[3 * (size_t)b + k]);
		nt++;
		if (dc_optimal(Q, tries[nt]))
			nt++;
	}
	for (i = 0; i < nt; i++)
	{
		double e = dc_eval(Q, tries[i]);
		if (e < best)
		{
			best = e;
			memcpy(x, tries[i], 3 * sizeof(double));
		}
	}
	*cost = (best > 0 ? best : 0) / (wsum > 0 ? wsum : 1);
	return 1;
}

/* no face round a or b (but the edge's) turned over or flattened by a and b moving to x */
static int dc_no_fold(qaws_hm const* m, unsigned int a, unsigned int b, double const* x)
{
	unsigned int s, v, c;
	for (s = 0; s < 2; s++)
	{
		v = s ? b : a;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			unsigned int f = m->cf[c], e = m->fc[f] + m->fn[f], q, hasa = 0, hasb = 0;
			double n0[3] = { 0, 0, 0 }, n1[3] = { 0, 0, 0 }, l0, l1;
			for (q = m->fc[f]; q < e; q++)
			{
				hasa |= m->cv[q] == a;
				hasb |= m->cv[q] == b;
			}
			if (hasa && hasb)
				continue;   /* (the edge's own faces go) */
			for (q = m->fc[f]; q < e; q++)
			{
				unsigned int r = q + 1 < e ? q + 1 : m->fc[f];
				double const *p = m->x + 3 * (size_t)m->cv[q], *pr = m->x + 3 * (size_t)m->cv[r];
				double const *P = m->cv[q] == a || m->cv[q] == b ? x : p, *PR = m->cv[r] == a || m->cv[r] == b ? x : pr;
				n0[0] += (p[1] - pr[1]) * (p[2] + pr[2]);
				n0[1] += (p[2] - pr[2]) * (p[0] + pr[0]);
				n0[2] += (p[0] - pr[0]) * (p[1] + pr[1]);
				n1[0] += (P[1] - PR[1]) * (P[2] + PR[2]);
				n1[1] += (P[2] - PR[2]) * (P[0] + PR[0]);
				n1[2] += (P[0] - PR[0]) * (P[1] + PR[1]);
			}
			l0 = sqrt(n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
			l1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
			/* (turned by more than about 78 degrees, or shrunk to a sliver) */
			if (!(l1 > 1e-6 * l0) || n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < 0.2 * l0 * l1)
				return 0;
		}
	}
	return 1;
}

static int dc_push_edges(dc_work const* w, dc_heap* H, unsigned int v, unsigned int** ring, unsigned int* cap)
{
	unsigned int n = hm_ring(w->m, v, *ring, *cap), i;
	if (n > *cap)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(*ring, sizeof(unsigned int) * (size_t)n);
		if (!g)
			return 0;
		*ring = g;
		*cap = n;
		hm_ring(w->m, v, *ring, *cap);
	}
	for (i = 0; i < n; i++)
	{
		unsigned int u = (*ring)[i], s;
		for (s = 0; s < 2; s++)
		{
			dc_entry e;
			double x[3];
			e.a = s ? u : v;
			e.b = s ? v : u;
			if (!dc_candidate(w, e.a, e.b, x, &e.cost))
				continue;
			e.sa = w->stamp[e.a];
			e.sb = w->stamp[e.b];
			if (!dc_push(H, &e))
				return 0;
		}
	}
	return 1;
}

static unsigned int dc_faces(qaws_hm const* m)
{
	unsigned int f, n = 0;
	for (f = 0; f < m->nf; f++)
		n += !(m->fflag[f] & HM_DELETED);
	return n;
}

qaws_status qaws_mesh_3d_decimate(qaws_decimate_report* report, qaws_mesh_3d* mesh, qaws_decimate_desc const* desc)
{
	qaws_decimate_report local;
	qaws_mesh_reference* ref = NULL;
	qaws_hm* m;
	dc_work w;
	dc_heap H;
	unsigned int *ring = NULL, rcap = 64, f, v, k, c, live, target;
	double bound2;
	int failed = 0;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!mesh || (desc && !(desc->max_error >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	memset(&w, 0, sizeof w);
	memset(&H, 0, sizeof H);
	live = dc_faces(m);
	report->faces_before = report->faces_after = live;
	target = desc ? desc->target_faces : live / 2;
	bound2 = desc && desc->max_error > 0 ? (double)desc->max_error * (double)desc->max_error : HUGE_VAL;
	if (desc && (desc->flags & QAWS_DECIMATE_MEASURE))
	{
		/* the input kept as a LINEAR reference for the deviation */
		qaws_mesh_3d* copy = NULL;
		qaws_mesh_reference_desc rd;
		memset(&rd, 0, sizeof rd);
		rd.interpolant = QAWS_INTERPOLANT_LINEAR;
		st = qaws_mesh_3d_clone(&copy, mesh);
		if (st == QAWS_STATUS_OK)
			st = qaws_mesh_reference_create(&ref, copy, &rd);
		qaws_mesh_3d_destroy(copy);
		if (st != QAWS_STATUS_OK)
			return st;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	}
	w.m = m;
	w.lock = desc && (desc->flags & QAWS_DECIMATE_LOCK_BORDER);
	w.Q = (double*)qaws_internal_zalloc(qaws_internal_mul2(10, (size_t)m->nv + 1), sizeof(double));
	w.W = (double*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(double));
	w.stamp = (unsigned int*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(unsigned int));
	w.border = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	ring = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * rcap);
	if (!w.Q || !w.W || !w.stamp || !w.border || !ring)
		goto done;
	for (v = 0; v < m->nv; v++)
		w.border[v] = !(m->vflag[v] & HM_DELETED) && hm_is_border_vertex(m, v);
	/* each face's plane on its corners, by area; planes standing on the border's and the features' edges, heavier */
	for (f = 0; f < m->nf; f++)
	{
		double n[3], a, d;
		unsigned int e;
		if (m->fflag[f] & HM_DELETED)
			continue;
		hm_face_normal(m, f, n);
		a = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (!(a > 0))
			continue;
		for (k = 0; k < 3; k++)
			n[k] /= a;
		d = -(n[0] * m->x[3 * (size_t)m->cv[m->fc[f]]] + n[1] * m->x[3 * (size_t)m->cv[m->fc[f]] + 1] + n[2] * m->x[3 * (size_t)m->cv[m->fc[f]] + 2]);
		e = m->fc[f] + m->fn[f];
		for (c = m->fc[f]; c < e; c++)
		{
			dc_plane(w.Q + 10 * (size_t)m->cv[c], n, d, a);
			w.W[m->cv[c]] += a;
			if (m->co[c] == HM_NONE || (m->cflag[c] & HM_FEATURE))
			{
				double const *p = m->x + 3 * (size_t)m->cv[c], *q = m->x + 3 * (size_t)hm_to(m, c);
				double u[3] = { q[0] - p[0], q[1] - p[1], q[2] - p[2] }, s[3], l2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2], ls;
				s[0] = u[1] * n[2] - u[2] * n[1];
				s[1] = u[2] * n[0] - u[0] * n[2];
				s[2] = u[0] * n[1] - u[1] * n[0];
				ls = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
				if (ls > 0)
				{
					for (k = 0; k < 3; k++)
						s[k] /= ls;
					{
					double ds = -(s[0] * p[0] + s[1] * p[1] + s[2] * p[2]);
					dc_plane(w.Q + 10 * (size_t)m->cv[c], s, ds, 100 * l2);
					dc_plane(w.Q + 10 * (size_t)hm_to(m, c), s, ds, 100 * l2);
					}
				}
			}
		}
	}
	for (v = 0; v < m->nv; v++)
		if (!(m->vflag[v] & HM_DELETED) && !dc_push_edges(&w, &H, v, &ring, &rcap))
			goto done;
	while (H.n && live > target)
	{
		dc_entry e;
		double x[3], cost, t, l2;
		unsigned int a, b, gone = 0;
		dc_pop(&H, &e);
		a = e.a;
		b = e.b;
		if ((m->vflag[a] & HM_DELETED) || (m->vflag[b] & HM_DELETED) || e.sa != w.stamp[a] || e.sb != w.stamp[b])
			continue;   /* (stale) */
		if (e.cost > bound2)
			break;
		if (!dc_candidate(&w, a, b, x, &cost) || !dc_no_fold(m, a, b, x))
			continue;
		/* the faces of the edge that go: its triangles */
		for (c = m->vc[a]; c != HM_NONE; c = hm_round(m, a, c))
		{
			unsigned int g = m->cf[c], q, end = m->fc[g] + m->fn[g];
			for (q = m->fc[g]; q < end; q++)
				gone += m->fn[g] == 3 && m->cv[q] == b;
		}
		l2 = 0;
		t = 0;
		for (k = 0; k < 3; k++)
		{
			double ab = m->x[3 * (size_t)b + k] - m->x[3 * (size_t)a + k];
			t += (x[k] - m->x[3 * (size_t)a + k]) * ab;
			l2 += ab * ab;
		}
		t = l2 > 0 ? t / l2 : 1;
		t = t < 0 ? 0 : t > 1 ? 1 : t;
		if (!hm_collapse(m, hm_find_edge(m, a, b), a, t, x))
		{
			failed = 1;   /* (out of memory: what was done stays, compacted) */
			break;
		}
		for (k = 0; k < 10; k++)
			w.Q[10 * (size_t)b + k] += w.Q[10 * (size_t)a + k];
		w.W[b] += w.W[a];
		w.border[b] |= w.border[a];
		w.stamp[a]++;
		w.stamp[b]++;
		live -= gone;
		report->collapses++;
		if (sqrt(cost) > (double)report->max_error)
			report->max_error = (qaws_scalar)sqrt(cost);
		if (!dc_push_edges(&w, &H, b, &ring, &rcap))
		{
			failed = 1;
			break;
		}
	}
	hm_compact(m);
	report->faces_after = dc_faces(m);
	st = failed ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (ref)
		st = qaws_mesh_3d_measure_deviation(&report->deviation, ref, mesh, NULL);
done:
	qaws_mesh_reference_destroy(ref);
	qaws_internal_free(w.Q); qaws_internal_free(w.W); qaws_internal_free(w.stamp); qaws_internal_free(w.border);
	qaws_internal_free(ring); qaws_internal_free(H.e);
	return st;
}
