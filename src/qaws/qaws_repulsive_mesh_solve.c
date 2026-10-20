/*
 * Fractional Sobolev descent of the meshes' tangent-point energy
 * (qaws_repulsive.h, R3; Repulsive Surfaces, Sec. 3.2, 3.3 and 5.2), on
 * the core of qaws_repulsive_solve.c.
 *
 * The metric A = B + B0 is assembled densely over every ordered pair of
 * distinct faces. The preconditioner of the paper's iterative solve,
 * (-L)^+ L^(2 - s) (-L)^+ with L the cotan Laplacian and L^(2 - s) the
 * fractional operator of that order, is built here densely with a GMRES
 * using it, to validate it against the direct solve; the hierarchical
 * mode (qaws_repulsive_hsolve.c) applies it at scale.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include "internal/qaws_internal_broadphase.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define RM_CTOL 1e-5
#else
#define RM_CTOL 1e-8
#endif

typedef qaws_internal_rs_core rs_core;

struct qaws_repulsive_mesh_solver
{
	rs_core core;                   /* (first: the solver is its core) */
	qaws_repulsive_mesh const* m;
	unsigned int *pin, npin, *edges, ne;
	unsigned char* fixed;
	double *pin_target, bary0[3], A0, V0;
	double area_w;                  /* > 0: the area weighted into the energy */
	double vol_w, vol_t;            /* > 0: the volume's deviation from vol_t weighted into it */
	int bary, area, volume;
	/* translations added back (desc's translations 1 or 2): each vertex's
	   group (~0u: none, a component with a pin), the groups' sizes, 3 rows
	   each (their barycenters held while solving, never corrected); keep:
	   the added moves keep the overall vertex mean */
	unsigned int *grp, *gcnt, ngrp;
	int keep;
	qaws_internal_obstacle** obs;   /* nobs fixed obstacles' trees */
	unsigned int nobs;
	/* the obstacles for the collision check: their vertices and points (3
	   each), faces and edges numbered after the mesh's vertices */
	double* ox;
	unsigned int onv, *otri, onf, *oedge, one;
};

/* the barycenter-type rows: 3 per group, or 3 for the held barycenter */
static unsigned int rm_brows(qaws_repulsive_mesh_solver const* s)
{
	return s->grp ? 3 * s->ngrp : s->bary ? 3 : 0;
}

/* ------------------------------------------------------------------------ */
/*  Faces                                                                   */
/* ------------------------------------------------------------------------ */

typedef struct rm_face
{
	double a, N[3], X[3], g[3][3];   /* area, unit normal, barycenter, its corners' hat gradients */
} rm_face;

static void rm_cross(double const* a, double const* b, double* r)
{
	r[0] = a[1] * b[2] - a[2] * b[1]; r[1] = a[2] * b[0] - a[0] * b[2]; r[2] = a[0] * b[1] - a[1] * b[0];
}

static double rm_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* the faces at x; 0 when one has no area */
static int rm_faces(qaws_repulsive_mesh const* m, double const* x, rm_face* f)
{
	unsigned int s, c, q;
	for (s = 0; s < m->nf; s++)
	{
		unsigned int const* t = m->tris + 3 * s;
		double u[3], v[3], n[3], l;
		for (q = 0; q < 3; q++) { u[q] = x[3 * t[1] + q] - x[3 * t[0] + q]; v[q] = x[3 * t[2] + q] - x[3 * t[0] + q]; }
		rm_cross(u, v, n);
		l = sqrt(rm_dot(n, n));
		if (!(l > 0))
			return 0;
		f[s].a = 0.5 * l;
		for (q = 0; q < 3; q++)
		{
			f[s].N[q] = n[q] / l;
			f[s].X[q] = (x[3 * t[0] + q] + x[3 * t[1] + q] + x[3 * t[2] + q]) / 3;
		}
		/* the hat gradient of corner c: N x e / (2a), e the opposite edge in turn */
		for (c = 0; c < 3; c++)
		{
			double e[3];
			for (q = 0; q < 3; q++) e[q] = x[3 * t[(c + 2) % 3] + q] - x[3 * t[(c + 1) % 3] + q];
			rm_cross(f[s].N, e, f[s].g[c]);
			for (q = 0; q < 3; q++) f[s].g[c][q] /= 2 * f[s].a;
		}
	}
	return 1;
}

/* the pair (S, T)'s vertices: each with G (its hat gradient on S less on T)
   and h (1/3 if in S, less 1/3 if in T); merged where the faces share one */
static unsigned int rm_pair_vertices(qaws_repulsive_mesh const* m, rm_face const* f, unsigned int S, unsigned int T,
	unsigned int* id, double (*G)[3], double* h)
{
	unsigned int n = 0, side, c, k, q;
	for (side = 0; side < 2; side++)
	{
		unsigned int F = side ? T : S;
		double sg = side ? -1 : 1;
		for (c = 0; c < 3; c++)
		{
			unsigned int v = m->tris[3 * F + c];
			for (k = 0; k < n && id[k] != v; k++)
				;
			if (k == n) { id[n] = v; G[n][0] = G[n][1] = G[n][2] = 0; h[n] = 0; n++; }
			for (q = 0; q < 3; q++) G[k][q] += sg * f[F].g[c][q];
			h[k] += sg / 3;
		}
	}
	return n;
}

/* A = B + B0 (Eqs. 12, 13) at x, nv x nv */
static qaws_status rm_metric(qaws_repulsive_mesh const* m, double const* x, double* A)
{
	double s = 2 - 2 / m->p, sigma = s - 1, ex = 2 * sigma + 2;
	unsigned int S, T, a, b, nv = m->nv;
	rm_face* f = (rm_face*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(rm_face), (m->nf + 1)));
	if (!f)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (!rm_faces(m, x, f))
	{
		qaws_internal_free(f);
		return QAWS_STATUS_DEGENERATE_CURVE;
	}
	memset(A, 0, sizeof(double) * nv * nv);
	for (S = 0; S < m->nf; S++)
		for (T = S + 1; T < m->nf; T++)
		{
			double d[3], r2, w, K, G[6][3], h[6];
			unsigned int id[6], n, q;
			for (q = 0; q < 3; q++) d[q] = f[S].X[q] - f[T].X[q];
			r2 = rm_dot(d, d);
			if (!(r2 > 0))
			{
				qaws_internal_free(f);
				return QAWS_STATUS_NUMERICAL_FAILURE;
			}
			w = f[S].a * f[T].a * pow(r2, -0.5 * ex);
			K = (rm_dot(f[S].N, d) * rm_dot(f[S].N, d) + rm_dot(f[T].N, d) * rm_dot(f[T].N, d)) / (r2 * r2);
			n = rm_pair_vertices(m, f, S, T, id, G, h);
			for (a = 0; a < n; a++)
				for (b = 0; b < n; b++)
					A[id[a] * nv + id[b]] += 2 * w * rm_dot(G[a], G[b]) + w * K * h[a] * h[b];
		}
	qaws_internal_free(f);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Constraints                                                             */
/* ------------------------------------------------------------------------ */

static void rm_measures(qaws_repulsive_mesh const* m, double const* x, double* area, double* gA, double* volume, double* gV);

static qaws_status rm_energy(rs_core const* core, double const* x, double* E, double* grad)
{
	qaws_repulsive_mesh_solver const* s = (qaws_repulsive_mesh_solver const*)core;
	qaws_repulsive_mesh const* m = s->m;
	qaws_status st = qaws_internal_repulsive_mesh_bh(E, grad, m, x, core->theta);   /* (theta 0 opens every node: the direct terms, in closed form) */
	unsigned int o;
	for (o = 0; st == QAWS_STATUS_OK && o < s->nobs; o++) st = qaws_internal_obstacle_energy(E, grad, s->obs[o], m, x, core->theta);
	if (st == QAWS_STATUS_OK && (s->area_w > 0 || s->vol_w > 0))
	{
		/* the area weighted (the paper's minimize area); the volume's deviation
		   (its program's volume_deviation: w ((V - Vt) / Vt)^2) */
		double area, volume, dev, *gA = grad ? (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, m->nv)) : NULL;
		double* gV = gA ? gA + 3 * m->nv : NULL;
		unsigned int i;
		if (grad && !gA)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		rm_measures(m, x, &area, gA, &volume, gV);
		dev = s->vol_w > 0 ? (volume - s->vol_t) / s->vol_t : 0;
		*E += s->area_w * area + s->vol_w * dev * dev;
		for (i = 0; gA && i < 3 * m->nv; i++) grad[i] += s->area_w * gA[i] + 2 * s->vol_w * dev / s->vol_t * gV[i];
		qaws_internal_free(gA);
	}
	return st;
}

static qaws_status rm_metric_op(rs_core const* core, double const* x, double* A)
{
	return rm_metric(((qaws_repulsive_mesh_solver const*)core)->m, x, A);
}

/* the total area and enclosed volume at x, with their gradients (N each, may be NULL) */
static void rm_measures(qaws_repulsive_mesh const* m, double const* x, double* area, double* gA, double* volume, double* gV)
{
	unsigned int s, c, q;
	*area = *volume = 0;
	if (gA) memset(gA, 0, sizeof(double) * 3 * m->nv);
	if (gV) memset(gV, 0, sizeof(double) * 3 * m->nv);
	for (s = 0; s < m->nf; s++)
	{
		unsigned int const* t = m->tris + 3 * s;
		double u[3], v[3], n[3], l, jk[3];
		for (q = 0; q < 3; q++) { u[q] = x[3 * t[1] + q] - x[3 * t[0] + q]; v[q] = x[3 * t[2] + q] - x[3 * t[0] + q]; }
		rm_cross(u, v, n);
		l = sqrt(rm_dot(n, n));
		*area += 0.5 * l;
		rm_cross(x + 3 * t[1], x + 3 * t[2], jk);
		*volume += rm_dot(x + 3 * t[0], jk) / 6;
		if (l > 0)
			for (q = 0; q < 3; q++) n[q] /= l;
		for (c = 0; c < 3; c++)
		{
			unsigned int i = t[c], j = t[(c + 1) % 3], k = t[(c + 2) % 3];
			double e[3], r[3];
			if (gA && l > 0)
			{
				/* d area / d x_i = N x (x_k - x_j) / 2 */
				for (q = 0; q < 3; q++) e[q] = x[3 * k + q] - x[3 * j + q];
				rm_cross(n, e, r);
				for (q = 0; q < 3; q++) gA[3 * i + q] += 0.5 * r[q];
			}
			if (gV)
			{
				/* d volume / d x_i = x_j x x_k / 6 */
				rm_cross(x + 3 * j, x + 3 * k, r);
				for (q = 0; q < 3; q++) gV[3 * i + q] += r[q] / 6;
			}
		}
	}
}

/* rows: barycenter, pins, area / l, volume / l^2 (all in lengths) */
static void rm_constraints(rs_core const* core, double const* x, double* phi, double* C)
{
	qaws_repulsive_mesh_solver const* s = (qaws_repulsive_mesh_solver const*)core;
	unsigned int N = core->N, r = 0, i, q;
	double l = core->lmean, area, volume;
	if (C) memset(C, 0, sizeof(double) * core->k * N);
	if (s->grp)
	{
		/* (residuals 0: the groups move freely; only the direction keeps their means) */
		for (r = 0; r < 3 * s->ngrp; r++) phi[r] = 0;
		for (i = 0; C && i < core->nv; i++)
			for (q = 0; q < 3 && s->grp[i] != ~0u; q++) C[(3 * s->grp[i] + q) * N + 3 * i + q] = 1.0 / s->gcnt[s->grp[i]];
	}
	else if (s->bary)
		for (q = 0; q < 3; q++, r++)
		{
			double mm = 0;
			for (i = 0; i < core->nv; i++)
			{
				mm += x[3 * i + q] / core->nv;
				if (C) C[r * N + 3 * i + q] = 1.0 / core->nv;
			}
			phi[r] = mm - s->bary0[q];
		}
	for (i = 0; core->theta <= 0 && i < s->npin; i++)   /* (hierarchical: eliminated) */
		for (q = 0; q < 3; q++, r++)
		{
			phi[r] = x[3 * s->pin[i] + q] - s->pin_target[3 * i + q];
			if (C) C[r * N + 3 * s->pin[i] + q] = 1;
		}
	if (s->area || s->volume)
	{
		double* gA = C && s->area ? C + (size_t)r * N : NULL;
		double* gV = C && s->volume ? C + (size_t)(r + (s->area ? 1 : 0)) * N : NULL;
		rm_measures(s->m, x, &area, gA, &volume, gV);
		if (s->area)
		{
			phi[r] = (area - s->A0) / l;
			for (i = 0; gA && i < N; i++) gA[i] /= l;
			r++;
		}
		if (s->volume)
		{
			phi[r] = (volume - s->V0) / (l * l);
			for (i = 0; gV && i < N; i++) gV[i] /= l * l;
		}
	}
}

/* ------------------------------------------------------------------------ */
/*  Collision checks                                                        */
/* ------------------------------------------------------------------------ */

/* the distance from p to triangle a b c (Ericson's closest point) */
static double rm_pointtri(double const* p, double const* a, double const* b, double const* c)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, den, r[3];
	unsigned int q;
	for (q = 0; q < 3; q++) { ab[q] = b[q] - a[q]; ac[q] = c[q] - a[q]; ap[q] = p[q] - a[q]; bp[q] = p[q] - b[q]; cp[q] = p[q] - c[q]; }
	d1 = rm_dot(ab, ap); d2 = rm_dot(ac, ap);
	if (d1 <= 0 && d2 <= 0) { for (q = 0; q < 3; q++) r[q] = a[q]; goto out; }
	d3 = rm_dot(ab, bp); d4 = rm_dot(ac, bp);
	if (d3 >= 0 && d4 <= d3) { for (q = 0; q < 3; q++) r[q] = b[q]; goto out; }
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) { v = d1 / (d1 - d3); for (q = 0; q < 3; q++) r[q] = a[q] + v * ab[q]; goto out; }
	d5 = rm_dot(ab, cp); d6 = rm_dot(ac, cp);
	if (d6 >= 0 && d5 <= d6) { for (q = 0; q < 3; q++) r[q] = c[q]; goto out; }
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) { w = d2 / (d2 - d6); for (q = 0; q < 3; q++) r[q] = a[q] + w * ac[q]; goto out; }
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
	{
		w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		for (q = 0; q < 3; q++) r[q] = b[q] + w * (c[q] - b[q]);
		goto out;
	}
	den = 1 / (va + vb + vc);
	v = vb * den; w = vc * den;
	for (q = 0; q < 3; q++) r[q] = a[q] + ab[q] * v + ac[q] * w;
out:
	for (q = 0; q < 3; q++) r[q] -= p[q];
	return sqrt(rm_dot(r, r));
}

/* whether segment p q meets triangle a b c (Moller-Trumbore, ends included;
   a near-coplanar meeting counts too) */
static int rm_segtri(double const* p, double const* q, double const* a, double const* b, double const* c)
{
	double d[3], e1[3], e2[3], h[3], s[3], k[3], det, u, v, t, scale;
	unsigned int i;
	for (i = 0; i < 3; i++) { d[i] = q[i] - p[i]; e1[i] = b[i] - a[i]; e2[i] = c[i] - a[i]; s[i] = p[i] - a[i]; }
	rm_cross(d, e2, h);
	det = rm_dot(e1, h);
	scale = sqrt(rm_dot(d, d) * rm_dot(e1, e1) * rm_dot(e2, e2));
	if (fabs(det) <= 1e-14 * scale)
	{
		/* parallel: a meeting only if the segment lies in the plane and on
		   the triangle; the distances settle it */
		double n[3];
		rm_cross(e1, e2, n);
		return fabs(rm_dot(n, s)) <= 1e-14 * sqrt(rm_dot(n, n)) * sqrt(rm_dot(s, s) + 1e-300) &&
			(rm_pointtri(p, a, b, c) <= 0 || rm_pointtri(q, a, b, c) <= 0);
	}
	u = rm_dot(s, h) / det;
	if (u < 0 || u > 1)
		return 0;
	rm_cross(s, e1, k);
	v = rm_dot(d, k) / det;
	if (v < 0 || u + v > 1)
		return 0;
	t = rm_dot(e2, k) / det;
	return t >= 0 && t <= 1;
}

typedef struct rm_pair
{
	unsigned int v[4];   /* vertex-face: v[0] the vertex, v[1..3] the face; edge-edge: two edges */
	int vf;
	double const *x0, *x1;
} rm_pair;

static double rm_pair_at(void const* user, double t)
{
	rm_pair const* p = (rm_pair const*)user;
	double y[4][3];
	unsigned int k, q;
	for (k = 0; k < 4; k++)
		for (q = 0; q < 3; q++) y[k][q] = (1 - t) * p->x0[3 * p->v[k] + q] + t * p->x1[3 * p->v[k] + q];
	return p->vf ? rm_pointtri(y[0], y[1], y[2], y[3]) : qaws_internal_rs_segdist(y[0], y[1], y[2], y[3], 3);
}

static double rm_move(double const* x0, double const* x1, unsigned int v)
{
	double d[3];
	unsigned int q;
	for (q = 0; q < 3; q++) d[q] = x1[3 * v + q] - x0[3 * v + q];
	return sqrt(rm_dot(d, d));
}

/* every vertex against every face it is not a corner of, every two edges
   sharing no vertex: at least min(clearance, half the start distance);
   check = 0: whether no such pair touches at x0 (the input check) */
/* the vertex-face and edge-edge pairs (with check 0: the still input, and
   edge-face crossings) through the grid broad phase over one box list:
   the mesh's vertices [0, nv), edges [nv, nv + ne), faces; then the
   obstacles' vertices and points, edges, faces. Each box holds its
   element's straight motion x0 -> x1 (an obstacle's still), grown by pad;
   boxes farther apart than 2 pad hold elements that stay farther apart
   than the clearance. Positions: the mesh's, then the obstacles' (fixed),
   so one vertex number reads either. Obstacle pairs: a mesh vertex and an
   obstacle face, a mesh face and an obstacle vertex or point, two edges. */
typedef struct rm_bp
{
	qaws_repulsive_mesh_solver const* s;
	double const *x0, *x1;
	int check;
} rm_bp;

/* an element's kind (0 vertex, 1 edge, 2 face; 3, 4, 5 the obstacles'),
   and its vertices (numbers into the joint positions) */
static int rm_bp_kind(rm_bp const* b, unsigned int i, unsigned int const** vs)
{
	qaws_repulsive_mesh_solver const* s = b->s;
	unsigned int n0 = s->m->nv, n1 = n0 + s->ne, n2 = n1 + s->m->nf, n3 = n2 + s->onv, n4 = n3 + s->one;
	if (i < n0) { *vs = NULL; return 0; }
	if (i < n1) { *vs = s->edges + 2 * (i - n0); return 1; }
	if (i < n2) { *vs = s->m->tris + 3 * (i - n1); return 2; }
	if (i < n3) { *vs = NULL; return 3; }
	if (i < n4) { *vs = s->oedge + 2 * (i - n3); return 4; }
	*vs = s->otri + 3 * (i - n4);
	return 5;
}

/* a vertex element's number into the joint positions */
static unsigned int rm_bp_vertex(rm_bp const* b, unsigned int i, int kind)
{
	return kind == 0 ? i : b->s->m->nv + (i - b->s->m->nv - b->s->ne - b->s->m->nf);
}

/* the pair's two elements, the lower kind first */
static void rm_bp_order(rm_bp const* b, unsigned int* i, unsigned int* j, int* ki, int* kj, unsigned int const** vi, unsigned int const** vj)
{
	*ki = rm_bp_kind(b, *i, vi); *kj = rm_bp_kind(b, *j, vj);
	if (*ki > *kj)
	{
		unsigned int t = *i;
		unsigned int const* v = *vi;
		int k = *ki;
		*i = *j; *j = t; *vi = *vj; *vj = v; *ki = *kj; *kj = k;
	}
}

static int rm_bp_accept(void* user, unsigned int i, unsigned int j)
{
	rm_bp const* b = (rm_bp const*)user;
	unsigned int const *vi, *vj;
	int ki, kj;
	rm_bp_order(b, &i, &j, &ki, &kj, &vi, &vj);
	if (ki == 0 && kj == 2)
		return vj[0] != i && vj[1] != i && vj[2] != i;
	if (ki == 1 && (kj == 1 || (kj == 2 && !b->check)))
	{
		unsigned int n = kj == 1 ? 2 : 3, k;
		for (k = 0; k < n; k++)
			if (vj[k] == vi[0] || vj[k] == vi[1])
				return 0;
		return 1;
	}
	/* (an obstacle's own pairs are never checked: it does not move) */
	return (ki == 0 && kj == 5) || (ki == 2 && kj == 3) || (ki == 1 && kj == 4);
}

static qaws_status rm_bp_visit(void* user, unsigned int i, unsigned int j)
{
	rm_bp const* b = (rm_bp const*)user;
	qaws_repulsive_mesh_solver const* s = b->s;
	unsigned int const *vi, *vj;
	int ki, kj;
	double const *x0 = b->x0, *x1 = b->x1;
	double d0, L;
	rm_pair p;
	rm_bp_order(b, &i, &j, &ki, &kj, &vi, &vj);
	p.x0 = x0; p.x1 = x1;
	if (ki == 1 && kj == 2)
	{
		/* the input: an edge through a face (motion from a clear state cannot reach this without a vertex-face or edge-edge contact first) */
		return rm_segtri(x0 + 3 * vi[0], x0 + 3 * vi[1], x0 + 3 * vj[0], x0 + 3 * vj[1], x0 + 3 * vj[2]) ? QAWS_STATUS_NUMERICAL_FAILURE : QAWS_STATUS_OK;
	}
	if (kj == 2 || kj == 5 || (ki == 2 && kj == 3))
	{
		/* a vertex and a face (either one the obstacle's) */
		unsigned int v = ki == 2 ? rm_bp_vertex(b, j, kj) : rm_bp_vertex(b, i, ki);
		unsigned int const* t = ki == 2 ? vi : vj;
		p.vf = 1; p.v[0] = v; p.v[1] = t[0]; p.v[2] = t[1]; p.v[3] = t[2];
		L = rm_move(x0, x1, v) + fmax(rm_move(x0, x1, t[0]), fmax(rm_move(x0, x1, t[1]), rm_move(x0, x1, t[2])));
	}
	else
	{
		p.vf = 0; p.v[0] = vi[0]; p.v[1] = vi[1]; p.v[2] = vj[0]; p.v[3] = vj[1];
		L = fmax(rm_move(x0, x1, vi[0]), rm_move(x0, x1, vi[1])) + fmax(rm_move(x0, x1, vj[0]), rm_move(x0, x1, vj[1]));
	}
	d0 = rm_pair_at(&p, 0);
	if (!b->check)
		return d0 > 1e-12 * s->core.lmean ? QAWS_STATUS_OK : QAWS_STATUS_NUMERICAL_FAILURE;
	return qaws_internal_rs_advance(rm_pair_at, &p, d0, L, fmin(s->core.clearance, 0.5 * d0)) ? QAWS_STATUS_OK : QAWS_STATUS_NUMERICAL_FAILURE;
}

static int rm_pairs_ok(qaws_repulsive_mesh_solver const* s, double const* x0, double const* x1, int check)
{
	qaws_repulsive_mesh const* m = s->m;
	unsigned int nv = m->nv, ne = s->ne, n = nv + ne + m->nf + s->onv + s->one + s->onf, i, a, q;
	double pad = check ? 0.5 * s->core.clearance * (1 + 1e-9) : 1e-12 * s->core.lmean;
	qaws_bp_box* box = (qaws_bp_box*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), (n + 1)));
	double *y0 = NULL, *y1 = NULL;
	rm_bp b;
	qaws_status st;
	if (!box)
		return 0;
	if (s->onv)
	{
		/* the joint positions: the mesh's, then the obstacles' */
		y0 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (nv + s->onv)));
		y1 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (nv + s->onv)));
		if (!y0 || !y1)
		{
			qaws_internal_free(box); qaws_internal_free(y0); qaws_internal_free(y1);
			return 0;
		}
		memcpy(y0, x0, sizeof(double) * 3 * nv); memcpy(y0 + 3 * nv, s->ox, sizeof(double) * 3 * s->onv);
		memcpy(y1, x1, sizeof(double) * 3 * nv); memcpy(y1 + 3 * nv, s->ox, sizeof(double) * 3 * s->onv);
		x0 = y0; x1 = y1;
	}
	b.s = s; b.x0 = x0; b.x1 = x1; b.check = check;
	for (i = 0; i < n; i++)
	{
		unsigned int const* vs;
		unsigned int self[1];
		int kind = rm_bp_kind(&b, i, &vs);
		unsigned int k = kind % 3 == 0 ? 1 : kind % 3 == 1 ? 2 : 3;
		if (!vs) { self[0] = rm_bp_vertex(&b, i, kind); vs = self; }
		for (q = 0; q < 3; q++) { box[i].lo[q] = HUGE_VAL; box[i].hi[q] = -HUGE_VAL; }
		for (a = 0; a < k; a++)
			for (q = 0; q < 3; q++)
			{
				double u = x0[3 * vs[a] + q], v = x1[3 * vs[a] + q];
				box[i].lo[q] = fmin(box[i].lo[q], fmin(u, v) - pad);
				box[i].hi[q] = fmax(box[i].hi[q], fmax(u, v) + pad);
			}
	}
	st = qaws_internal_rbh_box_pairs(box, n, 3, rm_bp_accept, rm_bp_visit, &b);
	qaws_internal_free(box); qaws_internal_free(y0); qaws_internal_free(y1);
	return st == QAWS_STATUS_OK;
}

static int rm_path_ok(rs_core const* core, double const* x0, double const* x1)
{
	return rm_pairs_ok((qaws_repulsive_mesh_solver const*)core, x0, x1, 1);
}

/* the hierarchical operators: the metric; -L by the cotan weights (each
   face's corner weighs its opposite edge), masses a/3 per corner;
   L^(2 - s) by vertex quadrature, r^-(2 + 2 (2 - s)) */
static qaws_status rm_hier(rs_core const* core, double const* x, qaws_internal_rs_hier* h)
{
	qaws_repulsive_mesh const* m = ((qaws_repulsive_mesh_solver const*)core)->m;
	unsigned int f, c, q;
	double s = 2 - 2 / m->p;
	qaws_status st;
	memset(h, 0, sizeof(*h));
	h->lij = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, m->nf));
	h->lw = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, m->nf));
	h->mass = (double*)qaws_internal_zalloc(m->nv + 1, sizeof(double));
	if (!h->lij || !h->lw || !h->mass)
	{
		qaws_internal_rs_hier_free(h);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (f = 0; f < m->nf; f++)
		for (c = 0; c < 3; c++)
		{
			unsigned int o = m->tris[3 * f + c], i = m->tris[3 * f + (c + 1) % 3], j = m->tris[3 * f + (c + 2) % 3];
			double u[3], v[3], n[3], l;
			for (q = 0; q < 3; q++) { u[q] = x[3 * i + q] - x[3 * o + q]; v[q] = x[3 * j + q] - x[3 * o + q]; }
			rm_cross(u, v, n);
			l = sqrt(rm_dot(n, n));
			h->lij[2 * h->nl] = i; h->lij[2 * h->nl + 1] = j;
			h->lw[h->nl++] = l > 0 ? 0.5 * rm_dot(u, v) / l : 0;
			h->mass[o] += l / 6;
		}
	if ((st = qaws_internal_rbh_mesh_metric(&h->A, m, x, core->theta)) == QAWS_STATUS_OK)
		st = qaws_internal_rbh_fractional(&h->F, x, h->mass, m->nv, 3, 2 + 2 * (2 - s), core->theta, core->executor);
	if (st != QAWS_STATUS_OK)
		qaws_internal_rs_hier_free(h);
	return st;
}

static qaws_status rm_hvp(rs_core const* core, double const* x, double const* v, double* hv)
{
	return qaws_internal_repulsive_mesh_hvp(hv, ((qaws_repulsive_mesh_solver const*)core)->m, x, v);
}

/* hv += sum_r lambda_r Hess(c_r) v: the area row (over l) and the volume
   row (over l^2) are not linear. Per face, m = e1 x e2, a = |m| / 2:
   grad a . w = m . m' / (2 |m|) with m' = w1' x e2 + e1 x w2'; its gradient
   is g1 = e2 x u + w2' x z (on e1), g2 = u x e1 + z x w1' (on e2), with
   u = (m' - (m . m') m / |m|^2) / (2 |m|), z = m / (2 |m|). The volume's
   (1/6) x0 . (x1 x x2) gives (w1 x x2 + x1 x w2) / 6 at x0, cyclically. */
static void rm_constraints_hvp(rs_core const* core, double const* x, double const* lambda, double const* v, double* hv)
{
	qaws_repulsive_mesh_solver const* s = (qaws_repulsive_mesh_solver const*)core;
	unsigned int r = rm_brows(s) + (core->theta > 0 ? 0 : 3 * s->npin), f, c, q;
	double la = s->area ? lambda[r] / core->lmean : 0, lv = s->volume ? lambda[r + (s->area ? 1 : 0)] / (core->lmean * core->lmean) : 0;
	if (!s->area && !s->volume)
		return;
	for (f = 0; f < s->m->nf; f++)
	{
		unsigned int const* t = s->m->tris + 3 * f;
		double e1[3], e2[3], w1[3], w2[3], m[3], md[3], a1[3], a2[3], u[3], z[3], g1[3], g2[3], mm, mmd;
		if (la != 0)
		{
			for (q = 0; q < 3; q++)
			{
				e1[q] = x[3 * t[1] + q] - x[3 * t[0] + q]; e2[q] = x[3 * t[2] + q] - x[3 * t[0] + q];
				w1[q] = v[3 * t[1] + q] - v[3 * t[0] + q]; w2[q] = v[3 * t[2] + q] - v[3 * t[0] + q];
			}
			rm_cross(e1, e2, m);
			rm_cross(w1, e2, a1); rm_cross(e1, w2, a2);
			for (q = 0; q < 3; q++) md[q] = a1[q] + a2[q];
			mm = sqrt(rm_dot(m, m));
			if (mm > 0)
			{
				mmd = rm_dot(m, md);
				for (q = 0; q < 3; q++) { u[q] = (md[q] - mmd * m[q] / (mm * mm)) / (2 * mm); z[q] = m[q] / (2 * mm); }
				rm_cross(e2, u, g1); rm_cross(w2, z, a1);
				rm_cross(u, e1, g2); rm_cross(z, w1, a2);
				for (q = 0; q < 3; q++)
				{
					g1[q] += a1[q]; g2[q] += a2[q];
					hv[3 * t[1] + q] += la * g1[q]; hv[3 * t[2] + q] += la * g2[q]; hv[3 * t[0] + q] -= la * (g1[q] + g2[q]);
				}
			}
		}
		if (lv != 0)
			for (c = 0; c < 3; c++)
			{
				unsigned int i = t[c], j = t[(c + 1) % 3], k = t[(c + 2) % 3];
				rm_cross(v + 3 * j, x + 3 * k, a1); rm_cross(x + 3 * j, v + 3 * k, a2);
				for (q = 0; q < 3; q++) hv[3 * i + q] += lv * (a1[q] + a2[q]) / 6;
			}
	}
}

/* g += each group's dual-area mean of dE (the groups' rows kept g's
   vertex mean 0 in each); with keep, less their mean over the grouped
   vertices, so the overall vertex mean stays */
static void rm_translate(rs_core const* core, double const* x, double const* dE, double* g)
{
	qaws_repulsive_mesh_solver const* s = (qaws_repulsive_mesh_solver const*)core;
	unsigned int nv = core->nv, i, f, q, c, n = 0;
	double* w = (double*)qaws_internal_scratch_zalloc(nv + 4 * (size_t)s->ngrp, sizeof(double));
	double* t;
	double keep[3] = { 0, 0, 0 };
	if (!w)
		return;   /* shortcut: no memory, no translations added this step; the direction stays a descent direction */
	t = w + nv;
	for (f = 0; f < s->m->nf; f++)
	{
		unsigned int const* tri = s->m->tris + 3 * f;
		double e1[3], e2[3], m[3], a;
		for (q = 0; q < 3; q++) { e1[q] = x[3 * tri[1] + q] - x[3 * tri[0] + q]; e2[q] = x[3 * tri[2] + q] - x[3 * tri[0] + q]; }
		rm_cross(e1, e2, m);
		a = sqrt(rm_dot(m, m)) / 6;   /* (a third of the face's area) */
		for (c = 0; c < 3; c++) w[tri[c]] += a;
	}
	for (i = 0; i < nv; i++)
		if (s->grp[i] != ~0u)
		{
			double* tg = t + 4 * s->grp[i];
			for (q = 0; q < 3; q++) tg[q] += w[i] * dE[3 * i + q];
			tg[3] += w[i];
		}
	for (c = 0; c < s->ngrp; c++)
		for (q = 0; q < 3; q++)
		{
			t[4 * c + q] = t[4 * c + 3] > 0 ? t[4 * c + q] / t[4 * c + 3] : 0;
			keep[q] += s->gcnt[c] * t[4 * c + q];
		}
	for (c = 0; c < s->ngrp; c++) n += s->gcnt[c];
	for (i = 0; i < nv; i++)
		for (q = 0; q < 3 && s->grp[i] != ~0u; q++) g[3 * i + q] += t[4 * s->grp[i] + q] - (s->keep ? keep[q] / n : 0);
	qaws_internal_free(w);
}

static qaws_internal_rs_ops const rm_ops = { rm_energy, rm_metric_op, rm_constraints, rm_path_ok, rm_hier, NULL, rm_hvp, rm_constraints_hvp, NULL };
/* with the area weighted: no Hessian products, so no solution derivatives */
static qaws_internal_rs_ops const rm_ops_plain = { rm_energy, rm_metric_op, rm_constraints, rm_path_ok, rm_hier, NULL, NULL, NULL, NULL };
/* with translations added back: no solution derivatives either */
static qaws_internal_rs_ops const rm_ops_translate = { rm_energy, rm_metric_op, rm_constraints, rm_path_ok, rm_hier, NULL, NULL, NULL, rm_translate };

/* the rows' targets from t (to rows, dir 0) or back (dir 1) */
static void rm_targets(qaws_repulsive_mesh_solver const* s, qaws_repulsive_targets* t, double* rows, int dir)
{
	unsigned int r = 0, i, q;
	double l = s->core.lmean;
	if (dir)
	{
		memset(t->barycenter, 0, sizeof(t->barycenter));
		t->length = t->area = t->volume = 0;
	}
	if (s->bary)
		for (q = 0; q < 3; q++, r++)
		{
			if (dir) t->barycenter[q] = (qaws_scalar)rows[r]; else rows[r] = (double)t->barycenter[q];
		}
	for (i = 0; s->core.theta <= 0 && i < s->npin; i++)   /* (hierarchical: pins are not rows) */
		for (q = 0; q < 3; q++, r++)
		{
			if (dir) { if (t->pinned) t->pinned[3 * i + q] = (qaws_scalar)rows[r]; }
			else rows[r] = t->pinned ? (double)t->pinned[3 * i + q] : 0;
		}
	/* (the rows hold (A - A0) / l and (V - V0) / l^2) */
	if (s->area)
	{
		if (dir) t->area = (qaws_scalar)(rows[r] / l); else rows[r] = (double)t->area / l;
		r++;
	}
	if (s->volume)
	{
		if (dir) t->volume = (qaws_scalar)(rows[r] / (l * l)); else rows[r] = (double)t->volume / (l * l);
	}
}

qaws_status qaws_repulsive_mesh_solver_solution_jvp(qaws_scalar* out_dx, unsigned int* out_gauge, qaws_repulsive_mesh_solver const* s,
	qaws_repulsive_targets const* d)
{
	return qaws_repulsive_mesh_solver_solution_jvp_checked(out_dx, out_gauge, NULL, s, d);
}

qaws_status qaws_repulsive_mesh_solver_solution_jvp_checked(qaws_scalar* out_dx, unsigned int* out_gauge, qaws_repulsive_derivative_report* out_report,
	qaws_repulsive_mesh_solver const* s, qaws_repulsive_targets const* d)
{
	double *rows, *dx;
	unsigned int i;
	qaws_status st;
	qaws_internal_rs_scope sc;
	if (out_report) memset(out_report, 0, sizeof(*out_report));
	if (!s || !d || !out_dx)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_rs_enter(&sc, &s->core.alloc);
	rows = (double*)qaws_internal_scratch_zalloc(s->core.k + 1, sizeof(double));
	dx = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (s->core.N + 1)));
	if (!rows || !dx)
	{
		qaws_internal_free(rows); qaws_internal_free(dx);
		qaws_internal_rs_leave(&sc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	rm_targets(s, (qaws_repulsive_targets*)d, rows, 0);
	st = qaws_internal_rs_sensitivity(dx, out_gauge, out_report, &s->core, rows, NULL);
	for (i = 0; st == QAWS_STATUS_OK && i < s->core.N; i++) out_dx[i] = (qaws_scalar)dx[i];
	qaws_internal_free(rows); qaws_internal_free(dx);
	qaws_internal_rs_leave(&sc);
	return st;
}

qaws_status qaws_repulsive_mesh_solver_solution_vjp(qaws_repulsive_targets* out_bar, unsigned int* out_gauge, qaws_repulsive_mesh_solver const* s,
	qaws_scalar const* xbar)
{
	return qaws_repulsive_mesh_solver_solution_vjp_checked(out_bar, out_gauge, NULL, s, xbar);
}

qaws_status qaws_repulsive_mesh_solver_solution_vjp_checked(qaws_repulsive_targets* out_bar, unsigned int* out_gauge, qaws_repulsive_derivative_report* out_report,
	qaws_repulsive_mesh_solver const* s, qaws_scalar const* xbar)
{
	double *rows, *xb;
	unsigned int i;
	qaws_status st;
	qaws_internal_rs_scope sc;
	if (out_report) memset(out_report, 0, sizeof(*out_report));
	if (!s || !xbar || !out_bar)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_rs_enter(&sc, &s->core.alloc);
	rows = (double*)qaws_internal_scratch_zalloc(s->core.k + 1, sizeof(double));
	xb = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (s->core.N + 1)));
	if (!rows || !xb)
	{
		qaws_internal_free(rows); qaws_internal_free(xb);
		qaws_internal_rs_leave(&sc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < s->core.N; i++) xb[i] = (double)xbar[i];
	st = qaws_internal_rs_sensitivity(rows, out_gauge, out_report, &s->core, NULL, xb);
	if (st == QAWS_STATUS_OK)
		rm_targets(s, out_bar, rows, 1);
	qaws_internal_free(rows); qaws_internal_free(xb);
	qaws_internal_rs_leave(&sc);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  The solver                                                              */
/* ------------------------------------------------------------------------ */

void qaws_repulsive_mesh_solver_destroy(qaws_repulsive_mesh_solver* s)
{
	unsigned int i;
	if (!s)
		return;
	qaws_internal_free(s->pin); qaws_internal_free(s->pin_target); qaws_internal_free(s->fixed); qaws_internal_free(s->edges); qaws_internal_free(s->core.x); qaws_internal_free(s->grp); qaws_internal_free(s->gcnt);
	for (i = 0; i < s->nobs; i++) qaws_internal_obstacle_destroy(s->obs[i]);
	qaws_internal_free(s->ox); qaws_internal_free(s->otri); qaws_internal_free(s->oedge);
	qaws_internal_free(s->obs);
	qaws_internal_free(s);
}

/* the mesh's edges, each once (sorted ends); and whether every edge lies on
   exactly two faces, oppositely oriented (closed) */
static qaws_status rm_edges(qaws_repulsive_mesh const* m, unsigned int** out, unsigned int* out_n, int* closed)
{
	unsigned int* e = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 6, (m->nf + 1)));
	int* cnt = (int*)qaws_internal_scratch_zalloc(3 * (size_t)m->nf + 1, sizeof(int));
	unsigned int n = 0, f, c, k;
	if (!e || !cnt)
	{
		qaws_internal_free(e); qaws_internal_free(cnt);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	*closed = 1;
	for (f = 0; f < m->nf; f++)
		for (c = 0; c < 3; c++)
		{
			unsigned int a = m->tris[3 * f + c], b = m->tris[3 * f + (c + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
			for (k = 0; k < n && !(e[2 * k] == lo && e[2 * k + 1] == hi); k++)
				;
			if (k == n) { e[2 * n] = lo; e[2 * n + 1] = hi; n++; }
			cnt[k] += a < b ? 1 : 16;   /* (the two directions apart) */
		}
	for (k = 0; k < n; k++)
		if (cnt[k] != 17)
			*closed = 0;
	qaws_internal_free(cnt);
	*out = e;
	*out_n = n;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_mesh_solver_create(qaws_repulsive_mesh_solver** out_solver, qaws_repulsive_mesh const* m,
	qaws_scalar const* positions, qaws_repulsive_mesh_solve_desc const* desc)
{
	qaws_repulsive_mesh_solver* s;
	rs_core* k;
	qaws_repulsive_mesh_solve_desc d;
	unsigned int i, q;
	int closed = 0;
	double total = 0, E;
	qaws_status st;
	if (!m || !positions || !out_solver || m->nf < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_solver = NULL;
	memset(&d, 0, sizeof(d));
	if (desc) d = *desc;
	if ((d.pinned_count && !d.pinned) || !(d.area_weight >= 0) || !(d.initial_step >= 0) || !(d.max_step >= 0) || !(d.volume_weight >= 0) || d.translations < 0 || d.translations > 2 ||
		(d.translations == 1 && d.pinned_count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < d.pinned_count; i++) if (d.pinned[i] >= m->nv) return QAWS_STATUS_INVALID_ARGUMENT;
	s = (qaws_repulsive_mesh_solver*)qaws_internal_zalloc(1, sizeof(*s));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	k = &s->core;
	k->ops = &rm_ops;
	s->m = m; k->dim = 3; k->nv = m->nv; k->N = 3 * m->nv;
	s->npin = d.pinned_count;
	s->area = d.hold_area != 0; s->volume = d.hold_volume != 0;
	s->area_w = (double)d.area_weight;
	s->vol_w = (double)d.volume_weight;
	if (s->area_w > 0 || s->vol_w > 0 || d.obstacle_count)
		k->ops = &rm_ops_plain;
	if (d.obstacle_count)
	{
		if (!d.obstacles)
		{
			qaws_repulsive_mesh_solver_destroy(s);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		s->obs = (qaws_internal_obstacle**)qaws_internal_zalloc(d.obstacle_count, sizeof(qaws_internal_obstacle*));
		if (!s->obs)
		{
			qaws_repulsive_mesh_solver_destroy(s);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < d.obstacle_count; i++)
		{
			if ((st = qaws_internal_obstacle_create(s->obs + i, d.obstacles + i)) != QAWS_STATUS_OK)
			{
				qaws_repulsive_mesh_solver_destroy(s);
				return st;
			}
			s->nobs++;
		}
		if (d.collision_checked)
		{
			/* their geometry for the collision check, numbered after the mesh's vertices */
			unsigned int base = 0, fb = 0, o;
			for (o = 0; o < d.obstacle_count; o++) { s->onv += d.obstacles[o].vertex_count; s->onf += d.obstacles[o].triangles ? d.obstacles[o].triangle_count : 0; }
			s->ox = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (s->onv + 1)));
			s->otri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (s->onf + 1)));
			s->oedge = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (3 * s->onf + 1)));
			if (!s->ox || !s->otri || !s->oedge)
			{
				qaws_repulsive_mesh_solver_destroy(s);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			for (o = 0; o < d.obstacle_count; o++)
			{
				qaws_repulsive_obstacle const* ob = d.obstacles + o;
				for (i = 0; i < 3 * ob->vertex_count; i++) s->ox[3 * base + i] = (double)ob->positions[i];
				if (ob->triangles)
				{
					qaws_repulsive_mesh tm;
					unsigned int *e = NULL, ne = 0, j;
					int cl;
					memset(&tm, 0, sizeof(tm));
					tm.nv = ob->vertex_count; tm.nf = ob->triangle_count; tm.tris = (unsigned int*)ob->triangles;   /* (read only) */
					for (j = 0; j < 3 * ob->triangle_count; j++) s->otri[3 * fb + j] = m->nv + base + ob->triangles[j];
					if (rm_edges(&tm, &e, &ne, &cl) != QAWS_STATUS_OK)
					{
						qaws_repulsive_mesh_solver_destroy(s);
						return QAWS_STATUS_ALLOCATION_FAILURE;
					}
					for (j = 0; j < 2 * ne; j++) s->oedge[2 * s->one + j] = m->nv + base + e[j];
					s->one += ne;
					fb += ob->triangle_count;
					qaws_internal_free(e);
				}
				base += ob->vertex_count;
			}
		}
	}
	s->bary = d.hold_barycenter > 0 || (d.hold_barycenter == 0 && d.pinned_count == 0);
	k->checked = d.collision_checked != 0;
	k->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), k->N));
	s->pin = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (s->npin + 1)));
	s->pin_target = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * s->npin + 1)));
	s->fixed = (unsigned char*)qaws_internal_zalloc(m->nv + 1, 1);
	if (!k->x || !s->pin || !s->pin_target || !s->fixed || (st = rm_edges(m, &s->edges, &s->ne, &closed)) != QAWS_STATUS_OK)
	{
		qaws_repulsive_mesh_solver_destroy(s);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if ((s->volume || s->vol_w > 0) && !closed)
	{
		/* (an open mesh encloses no volume) */
		qaws_repulsive_mesh_solver_destroy(s);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	for (i = 0; i < k->N; i++) k->x[i] = (double)positions[i];
	for (q = 0; q < 3; q++)
	{
		s->bary0[q] = 0;
		for (i = 0; i < k->nv; i++) s->bary0[q] += k->x[3 * i + q] / k->nv;
	}
	for (i = 0; i < s->npin; i++)
	{
		s->pin[i] = d.pinned[i];
		for (q = 0; q < 3; q++) s->pin_target[3 * i + q] = k->x[3 * d.pinned[i] + q];
		s->fixed[d.pinned[i]] = 1;
	}
	if (d.translations)
	{
		/* the groups: one, or the connected components (union-find on the
		   edges, each rooted at its least vertex) without a pin */
		unsigned int* root = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv + 1)));
		unsigned char* pinned = (unsigned char*)qaws_internal_scratch_zalloc(m->nv + 1, 1);
		s->grp = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv + 1)));
		s->gcnt = (unsigned int*)qaws_internal_zalloc(m->nv + 1, sizeof(unsigned int));
		if (!root || !pinned || !s->grp || !s->gcnt)
		{
			qaws_internal_free(root); qaws_internal_free(pinned);
			qaws_repulsive_mesh_solver_destroy(s);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < m->nv; i++) root[i] = d.translations == 1 ? 0 : i;
		for (i = 0; d.translations == 2 && i < s->ne; i++)
		{
			unsigned int a = s->edges[2 * i], b = s->edges[2 * i + 1];
			while (root[a] != a) a = root[a];
			while (root[b] != b) b = root[b];
			if (a < b) root[b] = a; else if (b < a) root[a] = b;
		}
		for (i = 0; i < m->nv; i++) root[i] = root[root[i]];   /* (a root precedes its vertices: one pass flattens) */
		for (i = 0; i < s->npin; i++) pinned[root[d.pinned[i]]] = 1;
		for (i = 0, s->ngrp = 0; i < m->nv; i++)
			s->grp[i] = pinned[root[i]] ? ~0u : root[i] == i ? s->ngrp++ : s->grp[root[i]];
		for (i = 0; i < m->nv; i++) if (s->grp[i] != ~0u) s->gcnt[s->grp[i]]++;
		qaws_internal_free(root); qaws_internal_free(pinned);
		s->keep = d.translations == 2 && s->bary;
		s->bary = 0;   /* (the groups' rows replace the barycenter's) */
		k->ops = &rm_ops_translate;
	}
	for (i = 0; i < s->ne; i++)
	{
		double e[3];
		for (q = 0; q < 3; q++) e[q] = k->x[3 * s->edges[2 * i + 1] + q] - k->x[3 * s->edges[2 * i] + q];
		total += sqrt(rm_dot(e, e));
	}
	k->lmean = total / s->ne;
	rm_measures(m, k->x, &s->A0, NULL, &s->V0, NULL);
	s->vol_t = d.volume_target != 0 ? (double)d.volume_target : s->V0;
	k->theta = (double)d.theta;
	k->executor = m->executor;
	k->progress = d.progress;
	k->fixed = s->fixed;
	k->translates = s->npin == 0 && !s->bary && !s->grp;
	k->k = rm_brows(s) + (k->theta > 0 ? 0 : 3 * s->npin) + (s->area ? 1 : 0) + (s->volume ? 1 : 0);
	if (!(d.theta >= 0) || (k->theta > 0 && k->k > QAWS_INTERNAL_RS_HROWS) || (s->vol_w > 0 && !(s->vol_t != 0)) ||
		!qaws_internal_rs_settings(k, d.armijo, d.backtrack, d.max_trials, d.gmres_max_iterations, d.gmres_restart, d.gmres_tolerance))
	{
		qaws_repulsive_mesh_solver_destroy(s);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	k->clearance = d.clearance > 0 ? (double)d.clearance : 1e-3 * k->lmean;
	k->tol = d.tolerance > 0 ? (double)d.tolerance : 1e-4;
	k->ctol = (d.constraint_tolerance > 0 ? (double)d.constraint_tolerance : RM_CTOL) * k->lmean;
	k->h = (d.initial_step > 0 ? (double)d.initial_step : 0.25) * k->lmean;
	k->hmax = d.max_step > 0 ? (double)d.max_step * k->lmean : 0;
	k->rep.reason = QAWS_REPULSIVE_RUNNING;
	k->rep.stationarity = 1;
	if ((st = k->ops->energy(k, k->x, &E, NULL)) != QAWS_STATUS_OK)
	{
		qaws_repulsive_mesh_solver_destroy(s);
		return st;
	}
	k->rep.initial_energy = k->rep.energy = (qaws_scalar)E;
	if (k->checked && !rm_pairs_ok(s, k->x, k->x, 0))
	{
		/* (collision-checked: the input touches itself already) */
		qaws_repulsive_mesh_solver_destroy(s);
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
	*out_solver = s;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_mesh_solver_direction(qaws_scalar* out_direction, qaws_repulsive_mesh_solver const* s, qaws_scalar const* positions)
{
	return s ? qaws_internal_rs_direction(out_direction, &s->core, positions) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_mesh_solver_set_targets(qaws_repulsive_mesh_solver* s, qaws_scalar const* area, qaws_scalar const* volume)
{
	if (!s || (area && (!s->area || !((double)*area > 0))) || (volume && (!s->volume || !((double)*volume == (double)*volume))))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (area) s->A0 = (double)*area;
	if (volume) s->V0 = (double)*volume;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_mesh_solver_step(qaws_repulsive_report* out_report, qaws_scalar* positions, qaws_repulsive_mesh_solver* s)
{
	return s ? qaws_internal_rs_step(out_report, positions, &s->core) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_mesh_solver_objective(qaws_scalar* out_energy, qaws_scalar* out_gradient, qaws_repulsive_mesh_solver const* s,
	qaws_scalar const* positions)
{
	return s ? qaws_internal_rs_objective(out_energy, out_gradient, &s->core, positions) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_mesh_solve(qaws_repulsive_report* out_report, qaws_scalar* positions,
	qaws_repulsive_mesh const* mesh, qaws_repulsive_mesh_solve_desc const* desc)
{
	qaws_repulsive_mesh_solver* s = NULL;
	unsigned int max = desc && desc->max_iterations ? desc->max_iterations : 100;
	qaws_status st = qaws_repulsive_mesh_solver_create_ex(&s, mesh, positions, desc, mesh ? &mesh->alloc : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_rs_solve(out_report, positions, &s->core, max);
	qaws_repulsive_mesh_solver_destroy(s);
	return st;
}

qaws_status qaws_repulsive_mesh_solver_create_ex(qaws_repulsive_mesh_solver** out_solver, qaws_repulsive_mesh const* mesh,
	qaws_scalar const* positions, qaws_repulsive_mesh_solve_desc const* desc, qaws_allocator const* allocator)
{
	qaws_internal_rs_scope sc;
	qaws_status st;
	qaws_internal_rs_enter(&sc, allocator);
	st = qaws_repulsive_mesh_solver_create(out_solver, mesh, positions, desc);
	if (st == QAWS_STATUS_OK)
		qaws_internal_rs_keep(&(*out_solver)->core.alloc, allocator);
	qaws_internal_rs_leave(&sc);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  The preconditioner (Sec. 5.2) and a GMRES using it                      */
/* ------------------------------------------------------------------------ */

/* the fractional operator L^sg for the preconditioner, by vertex
   quadrature: sum over ordered pairs of distinct vertices of
   (u_i - u_j) (v_i - v_j) m_i m_j / r_ij^(2 sg + 2), m the lumped areas.
   Not the face-average form of Eq. 11: on a 3-colourable mesh (every
   subdivided octahedron or icosahedron) a vertex function differing
   between two colour classes averages to 0 on every face, so that form
   has a kernel beyond the constants and the preconditioned system is
   singular there. This form's kernel is the constants only. */
static qaws_status rm_fractional(qaws_repulsive_mesh const* m, rm_face const* f, double const* x, double sg, double* L)
{
	unsigned int i, j, s, c, nv = m->nv;
	double* mass = (double*)qaws_internal_scratch_zalloc(nv + 1, sizeof(double));
	if (!mass)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (s = 0; s < m->nf; s++)
		for (c = 0; c < 3; c++) mass[m->tris[3 * s + c]] += f[s].a / 3;
	memset(L, 0, sizeof(double) * nv * nv);
	for (i = 0; i < nv; i++)
		for (j = i + 1; j < nv; j++)
		{
			double d[3], r2, w;
			unsigned int q;
			for (q = 0; q < 3; q++) d[q] = x[3 * i + q] - x[3 * j + q];
			r2 = rm_dot(d, d);
			if (!(r2 > 0))
			{
				qaws_internal_free(mass);
				return QAWS_STATUS_NUMERICAL_FAILURE;
			}
			w = 2 * mass[i] * mass[j] * pow(r2, -(sg + 1));
			L[i * nv + i] += w; L[j * nv + j] += w;
			L[i * nv + j] -= w; L[j * nv + i] -= w;
		}
	qaws_internal_free(mass);
	return QAWS_STATUS_OK;
}

/* minus the cotan Laplacian (positive semidefinite, constants its kernel) */
static void rm_cotan(qaws_repulsive_mesh const* m, double const* x, double* K)
{
	unsigned int f, c, q, nv = m->nv;
	memset(K, 0, sizeof(double) * nv * nv);
	for (f = 0; f < m->nf; f++)
		for (c = 0; c < 3; c++)
		{
			/* the angle at corner c weighs the opposite edge (i, j) */
			unsigned int o = m->tris[3 * f + c], i = m->tris[3 * f + (c + 1) % 3], j = m->tris[3 * f + (c + 2) % 3];
			double u[3], v[3], n[3], w;
			for (q = 0; q < 3; q++) { u[q] = x[3 * i + q] - x[3 * o + q]; v[q] = x[3 * j + q] - x[3 * o + q]; }
			rm_cross(u, v, n);
			w = 0.5 * rm_dot(u, v) / sqrt(rm_dot(n, n));
			K[i * nv + j] -= w; K[j * nv + i] -= w;
			K[i * nv + i] += w; K[j * nv + j] += w;
		}
}

/* a factored n + 1 bordered system [K 1; 1^T 0]: its solve is K's
   pseudo-inverse on mean-zero right-hand sides */
typedef struct rm_pinv
{
	double* M;
	unsigned int* piv, n;
} rm_pinv;

static int rm_lu(double* m, unsigned int n, unsigned int* piv)
{
	unsigned int i, j, k;
	double big = 0;
	for (i = 0; i < n * n; i++) big = fmax(big, fabs(m[i]));
	for (k = 0; k < n; k++)
	{
		unsigned int p = k;
		double best = fabs(m[k * n + k]);
		for (i = k + 1; i < n; i++)
			if (fabs(m[i * n + k]) > best) { best = fabs(m[i * n + k]); p = i; }
		if (!(best > 1e-13 * big))
			return 0;
		piv[k] = p;
		if (p != k)
			for (j = 0; j < n; j++) { double t = m[k * n + j]; m[k * n + j] = m[p * n + j]; m[p * n + j] = t; }
		for (i = k + 1; i < n; i++)
		{
			double fct = m[i * n + k] / m[k * n + k];
			m[i * n + k] = fct;
			if (fct != 0)
				for (j = k + 1; j < n; j++) m[i * n + j] -= fct * m[k * n + j];
		}
	}
	return 1;
}

static void rm_lu_solve(double const* m, unsigned int n, unsigned int const* piv, double* b)
{
	unsigned int i, j, k;
	for (k = 0; k < n; k++)
		if (piv[k] != k) { double t = b[k]; b[k] = b[piv[k]]; b[piv[k]] = t; }
	for (i = 1; i < n; i++)
		for (j = 0; j < i; j++) b[i] -= m[i * n + j] * b[j];
	for (i = n; i-- > 0;)
	{
		for (j = i + 1; j < n; j++) b[i] -= m[i * n + j] * b[j];
		b[i] /= m[i * n + i];
	}
}

static int rm_pinv_make(rm_pinv* p, double const* K, unsigned int nv)
{
	unsigned int i, j, n = nv + 1;
	p->n = n;
	p->M = (double*)qaws_internal_zalloc(qaws_internal_mul2((size_t)n, n), sizeof(double));
	p->piv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), n));
	if (!p->M || !p->piv)
		return 0;
	for (i = 0; i < nv; i++)
	{
		for (j = 0; j < nv; j++) p->M[i * n + j] = K[i * nv + j];
		p->M[i * n + nv] = 1; p->M[nv * n + i] = 1;
	}
	return rm_lu(p->M, n, p->piv);
}

static void rm_pinv_apply(rm_pinv const* p, double* v, double* tmp)
{
	unsigned int i, nv = p->n - 1;
	double mean = 0;
	for (i = 0; i < nv; i++) mean += v[i] / nv;
	for (i = 0; i < nv; i++) tmp[i] = v[i] - mean;
	tmp[nv] = 0;
	rm_lu_solve(p->M, p->n, p->piv, tmp);
	for (i = 0; i < nv; i++) v[i] = tmp[i];
}

static void rm_matvec(double const* A, unsigned int n, double const* x, double* y)
{
	unsigned int i, j;
	for (i = 0; i < n; i++)
	{
		double s = 0;
		for (j = 0; j < n; j++) s += A[i * n + j] * x[j];
		y[i] = s;
	}
}

qaws_status qaws_internal_repulsive_mesh_metric(double* A, qaws_repulsive_mesh const* mesh, double const* x)
{
	return rm_metric(mesh, x, A);
}

/* A u = b (b mean-zero) by GMRES (restart 60), right-preconditioned by
   (-L)^+ L^(2 - s) (-L)^+ when asked; u mean-zero; the iterations taken */
qaws_status qaws_internal_repulsive_mesh_gmres(unsigned int* out_iterations, double* u, qaws_repulsive_mesh const* m, double const* x, double const* b,
	int precondition, double tol)
{
	enum { R = 60 };
	unsigned int nv = m->nv, it = 0, i, j, k, outer;
	double *A = NULL, *Lf = NULL, *K = NULL, *V = NULL, *Z = NULL, *H = NULL, *cs = NULL, *sn = NULL, *gq = NULL, *w = NULL, *t = NULL, *t2 = NULL, bn = 0;
	rm_face* f = NULL;
	rm_pinv P;
	qaws_status st = QAWS_STATUS_OK;
	memset(&P, 0, sizeof(P));
	A = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), nv, nv));
	Lf = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), nv, nv));
	K = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), nv, nv));
	V = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), nv, (R + 1)));
	Z = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), nv, (R + 1)));
	H = (double*)qaws_internal_scratch_zalloc(qaws_internal_mul2((R + 1), R), sizeof(double));
	cs = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), R)); sn = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), R));
	gq = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (R + 1)));
	w = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (nv + 1))); t = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (nv + 2))); t2 = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (nv + 2)));
	f = (rm_face*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(rm_face), (m->nf + 1)));
	if (!A || !Lf || !K || !V || !Z || !H || !cs || !sn || !gq || !w || !t || !t2 || !f)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if ((st = rm_metric(m, x, A)) != QAWS_STATUS_OK)
		goto done;
	if (!rm_faces(m, x, f))
	{
		st = QAWS_STATUS_DEGENERATE_CURVE;
		goto done;
	}
	if (precondition)
	{
		double s = 2 - 2 / m->p;
		if ((st = rm_fractional(m, f, x, 2 - s, Lf)) != QAWS_STATUS_OK)
			goto done;
		rm_cotan(m, x, K);
		if (!rm_pinv_make(&P, K, nv))
		{
			st = QAWS_STATUS_NUMERICAL_FAILURE;
			goto done;
		}
	}
	for (i = 0; i < nv; i++) { u[i] = 0; bn += b[i] * b[i]; }
	bn = sqrt(bn);
	if (!(bn > 0))
		goto done;
	for (outer = 0; outer < 50 && st == QAWS_STATUS_OK; outer++)
	{
		double beta = 0, res;
		/* r = b - A u */
		rm_matvec(A, nv, u, w);
		for (i = 0; i < nv; i++) { V[i] = b[i] - w[i]; beta += V[i] * V[i]; }
		beta = sqrt(beta);
		if (beta <= tol * bn)
			break;
		for (i = 0; i < nv; i++) V[i] /= beta;
		memset(gq, 0, sizeof(double) * (R + 1));
		gq[0] = beta;
		memset(H, 0, sizeof(double) * (R + 1) * R);
		for (j = 0; j < R; j++)
		{
			/* z = P v (right preconditioning), w = A z */
			double* z = Z + (size_t)j * nv;
			double* v = V + (size_t)j * nv;
			for (i = 0; i < nv; i++) z[i] = v[i];
			if (precondition)
			{
				rm_pinv_apply(&P, z, t);
				rm_matvec(Lf, nv, z, t2);
				for (i = 0; i < nv; i++) z[i] = t2[i];
				rm_pinv_apply(&P, z, t);
			}
			rm_matvec(A, nv, z, w);
			for (k = 0; k <= j; k++)
			{
				double hk = 0;
				for (i = 0; i < nv; i++) hk += w[i] * V[(size_t)k * nv + i];
				H[k * R + j] = hk;
				for (i = 0; i < nv; i++) w[i] -= hk * V[(size_t)k * nv + i];
			}
			H[(j + 1) * R + j] = 0;
			for (i = 0; i < nv; i++) H[(j + 1) * R + j] += w[i] * w[i];
			H[(j + 1) * R + j] = sqrt(H[(j + 1) * R + j]);
			if (H[(j + 1) * R + j] > 0)
				for (i = 0; i < nv; i++) V[(size_t)(j + 1) * nv + i] = w[i] / H[(j + 1) * R + j];
			/* Givens rotations */
			for (k = 0; k < j; k++)
			{
				double a = H[k * R + j], c2 = H[(k + 1) * R + j];
				H[k * R + j] = cs[k] * a + sn[k] * c2;
				H[(k + 1) * R + j] = -sn[k] * a + cs[k] * c2;
			}
			{
				double a = H[j * R + j], c2 = H[(j + 1) * R + j], r = hypot(a, c2);
				cs[j] = r > 0 ? a / r : 1; sn[j] = r > 0 ? c2 / r : 0;
				H[j * R + j] = r; H[(j + 1) * R + j] = 0;
				gq[j + 1] = -sn[j] * gq[j]; gq[j] = cs[j] * gq[j];
			}
			it++;
			res = fabs(gq[j + 1]);
			if (res <= tol * bn || j + 1 == R)
			{
				/* u += Z y, H y = g */
				double y[R];
				unsigned int n = j + 1;
				for (k = n; k-- > 0;)
				{
					double sm = gq[k];
					for (i = k + 1; i < n; i++) sm -= H[k * R + i] * y[i];
					y[k] = sm / H[k * R + k];
				}
				for (k = 0; k < n; k++)
					for (i = 0; i < nv; i++) u[i] += y[k] * Z[(size_t)k * nv + i];
				break;
			}
		}
	}
	{
		/* (the kernel's part dropped) */
		double mean = 0;
		for (i = 0; i < nv; i++) mean += u[i] / nv;
		for (i = 0; i < nv; i++) u[i] -= mean;
	}
done:
	if (out_iterations) *out_iterations = it;
	qaws_internal_free(A); qaws_internal_free(Lf); qaws_internal_free(K); qaws_internal_free(V); qaws_internal_free(Z); qaws_internal_free(H); qaws_internal_free(cs); qaws_internal_free(sn); qaws_internal_free(gq); qaws_internal_free(w); qaws_internal_free(t); qaws_internal_free(t2); qaws_internal_free(f);
	qaws_internal_free(P.M); qaws_internal_free(P.piv);
	return st;
}

int qaws_internal_repulsive_mesh_path_clear(qaws_repulsive_mesh const* mesh, double clearance, double const* x0, double const* x1)
{
	qaws_repulsive_mesh_solver s;
	int closed, ok;
	memset(&s, 0, sizeof(s));
	s.m = mesh; s.core.clearance = clearance; s.core.lmean = 1;
	if (rm_edges(mesh, &s.edges, &s.ne, &closed) != QAWS_STATUS_OK)
		return 0;
	ok = rm_pairs_ok(&s, x0, x1, 1);
	qaws_internal_free(s.edges);
	return ok;
}
