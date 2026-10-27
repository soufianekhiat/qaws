/*
 * Continuous collision detection (qaws_mesh_sweep.h): the first time a moving
 * mesh meets an obstacle, by slabs of time tested with each triangle's swept
 * prism against the obstacle's triangles, the first touching slab halved.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "../qaws_mesh_bvh.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct cc_ctx
{
	qaws_motion const* m;
	unsigned int nv, nt;
	unsigned int const* t;          /* the moving mesh's triangles */
	qaws_scalar const* body;        /* a rigid motion's points */
	qaws_scalar* buf;               /* scratch: nv points */
	double const* ox;               /* the obstacle's positions */
	unsigned int const* ot;         /* its triangles */
	unsigned int ont;
	qaws_mesh_bvh const* bvh;       /* over its triangles */
	unsigned int* cand;             /* scratch: box query results */
	unsigned int capc;
	double eps;
	unsigned long long tests;
} cc_ctx;

/* the moving mesh's vertices at t, in double */
static qaws_status cc_at(cc_ctx* c, double t, double* X)
{
	unsigned int i;
	qaws_status st = qaws_motion_is_per_vertex(c->m) ? qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)t, NULL, c->nv)
	                                                 : qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)t, c->body, c->nv);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < 3 * c->nv; i++)
		X[i] = (double)c->buf[i];
	return QAWS_STATUS_OK;
}

static double cc_orient(double const* a, double const* b, double const* c, double const* d)
{
	double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] }, w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
	return u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}

/* p inside the tetrahedron a b c d (on its border counts) */
static int cc_in_tet(double const* p, double const* a, double const* b, double const* c, double const* d)
{
	double s = cc_orient(a, b, c, d), s0, s1, s2, s3;
	if (s == 0)
		return 0;
	s0 = cc_orient(p, b, c, d);
	s1 = cc_orient(a, p, c, d);
	s2 = cc_orient(a, b, p, d);
	s3 = cc_orient(a, b, c, p);
	return s > 0 ? s0 >= 0 && s1 >= 0 && s2 >= 0 && s3 >= 0 : s0 <= 0 && s1 <= 0 && s2 <= 0 && s3 <= 0;
}

static int cc_degenerate(double const (*T)[3])
{
	double u[3], v[3], n[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = T[1][k] - T[0][k];
		v[k] = T[2][k] - T[0][k];
	}
	n[0] = u[1] * v[2] - u[2] * v[1];
	n[1] = u[2] * v[0] - u[0] * v[2];
	n[2] = u[0] * v[1] - u[1] * v[0];
	return n[0] == 0 && n[1] == 0 && n[2] == 0;
}

/* the obstacle's triangles whose boxes meet lo hi: into c->cand, their count; ~0u on allocation failure */
static unsigned int cc_query(cc_ctx* c, double const* lo, double const* hi)
{
	qaws_scalar l[3], h[3];
	unsigned int n = 0, k;
	for (k = 0; k < 3; k++)
	{
		l[k] = (qaws_scalar)(lo[k] - c->eps);
		h[k] = (qaws_scalar)(hi[k] + c->eps);
	}
	for (;;)
	{
		qaws_status st = qaws_mesh_bvh_query_box(c->cand, &n, c->bvh, l, h, c->capc);
		if (st == QAWS_STATUS_OK && n <= c->capc)
			return n;
		if (st != QAWS_STATUS_OK && st != QAWS_STATUS_BUFFER_TOO_SMALL)
			return ~0u;
		{
			unsigned int cap = n > 2 * c->capc ? n : 2 * c->capc + 64;
			unsigned int* g = (unsigned int*)qaws_internal_realloc(c->cand, sizeof(unsigned int) * cap);
			if (!g)
				return ~0u;
			c->cand = g;
			c->capc = cap;
		}
	}
}

/* the point of triangle a b c nearest p (Ericson's regions) */
static void cc_closest(double const* p, double const* a, double const* b, double const* c, double* o)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, den;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = b[k] - a[k]; ac[k] = c[k] - a[k]; ap[k] = p[k] - a[k]; bp[k] = p[k] - b[k]; cp[k] = p[k] - c[k];
	}
	d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
	d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
	d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
	d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
	d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
	d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
	if (d1 <= 0 && d2 <= 0) { memcpy(o, a, sizeof(double) * 3); return; }
	if (d3 >= 0 && d4 <= d3) { memcpy(o, b, sizeof(double) * 3); return; }
	if (d6 >= 0 && d5 <= d6) { memcpy(o, c, sizeof(double) * 3); return; }
	vc = d1 * d4 - d3 * d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) { v = d1 / (d1 - d3); for (k = 0; k < 3; k++) o[k] = a[k] + v * ab[k]; return; }
	vb = d5 * d2 - d1 * d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) { w = d2 / (d2 - d6); for (k = 0; k < 3; k++) o[k] = a[k] + w * ac[k]; return; }
	va = d3 * d6 - d5 * d4;
	if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); for (k = 0; k < 3; k++) o[k] = b[k] + w * (c[k] - b[k]); return; }
	den = 1 / (va + vb + vc);
	v = vb * den;
	w = vc * den;
	for (k = 0; k < 3; k++) o[k] = a[k] + v * ab[k] + w * ac[k];
}

/* does any triangle's swept prism between positions A and B (its caps, its sides split in two, its three tetrahedra)
   meet the obstacle? the first found into *f, *g, *p. -1 on allocation failure */
static int cc_touch(cc_ctx* c, double const* A, double const* B, unsigned int* f, unsigned int* g, double* p)
{
	unsigned int i;
	for (i = 0; i < c->nt; i++)
	{
		double P[6][3], lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, T[8][3][3];
		unsigned int k, j, q, nc, nT = 0;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
			{
				P[k][j] = A[3 * (size_t)c->t[3 * (size_t)i + k] + j];
				P[3 + k][j] = B[3 * (size_t)c->t[3 * (size_t)i + k] + j];
			}
		for (k = 0; k < 6; k++)
			for (j = 0; j < 3; j++)
			{
				lo[j] = fmin(lo[j], P[k][j]);
				hi[j] = fmax(hi[j], P[k][j]);
			}
		if ((nc = cc_query(c, lo, hi)) == ~0u)
			return -1;
		if (!nc)
			continue;
		/* the prism's border: the two caps and each side's two halves */
		memcpy(T[nT][0], P[0], sizeof P[0]); memcpy(T[nT][1], P[1], sizeof P[0]); memcpy(T[nT][2], P[2], sizeof P[0]); nT++;
		memcpy(T[nT][0], P[3], sizeof P[0]); memcpy(T[nT][1], P[4], sizeof P[0]); memcpy(T[nT][2], P[5], sizeof P[0]); nT++;
		for (k = 0; k < 3; k++)
		{
			unsigned int k1 = (k + 1) % 3;
			memcpy(T[nT][0], P[k], sizeof P[0]); memcpy(T[nT][1], P[k1], sizeof P[0]); memcpy(T[nT][2], P[3 + k1], sizeof P[0]);
			nT += !cc_degenerate((double const(*)[3])T[nT]);
			memcpy(T[nT][0], P[k], sizeof P[0]); memcpy(T[nT][1], P[3 + k1], sizeof P[0]); memcpy(T[nT][2], P[3 + k], sizeof P[0]);
			nT += !cc_degenerate((double const(*)[3])T[nT]);
		}
		for (q = 0; q < nc; q++)
		{
			unsigned int o = c->cand[q], z;
			double Q[3][3];
			static unsigned int const ip[3] = { 0, 1, 2 }, iq[3] = { 3, 4, 5 };
			for (k = 0; k < 3; k++)
				for (j = 0; j < 3; j++)
					Q[k][j] = c->ox[3 * (size_t)c->ot[3 * (size_t)o + k] + j];
			for (z = 0; z < nT; z++)
			{
				c->tests++;
				if (rp_tri_tri_at((double const(*)[3])T[z], (double const(*)[3])Q, ip, iq))
					break;
			}
			if (z == nT)
			{
				/* (inside the prism: an obstacle smaller than one step of a face) */
				for (k = 0; k < 3 && z == nT; k++)
					if (cc_in_tet(Q[k], P[0], P[1], P[2], P[3]) || cc_in_tet(Q[k], P[1], P[2], P[3], P[4]) || cc_in_tet(Q[k], P[2], P[3], P[4], P[5]))
						z = 0;
			}
			if (z < nT)
			{
				/* the contact's point: the nearest of the obstacle triangle to the moving one at the step's end (its corners), or of
				   it to the obstacle's corners, whichever pair is nearer */
				double best = HUGE_VAL, q2[3];
				*f = i;
				*g = o;
				for (k = 0; k < 3; k++)
				{
					double d2;
					cc_closest(P[3 + k], Q[0], Q[1], Q[2], q2);
					d2 = (q2[0] - P[3 + k][0]) * (q2[0] - P[3 + k][0]) + (q2[1] - P[3 + k][1]) * (q2[1] - P[3 + k][1]) + (q2[2] - P[3 + k][2]) * (q2[2] - P[3 + k][2]);
					if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + P[3 + k][j]); }
					cc_closest(Q[k], P[3], P[4], P[5], q2);
					d2 = (q2[0] - Q[k][0]) * (q2[0] - Q[k][0]) + (q2[1] - Q[k][1]) * (q2[1] - Q[k][1]) + (q2[2] - Q[k][2]) * (q2[2] - Q[k][2]);
					if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + Q[k][j]); }
				}
				return 1;
			}
		}
	}
	return 0;
}

/* the solid angle sum of the triangles T (vertices X) seen from p over 4 pi: 1 inside a closed surface, 0 outside */
static double cc_winding(double const* X, unsigned int const* T, unsigned int nt, double const* p)
{
	double w = 0;
	unsigned int i, k;
	for (i = 0; i < nt; i++)
	{
		double a[3], b[3], d[3], la, lb, ld, num, den;
		for (k = 0; k < 3; k++)
		{
			a[k] = X[3 * (size_t)T[3 * (size_t)i] + k] - p[k];
			b[k] = X[3 * (size_t)T[3 * (size_t)i + 1] + k] - p[k];
			d[k] = X[3 * (size_t)T[3 * (size_t)i + 2] + k] - p[k];
		}
		la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		lb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
		ld = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		num = a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) + a[2] * (b[0] * d[1] - b[1] * d[0]);
		den = la * lb * ld + (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) * ld + (a[0] * d[0] + a[1] * d[1] + a[2] * d[2]) * lb +
			(b[0] * d[0] + b[1] * d[1] + b[2] * d[2]) * la;
		w += 2 * atan2(num, den);
	}
	return w / (4 * 3.14159265358979323846);
}

/* the meshes meet at positions X: a pair of triangles crossing, or one inside the other */
static int cc_meet(cc_ctx* c, double const* X, unsigned int* f, unsigned int* g, double* p)
{
	int r = cc_touch(c, X, X, f, g, p);
	if (r)
		return r;
	if (c->nv && fabs(cc_winding(c->ox, c->ot, c->ont, X)) > 0.5)
	{
		*f = QAWS_MESH_NONE;
		*g = QAWS_MESH_NONE;
		memcpy(p, X, sizeof(double) * 3);
		return 1;
	}
	if (c->ont && fabs(cc_winding(X, c->t, c->nt, c->ox)) > 0.5)
	{
		*f = QAWS_MESH_NONE;
		*g = QAWS_MESH_NONE;
		memcpy(p, c->ox, sizeof(double) * 3);
		return 1;
	}
	return 0;
}

/* the first contact in [a, b] (positions XA, XB): halved while touching and wider than dt; its interval into *ta, *tb.
   0 none, 1 found, -1 a failure (*st) */
static int cc_first(cc_ctx* c, double a, double b, double const* XA, double const* XB, double dt, double* ta, double* tb, unsigned int* f, unsigned int* g,
	double* p, qaws_status* st)
{
	double* XM;
	double m;
	int r = cc_touch(c, XA, XB, f, g, p);
	if (r <= 0)
	{
		if (r < 0)
			*st = QAWS_STATUS_ALLOCATION_FAILURE;
		return r;
	}
	if (b - a <= dt)
	{
		*ta = a;
		*tb = b;
		return 1;
	}
	m = 0.5 * (a + b);
	XM = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c->nv + 1));
	if (!XM)
	{
		*st = QAWS_STATUS_ALLOCATION_FAILURE;
		return -1;
	}
	if ((*st = cc_at(c, m, XM)) != QAWS_STATUS_OK)
		r = -1;
	else
	{
		r = cc_first(c, a, m, XA, XM, dt, ta, tb, f, g, p, st);
		if (r == 0)
			r = cc_first(c, m, b, XM, XB, dt, ta, tb, f, g, p, st);
	}
	qaws_internal_free(XM);
	return r;
}

static int cc_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

qaws_status qaws_mesh_3d_first_contact(qaws_contact_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_mesh_3d const* obstacle,
	qaws_sweep_desc const* desc)
{
	qaws_contact_report local;
	cc_ctx c;
	qaws_mesh_bvh* bvh = NULL;
	double *X0 = NULL, *X1 = NULL, *Xm = NULL, *ox = NULL, *times = NULL, tol, diag = 0, lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
	double vmax = 0, dt, p[3] = { 0, 0, 0 };
	qaws_scalar* px = NULL;
	unsigned int *t = NULL, *ot = NULL, onv, i, k, f = QAWS_MESH_NONE, g = QAWS_MESH_NONE, ntimes = 0, captimes = 0, maxs;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	report->face = report->obstacle_face = QAWS_MESH_NONE;
	if (!mesh || !motion || !obstacle)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	c.m = motion;
	c.nv = qaws_mesh_3d_get_vertex_count(mesh);
	c.nt = qaws_mesh_3d_get_triangle_count(mesh);
	onv = qaws_mesh_3d_get_vertex_count(obstacle);
	c.ont = qaws_mesh_3d_get_triangle_count(obstacle);
	if (!c.nt || !c.ont || (qaws_motion_is_per_vertex(motion) && qaws_motion_get_vertex_count(motion) != c.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	px = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	c.buf = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
	t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.nt + 1));
	ot = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)c.ont + 1));
	ox = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)onv + 1));
	X0 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nv + 1));
	X1 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nv + 1));
	Xm = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.nv + 1));
	if (!px || !c.buf || !t || !ot || !ox || !X0 || !X1 || !Xm)
		goto done;
	qaws_mesh_3d_get_positions(px, mesh);
	qaws_mesh_3d_get_triangles(t, mesh);
	qaws_mesh_3d_get_triangles(ot, obstacle);
	{
		qaws_scalar* oxs = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)onv + 1));
		if (!oxs)
			goto done;
		qaws_mesh_3d_get_positions(oxs, obstacle);
		for (i = 0; i < 3 * onv; i++)
			ox[i] = (double)oxs[i];
		qaws_internal_free(oxs);
	}
	{
		/* (over the triangles, each its own element: polygons fanned as the triangles are read) */
		qaws_mesh_bvh_desc bd;
		qaws_scalar* oxs = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)onv + 1));
		if (!oxs)
			goto done;
		qaws_mesh_3d_get_positions(oxs, obstacle);
		memset(&bd, 0, sizeof bd);
		bd.positions = oxs;
		bd.vertex_count = onv;
		bd.dimension = 3;
		bd.simplices = ot;
		bd.simplex_size = 3;
		bd.simplex_count = c.ont;
		st = qaws_mesh_bvh_create(&bvh, &bd);
		qaws_internal_free(oxs);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	c.body = px;
	c.t = t;
	c.ox = ox;
	c.ot = ot;
	c.bvh = bvh;
	/* the tolerance from the moving mesh's box at 0 */
	if ((st = cc_at(&c, 0, X0)) != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < c.nv; i++)
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(lo[k], X0[3 * (size_t)i + k]);
			hi[k] = fmax(hi[k], X0[3 * (size_t)i + k]);
		}
	for (k = 0; k < 3; k++)
		diag += (hi[k] - lo[k]) * (hi[k] - lo[k]);
	diag = sqrt(diag);
	tol = desc && desc->tolerance > 0 ? (double)desc->tolerance : 1e-3 * diag;
	c.eps = 1e-9 * (diag > 0 ? diag : 1);
	maxs = desc && desc->max_slabs ? desc->max_slabs : 256;
	/* already meeting at 0 */
	st = QAWS_STATUS_OK;
	if (cc_meet(&c, X0, &f, &g, p))
	{
		report->contact = 1;
		report->time = report->time_upper = 0;
		goto found;
	}
	/* the fastest vertex (sampled), so the halving stops when a step moves under tol */
	for (i = 0; i <= 64; i++)
	{
		qaws_scalar* vel = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)c.nv + 1));
		if (!vel)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		st = qaws_motion_is_per_vertex(motion) ? qaws_motion_transform_points(c.buf, vel, motion, (qaws_scalar)(i / 64.0), NULL, c.nv)
		                                       : qaws_motion_transform_points(c.buf, vel, motion, (qaws_scalar)(i / 64.0), px, c.nv);
		for (k = 0; st == QAWS_STATUS_OK && k < c.nv; k++)
			vmax = fmax(vmax, sqrt((double)vel[3 * k] * vel[3 * k] + (double)vel[3 * k + 1] * vel[3 * k + 1] + (double)vel[3 * k + 2] * vel[3 * k + 2]));
		qaws_internal_free(vel);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	dt = vmax > 0 ? tol / (1.5 * vmax) : 1;
	/* the slabs: [0, 1] halved until the chords' deviation (at the middle and quarters) is within tol, as the sweep's */
	{
		double* stack = (double*)qaws_internal_malloc(sizeof(double) * 2 * (2 * (size_t)maxs + 2));
		unsigned int nsk = 0, ns = 1;
		if (!stack)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		stack[nsk++] = 0;
		stack[nsk++] = 1;
		while (nsk && st == QAWS_STATUS_OK)
		{
			double b = stack[--nsk], a = stack[--nsk], dev = 0;
			unsigned int q;
			if ((st = cc_at(&c, a, X0)) != QAWS_STATUS_OK || (st = cc_at(&c, b, X1)) != QAWS_STATUS_OK)
				break;
			for (q = 1; q <= 3 && st == QAWS_STATUS_OK; q++)
			{
				double fr = 0.25 * q;
				if ((st = cc_at(&c, a + fr * (b - a), Xm)) != QAWS_STATUS_OK)
					break;
				for (i = 0; i < c.nv; i++)
				{
					double d2 = 0;
					for (k = 0; k < 3; k++)
					{
						double e = Xm[3 * (size_t)i + k] - ((1 - fr) * X0[3 * (size_t)i + k] + fr * X1[3 * (size_t)i + k]);
						d2 += e * e;
					}
					dev = fmax(dev, sqrt(d2));
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
			if (dev > tol && ns < maxs)
			{
				stack[nsk++] = 0.5 * (a + b);
				stack[nsk++] = b;
				stack[nsk++] = a;
				stack[nsk++] = 0.5 * (a + b);
				ns++;
				continue;
			}
			if (ntimes + 2 > captimes)
			{
				unsigned int cap = captimes ? 2 * captimes : 64;
				double* gr = (double*)qaws_internal_realloc(times, sizeof(double) * cap);
				if (!gr)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					break;
				}
				times = gr;
				captimes = cap;
			}
			times[ntimes++] = a;
			times[ntimes++] = b;
		}
		qaws_internal_free(stack);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	qsort(times, ntimes / 2, 2 * sizeof(double), cc_cmp_double);
	report->slabs = ntimes / 2;
	/* the slabs in time order: the first that touches, halved */
	for (i = 0; i < ntimes / 2; i++)
	{
		double ta = 0, tb = 0;
		int r;
		if ((st = cc_at(&c, times[2 * (size_t)i], X0)) != QAWS_STATUS_OK || (st = cc_at(&c, times[2 * (size_t)i + 1], X1)) != QAWS_STATUS_OK)
			goto done;
		r = cc_first(&c, times[2 * (size_t)i], times[2 * (size_t)i + 1], X0, X1, dt, &ta, &tb, &f, &g, p, &st);
		if (r < 0)
			goto done;
		if (r > 0)
		{
			report->contact = 1;
			report->time = (qaws_scalar)ta;
			report->time_upper = (qaws_scalar)tb;
			goto found;
		}
	}
	report->time = report->time_upper = 1;
	goto done;
found:
	report->face = f;
	report->obstacle_face = g;
	report->point.x = (qaws_scalar)p[0];
	report->point.y = (qaws_scalar)p[1];
	report->point.z = (qaws_scalar)p[2];
done:
	report->tests = (unsigned int)(c.tests > 0xFFFFFFFFull ? 0xFFFFFFFFull : c.tests);
	qaws_mesh_bvh_destroy(bvh);
	qaws_internal_free(px); qaws_internal_free(c.buf); qaws_internal_free(t); qaws_internal_free(ot); qaws_internal_free(ox);
	qaws_internal_free(X0); qaws_internal_free(X1); qaws_internal_free(Xm); qaws_internal_free(times); qaws_internal_free(c.cand);
	return st;
}
