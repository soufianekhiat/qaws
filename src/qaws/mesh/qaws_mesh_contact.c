/*
 * Continuous collision detection (qaws_mesh_sweep.h): the first time a moving
 * mesh meets an obstacle (standing still or moving too), by slabs of time
 * tested with each triangle's swept prism against the obstacle's prisms, the
 * first touching slab halved; for straight moves the first contact time
 * solved exactly (a rational, on the exact kernel's integers).
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "../qaws_mesh_bvh.h"
#include "../exact/qaws_exact_int.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* a mesh and its motion (NULL: standing still) */
typedef struct cc_body
{
	qaws_motion const* m;
	unsigned int nv, nt;
	unsigned int* t;                /* its triangles */
	qaws_scalar* x;                 /* its positions as given */
	qaws_scalar* buf;               /* scratch: nv points */
} cc_body;

typedef struct cc_ctx
{
	cc_body a, o;                   /* the moving mesh, the obstacle */
	qaws_mesh_bvh* bvh;             /* over the obstacle's prisms: 8 triangles each, over its positions at a step's two ends */
	unsigned int* bsim;             /* (its simplices and their triangles) */
	unsigned int* bel;
	qaws_scalar* bpos;              /* 2 onv positions */
	double *OA, *OB;                /* the obstacle's positions at the step's ends */
	unsigned int* cand;             /* scratch: box query results */
	unsigned int capc;
	unsigned int* seen;             /* per obstacle triangle: the stamp it was last met with */
	unsigned int stamp;
	double eps;
	unsigned long long tests;
	unsigned long long* pairs;      /* with collect: the touching pairs (mesh triangle << 32 | obstacle triangle) */
	unsigned int npairs, cappairs;
	int collect;
} cc_ctx;

/* a body's vertices at t, in double */
static qaws_status cc_at(cc_body* b, double t, double* X)
{
	unsigned int i;
	qaws_status st = QAWS_STATUS_OK;
	if (!b->m)
	{
		for (i = 0; i < 3 * b->nv; i++)
			X[i] = (double)b->x[i];
		return QAWS_STATUS_OK;
	}
	st = qaws_motion_is_per_vertex(b->m) ? qaws_motion_transform_points(b->buf, NULL, b->m, (qaws_scalar)t, NULL, b->nv)
	                                     : qaws_motion_transform_points(b->buf, NULL, b->m, (qaws_scalar)t, b->x, b->nv);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < 3 * b->nv; i++)
		X[i] = (double)b->buf[i];
	return QAWS_STATUS_OK;
}

/* the obstacle at a and b: its positions, the BVH of its prisms refit (once for an obstacle standing still) */
static qaws_status cc_obstacle(cc_ctx* c, double a, double b)
{
	unsigned int i, n = c->o.nv;
	qaws_status st;
	if (!c->o.m && c->bvh)
		return QAWS_STATUS_OK;
	if ((st = cc_at(&c->o, a, c->OA)) != QAWS_STATUS_OK || (st = cc_at(&c->o, b, c->OB)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < 3 * n; i++)
	{
		c->bpos[i] = (qaws_scalar)c->OA[i];
		c->bpos[3 * (size_t)n + i] = (qaws_scalar)c->OB[i];
	}
	if (!c->bvh)
	{
		qaws_mesh_bvh_desc d;
		memset(&d, 0, sizeof d);
		d.positions = c->bpos;
		d.vertex_count = 2 * n;
		d.dimension = 3;
		d.simplices = c->bsim;
		d.simplex_size = 3;
		d.simplex_count = 8 * c->o.nt;
		d.elements = c->bel;
		return qaws_mesh_bvh_create(&c->bvh, &d);
	}
	return qaws_mesh_bvh_refit(c->bvh, c->bpos);
}

static double cc_orient(double const* a, double const* b, double const* c, double const* d)
{
	double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] }, w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
	return u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}

/* p inside the tetrahedron a b c d (on its border counts; flat ones hold nothing) */
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

/* triangle i's prism between positions X0 and X1: its 6 points into P, its border (caps, each side's halves; flat ones left
   out) into T; their count */
static unsigned int cc_prism(double const* X0, double const* X1, unsigned int const* t, unsigned int i, double (*P)[3], double (*T)[3][3])
{
	unsigned int k, j, n = 0;
	for (k = 0; k < 3; k++)
		for (j = 0; j < 3; j++)
		{
			P[k][j] = X0[3 * (size_t)t[3 * (size_t)i + k] + j];
			P[3 + k][j] = X1[3 * (size_t)t[3 * (size_t)i + k] + j];
		}
	memcpy(T[n][0], P[0], sizeof P[0]); memcpy(T[n][1], P[1], sizeof P[0]); memcpy(T[n][2], P[2], sizeof P[0]); n++;
	if (X0 == X1)
		return n;   /* (standing still: the triangle alone) */
	memcpy(T[n][0], P[3], sizeof P[0]); memcpy(T[n][1], P[4], sizeof P[0]); memcpy(T[n][2], P[5], sizeof P[0]); n++;
	for (k = 0; k < 3; k++)
	{
		unsigned int k1 = (k + 1) % 3;
		memcpy(T[n][0], P[k], sizeof P[0]); memcpy(T[n][1], P[k1], sizeof P[0]); memcpy(T[n][2], P[3 + k1], sizeof P[0]);
		n += !cc_degenerate((double const(*)[3])T[n]);
		memcpy(T[n][0], P[k], sizeof P[0]); memcpy(T[n][1], P[3 + k1], sizeof P[0]); memcpy(T[n][2], P[3 + k], sizeof P[0]);
		n += !cc_degenerate((double const(*)[3])T[n]);
	}
	return n;
}

static int cc_in_prism(double const* q, double const (*P)[3])
{
	return cc_in_tet(q, P[0], P[1], P[2], P[3]) || cc_in_tet(q, P[1], P[2], P[3], P[4]) || cc_in_tet(q, P[2], P[3], P[4], P[5]);
}

/* the obstacle's triangles whose prisms' boxes meet lo hi (each once): into c->cand, their count; ~0u on allocation failure */
static unsigned int cc_query(cc_ctx* c, double const* lo, double const* hi)
{
	qaws_scalar l[3], h[3];
	unsigned int n = 0, k, m = 0;
	for (k = 0; k < 3; k++)
	{
		l[k] = (qaws_scalar)(lo[k] - c->eps);
		h[k] = (qaws_scalar)(hi[k] + c->eps);
	}
	for (;;)
	{
		qaws_status st = qaws_mesh_bvh_query_box(c->cand, &n, c->bvh, l, h, c->capc);
		if (st == QAWS_STATUS_OK && n <= c->capc)
			break;
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
	/* (a prism's 8 triangles share its element: each once) */
	c->stamp++;
	for (k = 0; k < n; k++)
		if (c->seen[c->cand[k]] != c->stamp)
		{
			c->seen[c->cand[k]] = c->stamp;
			c->cand[m++] = c->cand[k];
		}
	return m;
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

/* does any triangle's swept prism between positions A and B meet one of the obstacle's (at c->OA, c->OB): their borders
   crossing, or a corner of one inside the other? the first found into *f, *g, *p (with collect: every touching pair into
   c->pairs, 1 when any). -1 on allocation failure */
static int cc_touch(cc_ctx* c, double const* A, double const* B, unsigned int* f, unsigned int* g, double* p)
{
	unsigned int i;
	int still = c->o.m == NULL, found = 0;
	for (i = 0; i < c->a.nt; i++)
	{
		double P[6][3], lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, T[8][3][3];
		unsigned int k, j, q, nc, nT;
		nT = cc_prism(A, B, c->a.t, i, P, T);
		for (k = 0; k < 6; k++)
			for (j = 0; j < 3; j++)
			{
				lo[j] = fmin(lo[j], P[k][j]);
				hi[j] = fmax(hi[j], P[k][j]);
			}
		if ((nc = cc_query(c, lo, hi)) == ~0u)
			return -1;
		for (q = 0; q < nc; q++)
		{
			unsigned int o = c->cand[q], z, w, nU, hit = 0;
			double Q[6][3], U[8][3][3];
			static unsigned int const ip[3] = { 0, 1, 2 }, iq[3] = { 3, 4, 5 };
			nU = cc_prism(c->OA, still ? c->OA : c->OB, c->o.t, o, Q, U);
			for (z = 0; z < nT && !hit; z++)
				for (w = 0; w < nU && !hit; w++)
				{
					c->tests++;
					hit = rp_tri_tri_at((double const(*)[3])T[z], (double const(*)[3])U[w], ip, iq);
				}
			/* (a corner inside the other's prism: one smaller than the other's step) */
			for (k = 0; k < (still ? 3u : 6u) && !hit; k++)
				hit = cc_in_prism(Q[k], (double const(*)[3])P);
			for (k = 0; k < 6 && !hit && !still; k++)
				hit = cc_in_prism(P[k], (double const(*)[3])Q);
			if (!hit)
				continue;
			if (c->collect)
			{
				if (c->npairs + 1 > c->cappairs)
				{
					unsigned int cap = 2 * c->cappairs + 64;
					unsigned long long* gp = (unsigned long long*)qaws_internal_realloc(c->pairs, sizeof(unsigned long long) * cap);
					if (!gp)
						return -1;
					c->pairs = gp;
					c->cappairs = cap;
				}
				c->pairs[c->npairs++] = (unsigned long long)i << 32 | o;
				found = 1;
				continue;
			}
			{
				/* the contact's point: between the nearest pair of the two triangles at the step's end (from their corners) */
				double best = HUGE_VAL, q2[3], *E = still ? Q[0] : Q[3];
				double const (*Qe)[3] = still ? (double const(*)[3])Q : (double const(*)[3])(Q + 3);
				(void)E;
				*f = i;
				*g = o;
				for (k = 0; k < 3; k++)
				{
					double d2;
					cc_closest(P[3 + k], Qe[0], Qe[1], Qe[2], q2);
					d2 = (q2[0] - P[3 + k][0]) * (q2[0] - P[3 + k][0]) + (q2[1] - P[3 + k][1]) * (q2[1] - P[3 + k][1]) + (q2[2] - P[3 + k][2]) * (q2[2] - P[3 + k][2]);
					if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + P[3 + k][j]); }
					cc_closest(Qe[k], P[3], P[4], P[5], q2);
					d2 = (q2[0] - Qe[k][0]) * (q2[0] - Qe[k][0]) + (q2[1] - Qe[k][1]) * (q2[1] - Qe[k][1]) + (q2[2] - Qe[k][2]) * (q2[2] - Qe[k][2]);
					if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + Qe[k][j]); }
				}
				return 1;
			}
		}
	}
	return found;
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

/* the meshes meet at 0 (the mesh at X, the obstacle at c->OA): a pair of triangles crossing, or one inside the other */
static int cc_meet(cc_ctx* c, double const* X, unsigned int* f, unsigned int* g, double* p)
{
	int r;
	double* keep = c->OB;
	c->OB = c->OA;
	r = cc_touch(c, X, X, f, g, p);
	c->OB = keep;
	if (r)
		return r;
	if (fabs(cc_winding(c->OA, c->o.t, c->o.nt, X)) > 0.5 || fabs(cc_winding(X, c->a.t, c->a.nt, c->OA)) > 0.5)
	{
		*f = QAWS_MESH_NONE;
		*g = QAWS_MESH_NONE;
		memcpy(p, fabs(cc_winding(c->OA, c->o.t, c->o.nt, X)) > 0.5 ? X : c->OA, sizeof(double) * 3);
		return 1;
	}
	return 0;
}

/* the first contact in [a, b] (the mesh at XA, XB): halved while touching and wider than dt; its interval into *ta, *tb.
   0 none, 1 found, -1 a failure (*st) */
static int cc_first(cc_ctx* c, double a, double b, double const* XA, double const* XB, double dt, double* ta, double* tb, unsigned int* f, unsigned int* g,
	double* p, qaws_status* st)
{
	double* XM;
	double m;
	int r;
	if ((*st = cc_obstacle(c, a, b)) != QAWS_STATUS_OK)
		return -1;
	r = cc_touch(c, XA, XB, f, g, p);
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
	XM = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c->a.nv + 1));
	if (!XM)
	{
		*st = QAWS_STATUS_ALLOCATION_FAILURE;
		return -1;
	}
	if ((*st = cc_at(&c->a, m, XM)) != QAWS_STATUS_OK)
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

/* ---- straight moves, exactly: every vertex at q + T(t), T(t) = T0 + t (T1 - T0) (q the body's points turned once, all
   doubles), so a vertex reaching a face's plane or two edges becoming coplanar is linear in t: the first contact time a
   rational, found on the exact kernel's integers (coordinates as integers times 2^E, E the least exponent) ---- */

typedef struct ccx_v
{
	qaws_exact_int c[3];
} ccx_v;

/* a body's straight move: its rotation fixed, its translation linear (positions at 0.25, 0.5, 0.75 the chords', its
   velocity its displacement); into R (by rows), T0, T1 */
static int ccx_straight(qaws_motion const* m, double* R, double* T0, double* T1)
{
	qaws_motion_state s[5];
	unsigned int i, k;
	if (!m)
	{
		for (k = 0; k < 9; k++)
			R[k] = k % 4 == 0;
		T0[0] = T0[1] = T0[2] = T1[0] = T1[1] = T1[2] = 0;
		return 1;
	}
	if (qaws_motion_is_per_vertex(m))
		return 0;
	for (i = 0; i < 5; i++)
		if (qaws_motion_evaluate(&s[i], m, (qaws_scalar)(0.25 * i)) != QAWS_STATUS_OK)
			return 0;
	for (i = 0; i < 5; i++)
	{
		double p[3] = { s[i].position.x, s[i].position.y, s[i].position.z }, v[3] = { s[i].velocity.x, s[i].velocity.y, s[i].velocity.z };
		double w = fabs((double)s[i].angular_velocity.x) + fabs((double)s[i].angular_velocity.y) + fabs((double)s[i].angular_velocity.z);
		double d[3] = { (double)s[4].position.x - (double)s[0].position.x, (double)s[4].position.y - (double)s[0].position.y, (double)s[4].position.z - (double)s[0].position.z };
		double sc = 1e-12 * (1 + fabs(d[0]) + fabs(d[1]) + fabs(d[2]) + fabs((double)s[0].position.x) + fabs((double)s[0].position.y) + fabs((double)s[0].position.z));
		if (w != 0)
			return 0;
		for (k = 0; k < 9; k++)
			if (s[i].rotation[k] != s[0].rotation[k])
				return 0;
		for (k = 0; k < 3; k++)
		{
			double p0 = k == 0 ? s[0].position.x : k == 1 ? s[0].position.y : s[0].position.z;
			if (fabs(p[k] - (p0 + 0.25 * i * d[k])) > sc || fabs(v[k] - d[k]) > sc)
				return 0;
		}
	}
	for (k = 0; k < 9; k++)
		R[k] = (double)s[0].rotation[k];
	T0[0] = s[0].position.x; T0[1] = s[0].position.y; T0[2] = s[0].position.z;
	T1[0] = s[4].position.x; T1[1] = s[4].position.y; T1[2] = s[4].position.z;
	return 1;
}

static int ccx_ok(qaws_status s)
{
	return s == QAWS_STATUS_OK;
}

/* d = m 2^(e - E) as an integer (E at most e) */
static int ccx_from(qaws_exact_int* r, double d, int E)
{
	int64_t m;
	int e;
	if (d == 0)
	{
		qaws_exact_int_zero(r);
		return 1;
	}
	if (!qaws_exact_split_double(d, &m, &e) || e < E)
		return 0;
	qaws_exact_int_from_i64(r, m);
	return ccx_ok(qaws_exact_int_shl(r, r, (unsigned int)(e - E)));
}

static int ccx_sub(ccx_v* r, ccx_v const* a, ccx_v const* b)
{
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (!ccx_ok(qaws_exact_int_sub(&r->c[k], &a->c[k], &b->c[k])))
			return 0;
	return 1;
}

static int ccx_cross(ccx_v* r, ccx_v const* a, ccx_v const* b)
{
	qaws_exact_int x, y;
	ccx_v o;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		unsigned int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
		if (!ccx_ok(qaws_exact_int_mul(&x, &a->c[k1], &b->c[k2])) || !ccx_ok(qaws_exact_int_mul(&y, &a->c[k2], &b->c[k1])) || !ccx_ok(qaws_exact_int_sub(&o.c[k], &x, &y)))
			return 0;
	}
	*r = o;
	return 1;
}

static int ccx_dot(qaws_exact_int* r, ccx_v const* a, ccx_v const* b)
{
	qaws_exact_int x;
	unsigned int k;
	qaws_exact_int_zero(r);
	for (k = 0; k < 3; k++)
		if (!ccx_ok(qaws_exact_int_mul(&x, &a->c[k], &b->c[k])) || !ccx_ok(qaws_exact_int_add(r, r, &x)))
			return 0;
	return 1;
}

/* r = a s + b u (vectors a, b; integers s, u) */
static int ccx_comb(ccx_v* r, ccx_v const* a, qaws_exact_int const* s, ccx_v const* b, qaws_exact_int const* u)
{
	qaws_exact_int x, y;
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (!ccx_ok(qaws_exact_int_mul(&x, &a->c[k], s)) || !ccx_ok(qaws_exact_int_mul(&y, &b->c[k], u)) || !ccx_ok(qaws_exact_int_add(&r->c[k], &x, &y)))
			return 0;
	return 1;
}

/* n . (e x w): the sign of w past edge e in the plane of normal n */
static int ccx_side(int* s, ccx_v const* n, ccx_v const* e, ccx_v const* w)
{
	ccx_v x;
	qaws_exact_int d;
	if (!ccx_cross(&x, e, w) || !ccx_dot(&d, n, &x))
		return 0;
	*s = qaws_exact_int_sign(&d);
	return 1;
}

/* num / den < bn / bd (dens positive) */
static int ccx_less(qaws_exact_int const* num, qaws_exact_int const* den, qaws_exact_int const* bn, qaws_exact_int const* bd, int* less)
{
	qaws_exact_int x, y;
	if (!ccx_ok(qaws_exact_int_mul(&x, num, bd)) || !ccx_ok(qaws_exact_int_mul(&y, bn, den)))
		return 0;
	*less = qaws_exact_int_cmp(&x, &y) < 0;
	return 1;
}

/* the time t = num / den (den > 0, 0 <= t <= 1) a point a moving by D reaches the triangle b0 b1 b2 standing still, inside
   it (its border counts); 0 none, 1 found, -1 out of range */
static int ccx_vf(ccx_v const* a, ccx_v const* D, ccx_v const* b0, ccx_v const* b1, ccx_v const* b2, qaws_exact_int* num, qaws_exact_int* den)
{
	ccx_v e1, e2, n, w, P, B[3], E[3];
	qaws_exact_int o0, one;
	int s[3], k;
	ccx_v const* bb[3] = { b0, b1, b2 };
	if (!ccx_sub(&e1, b1, b0) || !ccx_sub(&e2, b2, b0) || !ccx_cross(&n, &e1, &e2) || !ccx_sub(&w, a, b0) || !ccx_dot(&o0, &n, &w) || !ccx_dot(den, &n, D))
		return -1;
	if (qaws_exact_int_is_zero(den))
		return 0;   /* (moving along the plane: its edges meet the triangle's first) */
	qaws_exact_int_neg(num, &o0);
	if (qaws_exact_int_sign(den) < 0)
	{
		qaws_exact_int_neg(num, num);
		qaws_exact_int_neg(den, den);
	}
	if (qaws_exact_int_sign(num) < 0 || qaws_exact_int_cmp(num, den) > 0)
		return 0;
	/* inside at t: den p(t) = a den + num D against the triangle times den */
	if (!ccx_comb(&P, a, den, D, num))
		return -1;
	qaws_exact_int_from_i64(&one, 0);
	for (k = 0; k < 3; k++)
	{
		ccx_v z;
		unsigned int q;
		for (q = 0; q < 3; q++)
			z.c[q] = one;
		if (!ccx_comb(&B[k], bb[k], den, &z, &one) || !ccx_sub(&E[k], bb[(k + 1) % 3], bb[k]))
			return -1;
	}
	for (k = 0; k < 3; k++)
	{
		ccx_v r;
		if (!ccx_sub(&r, &P, &B[k]) || !ccx_side(&s[k], &n, &E[k], &r))
			return -1;
	}
	return (s[0] >= 0 && s[1] >= 0 && s[2] >= 0) || (s[0] <= 0 && s[1] <= 0 && s[2] <= 0) ? 1 : 0;
}

/* the time edge a0 a1 moving by D meets edge b0 b1 standing still (not parallel), within both; 0 none, 1 found, -1 out of
   range */
static int ccx_ee(ccx_v const* a0, ccx_v const* a1, ccx_v const* D, ccx_v const* b0, ccx_v const* b1, qaws_exact_int* num, qaws_exact_int* den)
{
	ccx_v ea, eb, n, w0, w1, x, A0, A1, B0, B1, r;
	qaws_exact_int zero;
	int s1, s2, s3, s4;
	unsigned int q;
	if (!ccx_sub(&ea, a1, a0) || !ccx_sub(&eb, b1, b0) || !ccx_cross(&n, &ea, &eb))
		return -1;
	if (qaws_exact_int_is_zero(&n.c[0]) && qaws_exact_int_is_zero(&n.c[1]) && qaws_exact_int_is_zero(&n.c[2]))
		return 0;   /* (parallel: their ends meet the other's faces first) */
	/* coplanar when (b0 - a(t)) . n = 0: (b0 - a0) . n - t D . n */
	if (!ccx_sub(&w0, b0, a0) || !ccx_dot(num, &w0, &n) || !ccx_dot(den, D, &n))
		return -1;
	if (qaws_exact_int_is_zero(den))
		return 0;
	if (qaws_exact_int_sign(den) < 0)
	{
		qaws_exact_int_neg(num, num);
		qaws_exact_int_neg(den, den);
	}
	if (qaws_exact_int_sign(num) < 0 || qaws_exact_int_cmp(num, den) > 0)
		return 0;
	/* at t, times den: crossing in their plane */
	qaws_exact_int_zero(&zero);
	for (q = 0; q < 3; q++)
		x.c[q] = zero;
	if (!ccx_comb(&A0, a0, den, D, num) || !ccx_comb(&A1, a1, den, D, num) || !ccx_comb(&B0, b0, den, &x, &zero) || !ccx_comb(&B1, b1, den, &x, &zero))
		return -1;
	if (!ccx_sub(&r, &B0, &A0) || !ccx_side(&s1, &n, &ea, &r) || !ccx_sub(&r, &B1, &A0) || !ccx_side(&s2, &n, &ea, &r))
		return -1;
	if (!ccx_sub(&r, &A0, &B0) || !ccx_side(&s3, &n, &eb, &r) || !ccx_sub(&r, &A1, &B0) || !ccx_side(&s4, &n, &eb, &r))
		return -1;
	(void)w1;
	return s1 * s2 <= 0 && s3 * s4 <= 0 ? 1 : 0;
}

/* the exact first contact among the pairs in c->pairs (the mesh's vertices at qa + Ta(t), the obstacle's at qo + To(t)):
   the earliest event into *num / *den and its pair; 0 none, 1 found, -1 out of range or allocation */
static int ccx_first(cc_ctx* c, double const* qa, double const* Ta0, double const* Ta1, double const* qo, double const* To0, double const* To1,
	qaws_exact_int* bnum, qaws_exact_int* bden, unsigned int* bf, unsigned int* bg)
{
	int E = 0x7FFFFFFF, have = 0;
	unsigned int i, k, z;
	ccx_v D, Da, Do;
	double const* lists[6];
	unsigned int counts[6];
	/* the least exponent of every double in play */
	lists[0] = qa; counts[0] = 3 * c->a.nv;
	lists[1] = Ta0; counts[1] = 3;
	lists[2] = Ta1; counts[2] = 3;
	lists[3] = qo; counts[3] = 3 * c->o.nv;
	lists[4] = To0; counts[4] = 3;
	lists[5] = To1; counts[5] = 3;
	for (z = 0; z < 6; z++)
		for (i = 0; i < counts[z]; i++)
		{
			int64_t m;
			int e;
			if (lists[z][i] != 0 && qaws_exact_split_double(lists[z][i], &m, &e) && e < E)
				E = e;
		}
	if (E == 0x7FFFFFFF)
		E = 0;
	{
		ccx_v t0, t1;
		for (k = 0; k < 3; k++)
			if (!ccx_from(&t0.c[k], Ta0[k], E) || !ccx_from(&t1.c[k], Ta1[k], E))
				return -1;
		if (!ccx_sub(&Da, &t1, &t0))
			return -1;
		for (k = 0; k < 3; k++)
			if (!ccx_from(&t0.c[k], To0[k], E) || !ccx_from(&t1.c[k], To1[k], E))
				return -1;
		if (!ccx_sub(&Do, &t1, &t0))
			return -1;
		if (!ccx_sub(&D, &Da, &Do))   /* (the mesh's motion against the obstacle's) */
			return -1;
	}
	for (i = 0; i < c->npairs; i++)
	{
		unsigned int f = (unsigned int)(c->pairs[i] >> 32), g = (unsigned int)c->pairs[i], q;
		ccx_v A[3], O[3], ND;
		/* the two triangles' corners at 0 */
		for (q = 0; q < 3; q++)
			for (k = 0; k < 3; k++)
			{
				qaws_exact_int x, y;
				unsigned int va = c->a.t[3 * (size_t)f + q], vo = c->o.t[3 * (size_t)g + q];
				if (!ccx_from(&x, qa[3 * (size_t)va + k], E) || !ccx_from(&y, Ta0[k], E) || !ccx_ok(qaws_exact_int_add(&A[q].c[k], &x, &y)))
					return -1;
				if (!ccx_from(&x, qo[3 * (size_t)vo + k], E) || !ccx_from(&y, To0[k], E) || !ccx_ok(qaws_exact_int_add(&O[q].c[k], &x, &y)))
					return -1;
			}
		for (k = 0; k < 3; k++)
			qaws_exact_int_neg(&ND.c[k], &D.c[k]);
		/* the mesh's corners against the obstacle's face, the obstacle's against the mesh's (moving back), the edges */
		for (q = 0; q < 15; q++)
		{
			qaws_exact_int num, den;
			int r, less = 1;
			if (q < 3)
				r = ccx_vf(&A[q], &D, &O[0], &O[1], &O[2], &num, &den);
			else if (q < 6)
				r = ccx_vf(&O[q - 3], &ND, &A[0], &A[1], &A[2], &num, &den);
			else
				r = ccx_ee(&A[(q - 6) / 3], &A[((q - 6) / 3 + 1) % 3], &D, &O[(q - 6) % 3], &O[((q - 6) % 3 + 1) % 3], &num, &den);
			if (r < 0)
				return -1;
			if (!r)
				continue;
			if (have && !ccx_less(&num, &den, bnum, bden, &less))
				return -1;
			if (less)
			{
				*bnum = num;
				*bden = den;
				*bf = f;
				*bg = g;
				have = 1;
			}
		}
	}
	return have;
}

static int cc_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

/* a body read: its triangles, positions, motion (NULL: still); INVALID_ARGUMENT for a vertex motion of another count */
static qaws_status cc_body_init(cc_body* b, qaws_mesh_3d const* m, qaws_motion const* mo)
{
	memset(b, 0, sizeof *b);
	b->m = mo;
	b->nv = qaws_mesh_3d_get_vertex_count(m);
	b->nt = qaws_mesh_3d_get_triangle_count(m);
	if (!b->nt || (mo && qaws_motion_is_per_vertex(mo) && qaws_motion_get_vertex_count(mo) != b->nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	b->t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (3 * (size_t)b->nt + 1));
	b->x = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)b->nv + 1));
	b->buf = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)b->nv + 1));
	if (!b->t || !b->x || !b->buf)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	qaws_mesh_3d_get_triangles(b->t, m);
	qaws_mesh_3d_get_positions(b->x, m);
	return QAWS_STATUS_OK;
}

static void cc_body_free(cc_body* b)
{
	qaws_internal_free(b->t);
	qaws_internal_free(b->x);
	qaws_internal_free(b->buf);
}

/* the fastest vertex of a body (sampled at 65 times) */
static qaws_status cc_speed(cc_body* b, double* vmax)
{
	unsigned int i, k;
	qaws_status st = QAWS_STATUS_OK;
	qaws_scalar* vel;
	if (!b->m)
		return QAWS_STATUS_OK;
	vel = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (3 * (size_t)b->nv + 1));
	if (!vel)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i <= 64 && st == QAWS_STATUS_OK; i++)
	{
		st = qaws_motion_is_per_vertex(b->m) ? qaws_motion_transform_points(b->buf, vel, b->m, (qaws_scalar)(i / 64.0), NULL, b->nv)
		                                     : qaws_motion_transform_points(b->buf, vel, b->m, (qaws_scalar)(i / 64.0), b->x, b->nv);
		for (k = 0; st == QAWS_STATUS_OK && k < b->nv; k++)
			*vmax = fmax(*vmax, sqrt((double)vel[3 * k] * vel[3 * k] + (double)vel[3 * k + 1] * vel[3 * k + 1] + (double)vel[3 * k + 2] * vel[3 * k + 2]));
	}
	qaws_internal_free(vel);
	return st;
}

qaws_status qaws_mesh_3d_first_contact(qaws_contact_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_mesh_3d const* obstacle,
	qaws_motion const* obstacle_motion, qaws_sweep_desc const* desc)
{
	qaws_contact_report local;
	cc_ctx c;
	double *X0 = NULL, *X1 = NULL, *Xm = NULL, *O0 = NULL, *O1 = NULL, *Om = NULL, *times = NULL, tol, diag = 0;
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, vmax = 0, dt, p[3] = { 0, 0, 0 };
	unsigned int i, k, f = QAWS_MESH_NONE, g = QAWS_MESH_NONE, ntimes = 0, captimes = 0, maxs;
	qaws_status st;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	report->face = report->obstacle_face = QAWS_MESH_NONE;
	report->time = report->time_upper = 1;
	if (!mesh || !motion || !obstacle)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	if ((st = cc_body_init(&c.a, mesh, motion)) != QAWS_STATUS_OK || (st = cc_body_init(&c.o, obstacle, obstacle_motion)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	X0 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.a.nv + 1));
	X1 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.a.nv + 1));
	Xm = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.a.nv + 1));
	O0 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
	O1 = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
	Om = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
	c.OA = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
	c.OB = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
	c.bpos = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * (6 * (size_t)c.o.nv + 1));
	c.bsim = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (24 * (size_t)c.o.nt + 1));
	c.bel = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (8 * (size_t)c.o.nt + 1));
	c.seen = (unsigned int*)qaws_internal_zalloc((size_t)c.o.nt + 1, sizeof(unsigned int));
	if (!X0 || !X1 || !Xm || !O0 || !O1 || !Om || !c.OA || !c.OB || !c.bpos || !c.bsim || !c.bel || !c.seen)
		goto done;
	/* the obstacle's prisms in its BVH: per triangle its two caps and its sides' halves over the positions at both ends */
	for (i = 0; i < c.o.nt; i++)
	{
		unsigned int const* v = c.o.t + 3 * (size_t)i;
		unsigned int n = c.o.nv, *s = c.bsim + 24 * (size_t)i;
		s[0] = v[0]; s[1] = v[1]; s[2] = v[2];
		s[3] = n + v[0]; s[4] = n + v[1]; s[5] = n + v[2];
		for (k = 0; k < 3; k++)
		{
			unsigned int x = v[k], y = v[(k + 1) % 3];
			s[6 + 6 * k] = x; s[7 + 6 * k] = y; s[8 + 6 * k] = n + y;
			s[9 + 6 * k] = x; s[10 + 6 * k] = n + y; s[11 + 6 * k] = n + x;
		}
		for (k = 0; k < 8; k++)
			c.bel[8 * (size_t)i + k] = i;
	}
	/* the tolerance from the moving mesh's box at 0 */
	if ((st = cc_at(&c.a, 0, X0)) != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < c.a.nv; i++)
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
	if ((st = cc_obstacle(&c, 0, 0)) != QAWS_STATUS_OK)
		goto done;
	{
		int r = cc_meet(&c, X0, &f, &g, p);
		if (r < 0)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		if (r)
		{
			report->contact = 1;
			report->exact = 1;
			report->time = report->time_upper = 0;
			goto found;
		}
	}
	/* the fastest vertex of either (sampled): the halving stops when a step moves none past tol */
	if ((st = cc_speed(&c.a, &vmax)) != QAWS_STATUS_OK || (st = cc_speed(&c.o, &vmax)) != QAWS_STATUS_OK)
		goto done;
	dt = vmax > 0 ? tol / (1.5 * vmax) : 1;
	/* the slabs: [0, 1] halved until both bodies' chords (at the middle and quarters) are within tol, as the sweep's */
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
			unsigned int q, bo;
			for (bo = 0; bo < 2 && st == QAWS_STATUS_OK; bo++)
			{
				cc_body* B = bo ? &c.o : &c.a;
				double *Y0 = bo ? O0 : X0, *Y1 = bo ? O1 : X1, *Ym = bo ? Om : Xm;
				if (!B->m)
					continue;
				if ((st = cc_at(B, a, Y0)) != QAWS_STATUS_OK || (st = cc_at(B, b, Y1)) != QAWS_STATUS_OK)
					break;
				for (q = 1; q <= 3 && st == QAWS_STATUS_OK; q++)
				{
					double fr = 0.25 * q;
					if ((st = cc_at(B, a + fr * (b - a), Ym)) != QAWS_STATUS_OK)
						break;
					for (i = 0; i < B->nv; i++)
					{
						double d2 = 0;
						for (k = 0; k < 3; k++)
						{
							double e = Ym[3 * (size_t)i + k] - ((1 - fr) * Y0[3 * (size_t)i + k] + fr * Y1[3 * (size_t)i + k]);
							d2 += e * e;
						}
						dev = fmax(dev, sqrt(d2));
					}
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
	{
		/* the slabs in time order: the first that touches, halved; for straight moves the exact time among the pairs touching
		   about it, the search going on past the interval when none of them meets */
		double Ra[9], Ta0[3], Ta1[3], Ro[9], To0[3], To1[3];
		int straight = ccx_straight(motion, Ra, Ta0, Ta1) && ccx_straight(obstacle_motion, Ro, To0, To1);
		double *qa = NULL, *qo = NULL;
		if (straight)
		{
			/* (the bodies' points turned once: every vertex at q + T(t)) */
			qa = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.a.nv + 1));
			qo = (double*)qaws_internal_malloc(sizeof(double) * (3 * (size_t)c.o.nv + 1));
			if (!qa || !qo)
			{
				qaws_internal_free(qa);
				qaws_internal_free(qo);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			for (i = 0; i < c.a.nv; i++)
				for (k = 0; k < 3; k++)
					qa[3 * (size_t)i + k] = Ra[3 * k] * (double)c.a.x[3 * (size_t)i] + Ra[3 * k + 1] * (double)c.a.x[3 * (size_t)i + 1] + Ra[3 * k + 2] * (double)c.a.x[3 * (size_t)i + 2];
			for (i = 0; i < c.o.nv; i++)
				for (k = 0; k < 3; k++)
					qo[3 * (size_t)i + k] = Ro[3 * k] * (double)c.o.x[3 * (size_t)i] + Ro[3 * k + 1] * (double)c.o.x[3 * (size_t)i + 1] + Ro[3 * k + 2] * (double)c.o.x[3 * (size_t)i + 2];
		}
		for (i = 0; i < ntimes / 2; i++)
		{
			double start = times[2 * (size_t)i], end = times[2 * (size_t)i + 1];
			while (start < end)
			{
				double ta = 0, tb = 0;
				int r;
				if ((st = cc_at(&c.a, start, X0)) != QAWS_STATUS_OK || (st = cc_at(&c.a, end, X1)) != QAWS_STATUS_OK)
					break;
				r = cc_first(&c, start, end, X0, X1, dt, &ta, &tb, &f, &g, p, &st);
				if (r < 0)
					break;
				if (!r)
					break;
				if (!straight)
				{
					report->contact = 1;
					report->time = (qaws_scalar)ta;
					report->time_upper = (qaws_scalar)tb;
					qaws_internal_free(qa);
					qaws_internal_free(qo);
					goto found;
				}
				{
					/* the pairs touching over the interval widened by a step each way, their exact events */
					double wa = fmax(0, ta - dt), wb = fmin(1, tb + dt);
					qaws_exact_int num, den;
					unsigned int ef = 0, eg = 0;
					int rx;
					if ((st = cc_at(&c.a, wa, X0)) != QAWS_STATUS_OK || (st = cc_at(&c.a, wb, X1)) != QAWS_STATUS_OK || (st = cc_obstacle(&c, wa, wb)) != QAWS_STATUS_OK)
						break;
					c.collect = 1;
					c.npairs = 0;
					rx = cc_touch(&c, X0, X1, &f, &g, p);
					c.collect = 0;
					if (rx < 0)
					{
						st = QAWS_STATUS_ALLOCATION_FAILURE;
						break;
					}
					rx = ccx_first(&c, qa, Ta0, Ta1, qo, To0, To1, &num, &den, &ef, &eg);
					if (rx < 0)
					{
						/* (past the integers' range: the bracket as found) */
						report->contact = 1;
						report->time = (qaws_scalar)ta;
						report->time_upper = (qaws_scalar)tb;
						qaws_internal_free(qa);
						qaws_internal_free(qo);
						goto found;
					}
					if (rx > 0)
					{
						double t = qaws_exact_ratio_to_double(&num, &den);
						report->contact = 1;
						report->exact = 1;
						report->time = report->time_upper = (qaws_scalar)t;
						f = ef;
						g = eg;
						/* (the point: the two triangles' nearest at t) */
						if ((st = cc_at(&c.a, t, X0)) == QAWS_STATUS_OK && (st = cc_obstacle(&c, t, t)) == QAWS_STATUS_OK)
						{
							double best = HUGE_VAL, q2[3], Pa[3][3], Qb[3][3];
							unsigned int j;
							for (k = 0; k < 3; k++)
								for (j = 0; j < 3; j++)
								{
									Pa[k][j] = X0[3 * (size_t)c.a.t[3 * (size_t)f + k] + j];
									Qb[k][j] = c.OA[3 * (size_t)c.o.t[3 * (size_t)g + k] + j];
								}
							for (k = 0; k < 3; k++)
							{
								double d2;
								cc_closest(Pa[k], Qb[0], Qb[1], Qb[2], q2);
								d2 = (q2[0] - Pa[k][0]) * (q2[0] - Pa[k][0]) + (q2[1] - Pa[k][1]) * (q2[1] - Pa[k][1]) + (q2[2] - Pa[k][2]) * (q2[2] - Pa[k][2]);
								if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + Pa[k][j]); }
								cc_closest(Qb[k], Pa[0], Pa[1], Pa[2], q2);
								d2 = (q2[0] - Qb[k][0]) * (q2[0] - Qb[k][0]) + (q2[1] - Qb[k][1]) * (q2[1] - Qb[k][1]) + (q2[2] - Qb[k][2]) * (q2[2] - Qb[k][2]);
								if (d2 < best) { best = d2; for (j = 0; j < 3; j++) p[j] = 0.5 * (q2[j] + Qb[k][j]); }
							}
						}
						qaws_internal_free(qa);
						qaws_internal_free(qo);
						goto found;
					}
					/* (no event there: they only came within the tolerance) */
					start = tb;
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
		}
		qaws_internal_free(qa);
		qaws_internal_free(qo);
		goto done;
	}
found:
	report->face = f;
	report->obstacle_face = g;
	report->point.x = (qaws_scalar)p[0];
	report->point.y = (qaws_scalar)p[1];
	report->point.z = (qaws_scalar)p[2];
done:
	report->tests = (unsigned int)(c.tests > 0xFFFFFFFFull ? 0xFFFFFFFFull : c.tests);
	qaws_mesh_bvh_destroy(c.bvh);
	cc_body_free(&c.a);
	cc_body_free(&c.o);
	qaws_internal_free(X0); qaws_internal_free(X1); qaws_internal_free(Xm); qaws_internal_free(O0); qaws_internal_free(O1); qaws_internal_free(Om);
	qaws_internal_free(c.OA); qaws_internal_free(c.OB); qaws_internal_free(c.bpos); qaws_internal_free(c.bsim); qaws_internal_free(c.bel);
	qaws_internal_free(c.seen); qaws_internal_free(c.cand); qaws_internal_free(c.pairs); qaws_internal_free(times);
	return st;
}
