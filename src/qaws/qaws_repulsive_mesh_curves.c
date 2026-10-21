/*
 * Repulsive curves on a triangle mesh (qaws_repulsive_mesh_curves.h): the
 * unknowns are each vertex's move in its triangle's tangent plane (two per
 * vertex, in an orthonormal frame), zero at the start of every step; an
 * accepted step walks each vertex along the triangles and the frames are
 * taken again where it lands. The descent runs on the core of
 * qaws_repulsive_solve.c with the curves' metric pulled back through the
 * frames, as qaws_repulsive_surface.c does through a surface's derivatives.
 */

#include "qaws_repulsive_mesh_curves.h"
#include "qaws_mesh_reference.h"
#include "qaws_mesh_geodesic.h"
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RMC_REG 1e-10   /* the diagonal blocks' lift, relative (as on surfaces) */

typedef qaws_internal_rs_core rs_core;

struct qaws_repulsive_mesh_curves_solver
{
	rs_core core;                   /* (first) x: the moves in the frames, two per vertex */
	qaws_repulsive_curves const* c;
	qaws_mesh_reference* ref;       /* the flat triangles: walks and positions */
	qaws_mesh_geodesic* geo;        /* NULL: the chords only */
	qaws_repulsive_curves_solver* inner;   /* the curves' constraints in 3D (pins and barycenter aside) */
	unsigned int nt, nmv, kin, npin;
	unsigned int *tri, *pin;
	unsigned char* fixed;
	double* P;                      /* the mesh's positions */
	qaws_mesh_location *loc, *wl;   /* the vertices' locations; walked ones (scratch) */
	double *F;                      /* frames at loc: e1, e2 (6 per vertex) */
	double *x, *x1, *A, *gx, *Cx;
	double *Dv, *Dg;                /* the distance from vertex j at vertex i (n x n) and its gradient at i (3 each) */
	qaws_scalar *xs, *ds, *q;
};

/* ------------------------------------------------------------------------ */
/*  Locations, frames, walks                                                */
/* ------------------------------------------------------------------------ */

static void rmc_corners(qaws_repulsive_mesh_curves_solver const* s, unsigned int t, double const** p)
{
	p[0] = s->P + 3 * (size_t)s->tri[3 * (size_t)t];
	p[1] = s->P + 3 * (size_t)s->tri[3 * (size_t)t + 1];
	p[2] = s->P + 3 * (size_t)s->tri[3 * (size_t)t + 2];
}

static void rmc_cross(double* o, double const* a, double const* b)
{
	o[0] = a[1] * b[2] - a[2] * b[1]; o[1] = a[2] * b[0] - a[0] * b[2]; o[2] = a[0] * b[1] - a[1] * b[0];
}

/* the frames of the vertices' triangles; 0 for a triangle without area */
static int rmc_frames(qaws_repulsive_mesh_curves_solver* s)
{
	unsigned int i, k;
	for (i = 0; i < s->core.nv; i++)
	{
		double const* p[3];
		double u[3], w[3], n[3], lu, ln;
		double* f = s->F + 6 * i;
		rmc_corners(s, s->loc[i].patch, p);
		for (k = 0; k < 3; k++) { u[k] = p[1][k] - p[0][k]; w[k] = p[2][k] - p[0][k]; }
		rmc_cross(n, u, w);
		lu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
		ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (!(lu > 0 && ln > 0))
			return 0;
		for (k = 0; k < 3; k++) { f[k] = u[k] / lu; n[k] /= ln; }
		rmc_cross(f + 3, n, f);
	}
	return 1;
}

/* the locations walked from loc by q (two per vertex, in the frames) into wl, their positions into x */
static qaws_status rmc_walk(qaws_repulsive_mesh_curves_solver const* s, double const* q, double* x)
{
	qaws_repulsive_mesh_curves_solver* m = (qaws_repulsive_mesh_curves_solver*)s;
	unsigned int n = s->core.nv, i, k;
	qaws_status st;
	for (i = 0; i < n; i++)
	{
		qaws_scalar t[3];
		double const* f = s->F + 6 * i;
		for (k = 0; k < 3; k++) t[k] = (qaws_scalar)(f[k] * q[2 * i] + f[3 + k] * q[2 * i + 1]);
		if (q[2 * i] == 0 && q[2 * i + 1] == 0)
			m->wl[i] = s->loc[i];
		else if ((st = qaws_mesh_reference_walk(m->wl + i, s->ref, t, s->loc[i])) != QAWS_STATUS_OK)
			return st;
	}
	if ((st = qaws_mesh_reference_evaluate(m->xs, NULL, NULL, s->ref, s->wl, n)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < 3 * n; i++) x[i] = (double)s->xs[i];
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  The energy                                                              */
/* ------------------------------------------------------------------------ */

/* Dv, Dg from the distance fields of the locations wl */
static qaws_status rmc_distances(qaws_repulsive_mesh_curves_solver const* s)
{
	qaws_repulsive_mesh_curves_solver* m = (qaws_repulsive_mesh_curves_solver*)s;
	unsigned int n = s->core.nv, i, j, k;
	qaws_status st;
	for (j = 0; j < n; j++)
	{
		if ((st = qaws_mesh_geodesic_distance(m->ds, s->geo, s->wl + j, 1)) != QAWS_STATUS_OK)
			return st;
		for (i = 0; i < n; i++)
		{
			/* the field's linear interpolation in i's triangle and its gradient there */
			unsigned int const* c = s->tri + 3 * (size_t)s->wl[i].patch;
			double const* p[3];
			double a = (double)s->wl[i].a, b = (double)s->wl[i].b, d[3], N[3], u[3], w[3], e[3], g[3] = { 0, 0, 0 }, A2;
			double* out = m->Dg + 3 * ((size_t)j * n + i);
			for (k = 0; k < 3; k++) d[k] = (double)s->ds[c[k]];
			m->Dv[(size_t)j * n + i] = (1 - a - b) * d[0] + a * d[1] + b * d[2];
			rmc_corners(s, s->wl[i].patch, p);
			for (k = 0; k < 3; k++) { u[k] = p[1][k] - p[0][k]; w[k] = p[2][k] - p[0][k]; }
			rmc_cross(N, u, w);
			A2 = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
			if (A2 > 0)
			{
				unsigned int r;
				for (k = 0; k < 3; k++) N[k] /= A2;
				/* grad phi_r = N x (the edge opposite r, counterclockwise) / 2A */
				for (r = 0; r < 3; r++)
				{
					double const *pa = p[(r + 1) % 3], *pb = p[(r + 2) % 3];
					double ed[3];
					for (k = 0; k < 3; k++) ed[k] = pb[k] - pa[k];
					rmc_cross(e, N, ed);
					for (k = 0; k < 3; k++) g[k] += d[r] * e[k] / A2;
				}
			}
			for (k = 0; k < 3; k++) out[k] = g[k];
		}
	}
	return QAWS_STATUS_OK;
}

/* E at positions x (locations wl) and, if g, g += dE/dx (3 per vertex) */
static qaws_status rmc_sum(qaws_repulsive_mesh_curves_solver const* s, double const* x, double* E, double* g)
{
	qaws_repulsive_curves const* c = s->c;
	unsigned int n = c->nv, I, J, k;
	double al = c->alpha, be = c->beta, sum = 0;
	qaws_status st;
	if (s->geo && (st = rmc_distances(s)) != QAWS_STATUS_OK)
		return st;
	for (I = 0; I < c->ne; I++)
	{
		unsigned int i0 = c->edges[2 * I], i1 = c->edges[2 * I + 1];
		double T[3], lI = 0;
		for (k = 0; k < 3; k++) { T[k] = x[3 * i1 + k] - x[3 * i0 + k]; lI += T[k] * T[k]; }
		lI = sqrt(lI);
		if (!(lI > 0))
			return QAWS_STATUS_DEGENERATE_CURVE;
		for (k = 0; k < 3; k++) T[k] /= lI;
		for (J = 0; J < c->ne; J++)
		{
			unsigned int j0 = c->edges[2 * J], j1 = c->edges[2 * J + 1], a2, b2;
			unsigned int const ii[2] = { i0, i1 }, jj[2] = { j0, j1 };
			double TJ[3], lJ = 0, w;
			if (J == I || j0 == i0 || j0 == i1 || j1 == i0 || j1 == i1)
				continue;
			for (k = 0; k < 3; k++) { TJ[k] = x[3 * j1 + k] - x[3 * j0 + k]; lJ += TJ[k] * TJ[k]; }
			lJ = sqrt(lJ);
			if (!(lJ > 0))
				return QAWS_STATUS_DEGENERATE_CURVE;
			for (k = 0; k < 3; k++) TJ[k] /= lJ;
			w = 0.25 * lI * lJ;
			for (a2 = 0; a2 < 2; a2++)
				for (b2 = 0; b2 < 2; b2++)
				{
					unsigned int i = ii[a2], j = jj[b2];
					double v[3], e = 0, r, a = 0, P, kk, kr, ka, ri[3], rj[3];
					for (k = 0; k < 3; k++) { v[k] = x[3 * i + k] - x[3 * j + k]; e += v[k] * v[k]; a += T[k] * v[k]; }
					e = sqrt(e);
					if (!(e > 0))
						return QAWS_STATUS_NUMERICAL_FAILURE;
					r = e;
					for (k = 0; k < 3; k++) { ri[k] = v[k] / e; rj[k] = -v[k] / e; }
					if (s->geo)
					{
						double D = 0.5 * (s->Dv[(size_t)j * n + i] + s->Dv[(size_t)i * n + j]);
						if (!(D < HUGE_VAL))
							continue;   /* (another part of the mesh: no interaction) */
						if (D > e)
						{
							r = D;
							for (k = 0; k < 3; k++) { ri[k] = s->Dg[3 * ((size_t)j * n + i) + k]; rj[k] = s->Dg[3 * ((size_t)i * n + j) + k]; }
						}
					}
					P = r * r - a * a;
					if (P < 0) P = 0;
					kk = pow(P, 0.5 * al) * pow(r, -be);
					sum += w * kk;
					if (!g)
						continue;
					kr = al * pow(P, 0.5 * al - 1) * pow(r, 1 - be) - be * kk / r;
					ka = -al * pow(P, 0.5 * al - 1) * a * pow(r, -be);
					for (k = 0; k < 3; k++)
					{
						double tv = (v[k] - a * T[k]) / lI;   /* d a / d x_i1 through T */
						g[3 * i + k] += w * (kr * ri[k] + ka * T[k]);
						g[3 * j + k] += w * (kr * rj[k] - ka * T[k]);
						g[3 * i1 + k] += w * ka * tv + kk * 0.25 * lJ * T[k];
						g[3 * i0 + k] -= w * ka * tv + kk * 0.25 * lJ * T[k];
						g[3 * j1 + k] += kk * 0.25 * lI * TJ[k];
						g[3 * j0 + k] -= kk * 0.25 * lI * TJ[k];
					}
				}
		}
	}
	*E = sum;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  The core's callbacks                                                    */
/* ------------------------------------------------------------------------ */

static qaws_status rmc_energy(rs_core const* core, double const* q, double* E, double* grad)
{
	qaws_repulsive_mesh_curves_solver const* s = (qaws_repulsive_mesh_curves_solver const*)core;
	unsigned int n = core->nv, i, a, k;
	qaws_status st;
	if ((st = rmc_walk(s, q, s->x)) != QAWS_STATUS_OK)
		return st;
	if (grad) memset(s->gx, 0, sizeof(double) * 3 * n);
	if ((st = rmc_sum(s, s->x, E, grad ? s->gx : NULL)) != QAWS_STATUS_OK || !grad)
		return st;
	for (i = 0; i < n; i++)
		for (a = 0; a < 2; a++)
			for (k = 0; k < 3; k++) grad[2 * i + a] += s->F[6 * i + 3 * a + k] * s->gx[3 * i + k];
	return QAWS_STATUS_OK;
}

/* M = blocks A_ij F_i^T F_j, the diagonal ones lifted by RMC_REG */
static qaws_status rmc_metric(rs_core const* core, double const* q, double* M)
{
	qaws_repulsive_mesh_curves_solver const* s = (qaws_repulsive_mesh_curves_solver const*)core;
	unsigned int n = core->nv, N = core->N, i, j, a, b, k;
	qaws_status st;
	if ((st = rmc_walk(s, q, s->x)) != QAWS_STATUS_OK || (st = qaws_internal_repulsive_curves_metric(s->A, s->c, s->x)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (a = 0; a < 2; a++)
				for (b = 0; b < 2; b++)
				{
					double d = 0, w = s->A[(size_t)i * n + j] * (i == j ? 1 + RMC_REG : 1);
					for (k = 0; k < 3; k++) d += s->F[6 * i + 3 * a + k] * s->F[6 * j + 3 * b + k];
					M[(size_t)(2 * i + a) * N + 2 * j + b] = w * d;
				}
	return QAWS_STATUS_OK;
}

/* rows: the pins (their moves held at 0), then the curves' rows times the frames */
static void rmc_constraints(rs_core const* core, double const* q, double* phi, double* C)
{
	qaws_repulsive_mesh_curves_solver const* s = (qaws_repulsive_mesh_curves_solver const*)core;
	unsigned int n = core->nv, N = core->N, r = 0, i, a, k, row;
	rs_core const* in = (rs_core const*)s->inner;
	if (C) memset(C, 0, sizeof(double) * core->k * N);
	for (i = 0; i < s->npin; i++)
		for (a = 0; a < 2; a++, r++)
		{
			phi[r] = q[2 * s->pin[i] + a];
			if (C) C[(size_t)r * N + 2 * s->pin[i] + a] = 1;
		}
	if (!s->kin)
		return;
	if (rmc_walk(s, q, s->x) != QAWS_STATUS_OK)
	{
		for (row = 0; row < s->kin; row++) phi[r + row] = HUGE_VAL;
		return;
	}
	in->ops->constraints(in, s->x, phi + r, C ? s->Cx : NULL);
	if (C)
		for (row = 0; row < s->kin; row++)
			for (i = 0; i < n; i++)
				for (a = 0; a < 2; a++)
				{
					double v = 0;
					for (k = 0; k < 3; k++) v += s->Cx[(size_t)row * 3 * n + 3 * i + k] * s->F[6 * i + 3 * a + k];
					C[(size_t)(r + row) * N + 2 * i + a] = v;
				}
}

static int rmc_path_ok(rs_core const* core, double const* q0, double const* q1)
{
	qaws_repulsive_mesh_curves_solver const* s = (qaws_repulsive_mesh_curves_solver const*)core;
	if (rmc_walk(s, q0, s->x) != QAWS_STATUS_OK)
		return 0;
	memcpy(s->x1, s->x, sizeof(double) * 3 * core->nv);
	return rmc_walk(s, q1, s->x) == QAWS_STATUS_OK && qaws_internal_repulsive_curves_path_clear(s->c, core->clearance, s->x1, s->x);
}

static qaws_internal_rs_ops const rmc_ops = { rmc_energy, rmc_metric, rmc_constraints, rmc_path_ok, NULL, NULL, NULL, NULL, NULL };

/* ------------------------------------------------------------------------ */
/*  The solver                                                              */
/* ------------------------------------------------------------------------ */

void qaws_repulsive_mesh_curves_solver_destroy(qaws_repulsive_mesh_curves_solver* s)
{
	if (!s)
		return;
	qaws_repulsive_curves_solver_destroy(s->inner);
	qaws_mesh_reference_destroy(s->ref);
	qaws_mesh_geodesic_destroy(s->geo);
	qaws_internal_free(s->core.x); qaws_internal_free(s->tri); qaws_internal_free(s->pin); qaws_internal_free(s->fixed); qaws_internal_free(s->P);
	qaws_internal_free(s->loc); qaws_internal_free(s->wl); qaws_internal_free(s->F); qaws_internal_free(s->x); qaws_internal_free(s->x1);
	qaws_internal_free(s->A); qaws_internal_free(s->gx); qaws_internal_free(s->Cx); qaws_internal_free(s->Dv); qaws_internal_free(s->Dg);
	qaws_internal_free(s->xs); qaws_internal_free(s->ds); qaws_internal_free(s->q);
	qaws_internal_free(s);
}

static int rmc_location_ok(qaws_repulsive_mesh_curves_solver const* s, qaws_mesh_location const* l)
{
	return l->patch < s->nt && (double)l->a >= 0 && (double)l->b >= 0 && (double)l->a + (double)l->b <= 1 + 1e-6;
}

qaws_status qaws_repulsive_mesh_curves_solver_create(qaws_repulsive_mesh_curves_solver** out_solver, qaws_repulsive_curves const* c,
	qaws_mesh_location const* locations, qaws_repulsive_mesh_curves_desc const* md, qaws_repulsive_curves_solve_desc const* desc)
{
	qaws_repulsive_mesh_curves_solver* s;
	qaws_repulsive_curves_solve_desc d, di;
	qaws_mesh_reference_desc rd;
	qaws_mesh_geodesic_desc gd;
	rs_core* k;
	rs_core const* in;
	qaws_scalar* ps = NULL;
	unsigned int n, i;
	double E;
	qaws_status st;
	if (!out_solver || !c || !locations || !md || !md->mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_solver = NULL;
	if (c->dim != 3)
		return QAWS_STATUS_INVALID_DIMENSION;
	memset(&d, 0, sizeof(d));
	if (desc) d = *desc;
	if (d.theta != 0 || d.hold_barycenter != 0 || (d.pinned_count && !d.pinned) || !(md->geodesic_time >= 0))
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = c->nv;
	for (i = 0; i < d.pinned_count; i++) if (d.pinned[i] >= n) return QAWS_STATUS_INVALID_ARGUMENT;
	s = (qaws_repulsive_mesh_curves_solver*)qaws_internal_zalloc(1, sizeof(*s));
	if (!s)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s->c = c;
	memset(&rd, 0, sizeof(rd));
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	if ((st = qaws_mesh_reference_create(&s->ref, md->mesh, &rd)) != QAWS_STATUS_OK)
		goto fail;
	if (!md->euclidean)
	{
		memset(&gd, 0, sizeof(gd));
		gd.time = md->geodesic_time;
		if ((st = qaws_mesh_geodesic_create(&s->geo, md->mesh, &gd)) != QAWS_STATUS_OK)
			goto fail;
	}
	s->nt = qaws_mesh_3d_get_triangle_count(md->mesh);
	s->nmv = qaws_mesh_3d_get_vertex_count(md->mesh);
	k = &s->core;
	k->ops = &rmc_ops; k->dim = 2; k->nv = n; k->N = 2 * n; k->full = 1;
	k->x = (double*)qaws_internal_zalloc(2 * (size_t)n + 1, sizeof(double));
	s->tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), 3 * (size_t)s->nt + 1));
	s->pin = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (d.pinned_count + 1)));
	s->fixed = (unsigned char*)qaws_internal_zalloc(n + 1, 1);
	s->P = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 3 * (size_t)s->nmv + 1));
	s->loc = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)n + 1));
	s->wl = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)n + 1));
	s->F = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 6 * (size_t)n + 1));
	s->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 3 * (size_t)n + 1));
	s->x1 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 3 * (size_t)n + 1));
	s->A = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)n * n + 1));
	s->gx = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 3 * (size_t)n + 1));
	s->xs = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 3 * (size_t)n + 1));
	s->q = (qaws_scalar*)qaws_internal_zalloc(2 * (size_t)n + 1, sizeof(qaws_scalar));
	ps = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), 3 * (size_t)s->nmv + 1));
	if (s->geo)
	{
		s->Dv = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)n * n + 1));
		s->Dg = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), 3 * (size_t)n * n + 1));
		s->ds = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)s->nmv + 1));
	}
	if (!k->x || !s->tri || !s->pin || !s->fixed || !s->P || !s->loc || !s->wl || !s->F || !s->x || !s->x1 || !s->A || !s->gx || !s->xs || !s->q || !ps ||
		(s->geo && (!s->Dv || !s->Dg || !s->ds)))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	qaws_mesh_3d_get_triangles(s->tri, md->mesh);
	qaws_mesh_3d_get_positions(ps, md->mesh);
	for (i = 0; i < 3 * s->nmv; i++) s->P[i] = (double)ps[i];
	for (i = 0; i < n; i++)
	{
		if (!rmc_location_ok(s, locations + i))
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto fail;
		}
		s->loc[i] = locations[i];
	}
	if (!rmc_frames(s))
	{
		st = QAWS_STATUS_INVALID_ARGUMENT;
		goto fail;
	}
	if ((st = rmc_walk(s, k->x, s->x)) != QAWS_STATUS_OK)
		goto fail;
	/* the 3D rows (targets at the start) from a curve solver on the positions */
	di = d;
	di.pinned = NULL; di.pinned_count = 0; di.hold_barycenter = -1; di.theta = 0;
	if ((st = qaws_repulsive_curves_solver_create(&s->inner, c, s->xs, &di)) != QAWS_STATUS_OK)
		goto fail;
	in = (rs_core const*)s->inner;
	s->kin = in->k;
	s->Cx = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)s->kin * 3 * n + 1)));
	if (!s->Cx)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	s->npin = d.pinned_count;
	for (i = 0; i < s->npin; i++)
	{
		s->pin[i] = d.pinned[i];
		s->fixed[d.pinned[i]] = 1;
	}
	k->fixed = s->fixed;
	k->k = 2 * s->npin + s->kin;
	k->checked = in->checked;
	k->lmean = in->lmean; k->clearance = in->clearance; k->tol = in->tol; k->ctol = in->ctol; k->h = in->h;
	k->armijo = in->armijo; k->backtrack = in->backtrack; k->trials = in->trials;
	k->gmres_max = in->gmres_max; k->gmres_restart = in->gmres_restart; k->gmres_tol = in->gmres_tol;
	k->progress = in->progress;
	k->rep.reason = QAWS_REPULSIVE_RUNNING;
	k->rep.stationarity = 1;
	if ((st = rmc_sum(s, s->x, &E, NULL)) != QAWS_STATUS_OK)
		goto fail;
	k->rep.initial_energy = k->rep.energy = (qaws_scalar)E;
	qaws_internal_free(ps);
	*out_solver = s;
	return QAWS_STATUS_OK;
fail:
	qaws_internal_free(ps);
	qaws_repulsive_mesh_curves_solver_destroy(s);
	return st;
}

static qaws_status rmc_step_run(qaws_repulsive_report* out_report, qaws_mesh_location* locations, qaws_repulsive_mesh_curves_solver* s)
{
	unsigned int n = s->core.nv, i, moved = 0;
	qaws_status st;
	for (i = 0; i < n; i++)
	{
		qaws_mesh_location const* l = locations + i;
		if (l->patch != s->loc[i].patch || l->a != s->loc[i].a || l->b != s->loc[i].b)
		{
			if (!rmc_location_ok(s, l))
				return QAWS_STATUS_INVALID_ARGUMENT;
			s->loc[i] = *l;
			moved = 1;
		}
	}
	if (moved && !rmc_frames(s))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* each step from the locations themselves: the moves start at 0 */
	memset(s->core.x, 0, sizeof(double) * 2 * n);
	memset(s->q, 0, sizeof(qaws_scalar) * 2 * n);
	st = qaws_internal_rs_step(out_report, s->q, &s->core);
	for (i = 0; i < 2 * n && s->core.x[i] == 0; i++) {}
	if (i < 2 * n)
	{
		/* accepted: walk there and take the frames again */
		if ((st = rmc_walk(s, s->core.x, s->x)) != QAWS_STATUS_OK)
			return st;
		memcpy(s->loc, s->wl, sizeof(qaws_mesh_location) * n);
		memset(s->core.x, 0, sizeof(double) * 2 * n);
		if (!rmc_frames(s))
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	memcpy(locations, s->loc, sizeof(qaws_mesh_location) * n);
	return st;
}

qaws_status qaws_repulsive_mesh_curves_solver_step(qaws_repulsive_report* out_report, qaws_mesh_location* locations, qaws_repulsive_mesh_curves_solver* s)
{
	qaws_internal_rs_scope sc;
	qaws_status st;
	if (!s || !locations)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_rs_enter(&sc, &s->core.alloc);
	st = rmc_step_run(out_report, locations, s);
	qaws_internal_rs_leave(&sc);
	return st;
}

qaws_status qaws_repulsive_mesh_curves_solver_energy(qaws_scalar* out_energy, qaws_scalar* out_gradient, qaws_repulsive_mesh_curves_solver const* s)
{
	qaws_internal_rs_scope sc;
	unsigned int n, i;
	double E = 0;
	qaws_status st;
	if (!s || !out_energy)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = s->core.nv;
	qaws_internal_rs_enter(&sc, &s->core.alloc);
	memset(s->q, 0, sizeof(qaws_scalar) * 2 * n);
	memset(s->gx, 0, sizeof(double) * 3 * n);
	memset(s->core.x, 0, sizeof(double) * 2 * n);
	st = rmc_walk(s, s->core.x, s->x);
	if (st == QAWS_STATUS_OK)
		st = rmc_sum(s, s->x, &E, out_gradient ? s->gx : NULL);
	if (st == QAWS_STATUS_OK)
	{
		*out_energy = (qaws_scalar)E;
		for (i = 0; out_gradient && i < 3 * n; i++) out_gradient[i] = (qaws_scalar)s->gx[i];
	}
	qaws_internal_rs_leave(&sc);
	return st;
}

qaws_status qaws_repulsive_mesh_curves_solver_get_positions(qaws_scalar* out_positions, qaws_repulsive_mesh_curves_solver const* s)
{
	return s && out_positions ? qaws_mesh_reference_evaluate(out_positions, NULL, NULL, s->ref, s->loc, s->core.nv) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_mesh_curves_solve(qaws_repulsive_report* out_report, qaws_mesh_location* locations, qaws_repulsive_curves const* c,
	qaws_repulsive_mesh_curves_desc const* md, qaws_repulsive_curves_solve_desc const* desc)
{
	qaws_repulsive_mesh_curves_solver* s = NULL;
	qaws_repulsive_report rep;
	unsigned int max = desc && desc->max_iterations ? desc->max_iterations : 100, it;
	qaws_status st = qaws_repulsive_mesh_curves_solver_create_ex(&s, c, locations, md, desc, c ? &c->alloc : NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(&rep, 0, sizeof(rep));
	for (it = 0; it < max && st == QAWS_STATUS_OK; it++)
	{
		unsigned int before = s->core.rep.iterations;
		st = qaws_repulsive_mesh_curves_solver_step(&rep, locations, s);
		if (st == QAWS_STATUS_OK && rep.iterations > before && s->core.progress.on_step &&
			qaws_repulsive_mesh_curves_solver_get_positions(s->xs, s) == QAWS_STATUS_OK)
			s->core.progress.on_step(&rep, s->xs, s->core.progress.user_data);
		if (rep.reason != QAWS_REPULSIVE_RUNNING)
			break;
	}
	if (st == QAWS_STATUS_OK && rep.reason == QAWS_REPULSIVE_RUNNING)
		rep.reason = QAWS_REPULSIVE_ITERATION_LIMIT;
	if (out_report) *out_report = rep;
	qaws_repulsive_mesh_curves_solver_destroy(s);
	return st;
}

qaws_status qaws_repulsive_mesh_curves_solver_create_ex(qaws_repulsive_mesh_curves_solver** out_solver, qaws_repulsive_curves const* c,
	qaws_mesh_location const* locations, qaws_repulsive_mesh_curves_desc const* md, qaws_repulsive_curves_solve_desc const* desc,
	qaws_allocator const* allocator)
{
	qaws_internal_rs_scope sc;
	qaws_status st;
	qaws_internal_rs_enter(&sc, allocator);
	st = qaws_repulsive_mesh_curves_solver_create(out_solver, c, locations, md, desc);
	if (st == QAWS_STATUS_OK)
		qaws_internal_rs_keep(&(*out_solver)->core.alloc, allocator);
	qaws_internal_rs_leave(&sc);
	return st;
}
