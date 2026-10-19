/*
 * Curves on a parametric surface (qaws_repulsive.h, R5): the vertices'
 * surface parameters are the unknowns, x_i = S(u_i, v_i), the energy the
 * curves' own in the embedding. The descent runs on the core of
 * qaws_repulsive_solve.c with the metric pulled back through S (an N x N
 * metric, N = 2 |V|); the curves' length and direction constraints come
 * from an x-space curve solver (its rows times J), the pins hold
 * parameters.
 */

#include "qaws_repulsive.h"
#include "internal/qaws_internal_memory.h"
#include "qaws_surface.h"
#include "qaws_surface_trim.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_solve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RSF_REG 1e-10   /* the diagonal blocks' lift, relative */
#if QAWS_SCALAR_IS_FLOAT
#define RSF_SINGULAR 1e-5   /* |S_u x S_v| over |S_u|^2 + |S_v|^2: a pole, at the precision of the evaluation */
#else
#define RSF_SINGULAR 1e-10
#endif

typedef qaws_internal_rs_core rs_core;

/* ------------------------------------------------------------------------ */
/*  The map                                                                 */
/* ------------------------------------------------------------------------ */

/* S at one vertex: x (3) and, if J, [S_u S_v] (J[c] = S_u, J[3 + c] = S_v),
   if H, [S_uu S_uv S_vv]; 0 outside the domain or at a singular point */
static int rsf_at(qaws_repulsive_surface_desc const* sd, double u, double v, double* x, double* J, double* H)
{
	qaws_range ru = qaws_surface_get_u_range(sd->surface), rv = qaws_surface_get_v_range(sd->surface);
	double ulo = (double)ru.min_value, uhi = (double)ru.max_value, vlo = (double)rv.min_value, vhi = (double)rv.max_value, du2 = 0, dv2 = 0, cr[3];
	qaws_surface_eval_result r;
	unsigned int c;
	if (sd->periodic_u) { u = ulo + fmod(u - ulo, uhi - ulo); if (u < ulo) u += uhi - ulo; }
	else if (!(u >= ulo && u <= uhi)) return 0;
	if (sd->periodic_v) { v = vlo + fmod(v - vlo, vhi - vlo); if (v < vlo) v += vhi - vlo; }
	else if (!(v >= vlo && v <= vhi)) return 0;
	if (qaws_surface_get_kind(sd->surface) == QAWS_SURFACE_KIND_TRIMMED)
	{
		int in = 0;
		if (qaws_surface_trim_contains(sd->surface, (qaws_scalar)u, (qaws_scalar)v, &in) != QAWS_STATUS_OK || !in)
			return 0;
	}
	if (qaws_surface_evaluate(sd->surface, (qaws_scalar)u, (qaws_scalar)v,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | (H ? QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV : 0),
		&r) != QAWS_STATUS_OK)
		return 0;
	x[0] = (double)r.position.x; x[1] = (double)r.position.y; x[2] = (double)r.position.z;
	{
		double a[3], b[3];
		a[0] = (double)r.du.x; a[1] = (double)r.du.y; a[2] = (double)r.du.z;
		b[0] = (double)r.dv.x; b[1] = (double)r.dv.y; b[2] = (double)r.dv.z;
		cr[0] = a[1] * b[2] - a[2] * b[1]; cr[1] = a[2] * b[0] - a[0] * b[2]; cr[2] = a[0] * b[1] - a[1] * b[0];
		for (c = 0; c < 3; c++) { du2 += a[c] * a[c]; dv2 += b[c] * b[c]; }
		/* (singular: the tangent plane degenerates, as at a sphere's pole) */
		if (!(sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]) > RSF_SINGULAR * (du2 + dv2)))
			return 0;
		if (J)
			for (c = 0; c < 3; c++) { J[c] = a[c]; J[3 + c] = b[c]; }
	}
	if (H)
	{
		H[0] = (double)r.duu.x; H[1] = (double)r.duu.y; H[2] = (double)r.duu.z;
		H[3] = (double)r.duv.x; H[4] = (double)r.duv.y; H[5] = (double)r.duv.z;
		H[6] = (double)r.dvv.x; H[7] = (double)r.dvv.y; H[8] = (double)r.dvv.z;
	}
	return 1;
}

/* all n vertices: x (3n), J (6n), H (9n) as rsf_at */
static int rsf_map(qaws_repulsive_surface_desc const* sd, double const* q, unsigned int n, double* x, double* J, double* H)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		if (!rsf_at(sd, q[2 * i], q[2 * i + 1], x + 3 * i, J ? J + 6 * i : NULL, H ? H + 9 * i : NULL))
			return 0;
	return 1;
}

/* ------------------------------------------------------------------------ */
/*  The energy in uv                                                        */
/* ------------------------------------------------------------------------ */

typedef struct rsf_eval
{
	double *x, *J, *H;
	qaws_scalar *xs, *dxs, *g;
} rsf_eval;

static void rsf_eval_free(rsf_eval* e)
{
	qaws_internal_free(e->x); qaws_internal_free(e->J); qaws_internal_free(e->H); qaws_internal_free(e->xs); qaws_internal_free(e->dxs); qaws_internal_free(e->g);
}

/* the positions (and their velocity along duv) at uv; the gradient g at x */
static qaws_status rsf_eval_make(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv, qaws_scalar const* duv,
	rsf_eval* e, double* E)
{
	unsigned int n = c->nv, i, k;
	double* q = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (2 * n + 1)));
	qaws_scalar Es = 0;
	qaws_status st;
	memset(e, 0, sizeof(*e));
	e->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * n + 1)));
	e->J = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (6 * n + 1)));
	e->H = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (9 * n + 1)));
	e->xs = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (3 * n + 1)));
	e->dxs = (qaws_scalar*)qaws_internal_zalloc(3 * (size_t)n + 1, sizeof(qaws_scalar));
	e->g = (qaws_scalar*)qaws_internal_zalloc(3 * (size_t)n + 1, sizeof(qaws_scalar));
	if (!q || !e->x || !e->J || !e->H || !e->xs || !e->dxs || !e->g)
	{
		qaws_internal_free(q);
		rsf_eval_free(e);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < 2 * n; i++) q[i] = (double)uv[i];
	if (!rsf_map(sd, q, n, e->x, e->J, e->H))
	{
		qaws_internal_free(q);
		rsf_eval_free(e);
		return QAWS_STATUS_OUT_OF_RANGE;
	}
	qaws_internal_free(q);
	for (i = 0; i < 3 * n; i++) e->xs[i] = (qaws_scalar)e->x[i];
	if (duv)
		for (i = 0; i < n; i++)
			for (k = 0; k < 3; k++)
				e->dxs[3 * i + k] = (qaws_scalar)(e->J[6 * i + k] * (double)duv[2 * i] + e->J[6 * i + 3 + k] * (double)duv[2 * i + 1]);
	st = qaws_repulsive_curves_gradient(c, e->xs, 1, e->g, &Es);
	if (st != QAWS_STATUS_OK)
		rsf_eval_free(e);
	*E = (double)Es;
	return st;
}

static qaws_status rsf_check(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv)
{
	return !c || !sd || !sd->surface || !uv ? QAWS_STATUS_INVALID_ARGUMENT : c->dim != 3 ? QAWS_STATUS_INVALID_DIMENSION : QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_surface_positions(qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv, unsigned int count, qaws_scalar* out_positions)
{
	unsigned int i;
	if (!sd || !sd->surface || !uv || !out_positions)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		double x[3];
		if (!rsf_at(sd, (double)uv[2 * i], (double)uv[2 * i + 1], x, NULL, NULL))
			return QAWS_STATUS_OUT_OF_RANGE;
		out_positions[3 * i] = (qaws_scalar)x[0]; out_positions[3 * i + 1] = (qaws_scalar)x[1]; out_positions[3 * i + 2] = (qaws_scalar)x[2];
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_curves_on_surface_eval(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv,
	qaws_scalar const* duv, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2)
{
	rsf_eval e;
	double E;
	qaws_scalar v = 0, t = 0, t2 = 0;
	qaws_status st = rsf_check(c, sd, uv);
	unsigned int i, k;
	if (st != QAWS_STATUS_OK || (st = rsf_eval_make(c, sd, uv, duv, &e, &E)) != QAWS_STATUS_OK)
		return st;
	if (duv)
	{
		/* t2 = d2E[dx, dx] + g . S''[duv, duv] */
		double acc = 0;
		st = qaws_repulsive_curves_eval(c, e.xs, e.dxs, &v, &t, &t2);
		for (i = 0; i < c->nv; i++)
		{
			double a = (double)duv[2 * i], b = (double)duv[2 * i + 1];
			for (k = 0; k < 3; k++) acc += (double)e.g[3 * i + k] * (e.H[9 * i + k] * a * a + 2 * e.H[9 * i + 3 + k] * a * b + e.H[9 * i + 6 + k] * b * b);
		}
		t2 = (qaws_scalar)((double)t2 + acc);
	}
	if (st == QAWS_STATUS_OK)
	{
		if (out_value) *out_value = (qaws_scalar)E;
		if (out_tangent) *out_tangent = t;
		if (out_tangent2) *out_tangent2 = t2;
	}
	rsf_eval_free(&e);
	return st;
}

qaws_status qaws_repulsive_curves_on_surface_gradient(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv,
	qaws_scalar seed, qaws_scalar* gradient, qaws_scalar* out_value)
{
	rsf_eval e;
	double E;
	unsigned int i, k, a;
	qaws_status st = rsf_check(c, sd, uv);
	if (st != QAWS_STATUS_OK || !gradient || (st = rsf_eval_make(c, sd, uv, NULL, &e, &E)) != QAWS_STATUS_OK)
		return st != QAWS_STATUS_OK ? st : QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < c->nv; i++)
		for (a = 0; a < 2; a++)
		{
			double s = 0;
			for (k = 0; k < 3; k++) s += e.J[6 * i + 3 * a + k] * (double)e.g[3 * i + k];
			gradient[2 * i + a] = (qaws_scalar)((double)gradient[2 * i + a] + (double)seed * s);
		}
	if (out_value) *out_value = (qaws_scalar)E;
	rsf_eval_free(&e);
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_curves_on_surface_hvp(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv,
	qaws_scalar const* duv, qaws_scalar* out_hv)
{
	rsf_eval e;
	double E;
	unsigned int i, k, a, n;
	qaws_scalar* hx;
	qaws_status st = rsf_check(c, sd, uv);
	if (st != QAWS_STATUS_OK || !duv || !out_hv)
		return st != QAWS_STATUS_OK ? st : QAWS_STATUS_INVALID_ARGUMENT;
	n = c->nv;
	if ((st = rsf_eval_make(c, sd, uv, duv, &e, &E)) != QAWS_STATUS_OK)
		return st;
	hx = (qaws_scalar*)qaws_internal_scratch_zalloc(3 * (size_t)n + 1, sizeof(qaws_scalar));
	if (!hx)
	{
		rsf_eval_free(&e);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* J^T H_x (J duv) + sum_c g_c S_c''[duv, .] */
	st = qaws_repulsive_curves_hvp(c, e.xs, e.dxs, hx);
	for (i = 0; st == QAWS_STATUS_OK && i < n; i++)
	{
		double du = (double)duv[2 * i], dv = (double)duv[2 * i + 1];
		for (a = 0; a < 2; a++)
		{
			double s = 0;
			for (k = 0; k < 3; k++)
			{
				/* S''[duv, e_a]: a = u: S_uu du + S_uv dv; a = v: S_uv du + S_vv dv */
				double sec = a == 0 ? e.H[9 * i + k] * du + e.H[9 * i + 3 + k] * dv : e.H[9 * i + 3 + k] * du + e.H[9 * i + 6 + k] * dv;
				s += e.J[6 * i + 3 * a + k] * (double)hx[3 * i + k] + (double)e.g[3 * i + k] * sec;
			}
			out_hv[2 * i + a] = (qaws_scalar)((double)out_hv[2 * i + a] + s);
		}
	}
	qaws_internal_free(hx);
	rsf_eval_free(&e);
	return st;
}

/* ------------------------------------------------------------------------ */
/*  The solver                                                              */
/* ------------------------------------------------------------------------ */

struct qaws_repulsive_surface_curves_solver
{
	rs_core core;                   /* (first) x: the parameters, two per vertex */
	qaws_repulsive_curves const* c;
	qaws_repulsive_surface_desc sd;
	qaws_repulsive_curves_solver* inner;   /* the curves' constraints in 3D (pins and barycenter aside) */
	unsigned int *pin, npin, kin;
	unsigned char* fixed;
	double *pin_target, *x, *x1, *J, *Cx, *A, *gx;
};

static qaws_status rsf_energy(rs_core const* core, double const* q, double* E, double* grad)
{
	qaws_repulsive_surface_curves_solver const* s = (qaws_repulsive_surface_curves_solver const*)core;
	unsigned int n = core->nv, i, a, k;
	qaws_status st;
	if (!rsf_map(&s->sd, q, n, s->x, grad ? s->J : NULL, NULL))
		return QAWS_STATUS_OUT_OF_RANGE;   /* (a step out of the domain: rejected) */
	if (grad) memset(s->gx, 0, sizeof(double) * 3 * n);
	if ((st = qaws_internal_repulsive_curves_energy(s->c, s->x, E, grad ? s->gx : NULL)) != QAWS_STATUS_OK || !grad)
		return st;
	for (i = 0; i < n; i++)
		for (a = 0; a < 2; a++)
			for (k = 0; k < 3; k++) grad[2 * i + a] += s->J[6 * i + 3 * a + k] * s->gx[3 * i + k];
	return QAWS_STATUS_OK;
}

/* M = blocks A_ij J_i^T J_j, the diagonal ones lifted by RSF_REG A_ii J_i^T J_i */
static qaws_status rsf_metric(rs_core const* core, double const* q, double* M)
{
	qaws_repulsive_surface_curves_solver const* s = (qaws_repulsive_surface_curves_solver const*)core;
	unsigned int n = core->nv, N = core->N, i, j, a, b, k;
	qaws_status st;
	if (!rsf_map(&s->sd, q, n, s->x, s->J, NULL))
		return QAWS_STATUS_OUT_OF_RANGE;
	if ((st = qaws_internal_repulsive_curves_metric(s->c, s->x, s->A)) != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (a = 0; a < 2; a++)
				for (b = 0; b < 2; b++)
				{
					double g = 0, w = s->A[(size_t)i * n + j] * (i == j ? 1 + RSF_REG : 1);
					for (k = 0; k < 3; k++) g += s->J[6 * i + 3 * a + k] * s->J[6 * j + 3 * b + k];
					M[(size_t)(2 * i + a) * N + 2 * j + b] = w * g;
				}
	return QAWS_STATUS_OK;
}

/* rows: the pins (parameters), then the curves' rows times J */
static void rsf_constraints(rs_core const* core, double const* q, double* phi, double* C)
{
	qaws_repulsive_surface_curves_solver const* s = (qaws_repulsive_surface_curves_solver const*)core;
	unsigned int n = core->nv, N = core->N, r = 0, i, a, k, row;
	rs_core const* in = (rs_core const*)s->inner;
	if (C) memset(C, 0, sizeof(double) * core->k * N);
	for (i = 0; i < s->npin; i++)
		for (a = 0; a < 2; a++, r++)
		{
			phi[r] = q[2 * s->pin[i] + a] - s->pin_target[2 * i + a];
			if (C) C[(size_t)r * N + 2 * s->pin[i] + a] = 1;
		}
	if (!s->kin)
		return;
	if (!rsf_map(&s->sd, q, n, s->x, s->J, NULL))
	{
		/* (outside the domain: unmeetable, so the projection fails) */
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
					for (k = 0; k < 3; k++) v += s->Cx[(size_t)row * 3 * n + 3 * i + k] * s->J[6 * i + 3 * a + k];
					C[(size_t)(r + row) * N + 2 * i + a] = v;
				}
}

static int rsf_path_ok(rs_core const* core, double const* q0, double const* q1)
{
	qaws_repulsive_surface_curves_solver const* s = (qaws_repulsive_surface_curves_solver const*)core;
	return rsf_map(&s->sd, q0, core->nv, s->x, NULL, NULL) && rsf_map(&s->sd, q1, core->nv, s->x1, NULL, NULL) &&
		qaws_internal_repulsive_curves_path_clear(s->c, core->clearance, s->x, s->x1);
}

/* |J_v d_v|: the vertex's speed in 3D */
static double rsf_move(rs_core const* core, double const* q, double const* d, unsigned int v)
{
	qaws_repulsive_surface_curves_solver const* s = (qaws_repulsive_surface_curves_solver const*)core;
	double x[3], J[6], m = 0;
	unsigned int k;
	if (!rsf_at(&s->sd, q[2 * v], q[2 * v + 1], x, J, NULL))
		return HUGE_VAL;
	for (k = 0; k < 3; k++)
	{
		double w = J[k] * d[2 * v] + J[3 + k] * d[2 * v + 1];
		m += w * w;
	}
	return sqrt(m);
}

static qaws_internal_rs_ops const rsf_ops = { rsf_energy, rsf_metric, rsf_constraints, rsf_path_ok, NULL, rsf_move, NULL, NULL, NULL };

void qaws_repulsive_surface_curves_solver_destroy(qaws_repulsive_surface_curves_solver* s)
{
	if (!s)
		return;
	qaws_repulsive_curves_solver_destroy(s->inner);
	qaws_internal_free(s->core.x); qaws_internal_free(s->pin); qaws_internal_free(s->fixed); qaws_internal_free(s->pin_target); qaws_internal_free(s->x); qaws_internal_free(s->x1); qaws_internal_free(s->J); qaws_internal_free(s->Cx); qaws_internal_free(s->A); qaws_internal_free(s->gx);
	qaws_internal_free(s);
}

qaws_status qaws_repulsive_surface_curves_solver_create(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd,
	qaws_repulsive_curves_solve_desc const* desc, qaws_scalar const* uv, qaws_repulsive_surface_curves_solver** out_solver)
{
	qaws_repulsive_surface_curves_solver* s;
	qaws_repulsive_curves_solve_desc d, di;
	rs_core* k;
	rs_core const* in;
	qaws_scalar* xs;
	unsigned int n, i;
	double E;
	qaws_status st = rsf_check(c, sd, uv);
	if (st != QAWS_STATUS_OK || !out_solver)
		return st != QAWS_STATUS_OK ? st : QAWS_STATUS_INVALID_ARGUMENT;
	*out_solver = NULL;
	memset(&d, 0, sizeof(d));
	if (desc) d = *desc;
	if (d.theta != 0 || d.hold_barycenter != 0 || (d.pinned_count && !d.pinned))
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (dense mode only; no barycenter on a surface) */
	n = c->nv;
	for (i = 0; i < d.pinned_count; i++) if (d.pinned[i] >= n) return QAWS_STATUS_INVALID_ARGUMENT;
	s = (qaws_repulsive_surface_curves_solver*)qaws_internal_zalloc(1, sizeof(*s));
	xs = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (3 * n + 1)));
	if (!s || !xs)
	{
		qaws_internal_free(s); qaws_internal_free(xs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	s->c = c; s->sd = *sd;
	k = &s->core;
	k->ops = &rsf_ops; k->dim = 2; k->nv = n; k->N = 2 * n; k->full = 1;
	k->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (2 * n + 1)));
	s->pin = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (d.pinned_count + 1)));
	s->fixed = (unsigned char*)qaws_internal_zalloc(n + 1, 1);
	s->pin_target = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (2 * d.pinned_count + 1)));
	s->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * n + 1))); s->x1 = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * n + 1)));
	s->J = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (6 * n + 1)));
	s->A = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)n * n + 1)));
	s->gx = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * n + 1)));
	if (!k->x || !s->pin || !s->fixed || !s->pin_target || !s->x || !s->x1 || !s->J || !s->A || !s->gx)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto fail;
	}
	for (i = 0; i < 2 * n; i++) k->x[i] = (double)uv[i];
	if (!rsf_map(sd, k->x, n, s->x, NULL, NULL))
	{
		st = QAWS_STATUS_OUT_OF_RANGE;   /* (a vertex outside the domain, or at a singular point) */
		goto fail;
	}
	/* the 3D rows (targets at the start) from a curve solver on the positions */
	for (i = 0; i < 3 * n; i++) xs[i] = (qaws_scalar)s->x[i];
	di = d;
	di.pinned = NULL; di.pinned_count = 0; di.hold_barycenter = -1; di.theta = 0;
	if ((st = qaws_repulsive_curves_solver_create(c, &di, xs, &s->inner)) != QAWS_STATUS_OK)
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
		s->pin_target[2 * i] = k->x[2 * d.pinned[i]]; s->pin_target[2 * i + 1] = k->x[2 * d.pinned[i] + 1];
	}
	k->fixed = s->fixed;
	k->k = 2 * s->npin + s->kin;
	k->checked = in->checked;
	k->lmean = in->lmean; k->clearance = in->clearance; k->tol = in->tol; k->ctol = in->ctol; k->h = in->h;
	k->progress = in->progress;
	k->rep.reason = QAWS_REPULSIVE_RUNNING;
	k->rep.stationarity = 1;
	if ((st = qaws_internal_repulsive_curves_energy(c, s->x, &E, NULL)) != QAWS_STATUS_OK)
		goto fail;
	k->rep.initial_energy = k->rep.energy = (qaws_scalar)E;
	qaws_internal_free(xs);
	*out_solver = s;
	return QAWS_STATUS_OK;
fail:
	qaws_internal_free(xs);
	qaws_repulsive_surface_curves_solver_destroy(s);
	return st;
}

qaws_status qaws_repulsive_surface_curves_solver_direction(qaws_repulsive_surface_curves_solver const* s, qaws_scalar const* uv,
	qaws_scalar* out_direction)
{
	return s ? qaws_internal_rs_direction(&s->core, uv, out_direction) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_surface_curves_solver_step(qaws_repulsive_surface_curves_solver* s, qaws_scalar* uv, qaws_repulsive_report* out_report)
{
	return s ? qaws_internal_rs_step(&s->core, uv, out_report) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_surface_curves_solver_objective(qaws_repulsive_surface_curves_solver const* s, qaws_scalar const* uv, qaws_scalar* out_energy,
	qaws_scalar* out_gradient)
{
	return s ? qaws_internal_rs_objective(&s->core, uv, out_energy, out_gradient) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_repulsive_surface_curves_solve(qaws_repulsive_curves const* c, qaws_repulsive_surface_desc const* sd,
	qaws_repulsive_curves_solve_desc const* desc, qaws_scalar* uv, qaws_repulsive_report* out_report)
{
	qaws_repulsive_surface_curves_solver* s = NULL;
	unsigned int max = desc && desc->max_iterations ? desc->max_iterations : 100;
	qaws_status st = qaws_repulsive_surface_curves_solver_create(c, sd, desc, uv, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_internal_rs_solve(&s->core, max, uv, out_report);
	qaws_repulsive_surface_curves_solver_destroy(s);
	return st;
}
