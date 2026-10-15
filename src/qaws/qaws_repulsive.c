/*
 * Tangent-point (repulsive) energies of polygonal curves and triangle
 * meshes: the direct all-pairs evaluator and its exact derivatives
 * (qaws_repulsive.h; the models in docs/repulsive_plan.md).
 *
 * Each pair term is a function of the positions of 4 vertices (curves) or
 * 6 (meshes). It is evaluated once in a second-order forward jet: the
 * value, its gradient in the term's own variables, its first and second
 * derivatives along the direction, and the derivative of the gradient
 * along the direction (the Hessian times the direction). The jet's rules
 * are hyper-dual numbers' with a vector first slot. The terms are scattered
 * to their vertices.
 */

#include "qaws_repulsive.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RP_MAXN 18      /* a pair's variables: 6 vertices of 3 coordinates */

/* ------------------------------------------------------------------------ */
/*  Second-order jets                                                       */
/* ------------------------------------------------------------------------ */

/* f and its derivatives: g along each variable, t and t2 along the
   direction, h the derivative of g along the direction */
typedef struct rp_jet
{
	double v, t, t2;
	double g[RP_MAXN], h[RP_MAXN];
} rp_jet;

static void rp_const(rp_jet* r, double c, unsigned int n)
{
	unsigned int k;
	r->v = c; r->t = r->t2 = 0;
	for (k = 0; k < n; k++) r->g[k] = r->h[k] = 0;
}

/* variable number slot (n = 0: no gradient wanted), moving at rate dx */
static void rp_var(rp_jet* r, double x, double dx, unsigned int slot, unsigned int n)
{
	rp_const(r, x, n);
	r->t = dx;
	if (slot < n) r->g[slot] = 1;
}

static void rp_add(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	unsigned int k;
	r->v = a->v + b->v; r->t = a->t + b->t; r->t2 = a->t2 + b->t2;
	for (k = 0; k < n; k++) { r->g[k] = a->g[k] + b->g[k]; r->h[k] = a->h[k] + b->h[k]; }
}

static void rp_sub(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	unsigned int k;
	r->v = a->v - b->v; r->t = a->t - b->t; r->t2 = a->t2 - b->t2;
	for (k = 0; k < n; k++) { r->g[k] = a->g[k] - b->g[k]; r->h[k] = a->h[k] - b->h[k]; }
}

static void rp_scale(rp_jet* r, rp_jet const* a, double c, unsigned int n)
{
	unsigned int k;
	r->v = c * a->v; r->t = c * a->t; r->t2 = c * a->t2;
	for (k = 0; k < n; k++) { r->g[k] = c * a->g[k]; r->h[k] = c * a->h[k]; }
}

static void rp_mul(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	rp_jet q;
	unsigned int k;
	q.v = a->v * b->v;
	q.t = a->t * b->v + a->v * b->t;
	q.t2 = a->t2 * b->v + 2 * a->t * b->t + a->v * b->t2;
	for (k = 0; k < n; k++)
	{
		q.g[k] = a->g[k] * b->v + a->v * b->g[k];
		q.h[k] = a->h[k] * b->v + a->g[k] * b->t + a->t * b->g[k] + a->v * b->h[k];
	}
	*r = q;
}

/* u^s for u >= 0. At u = 0 (u a square, so its derivatives vanish there
   too) the limits of the smooth formulas: u itself for s = 1, all zero for
   s > 1 (|w|^(2s) has zero value, gradient and Hessian at w = 0). 0 when
   undefined (u < 0, or u = 0 with s < 1). */
static int rp_pow(rp_jet* r, rp_jet const* u, double s, unsigned int n)
{
	double f0, f1, f2;
	unsigned int k;
	rp_jet q;
	if (u->v <= 0)
	{
		if (u->v < 0 || s < 1)
			return 0;
		if (s == 1)
			*r = *u;
		else
			rp_const(r, 0, n);
		return 1;
	}
	f0 = pow(u->v, s);
	f1 = s * f0 / u->v;
	f2 = (s - 1) * f1 / u->v;
	q.v = f0;
	q.t = f1 * u->t;
	q.t2 = f1 * u->t2 + f2 * u->t * u->t;
	for (k = 0; k < n; k++)
	{
		q.g[k] = f1 * u->g[k];
		q.h[k] = f1 * u->h[k] + f2 * u->g[k] * u->t;
	}
	*r = q;
	return 1;
}

static void rp_dot(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int dim, unsigned int n)
{
	rp_jet p;
	unsigned int c;
	rp_mul(r, &a[0], &b[0], n);
	for (c = 1; c < dim; c++)
	{
		rp_mul(&p, &a[c], &b[c], n);
		rp_add(r, r, &p, n);
	}
}

static void rp_cross3(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int n)
{
	rp_jet p, q;
	unsigned int c;
	for (c = 0; c < 3; c++)
	{
		rp_mul(&p, &a[(c + 1) % 3], &b[(c + 2) % 3], n);
		rp_mul(&q, &a[(c + 2) % 3], &b[(c + 1) % 3], n);
		rp_sub(&r[c], &p, &q, n);
	}
}

/* |a x b|^2 (2D: the scalar cross product squared) */
static void rp_cross2(rp_jet* r, rp_jet const* a, rp_jet const* b, unsigned int dim, unsigned int n)
{
	rp_jet c[3];
	if (dim == 3)
	{
		rp_cross3(c, a, b, n);
		rp_dot(r, c, c, 3, n);
	}
	else
	{
		rp_jet p, q;
		rp_mul(&p, &a[0], &b[1], n);
		rp_mul(&q, &a[1], &b[0], n);
		rp_sub(&c[0], &p, &q, n);
		rp_mul(r, &c[0], &c[0], n);
	}
}

/* ------------------------------------------------------------------------ */
/*  Pair terms                                                              */
/* ------------------------------------------------------------------------ */

#define RP_OK 0
#define RP_DEGENERATE 1
#define RP_COINCIDENT 2

/* both ordered terms of edges I = (x[0], x[1]) and J = (x[2], x[3]):
   l_P l_Q / 4 sum k(x_i, x_j, T_P) = l_P^(1 - alpha) l_Q / 4
   sum |e_P x (x_i - x_j)|^alpha |x_i - x_j|^-beta */
static int rp_curve_pair(rp_jet (*x)[3], unsigned int dim, double alpha, double beta, unsigned int n, rp_jet* out)
{
	rp_jet e[2][3], L2[2], w, a, b, sum, d[3], dd, cr, k;
	unsigned int s, c, i, j;
	for (s = 0; s < 2; s++)
	{
		for (c = 0; c < dim; c++) rp_sub(&e[s][c], &x[2 * s + 1][c], &x[2 * s][c], n);
		for (c = dim; c < 3; c++) rp_const(&e[s][c], 0, n);   /* (2D: unread) */
		rp_dot(&L2[s], e[s], e[s], dim, n);
		if (!(L2[s].v > 0))
			return RP_DEGENERATE;
	}
	rp_const(out, 0, n);
	for (s = 0; s < 2; s++)
	{
		/* P = edge s, Q = the other */
		rp_pow(&a, &L2[s], 0.5 * (1 - alpha), n);
		rp_pow(&b, &L2[1 - s], 0.5, n);
		rp_mul(&w, &a, &b, n);
		rp_const(&sum, 0, n);
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
			{
				rp_jet* p = x[2 * s + i];
				rp_jet* q = x[2 * (1 - s) + j];
				for (c = 0; c < dim; c++) rp_sub(&d[c], &p[c], &q[c], n);
				for (c = dim; c < 3; c++) rp_const(&d[c], 0, n);
				rp_dot(&dd, d, d, dim, n);
				if (!(dd.v > 0))
					return RP_COINCIDENT;
				rp_cross2(&cr, e[s], d, dim, n);
				rp_pow(&a, &cr, 0.5 * alpha, n);
				rp_pow(&b, &dd, -0.5 * beta, n);
				rp_mul(&k, &a, &b, n);
				rp_add(&sum, &sum, &k, n);
			}
		rp_mul(&k, &w, &sum, n);
		rp_scale(&k, &k, 0.25, n);
		rp_add(out, out, &k, n);
	}
	return RP_OK;
}

/* both ordered terms of faces S = (x[0..2]) and T = (x[3..5]):
   a_S a_T |N_S . d|^p |d|^-2p = |m_S|^(1 - p) |m_T| / 4 (m_S . d)^p |d|^-2p,
   m the face's (unnormalized) normal, d the barycenters' difference */
static int rp_mesh_pair(rp_jet (*x)[3], double p, unsigned int n, rp_jet* out)
{
	rp_jet m[2][3], M2[2], X[2][3], u[3], v[3], d[3], dd, w, a, b, k, pr;
	unsigned int s, c;
	for (s = 0; s < 2; s++)
	{
		for (c = 0; c < 3; c++)
		{
			rp_sub(&u[c], &x[3 * s + 1][c], &x[3 * s][c], n);
			rp_sub(&v[c], &x[3 * s + 2][c], &x[3 * s][c], n);
			rp_add(&X[s][c], &x[3 * s][c], &x[3 * s + 1][c], n);
			rp_add(&X[s][c], &X[s][c], &x[3 * s + 2][c], n);
			rp_scale(&X[s][c], &X[s][c], 1.0 / 3, n);
		}
		rp_cross3(m[s], u, v, n);
		rp_dot(&M2[s], m[s], m[s], 3, n);
		if (!(M2[s].v > 0))
			return RP_DEGENERATE;
	}
	for (c = 0; c < 3; c++) rp_sub(&d[c], &X[0][c], &X[1][c], n);
	rp_dot(&dd, d, d, 3, n);
	if (!(dd.v > 0))
		return RP_COINCIDENT;
	rp_pow(&b, &dd, -p, n);
	rp_const(out, 0, n);
	for (s = 0; s < 2; s++)
	{
		rp_pow(&a, &M2[s], 0.5 * (1 - p), n);
		rp_pow(&w, &M2[1 - s], 0.5, n);
		rp_mul(&w, &w, &a, n);
		rp_dot(&pr, m[s], d, 3, n);
		rp_mul(&pr, &pr, &pr, n);
		rp_pow(&k, &pr, 0.5 * p, n);
		rp_mul(&k, &k, &b, n);
		rp_mul(&k, &k, &w, n);
		rp_scale(&k, &k, 0.25, n);
		rp_add(out, out, &k, n);
	}
	return RP_OK;
}

/* ------------------------------------------------------------------------ */
/*  The sums                                                                */
/* ------------------------------------------------------------------------ */

typedef enum rp_want
{
	RP_EVAL = 0,        /* value, tangent, tangent2 */
	RP_GRAD,            /* value, gradient */
	RP_HVP              /* Hessian times the direction */
} rp_want;

struct qaws_repulsive_curves
{
	unsigned int dim, nv, ne;
	unsigned int* edges;
	double alpha, beta;
};

struct qaws_repulsive_mesh
{
	unsigned int nv, nf;
	unsigned int* tris;
	double p;
};

/* the pairs of elements (k vertices each, at idx), each with the term
   function; results in double, written to the outputs only on success */
static qaws_status rp_run(unsigned int dim, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k,
	int (*pair)(void const*, rp_jet (*)[3], unsigned int, rp_jet*), void const* user, int skip_shared,
	qaws_scalar const* positions, qaws_scalar const* direction, rp_want want, double seed,
	double* out_v, double* out_t, double* out_t2, qaws_scalar* out_vec)
{
	unsigned int i, j, a, c, n = want == RP_EVAL ? 0 : 2 * k * dim;
	double v = 0, t = 0, t2 = 0, *acc = NULL;
	rp_jet x[6][3], term;
	if (!positions || ((want == RP_HVP) && !direction))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (want != RP_EVAL)
	{
		if (!out_vec)
			return QAWS_STATUS_INVALID_ARGUMENT;
		acc = (double*)calloc((size_t)nv * dim + 1, sizeof(double));
		if (!acc)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < ne; i++)
		for (j = i + 1; j < ne; j++)
		{
			unsigned int const* ei = idx + k * i;
			unsigned int const* ej = idx + k * j;
			int r;
			if (skip_shared)
			{
				/* (curves: edges sharing a vertex do not interact) */
				int shared = 0;
				for (a = 0; a < k; a++)
					for (c = 0; c < k; c++)
						shared |= ei[a] == ej[c];
				if (shared)
					continue;
			}
			for (a = 0; a < 2 * k; a++)
			{
				unsigned int vtx = a < k ? ei[a] : ej[a - k];
				for (c = 0; c < dim; c++)
					rp_var(&x[a][c], (double)positions[dim * vtx + c], direction ? (double)direction[dim * vtx + c] : 0, a * dim + c, n);
			}
			r = pair(user, x, n, &term);
			if (r != RP_OK)
			{
				free(acc);
				return r == RP_DEGENERATE ? QAWS_STATUS_DEGENERATE_CURVE : QAWS_STATUS_NUMERICAL_FAILURE;
			}
			v += term.v; t += term.t; t2 += term.t2;
			for (a = 0; acc && a < 2 * k; a++)
			{
				unsigned int vtx = a < k ? ei[a] : ej[a - k];
				for (c = 0; c < dim; c++)
					acc[dim * vtx + c] += want == RP_GRAD ? seed * term.g[a * dim + c] : term.h[a * dim + c];
			}
		}
	if (out_v) *out_v = v;
	if (out_t) *out_t = direction ? t : 0;
	if (out_t2) *out_t2 = direction ? t2 : 0;
	for (i = 0; acc && i < nv * dim; i++)
		out_vec[i] = (qaws_scalar)((double)out_vec[i] + acc[i]);
	free(acc);
	return QAWS_STATUS_OK;
}

static int rp_curve_term(void const* user, rp_jet (*x)[3], unsigned int n, rp_jet* out)
{
	qaws_repulsive_curves const* c = (qaws_repulsive_curves const*)user;
	return rp_curve_pair(x, c->dim, c->alpha, c->beta, n, out);
}

static int rp_mesh_term(void const* user, rp_jet (*x)[3], unsigned int n, rp_jet* out)
{
	return rp_mesh_pair(x, ((qaws_repulsive_mesh const*)user)->p, n, out);
}

/* ------------------------------------------------------------------------ */
/*  Curves                                                                  */
/* ------------------------------------------------------------------------ */

qaws_status qaws_repulsive_curves_create(qaws_repulsive_curves_desc const* desc, qaws_repulsive_curves** out_curves)
{
	qaws_repulsive_curves* c;
	double alpha, beta;
	unsigned int i;
	if (!desc || !out_curves || (desc->dimension != 2 && desc->dimension != 3) || (desc->edge_count && !desc->edges))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curves = NULL;
	alpha = desc->alpha == 0 && desc->beta == 0 ? 3 : (double)desc->alpha;
	beta = desc->alpha == 0 && desc->beta == 0 ? 6 : (double)desc->beta;
	/* the energy's range (alpha > 1, beta in [alpha + 2, 2 alpha + 1]), and
	   alpha >= 2 for second derivatives where |T x d| = 0 */
	if (!(alpha >= 2 && beta >= alpha + 2 && beta <= 2 * alpha + 1))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->edge_count; i++)
		if (desc->edges[2 * i] >= desc->vertex_count || desc->edges[2 * i + 1] >= desc->vertex_count || desc->edges[2 * i] == desc->edges[2 * i + 1])
			return QAWS_STATUS_INVALID_ARGUMENT;
	c = (qaws_repulsive_curves*)calloc(1, sizeof(*c));
	if (!c)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	c->edges = (unsigned int*)malloc(sizeof(unsigned int) * 2 * desc->edge_count + 1);
	if (!c->edges)
	{
		free(c);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (desc->edge_count) memcpy(c->edges, desc->edges, sizeof(unsigned int) * 2 * desc->edge_count);
	c->dim = desc->dimension; c->nv = desc->vertex_count; c->ne = desc->edge_count;
	c->alpha = alpha; c->beta = beta;
	*out_curves = c;
	return QAWS_STATUS_OK;
}

void qaws_repulsive_curves_destroy(qaws_repulsive_curves* curves)
{
	if (!curves)
		return;
	free(curves->edges);
	free(curves);
}

qaws_status qaws_repulsive_curves_eval(qaws_repulsive_curves const* curves, qaws_scalar const* positions,
	qaws_scalar const* direction, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2)
{
	double v, t, t2;
	qaws_status st;
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, direction,
		RP_EVAL, 0, &v, &t, &t2, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)v;
	if (out_tangent) *out_tangent = (qaws_scalar)t;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)t2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_curves_gradient(qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value)
{
	double v;
	qaws_status st;
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, NULL,
		RP_GRAD, (double)seed, &v, NULL, NULL, gradient);
	if (st == QAWS_STATUS_OK && out_value) *out_value = (qaws_scalar)v;
	return st;
}

qaws_status qaws_repulsive_curves_hvp(qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar const* direction,
	qaws_scalar* out_hv)
{
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return rp_run(curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, direction,
		RP_HVP, 0, NULL, NULL, NULL, out_hv);
}

/* ------------------------------------------------------------------------ */
/*  Meshes                                                                  */
/* ------------------------------------------------------------------------ */

qaws_status qaws_repulsive_mesh_create(qaws_repulsive_mesh_desc const* desc, qaws_repulsive_mesh** out_mesh)
{
	qaws_repulsive_mesh* m;
	double p;
	unsigned int i;
	if (!desc || !out_mesh || (desc->triangle_count && !desc->triangles))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	p = desc->p == 0 ? 6 : (double)desc->p;
	if (!(p > 4))
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (p > 2n: the energy truly repulsive) */
	for (i = 0; i < desc->triangle_count; i++)
	{
		unsigned int const* t = desc->triangles + 3 * i;
		if (t[0] >= desc->vertex_count || t[1] >= desc->vertex_count || t[2] >= desc->vertex_count || t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	m = (qaws_repulsive_mesh*)calloc(1, sizeof(*m));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->tris = (unsigned int*)malloc(sizeof(unsigned int) * 3 * desc->triangle_count + 1);
	if (!m->tris)
	{
		free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (desc->triangle_count) memcpy(m->tris, desc->triangles, sizeof(unsigned int) * 3 * desc->triangle_count);
	m->nv = desc->vertex_count; m->nf = desc->triangle_count; m->p = p;
	*out_mesh = m;
	return QAWS_STATUS_OK;
}

void qaws_repulsive_mesh_destroy(qaws_repulsive_mesh* mesh)
{
	if (!mesh)
		return;
	free(mesh->tris);
	free(mesh);
}

qaws_status qaws_repulsive_mesh_eval(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions,
	qaws_scalar const* direction, qaws_scalar* out_value, qaws_scalar* out_tangent, qaws_scalar* out_tangent2)
{
	double v, t, t2;
	qaws_status st;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, direction, RP_EVAL, 0, &v, &t, &t2, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)v;
	if (out_tangent) *out_tangent = (qaws_scalar)t;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)t2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_mesh_gradient(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar seed,
	qaws_scalar* gradient, qaws_scalar* out_value)
{
	double v;
	qaws_status st;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, NULL, RP_GRAD, (double)seed, &v, NULL, NULL, gradient);
	if (st == QAWS_STATUS_OK && out_value) *out_value = (qaws_scalar)v;
	return st;
}

qaws_status qaws_repulsive_mesh_hvp(qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar const* direction,
	qaws_scalar* out_hv)
{
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return rp_run(3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, direction, RP_HVP, 0, NULL, NULL, NULL, out_hv);
}
