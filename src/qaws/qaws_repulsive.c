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
#include "internal/qaws_internal_memory.h"
#include "internal/qaws_internal_repulsive.h"
#include "internal/qaws_internal_repulsive_jet.h"
#include "qaws_alloc.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>


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

/* the pairs of elements (k vertices each, at idx), each with the term
   function, in double: positions, direction (may be NULL); acc (dim per
   vertex; NULL: values only) accumulates seed * gradient or the HVP. On
   failure acc's contents are unspecified. */
static qaws_status rp_run_d(unsigned int dim, unsigned int const* idx, unsigned int ne, unsigned int k,
	int (*pair)(void const*, rp_jet (*)[3], unsigned int, rp_jet*), void const* user, int skip_shared,
	double const* positions, double const* direction, rp_want want, double seed,
	double* out_v, double* out_t, double* out_t2, double* acc)
{
	unsigned int i, j, a, c, n;
	double v = 0, t = 0, t2 = 0;
	rp_jet x[6][3], term;
	if (!acc)
		want = RP_EVAL;
	n = want == RP_EVAL ? 0 : 2 * k * dim;
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
				return r == RP_DEGENERATE ? QAWS_STATUS_DEGENERATE_CURVE : QAWS_STATUS_NUMERICAL_FAILURE;
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
	return QAWS_STATUS_OK;
}

/* rp_run_d from qaws_scalar inputs: out_vec += the result, only on success */
static qaws_status rp_run(qaws_allocator const* alloc, unsigned int dim, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k,
	int (*pair)(void const*, rp_jet (*)[3], unsigned int, rp_jet*), void const* user, int skip_shared,
	qaws_scalar const* positions, qaws_scalar const* direction, rp_want want, double seed,
	double* out_v, double* out_t, double* out_t2, qaws_scalar* out_vec)
{
	unsigned int i, m = nv * dim;
	double *x, *d = NULL, *acc = NULL;
	qaws_status st;
	qaws_internal_rs_scope sc;
	if (!positions || (want == RP_HVP && !direction) || (want != RP_EVAL && !out_vec))
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_rs_enter(&sc, alloc);
	x = (double*)qaws_internal_scratch_zalloc(3 * (size_t)m + 1, sizeof(double));
	if (!x)
	{
		qaws_internal_rs_leave(&sc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (direction) d = x + m;
	if (want != RP_EVAL) acc = x + 2 * m;
	for (i = 0; i < m; i++)
	{
		x[i] = (double)positions[i];
		if (d) d[i] = (double)direction[i];
		if (acc) acc[i] = 0;
	}
	st = rp_run_d(dim, idx, ne, k, pair, user, skip_shared, x, d, want, seed, out_v, out_t, out_t2, acc);
	for (i = 0; st == QAWS_STATUS_OK && acc && i < m; i++)
		out_vec[i] = (qaws_scalar)((double)out_vec[i] + acc[i]);
	qaws_internal_free(x);
	qaws_internal_rs_leave(&sc);
	return st;
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
/*  Objects' allocators                                                     */
/* ------------------------------------------------------------------------ */

void qaws_internal_rs_keep(qaws_allocator* to, qaws_allocator const* a)
{
	memset(to, 0, sizeof(*to));
	if (!a || !a->alloc)
		return;
	to->alloc = a->alloc; to->dealloc = a->dealloc; to->user_data = a->user_data;
	if (a->alloc == qaws_allocator_extended)
	{
		/* (a legacy descriptor's other members are never read) */
		to->alloc_sized = a->alloc_sized; to->dealloc_sized = a->dealloc_sized; to->resize = a->resize;
	}
}

void qaws_internal_rs_enter(qaws_internal_rs_scope* sc, qaws_allocator const* a)
{
	sc->active = a && a->alloc;
	if (sc->active)
		qaws_exchange_thread_allocator(&sc->previous, a);
}

void qaws_internal_rs_leave(qaws_internal_rs_scope const* sc)
{
	if (sc->active)
		qaws_set_thread_allocator(&sc->previous);
}

/* ------------------------------------------------------------------------ */
/*  Curves                                                                  */
/* ------------------------------------------------------------------------ */

qaws_status qaws_repulsive_curves_create_ex(qaws_repulsive_curves** out_curves, qaws_repulsive_curves_desc const* desc, qaws_allocator const* allocator)
{
	qaws_internal_rs_scope sc;
	qaws_status st;
	qaws_internal_rs_enter(&sc, allocator);
	st = qaws_repulsive_curves_create(out_curves, desc);
	if (st == QAWS_STATUS_OK)
		qaws_internal_rs_keep(&(*out_curves)->alloc, allocator);
	qaws_internal_rs_leave(&sc);
	return st;
}

qaws_status qaws_repulsive_curves_create(qaws_repulsive_curves** out_curves, qaws_repulsive_curves_desc const* desc)
{
	qaws_repulsive_curves* c;
	double alpha, beta;
	unsigned int i;
	if (!desc || !out_curves || (desc->dimension != 2 && desc->dimension != 3) || (desc->edge_count && !desc->edges))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* 2 edge_count indices and vertex_count * dimension coordinates count in unsigned int */
	if (desc->edge_count > 0xFFFFFFFFu / 2 || desc->vertex_count > 0xFFFFFFFFu / 3)
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
	c = (qaws_repulsive_curves*)qaws_internal_zalloc(1, sizeof(*c));
	if (!c)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	c->edges = (unsigned int*)qaws_internal_malloc(qaws_internal_add2(qaws_internal_mul3(sizeof(unsigned int), 2, desc->edge_count), 1));
	if (!c->edges)
	{
		qaws_internal_free(c);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (desc->edge_count) memcpy(c->edges, desc->edges, sizeof(unsigned int) * 2 * desc->edge_count);
	c->dim = desc->dimension; c->nv = desc->vertex_count; c->ne = desc->edge_count;
	c->alpha = alpha; c->beta = beta; c->executor = desc->executor;
	*out_curves = c;
	return QAWS_STATUS_OK;
}

void qaws_repulsive_curves_destroy(qaws_repulsive_curves* curves)
{
	if (!curves)
		return;
	qaws_internal_free(curves->edges);
	qaws_internal_free(curves);
}

qaws_status qaws_repulsive_curves_eval(qaws_scalar* out_value, qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2, qaws_repulsive_curves const* curves, qaws_scalar const* positions, qaws_scalar const* direction)
{
	double v, t, t2;
	qaws_status st;
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(&curves->alloc, curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, direction,
		RP_EVAL, 0, &v, &t, &t2, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)v;
	if (out_tangent) *out_tangent = (qaws_scalar)t;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)t2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_curves_gradient(qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_curves const* curves,
	qaws_scalar const* positions, qaws_scalar seed)
{
	double v;
	qaws_status st;
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(&curves->alloc, curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, NULL,
		RP_GRAD, (double)seed, &v, NULL, NULL, gradient);
	if (st == QAWS_STATUS_OK && out_value) *out_value = (qaws_scalar)v;
	return st;
}

qaws_status qaws_repulsive_curves_hvp(qaws_scalar* out_hv, qaws_repulsive_curves const* curves, qaws_scalar const* positions,
	qaws_scalar const* direction)
{
	if (!curves)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return rp_run(&curves->alloc, curves->dim, curves->nv, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, positions, direction,
		RP_HVP, 0, NULL, NULL, NULL, out_hv);
}

/* ------------------------------------------------------------------------ */
/*  Meshes                                                                  */
/* ------------------------------------------------------------------------ */

qaws_status qaws_repulsive_mesh_create_ex(qaws_repulsive_mesh** out_mesh, qaws_repulsive_mesh_desc const* desc, qaws_allocator const* allocator)
{
	qaws_internal_rs_scope sc;
	qaws_status st;
	qaws_internal_rs_enter(&sc, allocator);
	st = qaws_repulsive_mesh_create(out_mesh, desc);
	if (st == QAWS_STATUS_OK)
		qaws_internal_rs_keep(&(*out_mesh)->alloc, allocator);
	qaws_internal_rs_leave(&sc);
	return st;
}

qaws_status qaws_repulsive_mesh_create(qaws_repulsive_mesh** out_mesh, qaws_repulsive_mesh_desc const* desc)
{
	qaws_repulsive_mesh* m;
	double p;
	unsigned int i;
	if (!desc || !out_mesh || (desc->triangle_count && !desc->triangles))
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* 3 triangle_count indices and 3 vertex_count coordinates count in unsigned int */
	if (desc->triangle_count > 0xFFFFFFFFu / 3 || desc->vertex_count > 0xFFFFFFFFu / 3)
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
	m = (qaws_repulsive_mesh*)qaws_internal_zalloc(1, sizeof(*m));
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->tris = (unsigned int*)qaws_internal_malloc(qaws_internal_add2(qaws_internal_mul3(sizeof(unsigned int), 3, desc->triangle_count), 1));
	if (!m->tris)
	{
		qaws_internal_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (desc->triangle_count) memcpy(m->tris, desc->triangles, sizeof(unsigned int) * 3 * desc->triangle_count);
	m->nv = desc->vertex_count; m->nf = desc->triangle_count; m->p = p; m->executor = desc->executor;
	*out_mesh = m;
	return QAWS_STATUS_OK;
}

void qaws_repulsive_mesh_destroy(qaws_repulsive_mesh* mesh)
{
	if (!mesh)
		return;
	qaws_internal_free(mesh->tris);
	qaws_internal_free(mesh);
}

qaws_status qaws_repulsive_mesh_eval(qaws_scalar* out_value, qaws_scalar* out_tangent,
	qaws_scalar* out_tangent2, qaws_repulsive_mesh const* mesh, qaws_scalar const* positions, qaws_scalar const* direction)
{
	double v, t, t2;
	qaws_status st;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(&mesh->alloc, 3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, direction, RP_EVAL, 0, &v, &t, &t2, NULL);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_value) *out_value = (qaws_scalar)v;
	if (out_tangent) *out_tangent = (qaws_scalar)t;
	if (out_tangent2) *out_tangent2 = (qaws_scalar)t2;
	return QAWS_STATUS_OK;
}

qaws_status qaws_repulsive_mesh_gradient(qaws_scalar* gradient, qaws_scalar* out_value, qaws_repulsive_mesh const* mesh,
	qaws_scalar const* positions, qaws_scalar seed)
{
	double v;
	qaws_status st;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = rp_run(&mesh->alloc, 3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, NULL, RP_GRAD, (double)seed, &v, NULL, NULL, gradient);
	if (st == QAWS_STATUS_OK && out_value) *out_value = (qaws_scalar)v;
	return st;
}

qaws_status qaws_repulsive_mesh_hvp(qaws_scalar* out_hv, qaws_repulsive_mesh const* mesh, qaws_scalar const* positions,
	qaws_scalar const* direction)
{
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return rp_run(&mesh->alloc, 3, mesh->nv, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, positions, direction, RP_HVP, 0, NULL, NULL, NULL, out_hv);
}

/* ------------------------------------------------------------------------ */
/*  For the solver                                                          */
/* ------------------------------------------------------------------------ */

qaws_status qaws_internal_repulsive_curves_energy(double* out_energy, double* grad, qaws_repulsive_curves const* curves, double const* x)
{
	return rp_run_d(curves->dim, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, x, NULL, grad ? RP_GRAD : RP_EVAL, 1,
		out_energy, NULL, NULL, grad);
}

qaws_status qaws_internal_repulsive_mesh_energy(double* out_energy, double* grad, qaws_repulsive_mesh const* mesh, double const* x)
{
	return rp_run_d(3, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, x, NULL, grad ? RP_GRAD : RP_EVAL, 1, out_energy, NULL, NULL, grad);
}

qaws_status qaws_internal_repulsive_curves_hvp(double* hv, qaws_repulsive_curves const* curves, double const* x, double const* v)
{
	return rp_run_d(curves->dim, curves->edges, curves->ne, 2, rp_curve_term, curves, 1, x, v, RP_HVP, 0, NULL, NULL, NULL, hv);
}

qaws_status qaws_internal_repulsive_mesh_hvp(double* hv, qaws_repulsive_mesh const* mesh, double const* x, double const* v)
{
	return rp_run_d(3, mesh->tris, mesh->nf, 3, rp_mesh_term, mesh, 0, x, v, RP_HVP, 0, NULL, NULL, NULL, hv);
}
