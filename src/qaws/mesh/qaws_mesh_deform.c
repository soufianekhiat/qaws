/*
 * Deformation by handles (qaws_mesh_deform.h): the rest shape's Laplacian
 * (qaws_mesh_smooth.c's sparse helpers), the linear editing systems, and
 * ARAP's local-global rounds on one factored system; a deformer keeps the
 * factors for many targets, qaws_mesh_3d_deform uses one for one.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_deform.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the rotation nearest a (row-major 3 x 3) in the sense of the largest trace of R^T a, from the quaternion q (updated: the
   start of the next call), by Muller, Bender, Chentanez, Macklin (2016): robust where a is singular */
static void df_rotation(double const* a, double* q, double* R)
{
	unsigned int it, k;
	for (it = 0; it < 20; it++)
	{
		double w = q[0], x = q[1], y = q[2], z = q[3], om[3] = { 0, 0, 0 }, den = 0, l, s, c, h[4];
		/* R from q */
		R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - w * z); R[2] = 2 * (x * z + w * y);
		R[3] = 2 * (x * y + w * z); R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - w * x);
		R[6] = 2 * (x * z - w * y); R[7] = 2 * (y * z + w * x); R[8] = 1 - 2 * (x * x + y * y);
		/* sum over the columns of (R's column) x (a's column), over the sum of their dots */
		for (k = 0; k < 3; k++)
		{
			double r[3] = { R[k], R[3 + k], R[6 + k] }, b[3] = { a[k], a[3 + k], a[6 + k] };
			om[0] += r[1] * b[2] - r[2] * b[1];
			om[1] += r[2] * b[0] - r[0] * b[2];
			om[2] += r[0] * b[1] - r[1] * b[0];
			den += r[0] * b[0] + r[1] * b[1] + r[2] * b[2];
		}
		den = fabs(den) + 1e-12;
		om[0] /= den;
		om[1] /= den;
		om[2] /= den;
		l = sqrt(om[0] * om[0] + om[1] * om[1] + om[2] * om[2]);
		if (l < 1e-12)
			break;
		s = sin(0.5 * l) / l;
		c = cos(0.5 * l);
		/* q = (turn by omega) q */
		h[0] = c * w - s * (om[0] * x + om[1] * y + om[2] * z);
		h[1] = c * x + s * (om[0] * w + om[1] * z - om[2] * y);
		h[2] = c * y + s * (om[1] * w + om[2] * x - om[0] * z);
		h[3] = c * z + s * (om[2] * w + om[0] * y - om[1] * x);
		l = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3]);
		for (k = 0; k < 4; k++)
			q[k] = h[k] / l;
	}
	{
		double w = q[0], x = q[1], y = q[2], z = q[3];
		R[0] = 1 - 2 * (y * y + z * z); R[1] = 2 * (x * y - w * z); R[2] = 2 * (x * z + w * y);
		R[3] = 2 * (x * y + w * z); R[4] = 1 - 2 * (x * x + z * z); R[5] = 2 * (y * z - w * x);
		R[6] = 2 * (x * z - w * y); R[7] = 2 * (y * z + w * x); R[8] = 1 - 2 * (x * x + y * y);
	}
}

struct qaws_mesh_deformer
{
	qaws_deform_method method;
	unsigned int nv, nh, iters, solved;
	unsigned int* handles;
	double *x0, *x, *mass, *rhs0, *rhs, *q, *R, *prev, tol;
	unsigned char* move;
	sm_csr K, D;
	sm_system SK, SD;           /* K's system (LAPLACIAN, ARAP's rounds), K M^-1 K's (BIHARMONIC, ARAP's start) */
};

void qaws_mesh_deformer_destroy(qaws_mesh_deformer* d)
{
	if (!d)
		return;
	sm_system_free(&d->SK);
	sm_system_free(&d->SD);
	sm_csr_free(&d->K);
	sm_csr_free(&d->D);
	qaws_internal_free(d->handles); qaws_internal_free(d->x0); qaws_internal_free(d->x); qaws_internal_free(d->mass); qaws_internal_free(d->rhs0);
	qaws_internal_free(d->rhs); qaws_internal_free(d->q); qaws_internal_free(d->R); qaws_internal_free(d->prev); qaws_internal_free(d->move);
	qaws_internal_free(d);
}

qaws_status qaws_mesh_deformer_create(qaws_mesh_deformer** out, qaws_mesh_3d const* mesh, unsigned int const* handles, unsigned int nh,
	qaws_deform_desc const* desc)
{
	qaws_deform_method method = desc ? desc->method : QAWS_DEFORM_ARAP;
	int uniform = desc && desc->weights == QAWS_SMOOTH_WEIGHTS_UNIFORM;
	unsigned int nv, nt, i, k, p;
	unsigned int* t = NULL;
	qaws_scalar* pos = NULL;
	double* inv = NULL;
	qaws_mesh_deformer* d;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh || (nh && !handles) || (unsigned int)method > (unsigned int)QAWS_DEFORM_BIHARMONIC
		|| (desc && (unsigned int)desc->weights > (unsigned int)QAWS_SMOOTH_WEIGHTS_UNIFORM))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	for (i = 0; i < nh; i++)
		if (handles[i] >= nv)
			return QAWS_STATUS_INVALID_ARGUMENT;
	d = (qaws_mesh_deformer*)qaws_internal_zalloc(1, sizeof *d);
	if (!d)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	d->method = method;
	d->nv = nv;
	d->nh = nh;
	d->iters = desc && desc->iterations ? desc->iterations : 200;
	d->handles = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nh + 1));
	d->x0 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	d->x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	d->rhs0 = (double*)qaws_internal_zalloc(qaws_internal_mul2(3, (size_t)nv + 1), sizeof(double));
	d->rhs = (double*)qaws_internal_zalloc(qaws_internal_mul2(3, (size_t)nv + 1), sizeof(double));
	d->mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	d->move = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	if (!d->handles || !d->x0 || !d->x || !d->rhs0 || !d->rhs || !d->mass || !d->move || !pos || !t)
		goto done;
	if (nh)
		memcpy(d->handles, handles, sizeof(unsigned int) * nh);
	if ((st = qaws_mesh_3d_get_positions(pos, mesh)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_get_triangles(t, mesh)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * nv; i++)
		d->x0[i] = d->x[i] = (double)pos[i];
	for (i = 0; i < nv; i++)
		d->move[i] = !desc || !desc->region || desc->region[i];
	for (i = 0; i < nh; i++)
		d->move[handles[i]] = 0;
	/* (ARAP's rotations fit badly on negative weights, from stretched or obtuse triangles: clamped at 0 there) */
	if (!sm_laplacian(&d->K, d->mass, d->x0, nv, t, nt, uniform, method == QAWS_DEFORM_ARAP))
		goto done;
	if (method == QAWS_DEFORM_LAPLACIAN)
	{
		/* K x = K x0: the rest's Laplacian coordinates */
		for (i = 0; i < nv; i++)
			for (p = d->K.rp[i]; p < d->K.rp[i + 1]; p++)
				for (k = 0; k < 3; k++)
					d->rhs0[3 * (size_t)i + k] += d->K.v[p] * d->x0[3 * (size_t)d->K.ci[p] + k];
		st = sm_system_factor(&d->SK, &d->K, d->x0, d->move);
		goto done;
	}
	/* the biharmonic displacement: K M^-1 K d = 0, the handles' given (ARAP's start) */
	inv = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	if (!inv)
		goto done;
	for (i = 0; i < nv; i++)
		inv[i] = 1 / d->mass[i];
	if (!sm_mul(&d->D, &d->K, inv, &d->K))
		goto done;
	st = sm_system_factor(&d->SD, &d->D, d->x0, d->move);
	if (st != QAWS_STATUS_OK || method != QAWS_DEFORM_ARAP)
		goto done;
	{
		/* ARAP: K factored once for the rounds; they stop once no vertex moves by 1e-6 of the box's diagonal */
		double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		d->q = (double*)qaws_internal_zalloc(qaws_internal_mul2(4, (size_t)nv + 1), sizeof(double));
		d->R = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 9, (size_t)nv + 1));
		d->prev = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
		if (!d->q || !d->R || !d->prev)
			goto done;
		for (i = 0; i < nv; i++)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], d->x0[3 * (size_t)i + k]);
				hi[k] = fmax(hi[k], d->x0[3 * (size_t)i + k]);
			}
		d->tol = 1e-6 * sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
		for (i = 0; i < nv; i++)
			d->q[4 * (size_t)i] = 1;
		st = sm_system_factor(&d->SK, &d->K, d->x0, d->move);
	}
done:
	qaws_internal_free(pos);
	qaws_internal_free(t);
	qaws_internal_free(inv);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_deformer_destroy(d);
		return st;
	}
	*out = d;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_deformer_solve(qaws_scalar* out_positions, unsigned int* out_rounds, qaws_mesh_deformer* d, qaws_scalar const* targets)
{
	unsigned int nv, i, j, k, p, it, rounds = 0;
	double* x;
	qaws_status st = QAWS_STATUS_OK;
	if (!out_positions || !d || (d->nh && !targets))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * d->nh; i++)
		if (!isfinite((double)targets[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	nv = d->nv;
	x = d->x;
	if (d->method == QAWS_DEFORM_ARAP && d->solved)
	{
		/* ARAP again: from the last solution (its rotations kept) plus the biharmonic spread of the handles' moves since */
		memset(d->prev, 0, sizeof(double) * 3 * nv);
		for (i = 0; i < d->nh; i++)
			for (k = 0; k < 3; k++)
				d->prev[3 * (size_t)d->handles[i] + k] = (double)targets[3 * (size_t)i + k] - x[3 * (size_t)d->handles[i] + k];
		st = sm_system_solve(&d->SD, d->prev, d->rhs0, 3);
		for (i = 0; st == QAWS_STATUS_OK && i < 3 * nv; i++) x[i] += d->prev[i];
	}
	else
	{
		/* the held vertices at rest, the handles at their targets */
		for (i = 0; i < 3 * nv; i++) x[i] = d->x0[i];
		for (i = 0; i < d->nh; i++)
			for (k = 0; k < 3; k++) x[3 * (size_t)d->handles[i] + k] = (double)targets[3 * (size_t)i + k];
		if (d->method == QAWS_DEFORM_LAPLACIAN)
			st = sm_system_solve(&d->SK, x, d->rhs0, 3);
		else
		{
			for (i = 0; i < 3 * nv; i++) x[i] -= d->x0[i];
			st = sm_system_solve(&d->SD, x, d->rhs0, 3);
			for (i = 0; i < 3 * nv; i++) x[i] += d->x0[i];
		}
	}
	if (st == QAWS_STATUS_OK && d->method == QAWS_DEFORM_ARAP)
	{
		/* local: each ring's rotation; global: K x = b, b_i = sum_j w_ij (R_i + R_j) (x0_i - x0_j) / 2 */
		sm_csr const* K = &d->K;
		double const* x0 = d->x0;
		for (it = 0; it < d->iters; it++)
		{
			double moved = 0;
			for (i = 0; i < nv; i++)
			{
				double a[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
				for (p = K->rp[i]; p < K->rp[i + 1]; p++)
				{
					unsigned int jj = K->ci[p];
					double w = -K->v[p], e[3], e0[3];
					if (jj == i)
						continue;
					for (k = 0; k < 3; k++)
					{
						e[k] = x[3 * (size_t)jj + k] - x[3 * (size_t)i + k];
						e0[k] = x0[3 * (size_t)jj + k] - x0[3 * (size_t)i + k];
					}
					for (k = 0; k < 3; k++)
						for (j = 0; j < 3; j++)
							a[3 * k + j] += w * e[k] * e0[j];
				}
				df_rotation(a, d->q + 4 * (size_t)i, d->R + 9 * (size_t)i);
			}
			memset(d->rhs, 0, sizeof(double) * 3 * nv);
			for (i = 0; i < nv; i++)
				for (p = K->rp[i]; p < K->rp[i + 1]; p++)
				{
					unsigned int jj = K->ci[p];
					double w = -K->v[p], e0[3];
					if (jj == i)
						continue;
					for (k = 0; k < 3; k++)
						e0[k] = x0[3 * (size_t)i + k] - x0[3 * (size_t)jj + k];
					for (k = 0; k < 3; k++)
						for (j = 0; j < 3; j++)
							d->rhs[3 * (size_t)i + k] += 0.5 * w * (d->R[9 * (size_t)i + 3 * k + j] + d->R[9 * (size_t)jj + 3 * k + j]) * e0[j];
				}
			memcpy(d->prev, x, sizeof(double) * 3 * nv);
			st = sm_system_solve(&d->SK, x, d->rhs, 3);
			if (st != QAWS_STATUS_OK)
				break;
			rounds++;
			for (i = 0; i < 3 * nv; i++)
				moved = fmax(moved, fabs(x[i] - d->prev[i]));
			if (moved <= d->tol)
				break;
		}
	}
	if (st != QAWS_STATUS_OK)
		return st;
	d->solved = 1;
	for (i = 0; i < 3 * nv; i++)
		out_positions[i] = (qaws_scalar)x[i];
	if (out_rounds)
		*out_rounds = rounds;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_deform(qaws_mesh_3d* mesh, unsigned int const* handles, qaws_scalar const* targets, unsigned int nh,
	qaws_deform_desc const* desc)
{
	qaws_mesh_deformer* d = NULL;
	qaws_scalar* pos;
	qaws_status st;
	unsigned int i;
	if (!mesh || (nh && (!handles || !targets)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * nh; i++)
		if (!isfinite((double)targets[i]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	st = qaws_mesh_deformer_create(&d, mesh, handles, nh, desc);
	if (st != QAWS_STATUS_OK)
		return st;
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)d->nv + 1));
	st = pos ? qaws_mesh_deformer_solve(pos, NULL, d, targets) : QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_set_positions(mesh, pos);
	qaws_internal_free(pos);
	qaws_mesh_deformer_destroy(d);
	return st;
}
