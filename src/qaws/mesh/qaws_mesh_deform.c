/*
 * Deformation by handles (qaws_mesh_deform.h): the rest shape's Laplacian
 * (qaws_mesh_smooth.c's sparse helpers), the linear editing systems, and
 * ARAP's local-global rounds on one factored system.
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

qaws_status qaws_mesh_3d_deform(qaws_mesh_3d* mesh, unsigned int const* handles, qaws_scalar const* targets, unsigned int nh,
	qaws_deform_desc const* desc)
{
	qaws_deform_method method = desc ? desc->method : QAWS_DEFORM_ARAP;
	int uniform = desc && desc->weights == QAWS_SMOOTH_WEIGHTS_UNIFORM;
	unsigned int iters = desc && desc->iterations ? desc->iterations : 200, nv, nt, i, j, k, p, it;
	double *x0 = NULL, *x = NULL, *mass = NULL, *rhs = NULL, *inv = NULL, *q = NULL, *R = NULL, *prev = NULL;
	unsigned int* t = NULL;
	unsigned char* move = NULL;
	qaws_scalar* pos = NULL;
	sm_csr K, D;
	sm_system S;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&K, 0, sizeof K);
	memset(&D, 0, sizeof D);
	memset(&S, 0, sizeof S);
	if (!mesh || (nh && (!handles || !targets)) || (unsigned int)method > (unsigned int)QAWS_DEFORM_BIHARMONIC
		|| (desc && (unsigned int)desc->weights > (unsigned int)QAWS_SMOOTH_WEIGHTS_UNIFORM))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	for (i = 0; i < nh; i++)
		if (handles[i] >= nv || !isfinite((double)targets[3 * (size_t)i]) || !isfinite((double)targets[3 * (size_t)i + 1])
			|| !isfinite((double)targets[3 * (size_t)i + 2]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	x0 = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	rhs = (double*)qaws_internal_zalloc(qaws_internal_mul2(3, (size_t)nv + 1), sizeof(double));
	mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	t = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	move = (unsigned char*)qaws_internal_malloc((size_t)nv + 1);
	if (!pos || !x0 || !x || !rhs || !mass || !t || !move)
		goto done;
	if ((st = qaws_mesh_3d_get_positions(pos, mesh)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_get_triangles(t, mesh)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * nv; i++)
		x0[i] = x[i] = (double)pos[i];
	for (i = 0; i < nv; i++)
		move[i] = !desc || !desc->region || desc->region[i];
	for (i = 0; i < nh; i++)
	{
		move[handles[i]] = 0;
		for (k = 0; k < 3; k++)
			x[3 * (size_t)handles[i] + k] = (double)targets[3 * (size_t)i + k];
	}
	/* (ARAP's rotations fit badly on negative weights, from stretched or obtuse triangles: clamped at 0 there) */
	if (!sm_laplacian(&K, mass, x0, nv, t, nt, uniform, method == QAWS_DEFORM_ARAP))
		goto done;
	if (method == QAWS_DEFORM_LAPLACIAN)
	{
		/* K x = K x0: the rest's Laplacian coordinates */
		for (i = 0; i < nv; i++)
			for (p = K.rp[i]; p < K.rp[i + 1]; p++)
				for (k = 0; k < 3; k++)
					rhs[3 * (size_t)i + k] += K.v[p] * x0[3 * (size_t)K.ci[p] + k];
		st = sm_system_factor(&S, &K, x0, move);
		if (st == QAWS_STATUS_OK)
			st = sm_system_solve(&S, x, rhs, 3);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	else
	{
		/* the biharmonic displacement: K M^-1 K d = 0, the handles' given (ARAP's start) */
		inv = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
		if (!inv)
			goto done;
		for (i = 0; i < nv; i++)
			inv[i] = 1 / mass[i];
		if (!sm_mul(&D, &K, inv, &K))
			goto done;
		for (i = 0; i < 3 * nv; i++)
			x[i] -= x0[i];
		st = sm_system_factor(&S, &D, x0, move);
		if (st == QAWS_STATUS_OK)
			st = sm_system_solve(&S, x, rhs, 3);
		sm_system_free(&S);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < 3 * nv; i++)
			x[i] += x0[i];
	}
	if (method == QAWS_DEFORM_ARAP)
	{
		/* local: each ring's rotation; global: K x = b, b_i = sum_j w_ij (R_i + R_j) (x0_i - x0_j) / 2, on K factored once */
		double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, tol;
		q = (double*)qaws_internal_zalloc(qaws_internal_mul2(4, (size_t)nv + 1), sizeof(double));
		R = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 9, (size_t)nv + 1));
		prev = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (!q || !R || !prev)
			goto done;
		for (i = 0; i < nv; i++)
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], x0[3 * (size_t)i + k]);
				hi[k] = fmax(hi[k], x0[3 * (size_t)i + k]);
			}
		/* (rounds until no vertex moves by more than 1e-6 of the box's diagonal) */
		tol = 1e-6 * sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
		for (i = 0; i < nv; i++)
			q[4 * (size_t)i] = 1;
		st = sm_system_factor(&S, &K, x0, move);
		if (st != QAWS_STATUS_OK)
			goto done;
		for (it = 0; it < iters; it++)
		{
			for (i = 0; i < nv; i++)
			{
				double a[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
				for (p = K.rp[i]; p < K.rp[i + 1]; p++)
				{
					unsigned int jj = K.ci[p];
					double w = -K.v[p], e[3], e0[3];
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
				df_rotation(a, q + 4 * (size_t)i, R + 9 * (size_t)i);
			}
			memset(rhs, 0, sizeof(double) * 3 * nv);
			for (i = 0; i < nv; i++)
				for (p = K.rp[i]; p < K.rp[i + 1]; p++)
				{
					unsigned int jj = K.ci[p];
					double w = -K.v[p], d[3];
					if (jj == i)
						continue;
					for (k = 0; k < 3; k++)
						d[k] = x0[3 * (size_t)i + k] - x0[3 * (size_t)jj + k];
					for (k = 0; k < 3; k++)
						for (j = 0; j < 3; j++)
							rhs[3 * (size_t)i + k] += 0.5 * w * (R[9 * (size_t)i + 3 * k + j] + R[9 * (size_t)jj + 3 * k + j]) * d[j];
				}
			memcpy(prev, x, sizeof(double) * 3 * nv);
			st = sm_system_solve(&S, x, rhs, 3);
			if (st != QAWS_STATUS_OK)
				goto done;
			{
				double moved = 0;
				for (i = 0; i < 3 * nv; i++)
					moved = fmax(moved, fabs(x[i] - prev[i]));
				if (moved <= tol)
					break;
			}
		}
	}
	for (i = 0; i < 3 * nv; i++)
		pos[i] = (qaws_scalar)x[i];
	st = qaws_mesh_3d_set_positions(mesh, pos);
done:
	sm_system_free(&S);
	sm_csr_free(&K);
	sm_csr_free(&D);
	qaws_internal_free(pos); qaws_internal_free(x0); qaws_internal_free(x); qaws_internal_free(rhs); qaws_internal_free(mass);
	qaws_internal_free(t); qaws_internal_free(move); qaws_internal_free(inv); qaws_internal_free(q); qaws_internal_free(R); qaws_internal_free(prev);
	return st;
}
