#include "qaws_surface_biquadratic.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_biquadratic_impl
{
	qaws_vec3 cp[9];  /* 3x3 grid, row-major */
} qaws_surface_biquadratic_impl;

/* Biquadratic patch using Bernstein basis:
   B0(t) = (1-t)^2,  B1(t) = 2*t*(1-t),  B2(t) = t^2
   dB0(t) = -2(1-t), dB1(t) = 2-4t,      dB2(t) = 2t
   d2B0 = 2,         d2B1 = -4,           d2B2 = 2

   S(u,v) = sum_i sum_j cp[j*3+i] * Bi(u) * Bj(v)
   All derivatives computed analytically. */
static qaws_status biquadratic_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_biquadratic_impl const* impl =
		(qaws_surface_biquadratic_impl const*)surface->impl;
	qaws_scalar one_minus_u = QAWS_ONE - u;
	qaws_scalar one_minus_v = QAWS_ONE - v;

	/* Bernstein basis values */
	qaws_scalar Bu[3];
	qaws_scalar Bv[3];

	/* Bernstein first derivative values */
	qaws_scalar dBu[3];
	qaws_scalar dBv[3];

	/* Bernstein second derivative values (constants) */
	qaws_scalar d2Bu[3];
	qaws_scalar d2Bv[3];

	int need_pos = (eval_flags & QAWS_SURFACE_EVAL_POSITION) != 0;
	int need_du  = (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL)) != 0;
	int need_dv  = (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL)) != 0;
	int need_duu = (eval_flags & QAWS_SURFACE_EVAL_DUU) != 0;
	int need_dvv = (eval_flags & QAWS_SURFACE_EVAL_DVV) != 0;
	int need_duv = (eval_flags & QAWS_SURFACE_EVAL_DUV) != 0;

	Bu[0] = one_minus_u * one_minus_u;
	Bu[1] = QAWS_LITERAL(2.0) * u * one_minus_u;
	Bu[2] = u * u;

	Bv[0] = one_minus_v * one_minus_v;
	Bv[1] = QAWS_LITERAL(2.0) * v * one_minus_v;
	Bv[2] = v * v;

	dBu[0] = QAWS_LITERAL(-2.0) * one_minus_u;
	dBu[1] = QAWS_LITERAL(2.0) - QAWS_LITERAL(4.0) * u;
	dBu[2] = QAWS_LITERAL(2.0) * u;

	dBv[0] = QAWS_LITERAL(-2.0) * one_minus_v;
	dBv[1] = QAWS_LITERAL(2.0) - QAWS_LITERAL(4.0) * v;
	dBv[2] = QAWS_LITERAL(2.0) * v;

	d2Bu[0] = QAWS_LITERAL(2.0);
	d2Bu[1] = QAWS_LITERAL(-4.0);
	d2Bu[2] = QAWS_LITERAL(2.0);

	d2Bv[0] = QAWS_LITERAL(2.0);
	d2Bv[1] = QAWS_LITERAL(-4.0);
	d2Bv[2] = QAWS_LITERAL(2.0);

	/* Position: S(u,v) = sum_i sum_j cp[j*3+i] * Bi(u) * Bj(v) */
	if (need_pos)
	{
		qaws_scalar px = 0, py = 0, pz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = Bu[i] * Bv[j];
				px += impl->cp[j * 3 + i].x * w;
				py += impl->cp[j * 3 + i].y * w;
				pz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->position.x = px;
		out_result->position.y = py;
		out_result->position.z = pz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = sum_i sum_j cp[j*3+i] * dBi(u) * Bj(v) */
	if (need_du)
	{
		qaws_scalar dx = 0, dy = 0, dz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = dBu[i] * Bv[j];
				dx += impl->cp[j * 3 + i].x * w;
				dy += impl->cp[j * 3 + i].y * w;
				dz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->du.x = dx;
		out_result->du.y = dy;
		out_result->du.z = dz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = sum_i sum_j cp[j*3+i] * Bi(u) * dBj(v) */
	if (need_dv)
	{
		qaws_scalar dx = 0, dy = 0, dz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = Bu[i] * dBv[j];
				dx += impl->cp[j * 3 + i].x * w;
				dy += impl->cp[j * 3 + i].y * w;
				dz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->dv.x = dx;
		out_result->dv.y = dy;
		out_result->dv.z = dz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* d2S/du2 = sum_i sum_j cp[j*3+i] * d2Bi(u) * Bj(v) */
	if (need_duu)
	{
		qaws_scalar dx = 0, dy = 0, dz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = d2Bu[i] * Bv[j];
				dx += impl->cp[j * 3 + i].x * w;
				dy += impl->cp[j * 3 + i].y * w;
				dz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->duu.x = dx;
		out_result->duu.y = dy;
		out_result->duu.z = dz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* d2S/dv2 = sum_i sum_j cp[j*3+i] * Bi(u) * d2Bj(v) */
	if (need_dvv)
	{
		qaws_scalar dx = 0, dy = 0, dz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = Bu[i] * d2Bv[j];
				dx += impl->cp[j * 3 + i].x * w;
				dy += impl->cp[j * 3 + i].y * w;
				dz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->dvv.x = dx;
		out_result->dvv.y = dy;
		out_result->dvv.z = dz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* d2S/dudv = sum_i sum_j cp[j*3+i] * dBi(u) * dBj(v) */
	if (need_duv)
	{
		qaws_scalar dx = 0, dy = 0, dz = 0;
		int i, j;
		for (j = 0; j < 3; j++)
		{
			for (i = 0; i < 3; i++)
			{
				qaws_scalar w = dBu[i] * dBv[j];
				dx += impl->cp[j * 3 + i].x * w;
				dy += impl->cp[j * 3 + i].y * w;
				dz += impl->cp[j * 3 + i].z * w;
			}
		}
		out_result->duv.x = dx;
		out_result->duv.y = dy;
		out_result->duv.z = dz;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void biquadratic_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_internal_dealloc(allocator, impl);
}

static int biquadratic_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const biquadratic_surface_vtable = {
	biquadratic_surface_eval,
	biquadratic_surface_destroy,
	biquadratic_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_biquadratic(
	qaws_surface_biquadratic_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_biquadratic_impl* impl;
	qaws_range u_range, v_range;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_BIQUADRATIC,
		2, 2, u_range, v_range,
		&biquadratic_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_biquadratic_impl*)malloc(sizeof(qaws_surface_biquadratic_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	memcpy(impl->cp, desc->control_points, sizeof(qaws_vec3) * 9);

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
