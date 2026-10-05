#include "qaws_surface_loft.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_loft_impl
{
	qaws_curve const** sections;   /* copied array of N section curve pointers */
	unsigned int section_count;
	qaws_scalar* v_params;         /* copied array of N v-parameter values */
	qaws_range* section_ranges;    /* copied array of N parameter ranges */
} qaws_surface_loft_impl;

/* Catmull-Rom interpolation through points at non-uniform parameters.
   pts[n_pts], params[n_pts], evaluate at parameter t.
   Returns blended position. If out_deriv is non-NULL, also returns derivative w.r.t. t. */
/* Evaluate the loft at (u,v) by evaluating each section at u
   then blending in v with Catmull-Rom. */
static qaws_status loft_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_loft_impl const* impl =
		(qaws_surface_loft_impl const*)surface->impl;
	unsigned int n = impl->section_count;
	unsigned int i;
	qaws_vec3* sec_pos = NULL;
	qaws_vec3* sec_du = NULL;
	qaws_vec3 pos, dv_vec, du_vec;
	qaws_status status;
	int need_pos, need_du, need_dv;

	memset(&pos, 0, sizeof(pos));
	memset(&dv_vec, 0, sizeof(dv_vec));
	memset(&du_vec, 0, sizeof(du_vec));

	need_pos = (eval_flags & QAWS_SURFACE_EVAL_POSITION) != 0;
	need_du = (eval_flags & QAWS_SURFACE_EVAL_DU) != 0;
	need_dv = (eval_flags & QAWS_SURFACE_EVAL_DV) != 0;

	/* For second derivatives and normal, we need du and dv */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		need_pos = 1;
		need_du = 1;
		need_dv = 1;
	}

	/* Allocate temporary arrays for section evaluations */
	sec_pos = (qaws_vec3*)malloc(n * sizeof(qaws_vec3));
	if (!sec_pos) return QAWS_STATUS_ALLOCATION_FAILURE;

	if (need_du)
	{
		sec_du = (qaws_vec3*)malloc(n * sizeof(qaws_vec3));
		if (!sec_du)
		{
			free(sec_pos);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	/* Evaluate each section curve at parameter u */
	for (i = 0; i < n; i++)
	{
		qaws_scalar s_range = impl->section_ranges[i].max_value
			- impl->section_ranges[i].min_value;
		qaws_scalar t_curve = impl->section_ranges[i].min_value + u * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_du)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->sections[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(sec_pos);
			if (sec_du) free(sec_du);
			return status;
		}

		sec_pos[i] = cr.position;
		if (need_du)
		{
			sec_du[i].x = cr.d1.x * s_range;
			sec_du[i].y = cr.d1.y * s_range;
			sec_du[i].z = cr.d1.z * s_range;
		}
	}

	/* Blend section positions in v via Catmull-Rom */
	if (need_pos || need_dv)
	{
		qaws_internal_surface_catmull_rom_blend(sec_pos, impl->v_params, n, v,
			&pos, need_dv ? &dv_vec : NULL);
	}

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position = pos;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv = dv_vec;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* dS/du: blend section derivatives in v with same weights */
	if (need_du)
	{
		qaws_internal_surface_catmull_rom_blend(sec_du, impl->v_params, n, v,
			&du_vec, NULL);
		if (eval_flags & QAWS_SURFACE_EVAL_DU)
		{
			out_result->du = du_vec;
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
		}
	}

	/* Second derivatives via central finite differences */
	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar h = QAWS_LITERAL(1e-5);
		qaws_surface_eval_result r_lo, r_hi;

		if (eval_flags & QAWS_SURFACE_EVAL_DUU)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			loft_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duu.x = (r_hi.du.x - r_lo.du.x) * inv2h;
				out_result->duu.y = (r_hi.du.y - r_lo.du.y) * inv2h;
				out_result->duu.z = (r_hi.du.z - r_lo.du.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_DVV)
		{
			qaws_scalar v_lo = v - h, v_hi = v + h;
			qaws_scalar hv;
			if (v_lo < 0) v_lo = 0;
			if (v_hi > 1) v_hi = 1;
			hv = (v_hi - v_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			loft_surface_eval(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hv);
				out_result->dvv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->dvv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->dvv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_DUV)
		{
			qaws_scalar u_lo = u - h, u_hi = u + h;
			qaws_scalar hu;
			if (u_lo < 0) u_lo = 0;
			if (u_hi > 1) u_hi = 1;
			hu = (u_hi - u_lo) * QAWS_LITERAL(0.5);

			memset(&r_lo, 0, sizeof(r_lo));
			memset(&r_hi, 0, sizeof(r_hi));
			loft_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			loft_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

			{
				qaws_scalar inv2h = QAWS_ONE / (QAWS_LITERAL(2.0) * hu);
				out_result->duv.x = (r_hi.dv.x - r_lo.dv.x) * inv2h;
				out_result->duv.y = (r_hi.dv.y - r_lo.dv.y) * inv2h;
				out_result->duv.z = (r_hi.dv.z - r_lo.dv.z) * inv2h;
			}
			out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
		}

		if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
		{
			qaws_internal_surface_normal(out_result->du, out_result->dv, &out_result->normal);
			out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
		}
	}

	free(sec_pos);
	if (sec_du) free(sec_du);
	return QAWS_STATUS_OK;
}

static void loft_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_loft_impl* li = (qaws_surface_loft_impl*)impl;
	if (li)
	{
		if (li->sections) free(li->sections);
		if (li->v_params) free(li->v_params);
		if (li->section_ranges) free(li->section_ranges);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int loft_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const loft_surface_vtable = {
	loft_surface_eval,
	loft_surface_destroy,
	loft_surface_is_rational,
	NULL /* diff */
};

qaws_status qaws_surface_create_loft(
	qaws_surface_loft_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_loft_impl* impl;
	qaws_range u_range, v_range;
	unsigned int i;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->sections) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->section_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate all sections are 3D */
	for (i = 0; i < desc->section_count; i++)
	{
		if (!desc->sections[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->sections[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_LOFT,
		0, 0, u_range, v_range,
		&loft_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_loft_impl*)malloc(sizeof(qaws_surface_loft_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_surface_loft_impl));

	impl->section_count = desc->section_count;

	/* Copy sections array */
	impl->sections = (qaws_curve const**)malloc(
		desc->section_count * sizeof(qaws_curve const*));
	if (!impl->sections)
	{
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->sections, desc->sections,
		desc->section_count * sizeof(qaws_curve const*));

	/* Copy or generate v_parameters */
	impl->v_params = (qaws_scalar*)malloc(
		desc->section_count * sizeof(qaws_scalar));
	if (!impl->v_params)
	{
		free(impl->sections);
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	if (desc->v_parameters)
	{
		memcpy(impl->v_params, desc->v_parameters,
			desc->section_count * sizeof(qaws_scalar));
	}
	else
	{
		/* Uniform spacing */
		for (i = 0; i < desc->section_count; i++)
		{
			impl->v_params[i] = (qaws_scalar)i
				/ (qaws_scalar)(desc->section_count - 1);
		}
	}

	/* Cache section parameter ranges */
	impl->section_ranges = (qaws_range*)malloc(
		desc->section_count * sizeof(qaws_range));
	if (!impl->section_ranges)
	{
		free(impl->v_params);
		free(impl->sections);
		free(impl);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < desc->section_count; i++)
	{
		impl->section_ranges[i] = qaws_curve_get_parameter_range(desc->sections[i]);
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}
