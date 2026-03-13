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

static void compute_normal(qaws_vec3 du, qaws_vec3 dv, qaws_vec3* out)
{
	qaws_scalar nx = du.y * dv.z - du.z * dv.y;
	qaws_scalar ny = du.z * dv.x - du.x * dv.z;
	qaws_scalar nz = du.x * dv.y - du.y * dv.x;
	qaws_scalar len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
	if (len > QAWS_LITERAL(1e-12))
	{
		out->x = nx / len; out->y = ny / len; out->z = nz / len;
	}
	else
	{
		out->x = 0; out->y = 0; out->z = 1;
	}
}

/* Catmull-Rom interpolation through points at non-uniform parameters.
   pts[n_pts], params[n_pts], evaluate at parameter t.
   Returns blended position. If out_deriv is non-NULL, also returns derivative w.r.t. t. */
static void catmull_rom_blend(
	qaws_vec3 const* pts, qaws_scalar const* params,
	unsigned int n_pts, qaws_scalar t,
	qaws_vec3* out_pos, qaws_vec3* out_deriv)
{
	unsigned int seg;
	unsigned int i0, i1, i2, i3;
	qaws_scalar s, dt;
	qaws_vec3 m0, m1;
	qaws_scalar s2, s3;
	qaws_scalar h00, h10, h01, h11;
	qaws_scalar dh00, dh10, dh01, dh11;
	qaws_scalar inv_dt;

	/* Special case: 2 points = linear interpolation */
	if (n_pts == 2)
	{
		qaws_scalar denom = params[1] - params[0];
		qaws_scalar alpha;
		if (QAWS_FABS(denom) < QAWS_LITERAL(1e-12))
			alpha = QAWS_LITERAL(0.5);
		else
			alpha = (t - params[0]) / denom;
		if (alpha < QAWS_ZERO) alpha = QAWS_ZERO;
		if (alpha > QAWS_ONE) alpha = QAWS_ONE;
		out_pos->x = (QAWS_ONE - alpha) * pts[0].x + alpha * pts[1].x;
		out_pos->y = (QAWS_ONE - alpha) * pts[0].y + alpha * pts[1].y;
		out_pos->z = (QAWS_ONE - alpha) * pts[0].z + alpha * pts[1].z;
		if (out_deriv)
		{
			qaws_scalar inv = (QAWS_FABS(denom) < QAWS_LITERAL(1e-12))
				? QAWS_ZERO : QAWS_ONE / denom;
			out_deriv->x = (pts[1].x - pts[0].x) * inv;
			out_deriv->y = (pts[1].y - pts[0].y) * inv;
			out_deriv->z = (pts[1].z - pts[0].z) * inv;
		}
		return;
	}

	/* Find segment: params[seg] <= t < params[seg+1] */
	seg = 0;
	{
		unsigned int k;
		for (k = 0; k < n_pts - 2; k++)
		{
			if (t < params[k + 1])
			{
				seg = k;
				break;
			}
			seg = k;
		}
		if (t >= params[n_pts - 2])
			seg = n_pts - 2;
	}

	/* Indices for the 4 surrounding points (clamped) */
	i1 = seg;
	i2 = seg + 1;
	i0 = (seg > 0) ? seg - 1 : 0;
	i3 = (seg + 2 < n_pts) ? seg + 2 : n_pts - 1;

	/* Local parameter s in [0,1] within the segment */
	dt = params[i2] - params[i1];
	if (QAWS_FABS(dt) < QAWS_LITERAL(1e-12))
	{
		*out_pos = pts[i1];
		if (out_deriv)
		{
			out_deriv->x = QAWS_ZERO;
			out_deriv->y = QAWS_ZERO;
			out_deriv->z = QAWS_ZERO;
		}
		return;
	}
	s = (t - params[i1]) / dt;
	if (s < QAWS_ZERO) s = QAWS_ZERO;
	if (s > QAWS_ONE) s = QAWS_ONE;

	/* Compute tangents at i1 and i2 using Catmull-Rom (non-uniform) */
	{
		qaws_scalar dp_prev = params[i2] - params[i0];
		if (QAWS_FABS(dp_prev) < QAWS_LITERAL(1e-12))
			dp_prev = QAWS_ONE;
		m0.x = (pts[i2].x - pts[i0].x) / dp_prev * dt;
		m0.y = (pts[i2].y - pts[i0].y) / dp_prev * dt;
		m0.z = (pts[i2].z - pts[i0].z) / dp_prev * dt;
	}
	{
		qaws_scalar dp_next = params[i3] - params[i1];
		if (QAWS_FABS(dp_next) < QAWS_LITERAL(1e-12))
			dp_next = QAWS_ONE;
		m1.x = (pts[i3].x - pts[i1].x) / dp_next * dt;
		m1.y = (pts[i3].y - pts[i1].y) / dp_next * dt;
		m1.z = (pts[i3].z - pts[i1].z) / dp_next * dt;
	}

	/* Hermite basis functions */
	s2 = s * s;
	s3 = s2 * s;
	h00 = QAWS_LITERAL(2.0) * s3 - QAWS_LITERAL(3.0) * s2 + QAWS_ONE;
	h10 = s3 - QAWS_LITERAL(2.0) * s2 + s;
	h01 = -QAWS_LITERAL(2.0) * s3 + QAWS_LITERAL(3.0) * s2;
	h11 = s3 - s2;

	out_pos->x = h00 * pts[i1].x + h10 * m0.x + h01 * pts[i2].x + h11 * m1.x;
	out_pos->y = h00 * pts[i1].y + h10 * m0.y + h01 * pts[i2].y + h11 * m1.y;
	out_pos->z = h00 * pts[i1].z + h10 * m0.z + h01 * pts[i2].z + h11 * m1.z;

	if (out_deriv)
	{
		/* Derivatives of Hermite basis w.r.t. s */
		dh00 = QAWS_LITERAL(6.0) * s2 - QAWS_LITERAL(6.0) * s;
		dh10 = QAWS_LITERAL(3.0) * s2 - QAWS_LITERAL(4.0) * s + QAWS_ONE;
		dh01 = -QAWS_LITERAL(6.0) * s2 + QAWS_LITERAL(6.0) * s;
		dh11 = QAWS_LITERAL(3.0) * s2 - QAWS_LITERAL(2.0) * s;

		/* ds/dt = 1/dt, so d/dt = d/ds * (1/dt) */
		inv_dt = QAWS_ONE / dt;
		out_deriv->x = (dh00 * pts[i1].x + dh10 * m0.x
			+ dh01 * pts[i2].x + dh11 * m1.x) * inv_dt;
		out_deriv->y = (dh00 * pts[i1].y + dh10 * m0.y
			+ dh01 * pts[i2].y + dh11 * m1.y) * inv_dt;
		out_deriv->z = (dh00 * pts[i1].z + dh10 * m0.z
			+ dh01 * pts[i2].z + dh11 * m1.z) * inv_dt;
	}
}

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
		catmull_rom_blend(sec_pos, impl->v_params, n, v,
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
		catmull_rom_blend(sec_du, impl->v_params, n, v,
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
			compute_normal(out_result->du, out_result->dv, &out_result->normal);
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
	loft_surface_is_rational
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
