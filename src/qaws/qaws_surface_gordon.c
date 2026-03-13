#include "qaws_surface_gordon.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

typedef struct qaws_surface_gordon_impl
{
	qaws_curve const** u_curves;    /* copied array of M u-curve pointers */
	unsigned int u_curve_count;      /* M */
	qaws_scalar* v_params;           /* copied array of M v-parameter values */
	qaws_range* u_curve_ranges;      /* copied array of M parameter ranges */

	qaws_curve const** v_curves;    /* copied array of N v-curve pointers */
	unsigned int v_curve_count;      /* N */
	qaws_scalar* u_params;           /* copied array of N u-parameter values */
	qaws_range* v_curve_ranges;      /* copied array of N parameter ranges */

	qaws_vec3* Q;                    /* precomputed intersection grid Q[M*N], row-major */
} qaws_surface_gordon_impl;

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

/* Evaluate the tensor product term T(u,v) by blending Q[M*N] grid.
   First blend each row (fixed i) in u, then blend those results in v. */
static void tensor_product_eval(
	qaws_surface_gordon_impl const* impl,
	qaws_scalar u, qaws_scalar v,
	qaws_vec3* out_pos, qaws_vec3* out_du, qaws_vec3* out_dv)
{
	unsigned int M = impl->u_curve_count;
	unsigned int N = impl->v_curve_count;
	unsigned int i;
	qaws_vec3* row_blended;
	qaws_vec3* row_du;

	row_blended = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
	if (!row_blended) return;

	row_du = NULL;
	if (out_du)
	{
		row_du = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
		if (!row_du)
		{
			free(row_blended);
			return;
		}
	}

	/* For each u-curve index i, blend the N intersection points Q[i][0..N-1] at u */
	for (i = 0; i < M; i++)
	{
		qaws_vec3 const* row = &impl->Q[i * N];
		catmull_rom_blend(row, impl->u_params, N, u,
			&row_blended[i], row_du ? &row_du[i] : NULL);
	}

	/* Blend the M results in v */
	catmull_rom_blend(row_blended, impl->v_params, M, v,
		out_pos, out_dv);

	/* For dT/du: blend the M du-results in v */
	if (out_du)
	{
		catmull_rom_blend(row_du, impl->v_params, M, v,
			out_du, NULL);
	}

	free(row_blended);
	if (row_du) free(row_du);
}

/* Gordon surface: S(u,v) = L_u(u,v) + L_v(u,v) - T(u,v)
   L_u: evaluate each u-curve at u, blend in v
   L_v: evaluate each v-curve at v, blend in u
   T: tensor product of intersection grid */
static qaws_status gordon_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_gordon_impl const* impl =
		(qaws_surface_gordon_impl const*)surface->impl;
	unsigned int M = impl->u_curve_count;
	unsigned int N = impl->v_curve_count;
	unsigned int i;
	qaws_vec3* u_pts = NULL;
	qaws_vec3* u_du = NULL;
	qaws_vec3* v_pts = NULL;
	qaws_vec3* v_dv = NULL;
	qaws_vec3 Lu_pos, Lu_dv, Lu_du;
	qaws_vec3 Lv_pos, Lv_du, Lv_dv;
	qaws_vec3 T_pos, T_du, T_dv;
	qaws_status status;
	int need_pos, need_du, need_dv;

	need_pos = (eval_flags & QAWS_SURFACE_EVAL_POSITION) != 0;
	need_du = (eval_flags & QAWS_SURFACE_EVAL_DU) != 0;
	need_dv = (eval_flags & QAWS_SURFACE_EVAL_DV) != 0;

	if (eval_flags & (QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DVV
		| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_NORMAL))
	{
		need_pos = 1;
		need_du = 1;
		need_dv = 1;
	}

	/* Initialize vectors */
	memset(&Lu_pos, 0, sizeof(Lu_pos));
	memset(&Lu_dv, 0, sizeof(Lu_dv));
	memset(&Lu_du, 0, sizeof(Lu_du));
	memset(&Lv_pos, 0, sizeof(Lv_pos));
	memset(&Lv_du, 0, sizeof(Lv_du));
	memset(&Lv_dv, 0, sizeof(Lv_dv));
	memset(&T_pos, 0, sizeof(T_pos));
	memset(&T_du, 0, sizeof(T_du));
	memset(&T_dv, 0, sizeof(T_dv));

	/* === L_u(u,v): evaluate each u-curve at u, blend in v === */
	u_pts = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
	if (!u_pts) return QAWS_STATUS_ALLOCATION_FAILURE;

	u_du = NULL;
	if (need_du)
	{
		u_du = (qaws_vec3*)malloc(M * sizeof(qaws_vec3));
		if (!u_du)
		{
			free(u_pts);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	for (i = 0; i < M; i++)
	{
		qaws_scalar s_range = impl->u_curve_ranges[i].max_value
			- impl->u_curve_ranges[i].min_value;
		qaws_scalar t_curve = impl->u_curve_ranges[i].min_value + u * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_du)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->u_curves[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(u_pts);
			if (u_du) free(u_du);
			return status;
		}

		u_pts[i] = cr.position;
		if (need_du)
		{
			u_du[i].x = cr.d1.x * s_range;
			u_du[i].y = cr.d1.y * s_range;
			u_du[i].z = cr.d1.z * s_range;
		}
	}

	catmull_rom_blend(u_pts, impl->v_params, M, v,
		&Lu_pos, need_dv ? &Lu_dv : NULL);

	if (need_du)
	{
		catmull_rom_blend(u_du, impl->v_params, M, v,
			&Lu_du, NULL);
	}

	free(u_pts);
	if (u_du) free(u_du);

	/* === L_v(u,v): evaluate each v-curve at v, blend in u === */
	v_pts = (qaws_vec3*)malloc(N * sizeof(qaws_vec3));
	if (!v_pts) return QAWS_STATUS_ALLOCATION_FAILURE;

	v_dv = NULL;
	if (need_dv)
	{
		v_dv = (qaws_vec3*)malloc(N * sizeof(qaws_vec3));
		if (!v_dv)
		{
			free(v_pts);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	for (i = 0; i < N; i++)
	{
		qaws_scalar s_range = impl->v_curve_ranges[i].max_value
			- impl->v_curve_ranges[i].min_value;
		qaws_scalar t_curve = impl->v_curve_ranges[i].min_value + v * s_range;
		unsigned int cf = QAWS_EVAL_FLAG_POSITION;
		qaws_eval_result_3d cr;

		if (need_dv)
			cf |= QAWS_EVAL_FLAG_D1;

		memset(&cr, 0, sizeof(cr));
		status = qaws_curve_evaluate_3d(impl->v_curves[i], t_curve, cf, &cr);
		if (status != QAWS_STATUS_OK)
		{
			free(v_pts);
			if (v_dv) free(v_dv);
			return status;
		}

		v_pts[i] = cr.position;
		if (need_dv)
		{
			v_dv[i].x = cr.d1.x * s_range;
			v_dv[i].y = cr.d1.y * s_range;
			v_dv[i].z = cr.d1.z * s_range;
		}
	}

	catmull_rom_blend(v_pts, impl->u_params, N, u,
		&Lv_pos, need_du ? &Lv_du : NULL);

	if (need_dv)
	{
		catmull_rom_blend(v_dv, impl->u_params, N, u,
			&Lv_dv, NULL);
	}

	free(v_pts);
	if (v_dv) free(v_dv);

	/* === T(u,v): tensor product of intersection grid === */
	tensor_product_eval(impl, u, v, &T_pos,
		need_du ? &T_du : NULL,
		need_dv ? &T_dv : NULL);

	/* === Combine: S = L_u + L_v - T === */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = Lu_pos.x + Lv_pos.x - T_pos.x;
		out_result->position.y = Lu_pos.y + Lv_pos.y - T_pos.y;
		out_result->position.z = Lu_pos.z + Lv_pos.z - T_pos.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DU)
	{
		out_result->du.x = Lu_du.x + Lv_du.x - T_du.x;
		out_result->du.y = Lu_du.y + Lv_du.y - T_du.y;
		out_result->du.z = Lu_du.z + Lv_du.z - T_du.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DV)
	{
		out_result->dv.x = Lu_dv.x + Lv_dv.x - T_dv.x;
		out_result->dv.y = Lu_dv.y + Lv_dv.y - T_dv.y;
		out_result->dv.z = Lu_dv.z + Lv_dv.z - T_dv.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
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
			gordon_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DU, &r_lo);
			gordon_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DU, &r_hi);

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
			gordon_surface_eval(surface, u, v_lo, QAWS_SURFACE_EVAL_DV, &r_lo);
			gordon_surface_eval(surface, u, v_hi, QAWS_SURFACE_EVAL_DV, &r_hi);

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
			gordon_surface_eval(surface, u_lo, v, QAWS_SURFACE_EVAL_DV, &r_lo);
			gordon_surface_eval(surface, u_hi, v, QAWS_SURFACE_EVAL_DV, &r_hi);

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

	return QAWS_STATUS_OK;
}

static void gordon_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_gordon_impl* gi = (qaws_surface_gordon_impl*)impl;
	if (gi)
	{
		if (gi->u_curves) free(gi->u_curves);
		if (gi->v_params) free(gi->v_params);
		if (gi->u_curve_ranges) free(gi->u_curve_ranges);
		if (gi->v_curves) free(gi->v_curves);
		if (gi->u_params) free(gi->u_params);
		if (gi->v_curve_ranges) free(gi->v_curve_ranges);
		if (gi->Q) free(gi->Q);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int gordon_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const gordon_surface_vtable = {
	gordon_surface_eval,
	gordon_surface_destroy,
	gordon_surface_is_rational
};

qaws_status qaws_surface_create_gordon(
	qaws_surface_gordon_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_gordon_impl* impl;
	qaws_range u_range, v_range;
	unsigned int M, N, i, j;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->u_curves || !desc->v_curves) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->v_params || !desc->u_params) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->u_curve_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->v_curve_count < 2) return QAWS_STATUS_INVALID_ARGUMENT;

	M = desc->u_curve_count;
	N = desc->v_curve_count;

	/* Validate all curves are 3D */
	for (i = 0; i < M; i++)
	{
		if (!desc->u_curves[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->u_curves[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	for (i = 0; i < N; i++)
	{
		if (!desc->v_curves[i]) return QAWS_STATUS_INVALID_ARGUMENT;
		if (qaws_curve_get_dimension(desc->v_curves[i]) != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}

	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_GORDON,
		0, 0, u_range, v_range,
		&gordon_surface_vtable);
	if (!surface) return QAWS_STATUS_ALLOCATION_FAILURE;

	impl = (qaws_surface_gordon_impl*)malloc(sizeof(qaws_surface_gordon_impl));
	if (!impl)
	{
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(impl, 0, sizeof(qaws_surface_gordon_impl));

	impl->u_curve_count = M;
	impl->v_curve_count = N;

	/* Copy u_curves array */
	impl->u_curves = (qaws_curve const**)malloc(M * sizeof(qaws_curve const*));
	if (!impl->u_curves) goto fail;
	memcpy(impl->u_curves, desc->u_curves, M * sizeof(qaws_curve const*));

	/* Copy v_params */
	impl->v_params = (qaws_scalar*)malloc(M * sizeof(qaws_scalar));
	if (!impl->v_params) goto fail;
	memcpy(impl->v_params, desc->v_params, M * sizeof(qaws_scalar));

	/* Cache u-curve parameter ranges */
	impl->u_curve_ranges = (qaws_range*)malloc(M * sizeof(qaws_range));
	if (!impl->u_curve_ranges) goto fail;
	for (i = 0; i < M; i++)
		impl->u_curve_ranges[i] = qaws_curve_get_parameter_range(desc->u_curves[i]);

	/* Copy v_curves array */
	impl->v_curves = (qaws_curve const**)malloc(N * sizeof(qaws_curve const*));
	if (!impl->v_curves) goto fail;
	memcpy(impl->v_curves, desc->v_curves, N * sizeof(qaws_curve const*));

	/* Copy u_params */
	impl->u_params = (qaws_scalar*)malloc(N * sizeof(qaws_scalar));
	if (!impl->u_params) goto fail;
	memcpy(impl->u_params, desc->u_params, N * sizeof(qaws_scalar));

	/* Cache v-curve parameter ranges */
	impl->v_curve_ranges = (qaws_range*)malloc(N * sizeof(qaws_range));
	if (!impl->v_curve_ranges) goto fail;
	for (i = 0; i < N; i++)
		impl->v_curve_ranges[i] = qaws_curve_get_parameter_range(desc->v_curves[i]);

	/* Precompute intersection grid Q[M*N].
	   Q[i][j] = u_curve[i] evaluated at the u-parameter corresponding to u_params[j].
	   u_curve[i] runs along u, so evaluate at:
	     t = range.min + u_params[j] * (range.max - range.min) */
	impl->Q = (qaws_vec3*)malloc(M * N * sizeof(qaws_vec3));
	if (!impl->Q) goto fail;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			qaws_scalar s_range = impl->u_curve_ranges[i].max_value
				- impl->u_curve_ranges[i].min_value;
			qaws_scalar t_curve = impl->u_curve_ranges[i].min_value
				+ desc->u_params[j] * s_range;
			qaws_eval_result_3d cr;
			qaws_status status;

			memset(&cr, 0, sizeof(cr));
			status = qaws_curve_evaluate_3d(desc->u_curves[i], t_curve,
				QAWS_EVAL_FLAG_POSITION, &cr);
			if (status != QAWS_STATUS_OK)
			{
				gordon_surface_destroy(impl, NULL);
				qaws_internal_surface_free(surface);
				return status;
			}
			impl->Q[i * N + j] = cr.position;
		}
	}

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;

fail:
	gordon_surface_destroy(impl, NULL);
	qaws_internal_surface_free(surface);
	return QAWS_STATUS_ALLOCATION_FAILURE;
}
