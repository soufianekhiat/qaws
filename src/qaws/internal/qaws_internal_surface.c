#include "qaws_internal_surface.h"
#include "qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>

qaws_surface* qaws_internal_surface_alloc_ex(
	qaws_surface_kind kind,
	unsigned int u_degree,
	unsigned int v_degree,
	qaws_range u_range,
	qaws_range v_range,
	qaws_surface_vtable const* vtable,
	qaws_allocator const* allocator)
{
	qaws_surface* s = (qaws_surface*)qaws_internal_alloc(allocator, sizeof(qaws_surface));
	if (!s) return NULL;
	s->kind = kind;
	s->u_degree = u_degree;
	s->v_degree = v_degree;
	s->u_range = u_range;
	s->v_range = v_range;
	s->vtable = vtable;
	s->impl = NULL;
	s->allocator = allocator;
	return s;
}

qaws_surface* qaws_internal_surface_alloc(
	qaws_surface_kind kind,
	unsigned int u_degree,
	unsigned int v_degree,
	qaws_range u_range,
	qaws_range v_range,
	qaws_surface_vtable const* vtable)
{
	return qaws_internal_surface_alloc_ex(
		kind, u_degree, v_degree, u_range, v_range, vtable, NULL);
}

void qaws_internal_surface_free(qaws_surface* surface)
{
	if (surface)
	{
		if (surface->vtable && surface->vtable->destroy_impl && surface->impl)
			surface->vtable->destroy_impl(surface->impl, surface->allocator);
		qaws_internal_dealloc(surface->allocator, surface);
	}
}

unsigned int qaws_internal_surface_uniform_knots(
	unsigned int degree,
	unsigned int num_cp,
	qaws_scalar* out_knots,
	unsigned int capacity)
{
	unsigned int knot_count = num_cp + degree + 1;
	unsigned int i;
	unsigned int n_internal;
	if (knot_count > capacity) return 0;

	/* Clamped: first (degree+1) knots = 0, last (degree+1) knots = 1 */
	for (i = 0; i <= degree; i++)
		out_knots[i] = (qaws_scalar)0;
	n_internal = knot_count - 2 * (degree + 1);
	for (i = 0; i < n_internal; i++)
		out_knots[degree + 1 + i] = (qaws_scalar)(i + 1) / (qaws_scalar)(n_internal + 1);
	for (i = 0; i <= degree; i++)
		out_knots[knot_count - 1 - i] = (qaws_scalar)1;
	return knot_count;
}

void qaws_internal_surface_catmull_rom_blend(
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
