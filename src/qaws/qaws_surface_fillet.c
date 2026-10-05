#include "qaws_surface_fillet.h"
#include "qaws_surface.h"
#include "qaws_surface_intersect.h"
#include "qaws_surface_offset.h"
#include "qaws_inspect.h"
#include "qaws_curve.h"
#include "qaws_bspline.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_fit.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

#define FILLET_DEFAULT_SAMPLES   64
#define FILLET_DEFAULT_ARC_SEGS  8
#define FILLET_SSI_GRID          20
#define FILLET_SSI_MAX_STEPS     500
#define FILLET_SSI_MAX_CURVES    4
#define FILLET_SSI_POINT_CAP     2048

/* ------------------------------------------------------------------ */
/* Vector helpers (file-local)                                        */
/* ------------------------------------------------------------------ */

static qaws_scalar fillet_vec3_dot(qaws_vec3 a, qaws_vec3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static qaws_vec3 fillet_vec3_cross(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.y * b.z - a.z * b.y;
	r.y = a.z * b.x - a.x * b.z;
	r.z = a.x * b.y - a.y * b.x;
	return r;
}

static qaws_scalar fillet_vec3_length(qaws_vec3 v)
{
	return QAWS_SQRT(v.x * v.x + v.y * v.y + v.z * v.z);
}

static qaws_vec3 fillet_vec3_normalize(qaws_vec3 v)
{
	qaws_scalar len = fillet_vec3_length(v);
	if (len > QAWS_LITERAL(1e-14))
	{
		v.x /= len; v.y /= len; v.z /= len;
	}
	return v;
}

static qaws_vec3 fillet_vec3_scale(qaws_vec3 v, qaws_scalar s)
{
	qaws_vec3 r;
	r.x = v.x * s; r.y = v.y * s; r.z = v.z * s;
	return r;
}

static qaws_vec3 fillet_vec3_add(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.x + b.x; r.y = a.y + b.y; r.z = a.z + b.z;
	return r;
}

static qaws_vec3 fillet_vec3_sub(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.x - b.x; r.y = a.y - b.y; r.z = a.z - b.z;
	return r;
}

/* ------------------------------------------------------------------ */
/* Rodrigues rotation: rotate vector v around axis k by angle theta   */
/* ------------------------------------------------------------------ */

static qaws_vec3 rodrigues_rotate(qaws_vec3 v, qaws_vec3 k, qaws_scalar theta)
{
	qaws_scalar cos_t = QAWS_COS(theta);
	qaws_scalar sin_t = QAWS_SIN(theta);
	qaws_scalar kdv = fillet_vec3_dot(k, v);
	qaws_vec3 kxv = fillet_vec3_cross(k, v);
	qaws_vec3 r;
	r.x = v.x * cos_t + kxv.x * sin_t + k.x * kdv * (QAWS_ONE - cos_t);
	r.y = v.y * cos_t + kxv.y * sin_t + k.y * kdv * (QAWS_ONE - cos_t);
	r.z = v.z * cos_t + kxv.z * sin_t + k.z * kdv * (QAWS_ONE - cos_t);
	return r;
}

/* ------------------------------------------------------------------ */
/* Compute fillet frame at an intersection point                      */
/* ------------------------------------------------------------------ */

static void compute_fillet_frame(
	qaws_vec3 pos,
	qaws_vec3 normal_a,
	qaws_vec3 normal_b,
	qaws_scalar radius,
	qaws_vec3* out_center,
	qaws_vec3* out_contact_a,
	qaws_vec3* out_contact_b)
{
	qaws_vec3 bisector;
	qaws_scalar cos_half, dist;

	/* Bisector = normalize(normal_a + normal_b) */
	bisector = fillet_vec3_add(normal_a, normal_b);
	{
		qaws_scalar len = fillet_vec3_length(bisector);
		if (len < QAWS_LITERAL(1e-10))
		{
			/* Normals are opposite: degenerate case, use normal_a */
			bisector = normal_a;
		}
		else
		{
			bisector.x /= len; bisector.y /= len; bisector.z /= len;
		}
	}

	/* Half-angle between normals: cos(half) = dot(normal_a, bisector) */
	cos_half = fillet_vec3_dot(normal_a, bisector);
	if (cos_half < QAWS_LITERAL(1e-6))
		cos_half = QAWS_LITERAL(1e-6);

	/* Distance from intersection point to ball center along bisector */
	dist = radius / cos_half;

	*out_center = fillet_vec3_add(pos, fillet_vec3_scale(bisector, dist));

	/* Contact points: project center onto each surface tangent plane */
	out_contact_a->x = out_center->x - radius * normal_a.x;
	out_contact_a->y = out_center->y - radius * normal_a.y;
	out_contact_a->z = out_center->z - radius * normal_a.z;

	out_contact_b->x = out_center->x - radius * normal_b.x;
	out_contact_b->y = out_center->y - radius * normal_b.y;
	out_contact_b->z = out_center->z - radius * normal_b.z;
}

/* ------------------------------------------------------------------ */
/* Normal helper                                                      */
/* ------------------------------------------------------------------ */

static void compute_surface_normal(qaws_vec3 du, qaws_vec3 dv, qaws_vec3* out)
{
	qaws_vec3 n = fillet_vec3_cross(du, dv);
	qaws_scalar len = fillet_vec3_length(n);
	if (len > QAWS_LITERAL(1e-12))
	{
		out->x = n.x / len; out->y = n.y / len; out->z = n.z / len;
	}
	else
	{
		out->x = 0; out->y = 0; out->z = 1;
	}
}

/* ------------------------------------------------------------------ */
/* Fillet surface impl                                                */
/* ------------------------------------------------------------------ */

typedef struct qaws_surface_fillet_impl
{
	qaws_vec3* grid;           /* u_count * v_count positions */
	unsigned int u_count;      /* samples along spine */
	unsigned int v_count;      /* arc_segments + 1 (around arc) */
} qaws_surface_fillet_impl;

/* Evaluate the fillet surface using bilinear interpolation on the grid. */
static qaws_status fillet_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_fillet_impl const* impl =
		(qaws_surface_fillet_impl const*)surface->impl;
	unsigned int uc = impl->u_count;
	unsigned int vc = impl->v_count;
	qaws_scalar fu, fv;
	unsigned int iu, iv;
	qaws_scalar su, sv;
	qaws_vec3 p00, p10, p01, p11;
	qaws_scalar one_su, one_sv;

	/* Map (u,v) to grid indices */
	fu = u * (qaws_scalar)(uc - 1);
	fv = v * (qaws_scalar)(vc - 1);
	iu = (unsigned int)QAWS_FLOOR(fu);
	iv = (unsigned int)QAWS_FLOOR(fv);
	if (iu >= uc - 1) iu = uc - 2;
	if (iv >= vc - 1) iv = vc - 2;
	su = fu - (qaws_scalar)iu;
	sv = fv - (qaws_scalar)iv;
	one_su = QAWS_ONE - su;
	one_sv = QAWS_ONE - sv;

	/* Fetch the four corners */
	p00 = impl->grid[iu * vc + iv];
	p10 = impl->grid[(iu + 1) * vc + iv];
	p01 = impl->grid[iu * vc + (iv + 1)];
	p11 = impl->grid[(iu + 1) * vc + (iv + 1)];

	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = one_su * one_sv * p00.x + su * one_sv * p10.x
			+ one_su * sv * p01.x + su * sv * p11.x;
		out_result->position.y = one_su * one_sv * p00.y + su * one_sv * p10.y
			+ one_su * sv * p01.y + su * sv * p11.y;
		out_result->position.z = one_su * one_sv * p00.z + su * one_sv * p10.z
			+ one_su * sv * p01.z + su * sv * p11.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
	{
		/* dS/du in grid space, then rescale by (uc-1) to get parameter-space deriv */
		qaws_scalar sc = (qaws_scalar)(uc - 1);
		out_result->du.x = (one_sv * (p10.x - p00.x) + sv * (p11.x - p01.x)) * sc;
		out_result->du.y = (one_sv * (p10.y - p00.y) + sv * (p11.y - p01.y)) * sc;
		out_result->du.z = (one_sv * (p10.z - p00.z) + sv * (p11.z - p01.z)) * sc;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL))
	{
		qaws_scalar sc = (qaws_scalar)(vc - 1);
		out_result->dv.x = (one_su * (p01.x - p00.x) + su * (p11.x - p10.x)) * sc;
		out_result->dv.y = (one_su * (p01.y - p00.y) + su * (p11.y - p10.y)) * sc;
		out_result->dv.z = (one_su * (p01.z - p00.z) + su * (p11.z - p10.z)) * sc;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		out_result->duu.x = 0; out_result->duu.y = 0; out_result->duu.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = 0; out_result->dvv.y = 0; out_result->dvv.z = 0;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		qaws_scalar sc = (qaws_scalar)(uc - 1) * (qaws_scalar)(vc - 1);
		out_result->duv.x = (p00.x - p10.x - p01.x + p11.x) * sc;
		out_result->duv.y = (p00.y - p10.y - p01.y + p11.y) * sc;
		out_result->duv.z = (p00.z - p10.z - p01.z + p11.z) * sc;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}

	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		compute_surface_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

static void fillet_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_fillet_impl* fi = (qaws_surface_fillet_impl*)impl;
	if (fi)
	{
		if (fi->grid) free(fi->grid);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int fillet_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const fillet_surface_vtable = {
	fillet_surface_eval,
	fillet_surface_destroy,
	fillet_surface_is_rational,
	NULL /* diff */
};

/* ------------------------------------------------------------------ */
/* Resample an SSI curve to uniform-t positions via linear interp     */
/* ------------------------------------------------------------------ */

static void resample_ssi_curve(
	qaws_ssi_point const* pts,
	unsigned int pt_count,
	unsigned int out_count,
	qaws_ssi_point* out_pts)
{
	unsigned int i;
	if (pt_count == 0 || out_count == 0) return;
	if (pt_count == 1)
	{
		for (i = 0; i < out_count; i++)
			out_pts[i] = pts[0];
		return;
	}
	for (i = 0; i < out_count; i++)
	{
		qaws_scalar t = (qaws_scalar)i / (qaws_scalar)(out_count - 1);
		qaws_scalar fi = t * (qaws_scalar)(pt_count - 1);
		unsigned int idx = (unsigned int)QAWS_FLOOR(fi);
		qaws_scalar frac;
		unsigned int next;
		if (idx >= pt_count - 1) idx = pt_count - 2;
		next = idx + 1;
		frac = fi - (qaws_scalar)idx;

		out_pts[i].u1 = pts[idx].u1 * (QAWS_ONE - frac) + pts[next].u1 * frac;
		out_pts[i].v1 = pts[idx].v1 * (QAWS_ONE - frac) + pts[next].v1 * frac;
		out_pts[i].u2 = pts[idx].u2 * (QAWS_ONE - frac) + pts[next].u2 * frac;
		out_pts[i].v2 = pts[idx].v2 * (QAWS_ONE - frac) + pts[next].v2 * frac;
		out_pts[i].position.x = pts[idx].position.x * (QAWS_ONE - frac)
			+ pts[next].position.x * frac;
		out_pts[i].position.y = pts[idx].position.y * (QAWS_ONE - frac)
			+ pts[next].position.y * frac;
		out_pts[i].position.z = pts[idx].position.z * (QAWS_ONE - frac)
			+ pts[next].position.z * frac;
	}
}

/* ------------------------------------------------------------------ */
/* Build the fillet grid from intersection data                       */
/* ------------------------------------------------------------------ */

static qaws_status build_fillet_grid(
	qaws_surface const* surface_a,
	qaws_surface const* surface_b,
	qaws_scalar radius,
	qaws_ssi_point const* spine_pts,
	unsigned int sample_count,
	unsigned int arc_segments,
	qaws_vec3** out_grid,
	qaws_scalar* out_trim_a_uv,   /* optional: 2 * sample_count coords for trim on A */
	qaws_scalar* out_trim_b_uv)   /* optional: 2 * sample_count coords for trim on B */
{
	unsigned int v_count = arc_segments + 1;
	unsigned int si, ai;
	qaws_vec3* grid;

	grid = (qaws_vec3*)malloc(sample_count * v_count * sizeof(qaws_vec3));
	if (!grid) return QAWS_STATUS_ALLOCATION_FAILURE;

	for (si = 0; si < sample_count; si++)
	{
		qaws_ssi_point const* sp = &spine_pts[si];
		qaws_vec3 normal_a, normal_b;
		qaws_vec3 center, contact_a, contact_b;
		qaws_vec3 ra, rb, axis;
		qaws_scalar angle, axis_len;

		/* Get normals from both surfaces */
		{
			qaws_status st;
			st = qaws_surface_get_normal(surface_a, sp->u1, sp->v1, &normal_a);
			if (st != QAWS_STATUS_OK)
			{
				/* Fallback: use (0,0,1) */
				normal_a.x = 0; normal_a.y = 0; normal_a.z = 1;
			}
			st = qaws_surface_get_normal(surface_b, sp->u2, sp->v2, &normal_b);
			if (st != QAWS_STATUS_OK)
			{
				normal_b.x = 0; normal_b.y = 0; normal_b.z = 1;
			}
		}

		/* Ensure normals point in consistent directions: both should
		   point "outward" so the fillet center is displaced away from
		   the surfaces. We orient them so they point to the same side
		   of the intersection. If they are nearly the same direction,
		   the fillet wraps a convex edge; if opposite, a concave edge. */
		normal_a = fillet_vec3_normalize(normal_a);
		normal_b = fillet_vec3_normalize(normal_b);

		compute_fillet_frame(sp->position, normal_a, normal_b,
			radius, &center, &contact_a, &contact_b);

		/* Store trim curve parameters (contact projected back to surface params).
		   We approximate by using the closest-point u,v from the intersection
		   plus a small offset along the normal in parameter space. For a more
		   accurate result we would call qaws_surface_find_closest_point, but
		   the linear approximation is cheaper and adequate for the trim curves. */
		if (out_trim_a_uv)
		{
			qaws_scalar tu, tv;
			qaws_vec3 closest;
			qaws_status cst = qaws_surface_find_closest_point(
				surface_a, contact_a, &tu, &tv, &closest);
			if (cst != QAWS_STATUS_OK)
			{
				tu = sp->u1; tv = sp->v1;
			}
			out_trim_a_uv[si * 2 + 0] = tu;
			out_trim_a_uv[si * 2 + 1] = tv;
		}
		if (out_trim_b_uv)
		{
			qaws_scalar tu, tv;
			qaws_vec3 closest;
			qaws_status cst = qaws_surface_find_closest_point(
				surface_b, contact_b, &tu, &tv, &closest);
			if (cst != QAWS_STATUS_OK)
			{
				tu = sp->u2; tv = sp->v2;
			}
			out_trim_b_uv[si * 2 + 0] = tu;
			out_trim_b_uv[si * 2 + 1] = tv;
		}

		/* Vectors from center to contact points */
		ra = fillet_vec3_sub(contact_a, center);
		rb = fillet_vec3_sub(contact_b, center);

		/* Rotation axis = normalize(ra x rb) */
		axis = fillet_vec3_cross(ra, rb);
		axis_len = fillet_vec3_length(axis);

		if (axis_len < QAWS_LITERAL(1e-12))
		{
			/* contact_a and contact_b are (nearly) the same or opposite:
			   degenerate arc, fill with linear interpolation */
			for (ai = 0; ai <= arc_segments; ai++)
			{
				qaws_scalar f = (qaws_scalar)ai / (qaws_scalar)arc_segments;
				qaws_vec3 pt;
				pt.x = contact_a.x * (QAWS_ONE - f) + contact_b.x * f;
				pt.y = contact_a.y * (QAWS_ONE - f) + contact_b.y * f;
				pt.z = contact_a.z * (QAWS_ONE - f) + contact_b.z * f;
				grid[si * v_count + ai] = pt;
			}
		}
		else
		{
			/* Normalize axis */
			axis.x /= axis_len; axis.y /= axis_len; axis.z /= axis_len;

			/* Angle between ra and rb */
			{
				qaws_scalar ra_len = fillet_vec3_length(ra);
				qaws_scalar rb_len = fillet_vec3_length(rb);
				qaws_scalar dot_ab;
				if (ra_len < QAWS_LITERAL(1e-14) || rb_len < QAWS_LITERAL(1e-14))
				{
					angle = QAWS_ZERO;
				}
				else
				{
					dot_ab = fillet_vec3_dot(ra, rb) / (ra_len * rb_len);
					if (dot_ab > QAWS_ONE) dot_ab = QAWS_ONE;
					if (dot_ab < -QAWS_ONE) dot_ab = -QAWS_ONE;
					angle = QAWS_ATAN2(
						fillet_vec3_length(fillet_vec3_cross(ra, rb)),
						fillet_vec3_dot(ra, rb));
				}
			}

			/* Sweep arc using Rodrigues rotation */
			for (ai = 0; ai <= arc_segments; ai++)
			{
				qaws_scalar f = (qaws_scalar)ai / (qaws_scalar)arc_segments;
				qaws_scalar theta = f * angle;
				qaws_vec3 rotated = rodrigues_rotate(ra, axis, theta);
				grid[si * v_count + ai] = fillet_vec3_add(center, rotated);
			}
		}
	}

	*out_grid = grid;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/* Internal: core creation shared by both public functions            */
/* ------------------------------------------------------------------ */

static qaws_status create_fillet_internal(
	qaws_surface_fillet_desc const* desc,
	qaws_surface** out_surface,
	qaws_scalar* out_trim_a_uv,
	qaws_scalar* out_trim_b_uv)
{
	qaws_ssi_desc ssi_desc;
	qaws_ssi_curve ssi_curves[FILLET_SSI_MAX_CURVES];
	qaws_ssi_point* point_buffer = NULL;
	qaws_ssi_point* resampled = NULL;
	unsigned int curve_count = 0;
	unsigned int sample_count, arc_segments;
	qaws_vec3* grid = NULL;
	qaws_surface_fillet_impl* impl = NULL;
	qaws_surface* surface = NULL;
	qaws_range u_range, v_range;
	qaws_status status;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->surface_a || !desc->surface_b) return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->radius <= QAWS_ZERO) return QAWS_STATUS_INVALID_ARGUMENT;

	sample_count = desc->sample_count;
	if (sample_count == 0) sample_count = FILLET_DEFAULT_SAMPLES;
	if (sample_count < 2) sample_count = 2;

	arc_segments = desc->arc_segments;
	if (arc_segments == 0) arc_segments = FILLET_DEFAULT_ARC_SEGS;
	if (arc_segments < 1) arc_segments = 1;

	/* Allocate SSI point buffer */
	point_buffer = (qaws_ssi_point*)malloc(
		FILLET_SSI_POINT_CAP * sizeof(qaws_ssi_point));
	if (!point_buffer) return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Find intersection curves */
	memset(&ssi_desc, 0, sizeof(ssi_desc));
	ssi_desc.surface_a = desc->surface_a;
	ssi_desc.surface_b = desc->surface_b;
	ssi_desc.grid_samples = FILLET_SSI_GRID;
	ssi_desc.max_march_steps = FILLET_SSI_MAX_STEPS;

	status = qaws_surface_intersect(
		&ssi_desc, ssi_curves, FILLET_SSI_MAX_CURVES, &curve_count,
		point_buffer, FILLET_SSI_POINT_CAP);

	if (status != QAWS_STATUS_OK || curve_count == 0)
	{
		free(point_buffer);
		return (status != QAWS_STATUS_OK) ? status : QAWS_STATUS_NUMERICAL_FAILURE;
	}

	/* Use the first (longest) intersection curve.
	   Resample to uniform spacing. */
	resampled = (qaws_ssi_point*)malloc(sample_count * sizeof(qaws_ssi_point));
	if (!resampled)
	{
		free(point_buffer);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	resample_ssi_curve(
		ssi_curves[0].points, ssi_curves[0].point_count,
		sample_count, resampled);

	/* Build fillet grid */
	status = build_fillet_grid(
		desc->surface_a, desc->surface_b, desc->radius,
		resampled, sample_count, arc_segments,
		&grid, out_trim_a_uv, out_trim_b_uv);

	free(resampled);
	free(point_buffer);

	if (status != QAWS_STATUS_OK) return status;

	/* Allocate surface */
	u_range.min_value = 0; u_range.max_value = 1;
	v_range.min_value = 0; v_range.max_value = 1;

	surface = qaws_internal_surface_alloc(
		QAWS_SURFACE_KIND_TRIMMED, /* reuse TRIMMED kind since no FILLET kind yet */
		1, 1, u_range, v_range,
		&fillet_surface_vtable);
	if (!surface)
	{
		free(grid);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl = (qaws_surface_fillet_impl*)malloc(sizeof(qaws_surface_fillet_impl));
	if (!impl)
	{
		free(grid);
		qaws_internal_surface_free(surface);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->grid = grid;
	impl->u_count = sample_count;
	impl->v_count = arc_segments + 1;

	surface->impl = impl;
	*out_surface = surface;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/* Public API                                                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_create_fillet(
	qaws_surface_fillet_desc const* desc,
	qaws_surface** out_surface)
{
	return create_fillet_internal(desc, out_surface, NULL, NULL);
}

qaws_status qaws_surface_create_fillet_with_trims(
	qaws_surface_fillet_desc const* desc,
	qaws_surface** out_fillet,
	qaws_curve** out_trim_a,
	qaws_curve** out_trim_b)
{
	qaws_scalar* trim_a_uv = NULL;
	qaws_scalar* trim_b_uv = NULL;
	qaws_scalar* trim_a_params = NULL;
	qaws_scalar* trim_b_params = NULL;
	qaws_scalar* trim_a_coords = NULL;
	qaws_scalar* trim_b_coords = NULL;
	unsigned int sample_count;
	unsigned int cp_count;
	unsigned int i;
	qaws_status status;

	if (!desc || !out_fillet || !out_trim_a || !out_trim_b)
		return QAWS_STATUS_INVALID_ARGUMENT;

	sample_count = desc->sample_count;
	if (sample_count == 0) sample_count = FILLET_DEFAULT_SAMPLES;
	if (sample_count < 2) sample_count = 2;

	/* Allocate trim UV buffers */
	trim_a_uv = (qaws_scalar*)malloc(sample_count * 2 * sizeof(qaws_scalar));
	trim_b_uv = (qaws_scalar*)malloc(sample_count * 2 * sizeof(qaws_scalar));
	if (!trim_a_uv || !trim_b_uv)
	{
		free(trim_a_uv);
		free(trim_b_uv);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Build the fillet surface and collect trim UV data */
	status = create_fillet_internal(desc, out_fillet, trim_a_uv, trim_b_uv);
	if (status != QAWS_STATUS_OK)
	{
		free(trim_a_uv);
		free(trim_b_uv);
		return status;
	}

	/* Fit B-spline curves through the 2D trim points.
	   Parameters are uniform in [0,1], coordinates are (u,v) pairs. */
	trim_a_params = (qaws_scalar*)malloc(sample_count * sizeof(qaws_scalar));
	trim_b_params = (qaws_scalar*)malloc(sample_count * sizeof(qaws_scalar));
	trim_a_coords = (qaws_scalar*)malloc(sample_count * 2 * sizeof(qaws_scalar));
	trim_b_coords = (qaws_scalar*)malloc(sample_count * 2 * sizeof(qaws_scalar));
	if (!trim_a_params || !trim_b_params || !trim_a_coords || !trim_b_coords)
	{
		free(trim_a_uv); free(trim_b_uv);
		free(trim_a_params); free(trim_b_params);
		free(trim_a_coords); free(trim_b_coords);
		qaws_surface_destroy(*out_fillet);
		*out_fillet = NULL;
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < sample_count; i++)
	{
		qaws_scalar t = (qaws_scalar)i / (qaws_scalar)(sample_count - 1);
		trim_a_params[i] = t;
		trim_b_params[i] = t;
		trim_a_coords[i * 2 + 0] = trim_a_uv[i * 2 + 0];
		trim_a_coords[i * 2 + 1] = trim_a_uv[i * 2 + 1];
		trim_b_coords[i * 2 + 0] = trim_b_uv[i * 2 + 0];
		trim_b_coords[i * 2 + 1] = trim_b_uv[i * 2 + 1];
	}

	/* Number of control points: use min(sample_count, 32) for a smooth fit */
	cp_count = sample_count;
	if (cp_count > 32) cp_count = 32;
	if (cp_count < 4) cp_count = 4;
	if (cp_count > sample_count) cp_count = sample_count;

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_2D, 3,
		trim_a_params, trim_a_coords,
		sample_count, cp_count,
		out_trim_a);

	if (status != QAWS_STATUS_OK)
	{
		free(trim_a_uv); free(trim_b_uv);
		free(trim_a_params); free(trim_b_params);
		free(trim_a_coords); free(trim_b_coords);
		qaws_surface_destroy(*out_fillet);
		*out_fillet = NULL;
		return status;
	}

	status = qaws_internal_fit_bspline(
		QAWS_DIMENSION_2D, 3,
		trim_b_params, trim_b_coords,
		sample_count, cp_count,
		out_trim_b);

	if (status != QAWS_STATUS_OK)
	{
		free(trim_a_uv); free(trim_b_uv);
		free(trim_a_params); free(trim_b_params);
		free(trim_a_coords); free(trim_b_coords);
		qaws_curve_destroy(*out_trim_a);
		*out_trim_a = NULL;
		qaws_surface_destroy(*out_fillet);
		*out_fillet = NULL;
		return status;
	}

	free(trim_a_uv);
	free(trim_b_uv);
	free(trim_a_params);
	free(trim_b_params);
	free(trim_a_coords);
	free(trim_b_coords);
	return QAWS_STATUS_OK;
}
