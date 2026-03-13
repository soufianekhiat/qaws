#include "qaws_surface_intersect.h"
#include "qaws_surface.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* Default parameters */
#define SSI_DEFAULT_TOLERANCE    QAWS_LITERAL(1e-6)
#define SSI_DEFAULT_GRID         20
#define SSI_DEFAULT_MAX_STEPS    500
#define SSI_NEWTON_MAX_ITER      10
#define SSI_SEED_DIST_THRESHOLD  QAWS_LITERAL(0.05)
#define SSI_MERGE_THRESHOLD      QAWS_LITERAL(1e-4)

/* Helper: 3D vector length */
static qaws_scalar vec3_length(qaws_vec3 v)
{
	return QAWS_SQRT(v.x * v.x + v.y * v.y + v.z * v.z);
}

/* Helper: 3D vector distance */
static qaws_scalar vec3_dist(qaws_vec3 a, qaws_vec3 b)
{
	qaws_scalar dx = a.x - b.x;
	qaws_scalar dy = a.y - b.y;
	qaws_scalar dz = a.z - b.z;
	return QAWS_SQRT(dx * dx + dy * dy + dz * dz);
}

/* Helper: cross product */
static qaws_vec3 vec3_cross(qaws_vec3 a, qaws_vec3 b)
{
	qaws_vec3 r;
	r.x = a.y * b.z - a.z * b.y;
	r.y = a.z * b.x - a.x * b.z;
	r.z = a.x * b.y - a.y * b.x;
	return r;
}

/* Helper: dot product */
static qaws_scalar vec3_dot(qaws_vec3 a, qaws_vec3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* Helper: clamp to [0,1] */
static qaws_scalar clamp01(qaws_scalar x)
{
	if (x < QAWS_ZERO) return QAWS_ZERO;
	if (x > QAWS_ONE) return QAWS_ONE;
	return x;
}

/* Solve 2x2 linear system: [a b; c d] * [x; y] = [e; f]
   Returns 0 on singular, 1 on success. */
static int solve_2x2(
	qaws_scalar a, qaws_scalar b,
	qaws_scalar c, qaws_scalar d,
	qaws_scalar e, qaws_scalar f,
	qaws_scalar* out_x, qaws_scalar* out_y)
{
	qaws_scalar det = a * d - b * c;
	if (QAWS_FABS(det) < QAWS_LITERAL(1e-30))
		return 0;
	*out_x = (e * d - b * f) / det;
	*out_y = (a * f - e * c) / det;
	return 1;
}

/* Newton refinement: given (u1,v1) on surface_a and (u2,v2) on surface_b,
   refine so that S1(u1,v1) = S2(u2,v2).
   Uses midpoint projection approach: each surface projects independently. */
static int newton_refine(
	qaws_surface const* sa,
	qaws_surface const* sb,
	qaws_scalar* u1, qaws_scalar* v1,
	qaws_scalar* u2, qaws_scalar* v2,
	qaws_vec3* out_pos,
	qaws_scalar tolerance)
{
	int iter;
	unsigned int eval_flags = QAWS_SURFACE_EVAL_POSITION
		| QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;

	for (iter = 0; iter < SSI_NEWTON_MAX_ITER; ++iter)
	{
		qaws_surface_eval_result r1, r2;
		qaws_vec3 target;
		qaws_vec3 diff1, diff2;
		qaws_scalar du1_val, dv1_val, du2_val, dv2_val;

		/* Evaluate both surfaces */
		memset(&r1, 0, sizeof(r1));
		memset(&r2, 0, sizeof(r2));
		if (qaws_surface_evaluate(sa, *u1, *v1, eval_flags, &r1) != QAWS_STATUS_OK)
			return 0;
		if (qaws_surface_evaluate(sb, *u2, *v2, eval_flags, &r2) != QAWS_STATUS_OK)
			return 0;

		/* Check convergence */
		{
			qaws_scalar dist = vec3_dist(r1.position, r2.position);
			if (dist < tolerance)
			{
				out_pos->x = (r1.position.x + r2.position.x) * QAWS_LITERAL(0.5);
				out_pos->y = (r1.position.y + r2.position.y) * QAWS_LITERAL(0.5);
				out_pos->z = (r1.position.z + r2.position.z) * QAWS_LITERAL(0.5);
				return 1;
			}
		}

		/* Target = midpoint of current evaluations */
		target.x = (r1.position.x + r2.position.x) * QAWS_LITERAL(0.5);
		target.y = (r1.position.y + r2.position.y) * QAWS_LITERAL(0.5);
		target.z = (r1.position.z + r2.position.z) * QAWS_LITERAL(0.5);

		/* Solve for surface 1: S1 + S1u*du1 + S1v*dv1 = target
		   diff1 = target - S1 */
		diff1.x = target.x - r1.position.x;
		diff1.y = target.y - r1.position.y;
		diff1.z = target.z - r1.position.z;

		/* Project diff1 onto (S1u, S1v) using normal equations:
		   [S1u.S1u  S1u.S1v] [du1]   [S1u.diff1]
		   [S1v.S1u  S1v.S1v] [dv1] = [S1v.diff1] */
		{
			qaws_scalar a11 = vec3_dot(r1.du, r1.du);
			qaws_scalar a12 = vec3_dot(r1.du, r1.dv);
			qaws_scalar a22 = vec3_dot(r1.dv, r1.dv);
			qaws_scalar b1 = vec3_dot(r1.du, diff1);
			qaws_scalar b2 = vec3_dot(r1.dv, diff1);
			if (!solve_2x2(a11, a12, a12, a22, b1, b2, &du1_val, &dv1_val))
				return 0;
		}

		/* Solve for surface 2: S2 + S2u*du2 + S2v*dv2 = target */
		diff2.x = target.x - r2.position.x;
		diff2.y = target.y - r2.position.y;
		diff2.z = target.z - r2.position.z;

		{
			qaws_scalar a11 = vec3_dot(r2.du, r2.du);
			qaws_scalar a12 = vec3_dot(r2.du, r2.dv);
			qaws_scalar a22 = vec3_dot(r2.dv, r2.dv);
			qaws_scalar b1 = vec3_dot(r2.du, diff2);
			qaws_scalar b2 = vec3_dot(r2.dv, diff2);
			if (!solve_2x2(a11, a12, a12, a22, b1, b2, &du2_val, &dv2_val))
				return 0;
		}

		*u1 = clamp01(*u1 + du1_val);
		*v1 = clamp01(*v1 + dv1_val);
		*u2 = clamp01(*u2 + du2_val);
		*v2 = clamp01(*v2 + dv2_val);
	}

	/* Did not converge within max iterations; check final distance */
	{
		qaws_surface_eval_result r1, r2;
		memset(&r1, 0, sizeof(r1));
		memset(&r2, 0, sizeof(r2));
		qaws_surface_evaluate(sa, *u1, *v1, QAWS_SURFACE_EVAL_POSITION, &r1);
		qaws_surface_evaluate(sb, *u2, *v2, QAWS_SURFACE_EVAL_POSITION, &r2);
		if (vec3_dist(r1.position, r2.position) < tolerance * QAWS_LITERAL(10.0))
		{
			out_pos->x = (r1.position.x + r2.position.x) * QAWS_LITERAL(0.5);
			out_pos->y = (r1.position.y + r2.position.y) * QAWS_LITERAL(0.5);
			out_pos->z = (r1.position.z + r2.position.z) * QAWS_LITERAL(0.5);
			return 1;
		}
	}
	return 0;
}

/* Check if a seed point is a duplicate of an existing seed */
static int is_duplicate_seed(
	qaws_ssi_point const* seeds, unsigned int count,
	qaws_scalar u1, qaws_scalar v1,
	qaws_scalar u2, qaws_scalar v2,
	qaws_scalar threshold)
{
	unsigned int i;
	for (i = 0; i < count; ++i)
	{
		qaws_scalar du1 = seeds[i].u1 - u1;
		qaws_scalar dv1 = seeds[i].v1 - v1;
		qaws_scalar du2 = seeds[i].u2 - u2;
		qaws_scalar dv2 = seeds[i].v2 - v2;
		qaws_scalar d = QAWS_SQRT(du1 * du1 + dv1 * dv1 + du2 * du2 + dv2 * dv2);
		if (d < threshold) return 1;
	}
	return 0;
}

/* March along the intersection curve in one direction from a seed point.
   direction: +1 or -1.
   Returns number of points written. */
static unsigned int march_direction(
	qaws_surface const* sa,
	qaws_surface const* sb,
	qaws_scalar u1_start, qaws_scalar v1_start,
	qaws_scalar u2_start, qaws_scalar v2_start,
	qaws_vec3 pos_start,
	int direction,
	qaws_scalar step_size,
	qaws_scalar tolerance,
	unsigned int max_steps,
	qaws_ssi_point* out_points,
	unsigned int point_capacity)
{
	unsigned int count = 0;
	qaws_scalar u1 = u1_start;
	qaws_scalar v1 = v1_start;
	qaws_scalar u2 = u2_start;
	qaws_scalar v2 = v2_start;
	qaws_vec3 pos = pos_start;
	unsigned int step;
	unsigned int boundary_streak = 0;

	for (step = 0; step < max_steps && count < point_capacity; ++step)
	{
		qaws_surface_eval_result r1, r2;
		qaws_vec3 T, n1, n2;
		qaws_scalar t_len, actual_step;
		qaws_scalar pred_u1 = u1, pred_v1 = v1, pred_u2 = u2, pred_v2 = v2;
		qaws_vec3 new_pos = pos;
		int on_boundary;
		unsigned int eval_flags = QAWS_SURFACE_EVAL_POSITION
			| QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV
			| QAWS_SURFACE_EVAL_NORMAL;

		/* Evaluate normals at current point */
		memset(&r1, 0, sizeof(r1));
		memset(&r2, 0, sizeof(r2));
		if (qaws_surface_evaluate(sa, u1, v1, eval_flags, &r1) != QAWS_STATUS_OK)
			break;
		if (qaws_surface_evaluate(sb, u2, v2, eval_flags, &r2) != QAWS_STATUS_OK)
			break;

		n1 = r1.normal;
		n2 = r2.normal;

		/* Intersection tangent: T = N1 x N2 */
		T = vec3_cross(n1, n2);
		t_len = vec3_length(T);
		if (t_len < QAWS_LITERAL(1e-14))
			break; /* surfaces truly tangent; can't march */

		/* Normalize and apply direction */
		{
			qaws_scalar inv_len = (qaws_scalar)direction / t_len;
			T.x *= inv_len;
			T.y *= inv_len;
			T.z *= inv_len;
		}

		/* Adaptive step: reduce near tangency where the intersection
		   curve has high curvature and needs finer tracking. */
		actual_step = step_size;
		if (t_len < QAWS_LITERAL(0.1))
		{
			actual_step = step_size * (t_len / QAWS_LITERAL(0.1));
			if (actual_step < step_size * QAWS_LITERAL(0.02))
				actual_step = step_size * QAWS_LITERAL(0.02);
		}

		/* Predict + project + Newton-refine, with retry on failure.
		   If Newton doesn't converge, halve the step and try again. */
		{
			int converged = 0;
			int retry;
			qaws_scalar try_step = actual_step;

			for (retry = 0; retry < 4 && !converged; ++retry)
			{
				qaws_vec3 pred_pos, dp;
				qaws_scalar a11, a12, a22, b1_val, b2_val;
				qaws_scalar dparam_u, dparam_v;

				pred_pos.x = pos.x + try_step * T.x;
				pred_pos.y = pos.y + try_step * T.y;
				pred_pos.z = pos.z + try_step * T.z;

				/* Project onto surface 1 tangent plane */
				dp.x = pred_pos.x - r1.position.x;
				dp.y = pred_pos.y - r1.position.y;
				dp.z = pred_pos.z - r1.position.z;
				a11 = vec3_dot(r1.du, r1.du);
				a12 = vec3_dot(r1.du, r1.dv);
				a22 = vec3_dot(r1.dv, r1.dv);
				b1_val = vec3_dot(r1.du, dp);
				b2_val = vec3_dot(r1.dv, dp);
				if (!solve_2x2(a11, a12, a12, a22, b1_val, b2_val,
					&dparam_u, &dparam_v))
				{
					try_step *= QAWS_LITERAL(0.5);
					continue;
				}
				pred_u1 = clamp01(u1 + dparam_u);
				pred_v1 = clamp01(v1 + dparam_v);

				/* Project onto surface 2 tangent plane */
				dp.x = pred_pos.x - r2.position.x;
				dp.y = pred_pos.y - r2.position.y;
				dp.z = pred_pos.z - r2.position.z;
				a11 = vec3_dot(r2.du, r2.du);
				a12 = vec3_dot(r2.du, r2.dv);
				a22 = vec3_dot(r2.dv, r2.dv);
				b1_val = vec3_dot(r2.du, dp);
				b2_val = vec3_dot(r2.dv, dp);
				if (!solve_2x2(a11, a12, a12, a22, b1_val, b2_val,
					&dparam_u, &dparam_v))
				{
					try_step *= QAWS_LITERAL(0.5);
					continue;
				}
				pred_u2 = clamp01(u2 + dparam_u);
				pred_v2 = clamp01(v2 + dparam_v);

				/* Newton-refine */
				if (newton_refine(sa, sb, &pred_u1, &pred_v1,
					&pred_u2, &pred_v2, &new_pos, tolerance))
					converged = 1;
				else
					try_step *= QAWS_LITERAL(0.5);
			}
			if (!converged)
				break;
		}

		/* Check if any parameter is stuck at the domain boundary.
		   Use a streak counter: near-tangent intersection curves may
		   run along a boundary for several steps before turning back. */
		on_boundary = 0;
		if ((pred_u1 <= QAWS_LITERAL(1e-8) && u1 <= QAWS_LITERAL(1e-8)) ||
			(pred_u1 >= QAWS_ONE - QAWS_LITERAL(1e-8) && u1 >= QAWS_ONE - QAWS_LITERAL(1e-8)) ||
			(pred_v1 <= QAWS_LITERAL(1e-8) && v1 <= QAWS_LITERAL(1e-8)) ||
			(pred_v1 >= QAWS_ONE - QAWS_LITERAL(1e-8) && v1 >= QAWS_ONE - QAWS_LITERAL(1e-8)) ||
			(pred_u2 <= QAWS_LITERAL(1e-8) && u2 <= QAWS_LITERAL(1e-8)) ||
			(pred_u2 >= QAWS_ONE - QAWS_LITERAL(1e-8) && u2 >= QAWS_ONE - QAWS_LITERAL(1e-8)) ||
			(pred_v2 <= QAWS_LITERAL(1e-8) && v2 <= QAWS_LITERAL(1e-8)) ||
			(pred_v2 >= QAWS_ONE - QAWS_LITERAL(1e-8) && v2 >= QAWS_ONE - QAWS_LITERAL(1e-8)))
			on_boundary = 1;

		if (on_boundary)
		{
			if (++boundary_streak > 5)
				break;
		}
		else
		{
			boundary_streak = 0;
		}

		/* Store point */
		out_points[count].u1 = pred_u1;
		out_points[count].v1 = pred_v1;
		out_points[count].u2 = pred_u2;
		out_points[count].v2 = pred_v2;
		out_points[count].position = new_pos;
		++count;

		/* Update for next step */
		u1 = pred_u1;
		v1 = pred_v1;
		u2 = pred_u2;
		v2 = pred_v2;
		pos = new_pos;

		/* Check if we've returned to start (closed curve) */
		if (count > 2)
		{
			qaws_scalar d_ret = vec3_dist(pos, pos_start);
			if (d_ret < step_size * QAWS_LITERAL(1.5))
			{
				qaws_scalar du1_ret = u1 - u1_start;
				qaws_scalar dv1_ret = v1 - v1_start;
				qaws_scalar dp_ret = QAWS_SQRT(du1_ret * du1_ret + dv1_ret * dv1_ret);
				if (dp_ret < step_size * QAWS_LITERAL(2.0))
					break;
			}
		}
	}

	return count;
}

qaws_status qaws_surface_intersect(
	qaws_ssi_desc const* desc,
	qaws_ssi_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* point_buffer,
	unsigned int point_capacity)
{
	qaws_surface const* sa;
	qaws_surface const* sb;
	qaws_scalar tolerance;
	unsigned int grid_n;
	unsigned int max_march;
	qaws_scalar step;
	unsigned int seed_capacity;
	qaws_ssi_point* seeds;
	unsigned int seed_count = 0;
	unsigned int total_points = 0;
	unsigned int curve_count = 0;
	unsigned int gi, gj;
	int* seed_used;
	unsigned int si;

	if (!desc || !out_curves || !out_curve_count || !point_buffer)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->surface_a || !desc->surface_b)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve_capacity == 0 || point_capacity == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	sa = desc->surface_a;
	sb = desc->surface_b;
	tolerance = (desc->tolerance > QAWS_ZERO) ? desc->tolerance : SSI_DEFAULT_TOLERANCE;
	grid_n = (desc->grid_samples > 0) ? desc->grid_samples : SSI_DEFAULT_GRID;
	max_march = (desc->max_march_steps > 0) ? desc->max_march_steps : SSI_DEFAULT_MAX_STEPS;
	step = desc->step_size;

	/* Auto step size: approximate diagonal / grid_samples */
	if (step <= QAWS_ZERO)
		step = QAWS_LITERAL(1.0) / (qaws_scalar)grid_n;

	/* Phase 1: Find seed points via grid sampling.
	   Compute dynamic proximity threshold from surface extents. */
	{
		qaws_scalar seed_threshold;
		qaws_surface_eval_result rc0, rc1;
		qaws_scalar diag_a, diag_b, max_diag;

		memset(&rc0, 0, sizeof(rc0));
		memset(&rc1, 0, sizeof(rc1));
		qaws_surface_evaluate(sa, 0, 0, QAWS_SURFACE_EVAL_POSITION, &rc0);
		qaws_surface_evaluate(sa, 1, 1, QAWS_SURFACE_EVAL_POSITION, &rc1);
		diag_a = vec3_dist(rc0.position, rc1.position);

		memset(&rc0, 0, sizeof(rc0));
		memset(&rc1, 0, sizeof(rc1));
		qaws_surface_evaluate(sb, 0, 0, QAWS_SURFACE_EVAL_POSITION, &rc0);
		qaws_surface_evaluate(sb, 1, 1, QAWS_SURFACE_EVAL_POSITION, &rc1);
		diag_b = vec3_dist(rc0.position, rc1.position);

		max_diag = (diag_a > diag_b) ? diag_a : diag_b;
		if (max_diag < QAWS_LITERAL(1.0)) max_diag = QAWS_LITERAL(1.0);
		seed_threshold = max_diag / (qaws_scalar)grid_n;

		seed_capacity = grid_n * grid_n;
		if (seed_capacity > 1000) seed_capacity = 1000;
		seeds = (qaws_ssi_point*)malloc(seed_capacity * sizeof(qaws_ssi_point));
		if (!seeds) return QAWS_STATUS_ALLOCATION_FAILURE;

		for (gi = 0; gi < grid_n && seed_count < seed_capacity; ++gi)
		{
			for (gj = 0; gj < grid_n && seed_count < seed_capacity; ++gj)
			{
				qaws_scalar u1_c = ((qaws_scalar)gi + QAWS_LITERAL(0.5)) / (qaws_scalar)grid_n;
				qaws_scalar v1_c = ((qaws_scalar)gj + QAWS_LITERAL(0.5)) / (qaws_scalar)grid_n;
				qaws_surface_eval_result r1;
				unsigned int gi2, gj2;

				memset(&r1, 0, sizeof(r1));
				if (qaws_surface_evaluate(sa, u1_c, v1_c, QAWS_SURFACE_EVAL_POSITION, &r1)
					!= QAWS_STATUS_OK)
					continue;

				for (gi2 = 0; gi2 < grid_n && seed_count < seed_capacity; ++gi2)
				{
					for (gj2 = 0; gj2 < grid_n && seed_count < seed_capacity; ++gj2)
					{
						qaws_scalar u2_c = ((qaws_scalar)gi2 + QAWS_LITERAL(0.5)) / (qaws_scalar)grid_n;
						qaws_scalar v2_c = ((qaws_scalar)gj2 + QAWS_LITERAL(0.5)) / (qaws_scalar)grid_n;
						qaws_surface_eval_result r2;
						qaws_scalar dist;

						memset(&r2, 0, sizeof(r2));
						if (qaws_surface_evaluate(sb, u2_c, v2_c, QAWS_SURFACE_EVAL_POSITION, &r2)
							!= QAWS_STATUS_OK)
							continue;

						dist = vec3_dist(r1.position, r2.position);
						if (dist < seed_threshold)
						{
							qaws_scalar su1 = u1_c, sv1 = v1_c, su2 = u2_c, sv2 = v2_c;
							qaws_vec3 seed_pos;

							if (newton_refine(sa, sb, &su1, &sv1, &su2, &sv2, &seed_pos, tolerance))
							{
								if (!is_duplicate_seed(seeds, seed_count, su1, sv1, su2, sv2,
									step * QAWS_LITERAL(0.5)))
								{
									seeds[seed_count].u1 = su1;
									seeds[seed_count].v1 = sv1;
									seeds[seed_count].u2 = su2;
									seeds[seed_count].v2 = sv2;
									seeds[seed_count].position = seed_pos;
									++seed_count;
								}
							}
						}
					}
				}
			}
		}
	}

	if (seed_count == 0)
	{
		free(seeds);
		*out_curve_count = 0;
		return QAWS_STATUS_OK;
	}

	/* Phase 2: March from each unused seed to form curves */
	seed_used = (int*)calloc(seed_count, sizeof(int));
	if (!seed_used)
	{
		free(seeds);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (si = 0; si < seed_count && curve_count < curve_capacity; ++si)
	{
		unsigned int fwd_count, bwd_count, total_curve_pts;
		unsigned int remaining;
		qaws_ssi_point* fwd_buf;
		qaws_ssi_point* bwd_buf;
		unsigned int fwd_cap, bwd_cap;
		unsigned int pi_idx;

		if (seed_used[si]) continue;
		seed_used[si] = 1;

		remaining = point_capacity - total_points;
		if (remaining < 3) break;

		/* Allocate temporary buffers for forward and backward marching */
		fwd_cap = remaining / 2;
		bwd_cap = remaining - fwd_cap;
		if (fwd_cap > max_march) fwd_cap = max_march;
		if (bwd_cap > max_march) bwd_cap = max_march;

		fwd_buf = (qaws_ssi_point*)malloc(fwd_cap * sizeof(qaws_ssi_point));
		bwd_buf = (qaws_ssi_point*)malloc(bwd_cap * sizeof(qaws_ssi_point));
		if (!fwd_buf || !bwd_buf)
		{
			free(fwd_buf);
			free(bwd_buf);
			break;
		}

		/* March forward */
		fwd_count = march_direction(sa, sb,
			seeds[si].u1, seeds[si].v1,
			seeds[si].u2, seeds[si].v2,
			seeds[si].position,
			1, step, tolerance, max_march,
			fwd_buf, fwd_cap);

		/* March backward */
		bwd_count = march_direction(sa, sb,
			seeds[si].u1, seeds[si].v1,
			seeds[si].u2, seeds[si].v2,
			seeds[si].position,
			-1, step, tolerance, max_march,
			bwd_buf, bwd_cap);

		/* Assemble curve: backward (reversed) + seed + forward */
		total_curve_pts = bwd_count + 1 + fwd_count;
		if (total_points + total_curve_pts > point_capacity)
		{
			/* Truncate */
			if (total_points + 1 > point_capacity)
			{
				free(fwd_buf);
				free(bwd_buf);
				break;
			}
			total_curve_pts = point_capacity - total_points;
		}

		/* Write backward points in reverse order */
		{
			unsigned int write_bwd = bwd_count;
			if (write_bwd > total_curve_pts - 1) write_bwd = total_curve_pts - 1;
			for (pi_idx = 0; pi_idx < write_bwd; ++pi_idx)
			{
				point_buffer[total_points + pi_idx] = bwd_buf[write_bwd - 1 - pi_idx];
			}

			/* Write seed point */
			if (write_bwd < total_curve_pts)
			{
				point_buffer[total_points + write_bwd] = seeds[si];
			}

			/* Write forward points */
			{
				unsigned int fwd_start = write_bwd + 1;
				unsigned int write_fwd = 0;
				if (fwd_start < total_curve_pts)
					write_fwd = total_curve_pts - fwd_start;
				if (write_fwd > fwd_count) write_fwd = fwd_count;
				for (pi_idx = 0; pi_idx < write_fwd; ++pi_idx)
				{
					point_buffer[total_points + fwd_start + pi_idx] = fwd_buf[pi_idx];
				}
				total_curve_pts = fwd_start + write_fwd;
			}
		}

		/* Mark nearby seeds as used */
		{
			unsigned int sj;
			for (sj = si + 1; sj < seed_count; ++sj)
			{
				unsigned int pk;
				if (seed_used[sj]) continue;
				for (pk = 0; pk < total_curve_pts; ++pk)
				{
					qaws_scalar du1 = seeds[sj].u1 - point_buffer[total_points + pk].u1;
					qaws_scalar dv1 = seeds[sj].v1 - point_buffer[total_points + pk].v1;
					qaws_scalar d = QAWS_SQRT(du1 * du1 + dv1 * dv1);
					if (d < step * QAWS_LITERAL(2.0))
					{
						seed_used[sj] = 1;
						break;
					}
				}
			}
		}

		/* Record curve */
		if (total_curve_pts > 0)
		{
			out_curves[curve_count].points = &point_buffer[total_points];
			out_curves[curve_count].point_count = total_curve_pts;
			++curve_count;
			total_points += total_curve_pts;
		}

		free(fwd_buf);
		free(bwd_buf);
	}

	free(seed_used);
	free(seeds);

	*out_curve_count = curve_count;
	return QAWS_STATUS_OK;
}
