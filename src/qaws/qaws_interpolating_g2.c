/*
 * Curvature-continuous cubic spline through points, solved for the tangent angle at each knot.
 *
 * Follows aurimasg/splines:
 *
 * MIT License
 *
 * Copyright (c) 2026 Aurimas Gasiulis
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * Each segment is a cubic set by two angles, alpha at its start and beta at its end, measured against its chord
 * (https://spline.technology/paper1.pdf). The knot angles start from a local guess moved to the middle of each knot's
 * rotation range (past its ends a segment's arm reaches zero length and the curve flips). Sweeps then solve one knot at
 * a time with the neighbouring knots fixed, by the ITP root finder, until no knot moves by more than a tolerance. Newton
 * iterations on the tridiagonal system of curvature differences finish the solution. A Newton phase that diverges or
 * runs out of iterations restores its starting angles and tightens the sweep tolerance, within the limit on all sweeps.
 */

#include "internal/qaws_internal_bezier_chain.h"
#include "internal/qaws_internal_memory.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define G2_PI 3.14159265358979323846
#define G2_ITP_HALVINGS 32U                /* the root finder's result is within the range / 2^32 of the root */
#define G2_ITP_TRUNCATION 0.2              /* weight of the truncation step, relative to the range's width */
#define G2_MAX_CENTERING 3U
#define G2_MAX_SWEEPS 20U                  /* all sweeps of one solve */
#define G2_MAX_NEWTON 12U
#define G2_SWEEP_ERROR 1e-2                /* a sweep moving no knot by this much (radians) ends the sweeps */
#define G2_FALLBACK_FACTOR 0.01            /* on the sweep tolerance after a failed Newton phase */
#define G2_NEWTON_ERROR 1e-15

typedef struct g2_range
{
	double lo, hi;
} g2_range;

/* Curvatures at both ends of one segment (up to a factor that cancels in the difference at a knot) and their partial
   derivatives by the segment's two angles. */
typedef struct g2_curvatures
{
	double start, end, start_by_alpha, start_by_beta, end_by_alpha, end_by_beta;
} g2_curvatures;

typedef struct g2_solver
{
	double* lengths;                   /* per segment */
	double* angles;                    /* per segment, absolute, in -pi..pi */
	double* knot_angles;               /* per knot, absolute */
	double* saved;                     /* the knot angles at the start of a Newton phase */
	double* a;                         /* tridiagonal system: subdiagonal, diagonal, superdiagonal, right side */
	double* b;
	double* c;
	double* x;
	double* scratch;                   /* 4 n doubles for the tridiagonal solvers */
	g2_curvatures* curvatures;         /* per segment */
	unsigned int n, segments;
} g2_solver;

/* a NaN in b is passed over and a NaN in a stays, on which the sweeps' stop depends */
static double g2_max(double a, double b)
{
	return a < b ? b : a;
}

/* a vector's angle, positive counter-clockwise with y pointing down */
static double g2_vector_angle(double x, double y)
{
	return atan2(-y, x);
}

static double g2_normalize_angle(double radians)
{
	double const r = radians * (1.0 / (G2_PI + G2_PI));
	return (G2_PI + G2_PI) * (r - round(r));
}

static double g2_absolute_to_alpha(double knot_angle, double segment_angle)
{
	return knot_angle - segment_angle;
}

static double g2_alpha_to_absolute(double alpha, double segment_angle)
{
	return alpha + segment_angle;
}

static double g2_absolute_to_beta(double knot_angle, double segment_angle)
{
	return segment_angle - knot_angle;
}

static double g2_beta_to_absolute(double beta, double segment_angle)
{
	return -beta + segment_angle;
}

/* the arm at one end of the unit cubic, as a fraction of the chord: (5 / 12) (cos(x) - 0.2 cos(3x)), which
   cos(3x) = 4 cos(x)^3 - 3 cos(x) turns into cos(x) (2 - cos(x)^2) / 3 */
static double g2_arm_length(double a, double b)
{
	double const x = a - 0.3 * sin(2.0 * b - 0.4 * sin(2.0 * b));
	double const c = cos(x);
	return c * (2.0 - c * c) / 3.0;
}

/*
 * The curvatures at both ends of a segment and their derivatives by its two angles. Each arm depends on its own angle
 * directly and on the opposite one through a turn term:
 *
 *     arm = L(own - t(opposite))
 *     t(b) = 0.3 sin(2b - 0.4 sin(2b))              t'(b) = 0.3 cos(2b - 0.4 sin(2b)) (2 - 0.8 cos(2b))
 *     L(u) = cos(u) (2 - cos(u)^2) / 3              L'(u) = -sin(u) (2 - 3 cos(u)^2) / 3
 *
 * With arms A and B, start curvature ~ (B sin(alpha + beta) - sin(alpha)) / (A^2 length) and
 * end curvature ~ (A sin(alpha + beta) - sin(beta)) / (B^2 length). The derivative of n / (arm^2 length) by an angle is
 * (n' - 2 n arm' / arm) / (arm^2 length).
 */
static g2_curvatures g2_segment_curvatures(double alpha, double beta, double length)
{
	double const turn_alpha_argument = 2.0 * alpha - 0.4 * sin(2.0 * alpha);
	double const turn_beta_argument = 2.0 * beta - 0.4 * sin(2.0 * beta);
	double const turn_alpha = 0.3 * sin(turn_alpha_argument);
	double const turn_beta = 0.3 * sin(turn_beta_argument);
	double const turn_alpha_slope = 0.3 * cos(turn_alpha_argument) * (2.0 - 0.8 * cos(2.0 * alpha));
	double const turn_beta_slope = 0.3 * cos(turn_beta_argument) * (2.0 - 0.8 * cos(2.0 * beta));

	/* an arm's derivative by its own angle is L', by the opposite angle -L' t' */
	double const start_cosine = cos(alpha - turn_beta);
	double const end_cosine = cos(beta - turn_alpha);
	double const start_arm = start_cosine * (2.0 - start_cosine * start_cosine) / 3.0;
	double const end_arm = end_cosine * (2.0 - end_cosine * end_cosine) / 3.0;
	double const start_arm_slope = -sin(alpha - turn_beta) * (2.0 - 3.0 * start_cosine * start_cosine) / 3.0;
	double const end_arm_slope = -sin(beta - turn_alpha) * (2.0 - 3.0 * end_cosine * end_cosine) / 3.0;

	double const sine_sum = sin(alpha + beta);
	double const cosine_sum = cos(alpha + beta);
	double const start_numerator = end_arm * sine_sum - sin(alpha);
	double const end_numerator = start_arm * sine_sum - sin(beta);
	double const start_scale = 1.0 / (start_arm * start_arm * length);
	double const end_scale = 1.0 / (end_arm * end_arm * length);

	double const start_quotient = 2.0 * start_numerator * start_arm_slope / start_arm;
	double const end_quotient = 2.0 * end_numerator * end_arm_slope / end_arm;

	g2_curvatures r;
	r.start = start_numerator * start_scale;
	r.end = end_numerator * end_scale;
	r.start_by_alpha = (end_arm * cosine_sum - cos(alpha) - end_arm_slope * turn_alpha_slope * sine_sum - start_quotient) * start_scale;
	r.start_by_beta = (end_arm * cosine_sum + end_arm_slope * sine_sum + start_quotient * turn_beta_slope) * start_scale;
	r.end_by_alpha = (start_arm * cosine_sum + start_arm_slope * sine_sum + end_quotient * turn_alpha_slope) * end_scale;
	r.end_by_beta = (start_arm * cosine_sum - cos(beta) - start_arm_slope * turn_beta_slope * sine_sum - end_quotient) * end_scale;
	return r;
}

/* The angle beyond which the curve flips, given the opposite angle: where the arm becomes zero. */
static double g2_rotation_limit_axis(double beta)
{
	double const beta2 = beta + beta;
	return 0.3 * sin(beta2 - (0.4 * sin(beta2)));
}

static g2_range g2_rotation_range(double prev_knot_angle, double next_knot_angle, double prev_segment_angle, double next_segment_angle)
{
	double const beta_limit = g2_beta_to_absolute(g2_rotation_limit_axis(g2_absolute_to_alpha(prev_knot_angle, prev_segment_angle)), prev_segment_angle);
	double const alpha_limit = g2_alpha_to_absolute(g2_rotation_limit_axis(g2_absolute_to_beta(next_knot_angle, next_segment_angle)), next_segment_angle);
	double const axis = (beta_limit + alpha_limit) / 2.0;
	double const range = (G2_PI - fabs(g2_normalize_angle(beta_limit - alpha_limit))) / 2.0;
	g2_range r;
	r.lo = axis - range;
	r.hi = axis + range;
	return r;
}

/* The rotation range of knot j, whose neighbours are prev and next. A segment has its start knot's index, so the
   segment before j is prev and the one after it is j. */
static g2_range g2_knot_range(g2_solver const* s, unsigned int prev, unsigned int j, unsigned int next)
{
	return g2_rotation_range(s->knot_angles[prev], s->knot_angles[next], s->angles[prev], s->angles[j]);
}

static double g2_center_knot(g2_solver const* s, unsigned int prev, unsigned int j, unsigned int next)
{
	g2_range const r = g2_knot_range(s, prev, j, next);
	return (r.lo + r.hi) * 0.5;
}

/*
 * The curvature difference at the knot between two segments, multiplied by both curvatures' denominators: n1 d0 - n0 d1
 * for curvatures n1 / d1 after the knot and n0 / d0 before it. Inside the knot's rotation range d0 and d1 are positive,
 * so this has the sign and the root of the difference without its poles at the range's ends, which keeps the root finder
 * effective. x, y and z are the absolute angles at the previous, this and the next knot; a0, h0 and a1, h1 the angles
 * and lengths of the segments before and after.
 */
static double g2_scaled_curvature_delta(double x, double y, double z, double a0, double a1, double h0, double h1)
{
	double const alpha0 = g2_absolute_to_alpha(x, a0);
	double const beta0 = g2_absolute_to_beta(y, a0);
	double const alpha1 = g2_absolute_to_alpha(y, a1);
	double const beta1 = g2_absolute_to_beta(z, a1);

	double const alpha_arm0 = g2_arm_length(alpha0, beta0);
	double const beta_arm0 = g2_arm_length(beta0, alpha0);
	double const alpha_arm1 = g2_arm_length(alpha1, beta1);
	double const beta_arm1 = g2_arm_length(beta1, alpha1);

	double const numerator1 = (beta_arm1 * h1) * sin(alpha1 + beta1) - h1 * sin(alpha1);
	double const root1 = alpha_arm1 * h1;
	double const denominator1 = root1 * root1;

	double const numerator0 = (alpha_arm0 * h0) * sin(alpha0 + beta0) - h0 * sin(beta0);
	double const root0 = beta_arm0 * h0;
	double const denominator0 = root0 * root0;

	return numerator1 * denominator0 - numerator0 * denominator1;
}

/*
 * The angle at knot j at which the curvatures of its two segments match, its neighbours fixed, by the ITP method of
 * Oliveira and Takahashi (2020). The scaled difference decreases strictly over the rotation range and has poles of
 * opposite sign at its ends, so the range brackets one root. When it does not change sign, the range's upper end.
 */
static double g2_solve_knot_angle(g2_solver const* s, unsigned int prev, unsigned int j, unsigned int next)
{
	double const prev_knot = s->knot_angles[prev];
	double const next_knot = s->knot_angles[next];
	double const a0 = s->angles[prev];
	double const a1 = s->angles[j];
	double const h0 = s->lengths[prev];
	double const h1 = s->lengths[j];
	g2_range const r = g2_rotation_range(prev_knot, next_knot, a0, a1);
	double lo = r.lo, hi = r.hi;
	double lo_value = g2_scaled_curvature_delta(prev_knot, lo, next_knot, a0, a1, h0, h1);
	double hi_value = g2_scaled_curvature_delta(prev_knot, hi, next_knot, a0, a1, h0, h1);
	double sign = 1.0;
	double range, epsilon, truncation;
	unsigned int const steps = G2_ITP_HALVINGS + 1U;
	unsigned int step;

	if (lo_value == 0.0)
		return lo;
	if (hi_value == 0.0)
		return hi;
	if ((lo_value < 0.0) == (hi_value < 0.0))
		return hi;

	/* the steps below expect a negative value at lo and a positive one at hi */
	if (lo_value > 0.0)
		sign = -1.0;
	lo_value *= sign;
	hi_value *= sign;

	range = hi - lo;
	epsilon = ldexp(range, -(int)(G2_ITP_HALVINGS + 1U));
	truncation = G2_ITP_TRUNCATION / range;

	for (step = 0; step < steps; step++)
	{
		double const width = hi - lo;
		double middle, secant, offset, direction, truncated, radius, probe, value;
		if (width <= 2.0 * epsilon)
			break;
		middle = 0.5 * (lo + hi);

		/* regula falsi: inside the bracket because the values at its ends have opposite signs */
		secant = (hi_value * lo - lo_value * hi) / (hi_value - lo_value);

		/* truncation: towards the middle by an offset shrinking with the square of the width, never past the middle */
		offset = truncation * width * width;
		direction = secant > middle ? -1.0 : 1.0;
		truncated = middle;
		if (offset <= fabs(middle - secant))
			truncated = secant + direction * offset;

		/* projection: into a window around the middle whose radius halves with each step */
		radius = ldexp(epsilon, (int)(steps - step)) - 0.5 * width;
		if (radius < 0.0)
			radius = 0.0;
		probe = middle - direction * radius;
		if (fabs(truncated - middle) <= radius)
			probe = truncated;

		value = sign * g2_scaled_curvature_delta(prev_knot, probe, next_knot, a0, a1, h0, h1);
		if (value > 0.0)
		{
			hi = probe;
			hi_value = value;
		}
		else if (value < 0.0)
		{
			lo = probe;
			lo_value = value;
		}
		else
			return probe;
	}
	return 0.5 * (lo + hi);
}

/* The end rules of an open spline: the first and the last knot follow their neighbours. */
static double g2_open_start_angle(double first_segment_angle, double second_knot_angle)
{
	double const beta = g2_absolute_to_beta(second_knot_angle, first_segment_angle);
	double const alpha = 0.5 * sin(2.0 * beta);
	return g2_alpha_to_absolute(alpha, first_segment_angle);
}

static double g2_open_end_angle(double final_segment_angle, double knot_angle)
{
	double const alpha = g2_absolute_to_alpha(knot_angle, final_segment_angle);
	double const beta = 0.5 * sin(2.0 * alpha);
	return g2_beta_to_absolute(beta, final_segment_angle);
}

/* superdiagonal and right_side are count doubles each; x may be d */
static void g2_solve_tridiagonal(unsigned int count, double const* a, double const* b, double const* c, double const* d, double* x, double* scratch)
{
	double* superdiagonal = scratch;
	double* right_side = scratch + count;
	double previous_superdiagonal = c[0] / b[0];
	double previous_right_side = d[0] / b[0];
	double previous_solution;
	unsigned int index;

	superdiagonal[0] = previous_superdiagonal;
	right_side[0] = previous_right_side;
	for (index = 1; index < count; index++)
	{
		double const denominator = b[index] - a[index] * previous_superdiagonal;
		previous_superdiagonal = c[index] / denominator;
		previous_right_side = (d[index] - a[index] * previous_right_side) / denominator;
		superdiagonal[index] = previous_superdiagonal;
		right_side[index] = previous_right_side;
	}

	previous_solution = previous_right_side;
	x[count - 1] = previous_solution;
	for (index = count - 1; index > 0; index--)
	{
		unsigned int const previous = index - 1;
		previous_solution = right_side[previous] - superdiagonal[previous] * previous_solution;
		x[previous] = previous_solution;
	}
}

/* scratch: 4 count doubles; x may be d */
static void g2_solve_tridiagonal_cyclic(unsigned int count, double const* a, double const* b, double const* c, double const* d, double* x, double* scratch)
{
	unsigned int const last = count - 1;
	double* update = scratch + 2 * (size_t)count;
	double* diagonal = scratch + 3 * (size_t)count;
	double first_subdiagonal, last_superdiagonal, factor;
	unsigned int index;

	/* with two equations each cyclic link joins the same two unknowns as the diagonal next to it */
	if (count == 2)
	{
		double subdiagonal[2], superdiagonal[2];
		subdiagonal[0] = 0.0;
		subdiagonal[1] = a[1] + c[1];
		superdiagonal[0] = a[0] + c[0];
		superdiagonal[1] = 0.0;
		g2_solve_tridiagonal(2, subdiagonal, b, superdiagonal, d, x, scratch);
		return;
	}

	first_subdiagonal = a[0];
	last_superdiagonal = c[last];
	if (fabs(first_subdiagonal) < DBL_EPSILON && fabs(last_superdiagonal) < DBL_EPSILON)
	{
		g2_solve_tridiagonal(count, a, b, c, d, x, scratch);
		return;
	}

	/* Sherman-Morrison: the cyclic links are a rank one update of a tridiagonal matrix with a modified diagonal */
	memcpy(diagonal, b, sizeof(double) * count);
	diagonal[0] -= first_subdiagonal;
	diagonal[last] -= last_superdiagonal;
	memset(update, 0, sizeof(double) * count);
	update[0] = first_subdiagonal;
	update[last] = last_superdiagonal;

	g2_solve_tridiagonal(count, a, diagonal, c, update, update, scratch);
	g2_solve_tridiagonal(count, a, diagonal, c, d, x, scratch);

	factor = (x[0] + x[last]) / (1.0 + update[0] + update[last]);
	for (index = 0; index < count; index++)
		x[index] -= factor * update[index];
}

static void g2_update_curvatures(g2_solver* s)
{
	unsigned int j;
	for (j = 0; j < s->segments; j++)
	{
		double const segment_angle = s->angles[j];
		double const alpha = g2_absolute_to_alpha(s->knot_angles[j], segment_angle);
		double const beta = g2_absolute_to_beta(s->knot_angles[j + 1 < s->n ? j + 1 : 0], segment_angle);
		s->curvatures[j] = g2_segment_curvatures(alpha, beta, s->lengths[j]);
	}
}

/* Row i of the Newton system of the knot between the two segments: alpha moves with the knot at a segment's start and
   beta against the knot at its end, which sets the signs. The residual is the curvature difference. */
static void g2_knot_row(g2_solver* s, unsigned int i, g2_curvatures const* previous, g2_curvatures const* next)
{
	s->a[i] = -previous->end_by_alpha;
	s->b[i] = next->start_by_alpha + previous->end_by_beta;
	s->c[i] = -next->start_by_beta;
	s->x[i] = next->start - previous->end;
}

static void g2_solve_open(g2_solver* s)
{
	unsigned int const count = s->n;
	unsigned int const last = count - 1;           /* the last knot, also the number of segments */
	unsigned int const interior = count - 2;
	double* k = s->knot_angles;
	double const* angles = s->angles;
	double tolerance = interior == 0 ? G2_NEWTON_ERROR : G2_SWEEP_ERROR;
	unsigned int sweeps = 0, i, j;
	int converged = 0;

	for (i = 0; i < G2_MAX_CENTERING; i++)
	{
		k[0] = g2_open_start_angle(angles[0], k[1]);
		for (j = 1; j < last; j++)
			k[j] = g2_center_knot(s, j - 1, j, j + 1);
		k[last] = g2_open_end_angle(angles[last - 1], k[last - 1]);
	}

	while (!converged)
	{
		double sweep_error, previous_step;
		int diverged = 0;

		do
		{
			double const angle0 = g2_open_start_angle(angles[0], k[1]);
			double max_error = fabs(k[0] - angle0);
			double angle1;
			k[0] = angle0;
			for (j = 1; j < last; j++)
			{
				double const angle = g2_solve_knot_angle(s, j - 1, j, j + 1);
				max_error = g2_max(max_error, fabs(k[j] - angle));
				k[j] = angle;
			}
			angle1 = g2_open_end_angle(angles[last - 1], k[last - 1]);
			max_error = g2_max(max_error, fabs(k[last] - angle1));
			k[last] = angle1;
			sweeps++;
			sweep_error = max_error;
		} while (sweep_error >= tolerance && sweeps < G2_MAX_SWEEPS);

		if (interior == 0)
			break;

		/* a sweep sets the first knot before it solves the second, so the start rule is one sweep behind: apply it
		   again, for a Newton phase starting with both end rules satisfied */
		k[0] = g2_open_start_angle(angles[0], k[1]);
		memcpy(s->saved, k, sizeof(double) * count);

		/* the end knots follow their neighbours, so dy0 / dy1 = -cos(2 beta0) and dyn / dy(n - 1) = -cos(2 alpha(n - 1))
		   fold into the first and last rows; a growing step or a knot outside its range means the phase diverges */
		previous_step = G2_PI;
		for (i = 0; i < G2_MAX_NEWTON; i++)
		{
			double max_step = 0.0;
			double start_beta, end_alpha;

			g2_update_curvatures(s);
			for (j = 1; j < last; j++)
				g2_knot_row(s, j - 1, &s->curvatures[j - 1], &s->curvatures[j]);

			start_beta = g2_absolute_to_beta(k[1], angles[0]);
			end_alpha = g2_absolute_to_alpha(k[last - 1], angles[last - 1]);
			s->b[0] -= s->a[0] * cos(2.0 * start_beta);
			s->b[interior - 1] -= s->c[interior - 1] * cos(2.0 * end_alpha);

			if (interior == 1)
				s->x[0] = s->x[0] / s->b[0];
			else
				g2_solve_tridiagonal(interior, s->a, s->b, s->c, s->x, s->x, s->scratch);

			for (j = 1; j < last; j++)
				max_step = g2_max(max_step, fabs(s->x[j - 1]));
			if (max_step > previous_step)
			{
				diverged = 1;
				break;
			}
			previous_step = max_step;

			for (j = 1; j < last; j++)
				k[j] -= s->x[j - 1];
			k[0] = g2_open_start_angle(angles[0], k[1]);
			k[last] = g2_open_end_angle(angles[last - 1], k[last - 1]);

			for (j = 1; j < last; j++)
			{
				g2_range const r = g2_knot_range(s, j - 1, j, j + 1);
				if (k[j] < r.lo || k[j] > r.hi)
					diverged = 1;
			}
			if (diverged)
				break;
			if (max_step <= G2_NEWTON_ERROR)
			{
				converged = 1;
				break;
			}
		}
		if (converged)
			break;

		memcpy(k, s->saved, sizeof(double) * count);
		if (sweeps >= G2_MAX_SWEEPS)
			break;
		tolerance *= G2_FALLBACK_FACTOR;
	}
}

static void g2_solve_closed(g2_solver* s)
{
	unsigned int const count = s->n;
	unsigned int const last = count - 1;
	double* k = s->knot_angles;
	double tolerance = G2_SWEEP_ERROR;
	unsigned int sweeps = 0, i, j;
	int converged = 0;

	for (i = 0; i < G2_MAX_CENTERING; i++)
		for (j = 0; j < count; j++)
			k[j] = g2_center_knot(s, j ? j - 1 : last, j, j < last ? j + 1 : 0);

	while (!converged)
	{
		double sweep_error, previous_step;
		int diverged = 0;

		do
		{
			double max_error = 0.0;
			for (j = 0; j < count; j++)
			{
				double const angle = g2_solve_knot_angle(s, j ? j - 1 : last, j, j < last ? j + 1 : 0);
				double const error = fabs(k[j] - angle);
				max_error = j ? g2_max(max_error, error) : error;
				k[j] = angle;
			}
			sweeps++;
			sweep_error = max_error;
		} while (sweep_error >= tolerance && sweeps < G2_MAX_SWEEPS);

		memcpy(s->saved, k, sizeof(double) * count);

		previous_step = G2_PI;
		for (i = 0; i < G2_MAX_NEWTON; i++)
		{
			double max_step = 0.0;

			g2_update_curvatures(s);
			for (j = 0; j < count; j++)
				g2_knot_row(s, j, &s->curvatures[j ? j - 1 : last], &s->curvatures[j]);

			g2_solve_tridiagonal_cyclic(count, s->a, s->b, s->c, s->x, s->x, s->scratch);

			for (j = 0; j < count; j++)
				max_step = g2_max(max_step, fabs(s->x[j]));
			if (max_step > previous_step)
			{
				diverged = 1;
				break;
			}
			previous_step = max_step;

			for (j = 0; j < count; j++)
				k[j] -= s->x[j];

			for (j = 0; j < count; j++)
			{
				g2_range const r = g2_knot_range(s, j ? j - 1 : last, j, j < last ? j + 1 : 0);
				if (k[j] < r.lo || k[j] > r.hi)
					diverged = 1;
			}
			if (diverged)
				break;
			if (max_step <= G2_NEWTON_ERROR)
			{
				converged = 1;
				break;
			}
		}
		if (converged)
			break;

		memcpy(k, s->saved, sizeof(double) * count);
		if (sweeps >= G2_MAX_SWEEPS)
			break;
		tolerance *= G2_FALLBACK_FACTOR;
	}
}

/* The two control points and the end of one segment from p0 to p1, as 6 doubles. */
static void g2_write_cubic(double* out, double const* p0, double const* p1, double knot_angle0, double knot_angle1, double segment_angle, double length)
{
	double const a0 = g2_absolute_to_alpha(knot_angle0, segment_angle);
	double const a1 = g2_absolute_to_beta(knot_angle1, segment_angle);
	double const arm0 = g2_arm_length(a0, a1);
	double const arm1 = g2_arm_length(a1, a0);

	/* the control points of the unit cubic, from (0, 0) to (1, 0), scaled by the chord */
	double const c0x = cos(a0) * arm0 * length;
	double const c0y = -sin(a0) * arm0 * length;
	double const c1x = (1.0 - cos(a1) * arm1) * length;
	double const c1y = -sin(a1) * arm1 * length;

	double const s = sin(segment_angle);
	double const c = cos(segment_angle);

	out[0] = p0[0] + (c0x * c + c0y * s);
	out[1] = p0[1] + (-c0x * s + c0y * c);
	out[2] = p0[0] + (c1x * c + c1y * s);
	out[3] = p0[1] + (-c1x * s + c1y * c);
	out[4] = p1[0];
	out[5] = p1[1];
}

qaws_status qaws_internal_interpolate_cubic_g2(double** out_points, unsigned int* out_segments, double const* points, unsigned int n, int closed)
{
	unsigned int const segments = closed ? n : n - 1;
	/* per knot: curvatures (6), lengths, angles, knot angles, saved, a, b, c, x, 4 of solver scratch */
	double* work = (double*)qaws_internal_malloc(qaws_internal_mul3(18, n, sizeof(double)));
	double* out = (double*)qaws_internal_malloc(qaws_internal_mul3(2, qaws_internal_add2(qaws_internal_mul2(3, segments), 1), sizeof(double)));
	g2_solver s;
	size_t coordinate, coordinates;
	unsigned int i;

	*out_points = NULL;
	*out_segments = 0;

	if (!work || !out)
	{
		qaws_internal_free(work);
		qaws_internal_free(out);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	s.curvatures = (g2_curvatures*)work;
	s.lengths = work + 6 * (size_t)n;
	s.angles = s.lengths + n;
	s.knot_angles = s.angles + n;
	s.saved = s.knot_angles + n;
	s.a = s.saved + n;
	s.b = s.a + n;
	s.c = s.b + n;
	s.x = s.c + n;
	s.scratch = s.x + n;
	s.n = n;
	s.segments = segments;

	for (i = 0; i < segments; i++)
	{
		unsigned int const next = i + 1 < n ? i + 1 : 0;
		double const dx = points[2 * (size_t)next] - points[2 * (size_t)i];
		double const dy = points[2 * (size_t)next + 1] - points[2 * (size_t)i + 1];
		s.lengths[i] = sqrt(dx * dx + dy * dy);
		s.angles[i] = g2_vector_angle(dx, dy);
	}

	for (i = 1; i + 1 < n; i++)
	{
		double const a0 = s.angles[i - 1];
		double const bend = s.angles[i] - a0;
		double const l0 = s.lengths[i - 1];
		s.knot_angles[i] = a0 + ((bend * l0) / (l0 + s.lengths[i]));
	}
	s.knot_angles[0] = g2_alpha_to_absolute(0.0, s.angles[0]);
	s.knot_angles[n - 1] = g2_beta_to_absolute(0.0, s.angles[segments - 1]);

	if (closed)
		g2_solve_closed(&s);
	else
		g2_solve_open(&s);

	out[0] = points[0];
	out[1] = points[1];
	for (i = 0; i < segments; i++)
	{
		unsigned int const next = i + 1 < n ? i + 1 : 0;
		g2_write_cubic(out + 2 + 6 * (size_t)i, points + 2 * (size_t)i, points + 2 * (size_t)next, s.knot_angles[i], s.knot_angles[next], s.angles[i], s.lengths[i]);
	}
	qaws_internal_free(work);

	coordinates = 2 * (3 * (size_t)segments + 1);
	for (coordinate = 0; coordinate < coordinates; coordinate++)
		if (!isfinite(out[coordinate]))
		{
			qaws_internal_free(out);
			return QAWS_STATUS_NUMERICAL_FAILURE;
		}

	*out_points = out;
	*out_segments = segments;
	return QAWS_STATUS_OK;
}
