/*
 * Curvature-matching spline through points, built from S9 curves and returned as a chain of cubic Bezier segments.
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

#include "internal/qaws_internal_bezier_chain.h"
#include "internal/qaws_internal_memory.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/*
 * One S9 curve joins a unit chord from (0, 0) to (1, 0), y up, and a positive curvature turns counter-clockwise. With t
 * the normalized arc length, its turn rate is
 *
 *     f(t) = StartLine g0(t) + EndLine g1(t)
 *     g0(t) = (1 - t) + c0(t)
 *     g1(t) = t + c1(1 - t)
 *
 * and the corner term at the distance r from its end, with x = lambda r, is
 *
 *     c(r) = Q (p(x) + phi q(x))
 *     p(x) = 1 / (1 + x)^1.5
 *     q(x) = eta x / (1 + x)^3
 *     phi = lambda / (lambda + lambda_c)
 *
 * The scale Q makes the integral over the curve equal to lambda. Both terms and their integrals are polynomials in
 * u = 1 / sqrt(1 + x), so the tangent angle needs no quadrature.
 */

#define S9_PI 3.141592653589793238462643383279502884

#define S9_ROUNDING_LAMBDA 0.5             /* the corner lambda at which the rounding term has half of its weight */
#define S9_ROUNDING_ETA 1.5

/* a smaller weight magnitude is raised to this: the corner is then numerically infinite, its pole about 1e-8 away */
#define S9_SMALLEST_WEIGHT 1e-8

/* quadrature grid */
#define S9_BASE_PANEL_COUNT 8U
#define S9_COARSE_BASE_PANEL_COUNT 4U
#define S9_GREATEST_GRADED_EDGE_COUNT 30U
#define S9_GREATEST_EDGE_COUNT (S9_BASE_PANEL_COUNT + 1U + 2U * S9_GREATEST_GRADED_EDGE_COUNT)
#define S9_GREATEST_NODE_COUNT ((S9_GREATEST_EDGE_COUNT - 1U) * 16U)
#define S9_GREATEST_COARSE_NODE_COUNT ((S9_COARSE_BASE_PANEL_COUNT + 2U * S9_GREATEST_GRADED_EDGE_COUNT) * 8U)

/* root search */
#define S9_SCAN_SAMPLE_COUNT 1440U
#define S9_SCAN_REACH (6.0 * S9_PI)
#define S9_ROOT_TOLERANCE 1e-12
#define S9_REFINE_COUNT 64U
#define S9_SAME_TURN_TOLERANCE 1e-9
#define S9_SMALLEST_CHORD_LENGTH 1e-9
#define S9_SCAN_RESTART_INTERVAL 32U
#define S9_PREDICTED_ANGLE_LIMIT 2.1
#define S9_PREDICTED_REACH 3.0
#define S9_PREDICTED_STEP_COUNT 16U
#define S9_GREATEST_NEAR_STEP_COUNT 12U
#define S9_NEAR_REACH 2.0
#define S9_GREATEST_NEAR_ANGLE_CHANGE 0.5

/* cubic approximation of one curve */
#define S9_GREATEST_PANEL_TURN 0.25
#define S9_INTEGRATION_ALLOWANCE 1e-12
#define S9_GREATEST_DEPTH 34U
#define S9_GREATEST_CUBIC_COUNT 2048U
#define S9_GREATEST_PANEL_COUNT 60000U
#define S9_MEASUREMENTS_PER_RANGE 5U
#define S9_PIECE_POINT_COUNT (3U * S9_GREATEST_CUBIC_COUNT + 1U)

/* distance between the curve of a segment and its cubics, as a fraction of the chord (checked at three points per cubic) */
#define S9_MAXIMUM_CUBIC_ERROR 1e-5

/* distance from an end of the rotation range at which the solver reads the sign of the curvature delta: at the end
   itself the weight is zero, and rounding can make it negative, which reverses that end and can change the sign */
#define S9_SIGN_PROBE_STEP 1e-4

/* the root of one knot: the tolerance is the noise floor of the knot angles (with a smaller one angles still move by up
   to 3e-13 radians from sweep to sweep through the rounding of the curve solves) */
#define S9_KNOT_ANGLE_TOLERANCE 1e-13
#define S9_GREATEST_ROOT_ITERATION_COUNT 64U
#define S9_FIRST_BRACKET_STEP 1e-3
#define S9_SMALLEST_BRACKET_STEP 1e-9
#define S9_BRACKET_STEP_GROWTH 4.0

#define S9_MAXIMUM_CENTERING_ITERATION_COUNT 3U
#define S9_MAXIMUM_ITP_ITERATION_COUNT 20U
#define S9_MAXIMUM_ITP_ERROR 1e-11         /* sweeps stop when no knot angle moves more than this */
#define S9_OSCILLATION_DAMPING 0.7         /* the part of an update kept when it goes against the knot's previous update */

/* the coupled Newton iteration */
#define S9_KNOT_SWEEP_COUNT 12U
#define S9_COUPLED_START_ERROR 3e-4
#define S9_COUPLED_GREATEST_ITERATION_COUNT 40U
#define S9_COUPLED_GREATEST_HALVING_COUNT 5U
#define S9_COUPLED_DIFFERENCE_STEP 1e-6
#define S9_COUPLED_DECREASE 1e-4
#define S9_COUPLED_RANGE_FRACTION 0.25
#define S9_COUPLED_STALL_FACTOR 0.2
#define S9_COUPLED_GREATEST_STALL_COUNT 3U

/* a knot whose rotation range is narrower than the first is pinned to the range's middle (both its ends are within 0.3
   degrees of their curvature poles), and freed only when wider than the second, so a range moving about the pinning
   width does not change the mark at every iteration */
#define S9_PINNED_RANGE_WIDTH 1e-2
#define S9_UNPIN_RANGE_WIDTH 2e-2
#define S9_CENTRING_PASS_COUNT 3U

static double const s9_gauss_legendre_x_8[8] = {
	0.0198550717512319119251,
	0.1016667612931866360171,
	0.2372337950418355045912,
	0.4082826787521750544485,
	0.5917173212478249455515,
	0.7627662049581644954088,
	0.8983332387068133639829,
	0.9801449282487680880749,
};

static double const s9_gauss_legendre_scale_8[8] = {
	0.0506142681451880876997,
	0.1111905172266874491926,
	0.1568533229389436067791,
	0.1813418916891808563285,
	0.1813418916891808563285,
	0.1568533229389436067791,
	0.1111905172266874491926,
	0.0506142681451880876997,
};

static double const s9_gauss_legendre_x_16[16] = {
	0.0052995325041750307449,
	0.0277124884633836998660,
	0.0671843988060840668908,
	0.1222977958224985006730,
	0.1910618777986781147149,
	0.2709916111713863151600,
	0.3591982246103705422868,
	0.4524937450811812866824,
	0.5475062549188187688287,
	0.6408017753896294577132,
	0.7290083888286136293289,
	0.8089381222013218852851,
	0.8777022041775015548382,
	0.9328156011939159331092,
	0.9722875115366163001340,
	0.9947004674958249692551,
};

static double const s9_gauss_legendre_scale_16[16] = {
	0.0135762297058778080155,
	0.0311267619693236866074,
	0.0475792558412461707995,
	0.0623144856277668551803,
	0.0747979944082883402290,
	0.0845782596975012124219,
	0.0913017075224617363771,
	0.0947253052275341816957,
	0.0947253052275341816957,
	0.0913017075224617363771,
	0.0845782596975012124219,
	0.0747979944082883402290,
	0.0623144856277668551803,
	0.0475792558412461707995,
	0.0311267619693236866074,
	0.0135762297058778080155,
};

typedef struct s9_point
{
	double x, y;
} s9_point;

/* the inputs of one curve: angles in radians (alpha counter-clockwise from the chord, beta clockwise), signed weights
   from -1 to 1 whose magnitude sets the corner lambda and whose negative value reverses that end, and the etas */
typedef struct s9_curve_parameters
{
	double alpha, beta;
	double start_weight, end_weight;
	double start_eta, end_eta;
} s9_curve_parameters;

/* one curve for an implicit unit chord */
typedef struct s9_curve
{
	double start_angle, end_angle;     /* forward tangent angles, counter-clockwise from the chord */
	double start_lambda, end_lambda;
	double start_eta, end_eta;
	double turn;                       /* net turn: the difference of the angles, possibly plus or minus a full turn */
	double start_line, end_line;       /* turn rate of the line part at each end */
	double length;                     /* of the curve, for a chord of length one */
} s9_curve;

/* the quadrature nodes of one curve */
typedef struct s9_grid
{
	unsigned int count;
	double weight[S9_GREATEST_NODE_COUNT];
	double start_integral[S9_GREATEST_NODE_COUNT];   /* G0 */
	double end_integral[S9_GREATEST_NODE_COUNT];     /* G1 */
	double share_rate[S9_GREATEST_NODE_COUNT];       /* how the tangent angle changes with the start share */
	double start_rate[S9_GREATEST_NODE_COUNT];       /* g0, built on request */
	double end_rate[S9_GREATEST_NODE_COUNT];         /* g1, built on request */
} s9_grid;

/* the end point of a curve in units of its length, and its derivative with respect to the start share */
typedef struct s9_end_point
{
	s9_point point, slope;
} s9_end_point;

/* the best root of a scan so far: the one with the least total turning */
typedef struct s9_root_choice
{
	double share, distance, turning;
	s9_point end_point;
	int found;
} s9_root_choice;

/* the cubics of one curve (S9_PIECE_POINT_COUNT points are enough for S9_GREATEST_CUBIC_COUNT cubics) */
typedef struct s9_pieces
{
	s9_point* points;
	unsigned int count;
} s9_pieces;

typedef struct s9_approximation
{
	s9_curve const* shape;
	s9_pieces* pieces;
	double maximum_error;
	unsigned int remaining_measurements, remaining_cubics;
	int failed;
} s9_approximation;

/* the result: a chain of cubics */
typedef struct s9_output
{
	double* points;                    /* x, y */
	size_t point_count, capacity;
} s9_output;

/* the last solved curve of a segment: the next solve of the segment continues it, which is much faster */
typedef struct s9_segment_curve
{
	s9_curve shape;
	int is_known;
} s9_segment_curve;

/* the work arrays, each of as many doubles as knots */
enum
{
	S9_WORK_SEGMENT_LENGTHS,
	S9_WORK_SEGMENT_ANGLES,
	S9_WORK_KNOT_ANGLES,
	S9_WORK_PREVIOUS_UPDATES,
	S9_WORK_TRIAL,
	S9_WORK_START_CURVATURES,
	S9_WORK_END_CURVATURES,
	S9_WORK_TRIAL_START_CURVATURES,
	S9_WORK_TRIAL_END_CURVATURES,
	S9_WORK_DELTAS,
	S9_WORK_TRIAL_DELTAS,
	S9_WORK_SUBDIAGONAL,
	S9_WORK_DIAGONAL,
	S9_WORK_SUPERDIAGONAL,
	S9_WORK_RIGHT_SIDE,
	S9_WORK_STEP,
	S9_WORK_LOWS,
	S9_WORK_HIGHS,
	S9_WORK_TRIAL_ANGLES,
	S9_WORK_SOLVE_SUPERDIAGONAL,
	S9_WORK_SOLVE_RIGHT_SIDE,
	S9_WORK_SOLVE_UPDATE,
	S9_WORK_SOLVE_DIAGONAL,
	S9_WORK_COUNT
};

/* all knots and segments of one spline */
typedef struct s9_knot_system
{
	unsigned int count;                /* of knots */
	unsigned int segment_count;        /* the knot count when closed (the last segment ends at knot zero), else one less */
	int is_cyclic;
	double const* segment_angles;      /* absolute, radians */
	double const* segment_lengths;     /* chords */
	double* knot_angles;               /* absolute, radians; the solver updates them */
	s9_segment_curve* segment_curves;
	signed char* pinned;               /* per knot, see s9_classify_knots */
	double* work[S9_WORK_COUNT];
} s9_knot_system;

/* written as comparisons rather than fabs and fmax, whose treatment of NaN and of the sign of zero would differ */
static double s9_abs(double t) { return t >= 0 ? t : -t; }
static double s9_max(double a, double b) { return a < b ? b : a; }
static double s9_min(double a, double b) { return a < b ? a : b; }
static double s9_clamp(double value, double low, double high) { return value > high ? high : value < low ? low : value; }

static s9_point s9_xy(double x, double y)
{
	s9_point p;
	p.x = x;
	p.y = y;
	return p;
}

static s9_point s9_sub(s9_point a, s9_point b) { return s9_xy(a.x - b.x, a.y - b.y); }

/* the angle of a vector; positive is counter-clockwise */
static double s9_vector_angle(s9_point v) { return atan2(-v.y, v.x); }

static s9_curve s9_curve_default(void)
{
	s9_curve c;
	c.start_angle = 0.0;
	c.end_angle = 0.0;
	c.start_lambda = 0.0;
	c.end_lambda = 0.0;
	c.start_eta = S9_ROUNDING_ETA;
	c.end_eta = S9_ROUNDING_ETA;
	c.turn = 0.0;
	c.start_line = 0.0;
	c.end_line = 0.0;
	c.length = 1.0;
	return c;
}

static double s9_corner_rounding(double lambda)
{
	return lambda / (lambda + S9_ROUNDING_LAMBDA);
}

/* The mean of the unit corner integrand, the corner term without its scale Q, over an interval of scaled distances:
       u1 = 1 / sqrt(1 + x1), u2 = 1 / sqrt(1 + x2)
       mean = u1^2 u2^2 (2 / (u1 + u2) + eta phi / 2 (x1 u1^2 + x2 u2^2)) */
static double s9_unit_corner_mean(double near_distance, double far_distance, double rounding, double eta)
{
	double const near_root = 1.0 / sqrt(1.0 + near_distance);
	double const far_root = 1.0 / sqrt(1.0 + far_distance);
	double const product = near_root * far_root;
	double const near_square = near_root * near_root;
	double const far_square = far_root * far_root;
	double const square_product = product * product;
	double const pole_mean = 2.0 / (near_root + far_root);
	double const rounding_mean = 0.5 * eta * rounding * (near_distance * near_square + far_distance * far_square);

	return square_product * (pole_mean + rounding_mean);
}

/* the scale Q of the unit corner term; a corner lambda of zero gives zero, an end without a corner */
static double s9_unit_corner_scale(double lambda, double eta)
{
	double const rounding = s9_corner_rounding(lambda);
	double const mean = s9_unit_corner_mean(0.0, lambda, rounding, eta);

	return lambda / mean;
}

/* the unit corner term at a distance from its end */
static double s9_unit_corner_value_with_scale(double lambda, double eta, double scale, double distance)
{
	double const scaled = distance * lambda;
	double const rounding = s9_corner_rounding(lambda);
	double const unit_root = 1.0 / sqrt(1.0 + scaled);
	double const pole = unit_root * unit_root * unit_root;
	double const rounding_term = eta * rounding * scaled * pole * pole;

	return scale * (pole + rounding_term);
}

static double s9_unit_corner_value(double lambda, double eta, double distance)
{
	double const scale = s9_unit_corner_scale(lambda, eta);

	return s9_unit_corner_value_with_scale(lambda, eta, scale, distance);
}

static double s9_unit_corner_turn_from_end_with_scale(double lambda, double eta, double scale, double distance)
{
	double const scaled_distance = distance * lambda;
	double const rounding = s9_corner_rounding(lambda);
	double const mean = s9_unit_corner_mean(0.0, scaled_distance, rounding, eta);

	return scale * distance * mean;
}

static double s9_unit_corner_turn_from_end(double lambda, double eta, double distance)
{
	double const scale = s9_unit_corner_scale(lambda, eta);

	return s9_unit_corner_turn_from_end_with_scale(lambda, eta, scale, distance);
}

static double s9_unit_corner_turn_to_far_end_with_scale(double lambda, double eta, double scale, double width)
{
	double const near_distance = (1.0 - width) * lambda;
	double const rounding = s9_corner_rounding(lambda);
	double const mean = s9_unit_corner_mean(near_distance, lambda, rounding, eta);

	return scale * width * mean;
}

static double s9_unit_corner_turn_to_far_end(double lambda, double eta, double width)
{
	double const scale = s9_unit_corner_scale(lambda, eta);

	return s9_unit_corner_turn_to_far_end_with_scale(lambda, eta, scale, width);
}

/* the turn J = 1/2 + lambda that a basis makes over the whole curve */
static double s9_start_basis_turn(s9_curve const* shape) { return 0.5 + shape->start_lambda; }
static double s9_end_basis_turn(s9_curve const* shape) { return 0.5 + shape->end_lambda; }

static double s9_start_basis(s9_curve const* shape, double t)
{
	double const corner = s9_unit_corner_value(shape->start_lambda, shape->start_eta, t);

	return (1.0 - t) + corner;
}

static double s9_end_basis(s9_curve const* shape, double t)
{
	double const corner = s9_unit_corner_value(shape->end_lambda, shape->end_eta, 1.0 - t);

	return t + corner;
}

static double s9_start_basis_integral(s9_curve const* shape, double t)
{
	double const corner = s9_unit_corner_turn_from_end(shape->start_lambda, shape->start_eta, t);

	return t - 0.5 * t * t + corner;
}

static double s9_end_basis_integral(s9_curve const* shape, double t)
{
	double const corner = s9_unit_corner_turn_to_far_end(shape->end_lambda, shape->end_eta, t);

	return 0.5 * t * t + corner;
}

static double s9_turn_rate(s9_curve const* shape, double t)
{
	double const start_basis = s9_start_basis(shape, t);
	double const end_basis = s9_end_basis(shape, t);

	return shape->start_line * start_basis + shape->end_line * end_basis;
}

/* the signed curvature, for a chord of length one */
static double s9_curvature(s9_curve const* shape, double t)
{
	double const turn_rate = s9_turn_rate(shape, t);

	return turn_rate / shape->length;
}

static double s9_turn_from_start(s9_curve const* shape, double t)
{
	double const start_integral = s9_start_basis_integral(shape, t);
	double const end_integral = s9_end_basis_integral(shape, t);

	return shape->start_line * start_integral + shape->end_line * end_integral;
}

static double s9_tangent_angle(s9_curve const* shape, double t)
{
	double const turn = s9_turn_from_start(shape, t);

	return shape->start_angle + turn;
}

static double s9_reduce_angle(double angle)
{
	return remainder(angle, 2.0 * S9_PI);
}

static double s9_heading(double forward_angle, double weight)
{
	double heading = forward_angle;

	if (weight < 0.0) {
		heading += S9_PI;
	}

	return s9_reduce_angle(heading);
}

static double s9_weight_lambda(double weight)
{
	double const magnitude = s9_abs(weight);
	double const limited = s9_clamp(magnitude, S9_SMALLEST_WEIGHT, 1.0);

	return (1.0 - limited) / limited;
}

/* sets the angles, the corner lambdas, the etas and the net turn of a shape from its parameters */
static void s9_prepare_shape(s9_curve_parameters const* parameters, s9_curve* shape)
{
	shape->start_angle = s9_heading(parameters->alpha, parameters->start_weight);
	shape->end_angle = s9_heading(-parameters->beta, parameters->end_weight);
	shape->start_lambda = s9_weight_lambda(parameters->start_weight);
	shape->end_lambda = s9_weight_lambda(parameters->end_weight);
	shape->start_eta = parameters->start_eta;
	shape->end_eta = parameters->end_eta;
	shape->turn = s9_reduce_angle(shape->end_angle - shape->start_angle);
}

static double s9_raw_turn(s9_curve const* shape)
{
	return shape->end_angle - shape->start_angle;
}

static void s9_set_lines(double share, s9_curve* shape)
{
	double const start_basis_turn = s9_start_basis_turn(shape);
	double const end_basis_turn = s9_end_basis_turn(shape);

	shape->start_line = share / start_basis_turn;
	shape->end_line = (shape->turn - share) / end_basis_turn;
}

/* the start share of the net turn of a solved shape */
static double s9_start_share(s9_curve const* shape)
{
	double const start_basis_turn = s9_start_basis_turn(shape);

	return shape->start_line * start_basis_turn;
}

static void s9_add_graded_edges(double lambda, int from_end, double* edges, unsigned int* count)
{
	double distance;
	unsigned int added = 0U;

	if (lambda <= 0.0) {
		return;
	}

	distance = 0.125 / lambda;

	while (distance < 0.5 && added < S9_GREATEST_GRADED_EDGE_COUNT && *count < S9_GREATEST_EDGE_COUNT) {
		double edge = distance;

		if (from_end) {
			edge = 1.0 - distance;
		}

		edges[*count] = edge;
		(*count)++;
		added++;
		distance *= 2.0;
	}
}

static void s9_build_grid(s9_curve const* shape, unsigned int base_panel_count, unsigned int nodes_per_panel, int with_rates,
	s9_grid* grid)
{
	double edges[S9_GREATEST_EDGE_COUNT];
	unsigned int edge_count = 0U;
	unsigned int index, node, count, slot;
	double const* node_x = s9_gauss_legendre_x_16;
	double const* node_scale = s9_gauss_legendre_scale_16;
	double const start_lambda = shape->start_lambda;
	double const end_lambda = shape->end_lambda;
	double const start_eta = shape->start_eta;
	double const end_eta = shape->end_eta;
	double const start_scale = s9_unit_corner_scale(start_lambda, start_eta);
	double const end_scale = s9_unit_corner_scale(end_lambda, end_eta);
	double const start_basis_turn = s9_start_basis_turn(shape);
	double const end_basis_turn = s9_end_basis_turn(shape);

	for (index = 0U; index <= base_panel_count; index++) {
		edges[edge_count] = (double)index / (double)base_panel_count;
		edge_count++;
	}

	s9_add_graded_edges(start_lambda, 0, edges, &edge_count);
	s9_add_graded_edges(end_lambda, 1, edges, &edge_count);

	/* insertion sort: the list has at most S9_GREATEST_EDGE_COUNT edges */
	for (index = 1U; index < edge_count; index++) {
		double const edge = edges[index];
		unsigned int position = index;

		while (position > 0U && edges[position - 1U] > edge) {
			edges[position] = edges[position - 1U];
			position--;
		}

		edges[position] = edge;
	}

	if (nodes_per_panel == 8U) {
		node_x = s9_gauss_legendre_x_8;
		node_scale = s9_gauss_legendre_scale_8;
	}

	count = 0U;

	for (index = 1U; index < edge_count; index++) {
		double const begin = edges[index - 1U];
		double const end = edges[index];
		double const width = end - begin;

		/* a duplicate edge gives an empty panel, which adds nothing */
		if (width <= 0.0) {
			continue;
		}

		for (node = 0U; node < nodes_per_panel; node++) {
			double const t = begin + width * node_x[node];
			double const start_corner = s9_unit_corner_turn_from_end_with_scale(start_lambda, start_eta, start_scale, t);
			double const end_corner = s9_unit_corner_turn_to_far_end_with_scale(end_lambda, end_eta, end_scale, t);
			double const start_integral = t - 0.5 * t * t + start_corner;
			double const end_integral = 0.5 * t * t + end_corner;

			grid->weight[count] = width * node_scale[node];
			grid->start_integral[count] = start_integral;
			grid->end_integral[count] = end_integral;
			grid->share_rate[count] = start_integral / start_basis_turn - end_integral / end_basis_turn;
			count++;
		}
	}

	grid->count = count;

	if (!with_rates) {
		return;
	}

	/* a second pass for the basis curvatures, at the same nodes */
	slot = 0U;

	for (index = 1U; index < edge_count; index++) {
		double const begin = edges[index - 1U];
		double const end = edges[index];
		double const width = end - begin;

		if (width <= 0.0) {
			continue;
		}

		for (node = 0U; node < nodes_per_panel; node++) {
			double const t = begin + width * node_x[node];
			double const start_corner = s9_unit_corner_value_with_scale(start_lambda, start_eta, start_scale, t);
			double const end_corner = s9_unit_corner_value_with_scale(end_lambda, end_eta, end_scale, 1.0 - t);

			grid->start_rate[slot] = (1.0 - t) + start_corner;
			grid->end_rate[slot] = t + end_corner;
			slot++;
		}
	}
}

static s9_end_point s9_measure_end_point(s9_grid const* grid, s9_curve* shape, double share)
{
	s9_end_point result;
	double start_angle, start_line, end_line;
	double horizontal = 0.0;
	double vertical = 0.0;
	double horizontal_slope = 0.0;
	double vertical_slope = 0.0;
	unsigned int index;

	s9_set_lines(share, shape);

	start_angle = shape->start_angle;
	start_line = shape->start_line;
	end_line = shape->end_line;

	for (index = 0U; index < grid->count; index++) {
		double const start_integral = grid->start_integral[index];
		double const end_integral = grid->end_integral[index];
		double const angle = start_angle + start_line * start_integral + end_line * end_integral;
		double const rate = grid->share_rate[index];
		double const weight = grid->weight[index];
		double const cosine = cos(angle);
		double const sine = sin(angle);

		horizontal += weight * cosine;
		vertical += weight * sine;
		horizontal_slope -= weight * sine * rate;
		vertical_slope += weight * cosine * rate;
	}

	result.point = s9_xy(horizontal, vertical);
	result.slope = s9_xy(horizontal_slope, vertical_slope);

	return result;
}

static double s9_measure_turning(s9_grid const* grid, s9_curve* shape, double share)
{
	double start_line, end_line;
	double turning = 0.0;
	unsigned int index;

	s9_set_lines(share, shape);

	start_line = shape->start_line;
	end_line = shape->end_line;

	for (index = 0U; index < grid->count; index++) {
		double const turn_rate = start_line * grid->start_rate[index] + end_line * grid->end_rate[index];

		turning += grid->weight[index] * s9_abs(turn_rate);
	}

	return turning;
}

/* sets the length of a shape from the end point of its solved share */
static void s9_complete_shape(s9_point end_point, s9_curve* shape)
{
	double const distance = hypot(end_point.x, end_point.y);

	shape->length = 1.0 / distance;
}

static int s9_parts_bracket_root(double first, double second)
{
	int const first_is_negative = first < 0.0;
	int const second_is_negative = second < 0.0;

	return first_is_negative != second_is_negative;
}

static int s9_end_point_is_root(s9_point end_point)
{
	int const is_on_chord_line = isfinite(end_point.y) && s9_abs(end_point.y) <= S9_ROOT_TOLERANCE;

	return is_on_chord_line && end_point.x > 0.0;
}

static void s9_offer_root(s9_grid const* fine, s9_curve* shape, double share, s9_point end_point, s9_root_choice* choice)
{
	double turning;

	if (!(end_point.x > 0.0)) {
		return;
	}

	turning = s9_measure_turning(fine, shape, share);

	if (!choice->found || turning < choice->turning) {
		choice->share = share;
		choice->distance = hypot(end_point.x, end_point.y);
		choice->turning = turning;
		choice->end_point = end_point;
		choice->found = 1;
	}
}

static int s9_try_refine_root(s9_grid const* grid, s9_curve* shape, double low, double high, double low_part, double high_part,
	double* root, s9_point* end_point)
{
	double low_share = low;
	double high_share = high;
	double low_value = low_part;
	double high_value = high_part;
	double guess, value_span;
	unsigned int index;

	if (s9_abs(low_part) <= S9_ROOT_TOLERANCE) {
		s9_end_point const measured = s9_measure_end_point(grid, shape, low);

		*root = low;
		*end_point = measured.point;

		return 1;
	}

	if (s9_abs(high_part) <= S9_ROOT_TOLERANCE) {
		s9_end_point const measured = s9_measure_end_point(grid, shape, high);

		*root = high;
		*end_point = measured.point;

		return 1;
	}

	guess = low_share + 0.5 * (high_share - low_share);
	value_span = low_value - high_value;

	if (value_span != 0.0) {
		double const fraction = low_value / value_span;

		if (fraction > 0.02 && fraction < 0.98) {
			guess = low_share + (high_share - low_share) * fraction;
		}
	}

	for (index = 0U; index < S9_REFINE_COUNT; index++) {
		s9_end_point measured;
		double value, slope, next;
		int same_side_as_low;

		if (guess == low_share || guess == high_share) {
			break;
		}

		measured = s9_measure_end_point(grid, shape, guess);
		value = measured.point.y;

		if (!isfinite(value)) {
			break;
		}

		if (s9_abs(value) <= S9_ROOT_TOLERANCE) {
			*root = guess;
			*end_point = measured.point;

			return 1;
		}

		same_side_as_low = (value < 0.0) == (low_value < 0.0);

		if (same_side_as_low) {
			low_share = guess;
			low_value = value;
		} else {
			high_share = guess;
			high_value = value;
		}

		/* a Newton step from this measurement, kept inside the bracket */
		slope = measured.slope.y;
		next = low_share + 0.5 * (high_share - low_share);

		if (slope != 0.0 && isfinite(slope)) {
			double const newton = guess - value / slope;

			if (newton > low_share && newton < high_share) {
				next = newton;
			}
		}

		guess = next;
	}

	return 0;
}

static void s9_refine_bracket(s9_grid const* fine, s9_curve* shape, double low, double high, s9_point low_point, s9_point high_point,
	s9_root_choice* choice)
{
	double root = 0.0;
	s9_point end_point = s9_xy(0.0, 0.0);

	if (!(isfinite(low_point.y) && isfinite(high_point.y))) {
		return;
	}

	if (!s9_parts_bracket_root(low_point.y, high_point.y)) {
		return;
	}

	if (s9_try_refine_root(fine, shape, low, high, low_point.y, high_point.y, &root, &end_point)) {
		s9_offer_root(fine, shape, root, end_point, choice);
	}
}

static void s9_examine_bracket(s9_grid const* fine, s9_curve* shape, double low, double high, double step, s9_root_choice* choice)
{
	s9_end_point const fine_low = s9_measure_end_point(fine, shape, low);
	s9_end_point const fine_high = s9_measure_end_point(fine, shape, high);
	int const low_is_root = s9_end_point_is_root(fine_low.point);
	int const high_is_root = s9_end_point_is_root(fine_high.point);
	double lower_low, upper_high;
	s9_end_point fine_lower_low, fine_upper_high;
	int lower_brackets;

	if (low_is_root) {
		s9_offer_root(fine, shape, low, fine_low.point, choice);
	}

	if (high_is_root) {
		s9_offer_root(fine, shape, high, fine_high.point, choice);
	}

	if (low_is_root || high_is_root) {
		return;
	}

	if (!(isfinite(fine_low.point.y) && isfinite(fine_high.point.y))) {
		return;
	}

	if (s9_parts_bracket_root(fine_low.point.y, fine_high.point.y)) {
		s9_refine_bracket(fine, shape, low, high, fine_low.point, fine_high.point, choice);

		return;
	}

	/* the root is near one end, in the interval on the other side of that end, which the coarse grid did not see as a bracket */
	lower_low = low - step;
	fine_lower_low = s9_measure_end_point(fine, shape, lower_low);
	lower_brackets = isfinite(fine_lower_low.point.y) && s9_parts_bracket_root(fine_lower_low.point.y, fine_high.point.y);

	if (lower_brackets) {
		s9_refine_bracket(fine, shape, lower_low, high, fine_lower_low.point, fine_high.point, choice);

		return;
	}

	upper_high = high + step;
	fine_upper_high = s9_measure_end_point(fine, shape, upper_high);

	s9_refine_bracket(fine, shape, low, upper_high, fine_low.point, fine_upper_high.point, choice);
}

/* The end points of the scan's shares. The angle of a node is linear in the share, so a share step adds a fixed delta per
   node: the angle addition identities update each node's sine and cosine from those of its delta, instead of calling
   sin and cos for every share. */
static void s9_sample_share_scan(s9_grid const* grid, s9_curve* shape, s9_point* points)
{
	double const start_basis_turn = s9_start_basis_turn(shape);
	double const end_basis_turn = s9_end_basis_turn(shape);
	double const first_share = 0.5 * shape->turn - S9_SCAN_REACH;
	double const step = 2.0 * S9_SCAN_REACH / (double)(S9_SCAN_SAMPLE_COUNT - 1U);
	double cosines[S9_GREATEST_COARSE_NODE_COUNT];
	double sines[S9_GREATEST_COARSE_NODE_COUNT];
	double step_cosines[S9_GREATEST_COARSE_NODE_COUNT];
	double step_sines[S9_GREATEST_COARSE_NODE_COUNT];
	double weight_magnitude = 0.0;
	double greatest_integral_magnitude = 0.0;
	double last_share, absolute_first_share, absolute_last_share, greatest_share, absolute_start_angle, absolute_turn;
	double phase_scale, weight_scale, sign_allowance;
	unsigned int index, sample;

	for (index = 0U; index < grid->count; index++) {
		double const start_integral = grid->start_integral[index] / start_basis_turn;
		double const end_integral = grid->end_integral[index] / end_basis_turn;
		double const rate = start_integral - end_integral;
		double const increment = step * rate;
		double const absolute_start_integral = s9_abs(start_integral);
		double const absolute_end_integral = s9_abs(end_integral);
		double const integral_magnitude = absolute_start_integral + absolute_end_integral;
		double const absolute_weight = s9_abs(grid->weight[index]);

		step_cosines[index] = cos(increment);
		step_sines[index] = sin(increment);
		weight_magnitude += absolute_weight;
		greatest_integral_magnitude = s9_max(greatest_integral_magnitude, integral_magnitude);
	}

	last_share = first_share + step * (double)(S9_SCAN_SAMPLE_COUNT - 1U);
	absolute_first_share = s9_abs(first_share);
	absolute_last_share = s9_abs(last_share);
	greatest_share = s9_max(absolute_first_share, absolute_last_share);
	absolute_start_angle = s9_abs(shape->start_angle);
	absolute_turn = s9_abs(shape->turn);
	phase_scale = 1.0 + absolute_start_angle + (absolute_turn + 2.0 * greatest_share) * greatest_integral_magnitude;
	weight_scale = s9_max(1.0, weight_magnitude);
	sign_allowance = 64.0 * DBL_EPSILON * (double)(grid->count + S9_SCAN_RESTART_INTERVAL + 8U) * phase_scale * weight_scale;

	for (sample = 0U; sample < S9_SCAN_SAMPLE_COUNT; sample++) {
		double const share = first_share + step * (double)sample;
		double horizontal = 0.0;
		double vertical = 0.0;
		double absolute_horizontal, absolute_vertical;

		if (sample % S9_SCAN_RESTART_INTERVAL == 0U) {
			s9_set_lines(share, shape);

			for (index = 0U; index < grid->count; index++) {
				double const angle = shape->start_angle + shape->start_line * grid->start_integral[index] +
					shape->end_line * grid->end_integral[index];
				double const cosine = cos(angle);
				double const sine = sin(angle);
				double const weight = grid->weight[index];

				cosines[index] = cosine;
				sines[index] = sine;
				horizontal += weight * cosine;
				vertical += weight * sine;
			}
		} else {
			for (index = 0U; index < grid->count; index++) {
				double const cosine = cosines[index] * step_cosines[index] - sines[index] * step_sines[index];
				double const sine = sines[index] * step_cosines[index] + cosines[index] * step_sines[index];
				double const weight = grid->weight[index];

				cosines[index] = cosine;
				sines[index] = sine;
				horizontal += weight * cosine;
				vertical += weight * sine;
			}
		}

		absolute_horizontal = s9_abs(horizontal);
		absolute_vertical = s9_abs(vertical);

		if (!isfinite(sign_allowance) || !isfinite(horizontal) || !isfinite(vertical) ||
			absolute_horizontal <= sign_allowance || absolute_vertical <= sign_allowance) {
			s9_end_point const direct = s9_measure_end_point(grid, shape, share);

			points[sample] = direct.point;
		} else {
			points[sample] = s9_xy(horizontal, vertical);
		}
	}
}

static void s9_scan_share(s9_grid const* coarse, s9_grid const* fine, s9_curve* shape, s9_root_choice* choice)
{
	double const center = 0.5 * shape->turn;
	double const step = 2.0 * S9_SCAN_REACH / (double)(S9_SCAN_SAMPLE_COUNT - 1U);
	double previous_share = center - S9_SCAN_REACH;
	s9_point points[S9_SCAN_SAMPLE_COUNT];
	s9_point previous;
	unsigned int sample;

	s9_sample_share_scan(coarse, shape, points);

	previous = points[0U];

	for (sample = 1U; sample < S9_SCAN_SAMPLE_COUNT; sample++) {
		double const share = center - S9_SCAN_REACH + step * (double)sample;
		s9_point const current = points[sample];
		int const is_finite = isfinite(previous.y) && isfinite(current.y);
		int const brackets = is_finite && s9_parts_bracket_root(previous.y, current.y);

		if (brackets) {
			/* a bracket of the rotated samples is checked with direct measurements */
			s9_end_point const direct_previous = s9_measure_end_point(coarse, shape, previous_share);
			s9_end_point const direct_current = s9_measure_end_point(coarse, shape, share);
			int const direct_is_finite = isfinite(direct_previous.point.y) && isfinite(direct_current.point.y);
			int const faces_forward = direct_previous.point.x > 0.0 || direct_current.point.x > 0.0;
			int const direct_brackets = direct_is_finite && s9_parts_bracket_root(direct_previous.point.y, direct_current.point.y);

			if (faces_forward && direct_brackets) {
				s9_examine_bracket(fine, shape, previous_share, share, step, choice);
			}
		}

		previous_share = share;
		previous = current;
	}
}

static int s9_try_predicted_root(s9_grid const* fine, s9_curve* shape, double* share, s9_point* end_point)
{
	double const start_basis_turn = s9_start_basis_turn(shape);
	double const end_basis_turn = s9_end_basis_turn(shape);
	double start_sum = 0.0;
	double end_sum = 0.0;
	double rate, predicted, current;
	unsigned int index;

	for (index = 0U; index < fine->count; index++) {
		start_sum += fine->weight[index] * fine->start_integral[index];
		end_sum += fine->weight[index] * fine->end_integral[index];
	}

	rate = start_sum / start_basis_turn - end_sum / end_basis_turn;
	predicted = 0.5 * shape->turn;

	if (rate != 0.0) {
		predicted = -(shape->start_angle + shape->turn * end_sum / end_basis_turn) / rate;
	}

	if (!isfinite(predicted)) {
		return 0;
	}

	current = predicted;

	for (index = 0U; index < S9_PREDICTED_STEP_COUNT; index++) {
		s9_end_point const measured = s9_measure_end_point(fine, shape, current);
		double const value = measured.point.y;
		double slope, next, travel;

		if (!isfinite(value)) {
			return 0;
		}

		if (s9_abs(value) <= S9_ROOT_TOLERANCE) {
			double const distance = hypot(measured.point.x, measured.point.y);
			int const is_ahead = measured.point.x > 0.0 && distance > S9_SMALLEST_CHORD_LENGTH;
			double const offset = s9_abs(current - 0.5 * shape->turn);
			int const is_central = offset <= S9_PREDICTED_REACH;

			if (is_ahead && is_central) {
				*share = current;
				*end_point = measured.point;

				return 1;
			}

			return 0;
		}

		slope = measured.slope.y;

		if (slope == 0.0 || !isfinite(slope)) {
			return 0;
		}

		next = current - value / slope;
		travel = s9_abs(next - predicted);

		if (!isfinite(next) || travel > S9_PREDICTED_REACH) {
			return 0;
		}

		current = next;
	}

	return 0;
}

/* solves the curve for given weights; true when a curve that meets both tangent angles was found */
static int s9_try_solve_with_weights(s9_curve* self, s9_curve_parameters const* parameters)
{
	s9_curve candidate = s9_curve_default();
	s9_grid coarse;
	s9_grid fine;
	double const full_turn = 2.0 * S9_PI;
	double principal_turn, raw_turn, best_turn;
	double turns[4];
	s9_root_choice choice;
	unsigned int index;

	s9_prepare_shape(parameters, &candidate);

	s9_build_grid(&candidate, S9_COARSE_BASE_PANEL_COUNT, 8U, 0, &coarse);
	s9_build_grid(&candidate, S9_BASE_PANEL_COUNT, 16U, 1, &fine);

	principal_turn = candidate.turn;
	raw_turn = s9_raw_turn(&candidate);

	if (s9_abs(candidate.start_angle) <= S9_PREDICTED_ANGLE_LIMIT && s9_abs(candidate.end_angle) <= S9_PREDICTED_ANGLE_LIMIT) {
		double predicted_share = 0.0;
		s9_point predicted_end = s9_xy(0.0, 0.0);

		candidate.turn = raw_turn;

		if (s9_try_predicted_root(&fine, &candidate, &predicted_share, &predicted_end)) {
			s9_set_lines(predicted_share, &candidate);
			s9_complete_shape(predicted_end, &candidate);
			*self = candidate;

			return 1;
		}

		candidate.turn = principal_turn;
	}

	turns[0] = raw_turn;
	turns[1] = principal_turn;
	turns[2] = principal_turn - full_turn;
	turns[3] = principal_turn + full_turn;

	choice.share = 0.0;
	choice.distance = 0.0;
	choice.turning = 0.0;
	choice.end_point = s9_xy(0.0, 0.0);
	choice.found = 0;

	best_turn = raw_turn;

	for (index = 0U; index < 4U && !choice.found; index++) {
		double const turn = turns[index];
		double const difference = s9_abs(turn - raw_turn);

		if (index > 0U && difference <= S9_SAME_TURN_TOLERANCE) {
			continue;
		}

		candidate.turn = turn;
		s9_scan_share(&coarse, &fine, &candidate, &choice);
		best_turn = turn;
	}

	if (!choice.found || choice.distance <= S9_SMALLEST_CHORD_LENGTH) {
		return 0;
	}

	candidate.turn = best_turn;
	s9_set_lines(choice.share, &candidate);
	s9_complete_shape(choice.end_point, &candidate);
	*self = candidate;

	return 1;
}

/* solves the curve for given weights from a nearby solved curve */
static int s9_try_solve_with_weights_near(s9_curve* self, s9_curve_parameters const* parameters, s9_curve const* nearby)
{
	s9_curve candidate = s9_curve_default();
	double start_change, end_change, raw_turn, principal_difference, raw_difference;
	double branch_turn, branch_difference, nearby_share;
	int angles_are_near;
	int is_solved = 0;

	s9_prepare_shape(parameters, &candidate);

	start_change = candidate.start_angle - nearby->start_angle;
	end_change = candidate.end_angle - nearby->end_angle;
	raw_turn = s9_raw_turn(&candidate);
	principal_difference = s9_abs(nearby->turn - candidate.turn);
	raw_difference = s9_abs(nearby->turn - raw_turn);
	branch_turn = candidate.turn;
	branch_difference = principal_difference;

	if (raw_difference < principal_difference) {
		branch_turn = raw_turn;
		branch_difference = raw_difference;
	}

	angles_are_near = s9_abs(start_change) <= S9_GREATEST_NEAR_ANGLE_CHANGE &&
		s9_abs(end_change) <= S9_GREATEST_NEAR_ANGLE_CHANGE &&
		branch_difference < S9_PI;
	nearby_share = s9_start_share(nearby);

	if (angles_are_near && isfinite(nearby_share)) {
		s9_grid fine;
		double share = nearby_share;
		unsigned int index;

		s9_build_grid(&candidate, S9_BASE_PANEL_COUNT, 16U, 0, &fine);
		candidate.turn = branch_turn;

		for (index = 0U; index < S9_GREATEST_NEAR_STEP_COUNT; index++) {
			s9_end_point const measured = s9_measure_end_point(&fine, &candidate, share);
			double const value = measured.point.y;
			double slope, next, distance;

			if (!isfinite(value)) {
				break;
			}

			if (s9_abs(value) <= S9_ROOT_TOLERANCE) {
				double const chord = hypot(measured.point.x, measured.point.y);

				if (measured.point.x > 0.0 && chord > S9_SMALLEST_CHORD_LENGTH) {
					s9_complete_shape(measured.point, &candidate);
					is_solved = 1;
				}

				break;
			}

			slope = measured.slope.y;

			if (slope == 0.0 || !isfinite(slope)) {
				break;
			}

			next = share - value / slope;
			distance = s9_abs(next - nearby_share);

			if (!isfinite(next) || distance > S9_NEAR_REACH) {
				break;
			}

			share = next;
		}
	}

	if (is_solved) {
		*self = candidate;

		return 1;
	}

	return s9_try_solve_with_weights(self, parameters);
}

static s9_point s9_measure_panel(s9_curve const* shape, double begin, double width)
{
	double horizontal = 0.0;
	double vertical = 0.0;
	unsigned int index;

	for (index = 0U; index < 16U; index++) {
		double const t = s9_min(begin + width * s9_gauss_legendre_x_16[index], 1.0);
		double const angle = s9_tangent_angle(shape, t);
		double const weight = width * s9_gauss_legendre_scale_16[index];

		horizontal += weight * cos(angle);
		vertical += weight * sin(angle);
	}

	return s9_xy(horizontal * shape->length, vertical * shape->length);
}

static void s9_pieces_curve_to(s9_pieces* pieces, s9_point first, s9_point second, s9_point last)
{
	pieces->points[pieces->count] = first;
	pieces->points[pieces->count + 1U] = second;
	pieces->points[pieces->count + 2U] = last;
	pieces->count += 3U;
}

static void s9_add_cubic(s9_approximation* approximation, s9_point start_point, s9_point end_point, double start_angle, double end_angle,
	double arc_length)
{
	double const handle = arc_length / 3.0;
	s9_point const start_direction = s9_xy(cos(start_angle), sin(start_angle));
	s9_point const end_direction = s9_xy(cos(end_angle), sin(end_angle));
	s9_point const first_control = s9_xy(start_point.x + start_direction.x * handle, start_point.y + start_direction.y * handle);
	s9_point const second_control = s9_xy(end_point.x - end_direction.x * handle, end_point.y - end_direction.y * handle);

	s9_pieces_curve_to(approximation->pieces, first_control, second_control, end_point);
	approximation->remaining_cubics--;
}

static void s9_approximate_range(s9_approximation* approximation, double begin, double end, unsigned int depth, s9_point start_point,
	s9_point* end_point)
{
	s9_curve const* const shape = approximation->shape;
	double const width = end - begin;
	double const middle = begin + 0.5 * width;
	double const quarter_width = 0.25 * width;
	double const start_angle = s9_tangent_angle(shape, begin);
	s9_point piece_point = start_point;
	s9_point sample_points[3];
	double sample_angles[3];
	double variation = 0.0;
	double piece_angle = start_angle;
	double end_angle, split_turn, turn_magnitude;
	s9_point split_displacement;
	int needs_split = 0;
	unsigned int index;

	sample_points[0] = start_point;
	sample_points[1] = start_point;
	sample_points[2] = start_point;
	sample_angles[0] = start_angle;
	sample_angles[1] = start_angle;
	sample_angles[2] = start_angle;

	for (index = 0U; index < 4U; index++) {
		double const piece_begin = begin + quarter_width * (double)index;
		double const piece_end = s9_min(piece_begin + quarter_width, 1.0);
		s9_point const displacement = s9_measure_panel(shape, piece_begin, quarter_width);
		double const next_angle = s9_tangent_angle(shape, piece_end);

		piece_point.x += displacement.x;
		piece_point.y += displacement.y;
		variation += s9_abs(next_angle - piece_angle);
		piece_angle = next_angle;

		if (index < 3U) {
			sample_points[index] = piece_point;
			sample_angles[index] = piece_angle;
		}
	}

	end_angle = piece_angle;
	split_turn = end_angle - start_angle;
	split_displacement = s9_sub(piece_point, start_point);
	turn_magnitude = s9_max(variation, s9_abs(split_turn));

	if (turn_magnitude > S9_GREATEST_PANEL_TURN) {
		needs_split = 1;
	}

	if (!needs_split) {
		s9_point const whole_displacement = s9_measure_panel(shape, begin, width);
		s9_point const difference = s9_sub(whole_displacement, split_displacement);
		double const error = hypot(difference.x, difference.y);
		double const share = S9_INTEGRATION_ALLOWANCE * width;
		double const scale = width + s9_abs(start_point.x) + s9_abs(start_point.y) + s9_abs(piece_point.x) + s9_abs(piece_point.y);
		double const floor_allowance = 32.0 * DBL_EPSILON * scale;
		double const panel_allowance = s9_max(share, floor_allowance);

		if (error > panel_allowance) {
			needs_split = 1;
		}
	}

	/* the cubic is compared at three interior points, because it can pass through the middle point and still leave the curve */
	if (!needs_split) {
		double const arc_length = width * shape->length;
		double const handle = arc_length / 3.0;
		s9_point const start_direction = s9_xy(cos(start_angle), sin(start_angle));
		s9_point const end_direction = s9_xy(cos(end_angle), sin(end_angle));
		s9_point const first_control = s9_xy(start_point.x + start_direction.x * handle, start_point.y + start_direction.y * handle);
		s9_point const second_control = s9_xy(piece_point.x - end_direction.x * handle, piece_point.y - end_direction.y * handle);
		double parameters[3];

		parameters[0] = 0.25;
		parameters[1] = 0.5;
		parameters[2] = 0.75;

		for (index = 0U; index < 3U; index++) {
			double const t = parameters[index];
			double const s = 1.0 - t;
			double const start_share = s * s * s;
			double const first_share = 3.0 * s * s * t;
			double const second_share = 3.0 * s * t * t;
			double const end_share = t * t * t;
			s9_point const on_cubic = s9_xy(
				start_point.x * start_share + first_control.x * first_share + second_control.x * second_share + piece_point.x * end_share,
				start_point.y * start_share + first_control.y * first_share + second_control.y * second_share + piece_point.y * end_share);
			double const sample_angle = sample_angles[index];
			s9_point const across = s9_xy(-sin(sample_angle), cos(sample_angle));
			s9_point const cubic_error = s9_sub(on_cubic, sample_points[index]);
			double const cross_distance = cubic_error.x * across.x + cubic_error.y * across.y;

			if (s9_abs(cross_distance) > approximation->maximum_error) {
				needs_split = 1;
			}
		}
	}

	if (approximation->remaining_measurements > S9_MEASUREMENTS_PER_RANGE) {
		approximation->remaining_measurements -= S9_MEASUREMENTS_PER_RANGE;
	} else {
		approximation->remaining_measurements = 0U;
		approximation->failed = 1;
	}

	if (needs_split && depth < S9_GREATEST_DEPTH && !approximation->failed) {
		s9_point inner_point = s9_xy(0.0, 0.0);

		s9_approximate_range(approximation, begin, middle, depth + 1U, start_point, &inner_point);
		s9_approximate_range(approximation, middle, end, depth + 1U, inner_point, end_point);

		return;
	}

	*end_point = piece_point;

	if (approximation->failed) {
		return;
	}

	if (approximation->remaining_cubics == 0U) {
		approximation->failed = 1;

		return;
	}

	s9_add_cubic(approximation, start_point, *end_point, start_angle, end_angle, width * shape->length);
}

/* the cubics of a solved curve, in unit chord coordinates */
static int s9_try_approximate_with_cubics(s9_curve const* solved, double maximum_error, s9_pieces* pieces)
{
	s9_point const start_point = s9_xy(0.0, 0.0);
	s9_point const last_point = s9_xy(1.0, 0.0);
	s9_approximation approximation;
	s9_point final_point = s9_xy(0.0, 0.0);
	unsigned int last_index;
	s9_point first_edge, second_edge, start_direction, end_direction;
	double first_handle, second_handle;

	if (!(solved->length > 0.0)) {
		return 0;
	}

	approximation.shape = solved;
	approximation.pieces = pieces;
	approximation.maximum_error = maximum_error;
	approximation.remaining_measurements = S9_GREATEST_PANEL_COUNT;
	approximation.remaining_cubics = S9_GREATEST_CUBIC_COUNT;
	approximation.failed = 0;

	pieces->points[0] = start_point;
	pieces->count = 1U;

	s9_approximate_range(&approximation, 0.0, 1.0, 0U, start_point, &final_point);

	if (approximation.failed || pieces->count <= 1U) {
		return 0;
	}

	/* the ends take the curve's exact tangents and the chord's exact end point */
	last_index = pieces->count - 1U;
	first_edge = s9_sub(pieces->points[1U], start_point);
	second_edge = s9_sub(pieces->points[last_index - 1U], last_point);
	first_handle = hypot(first_edge.x, first_edge.y);
	second_handle = hypot(second_edge.x, second_edge.y);
	start_direction = s9_xy(cos(solved->start_angle), sin(solved->start_angle));
	end_direction = s9_xy(cos(solved->end_angle), sin(solved->end_angle));

	pieces->points[1U] = s9_xy(start_point.x + start_direction.x * first_handle, start_point.y + start_direction.y * first_handle);
	pieces->points[last_index - 1U] = s9_xy(last_point.x - end_direction.x * second_handle, last_point.y - end_direction.y * second_handle);
	pieces->points[last_index] = last_point;

	return 1;
}

/* the knot angle conversions: alpha is relative to the segment leaving the knot, beta to the one arriving */
static double s9_absolute_to_alpha(double knot_angle, double segment_angle) { return knot_angle - segment_angle; }
static double s9_alpha_to_absolute(double knot_angle, double segment_angle) { return knot_angle + segment_angle; }
static double s9_absolute_to_beta(double knot_angle, double segment_angle) { return segment_angle - knot_angle; }
static double s9_beta_to_absolute(double knot_angle, double segment_angle) { return -knot_angle + segment_angle; }

/* to -PI..PI */
static double s9_normalize_angle(double radians)
{
	double const pi2 = S9_PI + S9_PI;
	double const pi2_r = 1.0 / pi2;
	double const r = radians * pi2_r;

	return pi2 * (r - round(r));
}

static s9_point s9_rotate_point_around_origin(s9_point v, double radians)
{
	double const s = sin(radians);
	double const c = cos(radians);
	double const x = v.x * c + v.y * s;
	double const y = -v.x * s + v.y * c;

	return s9_xy(x, y);
}

static double s9_pole_axis(double beta)
{
	double const k = 0.3;
	double const beta2 = beta * 2.0;
	double const v = sinh((1.0 + k) / (2.0 * k) * atan(k * sin(beta2) / (1.0 + k * cos(beta2))));

	return v / 2.0;
}

/* the rotation range of a knot: the absolute angles between the curvature pole of the segment before it and that of the segment after it */
static void s9_calculate_rotation_range(double prev_knot_angle, double next_knot_angle, double prev_segment_angle, double next_segment_angle,
	double* min_angle, double* max_angle)
{
	double const beta_limit = s9_beta_to_absolute(s9_pole_axis(s9_absolute_to_alpha(prev_knot_angle, prev_segment_angle)), prev_segment_angle);
	double const alpha_limit = s9_alpha_to_absolute(s9_pole_axis(s9_absolute_to_beta(next_knot_angle, next_segment_angle)), next_segment_angle);
	double const axis = (beta_limit + alpha_limit) / 2.0;
	double const range = (S9_PI - s9_abs(s9_normalize_angle(beta_limit - alpha_limit))) / 2.0;

	*min_angle = axis - range;
	*max_angle = axis + range;
}

static s9_curve const* s9_known_shape(s9_segment_curve const* curve)
{
	if (curve == NULL || !curve->is_known) {
		return NULL;
	}

	return &curve->shape;
}

static void s9_keep_shape(s9_curve const* shape, s9_segment_curve* curve)
{
	if (curve == NULL) {
		return;
	}

	curve->shape = *shape;
	curve->is_known = 1;
}

/* the signed alignment of one endpoint; at the right end alpha and beta are clockwise from the chord */
static double s9_end_alignment(double alpha, double beta)
{
	double const t = 0.9;
	double const axis = s9_pole_axis(beta);
	double const phase = alpha - axis;

	return tanh(t * cos(phase)) / tanh(t);
}

static double s9_alignment_weight(double alignment)
{
	return s9_max(pow(s9_abs(alignment), 2.8), 1.7e-10);
}

static double s9_weight_eta(double weight)
{
	double const ww = weight * weight;
	double const t = S9_PI * sqrt(5.0) * weight / sqrt(1.0 + 4.0 * ww);
	double const s = sin(t);

	return s * s + 1.5;
}

static s9_curve_parameters s9_parameters_from_angles(double alpha, double beta)
{
	double const start_alignment = s9_end_alignment(alpha, beta);
	double const end_alignment = s9_end_alignment(beta, alpha);
	double const start_alignment_weight = s9_alignment_weight(start_alignment);
	double const end_alignment_weight = s9_alignment_weight(end_alignment);
	/* one for two flat ends, near two when either end is a corner */
	double const middle = 2.0 - start_alignment_weight * end_alignment_weight;
	s9_curve_parameters parameters;

	parameters.alpha = alpha;
	parameters.beta = beta;
	parameters.start_weight = start_alignment_weight / middle;
	parameters.end_weight = end_alignment_weight / middle;
	parameters.start_eta = s9_weight_eta(parameters.start_weight);
	parameters.end_eta = s9_weight_eta(parameters.end_weight);

	if (start_alignment < 0.0) {
		parameters.start_weight = -parameters.start_weight;
	}

	if (end_alignment < 0.0) {
		parameters.end_weight = -parameters.end_weight;
	}

	return parameters;
}

/* solves a curve from nothing, or continuing nearby when it is given */
static int s9_try_solve_curve(double alpha, double beta, s9_curve const* nearby, s9_curve* curve)
{
	s9_curve_parameters const parameters = s9_parameters_from_angles(alpha, beta);

	if (nearby != NULL) {
		return s9_try_solve_with_weights_near(curve, &parameters, nearby);
	}

	return s9_try_solve_with_weights(curve, &parameters);
}

static int s9_try_calculate_end_curvatures(double alpha, double beta, s9_curve const* nearby, s9_curve* curve, double* start_curvature,
	double* end_curvature)
{
	s9_curve solved = s9_curve_default();

	if (!s9_try_solve_curve(alpha, beta, nearby, &solved)) {
		return 0;
	}

	*start_curvature = s9_curvature(&solved, 0.0);
	*end_curvature = s9_curvature(&solved, 1.0);
	*curve = solved;

	return 1;
}

/* The curvature delta at the center knot of two neighbouring segments: curvature at the start of the next minus curvature
   at the end of the previous, each for its own chord. It becomes very large when a weight at the center knot goes to zero.
   x, y and z are the absolute angles of the previous, center and next knots, a0 and a1 of the two segments, h0 and h1 the
   chords. The optional last curves of the two segments are continued when known, and both are replaced when both solves
   succeed. NaN when a solve fails. */
static double s9_calculate_curvature_delta(double x, double y, double z, double a0, double a1, double h0, double h1,
	s9_segment_curve* prev_curve, s9_segment_curve* next_curve)
{
	double const alpha0 = s9_absolute_to_alpha(x, a0);
	double const beta0 = s9_absolute_to_beta(y, a0);
	double const alpha1 = s9_absolute_to_alpha(y, a1);
	double const beta1 = s9_absolute_to_beta(z, a1);
	double start_curvature0 = 0.0;
	double unit_curvature0 = 0.0;
	double unit_curvature1 = 0.0;
	double end_curvature1 = 0.0;
	s9_curve const* const nearby_shape0 = s9_known_shape(prev_curve);
	s9_curve const* const nearby_shape1 = s9_known_shape(next_curve);
	s9_curve shape0 = s9_curve_default();
	s9_curve shape1 = s9_curve_default();
	double curvature0, curvature1;

	if (!s9_try_calculate_end_curvatures(alpha0, beta0, nearby_shape0, &shape0, &start_curvature0, &unit_curvature0)) {
		return NAN;
	}

	if (!s9_try_calculate_end_curvatures(alpha1, beta1, nearby_shape1, &shape1, &unit_curvature1, &end_curvature1)) {
		return NAN;
	}

	s9_keep_shape(&shape0, prev_curve);
	s9_keep_shape(&shape1, next_curve);

	/* curvature is inversely proportional to the chord length */
	curvature0 = unit_curvature0 / h0;
	curvature1 = unit_curvature1 / h1;

	return curvature1 - curvature0;
}

/* the middle of a knot's rotation range, which is the start angle of its solver */
static double s9_clamp_knot_angle(double prev_knot_angle, double next_knot_angle, double prev_segment_angle, double next_segment_angle)
{
	double min_angle, max_angle, sum;

	s9_calculate_rotation_range(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, &min_angle, &max_angle);

	sum = min_angle + max_angle;

	return sum * 0.5;
}

/* The fallback of s9_solve_knot_angle: the knot angle that makes the curvature delta zero, with 32 bisection steps over the
   whole rotation range. The sign near the start of the range is read with solves from nothing; the steps continue the
   last curves of the two segments when given. */
static double s9_bisection_angle(double prev_knot_angle, double next_knot_angle, double prev_segment_angle, double next_segment_angle,
	double prev_segment_length, double next_segment_length, s9_segment_curve* prev_curve, s9_segment_curve* next_curve)
{
	double lo, hi, sum, angle, range, probe_step, probe, probe_value;
	int sa, step;

	s9_calculate_rotation_range(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, &lo, &hi);

	sum = lo + hi;
	angle = sum / 2.0;

	range = hi - lo;
	probe_step = s9_min(S9_SIGN_PROBE_STEP, 0.25 * range);
	probe = lo + probe_step;
	probe_value = s9_calculate_curvature_delta(prev_knot_angle, probe, next_knot_angle, prev_segment_angle, next_segment_angle,
		prev_segment_length, next_segment_length, NULL, NULL);

	sa = (int)copysign(1.0, probe_value);

	for (step = 0; step < 32; step++) {
		double const y = s9_calculate_curvature_delta(prev_knot_angle, angle, next_knot_angle, prev_segment_angle, next_segment_angle,
			prev_segment_length, next_segment_length, prev_curve, next_curve);
		int const s = (int)copysign(1.0, y);

		if (sa == s) {
			lo = angle;
		} else {
			hi = angle;
		}

		sum = lo + hi;
		angle = sum / 2.0;
	}

	return angle;
}

/* true when a finite curvature delta is on the low side of the root of a knot, for the sign of the deltas near the start of the range */
static int s9_delta_is_below_root(double value, int low_is_negative)
{
	int const is_negative = signbit(value) != 0;

	return is_negative == low_is_negative;
}

/* The knot angle that makes the curvature delta zero, starting at the knot's current angle. The first call for a knot
   reads the sign of the delta near the start of the rotation range with solves from nothing and keeps it in low_sign. It
   then searches for a bracket of the root, with steps that start at twice the knot's previous update and grow 4 times per
   step, and refines the bracket with the Illinois form of the false position method, the middle of the bracket every 4th
   step. All deltas after the sign continue the last curves of the two segments. */
static double s9_solve_knot_angle(double prev_knot_angle, double current_knot_angle, double next_knot_angle, double prev_segment_angle,
	double next_segment_angle, double prev_segment_length, double next_segment_length, double previous_update,
	s9_segment_curve* prev_curve, s9_segment_curve* next_curve, signed char* low_sign)
{
	double min_angle, max_angle, range, probe_step, low, high;
	double guess, guess_value, first_step, step, lower, lower_value, upper, upper_value;
	int low_is_negative, guess_is_below;
	int has_bracket = 0;
	int last_side = 0;
	unsigned int index;

	s9_calculate_rotation_range(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, &min_angle, &max_angle);

	range = max_angle - min_angle;
	probe_step = s9_min(S9_SIGN_PROBE_STEP, 0.25 * range);
	low = min_angle + probe_step;
	high = max_angle - probe_step;

	if (!(high > low)) {
		*low_sign = 0;

		return s9_bisection_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, prev_segment_length,
			next_segment_length, prev_curve, next_curve);
	}

	/* solves from nothing: the last curves are no good start at a pole */
	if (*low_sign == 0) {
		double const low_value = s9_calculate_curvature_delta(prev_knot_angle, low, next_knot_angle, prev_segment_angle, next_segment_angle,
			prev_segment_length, next_segment_length, NULL, NULL);

		if (!isfinite(low_value)) {
			return s9_bisection_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, prev_segment_length,
				next_segment_length, prev_curve, next_curve);
		}

		*low_sign = 1;

		if (signbit(low_value)) {
			*low_sign = -1;
		}
	}

	low_is_negative = *low_sign < 0;
	guess = s9_clamp(current_knot_angle, low, high);
	guess_value = s9_calculate_curvature_delta(prev_knot_angle, guess, next_knot_angle, prev_segment_angle, next_segment_angle,
		prev_segment_length, next_segment_length, prev_curve, next_curve);

	if (!isfinite(guess_value)) {
		*low_sign = 0;

		return s9_bisection_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, prev_segment_length,
			next_segment_length, prev_curve, next_curve);
	}

	/* search for a bracket: lower is below the root and upper above it */
	first_step = S9_FIRST_BRACKET_STEP;

	if (previous_update != 0.0) {
		double const double_update = 2.0 * s9_abs(previous_update);

		first_step = s9_max(double_update, S9_SMALLEST_BRACKET_STEP);
	}

	step = s9_min(first_step, 0.25 * range);
	lower = guess;
	lower_value = guess_value;
	upper = guess;
	upper_value = guess_value;
	guess_is_below = s9_delta_is_below_root(guess_value, low_is_negative);

	for (index = 0U; index < S9_GREATEST_ROOT_ITERATION_COUNT; index++) {
		double candidate = s9_max(upper - step, low);
		double candidate_value;
		int is_below, is_at_end;

		if (guess_is_below) {
			candidate = s9_min(lower + step, high);
		}

		candidate_value = s9_calculate_curvature_delta(prev_knot_angle, candidate, next_knot_angle, prev_segment_angle, next_segment_angle,
			prev_segment_length, next_segment_length, prev_curve, next_curve);

		if (!isfinite(candidate_value)) {
			break;
		}

		is_below = s9_delta_is_below_root(candidate_value, low_is_negative);

		if (is_below) {
			lower = candidate;
			lower_value = candidate_value;
		} else {
			upper = candidate;
			upper_value = candidate_value;
		}

		if (is_below != guess_is_below) {
			has_bracket = 1;
			break;
		}

		is_at_end = candidate <= low || candidate >= high;

		if (is_at_end) {
			break;
		}

		step *= S9_BRACKET_STEP_GROWTH;
	}

	if (!has_bracket) {
		*low_sign = 0;

		return s9_bisection_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, prev_segment_length,
			next_segment_length, prev_curve, next_curve);
	}

	/* Illinois steps: the value at an end that holds for two steps is halved */
	for (index = 0U; index < S9_GREATEST_ROOT_ITERATION_COUNT; index++) {
		double const width = upper - lower;
		double middle, candidate, value_span, candidate_value;
		int takes_middle, is_below;

		if (width <= S9_KNOT_ANGLE_TOLERANCE) {
			break;
		}

		middle = lower + 0.5 * width;
		candidate = middle;
		takes_middle = (index % 4U) == 3U;
		value_span = upper_value - lower_value;

		if (!takes_middle && value_span != 0.0) {
			double const share = -lower_value / value_span;
			double const linear = lower + share * width;

			if (linear > lower && linear < upper) {
				candidate = linear;
			}
		}

		candidate_value = s9_calculate_curvature_delta(prev_knot_angle, candidate, next_knot_angle, prev_segment_angle, next_segment_angle,
			prev_segment_length, next_segment_length, prev_curve, next_curve);

		if (!isfinite(candidate_value)) {
			*low_sign = 0;

			return s9_bisection_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle, prev_segment_length,
				next_segment_length, prev_curve, next_curve);
		}

		if (candidate_value == 0.0) {
			return candidate;
		}

		is_below = s9_delta_is_below_root(candidate_value, low_is_negative);

		if (is_below) {
			lower = candidate;
			lower_value = candidate_value;

			if (last_side < 0) {
				upper_value *= 0.5;
			}

			last_side = -1;
		} else {
			upper = candidate;
			upper_value = candidate_value;

			if (last_side > 0) {
				lower_value *= 0.5;
			}

			last_side = 1;
		}
	}

	return lower + 0.5 * (upper - lower);
}

static double s9_acyclic_start_angle(double first_segment_angle, double second_knot_angle)
{
	double const beta = s9_absolute_to_beta(second_knot_angle, first_segment_angle);
	double const alpha = 0.5 * sin(2.0 * beta);

	return s9_alpha_to_absolute(alpha, first_segment_angle);
}

static double s9_acyclic_end_angle(double final_segment_angle, double first_knot_angle)
{
	double const alpha = s9_absolute_to_alpha(first_knot_angle, final_segment_angle);
	double const beta = 0.5 * sin(2.0 * alpha);

	return s9_beta_to_absolute(beta, final_segment_angle);
}

/* The next angle of a knot from the solved one. An update against the knot's previous update keeps only
   S9_OSCILLATION_DAMPING of it, so a knot does not move back and forth between two angles. An update of half a turn or more
   is kept whole: a knot angle and the same angle plus PI give the same curves. */
static double s9_calculate_damped_knot_angle(double current_angle, double solved_angle, double* previous_update)
{
	double const update = solved_angle - current_angle;
	int const turns_back = update * *previous_update < 0.0;
	int const is_small = s9_abs(update) < 0.5 * S9_PI;
	double applied_update = update;

	if (turns_back && is_small) {
		applied_update = update * S9_OSCILLATION_DAMPING;
	}

	*previous_update = applied_update;

	return current_angle + applied_update;
}

/* appends one cubic, from the current end point, growing the output when it is full; 0 on allocation failure */
static int s9_output_curve_to(s9_output* output, s9_point first, s9_point second, s9_point last)
{
	size_t const point_count = output->point_count;

	if (point_count + 3U > output->capacity) {
		size_t const capacity = output->capacity * 2U;
		double* const points = (double*)qaws_internal_malloc(qaws_internal_mul3(capacity, 2U, sizeof(double)));

		if (points == NULL) {
			return 0;
		}

		memcpy(points, output->points, point_count * 2U * sizeof(double));
		qaws_internal_free(output->points);
		output->points = points;
		output->capacity = capacity;
	}

	output->points[2U * point_count] = first.x;
	output->points[2U * point_count + 1U] = first.y;
	output->points[2U * point_count + 2U] = second.x;
	output->points[2U * point_count + 3U] = second.y;
	output->points[2U * point_count + 4U] = last.x;
	output->points[2U * point_count + 5U] = last.y;
	output->point_count = point_count + 3U;

	return 1;
}

static qaws_status s9_append_cubics(double alpha, double beta, s9_point start_point, s9_point end_point, double segment_angle,
	double segment_length, s9_segment_curve const* cache, s9_pieces* pieces, s9_output* output)
{
	s9_curve shape = s9_curve_default();
	unsigned int index, last_index, curve_count;

	if (!s9_try_solve_curve(alpha, beta, s9_known_shape(cache), &shape)) {
		return QAWS_STATUS_NUMERICAL_FAILURE;
	}

	if (!s9_try_approximate_with_cubics(&shape, S9_MAXIMUM_CUBIC_ERROR, pieces)) {
		return QAWS_STATUS_NUMERICAL_FAILURE;
	}

	last_index = pieces->count - 1U;

	for (index = 1U; index < last_index; index++) {
		s9_point const point = pieces->points[index];
		s9_point const mirrored = s9_xy(point.x * segment_length, -point.y * segment_length);
		s9_point const rotated = s9_rotate_point_around_origin(mirrored, segment_angle);

		pieces->points[index] = s9_xy(start_point.x + rotated.x, start_point.y + rotated.y);
	}

	/* the spline must go through its knots exactly */
	pieces->points[last_index] = end_point;

	curve_count = last_index / 3U;

	for (index = 0U; index < curve_count; index++) {
		s9_point const* const cubic = pieces->points + index * 3U;

		if (!s9_output_curve_to(output, cubic[1], cubic[2], cubic[3])) {
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	}

	return QAWS_STATUS_OK;
}

static unsigned int s9_unknown_count(s9_knot_system const* system)
{
	if (system->is_cyclic) {
		return system->count;
	}

	return system->count - 2U;
}

static unsigned int s9_knot_of_unknown(s9_knot_system const* system, unsigned int unknown)
{
	if (system->is_cyclic) {
		return unknown;
	}

	return unknown + 1U;
}

/* the segment that ends at a knot; for an open spline the knot must not be zero */
static unsigned int s9_segment_before_knot(s9_knot_system const* system, unsigned int knot)
{
	if (knot == 0U) {
		return system->count - 1U;
	}

	return knot - 1U;
}

static unsigned int s9_end_knot_of_segment(s9_knot_system const* system, unsigned int segment)
{
	unsigned int const next = segment + 1U;

	if (system->is_cyclic && next == system->count) {
		return 0U;
	}

	return next;
}

static double s9_start_knot_angle(s9_knot_system const* system, double const* angles)
{
	return s9_acyclic_start_angle(system->segment_angles[0], angles[1]);
}

static double s9_end_knot_angle(s9_knot_system const* system, double const* angles)
{
	unsigned int const last = system->count - 1U;

	return s9_acyclic_end_angle(system->segment_angles[last - 1U], angles[last - 1U]);
}

/* sets the two end knots of an open spline from the end rules */
static void s9_apply_end_rules(s9_knot_system const* system, double* angles)
{
	unsigned int last;

	if (system->is_cyclic) {
		return;
	}

	last = system->count - 1U;

	angles[0] = s9_start_knot_angle(system, angles);
	angles[last] = s9_end_knot_angle(system, angles);
}

/* solves one segment for a set of knot angles and returns its end curvatures for its chord */
static int s9_solve_segment(s9_knot_system const* system, double const* angles, unsigned int segment, s9_segment_curve* cache, int keep,
	double* start_curvature, double* end_curvature)
{
	unsigned int const start_knot = segment;
	unsigned int const end_knot = s9_end_knot_of_segment(system, segment);
	double const segment_angle = system->segment_angles[segment];
	double const length = system->segment_lengths[segment];
	double const alpha = s9_absolute_to_alpha(angles[start_knot], segment_angle);
	double const beta = s9_absolute_to_beta(angles[end_knot], segment_angle);
	s9_curve const* const nearby = s9_known_shape(cache);
	s9_curve shape = s9_curve_default();
	double unit_start = 0.0;
	double unit_end = 0.0;

	if (!s9_try_calculate_end_curvatures(alpha, beta, nearby, &shape, &unit_start, &unit_end)) {
		return 0;
	}

	if (keep) {
		s9_keep_shape(&shape, cache);
	}

	/* curvature is inversely proportional to the chord length */
	*start_curvature = unit_start / length;
	*end_curvature = unit_end / length;

	return 1;
}

/* solves every segment that has no stored curve yet, from nothing, and keeps the shapes */
static void s9_solve_unknown_segments(s9_knot_system const* system)
{
	unsigned int segment;

	for (segment = 0U; segment < system->segment_count; segment++) {
		s9_segment_curve* const cache = &system->segment_curves[segment];
		double start_curvature = 0.0;
		double end_curvature = 0.0;

		if (cache->is_known) {
			continue;
		}

		/* a failed solve leaves the segment unknown */
		s9_solve_segment(system, system->knot_angles, segment, cache, 1, &start_curvature, &end_curvature);
	}
}

/* solves all segments for a set of knot angles and calculates the curvature delta of every unknown knot */
static int s9_evaluate_knot_system(s9_knot_system const* system, double const* angles, int keep, double* start_curvatures,
	double* end_curvatures, double* deltas)
{
	unsigned int const unknown_count = s9_unknown_count(system);
	unsigned int segment, unknown;

	for (segment = 0U; segment < system->segment_count; segment++) {
		if (!s9_solve_segment(system, angles, segment, &system->segment_curves[segment], keep, &start_curvatures[segment], &end_curvatures[segment])) {
			return 0;
		}
	}

	for (unknown = 0U; unknown < unknown_count; unknown++) {
		unsigned int const knot = s9_knot_of_unknown(system, unknown);
		unsigned int const before = s9_segment_before_knot(system, knot);

		deltas[unknown] = start_curvatures[knot] - end_curvatures[before];
	}

	return 1;
}

/* the rotation range of one unknown knot, a probe step inside its ends; true when the range is not empty */
static int s9_calculate_knot_range(s9_knot_system const* system, double const* angles, unsigned int knot, double* low, double* high)
{
	unsigned int const before = s9_segment_before_knot(system, knot);
	unsigned int const after = knot;
	unsigned int const next_knot = s9_end_knot_of_segment(system, after);
	double min_angle, max_angle, width, probe_step;

	s9_calculate_rotation_range(angles[before], angles[next_knot], system->segment_angles[before], system->segment_angles[after],
		&min_angle, &max_angle);

	width = max_angle - min_angle;
	probe_step = s9_min(S9_SIGN_PROBE_STEP, 0.25 * width);

	*low = min_angle + probe_step;
	*high = max_angle - probe_step;

	return *high > *low;
}

static double s9_knot_range_width(s9_knot_system const* system, double const* angles, unsigned int knot)
{
	unsigned int const before = s9_segment_before_knot(system, knot);
	unsigned int const next_knot = s9_end_knot_of_segment(system, knot);
	double min_angle, max_angle;

	s9_calculate_rotation_range(angles[before], angles[next_knot], system->segment_angles[before], system->segment_angles[knot],
		&min_angle, &max_angle);

	return max_angle - min_angle;
}

/* Marks the unknown knots whose rotation range is narrower than S9_PINNED_RANGE_WIDTH as pinned, the others as free. Both
   ends of a pinned knot are so near their curvature poles that the end weights sit at their floor, where the curvature
   delta no longer changes with the knot angle and can have no root; such a knot stays in the middle of its range and
   follows its neighbours. True when a mark changed. */
static int s9_classify_knots(s9_knot_system const* system, double const* angles)
{
	unsigned int const unknown_count = s9_unknown_count(system);
	unsigned int unknown;
	int changed = 0;

	for (unknown = 0U; unknown < unknown_count; unknown++) {
		unsigned int const knot = s9_knot_of_unknown(system, unknown);
		double const width = s9_knot_range_width(system, angles, knot);
		int is_pinned = width < S9_PINNED_RANGE_WIDTH;

		/* a pinned knot stays pinned until its range widens clearly */
		if (system->pinned[knot] && width < S9_UNPIN_RANGE_WIDTH) {
			is_pinned = 1;
		}

		if (system->pinned[knot] != is_pinned) {
			changed = 1;
			system->pinned[knot] = (signed char)is_pinned;
		}
	}

	return changed;
}

/* Moves every pinned knot to the middle of its range, in knot order, the end rules applied before each pass and after the
   last. The middle depends on the neighbours, which can be pinned too, so several passes settle chains of pinned knots. */
static void s9_centre_pinned_knots(s9_knot_system const* system, double* angles)
{
	unsigned int const unknown_count = s9_unknown_count(system);
	unsigned int pass, unknown;

	for (pass = 0U; pass < S9_CENTRING_PASS_COUNT; pass++) {
		s9_apply_end_rules(system, angles);

		for (unknown = 0U; unknown < unknown_count; unknown++) {
			unsigned int const knot = s9_knot_of_unknown(system, unknown);
			unsigned int before, next_knot;

			if (!system->pinned[knot]) {
				continue;
			}

			before = s9_segment_before_knot(system, knot);
			next_knot = s9_end_knot_of_segment(system, knot);

			angles[knot] = s9_clamp_knot_angle(angles[before], angles[next_knot], system->segment_angles[before], system->segment_angles[knot]);
		}
	}

	s9_apply_end_rules(system, angles);
}

/* clamps every free knot into its rotation range and applies the end rules; a knot leaves its range when the ranges move with its neighbours */
static void s9_clamp_free_knots(s9_knot_system const* system, double* angles)
{
	unsigned int const unknown_count = s9_unknown_count(system);
	unsigned int unknown;

	for (unknown = 0U; unknown < unknown_count; unknown++) {
		unsigned int const knot = s9_knot_of_unknown(system, unknown);
		double low = 0.0;
		double high = 0.0;

		if (system->pinned[knot]) {
			continue;
		}

		if (s9_calculate_knot_range(system, angles, knot, &low, &high)) {
			angles[knot] = s9_clamp(angles[knot], low, high);
		}
	}

	s9_apply_end_rules(system, angles);
}

/* the greatest magnitude of the curvature deltas of the free knots: pinned knots have no root to reach */
static double s9_free_delta_norm(s9_knot_system const* system, unsigned int unknown_count, double const* deltas)
{
	double greatest = 0.0;
	unsigned int unknown;

	for (unknown = 0U; unknown < unknown_count; unknown++) {
		unsigned int const knot = s9_knot_of_unknown(system, unknown);

		if (system->pinned[knot]) {
			continue;
		}

		greatest = s9_max(greatest, s9_abs(deltas[unknown]));
	}

	return greatest;
}

/* The row of a pinned knot in the Newton system: the linearized constraint that the knot angle equals the middle of its
   rotation range. The middle is differentiated by finite differences with respect to the unknowns before and after the
   knot and the knot itself, with the end rules applied, because the end knots of an open spline follow their neighbours. */
static void s9_write_pinned_row(s9_knot_system const* system, unsigned int unknown_count, unsigned int unknown, double const* angles,
	double* trial_angles, double* subdiagonal, double* diagonal, double* superdiagonal, double* right_side)
{
	unsigned int const count = system->count;
	unsigned int const knot = s9_knot_of_unknown(system, unknown);
	unsigned int const before = s9_segment_before_knot(system, knot);
	unsigned int const next_knot = s9_end_knot_of_segment(system, knot);
	double const middle = s9_clamp_knot_angle(angles[before], angles[next_knot], system->segment_angles[before], system->segment_angles[knot]);
	int const has_previous = unknown > 0U || system->is_cyclic;
	int const has_next = unknown + 1U < unknown_count || system->is_cyclic;
	unsigned int const previous = (unknown + unknown_count - 1U) % unknown_count;
	unsigned int const next = (unknown + 1U) % unknown_count;
	unsigned int columns[3];
	int present[3];
	double slopes[3];
	unsigned int side;

	columns[0] = previous;
	columns[1] = unknown;
	columns[2] = next;
	present[0] = has_previous;
	present[1] = 1;
	present[2] = has_next;
	slopes[0] = 0.0;
	slopes[1] = 0.0;
	slopes[2] = 0.0;

	for (side = 0U; side < 3U; side++) {
		unsigned int column_knot, index;
		double moved;

		if (!present[side]) {
			continue;
		}

		column_knot = s9_knot_of_unknown(system, columns[side]);

		for (index = 0U; index < count; index++) {
			trial_angles[index] = angles[index];
		}

		trial_angles[column_knot] += S9_COUPLED_DIFFERENCE_STEP;

		s9_apply_end_rules(system, trial_angles);

		moved = s9_clamp_knot_angle(trial_angles[before], trial_angles[next_knot], system->segment_angles[before], system->segment_angles[knot]);

		slopes[side] = (moved - middle) / S9_COUPLED_DIFFERENCE_STEP;
	}

	diagonal[unknown] = 1.0 - slopes[1];
	right_side[unknown] = middle - angles[knot];

	if (has_previous) {
		subdiagonal[unknown] = -slopes[0];
	}

	if (has_next) {
		superdiagonal[unknown] = -slopes[2];
	}
}

/* solves the two segments at a knot for a set of knot angles without keeping the shapes */
static int s9_solve_segments_at_knot(s9_knot_system const* system, double const* angles, unsigned int knot, double* before_start,
	double* before_end, double* after_start, double* after_end)
{
	unsigned int const before = s9_segment_before_knot(system, knot);
	unsigned int const after = knot;

	if (!s9_solve_segment(system, angles, before, &system->segment_curves[before], 0, before_start, before_end)) {
		return 0;
	}

	return s9_solve_segment(system, angles, after, &system->segment_curves[after], 0, after_start, after_end);
}

/* One sweep over the knots in which each free knot takes one Newton step on its curvature delta, the slope from a finite
   difference, the step limited to a fraction of the knot's range and the angle clamped into it. Pinned knots are moved
   to the middle of their ranges first. The greatest update of a knot angle, or a negative value when a segment could not
   be solved or a slope was zero. */
static double s9_sweep_knot_newton_steps(s9_knot_system const* system, double* trial)
{
	unsigned int const count = system->count;
	unsigned int const unknown_count = s9_unknown_count(system);
	double* const angles = system->knot_angles;
	double greatest_update = 0.0;
	unsigned int unknown;

	s9_classify_knots(system, angles);
	s9_centre_pinned_knots(system, angles);

	for (unknown = 0U; unknown < unknown_count; unknown++) {
		unsigned int const knot = s9_knot_of_unknown(system, unknown);
		unsigned int before, after, index;
		double low = 0.0;
		double high = 0.0;
		double entering, current, delta, difference, trial_delta, slope, cap, update, next, wanted, moved, magnitude;
		double before_start = 0.0;
		double before_end = 0.0;
		double after_start = 0.0;
		double after_end = 0.0;
		double trial_before_start = 0.0;
		double trial_before_end = 0.0;
		double trial_after_start = 0.0;
		double trial_after_end = 0.0;

		/* a pinned knot follows its neighbours */
		if (system->pinned[knot]) {
			continue;
		}

		before = s9_segment_before_knot(system, knot);
		after = knot;

		if (!s9_calculate_knot_range(system, angles, knot, &low, &high)) {
			return -1.0;
		}

		entering = angles[knot];
		current = s9_clamp(entering, low, high);

		angles[knot] = current;

		s9_apply_end_rules(system, angles);

		if (!s9_solve_segment(system, angles, before, &system->segment_curves[before], 1, &before_start, &before_end)) {
			return -1.0;
		}

		if (!s9_solve_segment(system, angles, after, &system->segment_curves[after], 1, &after_start, &after_end)) {
			return -1.0;
		}

		delta = after_start - before_end;

		/* the finite difference goes towards the inside of the range */
		difference = S9_COUPLED_DIFFERENCE_STEP;

		if (current + difference > high) {
			difference = -difference;
		}

		for (index = 0U; index < count; index++) {
			trial[index] = angles[index];
		}

		trial[knot] = current + difference;

		s9_apply_end_rules(system, trial);

		if (!s9_solve_segments_at_knot(system, trial, knot, &trial_before_start, &trial_before_end, &trial_after_start, &trial_after_end)) {
			return -1.0;
		}

		trial_delta = trial_after_start - trial_before_end;
		slope = (trial_delta - delta) / difference;

		if (slope == 0.0 || !isfinite(slope)) {
			return -1.0;
		}

		cap = S9_COUPLED_RANGE_FRACTION * (high - low);
		update = s9_clamp(-delta / slope, -cap, cap);
		next = s9_clamp(current + update, low, high);
		wanted = s9_abs(update);
		moved = s9_abs(next - current);

		/* the update includes the clamp into the range */
		magnitude = s9_abs(next - entering);

		/* a knot held at an end of its range by an outward step counts the step it wanted as its update, so it does not look converged */
		if (moved == 0.0) {
			magnitude = s9_max(magnitude, wanted);
		}

		greatest_update = s9_max(greatest_update, magnitude);
		angles[knot] = next;

		s9_apply_end_rules(system, angles);
	}

	return greatest_update;
}

/* The Thomas algorithm. a is the subdiagonal (its first value is not read), c the superdiagonal (its last is not read), d
   the right side; superdiagonal and right_side are scratch of count doubles. x can be any of the inputs. */
static void s9_solve_tridiagonal(unsigned int count, double const* a, double const* b, double const* c, double const* d, double* x,
	double* superdiagonal, double* right_side)
{
	double previous_superdiagonal = c[0] / b[0];
	double previous_right_side = d[0] / b[0];
	double previous_solution;
	unsigned int index;

	superdiagonal[0] = previous_superdiagonal;
	right_side[0] = previous_right_side;

	for (index = 1U; index < count; index++) {
		double const denominator = b[index] - a[index] * previous_superdiagonal;

		previous_superdiagonal = c[index] / denominator;
		previous_right_side = (d[index] - a[index] * previous_right_side) / denominator;

		superdiagonal[index] = previous_superdiagonal;
		right_side[index] = previous_right_side;
	}

	/* back substitution from the last equation */
	previous_solution = previous_right_side;

	x[count - 1U] = previous_solution;

	for (index = count - 1U; index > 0U; index--) {
		unsigned int const previous = index - 1U;

		previous_solution = right_side[previous] - superdiagonal[previous] * previous_solution;

		x[previous] = previous_solution;
	}
}

/* The cyclic system, by Sherman-Morrison. a[0] and c[count - 1] link the first and the last equations. superdiagonal,
   right_side, update and diagonal are scratch of count doubles. */
static void s9_solve_tridiagonal_cyclic(unsigned int count, double const* a, double const* b, double const* c, double const* d, double* x,
	double* superdiagonal, double* right_side, double* update, double* diagonal)
{
	unsigned int const last = count - 1U;
	double first_subdiagonal, last_superdiagonal, factor;
	unsigned int index;

	/* with two equations each cyclic link joins the same two unknowns as the diagonal next to it, so the links fold into the off-diagonals */
	if (count == 2U) {
		double folded_subdiagonal[2];
		double folded_superdiagonal[2];

		folded_subdiagonal[0] = 0.0;
		folded_subdiagonal[1] = a[1] + c[1];
		folded_superdiagonal[0] = a[0] + c[0];
		folded_superdiagonal[1] = 0.0;

		s9_solve_tridiagonal(2U, folded_subdiagonal, b, folded_superdiagonal, d, x, superdiagonal, right_side);

		return;
	}

	first_subdiagonal = a[0];
	last_superdiagonal = c[last];

	/* without cyclic links the system is plain tridiagonal */
	if (fabs(first_subdiagonal) < DBL_EPSILON && fabs(last_superdiagonal) < DBL_EPSILON) {
		s9_solve_tridiagonal(count, a, b, c, d, x, superdiagonal, right_side);

		return;
	}

	/* the cyclic links are a rank one update of a tridiagonal matrix with a modified diagonal: the update vector is solved in place of a right side, and the two solutions combine */
	memcpy(diagonal, b, sizeof(double) * count);

	diagonal[0] -= first_subdiagonal;
	diagonal[last] -= last_superdiagonal;

	memset(update, 0, sizeof(double) * count);

	update[0] = first_subdiagonal;
	update[last] = last_superdiagonal;

	s9_solve_tridiagonal(count, a, diagonal, c, update, update, superdiagonal, right_side);
	s9_solve_tridiagonal(count, a, diagonal, c, d, x, superdiagonal, right_side);

	factor = (x[0] + x[last]) / (1.0 + update[0] + update[last]);

	for (index = 0U; index < count; index++) {
		x[index] -= factor * update[index];
	}
}

/* Newton steps on all knot angles at once until the angles do not move more than S9_MAXIMUM_ITP_ERROR. The curvature delta
   of a knot depends only on its angle and its two neighbours', so the Jacobian is tridiagonal, cyclic for a closed spline;
   its columns come from finite differences. Each step is scaled so that no knot moves more than the range fraction of its
   distance to an end of its range, then halved until the greatest curvature delta goes down. The iteration gives up when
   a segment cannot be solved, when the step cannot be found, when no halving makes the deltas go down, or when the cap
   scales the steps down for several iterations in a row; the knot angles then keep the last accepted iterate and the
   segment curves match them. Convergence is judged on the Newton step before the cap and the clamp, so a knot that cannot
   move does not look converged. True when the iteration converged. */
static int s9_solve_knot_system_coupled(s9_knot_system const* system)
{
	unsigned int const count = system->count;
	unsigned int const segment_count = system->segment_count;
	unsigned int const unknown_count = s9_unknown_count(system);
	double* const start_curvatures = system->work[S9_WORK_START_CURVATURES];
	double* const end_curvatures = system->work[S9_WORK_END_CURVATURES];
	double* const trial_start_curvatures = system->work[S9_WORK_TRIAL_START_CURVATURES];
	double* const trial_end_curvatures = system->work[S9_WORK_TRIAL_END_CURVATURES];
	double* const deltas = system->work[S9_WORK_DELTAS];
	double* const trial_deltas = system->work[S9_WORK_TRIAL_DELTAS];
	double* const subdiagonal = system->work[S9_WORK_SUBDIAGONAL];
	double* const diagonal = system->work[S9_WORK_DIAGONAL];
	double* const superdiagonal = system->work[S9_WORK_SUPERDIAGONAL];
	double* const right_side = system->work[S9_WORK_RIGHT_SIDE];
	double* const step = system->work[S9_WORK_STEP];
	double* const lows = system->work[S9_WORK_LOWS];
	double* const highs = system->work[S9_WORK_HIGHS];
	double* const trial_angles = system->work[S9_WORK_TRIAL_ANGLES];
	double* const angles = system->knot_angles;
	int is_complete = 0;
	int is_evaluated;
	unsigned int stall_count = 0U;
	unsigned int iteration = 0U;

	if (unknown_count == 0U) {
		return 1;
	}

	s9_classify_knots(system, angles);
	s9_clamp_free_knots(system, angles);
	s9_centre_pinned_knots(system, angles);

	is_evaluated = s9_evaluate_knot_system(system, angles, 1, start_curvatures, end_curvatures, deltas);

	while (is_evaluated && iteration < S9_COUPLED_GREATEST_ITERATION_COUNT) {
		double delta_norm, scale_factor, greatest_step, factor, greatest_update;
		int columns_are_valid, step_is_finite, is_accepted;
		unsigned int index, unknown, halving, segment;

		iteration++;

		/* the ranges move with the knots, so the marks are refreshed */
		if (iteration > 1U) {
			if (s9_classify_knots(system, angles)) {
				s9_clamp_free_knots(system, angles);
				s9_centre_pinned_knots(system, angles);

				is_evaluated = s9_evaluate_knot_system(system, angles, 1, start_curvatures, end_curvatures, deltas);

				if (!is_evaluated) {
					break;
				}
			}
		}

		delta_norm = s9_free_delta_norm(system, unknown_count, deltas);

		for (index = 0U; index < unknown_count; index++) {
			subdiagonal[index] = 0.0;
			diagonal[index] = 0.0;
			superdiagonal[index] = 0.0;
		}

		/* finite difference columns, one unknown at a time */
		columns_are_valid = 1;

		for (unknown = 0U; unknown < unknown_count; unknown++) {
			unsigned int const knot = s9_knot_of_unknown(system, unknown);
			unsigned int const before = s9_segment_before_knot(system, knot);
			unsigned int const after = knot;
			double difference, scale, before_start_slope, before_end_slope, after_start_slope, after_end_slope;
			double before_start = 0.0;
			double before_end = 0.0;
			double after_start = 0.0;
			double after_end = 0.0;

			if (!s9_calculate_knot_range(system, angles, knot, &lows[unknown], &highs[unknown])) {
				columns_are_valid = 0;
				break;
			}

			for (index = 0U; index < count; index++) {
				trial_angles[index] = angles[index];
			}

			difference = S9_COUPLED_DIFFERENCE_STEP;

			if (angles[knot] + difference > highs[unknown]) {
				difference = -difference;
			}

			trial_angles[knot] = angles[knot] + difference;

			s9_apply_end_rules(system, trial_angles);

			if (!s9_solve_segments_at_knot(system, trial_angles, knot, &before_start, &before_end, &after_start, &after_end)) {
				columns_are_valid = 0;
				break;
			}

			scale = 1.0 / difference;
			before_start_slope = (before_start - start_curvatures[before]) * scale;
			before_end_slope = (before_end - end_curvatures[before]) * scale;
			after_start_slope = (after_start - start_curvatures[after]) * scale;
			after_end_slope = (after_end - end_curvatures[after]) * scale;

			/* the delta of this knot: start of the segment after it minus end of the segment before it */
			diagonal[unknown] = after_start_slope - before_end_slope;

			/* the delta of the previous knot: the start of the segment before this knot changes */
			if (unknown > 0U) {
				superdiagonal[unknown - 1U] = before_start_slope;
			} else if (system->is_cyclic) {
				superdiagonal[unknown_count - 1U] = before_start_slope;
			}

			/* the delta of the next knot: the end of the segment after this knot changes */
			if (unknown + 1U < unknown_count) {
				subdiagonal[unknown + 1U] = -after_end_slope;
			} else if (system->is_cyclic) {
				subdiagonal[0] = -after_end_slope;
			}
		}

		if (!columns_are_valid) {
			break;
		}

		for (index = 0U; index < unknown_count; index++) {
			right_side[index] = -deltas[index];
		}

		/* a pinned knot has no delta to zero: its row keeps it in the middle of its range as the neighbours move */
		for (unknown = 0U; unknown < unknown_count; unknown++) {
			unsigned int const knot = s9_knot_of_unknown(system, unknown);

			if (system->pinned[knot]) {
				s9_write_pinned_row(system, unknown_count, unknown, angles, trial_angles, subdiagonal, diagonal, superdiagonal, right_side);
			}
		}

		if (unknown_count == 1U) {
			if (diagonal[0] == 0.0) {
				break;
			}

			step[0] = right_side[0] / diagonal[0];
		} else if (system->is_cyclic) {
			s9_solve_tridiagonal_cyclic(unknown_count, subdiagonal, diagonal, superdiagonal, right_side, step,
				system->work[S9_WORK_SOLVE_SUPERDIAGONAL], system->work[S9_WORK_SOLVE_RIGHT_SIDE], system->work[S9_WORK_SOLVE_UPDATE],
				system->work[S9_WORK_SOLVE_DIAGONAL]);
		} else {
			s9_solve_tridiagonal(unknown_count, subdiagonal, diagonal, superdiagonal, right_side, step,
				system->work[S9_WORK_SOLVE_SUPERDIAGONAL], system->work[S9_WORK_SOLVE_RIGHT_SIDE]);
		}

		/* cap the step of the free knots: the whole step is scaled, so its direction is kept */
		scale_factor = 1.0;
		greatest_step = 0.0;
		step_is_finite = 1;

		for (unknown = 0U; unknown < unknown_count; unknown++) {
			unsigned int const knot = s9_knot_of_unknown(system, unknown);
			double const value = step[unknown];
			double step_magnitude, room, positive_room, allowed, magnitude;

			if (!isfinite(value)) {
				step_is_finite = 0;
				break;
			}

			/* a pinned knot follows, so its step is not capped */
			if (system->pinned[knot]) {
				continue;
			}

			step_magnitude = s9_abs(value);

			greatest_step = s9_max(greatest_step, step_magnitude);

			room = highs[unknown] - angles[knot];

			if (value < 0.0) {
				room = angles[knot] - lows[unknown];
			}

			positive_room = s9_max(room, 0.0);
			allowed = S9_COUPLED_RANGE_FRACTION * positive_room + S9_MAXIMUM_ITP_ERROR;
			magnitude = s9_abs(value);

			if (magnitude > allowed) {
				scale_factor = s9_min(scale_factor, allowed / magnitude);
			}
		}

		if (!step_is_finite) {
			break;
		}

		if (scale_factor < S9_COUPLED_STALL_FACTOR) {
			stall_count++;
		} else {
			stall_count = 0U;
		}

		if (stall_count >= S9_COUPLED_GREATEST_STALL_COUNT) {
			break;
		}

		/* line search on the greatest delta */
		factor = scale_factor;
		is_accepted = 0;
		greatest_update = 0.0;

		for (halving = 0U; halving <= S9_COUPLED_GREATEST_HALVING_COUNT; halving++) {
			for (index = 0U; index < count; index++) {
				trial_angles[index] = angles[index];
			}

			for (unknown = 0U; unknown < unknown_count; unknown++) {
				unsigned int const knot = s9_knot_of_unknown(system, unknown);
				double const moved = angles[knot] + factor * step[unknown];

				trial_angles[knot] = moved;

				if (!system->pinned[knot]) {
					trial_angles[knot] = s9_clamp(moved, lows[unknown], highs[unknown]);
				}
			}

			/* pinned knots take their exact middles, not the linear ones */
			s9_centre_pinned_knots(system, trial_angles);

			greatest_update = 0.0;

			for (index = 0U; index < count; index++) {
				double const magnitude = s9_abs(trial_angles[index] - angles[index]);

				greatest_update = s9_max(greatest_update, magnitude);
			}

			if (s9_evaluate_knot_system(system, trial_angles, 1, trial_start_curvatures, trial_end_curvatures, trial_deltas)) {
				double const trial_norm = s9_free_delta_norm(system, unknown_count, trial_deltas);
				int const decreases = trial_norm <= delta_norm * (1.0 - S9_COUPLED_DECREASE * factor);

				/* the Newton step itself must be tiny, not the movement after the cap and the clamp, or a held knot looks converged */
				int const is_tiny = greatest_step < S9_MAXIMUM_ITP_ERROR;

				if (decreases || is_tiny) {
					is_accepted = 1;
					break;
				}
			}

			factor *= 0.5;
		}

		if (!is_accepted) {
			break;
		}

		for (index = 0U; index < count; index++) {
			angles[index] = trial_angles[index];
		}

		for (segment = 0U; segment < segment_count; segment++) {
			start_curvatures[segment] = trial_start_curvatures[segment];
			end_curvatures[segment] = trial_end_curvatures[segment];
		}

		for (unknown = 0U; unknown < unknown_count; unknown++) {
			deltas[unknown] = trial_deltas[unknown];
		}

		if (greatest_step < S9_MAXIMUM_ITP_ERROR) {
			is_complete = 1;
			break;
		}
	}

	/* the segment curves can hold a rejected trial: match them to the accepted angles */
	if (is_evaluated) {
		s9_evaluate_knot_system(system, angles, 1, start_curvatures, end_curvatures, deltas);
	}

	return is_complete;
}

/* Solves the knot angles so that the curvature is the same on both sides of every knot: sweeps of Newton steps until the
   greatest update of a sweep is under S9_COUPLED_START_ERROR, then the coupled iteration. True when the angles converged.
   False when a step failed: the angles then hold the last accepted iterate, a better start for the sweeps of
   s9_solve_knot_angle than the initial angles. */
static int s9_solve_knot_system(s9_knot_system const* system)
{
	double* const trial = system->work[S9_WORK_TRIAL];
	double greatest_update = 1e300;
	unsigned int sweep = 0U;
	int is_complete = 0;
	int is_failed = 0;

	s9_solve_unknown_segments(system);

	while (sweep < S9_KNOT_SWEEP_COUNT && greatest_update >= S9_COUPLED_START_ERROR) {
		greatest_update = s9_sweep_knot_newton_steps(system, trial);
		sweep++;

		if (greatest_update < 0.0) {
			is_failed = 1;
			break;
		}

		if (greatest_update < S9_MAXIMUM_ITP_ERROR) {
			is_complete = 1;
			break;
		}
	}

	if (is_failed) {
		return 0;
	}

	if (is_complete) {
		return 1;
	}

	return s9_solve_knot_system_coupled(system);
}

static s9_point s9_knot(double const* points, unsigned int index)
{
	return s9_xy(points[2U * index], points[2U * index + 1U]);
}

/* the sweeps of s9_solve_knot_angle over an open spline, when the Newton steps did not converge */
static void s9_sweep_open(s9_knot_system const* system, signed char* low_signs)
{
	unsigned int const n = system->count - 1U;
	double* const knot_angles = system->knot_angles;
	double const* const segment_angles = system->segment_angles;
	double const* const segment_lengths = system->segment_lengths;
	double* const previous_updates = system->work[S9_WORK_PREVIOUS_UPDATES];
	s9_segment_curve* const segment_curves = system->segment_curves;
	unsigned int i, j;

	for (i = 0U; i < S9_MAXIMUM_ITP_ITERATION_COUNT; i++) {
		double angle0, angle1, max_error;

		s9_classify_knots(system, knot_angles);

		angle0 = s9_start_knot_angle(system, knot_angles);
		max_error = s9_abs(knot_angles[0] - angle0);

		knot_angles[0] = angle0;

		for (j = 1U; j < n; j++) {
			double const prev_knot_angle = knot_angles[j - 1U];
			double const next_knot_angle = knot_angles[j + 1U];
			double const prev_segment_angle = segment_angles[j - 1U];
			double const next_segment_angle = segment_angles[j];
			double const prev_segment_length = segment_lengths[j - 1U];
			double const next_segment_length = segment_lengths[j];
			double solved_angle, angle;

			if (system->pinned[j]) {
				solved_angle = s9_clamp_knot_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle);
			} else {
				solved_angle = s9_solve_knot_angle(prev_knot_angle, knot_angles[j], next_knot_angle, prev_segment_angle, next_segment_angle,
					prev_segment_length, next_segment_length, previous_updates[j], &segment_curves[j - 1U], &segment_curves[j], &low_signs[j]);
			}

			angle = s9_calculate_damped_knot_angle(knot_angles[j], solved_angle, &previous_updates[j]);

			max_error = s9_max(max_error, s9_abs(knot_angles[j] - angle));

			knot_angles[j] = angle;
		}

		angle1 = s9_end_knot_angle(system, knot_angles);

		max_error = s9_max(max_error, s9_abs(knot_angles[n] - angle1));

		knot_angles[n] = angle1;

		if (max_error < S9_MAXIMUM_ITP_ERROR) {
			break;
		}
	}
}

/* the sweeps of s9_solve_knot_angle over a closed spline, when the Newton steps did not converge */
static void s9_sweep_closed(s9_knot_system const* system, signed char* low_signs)
{
	unsigned int const count = system->count;
	unsigned int const final_knot = count - 1U;
	unsigned int const final_segment = count - 1U;
	double* const knot_angles = system->knot_angles;
	double const* const segment_angles = system->segment_angles;
	double const* const segment_lengths = system->segment_lengths;
	double* const previous_updates = system->work[S9_WORK_PREVIOUS_UPDATES];
	s9_segment_curve* const segment_curves = system->segment_curves;
	unsigned int i, j;

	for (i = 0U; i < S9_MAXIMUM_ITP_ITERATION_COUNT; i++) {
		double solved_angle0, angle0, solved_angle1, angle1, max_error;

		s9_classify_knots(system, knot_angles);

		solved_angle0 = s9_solve_knot_angle(knot_angles[final_knot], knot_angles[0], knot_angles[1], segment_angles[final_segment],
			segment_angles[0], segment_lengths[final_segment], segment_lengths[0], previous_updates[0], &segment_curves[final_segment],
			&segment_curves[0], &low_signs[0]);

		angle0 = s9_calculate_damped_knot_angle(knot_angles[0], solved_angle0, &previous_updates[0]);

		max_error = s9_abs(knot_angles[0] - angle0);

		knot_angles[0] = angle0;

		for (j = 1U; j < count - 1U; j++) {
			double const prev_knot_angle = knot_angles[j - 1U];
			double const next_knot_angle = knot_angles[j + 1U];
			double const prev_segment_angle = segment_angles[j - 1U];
			double const next_segment_angle = segment_angles[j];
			double const prev_segment_length = segment_lengths[j - 1U];
			double const next_segment_length = segment_lengths[j];
			double solved_angle, angle;

			if (system->pinned[j]) {
				solved_angle = s9_clamp_knot_angle(prev_knot_angle, next_knot_angle, prev_segment_angle, next_segment_angle);
			} else {
				solved_angle = s9_solve_knot_angle(prev_knot_angle, knot_angles[j], next_knot_angle, prev_segment_angle, next_segment_angle,
					prev_segment_length, next_segment_length, previous_updates[j], &segment_curves[j - 1U], &segment_curves[j], &low_signs[j]);
			}

			angle = s9_calculate_damped_knot_angle(knot_angles[j], solved_angle, &previous_updates[j]);

			max_error = s9_max(max_error, s9_abs(knot_angles[j] - angle));

			knot_angles[j] = angle;
		}

		solved_angle1 = s9_solve_knot_angle(knot_angles[final_knot - 1U], knot_angles[final_knot], knot_angles[0], segment_angles[final_segment - 1U],
			segment_angles[final_segment], segment_lengths[final_segment - 1U], segment_lengths[final_segment], previous_updates[final_knot],
			&segment_curves[final_segment - 1U], &segment_curves[final_segment], &low_signs[final_knot]);

		angle1 = s9_calculate_damped_knot_angle(knot_angles[final_knot], solved_angle1, &previous_updates[final_knot]);

		max_error = s9_max(max_error, s9_abs(knot_angles[final_knot] - angle1));

		knot_angles[final_knot] = angle1;

		if (max_error < S9_MAXIMUM_ITP_ERROR) {
			break;
		}
	}
}

/* the starting angles of the knots: each interior knot leans by the chord lengths, the ends follow, then every knot is moved to the middle of its range */
static void s9_initial_knot_angles(s9_knot_system const* system)
{
	unsigned int const count = system->count;
	unsigned int const n = count - 1U;
	double* const knot_angles = system->knot_angles;
	double const* const segment_angles = system->segment_angles;
	double const* const segment_lengths = system->segment_lengths;
	unsigned int i, j;

	for (i = 1U; i < n; i++) {
		double const a0 = segment_angles[i - 1U];
		double const a1 = segment_angles[i];
		double const l0 = segment_lengths[i - 1U];
		double const l1 = segment_lengths[i];
		double const combined_lengths = l0 + l1;
		double const bend = a1 - a0;

		knot_angles[i] = a0 + ((bend * l0) / combined_lengths);
	}

	knot_angles[0] = s9_alpha_to_absolute(0.0, segment_angles[0]);
	knot_angles[n] = s9_beta_to_absolute(0.0, segment_angles[system->segment_count - 1U]);

	for (i = 0U; i < S9_MAXIMUM_CENTERING_ITERATION_COUNT; i++) {
		if (system->is_cyclic) {
			knot_angles[0] = s9_clamp_knot_angle(knot_angles[n], knot_angles[1], segment_angles[n], segment_angles[0]);
		} else {
			knot_angles[0] = s9_start_knot_angle(system, knot_angles);
		}

		for (j = 1U; j < n; j++) {
			knot_angles[j] = s9_clamp_knot_angle(knot_angles[j - 1U], knot_angles[j + 1U], segment_angles[j - 1U], segment_angles[j]);
		}

		if (system->is_cyclic) {
			knot_angles[n] = s9_clamp_knot_angle(knot_angles[n - 1U], knot_angles[0], segment_angles[n - 1U], segment_angles[n]);
		} else {
			knot_angles[n] = s9_end_knot_angle(system, knot_angles);
		}
	}
}

qaws_status qaws_internal_interpolate_s9(double** out_points, unsigned int* out_segments, double const* points, unsigned int n, int closed)
{
	unsigned int const segment_count = closed ? n : n - 1U;
	double* work = NULL;
	s9_segment_curve* segment_curves = NULL;
	signed char* flags = NULL;
	s9_point* piece_points = NULL;
	s9_output output;
	s9_knot_system system;
	s9_pieces pieces;
	qaws_status status = QAWS_STATUS_ALLOCATION_FAILURE;
	unsigned int i;
	size_t segments, total, k;

	output.points = NULL;
	output.point_count = 0U;
	output.capacity = 64U;

	work = (double*)qaws_internal_malloc(qaws_internal_mul3(n, S9_WORK_COUNT, sizeof(double)));
	segment_curves = (s9_segment_curve*)qaws_internal_malloc(qaws_internal_mul2(n, sizeof(s9_segment_curve)));
	flags = (signed char*)qaws_internal_malloc(qaws_internal_mul2(n, 2U));
	piece_points = (s9_point*)qaws_internal_malloc(qaws_internal_mul2(S9_PIECE_POINT_COUNT, sizeof(s9_point)));
	output.points = (double*)qaws_internal_malloc(qaws_internal_mul3(output.capacity, 2U, sizeof(double)));

	if (work == NULL || segment_curves == NULL || flags == NULL || piece_points == NULL || output.points == NULL) {
		goto done;
	}

	for (i = 0U; i < (unsigned int)S9_WORK_COUNT; i++) {
		system.work[i] = work + (size_t)i * n;
	}

	system.count = n;
	system.segment_count = segment_count;
	system.is_cyclic = closed;
	system.segment_angles = system.work[S9_WORK_SEGMENT_ANGLES];
	system.segment_lengths = system.work[S9_WORK_SEGMENT_LENGTHS];
	system.knot_angles = system.work[S9_WORK_KNOT_ANGLES];
	system.segment_curves = segment_curves;
	system.pinned = flags;

	memset(flags, 0, (size_t)n * 2U);

	for (i = 0U; i < n; i++) {
		segment_curves[i].shape = s9_curve_default();
		segment_curves[i].is_known = 0;
		system.work[S9_WORK_PREVIOUS_UPDATES][i] = 0.0;
	}

	for (i = 0U; i < segment_count; i++) {
		s9_point const d = s9_sub(s9_knot(points, s9_end_knot_of_segment(&system, i)), s9_knot(points, i));

		system.work[S9_WORK_SEGMENT_LENGTHS][i] = sqrt(d.x * d.x + d.y * d.y);
		system.work[S9_WORK_SEGMENT_ANGLES][i] = s9_vector_angle(d);
	}

	s9_initial_knot_angles(&system);

	if (!s9_solve_knot_system(&system)) {
		if (closed) {
			s9_sweep_closed(&system, flags + n);
		} else {
			s9_sweep_open(&system, flags + n);
		}
	}

	/* the chain starts at the first point, then each segment's cubics */
	output.points[0] = points[0];
	output.points[1] = points[1];
	output.point_count = 1U;

	pieces.points = piece_points;
	pieces.count = 0U;

	for (i = 0U; i < segment_count; i++) {
		unsigned int const end_knot = s9_end_knot_of_segment(&system, i);
		double const segment_angle = system.segment_angles[i];
		double const a0 = s9_absolute_to_alpha(system.knot_angles[i], segment_angle);
		double const a1 = s9_absolute_to_beta(system.knot_angles[end_knot], segment_angle);

		status = s9_append_cubics(a0, a1, s9_knot(points, i), s9_knot(points, end_knot), segment_angle, system.segment_lengths[i],
			&segment_curves[i], &pieces, &output);

		if (status != QAWS_STATUS_OK) {
			goto done;
		}
	}

	/* the result must have at least one cubic and only finite coordinates */
	status = QAWS_STATUS_NUMERICAL_FAILURE;

	if (output.point_count < 4U) {
		goto done;
	}

	segments = (output.point_count - 1U) / 3U;

	if (segments > (size_t)UINT_MAX) {
		status = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}

	total = 2U * output.point_count;

	for (k = 0U; k < total; k++) {
		if (!isfinite(output.points[k])) {
			goto done;
		}
	}

	*out_points = output.points;
	*out_segments = (unsigned int)segments;
	output.points = NULL;
	status = QAWS_STATUS_OK;

done:
	qaws_internal_free(output.points);
	qaws_internal_free(piece_points);
	qaws_internal_free(flags);
	qaws_internal_free(segment_curves);
	qaws_internal_free(work);

	return status;
}
