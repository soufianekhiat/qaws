#include "../qaws_platform.h"
#include "../qaws_types.h"
#include "../qaws_core_types.h"
#include "../internal/qaws_internal_basis.h"

/* ===================================================================
 * Batch NURBS evaluation (2D)
 *
 * Evaluates a NURBS curve at N parameter values using the quotient
 * rule for rational B-splines. For each parameter value the knot
 * span is located, basis functions are computed, and the weighted
 * rational point is evaluated as:
 *
 *   C(t) = sum( N_i(t) * w_i * P_i ) / sum( N_i(t) * w_i )
 *
 * knots:           full knot vector
 * knot_count:      number of knots
 * control_points:  flat array [x0,y0, x1,y1, ...]
 * weights:         per-control-point weights
 * cp_count:        number of control points
 * degree:          spline degree
 * t_values:        array of count parameter values
 * count:           number of parameter values
 * out_positions:   receives count * 2 scalars [x0,y0, x1,y1, ...]
 * =================================================================== */

void qaws_batch_nurbs_eval_2d(
    qaws_scalar const* knots,
    unsigned int knot_count,
    qaws_scalar const* control_points,
    qaws_scalar const* weights,
    unsigned int cp_count,
    unsigned int degree,
    qaws_scalar const* t_values,
    unsigned int count,
    qaws_scalar* out_positions)
{
    unsigned int i, j;
    unsigned int order = degree + 1;

    for (i = 0; i < count; i++) {
        qaws_scalar t = t_values[i];
        unsigned int span;
        qaws_scalar basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar ax, ay, w_sum;

        span = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);
        qaws_internal_bspline_basis(basis, knots, knot_count, degree, span, t);

        ax = QAWS_ZERO;
        ay = QAWS_ZERO;
        w_sum = QAWS_ZERO;

        for (j = 0; j < order; j++) {
            unsigned int idx = span - degree + j;
            unsigned int ci = idx * 2;
            qaws_scalar nw = basis[j] * weights[idx];
            w_sum += nw;
            ax += nw * control_points[ci];
            ay += nw * control_points[ci + 1];
        }

        {
            qaws_scalar inv_w = QAWS_ONE / w_sum;
            out_positions[i * 2]     = ax * inv_w;
            out_positions[i * 2 + 1] = ay * inv_w;
        }
    }
}

/* ===================================================================
 * Batch NURBS evaluation (3D)
 *
 * control_points: flat array [x0,y0,z0, x1,y1,z1, ...]
 * out_positions:  receives count * 3 scalars [x0,y0,z0, ...]
 * =================================================================== */

void qaws_batch_nurbs_eval_3d(
    qaws_scalar const* knots,
    unsigned int knot_count,
    qaws_scalar const* control_points,
    qaws_scalar const* weights,
    unsigned int cp_count,
    unsigned int degree,
    qaws_scalar const* t_values,
    unsigned int count,
    qaws_scalar* out_positions)
{
    unsigned int i, j;
    unsigned int order = degree + 1;

    for (i = 0; i < count; i++) {
        qaws_scalar t = t_values[i];
        unsigned int span;
        qaws_scalar basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar ax, ay, az, w_sum;

        span = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);
        qaws_internal_bspline_basis(basis, knots, knot_count, degree, span, t);

        ax = QAWS_ZERO;
        ay = QAWS_ZERO;
        az = QAWS_ZERO;
        w_sum = QAWS_ZERO;

        for (j = 0; j < order; j++) {
            unsigned int idx = span - degree + j;
            unsigned int ci = idx * 3;
            qaws_scalar nw = basis[j] * weights[idx];
            w_sum += nw;
            ax += nw * control_points[ci];
            ay += nw * control_points[ci + 1];
            az += nw * control_points[ci + 2];
        }

        {
            qaws_scalar inv_w = QAWS_ONE / w_sum;
            out_positions[i * 3]     = ax * inv_w;
            out_positions[i * 3 + 1] = ay * inv_w;
            out_positions[i * 3 + 2] = az * inv_w;
        }
    }
}
