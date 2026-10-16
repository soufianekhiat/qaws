#include "../qaws_platform.h"
#include "../qaws_types.h"
#include "../qaws_core_types.h"
#include "../internal/qaws_internal_basis.h"

/* ===================================================================
 * Batch B-spline evaluation (2D)
 *
 * Evaluates a B-spline curve at N parameter values. For each value
 * the knot span is located and basis functions computed, then the
 * position is accumulated as a weighted sum of local control points.
 *
 * knots:           full knot vector
 * knot_count:      number of knots
 * control_points:  flat array [x0,y0, x1,y1, ...]
 * cp_count:        number of control points
 * degree:          spline degree
 * t_values:        array of count parameter values
 * count:           number of parameter values
 * out_positions:   receives count * 2 scalars [x0,y0, x1,y1, ...]
 * =================================================================== */

void qaws_batch_bspline_eval_2d(
    qaws_scalar const* knots,
    unsigned int knot_count,
    qaws_scalar const* control_points,
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
        qaws_scalar vx, vy;

        span = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);
        qaws_internal_bspline_basis(basis, knots, knot_count, degree, span, t);

        vx = QAWS_ZERO;
        vy = QAWS_ZERO;

        for (j = 0; j < order; j++) {
            unsigned int ci = (span - degree + j) * 2;
            qaws_scalar b = basis[j];
            vx += b * control_points[ci];
            vy += b * control_points[ci + 1];
        }

        out_positions[i * 2]     = vx;
        out_positions[i * 2 + 1] = vy;
    }
}

/* ===================================================================
 * Batch B-spline evaluation (3D)
 *
 * control_points: flat array [x0,y0,z0, x1,y1,z1, ...]
 * out_positions:  receives count * 3 scalars [x0,y0,z0, ...]
 * =================================================================== */

void qaws_batch_bspline_eval_3d(
    qaws_scalar const* knots,
    unsigned int knot_count,
    qaws_scalar const* control_points,
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
        qaws_scalar vx, vy, vz;

        span = qaws_internal_find_knot_span(knots, knot_count, degree, cp_count, t);
        qaws_internal_bspline_basis(basis, knots, knot_count, degree, span, t);

        vx = QAWS_ZERO;
        vy = QAWS_ZERO;
        vz = QAWS_ZERO;

        for (j = 0; j < order; j++) {
            unsigned int ci = (span - degree + j) * 3;
            qaws_scalar b = basis[j];
            vx += b * control_points[ci];
            vy += b * control_points[ci + 1];
            vz += b * control_points[ci + 2];
        }

        out_positions[i * 3]     = vx;
        out_positions[i * 3 + 1] = vy;
        out_positions[i * 3 + 2] = vz;
    }
}
