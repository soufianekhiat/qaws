#ifndef QAWS_BSPLINE_EVAL_CORE_H
#define QAWS_BSPLINE_EVAL_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"

/* ===================================================================
 * B-spline evaluation (2D)
 *
 * local_cp: (degree+1) control points, flat [x0,y0, x1,y1, ...]
 * local_knots: relevant knot window (2*(degree+1) values)
 * span: knot span index within local_knots
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_bspline_eval_2d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int max_deriv = 0;
    int stride = degree + 1;
    int j;
    qaws_eval_2d result;
    qaws_scalar vx, vy, basis;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO;

    if ((eval_flags & 0x8) != 0) max_deriv = 3;
    else if ((eval_flags & 0x4) != 0) max_deriv = 2;
    else if ((eval_flags & 0x2) != 0) max_deriv = 1;
    if (max_deriv > degree) max_deriv = degree;

    qaws_bspline_basis_derivs(local_knots, degree, span, t, max_deriv, ders);

    /* Position */
    vx = QAWS_ZERO; vy = QAWS_ZERO;
    for (j = 0; j <= degree; j++) {
        basis = ders[0 * stride + j];
        vx += basis * local_cp[j * 2 + 0];
        vy += basis * local_cp[j * 2 + 1];
    }
    result.position.x = vx;
    result.position.y = vy;

    /* D1 */
    if (max_deriv >= 1) {
        vx = QAWS_ZERO; vy = QAWS_ZERO;
        for (j = 0; j <= degree; j++) {
            basis = ders[1 * stride + j];
            vx += basis * local_cp[j * 2 + 0];
            vy += basis * local_cp[j * 2 + 1];
        }
        result.d1.x = vx;
        result.d1.y = vy;
    }

    /* D2 */
    if (max_deriv >= 2) {
        vx = QAWS_ZERO; vy = QAWS_ZERO;
        for (j = 0; j <= degree; j++) {
            basis = ders[2 * stride + j];
            vx += basis * local_cp[j * 2 + 0];
            vy += basis * local_cp[j * 2 + 1];
        }
        result.d2.x = vx;
        result.d2.y = vy;
    }

    return result;
}

/* ===================================================================
 * B-spline evaluation (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_bspline_eval_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int max_deriv = 0;
    int stride = degree + 1;
    int j;
    qaws_eval_3d result;
    qaws_scalar vx, vy, vz, basis;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO; result.position.z = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO; result.d1.z = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO; result.d2.z = QAWS_ZERO;

    if ((eval_flags & 0x8) != 0) max_deriv = 3;
    else if ((eval_flags & 0x4) != 0) max_deriv = 2;
    else if ((eval_flags & 0x2) != 0) max_deriv = 1;
    if (max_deriv > degree) max_deriv = degree;

    qaws_bspline_basis_derivs(local_knots, degree, span, t, max_deriv, ders);

    /* Position */
    vx = QAWS_ZERO; vy = QAWS_ZERO; vz = QAWS_ZERO;
    for (j = 0; j <= degree; j++) {
        basis = ders[0 * stride + j];
        vx += basis * local_cp[j * 3 + 0];
        vy += basis * local_cp[j * 3 + 1];
        vz += basis * local_cp[j * 3 + 2];
    }
    result.position.x = vx;
    result.position.y = vy;
    result.position.z = vz;

    /* D1 */
    if (max_deriv >= 1) {
        vx = QAWS_ZERO; vy = QAWS_ZERO; vz = QAWS_ZERO;
        for (j = 0; j <= degree; j++) {
            basis = ders[1 * stride + j];
            vx += basis * local_cp[j * 3 + 0];
            vy += basis * local_cp[j * 3 + 1];
            vz += basis * local_cp[j * 3 + 2];
        }
        result.d1.x = vx;
        result.d1.y = vy;
        result.d1.z = vz;
    }

    /* D2 */
    if (max_deriv >= 2) {
        vx = QAWS_ZERO; vy = QAWS_ZERO; vz = QAWS_ZERO;
        for (j = 0; j <= degree; j++) {
            basis = ders[2 * stride + j];
            vx += basis * local_cp[j * 3 + 0];
            vy += basis * local_cp[j * 3 + 1];
            vz += basis * local_cp[j * 3 + 2];
        }
        result.d2.x = vx;
        result.d2.y = vy;
        result.d2.z = vz;
    }

    return result;
}

#endif /* QAWS_BSPLINE_EVAL_CORE_H */
